from __future__ import print_function

import copy

from smpspy.smps_read import StochasticModel

__author__ = 'vujanicr'
"""
Collection of inner problems, coded to fit the required oracle interface to function properly with the dual methods.

Oracles for the dual function should:
INPUT
    - lambda_k: (array of floats) the point at which the oracle is invoked
OUTPUT
    - x_k: (array of floats) the current inner solution, that is, x(lambda_k)
    - d_k: (float) dual function value
    - subgrad_lambda_k: (array of floats, same dimension as lambda_k) a valid subgradient at lambda_k
"""

# TODO: Move documentation related to interface in the abstract class definition.

import numpy as np
import gurobipy as gb
# from smps_read import *
from collections import defaultdict, OrderedDict
import multiprocessing as mp
# import pathos.multiprocessing as mp
import time
from os.path import expanduser
import sys
try:
    from ipyparallel import Client
except ImportError:
    print('ipyparallel needed for Elimination Based 2-Stage Programming solver.')


####################
# -- Parameters -- #
####################

if sys.platform == 'win32':
    HOME_PATH = expanduser("~")
else:
    HOME_PATH = '/tmp'

LOCAL_OPTIMIZATION_MIPGAP = 0.01  # default is 0.0001.
NUM_CORES = mp.cpu_count()
MINIMUM_N_VOTES = 2


# Abstract class that defines the interface
class InnerProblemInterface(object):
    def oracle(self, lambda_k):
        raise NotImplementedError()

    def projection_function(self, lambda_k):
        raise NotImplementedError()


class TwoStage_SMPS_InnerProblem(InnerProblemInterface):
    def __init__(self, filename='./data/2_sizes/sizes10', R=None):
        # generate stochastic model
        self.stoch_model = StochasticModel(filename)
        if self.stoch_model.nominal_model.ModelSense != 1:
            # TODO fix this.
            print('WARNING: local model is maximization. Should be transformed into minimization.')
        self.n_scenarios = self.stoch_model.scenarios.__len__()
        self.n_x = self.stoch_model.scenarios[self.stoch_model.scenarios.keys()[0]].c[1].__len__()
        self.n_y = self.stoch_model.scenarios[self.stoch_model.scenarios.keys()[0]].c[2].__len__()
        self.dimension = self.n_x * self.n_scenarios
        self.n_stages = self.stoch_model.scenarios[self.stoch_model.scenarios.keys()[0]].c.__len__()
        if self.n_stages > 2:
            print('This implementation only supports up to two-stage problems.')
        self.SCENARIO_TIME_LIMIT = np.ceil(600/self.n_scenarios)  # 5 minutes max in total

        # for record keeping
        self.instance_type = ''
        self.instance_subtype = ''
        self.instance_name = ''

        self.R = R  # used in softmax computation

    def _construct_local_problem(self, scenario, lambda_k):
        # TODO I generate the entire local model every iteration, even though the only thing changing is the objective
        # It would be faster to store the model once it's generated the first time and then to just change the objective.
        ##################################################################
        # local optimization program (inner problem) for single scenario #
        ##################################################################

        scenario_program = gb.Model()

        # -- variables -- #
        x = {}
        for i in range(self.n_x):
            x[i] = scenario_program.addVar(obj=scenario.c[1][i]*scenario.probability - lambda_k[i],
                                           lb=scenario.lb[1][i],
                                           ub=scenario.ub[1][i],
                                           vtype=scenario.vtype[1][i],
                                           name='x{}'.format(i))
        y = {}
        for i in range(self.n_y):
            y[i] = scenario_program.addVar(obj=scenario.c[2][i]*scenario.probability,
                                           lb=scenario.lb[2][i],
                                           ub=scenario.ub[2][i],
                                           vtype=scenario.vtype[2][i],
                                           name='y{}'.format(i))
        scenario_program.update()

        # -- constraints -- #
        # first stage (sparse matrix)
        lhs = defaultdict(gb.LinExpr)
        for entry in scenario.A[1,1].items():
            row, column = entry[0]
            value = entry[1]
            lhs[row] = lhs[row] + value*x[column]

        for row, value in enumerate(scenario.b[1]):
            scenario_program.addConstr(lhs[row], scenario.b_sense[1][row], scenario.b[1][row],
                                       name='stage1_row{}'.format(row))

        # constraints, second-stage; contribution to lhs from A[2,1]
        lhs_1 = defaultdict(gb.LinExpr)
        for entry in scenario.A[2,1].items():
            row, column = entry[0]
            value = entry[1]
            lhs_1[row] = lhs_1[row] + value*x[column]
        # contribution from A[2,2]
        lhs_2 = defaultdict(gb.LinExpr)
        for entry in scenario.A[2,2].items():
            row, column = entry[0]
            value = entry[1]
            lhs_2[row] = lhs_2[row] + value*y[column]

        for row, value in enumerate(scenario.b[2]):
                scenario_program.addConstr(lhs_1[row] + lhs_2[row],
                                           scenario.b_sense[2][row],
                                           scenario.b[2][row],
                                           name='stage2_sc{}_row{}'.format(scenario.scenario_id,row))
        scenario_program.update()
        return scenario_program, x, y

    def _calculate_oracle_outputs(self, scenario_objectives, scenario_optimizers):
        x_k = [item for sublist in scenario_optimizers for item in sublist]  # flatten list
        x_k = np.array(x_k)  # return as np.array (NOT list!)
        d_k = sum(scenario_objectives)

        # subgradient: calculate average, subgradient is then difference between local optimizers and the avg
        average_optimizer = np.average(scenario_optimizers, axis=0)

        subgrad_lambda_k = []
        for s in range(self.n_scenarios):
            subgrad_lambda_k.append(average_optimizer - scenario_optimizers[s])
        # flatten and return as np.array (NOT list!)
        subgrad_lambda_k = [item for sublist in subgrad_lambda_k for item in sublist]
        subgrad_lambda_k = np.array(subgrad_lambda_k)

        return x_k, d_k, subgrad_lambda_k

    def oracle(self, lambda_k):
        scenario_objectives = []
        scenario_optimizers = []

        for scenario, scenario_label in enumerate(self.stoch_model.scenarios):
            # Attention here. Since we use an OrderedDict to contain the scenarios in StochasticModel, the  order in
            # which scenarios are iterated over is consistent. This is crucial to avoid bugs.
            scenario_lambda_k = lambda_k[scenario*self.n_x:(scenario+1)*self.n_x]
            scenario_program, x, y = self._construct_local_problem(self.stoch_model.scenarios[scenario_label],
                                                                   scenario_lambda_k)

            scenario_program.setParam('OutputFlag', False)
            scenario_program.setParam('MIPGap', LOCAL_OPTIMIZATION_MIPGAP)
            scenario_program.setParam('TimeLimit', self.SCENARIO_TIME_LIMIT)
            scenario_program.optimize()
            if scenario_program.Status != gb.GRB.OPTIMAL:
                if scenario_program.Status == gb.GRB.TIME_LIMIT:
                    # we have stopped because of time limit; check that a feasible solution was found
                    assert scenario_program.SolCount > 0, 'Scenario program reached SCENARIO_TIME_LIMIT={} without ' \
                                                          'finding a feasible solution.'.format(self.SCENARIO_TIME_LIMIT)
                else:
                    # something else happened that caused an error
                    raise RuntimeError('Solution of single scenario encountered an error. '
                                       'Error code gb.GRB.Status = {}'.format(scenario_program.Status))



            # record local objective and optimizers
            scenario_objectives.append(scenario_program.ObjVal)
            temp = []
            for i in range(self.n_x):
                # TODO Only saving x here, what about y?
                temp.append(x[i].X)
            scenario_optimizers.append(temp)

        return self._calculate_oracle_outputs(scenario_objectives, scenario_optimizers)

    def projection_function(self, lambda_k):
        # reorganize lambda_k into a list of lists, each entry related to a scenario
        scenario_lambda_k = []
        for scenario in range(self.n_scenarios):
            scenario_lambda_k.append(lambda_k[scenario*self.n_x:(scenario+1)*self.n_x])

        # determine average lambda_k vector
        average_lambda_k = np.average(scenario_lambda_k, axis=0)

        # then remove average from each entry, to obtain a vector that sums to 0
        for scenario in range(self.n_scenarios):
            scenario_lambda_k[scenario] = scenario_lambda_k[scenario] - average_lambda_k

        # flatten again to a single list
        lambda_k = [item for sublist in scenario_lambda_k for item in sublist]
        lambda_k = np.array(lambda_k)
        return lambda_k

    def softmax_projection(self, lambda_k):
        assert self.R, 'inner problem must know the parameter R >= || lambda* - lambda_0 || to perform softmax projection'

        lambda_k = copy.copy(lambda_k)
        lambda_k = self.n_scenarios*self.R*lambda_k

        # reorganize lambda_k into a list of lists, each entry related to a scenario
        # INPUT: np.array([0.1, 0.2, 1.3, 0.2, 1.6, 0.5])
        # OUTPUT: [array([  0.1,   0.2,  1.3]), array([  0.2,  1.6,  0.5])]
        scenario_lambda_k = []
        for scenario in range(self.n_scenarios):
            scenario_lambda_k.append(lambda_k[scenario*self.n_x:(scenario+1)*self.n_x])

        # determine max lambda_k vector
        max_lambda_k = np.max(scenario_lambda_k, axis=0)

        # then remove max from each entry
        for scenario in range(self.n_scenarios):
            scenario_lambda_k[scenario] = scenario_lambda_k[scenario] - max_lambda_k

        # compute exponentiation and normalizing factor
        normalizing_factor = np.zeros(len(scenario_lambda_k[0]))
        for scenario in range(self.n_scenarios):
            scenario_lambda_k[scenario] = np.exp(scenario_lambda_k[scenario])
            normalizing_factor += scenario_lambda_k[scenario]

        # then divide
        for scenario in range(self.n_scenarios):
            scenario_lambda_k[scenario] = scenario_lambda_k[scenario] / normalizing_factor

        # flatten again to a single list
        lambda_k = [item for sublist in scenario_lambda_k for item in sublist]
        lambda_k = self.R * (self.n_scenarios * np.array(lambda_k) - 1)

        return lambda_k


def _external_optimization_of_single_scenario(scenario_n, scenario_lambda_k, points_to_eliminate):
    """
    We cannot pickle Gurobi objects; declaring this as a function external from every class makes it work with the
    multiprocessing package (which doesn't like to work with instance methods for some reason).
    """
    # Read Model
    scenario_program = gb.read(HOME_PATH+'/sc_model_{}.mps'.format(scenario_n))
    # Reconstruct variables
    x = {}
    y = {}
    for variable_n, variable in enumerate(scenario_program.getVars()):
        if 'x' in variable.VarName:
            x[variable_n] = variable
        elif 'y' in variable.VarName:
            y[variable_n] = variable
        else:
            print('Cannot recognize variable name.')

    a = 2
    # Change objective according to lambda_k
    for i, slot in enumerate(x):
        x[slot].Obj = x[slot].Obj - scenario_lambda_k[i]

    # Add cutting planes
    for point in points_to_eliminate:
        balas_cut = map(lambda z: 2*z-1, point)  # hyperplane equation: [0,1,0] -> [-1,1,-1]
        scenario_program.addConstr(gb.quicksum([balas_cut[i]*x[i] for i in range(point.__len__())]) <= sum(point)-1,
                                   name='Balas cut of {0}'.format(str(point)))
    # Solve
    # scenario_program.setParam('OutputFlag', False)
    scenario_program.setParam('LogFile', '')
    scenario_program.setParam('MIPGap', LOCAL_OPTIMIZATION_MIPGAP)
    scenario_program.update()
    scenario_program.optimize()

    # TODO this assumes that variables appear in the correct order.
    x_sol = []
    y_sol = []
    for variable in scenario_program.getVars():
        if 'x' in variable.VarName:
            x_sol.append(variable.X)
        elif 'y' in variable.VarName:
            y_sol.append(variable.X)
        else:
            print('Cannot recognize variable name.')

    obj_value = scenario_program.ObjVal
    return tuple(x_sol), tuple(y_sol), obj_value


def _external_optimization_of_single_scenario_distributed(scenario_n, scenario_lambda_k, points_to_eliminate):
    """
    We cannot pickle Gurobi objects; declaring this as a function external from every class makes it work with the
    multiprocessing package (which doesn't like to work with instance methods for some reason).
    """
    # Read Model
    scenario_program = gb.read(HOME_PATH+'/sc_model_{}.mps'.format(scenario_n))
    # Reconstruct variables
    x = {}
    y = {}
    for variable_n, variable in enumerate(scenario_program.getVars()):
        if 'x' in variable.VarName:
            x[variable_n] = variable
        elif 'y' in variable.VarName:
            y[variable_n] = variable
        else:
            print('Cannot recognize variable name.')

    a = 2
    # Change objective according to lambda_k
    for i, slot in enumerate(x):
        x[slot].Obj = x[slot].Obj - scenario_lambda_k[i]

    # Add cutting planes
    for point in points_to_eliminate:
        balas_cut = map(lambda z: 2*z-1, point)  # hyperplane equation: [0,1,0] -> [-1,1,-1]
        scenario_program.addConstr(gb.quicksum([balas_cut[i]*x[i] for i in range(point.__len__())]) <= sum(point)-1,
                                   name='Balas cut of {0}'.format(str(point)))
    # Solve
    scenario_program.setParam('OutputFlag', False)
    scenario_program.setParam('LogFile', '')
    scenario_program.setParam('MIPGap', LOCAL_OPTIMIZATION_MIPGAP)
    scenario_program.update()
    scenario_program.optimize()

    # TODO this assumes that variables appear in the correct order.
    x_sol = []
    y_sol = []
    for variable in scenario_program.getVars():
        if 'x' in variable.VarName:
            x_sol.append(variable.X)
        elif 'y' in variable.VarName:
            y_sol.append(variable.X)
        else:
            print('Cannot recognize variable name.')

    obj_value = scenario_program.ObjVal
    return tuple(x_sol), tuple(y_sol), obj_value


class Elimination_TwoStage_SMPS_InnerProblem(TwoStage_SMPS_InnerProblem):
    # variation of two-stage with elimination of first stage solutions, if they are purely binary
    def __init__(self, filename='./data/2_sslp/sslp_5_25_3_mymod'):
        super(Elimination_TwoStage_SMPS_InnerProblem, self).__init__(filename)
        # TODO add a check that this is not run unless the problem is purely integer in the first stage.
        # introduce a list to record the points we want to cut
        self.points_to_eliminate = []  # [[0,1,1],[0,0,0]]
        self.best_current_primal_value = []
        self.with_elimination = 0
        self.stop_dual_iterations = 0
        self.seen_x_ks = []  # memory of inner solutions seen, used in the 'frequent' variant of elimination procedure

        self.time_generating = 0
        self.time_writing = 0

        # write local optimization problems for each scenario; from iteration to iteration, the only thing that will
        # change is the objective (and the cutting planes added to eliminate solutions).
        for scenario, scenario_label in enumerate(self.stoch_model.scenarios):
            scenario_program, x, y = self._construct_initial_local_problem(self.stoch_model.scenarios[scenario_label])
            scenario_program.write(HOME_PATH+'/sc_model_{}.mps'.format(str(scenario)))

        # For the distributed implementation, direct views of the engines; this requires ipcluster to be started
        try:
            self.rc = Client()
            self.dview = self.rc[:]
            self.dview.execute('import gurobipy as gb')

            @self.dview.parallel
            def temp(scenario_n, scenario_lambda_k, points_to_eliminate):
                return _external_optimization_of_single_scenario_distributed(scenario_n,
                                                                             scenario_lambda_k,
                                                                             points_to_eliminate)
            self._optimization_of_single_scenario_distributed = temp


        except:
            print('WARNING. Could not connect to the engines. This is likely due to ipcluster not being started.')
            print('Distributed computations will not be available.')

        self._optimization_of_single_scenario_distributed = 0

    def _construct_initial_local_problem(self, scenario):
        dummy_lambda_k = np.zeros(self.n_x)
        return self._construct_local_problem(scenario, dummy_lambda_k)

    def _elimination_procedure(self, scenario_objectives, scenario_optimizers, type='all'):
        """
        This procedure updates self.points_to_eliminate. Different variations are possible (type):
        * type='all': test the performance of all inner solutions, and remove all but the best.
        * type='popular': only test inner solutions that produced by at least 'MINIMUM_N_VOTES' subproblems.
        * type='frequent': only test inner solutions that have appeared at least once in some past iteration.
        """

        # We have all the inner solutions x_k (first stage part) inside scenario_optimizers; we evaluate how
        # they perform, and only keep the best one.
        unique_scenario_optimizers_to_evaluate = list(OrderedDict.fromkeys(scenario_optimizers))

        print('Number of unique inner solutions: ' + str(unique_scenario_optimizers_to_evaluate.__len__()))

        x_k_performance = []
        popular_x_ks = []
        frequent_x_ks = []

        # If the inner solutions x_k are all the same one point, we should stop: the method has found the optimizer.
        if unique_scenario_optimizers_to_evaluate.__len__() == 1:
            print('METHOD HAS CONVERGED, OPTIMIZER FOUND.')
            print('Optimal objective: '
                   + str(self._calculate_oracle_outputs(scenario_objectives, scenario_optimizers)[1]))
            print('Number of cutting planes added to the subproblems: {} '
                   '(over {})'.format(self.points_to_eliminate.__len__(), 2 ** self.n_x))
            self.stop_dual_iterations = 1

        elif self.with_elimination and type == 'all':
            for x_k in unique_scenario_optimizers_to_evaluate:
                x_k_performance.append(self._evaluate_inner_solution_performance(x_k))
            # find index of the best objective (min); this is the only point that we want to keep.
            self.best_current_primal_value.append(min(x_k_performance))
            best_x_k_index = x_k_performance.index(min(x_k_performance))
            # eliminate the corresponding point from list of unique scenarios
            del unique_scenario_optimizers_to_evaluate[best_x_k_index]
            # add list to points to be removed
            for x_k in unique_scenario_optimizers_to_evaluate:
                self.points_to_eliminate.append(x_k)

        elif self.with_elimination and type == 'popular':
            # This variation isn't very effective; from the experiments it seems that often we have one inner solution
            # that collects many votes, and all other only have one vote. So we do not eliminate anything for that
            # iteration.
            for x_k in unique_scenario_optimizers_to_evaluate:
                if scenario_optimizers.count(x_k) >= MINIMUM_N_VOTES:
                    popular_x_ks.append(x_k)
                    # this particular inner solution appears at least MINIMUM_N_VOTES times; we test it.
                    x_k_performance.append(self._evaluate_inner_solution_performance(x_k))
            # rest is the same as the 'all' case above
            print('[popular] Of which {} appear at least {} times in this iteration '
                   '(these are also the only ones tested).'.format(popular_x_ks.__len__(), MINIMUM_N_VOTES))
            self.best_current_primal_value.append(min(x_k_performance))
            best_x_k_index = x_k_performance.index(min(x_k_performance))
            del popular_x_ks[best_x_k_index]
            for x_k in popular_x_ks:
                self.points_to_eliminate.append(x_k)

        elif self.with_elimination and type == 'frequent':
            for x_k in unique_scenario_optimizers_to_evaluate:
                # is it has been seen, add it to the set to be tested
                if x_k in self.seen_x_ks:
                    frequent_x_ks.append(x_k)
                # if it hasnt been seen, mark it as seen now
                else:
                    self.seen_x_ks.append(x_k)

            print('[frequent] Of which {} have already appeared in some past iteration '
                   '(these are also the only ones tested).'.format(frequent_x_ks.__len__()))
            # remove all but the best from seen_x_ks AND from the inner problem (we need at least two to have a
            # meaningful elimination).
            if frequent_x_ks.__len__() >= 2:
                for x_k in frequent_x_ks:
                    x_k_performance.append(self._evaluate_inner_solution_performance(x_k))
                self.best_current_primal_value.append(min(x_k_performance))
                best_x_k_index = x_k_performance.index(min(x_k_performance))
                del frequent_x_ks[best_x_k_index]
                for x_k in frequent_x_ks:
                    del self.seen_x_ks[self.seen_x_ks.index(x_k)]
                    self.points_to_eliminate.append(x_k)

        if self.best_current_primal_value:
            print('Objective of the best feasible solution available: {}'.format(self.best_current_primal_value[-1]))
        print('Total number of points eliminated: ' + str(self.points_to_eliminate.__len__()))

    def oracle(self, lambda_k):
        print('# -- Oracle called. -- #')
        # return self._oracle_sequential(lam,lambda_k)
        return self._oracle_parallel(lambda_k)

    def _oracle_sequential(self, lambda_k):
        ###################################
        # -- Sequential Implementation -- #
        ###################################
        # same implementation as base class, only we add cutting planes
        start = time.time()

        scenario_objectives = []
        scenario_optimizers = []

        for scenario, scenario_label in enumerate(self.stoch_model.scenarios):
            scenario_program, x, y = self._construct_local_problem_and_add_cutting_planes(lambda_k, scenario, scenario_label)
            scenario_program.setParam('MIPGap', LOCAL_OPTIMIZATION_MIPGAP)
            scenario_program.optimize()

            # record local objective and optimizers
            scenario_objective = scenario_program.ObjVal
            temp = []
            for i in range(self.n_x):
                # TODO Only saving x here, what about y?
                temp.append(x[i].X)
            # CARE HERE, with respect to the base class implementation, I return a list of tuples, not a list of lists.
            scenario_optimizer = tuple(temp)
            scenario_objectives.append(scenario_objective)
            scenario_optimizers.append(scenario_optimizer)

        # update self.points_to_eliminate
        self._elimination_procedure(scenario_objectives, scenario_optimizers)

        print('Execution time of oracle: ' + str(time.time() - start))
        return self._calculate_oracle_outputs(scenario_objectives, scenario_optimizers)

    def _oracle_parallel(self, lambda_k):
        #################################
        # -- Parallel Implementation -- #
        #################################

        # This parallel implementation will first generate gurobi models, one for each scenario; then write then into
        # .mps files (this happens centrally). We then start multiprocesses to read the .mps files, solve the models,
        # and report back the solution. The reason why we have to go through this writing (found this after 2 days
        # trying to keep everything in main memory) is that Gurobi models cannot be pickled (which is required by pretty
        # much anything that attempts to do parallel computations).
        # See: https://groups.google.com/forum/#!topic/gurobi/fwLRrWLLJqo

        start = time.time()

        # TODO move this in init; try both std multiproc as well as pathos version.
        pool = mp.ProcessingPool(NUM_CORES)
        n_scen = self.stoch_model.scenarios.__len__()

        print('Solving {} individual models for the scenarios.'.format(self.stoch_model.scenarios.__len__()))

        scenario_lambda_k = []
        # we create |S| copies of the points to eliminate; generally, we may want to distribute a different set of
        # points to eliminate to each individual sub problem optimizer
        scenario_points_to_eliminate = []
        for scenario, scenario_label in enumerate(self.stoch_model.scenarios):
            scenario_lambda_k.append(lambda_k[scenario*self.n_x:(scenario+1)*self.n_x])
            scenario_points_to_eliminate.append(self.points_to_eliminate)

        # Distribute computations!
        results = pool.map(_external_optimization_of_single_scenario,
                           range(n_scen),
                           scenario_lambda_k,
                           scenario_points_to_eliminate
                           )

        scenario_optimizers = []
        scenario_objectives = []
        for result in results:
            scenario_optimizers.append(result[0])
            scenario_objectives.append(result[2])

        # Do elimination if activated
        if self.with_elimination:
            start_elim = time.time()
            print('Start elimination ...')
            self._elimination_procedure(scenario_objectives, scenario_optimizers, type='all')
            print('Execution time of elimination: ' + str(time.time() - start_elim))

        print('Execution time of oracle: ' + str(time.time() - start))
        x_k, d_k, subgrad_lambda_k = self._calculate_oracle_outputs(scenario_objectives, scenario_optimizers)
        if self.best_current_primal_value:
            print('Optimality gap % (certified): {}'.format(
                float(self.best_current_primal_value[-1] - d_k) / float(d_k) * 100))
        return x_k, d_k, subgrad_lambda_k

    def _oracle_distributed(self, lambda_k):
        pass
        # THIS SHOULD GO IN A TRY STATEMENT IN init
        # 1) if it doesnt exist, create view, then load Gurobi on the engines

            # create variable  LOCAL_OPTIMIZATION_MIPGAP in the engine



        # 2) call external function

    def _evaluate_inner_solution_performance(self, x_k):
        # construct evaluation problem for fixed x_k
        evaluation_program = gb.Model()

        # -- variables -- #
        y = {}
        for s, scenario_label in enumerate(self.stoch_model.scenarios):
            scenario = self.stoch_model.scenarios[scenario_label]
            for i in range(self.n_y):
                y[s,i] = evaluation_program.addVar(obj=scenario.c[2][i]*scenario.probability,
                                                   lb=scenario.lb[2][i],
                                                   ub=scenario.ub[2][i],
                                                   vtype=scenario.vtype[2][i])
        evaluation_program.update()

        # -- constraints -- #
        for s, scenario_label in enumerate(self.stoch_model.scenarios):
            scenario = self.stoch_model.scenarios[scenario_label]
            # constraints, second-stage; contribution to lhs from A[2,1]
            lhs_1 = defaultdict(gb.LinExpr)
            for entry in scenario.A[2,1].items():
                row, column = entry[0]
                value = entry[1]
                lhs_1[row] = lhs_1[row] + value*x_k[column]
            # contribution from A[2,2]
            lhs_2 = defaultdict(gb.LinExpr)
            for entry in scenario.A[2,2].items():
                row, column = entry[0]
                value = entry[1]
                lhs_2[row] = lhs_2[row] + value*y[s,column]

            for row, value in enumerate(scenario.b[2]):
                    evaluation_program.addConstr(lhs_1[row] + lhs_2[row],
                                                 scenario.b_sense[2][row],
                                                 scenario.b[2][row])
        evaluation_program.update()
        evaluation_program.setParam('OutputFlag', False)
        # evaluation_program.setParam('MIPGap', LOCAL_OPTIMIZATION_MIPGAP)
        evaluation_program.optimize()

        # contribution from x_k
        objective_contribution_from_x_k = np.dot(self.stoch_model.c[1],x_k[:self.n_x])

        return evaluation_program.ObjVal + objective_contribution_from_x_k

