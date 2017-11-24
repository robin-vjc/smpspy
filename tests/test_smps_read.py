from collections import OrderedDict
import numpy as np

from smpspy.oracles import TwoStage_SMPS_InnerProblem
from smpspy.smps_read import StochasticModel

BENCHMARK_PATH = '../smpspy/benchmark_problems/'


def test_stochastic_model_from_smps_read():
    smps_model_name = BENCHMARK_PATH+'2_sizes/sizes10'
    sm = StochasticModel(smps_model_name)

    assert sm.mode_of_modification == 'REPLACE'
    assert sm.model_name == 'SIZES'

    #####################################
    # Base Deterministic Model Matrices #
    #####################################

    assert type(sm.A) == dict
    assert len(sm.A) == 4, 'matrix A has four blocks for 2 stage stoch programs'
    # the four blocks are:
    # (1,1): 1st stage constraints on 1st stage decision
    # (1,2): 1st stage constraints on 2nd stage decision SHOULD BE ALWAYS EMPTY
    assert len(sm.A[1,2]) == 0, 'this sub-block of matrix A should ALWAYS be 0.'
    # (2,1), (2,2): 2nd stage constraints (include 1st and 2nd stage decisions)
    assert len(sm.A[2,2]) == 150, 'this spars matrix should contain 150 non-zero entries'

    #############
    # Scenarios #
    #############
    assert len(sm.scenarios) == 10, 'sizes10 model should have 10 scenarios'
    assert type(sm.scenarios) == OrderedDict
    assert sm.scenarios.keys()[0] == 'SCEN01', 'name of first scenario is SCEN01 (ordered dict!)'
    assert sm.scenarios['SCEN01'].probability == 0.1

    # Scenario matrices
    # for this model, all modifications are on RHS (b vector) so..
    assert (sm.scenarios['SCEN01'].A[1,1]!=sm.A[1,1]).nnz == 0, 'det. A, and A in first scenario should be equal.'
    for key in sm.scenarios['SCEN01'].data_modification:
        # print(key)  #shows what entries are being modified in this scenario
        assert key[1] == 'RHS'


def test_caroe_schultz_model():
    print('Testing Caroe Schultz Model with 10 Scenarios')
    smps_model_name = BENCHMARK_PATH+'2_caroe_schultz/caroe_schultz_10'
    sm = TwoStage_SMPS_InnerProblem(smps_model_name)

    # c_x
    assert (sm.stoch_model.c[1] == [-1.5, -4]).all()
    # c_y = (16, 19, 23, 28)/121 for 121 scenarios
    np.testing.assert_allclose([-16., -19., -23., -28.], np.array(sm.stoch_model.c[2]), atol = 0.1)

    # LHS for nominal model is psi_1 = psi_2 = 5
    np.testing.assert_allclose([5, 5], np.array(sm.stoch_model.b[2]), atol = 0.1)
    assert sm.stoch_model.A[2, 1][0,0] == 1.0
    assert sm.stoch_model.A[2, 1][1,1] == 1.0

    # First constraint
    assert sm.stoch_model.A[2, 2][0,0] == 2.0
    assert sm.stoch_model.A[2, 2][0,1] == 3.0
    assert sm.stoch_model.A[2, 2][0,2] == 4.0
    assert sm.stoch_model.A[2, 2][0,3] == 5.0

    # Second Constraint
    assert sm.stoch_model.A[2, 2][1,0] == 6.0
    assert sm.stoch_model.A[2, 2][1,1] == 1.0
    assert sm.stoch_model.A[2, 2][1,2] == 3.0
    assert sm.stoch_model.A[2, 2][1,3] == 2.0

    # Verify stoch model / det equivalent is OK
    assert sm.stoch_model.scenarios['SCEN1'].probability == 0.1

    det_eq = sm.stoch_model.generate_deterministic_equivalent()
    det_eq.optimize()
    assert abs(det_eq.ObjVal+54.40) <= 0.1


def test_caroe_schultz_paper_model():
    """ test the eact model we have in the paper """
    print('Testing Caroe Schultz Model, exact as in the paper, with 121 Scenarios')
    smps_model_name = BENCHMARK_PATH+'2_caroe_schultz/caroe_schultz_121'
    sm = TwoStage_SMPS_InnerProblem(smps_model_name)

    det_eq = sm.stoch_model.generate_deterministic_equivalent()
    det_eq.optimize()

    assert abs(det_eq.Objval+62.29) <= 0.1

def test_multistage_smps_read():
    print('#Test: smps_read should work on models with more than 2 stages too.')

    sm = StochasticModel(BENCHMARK_PATH+'posts/2_sg/sgpf5y3')  # this is 3 stages
    assert len(sm.scenarios) == 25
    assert len(sm.periods) == 3

    sm = StochasticModel(BENCHMARK_PATH+'posts/2_sg/sgpf5y6')  # this is 6 stages
    assert len(sm.periods) == 6


def test_deterministic_equivalent():
    smps_model_name = BENCHMARK_PATH+'2_sizes/sizes10'
    sm = StochasticModel(smps_model_name)

    # Generate deterministic equivalent, and solve it with Gurobi
    det_eq = sm.generate_deterministic_equivalent()
    det_eq.params.MIPGap = 0.01
    det_eq.optimize()

    # Try the plotting function
    # sm.plot_scenario_tree()



# Confirmed working datasets
# - SIZES we only have .tim for the version with 10 scenarios (largest). Solves pretty instant, optimal
# objective confirmed to be correct.
# - SEMI model doesn't work; seems at least .cor file is corrupt
# sm = StochasticModel('./data/2_semi/2_semi3')
# - SMKP model seems to be OK, all instances 1...30 solve in <1s to 1% MipGap
# sm = StochasticModel('./data/2_smkp/smkp_3')
# - SSLP
# sm = StochasticModel('./data/2_sslp/sslp_5_25_50')  # works, objective confirmed, solves instant
# sm = StochasticModel('./data/2_sslp/sslp_10_50_500')  # works, objective confirmed, 2 min to build, 1200sec to solve
# - Vaccine
# SOME INFEASIBILITY, DEBUG WITH IIS
# sm = StochasticModel('./data/2_vaccine/vac100b')
# - DCAP
# sm = StochasticModel('./data/2_dcap/2_dcap233_200')  # works, objective confirmed, solves instant
# sm = StochasticModel('./data/2_dcap/2_dcap233_300')  # works, objective confirmed, solves instant
# sm = StochasticModel('./data/2_dcap/2_dcap233_500')  # works, objective confirmed, solves instant
# sm = StochasticModel('./data/2_dcap/dcap342_500')  # works, objective confirmed, solves instant
