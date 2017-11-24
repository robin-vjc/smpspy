# Create SMPS files for the example stochastic model (16) in Caroe / Schultz.
# Procedure:
# 1. Run script below with `build_only_nominal_model = True`
# 2. Rename `caroe_schultz.mps` into `caroe_schultz.cor`
# 3. Write .tim and .sto files (by hand for now)

import gurobipy as gb
import numpy as np
import os


def generate_caroe_schultz_cor_file(output_file='caroe_schultz', build_only_nominal_model=True, n_S=121):
    """ Generate nominal model .mps. """
    m = gb.Model()
    # m.setObjective(0, sense=gb.GRB.MAXIMIZE)
    m.setObjective(0, sense=gb.GRB.MINIMIZE)

    p_S = float(1)/n_S

    # Variables
    x = {
        1: m.addVar(vtype=gb.GRB.INTEGER, lb=0, ub=5, obj=-1.5, name='x1'),
         2: m.addVar(vtype=gb.GRB.INTEGER, lb=0, ub=5, obj=-4, name='x2')
    }
    m.update()

    # Constraints
    cst = {}
    y = {}


    if build_only_nominal_model:
        cst[0] = m.addConstr(x[1] + x[2] >= 0, name='cstr0')

        y[1] = m.addVar(vtype=gb.GRB.BINARY, obj=-16, name='y1')
        y[2] = m.addVar(vtype=gb.GRB.BINARY, obj=-19, name='y2')
        y[3] = m.addVar(vtype=gb.GRB.BINARY, obj=-23, name='y3')
        y[4] = m.addVar(vtype=gb.GRB.BINARY, obj=-28, name='y4')

        m.update()
        nominal_psi_1 = 5
        nominal_psi_2 = 5

        # cst[0] = m.addConstr(x[1] + x[2] >= 0, name='cstr0')
        cst[1] = m.addConstr(2 * y[1] + 3 * y[2] + 4 * y[3] + 5 * y[4] <= nominal_psi_1 - x[1], name='cstr1')
        cst[2] = m.addConstr(6 * y[1] + 1 * y[2] + 3 * y[3] + 2 * y[4] <= nominal_psi_2 - x[2], name='cstr2')

    else:
        # for psi_1 in range(5,5+int(np.sqrt(n_S))):
        #     for psi_2 in range(5,5+int(np.sqrt(n_S))):
        for psi_1 in np.linspace(5,15,int(np.sqrt(n_S))):
            for psi_2 in np.linspace(5,15,int(np.sqrt(n_S))):
                y[psi_1,psi_2,1] = m.addVar(vtype=gb.GRB.BINARY, obj=-16 * p_S)
                y[psi_1,psi_2,2] = m.addVar(vtype=gb.GRB.BINARY, obj=-19 * p_S)
                y[psi_1,psi_2,3] = m.addVar(vtype=gb.GRB.BINARY, obj=-23 * p_S)
                y[psi_1,psi_2,4] = m.addVar(vtype=gb.GRB.BINARY, obj=-28 * p_S)
                m.update()

                cst[psi_1, psi_2, 1] = m.addConstr(2 * y[psi_1,psi_2,1] + 3 * y[psi_1,psi_2,2] + 4 * y[psi_1,psi_2,3] + 5 * y[psi_1,psi_2,4] <= psi_1 - x[1])
                cst[psi_1, psi_2, 2] = m.addConstr(6 * y[psi_1,psi_2,1] + 1 * y[psi_1,psi_2,2] + 3 * y[psi_1,psi_2,3] + 2 * y[psi_1,psi_2,4] <= psi_2 - x[2])

    m.update()
    # m.optimize()
    output_filename = output_file+'_'+str(n_S)
    m.write(output_filename+'.mps')
    os.rename(output_filename+'.mps', output_filename+'.cor')


def generate_caroe_schultz_sto_and_tim_files(output_file='caroe_schultz', n_S=121):
    p_S = float(1) / float(n_S)

    with open(output_file+"_"+str(n_S)+".sto", "a") as sto_file:
        preamble = "STOCH\nSCENARIOS\t DISCRETE\n"
        sto_file.write(preamble)

        scenario_index = 0
        for psi_1 in np.linspace(5,15,int(np.sqrt(n_S))):
            for psi_2 in np.linspace(5,15,int(np.sqrt(n_S))):
                scenario_index += 1
                line_1 = "SC\tSCEN{}\tROOT\t{}\tPERIOD2\n".format(str(scenario_index), str(p_S))
                line_2 = "\tRHS\tcstr1\t{}\n".format(psi_1)
                line_3 = "\tRHS\tcstr2\t{}\n".format(psi_2)
                sto_file.write(line_1)
                sto_file.write(line_2)
                sto_file.write(line_3)

    with open(output_file + "_" + str(n_S) + ".tim", "a") as tim_file:
        tim_file.write("TIME\n")
        tim_file.write("PERIODS\tLP\n")
        tim_file.write("\tx1\tcstr0\tPERIOD1\n")
        tim_file.write("\ty1\tcstr1\tPERIOD2\n")
        tim_file.write("ENDATA\n")
