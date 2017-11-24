# Run `python run_final_experiments.py -h` for help on usage.

import argparse
import time

import numpy as np
import pandas as pd
from nsopy.methods.method_loggers import SlimDualMethodLogger
from nsopy.methods.methods_factory import DualMethodsFactory
from nsopy.methods.utils import record_logger
from smpspy.inner_problems_factory import BenchmarkInnerProblemsFactory, STOCH_INSTANCES

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
Runs the desired battery of experiments. File `input.csv` can be specified in two ways:
    (a) Columns: subtype, method_name, selected_parameter
        Script will run all instances in 'subtype', using method 'method_name' with
        the parameter tuning 'selected_parameter'.

        Example content of `input.csv`:
            subtype,method_name,selected_parameter
            sslp,bundle,0.01
            sslp,CP,0.1
            sslp,DSA,0.1

    (b) Columns: subtype, method_name, selected_parameter, instance_n
        Script will only run the `instance_n` of `subtype` with specified method and parameter.

        Example content of `input.csv`:
            subtype,method_name,selected_parameter,instance_n
            sslp,CP,0.1,0
            semi,DSA,1,2
            sizes,DSA,10,0

        Row `sizes,DSA,10,0`, e.g., indicates to solve instance number 0 of class sizes, with DSA
        method with and parameter setting 10.
""")
    parser.add_argument('-i', metavar='input.csv', type=str, help='input.csv file', required=True)
    parser.add_argument('-o', metavar='output.csv', type=str, default='output.csv',
                        help='file in which results are appended (default: output.csv)')

    args = parser.parse_args()

    methods = []
    loggers = []

    # experiments_to_run = pd.read_csv(sys.argv[1])
    experiments_to_run = pd.read_csv(args.i)

    assert 'subtype' in experiments_to_run.columns, 'csv file not valid; need subtype column'
    assert 'method_name' in experiments_to_run.columns, 'csv file not valid; need method_name column'
    assert 'selected_parameter' in experiments_to_run.columns, 'csv file not valid; need selected_parameter column'
    assert 4 >= len(experiments_to_run.columns) >= 3, 'csv file not valid; needs either 3 or 4 columns'

    # Version in which only stated instances are run
    ####################################################
    if len(experiments_to_run.columns) == 4:
        for index, row in experiments_to_run.iterrows():
            ip = BenchmarkInnerProblemsFactory(type='2 stage stoch', subtype=row.subtype, instance_n=row.instance_n)
            methods.append(DualMethodsFactory(ip, method=row.method_name, param=row.selected_parameter))
            loggers.append(SlimDualMethodLogger(methods[-1]))

    # Version where you pass methods parameters, and all instances are run
    ####################################################
    elif len(experiments_to_run.columns) == 3:
        for subtype in STOCH_INSTANCES:
            for instance_n, instance_filepath in enumerate(STOCH_INSTANCES[subtype]):
                ip = BenchmarkInnerProblemsFactory(type='2 stage stoch', subtype=subtype, instance_n=instance_n)
                for index, row in experiments_to_run[experiments_to_run.subtype==subtype].iterrows():
                    methods.append(DualMethodsFactory(ip, method=row.method_name, param=row.selected_parameter))
                    loggers.append(SlimDualMethodLogger(methods[-1]))

    # Setting appropriate dual domain for bundle and cutting planes methods
    for method in methods:
        if 'Cutting' in method.desc or 'Bundle' in method.desc:
            method.set_dual_domain(type='2 stage smpspy')

    for method_index, method in enumerate(methods):
            print(method)

            # Determine N_ORACLE_CALLS allocated for experiment
            start_time = time.time()
            method.oracle(method.projection_function(np.zeros(method.dimension, dtype=float)))
            oracle_time = time.time() - start_time

            # Set value
            N_ORACLE_CALLS = max(min(1000,  30*60/oracle_time), 100)
            print('oracle_time  = {}, N_ORACLE_CALLS = {}'.format(oracle_time, N_ORACLE_CALLS))

            i = 0
            while method.oracle_calls < N_ORACLE_CALLS and i < N_ORACLE_CALLS:
                method.dual_step()
                i += 1  # we have to add this because cutting planes stop querying oracles when they find the solution

            record_logger(loggers[method_index], filename=args.o)
