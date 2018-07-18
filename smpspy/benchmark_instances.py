####################################################################################################
# 2-Stage Stochastic Instances
####################################################################################################

import os

_main_smps_path = os.path.dirname(__file__)+'/benchmark_problems/'

# OK, dual_step() time 2 sec
_dcap_path  = _main_smps_path + '2_dcap/'
_dcap_names = [
    '2_dcap233_200',
    '2_dcap233_300',
    '2_dcap233_500',
    'dcap243_200',
    'dcap243_300',
    'dcap243_500',
]
_dcap_instances = list(map(lambda name: _dcap_path + name, _dcap_names))

# OK ~160 sec per it (time limit in action, optimization of certain scenarios is heavy)
_semi_path = _main_smps_path + '2_semi/'
_semi_names = [
    '2_semi2',  # 2 scenarios, 153 sec per it. 1st scenario ~1.0 sec, 2nd scenario 150s (time limit);
    '2_semi3',  # 3 scenarios, 167 sec per it. 1st scenario: 1 sec, 2nd: 150sec, 3rd: 0.3 sec
    '2_semi4',  # 4 scenarios, 159 sec per it. 1st insta, 2nd time limit, 3rd insta, 4th time limit
]
_semi_instances = list(map(lambda name: _semi_path + name, _semi_names))

# OK, dual_step() time 0.3 sec
_sizes_path = _main_smps_path + '2_sizes/'
_sizes_names = [
    'sizes10'
]
_sizes_instances = list(map(lambda name: _sizes_path + name, _sizes_names))

# OK, 18 sec per it
_smkp_path = _main_smps_path + '2_smkp/'
_smkp_names = [
    'smkp_1',
    'smkp_2',
    'smkp_3',
    'smkp_4',
    'smkp_5',
    'smkp_6',
    'smkp_7',
    'smkp_8',
    'smkp_9',
    'smkp_10',
    'smkp_11',
    'smkp_12',
    'smkp_13',
    'smkp_14',
    'smkp_15',
    'smkp_16',
    'smkp_17',
    'smkp_18',
    'smkp_19',
    'smkp_20',
    'smkp_21',
    'smkp_22',
    'smkp_23',
    'smkp_24',
    'smkp_25',
    'smkp_26',
    'smkp_27',
    'smkp_28',
    'smkp_29',
    'smkp_30',
]
_smkp_instances = list(map(lambda name: _smkp_path + name, _smkp_names))

# OK SSLP, 0.1 sec per it
_sslp_path = _main_smps_path + '2_sslp/'
_sslp_names = [
    'sslp_5_25_3_mymod',
    'sslp_5_25_50',
    'sslp_5_25_100',
    'sslp_10_50_10_mymod',
    'sslp_10_50_50',
    'sslp_10_50_100',
    'sslp_10_50_500',
    'sslp_10_50_1000',
    'sslp_10_50_2000',
    'sslp_15_45_5',
    'sslp_15_45_10',
    'sslp_15_45_15',
]
_sslp_instances = list(map(lambda name: _sslp_path + name, _sslp_names))

STOCH_INSTANCES = {'dcap': _dcap_instances,
                   'semi': _semi_instances,
                   'sizes': _sizes_instances,
                   'smkp': _smkp_instances,
                   'sslp': _sslp_instances,
                   # 'all': _dcap_instances + _semi_instances + _sizes_instances + _smkp_instances + _sslp_instances
                   }

# # NOT OK VACCINE
# # WHY: This is chance constrained
# path = './applications/smpspy/benchmark_problems/2_vaccine/'
# instances = [
#     'vac100a'
# ]

# # NOT OK POSTS - 2SG
# # WHY: pure continuous; seems more than 2 stages
# path = './applications/smpspy/benchmark_problems/posts/2_sg/'
# instances = [
#     'sgpf5y3',
#     'sgpf5y4',
#     'sgpf5y5',
#     'sgpf5y6',
# ]