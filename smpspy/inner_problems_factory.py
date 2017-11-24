from smpspy.smps.benchmark_instances import STOCH_INSTANCES
from smpspy.smps.oracles import TwoStage_SMPS_InnerProblem


def BenchmarkInnerProblemsFactory(type, instance_n, subtype='all'):
    if type == '2 stage stoch':
        # 52 instances in 'all'. Subtypes available: dcap, semi, sizes, smkp, sslp
        ip = TwoStage_SMPS_InnerProblem(STOCH_INSTANCES[subtype][instance_n])
        ip.instance_type = 'smpspy'
        ip.instance_subtype = subtype
        ip.instance_name = STOCH_INSTANCES[subtype][instance_n].split('/')[-1]
        return ip