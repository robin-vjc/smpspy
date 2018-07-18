from smpspy.inner_problems_factory import BenchmarkInnerProblemsFactory
from smpspy.smps_read import StochasticModel


def test_instances_are_importable():
    """
    first check that instances are actually importable, and make verifications that basic data about
    the stochastic model is correcrt.
    """
    # print(STOCH_INSTANCES['sslp'][0])
    ip = BenchmarkInnerProblemsFactory(type='2 stage stoch', subtype='sslp', instance_n=0)

    assert ip.n_scenarios == 3
    assert ip.n_stages == 2
    assert ip.n_x == 5, 'number of first stage variables should be 5 for sslp_5_25_3'
    assert ip.n_y == 130, 'number of second stage variables should be 130'
    assert ip.dimension == ip.n_x*ip.n_scenarios, 'number of coupling constraints not right'

    assert type(ip.stoch_model) == StochasticModel
