from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='smpspy',
    version='1.0',
    description='Oracle for Stochastic Mixed Integer Programs',
    # url='http://github.com/storborg/funniest',
    author='Robin Vujanic',
    author_email='vjc.robin@gmail.com',
    install_requires = required,
    packages=find_packages(),
    # packages=['dmp'],
)

# Version 2.0: when we have implemented a working version of the transport planners