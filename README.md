# Parallelized Optimization of Stochastic Multistage Integer Programs

<p align="center">
  <img src="./notebooks/img/stoch_tree.png" alt="Scenarios Tree" width="60%" href="#"/>
</p>

### Installation

* install `nsopy`, see [instructions here](https://github.com/robin-vjc/nsopy).
* then (last step optional):
~~~~
>> git clone https://github.com/robin-vjc/smpspy
>> cd smpspy
>> pip install -e .
>> py.test
~~~~
* the package `gurobipy` is also required, instructions [here](http://www.gurobi.com/documentation/6.5/quickstart_mac/the_gurobi_python_interfac.html).

Now you can cd into the `notebooks` directory and run `jupyter notebook`; you should be able to 
execute all the notebooks without seeing any error.

### Introduction

Implementation of a dual decomposition system, applied to parallelize the solution 
of hard multistage (integer) stochastic programs.

* For an introduction to dual decomposition for multistage stochastic integer programming, 
please check [1].
* The contents of this repo can be used to reproduce the results in [2]. 

### Contents
* utilities to parse files in the SMPS format; see [basic use examples](./notebooks/Dual Methods Performance-0 Basic Usage.ipynb).
* solving duals of multistage programs, see [notebook 1](./notebooks/Dual Methods Performance-1 Parameter Selection.ipynb) and 
[notebook 2](./notebooks/Dual Methods Performance-2 Comparison.ipynb). These reproduce the contents of [2].

### Replicating results from the command line

#### Step 1: Generating the csv results files

To run individual experiments, create a csv file, for example named `text_experiments.csv` whose contents are:

~~~~
subtype,method_name,selected_parameter,instance_n
sslp,SG 1/k,0.1,0
sslp,SG const,10,3
~~~~
Then issuing the command
~~~~
>> python run_experiments.py -i text_experiments.csv -o text_experiments_output.csv
~~~~
will run the following experiments
 * the `instance=0` of class `sslp` with subgradient method 1/k, and initial stepsize `0.1`
 * instance `instance=3` of class `sslp` with constant subgradient method, stepsize `10`
 
The outputs of the experiments will be saved in the file `text_experiments_output.csv`. Settings available:
* Available problem classes: `[dcap, semi, sizes, smkp, sslp]`
* Available methods: `[SG 1/k, SG const, UPGM, UDGM, UFGM, DSA, TA 1, TA 2, CP, bundle]` 

**Note:** running
~~~~
>> python run_experiments.py -i parameter_selection.csv -o experiments_output.csv 
~~~~
will run the complete battery of experiments discussed in the paper. The results of these experiments is in `final_experiments_record.csv`.

#### Step 1: Generating the plots

Given the file `final_experiments_record.csv`, the plots in the paper can be reproduced by running
~~~~
jupyter notebook --no-browser --port=8889
~~~~
and either opening a web browser locally at `http://localhost:8889/`, or to connect remotely:
~~~~
ssh -N -f -L localhost:8888:localhost:8889 username@your_remote_host_name
~~~~
and then opening the browser at `http://localhost:8888/`. Running all the cells in the notebook `Dual Methods Performance-2 Comparison.ipynb` generates the desired plots. 

### References

[1] C. Caroe and R. Schultz, *Dual decomposition in stochastic integer programming*, Operations
Research Letters 24 (1999), no. 1, 37-45. [Publication link](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.15.5253).

[2] Robin Vujanic, Peyman Mohajerin Esfahani, *Title*, Journal, 2016.

### Contributing

Pull requests are very welcome. Check the [TODO](TODO.txt). 


## Cite

~~~~
@article{Vujanic2018,
	title={Dual Decomposition of Stochastic Integer Programs: New Results and Experimental Comparison of Solution Methods},
	author={Vujanic, Robin and Esfahani, Peyman Mohajerin},
	journal={TO BE COMPLETED},
	volume={TO BE COMPLETED},
	number={TO BE COMPLETED},
	pages={TO BE COMPLETED},
	year={2018},
	publisher={TO BE COMPLETED}
}
~~~~ 