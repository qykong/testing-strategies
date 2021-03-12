## Epidemic simulations and testing

This repository hosts source code for the paper submission "Contact Tracing: Computational Bounds, Limitations and Implications". The code can be used to simulate compartmental epidemic models over random networks. The simulation code was originally from this [work](https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks) with modifications for including tracing and testing of epidemics.


## Quick tour
We implement some parts of the simulation code in [Cython](https://cython.org/) in order to leverage its high performance. For this reason, before running the simualtions, you would need to compile the Cython scripts first.
```bash
python pyscripts/setup.py build_ext --inplace
```

After this, we will be able to run the `main.py`. Here is an example
```bash
python3 main.py -pf parameters/abstract/SIR_with_Ho.json -num_tracers 100 -nI 10 -oracle_type forward -rate_I2R 0.0475 -rate_I2Ho 0.0025 -pCT 1 -seed 190 -simIts 50 -nxIts 10 -pQ 1 -pRT 1 -N 10000 -parallel 1  -return_full_data True -save_limited_data False -nt er -np '{"p_without_N": 2.708333}' -rate_IS2II 0.6
```
where `p_without_N=R0/(beta/(beta+gamma))`. This command runs SIR simulations with standard forward contact tracing and given parameters. Here's a quick look at the output of this simulations
```
This file will be used: parameters/abstract/SIR_with_Ho.json
An abstract model file is provided.
The following parameters will be applied:
{'nI': 10, 'maxdays': 1000, 'N': 10000, 'pC2T': 1.0, 'pQ': 1.0, 'pCT': 1.0, 'pRT': 1.0, 'tlist_size': 0, 'type': 'abstract', 'p2bCTs': {'S': 1.0, 'I': 1.0, 'R': 1.0, 'Ho': 0}, 'p2bRTs': {'S': 1.0, 'I': 1.0, 'R': 1.0, 'Ho': 0}, 'positive_instantaneous': ['Ho'], 'positive_induced': ['I'], 'to_trace_instantaneous': ['I'], 'recovered': ['R', 'Ho'], 'infected': ['I'], 'return_statuses': ('S', 'I', 'R', 'Ho'), 'gH': <networkx.classes.digraph.DiGraph object at 0x7fb12860f2e0>, 'J': <networkx.classes.digraph.DiGraph object at 0x7fb12860f190>, 'required_rates': ['rate_I2R', 'rate_IS2II', 'rate_I2Ho'], '_comment1': 'From paper: Estimation of the Basic Reproduction Number of Novel Influenza A (H1N1) pdm09 in Elementary Schools Using the SIR Model', '_comment2': 'reported R0: 1.33', 'seed': 190, 'network_type': 'er', 'network_params': '{"p_without_N": 2.708333}', 'oracle_type': 'forward', 'num_tracers': 100, 'nxIts': 10, 'simIts': 50, 'return_full_data': True, 'parallel': True, 'num_nodes': 0, 'reintroduction_nums': 0, 'rate_I2R': 0.0475, 'rate_IS2II': 0.6, 'rate_I2Ho': 0.0025}
Using random seed: 190
Start simulating loops...

Average Ratio of total infected 0.9522613211219033
Average Ratio of max infected 0.7308121562051657
Average Days to Epidemic End 24.535519080909193
Average Ratio of quarantined 0.23671631519437236
No secondary infection counted
Average Ratio of Positive Contact Tracing 0.20283854098340812
Average Ratio of Positive Random Tracing 0.007462874786371902
```

A complete list of arguments that are avaialbe from `main.py` can be found with `python3 main.py -pf parameters/abstract/SIR_with_Ho.json -h`.
```
usage: main.py [-h] [-seed SEED] [-nt NETWORK_TYPE] [-np NETWORK_PARAMS] [-oracle_type ORACLE_TYPE] [-num_tracers NUM_TRACERS] [-nxIts NXITS] [-simIts SIMITS] [-of OF] [-return_full_data RETURN_FULL_DATA] [-parallel PARALLEL]
               [-num_nodes NUM_NODES] [-pC2T PC2T] [-pQ PQ] [-pCT PCT] [-pRT PRT] [-tlist_size TLIST_SIZE] [-nI NI] [-maxdays MAXDAYS] [-N N] -rate_I2R RATE_I2R -rate_IS2II RATE_IS2II -rate_I2Ho RATE_I2HO

optional arguments:
  -h, --help            show this help message and exit
  -seed SEED            seed (if not specified seed is random
  -nt NETWORK_TYPE, --network-type NETWORK_TYPE
                        specify the network simulator to use
  -np NETWORK_PARAMS, --network-params NETWORK_PARAMS
                        specify the network parameter to use
  -oracle_type ORACLE_TYPE
                        type of oracle to use
  -num_tracers NUM_TRACERS
                        number of tracers
  -nxIts NXITS          number of iterations of random networks
  -simIts SIMITS        number of simulations per random network
  -of OF                path to the file to save results
  -return_full_data RETURN_FULL_DATA
                        Determine if individual transimissions are tracked (This will introduce additional computational overhead)
  -parallel PARALLEL    Run simulation in parallel.
  -num_nodes NUM_NODES  number of first set of nodes selected for computing secondary infections
```