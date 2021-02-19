## Epidemic simulations and testing

This repository hosts source code for the paper submission "Cautionary Tales of Contact Tracing Across Multiple Epidemics and the Promise of Intelligent Testing Strategies for the Future". The code can be used to simulate compartmental epidemic models over random networks. The simulation code was originally from this [work](https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks) with modifications for including tracing and testing of epidemics.


## Quick tour
We implement some parts of the simulation code in [Cython](https://cython.org/) in order to leverage its high performance. For this reason, before running the simualtions, you would need to compile the Cython scripts first.
```bash
python pyscripts/setup.py build_ext --inplace
```

After this, we will be able to run the `main.py`. Here is an example
```bash
python3 main.py -pf parameters/abstract/SIR_influenza_A_with_Ho.json -num_tracers 1 -nI 10 -oracle_type none -rate_I2R 0.0475 -rate_I2Ho 0.0025 -pCT 1 -seed 190 -simIts 50 -nxIts 10 -pQ 1 -pRT 1 -N 10000 -parallel 1  -return_full_data True -save_limited_data False -nt poisson -np '{"p": 0.000377746107944128}' -rate_IS2II 0.6
```
This command runs SIR simulations with standard forward contact tracing and given parameters. A complete list of arguments that are avaialbe from `main.py` can found in the following table: TODO
