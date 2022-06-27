# ELECTRA
Open source code from the "GENESIS" project

The ELECTRA (Electron Transport) code solves the Boltzmann transport equation (BTE) for charge transport using the wave-vector and energy dependent momentum relaxation time approximation. The simulator computes the charge transport coefficients. It considers charge carrier scattering with phonons, ionized dopants, and alloy scattering. ELECTRA takes as input the electronic bandstructure and certain scattering parameters as detailed below, and returns the charge transport coefficients.

ELECTRA is written in MATLAB® and makes use of the Parallel Computing Toolbox™. A version ported in C language exists, but it is still under test.

The ELECTRA_v1 folder contains all the source coded needed to run ELECTRA plus the input file and some input data as examples. 

contact: patrizio.graziosi@cnr.it
