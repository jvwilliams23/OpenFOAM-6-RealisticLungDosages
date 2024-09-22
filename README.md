# OpenFOAM-6 with modifications for lung simulations with realistic dosages

This repository contains a modified version of the OpenFOAM `lagrangian/intermediate` library for running simulations with realistic dosages.
The main modification to the library is found in commit [`f31fe2f`](https://github.com/jvwilliams23/OpenFOAM-6-RealisticLungDosages/commit/f31fe2ffb74355abb4125b10ada437d7a66cf470). It is a simple if statement, which checks if a particle is active (floating) or inactive (deposited), to avoid deposited particles interfering with the surrounding flow and creating numerical issues.

We also provide two solvers in `applications/solvers`:
- `MPPICFoamRealDose` is a standard MPPIC solver that is linked to the modified library. This solver facilitates two-way and four-way coupling.
- `uncoupledMPPICFoamRealDose` modifies the standard MPPIC solver by explicitly setting $\alpha_p = 0$ and also setting the source term to 0 (enforcing one-way coupling).  
Both new solvers come several turbulence models linked to them, which may prove useful.

## Usage
To install the library and solver, please do the following
```bash
# clone repo
cd $WM_PROJECT_USER_DIR
git clone https://github.com/jvwilliams23/OpenFOAM-6-RealisticLungDosages.git DPMRealDosage
cd DPMRealDosage
# compile library
cd src/lagrangian/intermediate/
wclean
wmake
cd -
# compile solver
cd applications/solvers/MPPICFoamRealDose/
./Allwmake
```

