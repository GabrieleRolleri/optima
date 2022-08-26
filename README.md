# optima
In the context of the Optima project (https://optima-hpc.eu) and with EnginSoft we are trying to prove that FPGA cards can offer an advantage to HPC applications in terms of computing speed and power effiency. In this particular instance a lattice Boltzman method fluid dynamics simulation has been implemented for CPU execution, afterwards about 60% of the computation has been moved onto FPGA cards, the resulting speed up still has to be measured.
Most of the coding work has been done on MaxIDE, while all the computation took place on jumax machines at Jülich Supercomputing Centre.
## main Branch
This branch contains the original CPU code for the simulation.
## parallel Branch
This branch contains the modified CPU code integrated with the Sli-C interface and all the code necessary for proper DFE functioning.
