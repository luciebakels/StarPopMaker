# StarPopMaker
Creating a population of stars and following it as it evolves. Includes main-sequence stars, high-mass X-ray binaries, black holes, neutron stars, and white dwarfs.Written in cpp using MPI.

Compiling the code:
- Necessary software: c++11 and mpi
- "make" in "source" directory

Running the code:
- Example parameter file in "examples" directory.
- mpirun -n "Number_of_cores_to_run_on" mpipop.exe example.param
- It runs embarrassingly parallel, writing out as many files as there are cores and combining those at the end of the run.

Input:
- TotalTime: sets the total time the simulation runs in Myr
- TimeStep: sets how large each time-step is in Myr
- TimeRes: inactive in this version

- Survival Fraction: "fudge-factor" for HMXB formation (see Power et al. 2009).
- IMF: choice of initial mass function, either Kroupa (K) or Salpeter (S).
- MinimumMass: sets the minimum stellar mass.
- MaximumMass: sets the maximum stellar mass.

- SFRType: choice of star formation history, either Exponential (E), Constant (C), Starburst (S), Disk (D), or Logarithmic-Normal (L)
Note: if Disk (D) is chosen, the parameters: DiskRadius and DiskBins need to be set (needs updating)

- TotalMass: sets the total stellar mass that will be generated
- BurstWidth and PeakTime set the width and the peak of the Logarithmic-Normal star formation history (see Gladders et al. 2013)
- SFR can be set in case of a constant (C) star formation rate (in Msun/Myr).
- StromgrenOn can be chosen in case of a Distk (D)

- Profile, CNFW, DeltaVir, GasCoeff, Radius, GasFraction, GasType, RadiationOn, PhotonsOn, and LogBins are inactive in this version

Output:
- NumData.txt: lists the number of each object at each time step
- LumData.txt: lists the luminosity and ionising luminosity at each time step
- MassData.txt: lists the mass of each object at each time step
- HMXB.txt: lists each high-mass X-ray binary that has formed and its properties.
- BH.txt: lists each black hole that has formed and its properties.
- NS.txt: lists each neutron star that has formed and its properties.
