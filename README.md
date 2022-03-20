# Molecular Dynamics / Monte Carlo Simulation

This code was written during the course _Simulation Methods in Statistical
Physics_ which is simulating Argon atoms interacting with a Lennard-Jones
potential.

In the `DEMO` part of the main function it runs both Molecular Dynamics (MD) and
Monte Carlo (MD) to see an example usage. 

The other code involved specific tasks related to written reports 1 to 5.


## Compiling the source code
> This code uses the
> [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library in
> version 3.4.0 for handling arrays and matrixes. Please ensure that is in your
> include path.

```bash
# create directory structure
mkdir data
# for every report create the directory like
mkdir data/01
# create a directory for the build files
mkdir build
# build the project
cd build
cmake .. -G Ninja
ninja
# exectute the program
./simulation
```

**This repository does not include a license which is intended.**