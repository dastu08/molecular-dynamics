# Molecular Dynamics Simulation

> This code uses the
> [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library in
> version 3.4.0 for handling arrays and matrixes.

## Compiling the source code
```bash
# create directory structure
mkdir data
mkdir data/02
mkdir build
# build the project
cd build
cmake .. -G Ninja
ninja
# exectute the program
./simulation
```