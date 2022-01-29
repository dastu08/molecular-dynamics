# Molecular Dynamics Simulation

## Compiling the source code
```bash
# create directory structure
mkdir data
mkdir data/01
mkdir build
# build the project
cd build
cmake .. -G Ninja
ninja
# exectute the program
./simulation
```