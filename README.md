# CCQE
Standalone version of the quasi-elastic reaction modelling according to the Llewellyn-Smith formula. 
Currently is dependen on the vect.h file that is a part of the NuWro project.

# Installation
Create a build directory:
```
mkdir build
```
Enter this directory:
```
cd build
```
Execute cmake:
```
cmake ../path_to_CMakeLists.txt
```
Execute the created Makefile:
```
make
```

# Run
The code needs 4 bool arguments:
**ifCC** if the channel is charged current or neutral current
**ifproton** if the outgoing hadron is a proton or neutron
**ifmuon** if the outgoing lepton is a muon
**ifnu** neutrino or antineutrino reaction

To run the code, use:
```
./CCQE 1 1 1 1   
```
Where 1 1 1 1 should be exchanged to the desired reaction configuration
