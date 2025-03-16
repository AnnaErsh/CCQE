# CCQE
Standalone version of the quasi-elastic reaction modelling according to the Llewellyn-Smith formula.
Currently is dependent on the vect.h file that is a part of the NuWro project.

## License
This project is licensed under the **GNU General Public License v3.0**.
Portions of this code are derived from the NuWro project, which is also licensed under GPLv3.
See the [LICENSE](LICENSE) file for full details.

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
