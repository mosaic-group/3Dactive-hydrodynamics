# 3D active-hydrodynamics
Repository for OpenFPM code that solves the 3D active Ericksen-Leslie Model

The codes in this repository requires the OpenFPM library.

## Citing This Work
If you use this software or its components in your research, please cite the relevant work as described in our [CITATION.md](./CITATION.md) file. This file provides detailed citation information for the software and associated research papers.

The file anaEulerian.cpp is the C++ code that can be used to simulate the active Ericksen-Leslie hydrodynamics in an Eulerian frame of reference.

The file anaHyb.cpp is the C++ code that can be used to simulate the active Ericksen-Leslie hydrodynamics in a Lagrangian frame of reference with remeshing at every time step.

The file annulus.cpp is the C++ code that can be used to simulate the active Ericksen-Leslie hydrodynamics in an Eulerian frame of reference for an annular geometry.

The file complexGeometries.cpp is the C++ code that can be used to simulate the active Ericksen-Leslie hydrodynamics in an Eulerian frame of reference any complex 3D geometry given point cloud data (.csv-file). The particle positions are given in the .csv-files (argv[2]) for the symmetric and the asymmetric geometry (the files have 7 columns specifying: x-position, y-position, z-position, bulk/boundary (0/1), x-component normal, y-component normal, z-component normal). The parametersComplex.txt (argv[1]) specifies the model and simulation parameters.   
