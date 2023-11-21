# Citation for 3Dactive-hydrodynamics Repository

This file provides citation information for the 3Dactive-hydrodynamics repository hosted on GitHub.

## Repository Information
- **Repository:** [3Dactive-hydrodynamics](https://github.com/mosaic-group/3Dactive-hydrodynamics)
- **Description:** Repository for OpenFPM code that solves the 3D active Ericksen-Leslie Model.
- **License:** GPL-3.0.
- **Main Features:**
  - `anaEulerian.cpp`: C++ code for simulating the active Ericksen-Leslie hydrodynamics in an Eulerian frame of reference.
  - `anaHyb.cpp`: C++ code for simulating the active Ericksen-Leslie hydrodynamics in a Lagrangian frame with remeshing at every time step.
  - `annulus.cpp`: C++ code for simulating the active Ericksen-Leslie hydrodynamics in an Eulerian frame for annular geometry.
  - `complexGeometries.cpp`: C++ code for simulating the active Ericksen-Leslie hydrodynamics in Eulerian frame for complex 3D geometry using point cloud data.

## Please cite the following related publications if you use these codes.
1. **A numerical solver for active hydrodynamics in three dimensions and its application to active turbulence**
   - **Authors:** Abhinav Singh, Philipp H. Suhrcke, Pietro Incardona, Ivo F. Sbalzarini.
   - **Published in:** Physics of Fluids, 1 October 2023.
   - **Volume:** 35, **Issue:** 10, **Article:** 105155.
   - **DOI:** [10.1063/5.0169546](https://doi.org/10.1063/5.0169546).
   - **Related Information:** This publication presents a higher-order convergent numerical solver for active polar hydrodynamics in three-dimensional domains of arbitrary shape. It involves a scalable open-source software implementation for shared- and distributed-memory parallel computers, enabling the computational study of the nonlinear dynamics of out-of-equilibrium materials from first principles.
  
2. **A C++ expression system for partial differential equations enables generic simulations of biological hydrodynamics**
   - **Published in:** European Physical Journal E (EPJ E)
   - **Year:** 2021
   - **Volume:** 44
   - **Article Number:** 117
   - **DOI:** [10.1140/epje/s10189-021-00121-x](https://doi.org/10.1140/epje/s10189-021-00121-x).


3. **OpenFPM: A scalable open framework for particle and particle-mesh codes on parallel computers**
   - **Authors:** Pietro Incardona et al.
   - **Published in:** Computer Physics Communications, 2019.
   - **DOI:** [10.1016/j.cpc.2019.03.007](https://doi.org/10.1016/j.cpc.2019.03.007).
