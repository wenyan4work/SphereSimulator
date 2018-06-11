# SimToolBox
A toolbox for object-tracking simulation. This toolbox contains a collection of useful tools for different tasks involved in a particle-tracking type simulation.

At the base level the code relies on `sctl` and `Trilinos`. The folder `sctl` contains the Scientific Computing Template Library (`sctl`), developed by Dhairya Malhotra, and publicly available at `http://git.dhairyamalhotra.com/dmalhotra/SCTL`. The version included here is slightly modified. The folder `Trilinos` contains an interface to the huge C++ project for distributed linear algebra and some other stuff.

The folder `Collision` contains the routines for LCP based collision resolution algorithms, supporting both OpenMP and MPI through `Trilinos`.

The folder `MPI` contains the routines for radius-based neighbor searching and syncing, supporting both OpenMP and MPI through `sctl`.

The folder `Protein`, `Sphere`, and `Sylinder` contains geometric primitives for both collision and Boundary Integral applications. Parallel IO is also supported for each geometric primitive through the XML-based vtk routines.

The folder `Util` contains some utilities to facilitate the application, including an interface to linear algebra `Eigen`, a parallel RNG interface to `TRNG`, a command line parser `cmdparser` (https://github.com/FlorianRappl/CmdParser), a timer, and a Gauss_Legendre table (https://github.com/sivaramambikasaran/Quadrature).

New functions will be continuously added to this toolbox to facilitate quick development of HPC simulation code.