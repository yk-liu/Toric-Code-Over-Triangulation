# Fortran Code for Toric Code over Triangulations

This repository is for the Fortran Code for computing the lowest eigenvalues and eigenvectors over triangulations. This is the appendix for [A Computer Code for Topological Quantum Spin Systems over Triangulated Surfaces](https://arxiv.org/abs/1912.12964). 

# Contents

The `EToricCode.f90` file computes the full eigenvalues of the Hamiltonian. 

The `PToricCode.f90` file computes lowest `p=5` eigenvalues of the Hamiltonian. This code is paralellized with OpenMP directives.

The `FToricCode.f90` file computes lowest `p=5` eigenvalues of the Hamiltonian. The difference between `PToricCode.f90` and `FToricCode.f90` is that the `FToricCode.f90` does not allocate the large `gpsi`s, making  the tasks inside the loop over Γ explicitly independent.

In the `Sample Inputs` folder contains input files listed in the paper. The Tetrahedron, Octahedron and Icosahedron are triangulation data for the sphere. The Torus1 is the minimal triangulation for a genus 1 torus and Torus2 is a conventional triangulation for the torus. The DoubleTorus is the minimal triangulation of the genus 2 torus or a double-holed torus. The illustrations for these triangulations are in the paper.

The `InputCheck` is for checking your version of triangulation data. It does not compute the genus of your data, but checks if each edge has exactly two vertices and two neighboring triangles, and if each triangle has exactly three egdes. (For large triangulations the data can be painful to check). Note that even if you passed this test there is no garantee that the data is indeed correct. Also, if you modify the code for other calculations orther than triangulation this file might not apply.

# How to run the code

To run the code, install Intel Fortran from the [official website](https://software.intel.com/en-us/fortran-compilers). After [setting up the correct environment](https://software.intel.com/en-us/articles/setting-up-the-build-environment-for-using-intel-c-or-fortran-compilers), you can run with the following command:

```
ifort -xHost -parallel [inputfilename] -o [outputfilename] -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -qopenmp

./outputfilename
```

To use your own triangulation data, modify the `include [filename]` in the three main code files and change the parameter`p` for the number of lowest eigenvalues you wish you compute.

# License

This code is free to use, modify, and redistribute for academic purposes. Please cite the original paper https://arxiv.org/abs/1912.12964.
