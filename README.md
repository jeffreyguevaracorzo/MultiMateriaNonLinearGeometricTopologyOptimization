# MultiMaterialNonLinearGeometricTopologyOptimization-FORTRAN

This repository contains a series of FORTRAN routines to perform a multi-material topology optimization using FEM in non-linear geometric elasticity. Unlike the previous codes of this series, which handled either several materials in linear elasticity or a single material under large displacements, this one covers both at the same time: an arbitrary number of candidate materials, each with its own Young's modulus and its own Poisson ratio, distributed over a 2D or 3D domain whose equilibrium is solved with geometrically non-linear kinematics.

The code operates by inputting structure and load information into the "input" folder and defining the characteristics of the problem in the MainMultiMaterialNL.f90 file. It utilizes the MA86 library from HSL (with an adapted interface), as well as the LAPACK, BLAS, Metis and OpenMP libraries. The optimization is carried out with the method of moving asymptotes (MMA) presented by Professor Svanberg (1987), with one volume constraint per material phase. This code is for academic purposes only.

The material distribution follows the mapping-based interpolation function of Yi et al. (2023) and Zheng et al. (2024), in which the ratio between the p-norm and the 1-norm of the design variables of an element decides which material occupies it. The extension to large displacements is done at the constitutive level: the scalar SIMP factor is replaced by an interpolation of the constitutive tensors of the competing phases,

```
D_e = D_void + SUM_i psi_i^q ( D_i - D_void )
```

which is enough because in the Total Lagrangian / St. Venant-Kirchhoff setting both the internal force vector and the tangent stiffness are linear in D. Interpolating D is therefore equivalent to interpolating the stored energy of the phases, and the Newton-Raphson kernel, the sparse assembly and the solver are reused without changes. The sensitivities are obtained with the adjoint method, solving K_T*lambda = f_ext, since the geometrically non-linear compliance problem is not self-adjoint and the usual shortcut of replacing the adjoint field by the displacement field is not valid here. In the linear limit lambda equals u and the classical self-adjoint expression is recovered.

A .makefile is included to facilitate compilation and execution. During the post-processing stage, the code generates files in EnSight format (https://dav.lbl.gov/archive/NERSC/Software/ensight/), which can be read using ParaView, an open-source tool (https://www.paraview.org/). One density field is written per material phase, plus a material index field where 0 is void and i is the i-th material. The code was based on the work of Liu & Tovar (2014) and Andreassen et al. (2011), considering the FEM theory presented by Prof. Bathe (2006). Any suggestions for improving the code performance or reports of bugs are welcome.

### A note on HSL

The HSL library is free for academic use but is licensed individually and cannot be redistributed, so its sources are **not** included in this repository. To compile the code you need to request your own academic licence at https://www.hsl.rl.ac.uk/ and place the files hsl_ma86d.f90, hsl_ma86s.f90, common.f, common90.f90 and sdeps90.f90 into src/solver/. The wrapper Solver_MA86Module.f90, which is part of this project, is included. That wrapper exposes a single function, SparseSystemMA86Solver(rows, cols, values, rhs), so using MUMPS, PARDISO or any other sparse symmetric solver instead only requires reimplementing it.

### Notes on running the code

Compile with a plain `make`, not `make -j`: the HSL sources define modules used by sibling files in the same folder and the parallel build races on the .mod files.

The projection parameter beta is doubled every fixed number of design iterations. With the default schedule, running fewer than 50 iterations means beta never grows, the design stays grey and the post-processing filter removes every element. Use 150 to 300 iterations to get a converged 0/1 layout.

Watch the residual of the Newton-Raphson iteration. It should fall quadratically, for example 1e-1 -> 1e-5 -> 1e-10. If it halves at every iteration instead, something is inconsistent between the tangent stiffness and the internal force vector. If the message "Max FEM iterations reached" appears in many design iterations, the low-density elements are distorting and that run should be treated as under-converged.

A ready-to-run 3D case with eight-node hexahedra is provided in tests/cantilever3D_hexa8.

### Verification

The implementation was checked in two independent ways. Setting the number of materials to one reproduces the response of the single-material non-linear code to 15 significant digits. The adjoint sensitivities were compared against central finite differences over twelve randomly chosen (element, material) pairs, with relative errors between 3e-8 and 6e-6. So far the tests cover 2D quad4 and 3D hexa8 elements with the Total Lagrangian formulation; the Updated Lagrangian formulation and the remaining element types are implemented but not yet covered.

### References

* Andreassen, E., Clausen, A., Schevenels, M., Lazarov, B. S., & Sigmund, O. (2011). Efficient topology optimization in MATLAB using 88 lines of code. Structural and Multidisciplinary Optimization, 43, 1-16.
* Bathe, K. J. (2006). Finite element procedures. Klaus-Jurgen Bathe.
* Bendsoe, M. P., & Sigmund, O. (2013). Topology optimization: theory, methods, and applications. Springer Science & Business Media.
* Buhl, T., Pedersen, C. B. W., & Sigmund, O. (2000). Stiffness design of geometrically nonlinear structures using topology optimization. Structural and Multidisciplinary Optimization, 19, 93-104.
* Liu, K., & Tovar, A. (2014). An efficient 3D topology optimization code written in Matlab. Structural and Multidisciplinary Optimization, 50, 1175-1196.
* Svanberg, K. (1987). The method of moving asymptotes—a new method for structural optimization. International Journal for Numerical Methods in Engineering, 24(2), 359-373.
* Wang, F., Lazarov, B. S., Sigmund, O., & Jensen, J. S. (2014). Interpolation scheme for fictitious domain techniques and topology optimization of finite strain elastic problems. Computer Methods in Applied Mechanics and Engineering, 276, 453-472.
* Yi, B., Yoon, G. H., Zheng, R., Liu, L., Li, D., & Peng, X. (2023). A unified material interpolation for topology optimization of multi-materials. Computers & Structures, 282, 107041.
* Zheng, R., Yi, B., Peng, X., & Yoon, G. H. (2024). An efficient code for the multi-material topology optimization of 2D/3D continuum structures written in Matlab. Applied Sciences, 14(2), 657.
