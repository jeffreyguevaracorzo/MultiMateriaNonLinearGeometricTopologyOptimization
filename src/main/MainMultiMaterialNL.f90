program MainMultiMaterialNonLinear
    ! 1. Modules
    use MMNL_Optimization_Module
    use MMNL_Paraview_Module
    implicit none
    type (MMNLTop)                                       :: TopoModel
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !          MULTI-MATERIAL TOPOLOGY OPTIMIZATION WITH GEOMETRIC NONLINEARITY                                  !
    !                                     (TOP-MM-NL)                                                            !
    !                                                                                                            !
    !  Fusion of two previously separate codes:                                                                  !
    !    * multi-material topology optimization in linear elasticity, using the mapping-based                    !
    !      interpolation function of Zheng, Yi, Peng & Yoon (Appl. Sci. 2024, 14, 657)                           !
    !    * single-material topology optimization with geometric nonlinearity (Total/Updated                      !
    !      Lagrangian, St. Venant-Kirchhoff, Newton-Raphson)                                                     !
    !                                                                                                            !
    !  The fusion happens at the CONSTITUTIVE level: the scalar SIMP factor rho^p is replaced by                 !
    !         D_e = D_void + sum_i Psi_i(e)^n ( D_i - D_void )                                                   !
    !  which works unchanged inside the nonlinear kernel because both the tangent stiffness and                  !
    !  the internal force vector are linear in D. See InterpolatedD in Base_FEA_MMNL_Module.f90.                 !
    !                                                                                                            !
    !  Sensitivities are computed with the ADJOINT method (K_T lambda = f_ext), which is required                !
    !  because the geometrically nonlinear compliance problem is NOT self-adjoint.                               !
    !                                                                                                            !
    !  Requires: LAPACK, BLAS, METIS, OpenMP and the HSL MA86 library (bundled in src/solver).                   !
    !  Academic use only.                                                                                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(unit=*, fmt=*) '1. Structure properties'
    ! ---------------- 1.1 Mechanical / mesh parameters ----------------
    call SetAnalysisType(TopoModel,'PlaneStress')       ! PlaneStress(2D), PlaneStrain(2D), SolidIso(3D)
    call SetElementType(TopoModel,'quad4')              ! tria3, tria6, quad4, quad8, tetra4, tetra10, hexa8, hexa20
    call SetThickness(TopoModel,50.0d0)                 ! Only for 2D
    call SetGaussAprox(TopoModel,2)

    ! ---------------- 1.2 Nonlinear FE analysis ----------------
    ! Formulation and number of load increments are now INDEPENDENT switches.
    ! 'TL' is path-independent, so Incremental > 1 is purely a robustness aid for
    ! Newton-Raphson (very useful once large voids appear in the design domain).
    ! 'UL' requires Incremental > 1 by construction.
    call SetFormulation(TopoModel,'TL')                 ! 'TL' or 'UL'
    call SetLoadIncremental(TopoModel,1)                ! number of load steps
    call SetMaxFEAIteration(TopoModel,40)               ! max. Newton-Raphson iterations
    call SetConvergenceTolerance(TopoModel,1.0d-4)       ! residual tolerance

    ! ---------------- 1.3 Multi-material properties (Eq. 4) ----------------
    call SetNMaterial(TopoModel,2)
    call SetMaterialProperties(TopoModel,[30000.0d0, 50000.0d0])      ! E of each material
    call SetPoissonModulusMaterial(TopoModel,[0.30d0, 0.30d0])      ! nu of each material, SAME ORDER
    call SetVoidModulus(TopoModel,1.0d-3)                           ! E_void (very small)
    call SetPNorm(TopoModel,6.0d0)                                  ! p (p-norm exponent)
    call SetDeltaMap(TopoModel,1.0d-9)                              ! delta, avoids 0/0
    call SetPenalFactorTO(TopoModel,3.0d0)                          ! n (penalization exponent)
    call SetVolFractionMaterial(TopoModel,[1.0d0/6.0d0, 1.0d0/6.0d0])   ! volume fraction per material

    ! ---------------- 1.4 Filter and projection (Eq. 2-3) ----------------
    call SetFilterRadiusTO(TopoModel,10.0d0*1.5d0)      ! filter radius (1.5x element size)
    call SetBetaParameters(TopoModel,1.0d0,64.0d0,50)   ! beta0, betaMax, iterations before doubling beta

    ! ---------------- 1.5 Optimization control ----------------
    call SetMaxIterationsTO(TopoModel,100)
    call SetOptimizationToleranceTO(TopoModel,0.01d0)
    call SetPostProcesingFilterTO(TopoModel,0.5d0)
    call ReadFiles(TopoModel)

    ! 2. Running the optimization
    write(unit=*, fmt=*) '2. Multi-material nonlinear topology optimization'
    call TopologyOptimizationProcessMultiMaterial(TopoModel)

    ! 3. Post-processing
    write(unit=*, fmt=*) '3. Postprocessing (Paraview)'
    call MMNLParaviewPostProcessing(TopoModel,'output/paraview','MultiMaterialNLResult')
end program MainMultiMaterialNonLinear
