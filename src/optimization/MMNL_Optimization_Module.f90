module MMNL_Optimization_Module
    ! -------------------------------------------------------------------------------------- !
    !  MULTI-MATERIAL TOPOLOGY OPTIMIZATION MODULE                                           !
    !  Implements the mapping-based interpolation function for multi-material topology       !
    !  optimization, as described in:                                                        !
    !    Zheng, R.; Yi, B.; Peng, X.; Yoon, G.-H. "An Efficient Code for the Multi-Material  !
    !    Topology Optimization of 2D/3D Continuum Structures Written in Matlab."             !
    !    Appl. Sci. 2024, 14, 657. https://doi.org/10.3390/app14020657                       !
    !  (interpolation originally proposed by Yi, B. et al., Comput. Struct. 2023, 282)       !
    !                                                                                        !
    !  This module extends the existing FEM_MMNL type (linear elasticity FEA) and re-uses the     !
    !  HSL-MA87 solver and the generic MMA solver already present in the code base. It is    !
    !  a self-contained addition: the original single-material Optimization_module is        !
    !  untouched.                                                                            !
    !                                                                                        !
    !  Mapping between the paper's equations and the subroutines below:                      !
    !    Eq. (2)  Density filter                -> ForwardFilterField / BuildFilterMulti     !
    !    Eq. (3)  tanh/Heaviside projection      -> ProjectionStep                           !
    !    Eq. (4)  mapping-based interpolation    -> InterpolationStep (Psi, EEff)            !
    !    Eq. (5)  per-material volume/mass       -> InterpolationStep (MassMaterial)         !
    !    Eq. (6)  optimization model             -> TopologyOptimizationProcessMultiMaterial !
    !    Eq. (7)-(9)  compliance sensitivities   -> ComplianceSensitivity                    !
    !    Eq. (10)-(11) volume sensitivities      -> handled directly in the main loop        !
    !    Eq. (12) projection derivative          -> ProjectionDerivative                     !
    !    Eq. (13) filter derivative (adjoint)    -> AdjointFilterField                       !
    ! -------------------------------------------------------------------------------------- !
    use MMA_Module
    use FEA_MMNL_Module
    implicit none

    type, extends(FEM_MMNL)                                          :: MMNLTop
        ! ---------------- material data ----------------
        ! NMaterial, EMaterial, PoissonMaterial, EVoid, PoissonVoid and DMaterial now
        ! live in FEM_MMNL_Base (Base_FEA_MMNL_Module.f90), because the finite element
        ! kernel itself needs them every Newton iteration.
        ! The linear code's K0Material(Ne,NMaterial+1,ndof,ndof) is GONE: with
        ! geometric nonlinearity the element stiffness depends on u and cannot be
        ! precomputed, so the interpolation is applied to the small constitutive
        ! tensors D instead -- exact, and orders of magnitude cheaper in memory.
        ! ---------------- interpolation / SIMP-like parameters ----------------
        double precision                                        :: PNorm                ! p-norm exponent (Eq. 4), paper suggests p=6
        double precision                                        :: DeltaMap             ! delta, avoids 0/0 (Eq. 4), ~1e-9
        double precision                                        :: PenalFactor          ! penalization exponent n (Eq. 4)
        double precision, dimension(:), allocatable             :: VolFractionMaterial  ! allowed volume fraction of each material (NMaterial)
        ! ---------------- filter ----------------
        double precision                                        :: FilterRadius
        integer, dimension(:), allocatable                      :: FilterRow            ! sparse filter, COO format (built once)
        integer, dimension(:), allocatable                      :: FilterCol
        double precision, dimension(:), allocatable             :: FilterVal
        double precision, dimension(:), allocatable             :: FilterHs             ! row sums (normalization)
        ! ---------------- projection (Eq. 3) ----------------
        double precision                                        :: Beta
        double precision                                        :: BetaMax
        integer                                                 :: BetaUpdateIter
        integer                                                 :: LoopBeta
        double precision                                        :: Eta
        ! ---------------- optimization control ----------------
        integer                                                 :: Iteration
        integer                                                 :: MaxTOIteration
        double precision                                        :: OptimizationTolerance
        double precision                                        :: Change
        double precision                                        :: PostProcesingFilter
        ! ---------------- design fields, all (Ne,NMaterial) ----------------
        double precision, dimension(:,:), allocatable           :: xDes                 ! raw design variables x
        double precision, dimension(:,:), allocatable           :: xTilde               ! filtered field
        double precision, dimension(:,:), allocatable           :: xProj                ! projected (Heaviside) field, x~~ in the paper
        double precision, dimension(:,:), allocatable           :: Psi                  ! mapping function psi_i (Eq. 4)
        double precision, dimension(:), allocatable             :: VolumePerElement
        double precision                                        :: Compliance
        double precision                                        :: ComplianceRef   ! iter-1 compliance, used to scale f0val/df0dx for MMA
        double precision, dimension(:), allocatable             :: MassMaterial         ! current volume of each material
        double precision, dimension(:), allocatable             :: MassMaterialTarget
        ! ---------------- MMA arrays ----------------
        integer                                                 :: mm
        integer                                                 :: nn
        double precision, dimension(:,:), allocatable           :: xmin
        double precision, dimension(:,:), allocatable           :: xmax
        double precision, dimension(:,:), allocatable           :: XMMA
        double precision, dimension(:,:), allocatable           :: xold1
        double precision, dimension(:,:), allocatable           :: xold2
        double precision, dimension(:,:), allocatable           :: low
        double precision, dimension(:,:), allocatable           :: upp
        double precision                                        :: f0val
        double precision, dimension(:,:), allocatable           :: df0dx
        double precision, dimension(:,:), allocatable           :: fval
        double precision, dimension(:,:), allocatable           :: dfdx
        double precision                                        :: a0
        double precision, dimension(:,:), allocatable           :: a
        double precision, dimension(:,:), allocatable           :: c
        double precision, dimension(:,:), allocatable           :: d
        ! ---------------- post-processing ----------------
        integer, dimension(:), allocatable                      :: MaterialIndex        ! 0 = void, i = material i
    contains
        procedure                                               :: SetPNorm
        procedure                                               :: SetDeltaMap
        procedure                                               :: SetVolFractionMaterial
        procedure                                               :: SetFilterRadiusTO
        procedure                                               :: SetPenalFactorTO
        procedure                                               :: SetMaxIterationsTO
        procedure                                               :: SetOptimizationToleranceTO
        procedure                                               :: SetBetaParameters
        procedure                                               :: SetPostProcesingFilterTO
        procedure                                               :: TopologyOptimizationProcessMultiMaterial
    end type MMNLTop

contains

    ! ----------------------------------------------------------------- !
    !                              SETTERS                              !
    ! ----------------------------------------------------------------- !



    ! One Poisson's ratio per material, same shape/usage as SetMaterialProperties,
    ! e.g. call SetPoissonModulusMaterial(TopoModel,[0.3d0, 0.3d0, 0.3d0]).


    subroutine SetPNorm(Self,PNorm)
        implicit none
        class(MMNLTop), intent(inout)                                  :: Self
        double precision, intent(in)                                 :: PNorm
        Self%PNorm = PNorm
    end subroutine SetPNorm

    subroutine SetDeltaMap(Self,DeltaMap)
        implicit none
        class(MMNLTop), intent(inout)                                  :: Self
        double precision, intent(in)                                 :: DeltaMap
        Self%DeltaMap = DeltaMap
    end subroutine SetDeltaMap

    subroutine SetVolFractionMaterial(Self,VolFractionMaterial)
        implicit none
        class(MMNLTop), intent(inout)                                  :: Self
        double precision, dimension(:), intent(in)                   :: VolFractionMaterial
        allocate(Self%VolFractionMaterial(size(VolFractionMaterial)))
        Self%VolFractionMaterial = VolFractionMaterial
    end subroutine SetVolFractionMaterial

    subroutine SetFilterRadiusTO(Self,FilterRadius)
        implicit none
        class(MMNLTop), intent(inout)                                  :: Self
        double precision, intent(in)                                 :: FilterRadius
        Self%FilterRadius = FilterRadius
    end subroutine SetFilterRadiusTO

    subroutine SetPenalFactorTO(Self,PenalFactor)
        implicit none
        class(MMNLTop), intent(inout)                                  :: Self
        double precision, intent(in)                                 :: PenalFactor
        Self%PenalFactor = PenalFactor
    end subroutine SetPenalFactorTO

    subroutine SetMaxIterationsTO(Self,MaxTOIteration)
        implicit none
        class(MMNLTop), intent(inout)                                  :: Self
        integer, intent(in)                                          :: MaxTOIteration
        Self%MaxTOIteration = MaxTOIteration
    end subroutine SetMaxIterationsTO

    subroutine SetOptimizationToleranceTO(Self,OptimizationTolerance)
        implicit none
        class(MMNLTop), intent(inout)                                  :: Self
        double precision, intent(in)                                 :: OptimizationTolerance
        Self%OptimizationTolerance = OptimizationTolerance
    end subroutine SetOptimizationToleranceTO

    subroutine SetBetaParameters(Self,BetaInit,BetaMax,BetaUpdateIter)
        implicit none
        class(MMNLTop), intent(inout)                                  :: Self
        double precision, intent(in)                                 :: BetaInit
        double precision, intent(in)                                 :: BetaMax
        integer, intent(in)                                          :: BetaUpdateIter
        Self%Beta = BetaInit
        Self%BetaMax = BetaMax
        Self%BetaUpdateIter = BetaUpdateIter
        Self%Eta = 0.5d0
        Self%LoopBeta = 0
    end subroutine SetBetaParameters

    subroutine SetPostProcesingFilterTO(Self,PostProcesingFilter)
        implicit none
        class(MMNLTop), intent(inout)                                  :: Self
        double precision, intent(in)                                 :: PostProcesingFilter
        Self%PostProcesingFilter = PostProcesingFilter
    end subroutine SetPostProcesingFilterTO

    ! ----------------------------------------------------------------- !
    !                      VOLUME PER ELEMENT (as in Optimization_module) !
    ! ----------------------------------------------------------------- !
    subroutine VolumeAllElementTO(Self)
        implicit none
        class(MMNLTop), intent(inout)                                  :: Self
        integer                                                      :: i
        double precision                                             :: AreaAux
        double precision, dimension(:), allocatable                  :: v1,v2,v3,v4,Centroid
        double precision, dimension(:,:), allocatable                :: Coordinates
        allocate(Self%VolumePerElement(Self%Ne))
        Self%VolumePerElement = 0.0d0
        do i = 1, Self%Ne, 1
            AreaAux = 0.0d0
            if ((Self%ElementType.eq.'tria3').or.(Self%ElementType.eq.'tria6')) then
                Coordinates = Self%Coordinates(Self%ConnectivityN(i,:),:)
                v1 = [Coordinates(3,:) - Coordinates(2,:),0.0d0]
                v2 = [Coordinates(1,:) - Coordinates(2,:),0.0d0]
                AreaAux = norm2(CrossProduct(v1,v2))/2.0d0
                Self%VolumePerElement(i) = AreaAux*Self%Thickness
            elseif ((Self%ElementType.eq.'quad4').or.(Self%ElementType.eq.'quad8')) then
                Coordinates = Self%Coordinates(Self%ConnectivityN(i,:),:)
                v1 = [Coordinates(2,:) - Coordinates(1,:),0.0d0]
                v2 = [Coordinates(4,:) - Coordinates(1,:),0.0d0]
                v3 = [Coordinates(4,:) - Coordinates(3,:),0.0d0]
                v4 = [Coordinates(2,:) - Coordinates(3,:),0.0d0]
                AreaAux = norm2(CrossProduct(v1,v2))/2.0d0 + norm2(CrossProduct(v3,v4))/2.0d0
                Self%VolumePerElement(i) = AreaAux*Self%Thickness
            elseif ((Self%ElementType.eq.'tetra4').or.(Self%ElementType.eq.'tetra10')) then
                Coordinates = Self%Coordinates(Self%ConnectivityN(i,:),:)
                v1 = Coordinates(3,:) - Coordinates(1,:)
                v2 = Coordinates(2,:) - Coordinates(1,:)
                v3 = Coordinates(4,:) - Coordinates(1,:)
                Self%VolumePerElement(i) = abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
            elseif ((Self%ElementType.eq.'hexa8').or.(Self%ElementType.eq.'hexa20')) then
                Coordinates = Self%Coordinates(Self%ConnectivityN(i,:),:)
                Centroid = sum(Coordinates,1)/size(Coordinates(:,1))
                v1 = Coordinates(1,:) - Centroid; v2 = Coordinates(2,:) - Centroid; v3 = Coordinates(5,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                v1 = Coordinates(2,:) - Centroid; v2 = Coordinates(5,:) - Centroid; v3 = Coordinates(6,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                v1 = Coordinates(3,:) - Centroid; v2 = Coordinates(7,:) - Centroid; v3 = Coordinates(8,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                v1 = Coordinates(3,:) - Centroid; v2 = Coordinates(4,:) - Centroid; v3 = Coordinates(8,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                v1 = Coordinates(2,:) - Centroid; v2 = Coordinates(6,:) - Centroid; v3 = Coordinates(7,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                v1 = Coordinates(2,:) - Centroid; v2 = Coordinates(3,:) - Centroid; v3 = Coordinates(7,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                v1 = Coordinates(1,:) - Centroid; v2 = Coordinates(5,:) - Centroid; v3 = Coordinates(8,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                v1 = Coordinates(1,:) - Centroid; v2 = Coordinates(4,:) - Centroid; v3 = Coordinates(8,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                v1 = Coordinates(1,:) - Centroid; v2 = Coordinates(2,:) - Centroid; v3 = Coordinates(3,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                v1 = Coordinates(1,:) - Centroid; v2 = Coordinates(3,:) - Centroid; v3 = Coordinates(4,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                v1 = Coordinates(5,:) - Centroid; v2 = Coordinates(6,:) - Centroid; v3 = Coordinates(7,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                v1 = Coordinates(5,:) - Centroid; v2 = Coordinates(7,:) - Centroid; v3 = Coordinates(8,:) - Centroid
                Self%VolumePerElement(i) = Self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
            end if
        end do
    end subroutine VolumeAllElementTO





    ! ----------------------------------------------------------------- !
    !   PER-MATERIAL STIFFNESS BASES (geometry-only, built ONCE)         !
    ! ----------------------------------------------------------------- !
    ! Each material -- and the void phase -- can now have its own Poisson's
    ! ratio, so a single shared "K0" (built with one Young's modulus and one
    ! Poisson's ratio, then scaled by a scalar EEff per element) is no longer
    ! valid: for a fixed nu, D(E) = E*D(1) is linear in E, but D is NOT linear
    ! across different nu's. So instead we build one full stiffness base per
    ! material (its real E_i and nu_i), plus one for the void phase (EVoid and
    ! PoissonModulusVoid), by reusing the existing single-material GetKlocal
    ! with a uniform unit density field (density=1, PenalFactor=1 => no-op).
    ! These bases depend ONLY on the mesh geometry -- not on the design
    ! variables -- so they are computed once, here, before the optimization
    ! loop, and simply re-combined every iteration in GetKlocalMultiMaterial.





    ! ----------------------------------------------------------------- !
    !     SPARSE DENSITY FILTER (Eq. 2), built ONCE (COO format)         !
    ! ----------------------------------------------------------------- !
    subroutine BuildFilterMulti(Self)
        implicit none
        class(MMNLTop), intent(inout)                                  :: Self
        integer                                                      :: i,j,cnt,nnz
        double precision                                             :: dist,w
        double precision, dimension(:,:), allocatable                :: PosPromE
        allocate(PosPromE(Self%Ne,Self%DimAnalysis))
        do i = 1, Self%Ne, 1
            PosPromE(i,:) = Sum(Self%Coordinates(Self%ConnectivityN(i,:),:),1)/Self%Npe
        end do
        ! pass 1: count non-zero entries
        nnz = 0
        !$omp parallel do default(none) shared(Self,PosPromE) private(i,j,dist) reduction(+:nnz)
        do i = 1, Self%Ne, 1
            do j = 1, Self%Ne, 1
                dist = norm2(PosPromE(i,:) - PosPromE(j,:))
                if (dist.le.Self%FilterRadius) nnz = nnz + 1
            end do
        end do
        !$omp end parallel do
        allocate(Self%FilterRow(nnz),Self%FilterCol(nnz),Self%FilterVal(nnz))
        allocate(Self%FilterHs(Self%Ne));  Self%FilterHs = 0.0d0
        ! pass 2: fill (serial, to keep insertion order deterministic and simple)
        cnt = 0
        do i = 1, Self%Ne, 1
            do j = 1, Self%Ne, 1
                dist = norm2(PosPromE(i,:) - PosPromE(j,:))
                w = Self%FilterRadius - dist
                if (w.ge.0.0d0) then
                    cnt = cnt + 1
                    Self%FilterRow(cnt) = i
                    Self%FilterCol(cnt) = j
                    Self%FilterVal(cnt) = w
                    Self%FilterHs(i) = Self%FilterHs(i) + w
                end if
            end do
        end do
        deallocate(PosPromE)
    end subroutine BuildFilterMulti

    ! x~_e = sum_j H_ej x_j / Hs_e   (Eq. 2, forward filter of one material's field)
    subroutine ForwardFilterField(Self,FieldIn,FieldOut)
        implicit none
        class(MMNLTop), intent(in)                                     :: Self
        double precision, dimension(:), intent(in)                   :: FieldIn
        double precision, dimension(:), allocatable, intent(inout)   :: FieldOut
        integer                                                      :: k
        if (.not.allocated(FieldOut)) allocate(FieldOut(Self%Ne))
        FieldOut = 0.0d0
        do k = 1, size(Self%FilterRow), 1
            FieldOut(Self%FilterRow(k)) = FieldOut(Self%FilterRow(k)) + Self%FilterVal(k)*FieldIn(Self%FilterCol(k))
        end do
        FieldOut = FieldOut/Self%FilterHs
    end subroutine ForwardFilterField

    ! sens_j = sum_e H_ej * (v_e/Hs_e)   (Eq. 13, adjoint used in the sensitivity chain rule)
    subroutine AdjointFilterField(Self,FieldIn,FieldOut)
        implicit none
        class(MMNLTop), intent(in)                                     :: Self
        double precision, dimension(:), intent(in)                   :: FieldIn
        double precision, dimension(:), allocatable, intent(inout)   :: FieldOut
        double precision, dimension(:), allocatable                  :: w
        integer                                                      :: k
        allocate(w(Self%Ne))
        w = FieldIn/Self%FilterHs
        if (.not.allocated(FieldOut)) allocate(FieldOut(Self%Ne))
        FieldOut = 0.0d0
        do k = 1, size(Self%FilterRow), 1
            ! H is symmetric (H_ej depends only on distance), so row/col can be swapped
            FieldOut(Self%FilterCol(k)) = FieldOut(Self%FilterCol(k)) + Self%FilterVal(k)*w(Self%FilterRow(k))
        end do
        deallocate(w)
    end subroutine AdjointFilterField

    ! ----------------------------------------------------------------- !
    !         PROJECTION (Eq. 3) and its derivative (Eq. 12)             !
    ! ----------------------------------------------------------------- !
    elemental double precision function ProjectionFunction(xt,beta,eta) result(xp)
        implicit none
        double precision, intent(in)                                 :: xt,beta,eta
        xp = (tanh(beta*eta) + tanh(beta*(xt-eta)))/(tanh(beta*eta)+tanh(beta*(1.0d0-eta)))
    end function ProjectionFunction

    elemental double precision function ProjectionDerivative(xt,beta,eta) result(dxp)
        implicit none
        double precision, intent(in)                                 :: xt,beta,eta
        dxp = beta*(1.0d0-tanh(beta*(xt-eta))**2)/(tanh(beta*eta)+tanh(beta*(1.0d0-eta)))
    end function ProjectionDerivative

    ! ----------------------------------------------------------------- !
    !                             MMA PARAMETERS                        !
    ! ----------------------------------------------------------------- !
    subroutine MMAParametersMulti(Self)
        implicit none
        class(MMNLTop), intent(inout)                                  :: Self
        double precision                                             :: xminval,xmaxval
        integer                                                      :: i
        xminval = 1.0d-3
        xmaxval = 1.0d0
        Self%mm = Self%NMaterial
        Self%nn = Self%Ne*Self%NMaterial
        allocate(Self%xmin(Self%nn,1));    Self%xmin = xminval
        allocate(Self%xmax(Self%nn,1));    Self%xmax = xmaxval
        allocate(Self%XMMA(Self%nn,1))
        allocate(Self%xold1(Self%nn,1))
        allocate(Self%xold2(Self%nn,1))
        allocate(Self%low(Self%nn,1));     Self%low = xminval
        allocate(Self%upp(Self%nn,1));     Self%upp = xmaxval
        allocate(Self%fval(Self%mm,1))
        allocate(Self%dfdx(Self%mm,Self%nn))
        allocate(Self%df0dx(Self%nn,1))
        allocate(Self%a(Self%mm,1),Self%c(Self%mm,1),Self%d(Self%mm,1))
        do i = 1, Self%NMaterial, 1
            Self%xold1((i-1)*Self%Ne+1:i*Self%Ne,1) = Self%xDes(:,i)
            Self%xold2((i-1)*Self%Ne+1:i*Self%Ne,1) = Self%xDes(:,i)
        end do
        Self%a0 = 1.0d0   ! Svanberg standard (the linear multi-material code used 0.0d0)
        Self%a = 0.0d0
        Self%c = 1.0d6
        Self%d = 0.0d0
    end subroutine MMAParametersMulti

    ! ----------------------------------------------------------------- !
    !   INTERPOLATION (Eq. 4-5): psi_i and per-material mass            !
    ! ----------------------------------------------------------------- !
    ! Note: the linear code also assembled a scalar "effective Young modulus"
    ! EEff(e) here. That is gone: with per-material Poisson ratios a single
    ! scalar cannot represent the interpolated constitutive behaviour, and the
    ! nonlinear kernel builds the full interpolated tensor D_e on the fly
    ! (InterpolatedD in Base_FEA_MMNL_Module). Only Psi is needed here.
    subroutine InterpolationStep(Self)
        implicit none
        class(MMNLTop), intent(inout)                                  :: Self
        integer                                                      :: e,i
        double precision                                             :: npnorm,s1norm
        if (.not.allocated(Self%Psi)) allocate(Self%Psi(Self%Ne,Self%NMaterial))
        do e = 1, Self%Ne, 1
            npnorm = (sum(Self%xProj(e,:)**Self%PNorm))**(1.0d0/Self%PNorm)
            s1norm = sum(Self%xProj(e,:))
            Self%Psi(e,:) = Self%xProj(e,:)*npnorm/(s1norm+Self%DeltaMap)
        end do
        ! mass of each material, weighted by the physical volume of each element (Eq. 5,
        ! generalized to non-uniform meshes: sum_e xProj(e,i) * VolumePerElement(e))
        if (.not.allocated(Self%MassMaterial)) allocate(Self%MassMaterial(Self%NMaterial))
        do i = 1, Self%NMaterial, 1
            Self%MassMaterial(i) = sum(Self%xProj(:,i)*Self%VolumePerElement)
        end do
    end subroutine InterpolationStep

    ! ----------------------------------------------------------------- !
    !   COMPLIANCE SENSITIVITY w.r.t. the PROJECTED field (Eq. 7-9)     !
    ! ----------------------------------------------------------------- !
    ! ADJOINT sensitivity of the end-compliance C = f_ext^T u.
    !
    ! Equilibrium:            R(u,x) = f_int(u,x) - f_ext = 0
    ! Implicit derivative:    K_T du/dx = -d(f_int)/dx
    ! Chain rule:             dC/dx = f_ext^T du/dx = -lambda^T d(f_int)/dx ,  K_T lambda = f_ext
    !
    ! In LINEAR elasticity f_int = K u and K_T = K, so lambda = u and this collapses
    ! to the familiar self-adjoint expression -u^T (dK/dx) u. With geometric
    ! nonlinearity lambda /= u and the self-adjoint shortcut is simply WRONG -- this
    ! is the single most important difference with respect to the linear code.
    !
    ! Since f_int,e is linear in D_e and D_e = D_void + sum_m Psi_m^n (D_m - D_void):
    !     d(f_int,e)/d(Psi_m) = n Psi_m^(n-1) ( Fphase(e,:,m) - Fphase(e,:,void) )
    ! so, per element e and material n:
    !     dC/d(xProj_n) = - sum_m  n Psi_m^(n-1) (dPsi_m/dxProj_n) *
    !                              lambda_e^T ( Fphase(e,:,m) - Fphase(e,:,void) )
    !
    ! The dPsi_m/dxProj_n Jacobian (Eq. 9) is unchanged from the linear code: it is
    ! pure algebra of the mapping function and knows nothing about the physics. Only
    ! the physical factor in the last parenthesis changes.
    subroutine ComplianceSensitivity(Self,dCompdxProj)
        implicit none
        class(MMNLTop), intent(inout)                                :: Self
        double precision, dimension(:,:), allocatable, intent(inout) :: dCompdxProj
        integer                                                      :: e,m,n
        double precision                                             :: npnorm,s1norm,dpsi_mn,PhysFac
        double precision, dimension(:), allocatable                  :: LambdaE
        if (.not.allocated(dCompdxProj)) allocate(dCompdxProj(Self%Ne,Self%NMaterial))
        dCompdxProj = 0.0d0
        allocate(LambdaE(Self%Npe*Self%DimAnalysis))
        do e = 1, Self%Ne, 1
            LambdaE = Self%Lambda(Self%ConnectivityD(e,:))
            npnorm = (sum(Self%xProj(e,:)**Self%PNorm))**(1.0d0/Self%PNorm)
            s1norm = sum(Self%xProj(e,:))
            do n = 1, Self%NMaterial, 1
                do m = 1, Self%NMaterial, 1
                    ! d(psi_m)/d(xProj_n), full Jacobian (diagonal + cross terms), Eq. 9
                    dpsi_mn = Self%xProj(e,m)*(Self%xProj(e,n)**(Self%PNorm-1.0d0))*(npnorm**(1.0d0-Self%PNorm)) &
                                /(s1norm+Self%DeltaMap) &
                              - Self%xProj(e,m)*npnorm/((s1norm+Self%DeltaMap)**2)
                    if (m.eq.n) dpsi_mn = dpsi_mn + npnorm/(s1norm+Self%DeltaMap)
                    ! physical factor: lambda_e^T ( f_int,e^(m) - f_int,e^(void) )
                    PhysFac = dot_product(LambdaE, Self%FphaseE(e,:,m) - Self%FphaseE(e,:,Self%NMaterial+1))
                    dCompdxProj(e,n) = dCompdxProj(e,n) - Self%PenalFactor*(Self%Psi(e,m)**(Self%PenalFactor-1.0d0)) &
                                        *dpsi_mn*PhysFac
                end do
            end do
        end do
        deallocate(LambdaE)
    end subroutine ComplianceSensitivity

    ! ----------------------------------------------------------------- !
    !   per-material unit-density strain energy: Ue^T K0Material_i Ue,  !
    !   plus the ACTUAL per-element strain energy Ue^T KLocal_e Ue,     !
    !   whose sum over elements equals the compliance c = U^T K U.      !
    !   (Both are exact regardless of whether nu differs per material,  !
    !    since they use the real per-material stiffness bases and the   !
    !    real combined element matrix -- no shared-D0/scalar-EEff       !
    !    shortcut needed anymore.)                                      !
    ! ----------------------------------------------------------------- !

    ! ----------------------------------------------------------------- !
    !                        FINAL TOPOLOGY (post-processing)           !
    ! ----------------------------------------------------------------- !
    subroutine FinalTopologyMulti(Self)
        implicit none
        class(MMNLTop), intent(inout)                                  :: Self
        integer                                                      :: e,imax
        if (.not.allocated(Self%MaterialIndex)) allocate(Self%MaterialIndex(Self%Ne))
        do e = 1, Self%Ne, 1
            imax = maxloc(Self%xProj(e,:),1)
            if (Self%xProj(e,imax).ge.Self%PostProcesingFilter) then
                Self%MaterialIndex(e) = imax
            else
                Self%MaterialIndex(e) = 0
            end if
        end do
        !call FilePrinting(Self%MaterialIndex,'V','DataResults/MultiMaterialIndex.txt')
        !call FilePrinting(Self%xProj,'DataResults/MultiMaterialDensity.txt')
    end subroutine FinalTopologyMulti

    ! ----------------------------------------------------------------- !
    !                    PRINTING CONVERGENCE INFO                       !
    ! ----------------------------------------------------------------- !
    subroutine PrintingConvergenceMulti(Self)
        implicit none
        class(MMNLTop), intent(inout)                                  :: Self
        write(unit=*, fmt=*) 'Ite,',Self%Iteration,'Beta,',Self%Beta,'Obj(compliance),',Self%Compliance, &
                                'Change,',Self%Change,'MassFraction,', &
                                Self%MassMaterial/sum(Self%VolumePerElement)
    end subroutine PrintingConvergenceMulti

    ! ----------------------------------------------------------------- !
    !                     MAIN OPTIMIZATION LOOP                         !
    ! ----------------------------------------------------------------- !
    subroutine TopologyOptimizationProcessMultiMaterial(Self)
        implicit none
        class(MMNLTop), intent(inout)                                  :: Self
        integer                                                      :: i,e
        double precision                                             :: TotalVolume
        double precision, dimension(:,:), allocatable                :: dCompdxProj
        double precision, dimension(:), allocatable                  :: temp1,temp2
        write(unit=*, fmt=*) '1. Multi-material pre-assembly'
        call PreAssemblyRoutine(Self)
        call VolumeAllElementTO(Self)
        call BuildFilterMulti(Self)
        ! constitutive tensor of every phase (+ void): built ONCE, they do not
        ! depend on the design variables nor on the deformation
        call BuildMaterialDTensors(Self)
        TotalVolume = sum(Self%VolumePerElement)
        allocate(Self%MassMaterialTarget(Self%NMaterial))
        Self%MassMaterialTarget = Self%VolFractionMaterial*TotalVolume
        ! initial design variables (0.5 for every material, as in the reference paper)
        allocate(Self%xDes(Self%Ne,Self%NMaterial));   Self%xDes = 0.5d0
        allocate(Self%xTilde(Self%Ne,Self%NMaterial))
        allocate(Self%xProj(Self%Ne,Self%NMaterial))
        call MMAParametersMulti(Self)
        Self%Iteration = 1
        Self%Change = 1.0d0
        write(unit=*, fmt=*) '---------- multi-material topology optimization process ----------'
        do
            ! -------- Eq. 2: density filter (per material) --------
            do i = 1, Self%NMaterial, 1
                temp1 = Self%xDes(:,i)
                call ForwardFilterField(Self,temp1,temp2)
                Self%xTilde(:,i) = temp2
            end do
            ! -------- Eq. 3: Heaviside/tanh projection --------
            Self%xProj = ProjectionFunction(Self%xTilde,Self%Beta,Self%Eta)
            ! -------- Eq. 4-5: interpolation function and per-material mass --------
            call InterpolationStep(Self)
            ! -------- GEOMETRICALLY NONLINEAR FE-Analysis --------
            ! Replaces the single linear solve of the reference code. NewtonRaphson
            ! also leaves behind everything the sensitivity analysis needs at the
            ! converged state: the per-phase internal forces FphaseE and the adjoint
            ! field Lambda.
            write(unit=*, fmt=*) '-Solving the multi-material NONLINEAR FE system'
            call NewtonRaphson(Self,Self%Psi,Self%PenalFactor)
            ! end-compliance C = f_ext^T u  (in the linear limit this equals U^T K U)
            Self%Compliance = dot_product(Self%FGlobal_ext_total,Self%UGlobal)
            ! -------- Eq. 7-9: compliance sensitivity w.r.t. projected field --------
            call ComplianceSensitivity(Self,dCompdxProj)
            ! -------- Eq. 12-13: chain rule through projection and filter --------
            do i = 1, Self%NMaterial, 1
                temp1 = dCompdxProj(:,i)*ProjectionDerivative(Self%xTilde(:,i),Self%Beta,Self%Eta)
                call AdjointFilterField(Self,temp1,temp2)
                Self%df0dx((i-1)*Self%Ne+1:i*Self%Ne,1) = temp2
                Self%XMMA((i-1)*Self%Ne+1:i*Self%Ne,1)  = Self%xDes(:,i)
            end do
            ! -------- scale objective (Sec. 3.4 of the reference paper: "the objective
            !          function and its sensitivity are scaled to balance the objective
            !          and constraints"). Without this, MMA's internal dual subproblem
            !          solver can divide by (numerically) zero on the very first call
            !          whenever compliance and volume constraints differ by orders of
            !          magnitude, producing NaN in xDes/EEff and a spuriously "singular"
            !          stiffness matrix on the NEXT FE solve.
            if (Self%Iteration.eq.1) Self%ComplianceRef = max(Self%Compliance,1.0d-12)
            Self%f0val = 100.0d0*Self%Compliance/Self%ComplianceRef
            Self%df0dx = 100.0d0*Self%df0dx/Self%ComplianceRef
            ! -------- volume constraints g_i = Mass_i - target_i <= 0 (normalized) --------
            Self%dfdx = 0.0d0
            do i = 1, Self%NMaterial, 1
                Self%fval(i,1) = 100.0d0*(Self%MassMaterial(i)-Self%MassMaterialTarget(i))/Self%MassMaterialTarget(i)
                temp1 = Self%VolumePerElement*ProjectionDerivative(Self%xTilde(:,i),Self%Beta,Self%Eta)
                call AdjointFilterField(Self,temp1,temp2)
                Self%dfdx(i,(i-1)*Self%Ne+1:i*Self%Ne) = 100.0d0*temp2/Self%MassMaterialTarget(i)
            end do
            ! -------- MMA update --------
            write(unit=*, fmt=*) '-Getting new solution (MMA)'
            call MMA_Solver_Interface(Self%mm,Self%nn,Self%Iteration,Self%XMMA,Self%xmin,Self%xmax, &
                                       Self%xold1,Self%xold2,Self%f0val,Self%df0dx,Self%fval,Self%dfdx, &
                                       Self%low,Self%upp,Self%a0,Self%a,Self%c,Self%d)
            Self%xold2 = Self%xold1
            do i = 1, Self%NMaterial, 1
                Self%xold1((i-1)*Self%Ne+1:i*Self%Ne,1) = Self%xDes(:,i)
                Self%Change = max(Self%Change, maxval(abs(Self%XMMA((i-1)*Self%Ne+1:i*Self%Ne,1)-Self%xDes(:,i))))
                Self%xDes(:,i) = Self%XMMA((i-1)*Self%Ne+1:i*Self%Ne,1)
            end do
            call PrintingConvergenceMulti(Self)
            ! -------- Heaviside continuation --------
            Self%LoopBeta = Self%LoopBeta + 1
            if ((Self%Beta.lt.Self%BetaMax).and.(Self%LoopBeta.ge.Self%BetaUpdateIter)) then
                Self%Beta = 2.0d0*Self%Beta
                Self%LoopBeta = 0
                write(unit=*, fmt=*) 'Parameter beta increased to ', Self%Beta
            end if
            if (Self%Change.lt.Self%OptimizationTolerance) then; write(unit=*, fmt=*) 'Precision achieved'       ; exit; end if
            if (Self%Iteration.ge.Self%MaxTOIteration)     then; write(unit=*, fmt=*) 'Max TO iterations reached'; exit; end if
            Self%Iteration = Self%Iteration + 1
            Self%Change = 0.0d0
        end do
        call FinalTopologyMulti(Self)
    end subroutine TopologyOptimizationProcessMultiMaterial

end module MMNL_Optimization_Module