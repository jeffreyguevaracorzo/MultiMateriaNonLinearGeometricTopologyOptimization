module FEA_MMNL_Module
    use Base_FEA_MMNL_Module
    ! Type variable for all the information of the structure
    type, extends(FEM_MMNL_Base)                          :: FEM_MMNL
        double precision, dimension(:), allocatable     :: FGlobal_ext      ! Total Global load (external)
        double precision, dimension(:), allocatable     :: FGlobal_ext_inc  ! Incremental Global load (external)
        double precision, dimension(:), allocatable     :: FGlobal_int      ! Total Global load (internal)
        double precision, dimension(:), allocatable     :: FGlobal_res      ! Total Global load (residual)
        double precision, dimension(:), allocatable     :: UGlobalTotal     ! Final displacement vector
        double precision, dimension(:), allocatable     :: UGlobal          ! Final global displacement vector
        double precision, dimension(:), allocatable     :: dUGlobal         ! Delta global displacement vector
        double precision, dimension(:,:,:), allocatable :: KLocal_NonLinear ! Local stiffness matrix per element
        ! performing variables
        integer, dimension(:), allocatable              :: Rows_KGlobal     ! Indx-rows of global stiffness matrix (reduced)
        integer, dimension(:), allocatable              :: Cols_KGlobal     ! Indx-cols of global stiffness matrix (reduced)
        double precision, dimension(:), allocatable     :: value_KGlobal    ! global stiffness vector
        double precision, dimension(:), allocatable     :: value_FGlobal    ! Global load vector
        double precision, dimension(:), allocatable     :: value_UGlobal    ! Global displacement vector
        ! Results
        double precision, dimension(:,:), allocatable   :: Displacement     ! Global dispacement matrix
        double precision, dimension(:,:), allocatable   :: StrainE          ! Strain per element
        double precision, dimension(:,:), allocatable   :: StressE          ! Stress per element
        double precision, dimension(:), allocatable     :: StrainEnergyE    ! Strain energy per element (UNPENALIZED, see note)
        ! ---------------- adjoint sensitivity analysis ----------------
        ! In geometrically nonlinear compliance minimization the problem is NOT
        ! self-adjoint: the linear-elasticity shortcut lambda = u is invalid.
        ! Minimizing the end-compliance C = f_ext^T u subject to the equilibrium
        ! residual R(u,x) = f_int(u,x) - f_ext = 0 gives
        !       dC/dx = - lambda^T * d(f_int)/dx ,   with   K_T * lambda = f_ext
        ! and lambda = u only in the linear case. K_T is the tangent stiffness at
        ! the CONVERGED state, so the adjoint costs one extra assembly + solve per
        ! topology-optimization iteration (~7% overhead over a 15-iteration Newton).
        double precision, dimension(:), allocatable     :: FGlobal_ext_total ! full external load (adjoint RHS)
        double precision, dimension(:), allocatable     :: Lambda            ! adjoint field
        ! FphaseE(Ne,ndof,NMaterial+1): internal force vector each element WOULD have
        ! if it were made entirely of phase m, evaluated at the converged u. Since
        ! f_int is linear in D, d(f_int,e)/d(Psi_m) = n*Psi_m^(n-1)*(Fphase_m - Fphase_void).
        double precision, dimension(:,:,:), allocatable :: FphaseE
    contains
        procedure                                       :: Initialize
        procedure                                       :: GetLocalTangentK
        procedure                                       :: GetLocalInternalF
        procedure                                       :: AssemblyNonLinearSystem
        procedure                                       :: SolveNonLinearSystem
        procedure                                       :: ReleasingMemory
        procedure                                       :: ProcessingResults
        procedure                                       :: GetPhaseInternalForces
        procedure                                       :: SolveAdjointSystem
        procedure                                       :: NewtonRaphson
    end type FEM_MMNL

    contains
    ! ------------ FINITE ELEMENT ANALYSISS FUNCTIONS AND SUBROUTINES ------------
    ! 1. Initialize system
    subroutine Initialize(Self)
        implicit none
        class(FEM_MMNL), intent(inout)                             :: Self
        if (.not.allocated(Self%UGlobal))     allocate(Self%UGlobal(Self%N*Self%DimAnalysis))
        if (.not.allocated(Self%FGlobal_res)) allocate(Self%FGlobal_res(Self%N*Self%DimAnalysis))
        Self%UGlobal = 0.0d0
        Self%FGlobal_res = 0.0d0
        Self%FGlobal_ext = (Self%FGlobal_DL + Self%FGlobal_PL)/Self%Incremental
        Self%FGlobal_ext_inc = (Self%FGlobal_DL + Self%FGlobal_PL)/Self%Incremental
        ! full external load, kept for the adjoint right-hand side (FGlobal_ext is
        ! overwritten by the incremental loading loop and must not be used for this)
        if (.not.allocated(Self%FGlobal_ext_total)) allocate(Self%FGlobal_ext_total(Self%N*Self%DimAnalysis))
        Self%FGlobal_ext_total = Self%FGlobal_DL + Self%FGlobal_PL
        ! in UL the mesh is convected, so it must be reset to the reference configuration
        Self%MatCoordinates = Self%Coordinates
        !call FilePrinting(Self%FGlobal_ext,'V','DataResults/.InternalData/FGlobal_ext.txt')
    end subroutine Initialize
    ! 2. Local Stiffness matrices
    subroutine GetLocalTangentK(Self,Psi,PenalFactor)
        implicit none
        class(FEM_MMNL), intent(inout)                             :: Self
        double precision, dimension(:,:), allocatable, intent(in)     :: Psi
        double precision, intent(in)                                :: PenalFactor
        ! internal variables
        integer                                                     :: el,i,j,k,l,m
        double precision                                            :: e,n,z,w1,w2,w3,DetJacobian,FAC
        double precision, dimension(:), allocatable                 :: Strain_ea,Strain_gl,StressL,StressNL
        double precision, dimension(:,:), allocatable               :: Jacobian,InvJacobian,D,DiffN,DiffNXY,IXY,S,H_L,H_NL
        double precision, dimension(:,:), allocatable               :: F,Cright,Cleft,Egl,Eea,IMatrix,B,BN,BG,BO
        double precision, dimension(:,:), allocatable               :: ElementCoordinates,ElementDisplacement
        if (Self%DimAnalysis.eq.2) then
            allocate(Self%KLocal_NonLinear(Self%Ne,Self%Npe*2,Self%Npe*2)); Self%KLocal_NonLinear = 0.0d0
            allocate(B(3,2*Self%Npe));                                                          B = 0.0d0
            allocate(BO(3,2*Self%Npe));                                                        BO = 0.0d0
            allocate(BN(3,2*Self%Npe));                                                        BN = 0.0d0
            allocate(BG(Self%DimAnalysis**2,2*Self%Npe));                                      BG = 0.0d0
            allocate(ElementCoordinates(Self%Npe,2));                          ElementCoordinates = 0.0d0
            allocate(ElementDisplacement(Self%Npe,2));                        ElementDisplacement = 0.0d0
            allocate(Strain_ea(3));                                                     Strain_ea = 0.0d0
            allocate(Strain_gl(3));                                                     Strain_gl = 0.0d0
            allocate(S(2,2));                                                                   S = 0.0d0
            allocate(H_L(4,4));                                                               H_L = 0.0d0
            allocate(H_NL(4,4));                                                             H_NL = 0.0d0
            allocate(IMatrix(2,2));                                                       IMatrix = 0.0d0
            IMatrix(1,1) = 1.0d0; IMatrix(2,2) = 1.0d0;
            do el = 1, Self%Ne, 1
                ! Multi-material interpolated constitutive tensor (Eq. 4)
                call InterpolatedD(Self,Psi(el,:),PenalFactor,D)
                if (Self%Formulation.eq.'TL') then
                    ElementCoordinates = Self%Coordinates(Self%ConnectivityN(el,:),1:2)
                else
                    ElementCoordinates = Self%MatCoordinates(Self%ConnectivityN(el,:),1:2)
                end if
                ElementDisplacement = reshape(Self%UGlobal(Self%ConnectivityD(el,:)),[2,Self%Npe])
                do i = 1, Self%QuadGauss, 1
                    e = Self%GaussPoint(i)
                    w1 = Self%GaussWeights(i)
                    do j = 1, Self%QuadGauss, 1
                        n = Self%GaussPoint(j)
                        w2 = Self%GaussWeights(j)
                        call DiffFormFunction(Self,DiffN,e,n)
                        Jacobian = matmul(DiffN,ElementCoordinates)
                        InvJacobian = Inverse(Jacobian)
                        DetJacobian = Determinant(Jacobian)
                        DiffNXY = matmul(InvJacobian,DiffN)
                        !IXY = matmul(DiffNXY,transpose(ElementDisplacement))
                        IXY = matmul(ElementDisplacement,transpose(DiffNXY))
                        F = IXY + IMatrix
                        Cleft = matmul(F,transpose(F))           ! Left Green Tensor
                        cleft = inverse(cleft)
                        Cright = matmul(transpose(F),F)          ! Right Green Tensor
                        Eea = 0.5d0*(IMatrix - Cleft)            ! Euler-Almansi Strain
                        Strain_ea = [Eea(1,1), Eea(2,2), 2.0d0*Eea(1,2)]
                        StressL = matmul(D,Strain_ea)
                        Egl = 0.5d0*(Cright - IMatrix)           ! Green-Lagrange Strain
                        Strain_gl = [Egl(1,1), Egl(2,2), 2.0d0*Egl(1,2)]
                        StressNL = matmul(D,Strain_gl)
                        ! BO - lineal
                        do k = 1, size(DiffN,2), 1
                            BO(1,2*k-1) = DiffNXY(1,k)
                            BO(2,2*k) = DiffNXY(2,k)
                            BO(3,2*k-1) = DiffNXY(2,k)
                            BO(3,2*k) = DiffNXY(1,k)
                        end do
                        ! BN - no lineal
                        do k = 1, Self%Npe, 1
                            BN(1,2*k-1) = F(1,1)*DiffNXY(1,k)
                            BN(1,2*k) = F(2,1)*DiffNXY(1,k)
                            BN(2,2*k-1) = F(1,2)*DiffNXY(2,k)
                            BN(2,2*k) = F(2,2)*DiffNXY(2,k)
                            BN(3,2*k-1) = F(1,2)*DiffNXY(1,k) + F(1,1)*DiffNXY(2,k)
                            BN(3,2*k) = F(2,2)*DiffNXY(1,k) + F(2,1)*DiffNXY(2,k)
                        end do
                        ! ---------------------------------------------------------------- !
                        !  CORRECTION: in the Total Lagrangian formulation the variation of  !
                        !  the Green-Lagrange strain is  dE = B_L du  with                   !
                        !      dE_11 = F_11 dN/dx du_1 + F_21 dN/dx du_2   (etc.)            !
                        !  i.e. B_L is built ENTIRELY from the deformation gradient F.       !
                        !  That is exactly what BN contains, so BN ALONE is the complete     !
                        !  B_L. BO is the small-strain matrix B_0, and since F = I + grad(u) !
                        !  the matrix BN ALREADY CONTAINS B_0 (its F = I part).              !
                        !  Writing B = BO + BN therefore counts B_0 twice: in the undeformed !
                        !  limit it gives B = 2*B_0, hence f_int = 2*f_int_exact and         !
                        !  K = 4*K_exact. Newton then converges only LINEARLY with rate 1/2  !
                        !  (residual halving every iteration) towards 2*K_0*u = f, i.e. HALF !
                        !  the correct displacement. Bathe (FEP, 6.3.2) writes the same      !
                        !  matrix split as B_L = B_L0 + B_L1 using the displacement gradient !
                        !  l_ij = du_i/dx_j; expressed through F_ij = delta_ij + l_ij the    !
                        !  two terms collapse into the single matrix BN. Mixing the two      !
                        !  notations is what produced the double count.                      !
                        ! ---------------------------------------------------------------- !
                        B = BN
                        ! BG - sigma
                        do k = 1, Self%Npe, 1
                            BG(1,2*k-1) = DiffNXY(1,k)
                            BG(2,2*k-1) = DiffNXY(2,k)
                            BG(3,2*k) = DiffNXY(1,k)
                            BG(4,2*k) = DiffNXY(2,k)
                        end do
                        ! H
                        H_L = 0.0d0; H_NL = 0.0d0;
                        S(1,:) = [StressL(1), StressL(3)]
                        S(2,:) = [StressL(3), StressL(2)]
                        H_L(1:2,1:2) = S
                        H_L(3:4,3:4) = S
                        S(1,:) = [StressNL(1), StressNL(3)]
                        S(2,:) = [StressNL(3), StressNL(2)]
                        H_NL(1:2,1:2) = S
                        H_NL(3:4,3:4) = S
                        ! K Local
                        FAC = DetJacobian*w1*w2*(Self%Thickness)
                        if (Self%Formulation.eq.'TL') then ! for TL formulation (Total lagrangian)
                            Self%KLocal_NonLinear(el,:,:) = Self%KLocal_NonLinear(el,:,:) &                 ! base
                                                        + FAC*(matmul(transpose(B),matmul(D,B))) &          ! Kmat
                                                        + FAC*(matmul(transpose(BG),matmul(H_NL,BG)))       ! Kgeo
                        else ! for UL formulation (Incremental)
                            Self%KLocal_NonLinear(el,:,:) = Self%KLocal_NonLinear(el,:,:) &                 ! base
                                                        + FAC*(matmul(transpose(BO),matmul(D,BO))) &        ! Kmat
                                                        + FAC*(matmul(transpose(BG),matmul(H_L,BG)))        ! Kgeo
                        end if
                        deallocate(DiffN)
                    end do
                end do
                deallocate(D)
            end do
        elseif(Self%DimAnalysis.eq.3) then
            allocate(Self%KLocal_NonLinear(Self%Ne,Self%Npe*3,Self%Npe*3)); Self%KLocal_NonLinear = 0.0d0
            allocate(B(6,3*Self%Npe));                                                          B = 0.0d0
            allocate(BO(6,3*Self%Npe));                                                        BO = 0.0d0
            allocate(BN(6,3*Self%Npe));                                                        BN = 0.0d0
            allocate(BG(Self%DimAnalysis**2,3*Self%Npe));                                      BG = 0.0d0
            allocate(ElementCoordinates(Self%Npe,3));                          ElementCoordinates = 0.0d0
            allocate(ElementDisplacement(Self%Npe,3));                        ElementDisplacement = 0.0d0
            allocate(Strain_ea(6));                                                     Strain_ea = 0.0d0
            allocate(Strain_gl(6));                                                     Strain_gl = 0.0d0
            allocate(S(3,3));                                                                   S = 0.0d0
            allocate(H_L(9,9));                                                               H_L = 0.0d0
            allocate(H_NL(9,9));                                                             H_NL = 0.0d0
            allocate(IMatrix(3,3));                                                       IMatrix = 0.0d0
            IMatrix(1,1) = 1.0d0; IMatrix(2,2) = 1.0d0; IMatrix(3,3) = 1.0d0;
            do el = 1, Self%Ne, 1
                ! Multi-material interpolated constitutive tensor (Eq. 4)
                call InterpolatedD(Self,Psi(el,:),PenalFactor,D)
                if (Self%Formulation.eq.'TL') then
                    ElementCoordinates = Self%Coordinates(Self%ConnectivityN(el,:),1:3)
                else
                    ElementCoordinates = Self%MatCoordinates(Self%ConnectivityN(el,:),1:3)
                end if
                ElementDisplacement = reshape(Self%UGlobal(Self%ConnectivityD(el,:)),[3,Self%Npe])
                do i = 1, Self%QuadGauss, 1
                    e = Self%GaussPoint(i)
                    w1 = Self%GaussWeights(i)
                    do j = 1, Self%QuadGauss, 1
                        n = Self%GaussPoint(j)
                        w2 = Self%GaussWeights(j)
                        do k = 1, Self%QuadGauss, 1
                            z = Self%GaussPoint(k)
                            w3 = Self%GaussWeights(k)
                            call DiffFormFunction(Self,DiffN,e,n,z)
                            Jacobian = matmul(DiffN,ElementCoordinates)
                            InvJacobian = Inverse(Jacobian)
                            DetJacobian = Determinant(Jacobian)
                            DiffNXY = matmul(InvJacobian,DiffN)
                            !IXY = matmul(DiffNXY,transpose(ElementDisplacement))
                            IXY = matmul(ElementDisplacement,transpose(DiffNXY))
                            F = IXY + IMatrix
                            Cleft = matmul(F,transpose(F))           ! Left Green Tensor
                            cleft = inverse(cleft)
                            Cright = matmul(transpose(F),F)          ! Right Green Tensor
                            Eea = 0.5d0*(IMatrix - Cleft)            ! Euler-Almansi Strain
                            Strain_ea = [Eea(1,1), Eea(2,2), Eea(3,3), 2.0d0*Eea(1,2), 2.0d0*Eea(2,3), 2.0d0*Eea(1,3)]
                            StressL = matmul(D,Strain_ea)
                            Egl = 0.5d0*(Cright - IMatrix)           ! Green-Lagrange Strain
                            Strain_gl = [Egl(1,1), Egl(2,2), Egl(3,3), 2.0d0*Egl(1,2), 2.0d0*Egl(2,3), 2.0d0*Egl(1,3)]
                            StressNL = matmul(D,Strain_gl)
                            ! BO - lineal
                            do l = 1, size(DiffN,2), 1
                                BO(1,3*l-2) = DiffNxy(1,l)
                                BO(2,3*l-1) = DiffNxy(2,l)
                                BO(3,3*l) = DiffNxy(3,l)
                                BO(4,3*l-2) = DiffNxy(2,l)
                                BO(4,3*l-1) = DiffNxy(1,l)
                                BO(5,3*l-1) = DiffNxy(3,l)
                                BO(5,3*l) = DiffNxy(2,l)
                                BO(6,3*l-2) = DiffNxy(3,l)
                                BO(6,3*l) = DiffNxy(1,l)
                            end do
                            ! BN - no lineal
                            do l = 1, Self%Npe, 1
                                BN(1,3*l-2) = F(1,1)*DiffNXY(1,l)
                                BN(1,3*l-1) = F(2,1)*DiffNXY(1,l)
                                BN(1,3*l) = F(3,1)*DiffNXY(1,l)
                                BN(2,3*l-2) = F(1,2)*DiffNXY(2,l)
                                BN(2,3*l-1) = F(2,2)*DiffNXY(2,l)
                                BN(2,3*l) = F(2,3)*DiffNXY(2,l)
                                BN(3,3*l-2) = F(1,3)*DiffNXY(3,l)
                                BN(3,3*l-1) = F(2,3)*DiffNXY(3,l)
                                BN(3,3*l) = F(3,3)*DiffNXY(3,l)
                                BN(4,3*l-2) = F(1,2)*DiffNXY(1,l) + F(1,1)*DiffNXY(2,l)
                                BN(4,3*l-1) = F(2,2)*DiffNXY(1,l) + F(2,1)*DiffNXY(2,l)
                                BN(4,3*l) = F(3,2)*DiffNXY(1,l) + F(3,1)*DiffNXY(2,l)
                                BN(5,3*l-2) = F(1,3)*DiffNXY(2,l) + F(1,2)*DiffNXY(3,l)
                                BN(5,3*l-1) = F(2,3)*DiffNXY(2,l) + F(2,2)*DiffNXY(3,l)
                                BN(5,3*l) = F(3,3)*DiffNXY(2,l) + F(3,2)*DiffNXY(3,l)
                                BN(6,3*l-2) = F(1,3)*DiffNXY(1,l) + F(1,1)*DiffNXY(3,l)
                                BN(6,3*l-1) = F(2,3)*DiffNXY(1,l) + F(2,1)*DiffNXY(3,l)
                                BN(6,3*l) = F(3,3)*DiffNXY(1,l) + F(3,1)*DiffNXY(3,l)
                            end do
                            B = BN          ! TL: BN alone is the complete B_L (see 2D note)
                            ! BG - sigma
                            do l = 1, Self%Npe, 1
                                BG(1,3*l-2) = DiffNxy(1,l)
                                BG(2,3*l-2) = DiffNxy(2,l)
                                BG(3,3*l-2) = DiffNxy(3,l)
                                BG(4,3*l-1) = DiffNxy(1,l)
                                BG(5,3*l-1) = DiffNxy(2,l)
                                BG(6,3*l-1) = DiffNxy(3,l)
                                BG(7,3*l) = DiffNxy(1,l)
                                BG(8,3*l) = DiffNxy(2,l)
                                BG(9,3*l) = DiffNxy(3,l)
                            end do
                            ! H
                            H_L = 0.0d0; H_NL = 0.0d0
                            S(1,:) = [StressL(1), StressL(4), StressL(6)]
                            S(2,:) = [StressL(4), StressL(2), StressL(5)]
                            S(3,:) = [StressL(6), StressL(5), StressL(3)]
                            H_L(1:3,1:3) = S
                            H_L(4:6,4:6) = S
                            H_L(7:9,7:9) = S 
                            S(1,:) = [StressNL(1), StressNL(4), StressNL(6)]
                            S(2,:) = [StressNL(4), StressNL(2), StressNL(5)]
                            S(3,:) = [StressNL(6), StressNL(5), StressNL(3)]
                            H_NL(1:3,1:3) = S
                            H_NL(4:6,4:6) = S
                            H_NL(7:9,7:9) = S
                            ! K Local
                            FAC = DetJacobian*w1*w2*w3
                            if (Self%Formulation.eq.'TL') then ! for TL formulation (Total lagrangian)
                                Self%KLocal_NonLinear(el,:,:) = Self%KLocal_NonLinear(el,:,:) &                 ! base
                                                            + FAC*(matmul(transpose(B),matmul(D,B))) &          ! Kmat
                                                            + FAC*(matmul(transpose(BG),matmul(H_NL,BG)))       ! Kgeo
                            else ! for UL formulation (Incremental)
                                Self%KLocal_NonLinear(el,:,:) = Self%KLocal_NonLinear(el,:,:) &                 ! base
                                                            + FAC*(matmul(transpose(BO),matmul(D,BO))) &        ! Kmat
                                                            + FAC*(matmul(transpose(BG),matmul(H_L,BG)))        ! Kgeo
                            end if
                            deallocate(DiffN)
                        end do
                    end do
                end do
                deallocate(D)
            end do
        end if
        deallocate(Strain_ea,Strain_gl,S,H_L,H_NL,IMatrix,B,BN,BG,BO)
        deallocate(ElementCoordinates,ElementDisplacement)
    end subroutine GetLocalTangentK
    ! 3. Global Load Vector (global sparse form)
    subroutine GetLocalInternalF(Self,Psi,PenalFactor)
        implicit none
        class(FEM_MMNL), intent(inout)                             :: Self
        double precision, dimension(:,:), allocatable, intent(in)     :: Psi
        double precision, intent(in)                                :: PenalFactor
        ! internal variables
        integer                                                     :: el,i,j,k,l,m
        double precision                                            :: e,n,z,w1,w2,w3,DetJacobian,FAC
        double precision, dimension(:), allocatable                 :: Strain_ea,Strain_gl,StressL,StressNL
        double precision, dimension(:,:), allocatable               :: Jacobian,InvJacobian,D,DiffN,DiffNXY,IXY
        double precision, dimension(:,:), allocatable               :: F,Cright,Cleft,Egl,Eea,IMatrix,B,BN,BO
        double precision, dimension(:,:), allocatable               :: ElementCoordinates,ElementDisplacement
        if (Self%DimAnalysis.eq.2) then
            allocate(Self%FGlobal_int(Self%N*Self%DimAnalysis));                 Self%FGlobal_int = 0.0d0
            allocate(B(3,2*Self%Npe));                                                          B = 0.0d0
            allocate(BO(3,2*Self%Npe));                                                        BO = 0.0d0
            allocate(BN(3,2*Self%Npe));                                                        BN = 0.0d0
            allocate(ElementCoordinates(Self%Npe,2));                          ElementCoordinates = 0.0d0
            allocate(ElementDisplacement(Self%Npe,2));                        ElementDisplacement = 0.0d0
            allocate(Strain_ea(3));                                                     Strain_ea = 0.0d0
            allocate(Strain_gl(3));                                                     Strain_gl = 0.0d0
            allocate(IMatrix(2,2));                                                       IMatrix = 0.0d0
            IMatrix(1,1) = 1.0d0; IMatrix(2,2) = 1.0d0
            do el = 1, Self%Ne, 1
                ! Multi-material interpolated constitutive tensor (Eq. 4)
                call InterpolatedD(Self,Psi(el,:),PenalFactor,D)
                if (Self%Formulation.eq.'TL') then
                    ElementCoordinates = Self%Coordinates(Self%ConnectivityN(el,:),1:2)
                else
                    ElementCoordinates = Self%MatCoordinates(Self%ConnectivityN(el,:),1:2)
                end if
                ElementDisplacement = reshape(Self%UGlobal(Self%ConnectivityD(el,:)),[2,Self%Npe])
                do i = 1, Self%QuadGauss, 1
                    e = Self%GaussPoint(i)
                    w1 = Self%GaussWeights(i)
                    do j = 1, Self%QuadGauss, 1
                        n = Self%GaussPoint(j)
                        w2 = Self%GaussWeights(j)
                        call DiffFormFunction(Self,DiffN,e,n)
                        Jacobian = matmul(DiffN,ElementCoordinates)
                        InvJacobian = Inverse(Jacobian)
                        DetJacobian = Determinant(Jacobian)
                        DiffNXY = matmul(InvJacobian,DiffN)
                        !IXY = matmul(DiffNXY,transpose(ElementDisplacement))
                        IXY = matmul(ElementDisplacement,transpose(DiffNXY))
                        F = IXY + IMatrix
                        Cleft =  matmul(F,transpose(F))          ! Left Green Tensor
                        cleft = inverse(cleft)
                        Cright = matmul(transpose(F),F)          ! Right Green Tensor
                        Eea = 0.5d0*(IMatrix - Cleft)            ! Euler-Almansi Strain
                        Strain_ea = [Eea(1,1), Eea(2,2), 2.0d0*Eea(1,2)]
                        StressL = matmul(D,Strain_ea)
                        Egl = 0.5d0*(Cright - IMatrix)           ! Green-Lagrange Strain
                        Strain_gl = [Egl(1,1), Egl(2,2), 2.0d0*Egl(1,2)]
                        StressNL = matmul(D,Strain_gl)
                        ! BO - lineal
                        do k = 1, Self%Npe, 1
                            BO(1,2*k-1) = DiffNxy(1,k)
                            BO(2,2*k) = DiffNxy(2,k)
                            BO(3,2*k-1) = DiffNxy(2,k)
                            BO(3,2*k) = DiffNxy(1,k)
                        end do
                        ! BN - no lineal
                        do k = 1, Self%Npe, 1
                            BN(1,2*k-1) = F(1,1)*DiffNXY(1,k)
                            BN(1,2*k) = F(2,1)*DiffNXY(1,k)
                            BN(2,2*k-1) = F(1,2)*DiffNXY(2,k)
                            BN(2,2*k) = F(2,2)*DiffNXY(2,k)
                            BN(3,2*k-1) = F(1,2)*DiffNXY(1,k) + F(1,1)*DiffNXY(2,k)
                            BN(3,2*k) = F(2,2)*DiffNXY(1,k) + F(2,1)*DiffNXY(2,k)
                        end do
                        ! TL: B_L = BN (see the note in GetLocalTangentK).
                        ! UL: the tangent integrates B_0^T D B_0 on the CURRENT configuration,
                        !     so the internal force must use the SAME B_0 -- otherwise K_T is
                        !     not the derivative of f_int and Newton loses convergence. The
                        !     original code used BO+BN here while the UL tangent used BO.
                        if (Self%Formulation.eq.'TL') then
                            B = BN
                        else
                            B = BO
                        end if
                        FAC = DetJacobian*w1*w2*(Self%Thickness)
                        if (Self%Formulation.eq.'TL') then     ! for TL formulation (Total lagrangian)
                            Self%FGlobal_int(Self%ConnectivityD(el,:)) = Self%FGlobal_int(Self%ConnectivityD(el,:)) &
                                                                        + FAC*matmul(transpose(B),StressNL)
                        else ! for UL formulation (Incremental)
                            Self%FGlobal_int(Self%ConnectivityD(el,:)) = Self%FGlobal_int(Self%ConnectivityD(el,:)) &
                                                                        + FAC*matmul(transpose(B),StressL)
                        end if
                        deallocate(DiffN)
                    end do
                end do
                deallocate(D)
            end do
        elseif(Self%DimAnalysis.eq.3) then
            allocate(Self%FGlobal_int(Self%N*Self%DimAnalysis));                 Self%FGlobal_int = 0.0d0
            allocate(B(6,3*Self%Npe));                                                          B = 0.0d0
            allocate(BN(6,3*Self%Npe));                                                        BN = 0.0d0
            allocate(BO(6,3*Self%Npe));                                                        BO = 0.0d0
            allocate(ElementCoordinates(Self%Npe,3));                          ElementCoordinates = 0.0d0
            allocate(ElementDisplacement(Self%Npe,3));                        ElementDisplacement = 0.0d0
            allocate(Strain_ea(6));                                                     Strain_ea = 0.0d0
            allocate(Strain_gl(6));                                                     Strain_gl = 0.0d0
            allocate(IMatrix(3,3));                                                       IMatrix = 0.0d0
            IMatrix(1,1) = 1.0d0; IMatrix(2,2) = 1.0d0; IMatrix(3,3) = 1.0d0
            do el = 1, Self%Ne, 1
                ! Multi-material interpolated constitutive tensor (Eq. 4)
                call InterpolatedD(Self,Psi(el,:),PenalFactor,D)
                if (Self%Formulation.eq.'TL') then
                    ElementCoordinates = Self%Coordinates(Self%ConnectivityN(el,:),1:3)
                else
                    ElementCoordinates = Self%MatCoordinates(Self%ConnectivityN(el,:),1:3)
                end if
                ElementDisplacement = reshape(Self%UGlobal(Self%ConnectivityD(el,:)),[3,Self%Npe])
                do i = 1, Self%QuadGauss, 1
                    e = Self%GaussPoint(i)
                    w1 = Self%GaussWeights(i)
                    do j = 1, Self%QuadGauss, 1
                        n = Self%GaussPoint(j)
                        w2 = Self%GaussWeights(j)
                        do k = 1, Self%QuadGauss, 1
                            z = Self%GaussPoint(k)
                            w3 = Self%GaussWeights(k)
                            call DiffFormFunction(Self,DiffN,e,n,z)
                            Jacobian = matmul(DiffN,ElementCoordinates)
                            InvJacobian = Inverse(Jacobian)
                            DetJacobian = Determinant(Jacobian)
                            DiffNXY = matmul(InvJacobian,DiffN)
                            !IXY = matmul(DiffNXY,transpose(ElementDisplacement))
                            IXY = matmul(ElementDisplacement,transpose(DiffNXY))
                            F = IXY + IMatrix
                            Cleft =  matmul(F,transpose(F))          ! Left Green Tensor
                            cleft = inverse(cleft)
                            Cright = matmul(transpose(F),F)          ! Right Green Tensor
                            Eea = 0.5d0*(IMatrix - Cleft)            ! Euler-Almansi Strain
                            Strain_ea = [Eea(1,1), Eea(2,2), Eea(3,3), 2.0d0*Eea(1,2), 2.0d0*Eea(2,3), 2.0d0*Eea(1,3)]
                            StressL = matmul(D,Strain_ea)
                            Egl = 0.5d0*(Cright - IMatrix)           ! Green-Lagrange Strain
                            Strain_gl = [Egl(1,1), Egl(2,2), Egl(3,3), 2.0d0*Egl(1,2), 2.0d0*Egl(2,3), 2.0d0*Egl(1,3)]
                            StressNL = matmul(D,Strain_gl)
                            ! BO - lineal
                            do l = 1, size(DiffN,2), 1
                                BO(1,3*l-2) = DiffNxy(1,l)
                                BO(2,3*l-1) = DiffNxy(2,l)
                                BO(3,3*l) = DiffNxy(3,l)
                                BO(4,3*l-2) = DiffNxy(2,l)
                                BO(4,3*l-1) = DiffNxy(1,l)
                                BO(5,3*l-1) = DiffNxy(3,l)
                                BO(5,3*l) = DiffNxy(2,l)
                                BO(6,3*l-2) = DiffNxy(3,l)
                                BO(6,3*l) = DiffNxy(1,l)
                            end do
                            ! BN - no lineal
                            do l = 1, Self%Npe, 1
                                BN(1,3*l-2) = F(1,1)*DiffNXY(1,l)
                                BN(1,3*l-1) = F(2,1)*DiffNXY(1,l)
                                BN(1,3*l) = F(3,1)*DiffNXY(1,l)
                                BN(2,3*l-2) = F(1,2)*DiffNXY(2,l)
                                BN(2,3*l-1) = F(2,2)*DiffNXY(2,l)
                                BN(2,3*l) = F(3,2)*DiffNXY(2,l)
                                BN(3,3*l-2) = F(1,3)*DiffNXY(3,l)
                                BN(3,3*l-1) = F(2,3)*DiffNXY(3,l)
                                BN(3,3*l) = F(3,3)*DiffNXY(3,l)
                                BN(4,3*l-2) = F(1,2)*DiffNXY(1,l) + F(1,1)*DiffNXY(2,l)
                                BN(4,3*l-1) = F(2,2)*DiffNXY(1,l) + F(2,1)*DiffNXY(2,l)
                                BN(4,3*l) = F(3,2)*DiffNXY(1,l) + F(3,1)*DiffNXY(2,l)
                                BN(5,3*l-2) = F(1,3)*DiffNXY(2,l) + F(1,2)*DiffNXY(3,l)
                                BN(5,3*l-1) = F(2,3)*DiffNXY(2,l) + F(2,2)*DiffNXY(3,l)
                                BN(5,3*l) = F(3,3)*DiffNXY(2,l) + F(3,2)*DiffNXY(3,l)
                                BN(6,3*l-2) = F(1,3)*DiffNXY(1,l) + F(1,1)*DiffNXY(3,l)
                                BN(6,3*l-1) = F(2,3)*DiffNXY(1,l) + F(2,1)*DiffNXY(3,l)
                                BN(6,3*l) = F(3,3)*DiffNXY(1,l) + F(3,1)*DiffNXY(3,l)
                            end do
                            if (Self%Formulation.eq.'TL') then
                                B = BN
                            else
                                B = BO
                            end if
                            FAC = DetJacobian*w1*w2*w3
                            if (Self%Formulation.eq.'TL') then     ! for TL formulation (Total lagrangian)
                                Self%FGlobal_int(Self%ConnectivityD(el,:)) = Self%FGlobal_int(Self%ConnectivityD(el,:)) &
                                                                            + FAC*matmul(transpose(B),StressNL)
                            else ! for UL formulation (Incremental)
                                Self%FGlobal_int(Self%ConnectivityD(el,:)) = Self%FGlobal_int(Self%ConnectivityD(el,:)) &
                                                                            + FAC*matmul(transpose(B),StressL)
                            end if
                            deallocate(DiffN)
                        end do
                    end do
                end do
                deallocate(D)
            end do
        end if
        !call FilePrinting(Self%FGlobal_int,'V','DataResults/.InternalData/FGlobal_int.txt')
        deallocate(Strain_ea,Strain_gl,IMatrix,B,BN,BO)
        deallocate(ElementCoordinates,ElementDisplacement)
    end subroutine GetLocalInternalF
    ! 4. Global Stifness Matrix (global sparse form)
    subroutine AssemblyNonLinearSystem(Self)
        implicit none
        class(FEM_MMNL), intent(inout)                             :: Self
        integer                                                     :: n,i,j,k
        logical, dimension(:), allocatable                          :: C1,C2,InLogical
        ! the global stiffness matrix is assembled here 
        n = size(Self%index_i_KGlobal)
        allocate(Self%value_KGlobal(n))
        Self%value_KGlobal = 0.0d0

        do i = 1, Self%Ne, 1
            do j = 1, Self%DimAnalysis*Self%Npe, 1
                do k = 1, Self%DimAnalysis*Self%Npe, 1
                    n = Self%Location_KGlobal(i,j,k)    
                    if (n.eq.0) cycle   ! doesnt exist 
                    Self%value_KGlobal(n) = Self%value_KGlobal(n) + Self%KLocal_NonLinear(i,j,k)
                end do
            end do
        end do

        ! eliminating the zero elements
        C1 = Self%value_KGlobal.ne.0.0d0
        C2 = Self%index_i_KGlobal.eq.Self%index_j_KGlobal
        InLogical = C1.or.C2
        Self%Rows_KGlobal = pack(Self%index_i_KGlobal,InLogical)
        Self%Cols_KGlobal = pack(Self%index_j_KGlobal,InLogical)
        Self%value_KGlobal = pack(Self%value_KGlobal,InLogical)

        !call FilePrinting(Self%Rows_KGlobal,'V','DataResults/.InternalData/Rows_KGlobal.txt')
        !call FilePrinting(Self%Cols_KGlobal,'V','DataResults/.InternalData/Cols_KGlobal.txt') 
        !call FilePrinting(Self%value_KGlobal,'V','DataResults/.InternalData/Vals_KGlobal.txt') 
        Self%FGlobal_res = Self%FGlobal_ext - Self%FGlobal_int
        !call FilePrinting(Self%FGlobal_res,'V','DataResults/.InternalData/FGlobal_res.txt')
    end subroutine AssemblyNonLinearSystem
    ! 5. Mechanical unbalance vector
    subroutine SolveNonLinearSystem(Self)
        implicit none
        class(FEM_MMNL), intent(inout)                             :: Self
        ! Applying HSL-MA86 Solver
        Self%value_FGlobal = Self%FGlobal_res(Self%FreeD)
        Self%value_UGlobal = SparseSystemMA86Solver(Self%Rows_KGlobal,Self%Cols_KGlobal,Self%value_KGlobal,Self%value_FGlobal)
        !Linking Solution
        allocate(Self%dUGlobal(Self%N*Self%DimAnalysis))
        Self%dUGlobal = 0.0d0
        Self%dUGlobal(Self%FreeD) = Self%value_UGlobal
        Self%UGlobal = Self%UGlobal + Self%dUGlobal
    end subroutine SolveNonLinearSystem
    ! 6. Releasing memory
    subroutine ReleasingMemory(Self)
        implicit none
        class(FEM_MMNL), intent(inout)                             :: Self
        deallocate(Self%FGlobal_int)
        deallocate(Self%dUGlobal)
        deallocate(Self%KLocal_NonLinear)
        deallocate(Self%value_KGlobal)
    end subroutine ReleasingMemory
    ! 7. Processing Results
    subroutine ProcessingResults(Self,Psi,PenalFactor)
        implicit none
        class(FEM_MMNL), intent(inout)                             :: Self
        double precision, dimension(:,:), allocatable, intent(in)     :: Psi
        double precision, intent(in)                                :: PenalFactor
        ! internal variables
        integer                                                     :: el,i,j,k
        double precision                                            :: e,n,z,w1,w2,w3,DetJacobian,FAC
        double precision, dimension(:), allocatable                 :: Strain,StressNL
        double precision, dimension(:,:), allocatable               :: F,Cright,Egl,IMatrix
        double precision, dimension(:,:), allocatable               :: Jacobian,InvJacobian,D,DiffN,DiffNXY,IXY
        double precision, dimension(:,:), allocatable               :: ElementCoordinates,ElementDisplacement
        if (Self%DimAnalysis.eq.2) then
            allocate(ElementCoordinates(Self%Npe,2));                          ElementCoordinates = 0.0d0
            allocate(ElementDisplacement(Self%Npe,2));                        ElementDisplacement = 0.0d0
            allocate(Strain(3));                                                           Strain = 0.0d0
            allocate(IMatrix(2,2));                                                       IMatrix = 0.0d0
            IMatrix(1,1) = 1.0d0; IMatrix(2,2) = 1.0d0; 
            ! Results
            ! Allocate ONLY on the first call: ProcessingResults runs once per topology
            ! optimization iteration, and the original single-material code relied on
            ! an explicit ErasePreliminaryResults between iterations to avoid a
            ! double-allocation crash. Guarding here removes that hidden coupling.
            if (.not.allocated(Self%StrainE))       allocate(Self%StrainE(Self%Ne,3))
            if (.not.allocated(Self%StressE))       allocate(Self%StressE(Self%Ne,3))
            if (.not.allocated(Self%StrainEnergyE)) allocate(Self%StrainEnergyE(Self%Ne))
            Self%StrainE = 0.0d0; Self%StressE = 0.0d0; Self%StrainEnergyE = 0.0d0
            do el = 1, Self%Ne, 1
                ! Multi-material interpolated constitutive tensor (Eq. 4)
                call InterpolatedD(Self,Psi(el,:),PenalFactor,D)
                ElementCoordinates = Self%Coordinates(Self%ConnectivityN(el,:),1:2)
                ElementDisplacement = reshape(Self%UGlobal(Self%ConnectivityD(el,:)),[2,Self%Npe])
                do i = 1, Self%QuadGauss, 1
                    e = Self%GaussPoint(i)
                    w1 = Self%GaussWeights(i)
                    do j = 1, Self%QuadGauss, 1
                        n = Self%GaussPoint(j)
                        w2 = Self%GaussWeights(j)
                        call DiffFormFunction(Self,DiffN,e,n)
                        Jacobian = matmul(DiffN,ElementCoordinates)
                        InvJacobian = Inverse(Jacobian)
                        DetJacobian = Determinant(Jacobian)
                        DiffNXY = matmul(InvJacobian,DiffN)
                        IXY = matmul(ElementDisplacement,transpose(DiffNXY))
                        F = IXY + IMatrix
                        Cright = matmul(transpose(F),F)         ! Right Green Tensor
                        Egl = 0.5*(Cright - IMatrix)            ! Green-Lagrange Strain
                        Strain = [Egl(1,1), Egl(2,2), 2*Egl(1,2)]
                        StressNL = matmul(D,Strain)
                        FAC = DetJacobian*w1*w2*(Self%Thickness)
                        ! Strain (per element)
                        Self%StrainE(el,:) = Self%StrainE(el,:) + FAC*Strain
                        ! Stress (per element)
                        Self%StressE(el,:) = Self%StressE(el,:) + FAC*StressNL
                        ! Strain Energy (per element)
                        Self%StrainEnergyE(el) = Self%StrainEnergyE(el) + FAC*0.5d0*dot_product(Strain,StressNL)
                        deallocate(DiffN)
                    end do
                end do
                deallocate(D)
            end do
        elseif(Self%DimAnalysis.eq.3) then
            allocate(ElementCoordinates(Self%Npe,3));                          ElementCoordinates = 0.0d0
            allocate(ElementDisplacement(Self%Npe,3));                        ElementDisplacement = 0.0d0
            allocate(Strain(6));                                                           Strain = 0.0d0
            allocate(IMatrix(3,3));                                                       IMatrix = 0.0d0
            IMatrix(1,1) = 1.0d0; IMatrix(2,2) = 1.0d0; IMatrix(3,3) = 1.0d0; 
            ! Results
            ! Allocate ONLY on the first call: ProcessingResults runs once per topology
            ! optimization iteration, and the original single-material code relied on
            ! an explicit ErasePreliminaryResults between iterations to avoid a
            ! double-allocation crash. Guarding here removes that hidden coupling.
            if (.not.allocated(Self%StrainE))       allocate(Self%StrainE(Self%Ne,6))
            if (.not.allocated(Self%StressE))       allocate(Self%StressE(Self%Ne,6))
            if (.not.allocated(Self%StrainEnergyE)) allocate(Self%StrainEnergyE(Self%Ne))
            Self%StrainE = 0.0d0; Self%StressE = 0.0d0; Self%StrainEnergyE = 0.0d0
            do el = 1, Self%Ne, 1
                ! Multi-material interpolated constitutive tensor (Eq. 4)
                call InterpolatedD(Self,Psi(el,:),PenalFactor,D)
                ElementCoordinates = Self%Coordinates(Self%ConnectivityN(el,:),1:3)
                ElementDisplacement = reshape(Self%UGlobal(Self%ConnectivityD(el,:)),[3,Self%Npe])
                do i = 1, Self%QuadGauss, 1
                    e = Self%GaussPoint(i)
                    w1 = Self%GaussWeights(i)
                    do j = 1, Self%QuadGauss, 1
                        n = Self%GaussPoint(j)
                        w2 = Self%GaussWeights(j)
                        do k = 1, Self%QuadGauss, 1
                            z = Self%GaussPoint(k)
                            w3 = Self%GaussWeights(k)
                            call DiffFormFunction(Self,DiffN,e,n,z)
                            Jacobian = matmul(DiffN,ElementCoordinates)
                            InvJacobian = Inverse(Jacobian)
                            DetJacobian = Determinant(Jacobian)
                            DiffNXY = matmul(InvJacobian,DiffN)
                            IXY = matmul(ElementDisplacement,transpose(DiffNXY))
                            F = IXY + IMatrix
                            Cright = matmul(transpose(F),F)           ! Right Green Tensor
                            Egl = 0.5d0*(Cright - IMatrix)            ! Green-Lagrange Strain
                            Strain = [Egl(1,1), Egl(2,2), Egl(3,3), 2*Egl(1,2), 2*Egl(2,3), 2*Egl(1,3)]
                            StressNL = matmul(D,Strain)
                            FAC = DetJacobian*w1*w2*w3
                            ! Strain (per element)
                            Self%StrainE(el,:) = Self%StrainE(el,:) + FAC*Strain
                            ! Stress (per element)
                            Self%StressE(el,:) = Self%StressE(el,:) + FAC*StressNL
                            ! Strain Energy (per element)
                            Self%StrainEnergyE(el) = Self%StrainEnergyE(el) + FAC*0.5d0*dot_product(Strain,StressNL)
                            deallocate(DiffN)
                        end do
                    end do
                end do
                deallocate(D)
            end do
        end if
        deallocate(ElementCoordinates,ElementDisplacement)
        deallocate(Strain,IMatrix)
        Self%Displacement = transpose(reshape(Self%UGlobal,[Self%DimAnalysis,Self%N]))
    end subroutine ProcessingResults

    ! ----------------------------------------------------------------- !
    !  8. PER-PHASE INTERNAL FORCES (for the adjoint sensitivity)        !
    ! ----------------------------------------------------------------- !
    !  FphaseE(el,:,m) = INT B^T ( D_m : E ) dV  evaluated at the CONVERGED u,
    !  i.e. the internal force element `el` would carry if it were made entirely
    !  of phase m. Slot NMaterial+1 is the void phase.
    !
    !  Because f_int is linear in D and the interpolation is
    !       D_e = D_void + sum_m Psi_m^n ( D_m - D_void )
    !  the design derivative of the element internal force is simply
    !       d(f_int,e)/d(Psi_m) = n * Psi_m^(n-1) * ( FphaseE(el,:,m) - FphaseE(el,:,void) )
    !  with NO extra finite-element integration needed per design variable.
    !
    !  All NMaterial+1 phases are accumulated inside a SINGLE Gauss loop: the
    !  kinematics (F, E, B) are identical for every phase, only D changes.
    subroutine GetPhaseInternalForces(Self)
        implicit none
        class(FEM_MMNL), intent(inout)                              :: Self
        integer                                                     :: el,i,j,k,l,m,ndof,nph
        double precision                                            :: e,n,z,w1,w2,w3,DetJacobian,FAC
        double precision, dimension(:), allocatable                 :: Strain,Stress
        double precision, dimension(:,:), allocatable               :: Jacobian,InvJacobian,DiffN,DiffNXY,IXY
        double precision, dimension(:,:), allocatable               :: F,Cright,Cleft,Egl,Eea,IMatrix,B,BN,BO
        double precision, dimension(:,:), allocatable               :: ElementCoordinates,ElementDisplacement
        ndof = Self%Npe*Self%DimAnalysis
        nph  = Self%NMaterial + 1
        if (.not.allocated(Self%FphaseE)) allocate(Self%FphaseE(Self%Ne,ndof,nph))
        Self%FphaseE = 0.0d0
        if (Self%DimAnalysis.eq.2) then
            allocate(B(3,2*Self%Npe));                    B = 0.0d0
            allocate(BO(3,2*Self%Npe));                  BO = 0.0d0
            allocate(BN(3,2*Self%Npe));                  BN = 0.0d0
            allocate(Strain(3));                     Strain = 0.0d0
            allocate(IMatrix(2,2));                 IMatrix = 0.0d0
            allocate(ElementCoordinates(Self%Npe,2))
            allocate(ElementDisplacement(Self%Npe,2))
            IMatrix(1,1) = 1.0d0; IMatrix(2,2) = 1.0d0
            do el = 1, Self%Ne, 1
                if (Self%Formulation.eq.'TL') then
                    ElementCoordinates = Self%Coordinates(Self%ConnectivityN(el,:),1:2)
                else
                    ElementCoordinates = Self%MatCoordinates(Self%ConnectivityN(el,:),1:2)
                end if
                ElementDisplacement = reshape(Self%UGlobal(Self%ConnectivityD(el,:)),[2,Self%Npe])
                do i = 1, Self%QuadGauss, 1
                    e = Self%GaussPoint(i); w1 = Self%GaussWeights(i)
                    do j = 1, Self%QuadGauss, 1
                        n = Self%GaussPoint(j); w2 = Self%GaussWeights(j)
                        call DiffFormFunction(Self,DiffN,e,n)
                        Jacobian = matmul(DiffN,ElementCoordinates)
                        InvJacobian = Inverse(Jacobian)
                        DetJacobian = Determinant(Jacobian)
                        DiffNXY = matmul(InvJacobian,DiffN)
                        IXY = matmul(ElementDisplacement,transpose(DiffNXY))
                        F = IXY + IMatrix
                        do k = 1, Self%Npe, 1
                            BO(1,2*k-1) = DiffNXY(1,k)
                            BO(2,2*k)   = DiffNXY(2,k)
                            BO(3,2*k-1) = DiffNXY(2,k)
                            BO(3,2*k)   = DiffNXY(1,k)
                            BN(1,2*k-1) = F(1,1)*DiffNXY(1,k)
                            BN(1,2*k)   = F(2,1)*DiffNXY(1,k)
                            BN(2,2*k-1) = F(1,2)*DiffNXY(2,k)
                            BN(2,2*k)   = F(2,2)*DiffNXY(2,k)
                            BN(3,2*k-1) = F(1,2)*DiffNXY(1,k) + F(1,1)*DiffNXY(2,k)
                            BN(3,2*k)   = F(2,2)*DiffNXY(1,k) + F(2,1)*DiffNXY(2,k)
                        end do
                        if (Self%Formulation.eq.'TL') then
                            B = BN          ! TL: BN alone is the complete B_L
                            Cright = matmul(transpose(F),F)
                            Egl = 0.5d0*(Cright - IMatrix)                  ! Green-Lagrange
                            Strain = [Egl(1,1), Egl(2,2), 2.0d0*Egl(1,2)]
                        else
                            B = BO
                            Cleft = matmul(F,transpose(F))
                            Cleft = Inverse(Cleft)
                            Eea = 0.5d0*(IMatrix - Cleft)                   ! Euler-Almansi
                            Strain = [Eea(1,1), Eea(2,2), 2.0d0*Eea(1,2)]
                        end if
                        FAC = DetJacobian*w1*w2*(Self%Thickness)
                        do m = 1, nph, 1
                            Stress = matmul(Self%DMaterial(m,:,:),Strain)
                            Self%FphaseE(el,:,m) = Self%FphaseE(el,:,m) + FAC*matmul(transpose(B),Stress)
                        end do
                        deallocate(DiffN)
                    end do
                end do
            end do
        elseif (Self%DimAnalysis.eq.3) then
            allocate(B(6,3*Self%Npe));                    B = 0.0d0
            allocate(BO(6,3*Self%Npe));                  BO = 0.0d0
            allocate(BN(6,3*Self%Npe));                  BN = 0.0d0
            allocate(Strain(6));                     Strain = 0.0d0
            allocate(IMatrix(3,3));                 IMatrix = 0.0d0
            allocate(ElementCoordinates(Self%Npe,3))
            allocate(ElementDisplacement(Self%Npe,3))
            IMatrix(1,1) = 1.0d0; IMatrix(2,2) = 1.0d0; IMatrix(3,3) = 1.0d0
            do el = 1, Self%Ne, 1
                if (Self%Formulation.eq.'TL') then
                    ElementCoordinates = Self%Coordinates(Self%ConnectivityN(el,:),1:3)
                else
                    ElementCoordinates = Self%MatCoordinates(Self%ConnectivityN(el,:),1:3)
                end if
                ElementDisplacement = reshape(Self%UGlobal(Self%ConnectivityD(el,:)),[3,Self%Npe])
                do i = 1, Self%QuadGauss, 1
                    e = Self%GaussPoint(i); w1 = Self%GaussWeights(i)
                    do j = 1, Self%QuadGauss, 1
                        n = Self%GaussPoint(j); w2 = Self%GaussWeights(j)
                        do k = 1, Self%QuadGauss, 1
                            z = Self%GaussPoint(k); w3 = Self%GaussWeights(k)
                            call DiffFormFunction(Self,DiffN,e,n,z)
                            Jacobian = matmul(DiffN,ElementCoordinates)
                            InvJacobian = Inverse(Jacobian)
                            DetJacobian = Determinant(Jacobian)
                            DiffNXY = matmul(InvJacobian,DiffN)
                            IXY = matmul(ElementDisplacement,transpose(DiffNXY))
                            F = IXY + IMatrix
                            do l = 1, Self%Npe, 1
                                BO(1,3*l-2) = DiffNXY(1,l)
                                BO(2,3*l-1) = DiffNXY(2,l)
                                BO(3,3*l)   = DiffNXY(3,l)
                                BO(4,3*l-2) = DiffNXY(2,l)
                                BO(4,3*l-1) = DiffNXY(1,l)
                                BO(5,3*l-1) = DiffNXY(3,l)
                                BO(5,3*l)   = DiffNXY(2,l)
                                BO(6,3*l-2) = DiffNXY(3,l)
                                BO(6,3*l)   = DiffNXY(1,l)
                                BN(1,3*l-2) = F(1,1)*DiffNXY(1,l)
                                BN(1,3*l-1) = F(2,1)*DiffNXY(1,l)
                                BN(1,3*l)   = F(3,1)*DiffNXY(1,l)
                                BN(2,3*l-2) = F(1,2)*DiffNXY(2,l)
                                BN(2,3*l-1) = F(2,2)*DiffNXY(2,l)
                                BN(2,3*l)   = F(3,2)*DiffNXY(2,l)
                                BN(3,3*l-2) = F(1,3)*DiffNXY(3,l)
                                BN(3,3*l-1) = F(2,3)*DiffNXY(3,l)
                                BN(3,3*l)   = F(3,3)*DiffNXY(3,l)
                                BN(4,3*l-2) = F(1,2)*DiffNXY(1,l) + F(1,1)*DiffNXY(2,l)
                                BN(4,3*l-1) = F(2,2)*DiffNXY(1,l) + F(2,1)*DiffNXY(2,l)
                                BN(4,3*l)   = F(3,2)*DiffNXY(1,l) + F(3,1)*DiffNXY(2,l)
                                BN(5,3*l-2) = F(1,3)*DiffNXY(2,l) + F(1,2)*DiffNXY(3,l)
                                BN(5,3*l-1) = F(2,3)*DiffNXY(2,l) + F(2,2)*DiffNXY(3,l)
                                BN(5,3*l)   = F(3,3)*DiffNXY(2,l) + F(3,2)*DiffNXY(3,l)
                                BN(6,3*l-2) = F(1,3)*DiffNXY(1,l) + F(1,1)*DiffNXY(3,l)
                                BN(6,3*l-1) = F(2,3)*DiffNXY(1,l) + F(2,1)*DiffNXY(3,l)
                                BN(6,3*l)   = F(3,3)*DiffNXY(1,l) + F(3,1)*DiffNXY(3,l)
                            end do
                            if (Self%Formulation.eq.'TL') then
                                B = BN      ! TL: BN alone is the complete B_L
                                Cright = matmul(transpose(F),F)
                                Egl = 0.5d0*(Cright - IMatrix)
                                Strain = [Egl(1,1),Egl(2,2),Egl(3,3),2.0d0*Egl(1,2),2.0d0*Egl(2,3),2.0d0*Egl(1,3)]
                            else
                                B = BO
                                Cleft = matmul(F,transpose(F))
                                Cleft = Inverse(Cleft)
                                Eea = 0.5d0*(IMatrix - Cleft)
                                Strain = [Eea(1,1),Eea(2,2),Eea(3,3),2.0d0*Eea(1,2),2.0d0*Eea(2,3),2.0d0*Eea(1,3)]
                            end if
                            FAC = DetJacobian*w1*w2*w3
                            do m = 1, nph, 1
                                Stress = matmul(Self%DMaterial(m,:,:),Strain)
                                Self%FphaseE(el,:,m) = Self%FphaseE(el,:,m) + FAC*matmul(transpose(B),Stress)
                            end do
                            deallocate(DiffN)
                        end do
                    end do
                end do
            end do
        end if
        deallocate(B,BO,BN,Strain,IMatrix,ElementCoordinates,ElementDisplacement)
    end subroutine GetPhaseInternalForces

    ! ----------------------------------------------------------------- !
    !  9. ADJOINT SYSTEM     K_T * lambda = f_ext                        !
    ! ----------------------------------------------------------------- !
    !  Solved ONCE per topology-optimization iteration, at the converged
    !  displacement field. The tangent is re-assembled here at the exact converged
    !  u rather than reusing the last Newton assembly (which was built at
    !  u_final - du); the cost is one assembly + one factorization, i.e. roughly
    !  one extra Newton iteration.
    !
    !  Sanity check: in linear elasticity f_int = K u and K_T = K, so lambda = u
    !  and the classical self-adjoint result -u^T (dK/dx) u is recovered exactly.
    subroutine SolveAdjointSystem(Self,Psi,PenalFactor)
        implicit none
        class(FEM_MMNL), intent(inout)                              :: Self
        double precision, dimension(:,:), allocatable, intent(in)   :: Psi
        double precision, intent(in)                                :: PenalFactor
        double precision, dimension(:), allocatable                 :: RHS,SOL
        ! tangent stiffness at the converged state
        call GetLocalTangentK(Self,Psi,PenalFactor)
        call GetLocalInternalF(Self,Psi,PenalFactor)
        call AssemblyNonLinearSystem(Self)
        ! right-hand side: the FULL external load (not the incremental one)
        allocate(RHS(size(Self%FreeD)))
        RHS = Self%FGlobal_ext_total(Self%FreeD)
        SOL = SparseSystemMA86Solver(Self%Rows_KGlobal,Self%Cols_KGlobal,Self%value_KGlobal,RHS)
        if (.not.allocated(Self%Lambda)) allocate(Self%Lambda(Self%N*Self%DimAnalysis))
        Self%Lambda = 0.0d0
        Self%Lambda(Self%FreeD) = SOL
        deallocate(RHS,SOL)
        deallocate(Self%FGlobal_int,Self%KLocal_NonLinear,Self%value_KGlobal)
    end subroutine SolveAdjointSystem

    ! 10. modified incremental newton raphson
    subroutine NewtonRaphson(Self,Psi,PenalFactor)
        implicit none
        ! T.O. Variables
        double precision, dimension(:,:), allocatable, intent(in)     :: Psi
        double precision, intent(in)                                :: PenalFactor
        ! FEA Variables
        class(FEM_MMNL), intent(inout)                             :: Self
        integer                                                     :: LoadInc
        integer                                                     :: ite
        double precision                                            :: Residual
        call Initialize(Self)
        ! -------------------- NEWTON RAPHSON -------------------- !
        do LoadInc = 1, Self%Incremental, 1
            write(unit=*, fmt=*) 'Load incremental No.', LoadInc
            write(unit=*, fmt=*) 'Star Newton-Raphson method.'
            ite = 1
            do
                ! Getting all local tangent stiffness matrices
                call GetLocalTangentK(Self,Psi,PenalFactor)
                ! Getting all local internal forces
                call GetLocalInternalF(Self,Psi,PenalFactor)
                ! Assembly non-linear system (getting the GlobalTK and the residual vector)
                call AssemblyNonLinearSystem(Self)
                ! Solve Non-linear system
                call SolveNonLinearSystem(Self)
                ! Update the convected mesh. Only the Updated Lagrangian formulation
                ! integrates on it; in TL the reference configuration is used
                ! throughout, so this is kept purely for post-processing.
                Self%MatCoordinates = Self%MatCoordinates + transpose(reshape(Self%dUGlobal,[Self%DimAnalysis,Self%N]))
                ! Releasing memory (FGlobal,KLocal,duGlobal,etc..)
                call ReleasingMemory(Self)
                ! Output
                Residual = maxval(Abs(Self%FGlobal_res(Self%FreeD)))
                write(unit=*, fmt=*) 'Ite,',ite, 'Res. Value',Residual
                if (ite.gt.Self%MaxFEAIteration) then; write(unit=*, fmt=*) 'Max FEM iterations reached'; exit; end if
                if (Residual.lt.Self%Tolerance)  then; write(unit=*, fmt=*) 'Precision achieved'        ; exit; end if   
                ite = ite + 1
            end do
            ! Update load (the last increment overshoots by one step, which is why
            ! the adjoint uses FGlobal_ext_total and not FGlobal_ext)
            Self%FGlobal_ext = Self%FGlobal_ext + Self%FGlobal_ext_inc
        end do
        call ProcessingResults(Self,Psi,PenalFactor)
        ! ---- everything the sensitivity analysis needs, at the converged state ----
        call GetPhaseInternalForces(Self)
        call SolveAdjointSystem(Self,Psi,PenalFactor)
    end subroutine NewtonRaphson
end module FEA_MMNL_Module