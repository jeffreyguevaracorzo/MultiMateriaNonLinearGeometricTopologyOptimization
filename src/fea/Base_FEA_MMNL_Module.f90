module Base_FEA_MMNL_Module
    use Base_Module
    use Solver_MA86Module
    implicit none
    type                                                :: FEM_MMNL_Base
        character(len=20)                               :: AnalysisType     ! PlaneStress(2D), PlainStrain(2D), SolidIso(3D)
        character(len=20)                               :: ElementType      ! tri3 tri6, cuad4, cuad8, tetra4, tetra10, hexa8, hexa20
        ! ---------------- kinematic formulation ----------------
        ! 'TL' = Total Lagrangian (reference configuration, Green-Lagrange strain,
        !        2nd Piola-Kirchhoff stress). Path-independent: the converged state
        !        does not depend on how the load was applied, so Incremental > 1 can
        !        be used purely as a Newton-Raphson robustness aid.
        ! 'UL' = Updated Lagrangian (current configuration, Euler-Almansi strain).
        !        Path-dependent by construction; requires Incremental > 1.
        ! NOTE: this used to be inferred from Incremental (=1 -> TL, >1 -> UL). It is
        !       now an INDEPENDENT switch, so that TL can be run with several load
        !       increments -- which is what makes Newton-Raphson converge in the
        !       low-density regions of a topology optimization run.
        character(len=2)                                :: Formulation = 'TL'
        integer                                         :: Incremental      ! Incremental Load Steps
        ! ---------------- multi-material data ----------------
        ! Each candidate material has its own Young's modulus AND its own Poisson's
        ! ratio; the void phase is stored in the last slot (NMaterial+1) of DMaterial.
        integer                                         :: NMaterial = 1
        double precision, dimension(:), allocatable     :: EMaterial        ! E_i        (NMaterial)
        double precision, dimension(:), allocatable     :: PoissonMaterial  ! nu_i       (NMaterial)
        double precision                                :: EVoid = 1.0d-9   ! E_void (very small)
        double precision                                :: PoissonVoid      ! nu_void (defaults to nu_1)
        ! DMaterial(NMaterial+1,ncomp,ncomp): constitutive tensor of every phase,
        ! built ONCE by BuildMaterialDTensors. ncomp = 3 (2D) or 6 (3D).
        ! Unlike the linear multi-material code -- which precomputes a full element
        ! stiffness base K0Material(Ne,NMaterial+1,ndof,ndof) -- here only the small
        ! constitutive tensors are stored, because in a geometrically nonlinear
        ! problem K_T depends on u and cannot be precomputed anyway. This is both
        ! exact and dramatically cheaper in memory (288 B/element instead of
        ! ~18 kB/element for a 3D hexa8 with 3 materials).
        double precision, dimension(:,:,:), allocatable :: DMaterial
        integer                                         :: MaxFEAIteration     ! Maximum iterations
        integer                                         :: DimAnalysis      ! 2D or 3D
        integer                                         :: QuadGauss        ! Number of GaussPoints
        integer                                         :: N                ! Number of nodes
        integer                                         :: Ne               ! Number of elements
        integer                                         :: Npe              ! Number of nodes per element
        integer                                         :: NBc              ! Number of nodes with restrictions
        integer                                         :: Npl              ! Number of point loads
        integer                                         :: Ndl              ! Number of distributed loads
        integer, Dimension(:), allocatable              :: FreeD            ! Free degrees of freedom (with filter)
        integer, Dimension(:), allocatable              :: FixedD           ! Fixed degrees of freedom (with filter)
        integer, Dimension(:), allocatable              :: BaseFreeD        ! Free degrees of freedom (without filter)
        integer, Dimension(:), allocatable              :: BaseFixedD       ! Fixed degrees of freedom (without filter)
        integer, dimension(:,:), allocatable            :: ConnectivityN, ConnectivityD     ! Connectivity
        double precision                                :: YoungModulus, PoissonModulus, Thickness     ! Material prop.
        double precision                                :: Tolerance        ! Convergence Tolerance
        double precision, dimension(:), allocatable     :: GaussPoint       ! Gauss Approximation points
        double precision, dimension(:), allocatable     :: GaussWeights     ! Gauss Approximation weights
        double precision, dimension(:), allocatable     :: FGlobal_PL       ! Global vector of point loads
        double precision, dimension(:), allocatable     :: FGlobal_DL       ! Global vector of distributed loads
        double precision, dimension(:,:), allocatable   :: Coordinates      ! Coordinates of nodes(original)
        double precision, dimension(:,:), allocatable   :: MatCoordinates   ! Material Coordinates of nodes(deformed)
        ! System of equations [K]{u}={f} in sparse format
        integer, dimension(:,:,:), allocatable          :: Location_KGlobal 
        integer, dimension(:), allocatable              :: index_i_KGlobal  ! Indx-rows of global stiffness matrix (general)
        integer, dimension(:), allocatable              :: index_j_KGlobal  ! Indx-cols of global stiffness matrix (genereal)
    contains
        procedure                                       :: SetAnalysisType
        procedure                                       :: SetElementType
        procedure                                       :: SetThickness
        procedure                                       :: SetYoungModulus
        procedure                                       :: SetPoissonModulus
        procedure                                       :: SetGaussAprox
        procedure                                       :: SetMaxFEAIteration
        procedure                                       :: SetConvergenceTolerance
        procedure                                       :: ReadFiles
        procedure                                       :: SetLoadIncremental
        procedure                                       :: SetFormulation
        procedure                                       :: SetNMaterial
        procedure                                       :: SetMaterialProperties
        procedure                                       :: SetPoissonModulusMaterial
        procedure                                       :: SetVoidModulus
        procedure                                       :: BuildMaterialDTensors
        procedure                                       :: PreAssemblyRoutine
    end type FEM_MMNL_Base
    contains
    ! --------------- ADDITIONAL BASE FUNCTIONS AND SUBROUTINES ---------------
    ! 1. Reading files routines
    Subroutine ReadingfileInteger(Path,Nrow,Ncol,Matrix)
        implicit none
        character(len=*), intent(in)                         :: Path
        integer                                              :: ios,iounit,i,j
        integer, intent(inout)                               :: Nrow
        integer, intent(in)                                  :: Ncol
        integer, dimension(:,:), allocatable, intent(inout)  :: Matrix
        open(unit=iounit, file=Path, iostat=ios, status="old", action="read")
            if ( ios /= 0 ) stop "Error opening file name"
            read(unit=iounit, fmt=*) !title
            read(unit=iounit, fmt=*) Nrow ; allocate(Matrix(Nrow,Ncol))
            read(unit=iounit, fmt=*) !references
            do i = 1, Nrow, 1
                read(unit=iounit, fmt=*) (Matrix(i,j), j = 1, Ncol, 1)
            end do
        close(iounit)
    end subroutine ReadingfileInteger
    Subroutine ReadingfileDP(Path,Nrow,Ncol,Matrix)
        implicit none
        character(len=*), intent(in)                                  :: Path
        integer                                                       :: ios,iounit,i,j
        integer, intent(inout)                                        :: Nrow
        integer, intent(in)                                           :: Ncol
        double precision, dimension(:,:), allocatable, intent(inout)  :: Matrix
        open(unit=iounit, file=Path, iostat=ios, status="old", action="read")
            if ( ios /= 0 ) stop "Error opening file name"
            read(unit=iounit, fmt=*) !title
            read(unit=iounit, fmt=*) Nrow ; allocate(Matrix(Nrow,Ncol))
            read(unit=iounit, fmt=*) !references
            do i = 1, Nrow, 1
                read(unit=iounit, fmt=*) (Matrix(i,j), j = 1, Ncol, 1)
            end do
        close(iounit)
    end subroutine ReadingfileDP
    ! 2. Area
    function Area(Coordinates,Type) result(AreaAprox)
        character(len=*), intent(in)                                :: Type
        double precision                                            :: AreaAprox
        double precision, dimension(:), allocatable                 :: vector1, vector2, vector3, vector4
        double precision, dimension(:), allocatable                 :: Av1, Av2
        double precision, dimension(:,:), allocatable, intent(in)   :: Coordinates
        allocate(vector1(3),vector2(3),vector3(3),vector4(3))
        allocate(Av1(3), Av2(3))
        if ((Type.eq.'tetra4').or.(Type.eq.'tetra10')) then
            vector1 = Coordinates(2,:)-Coordinates(1,:)
            vector2 = Coordinates(3,:)-Coordinates(1,:)
            Av1 = CrossProduct(vector1,vector2)
            AreaAprox = abs(Norm(Av1))/2
        elseif ((Type.eq.'hexa8').or.(Type.eq.'hexa20')) then
            vector1 = Coordinates(1,:)-Coordinates(2,:)
            vector2 = Coordinates(3,:)-Coordinates(2,:)
            vector3 = Coordinates(1,:)-Coordinates(4,:)
            vector4 = Coordinates(3,:)-Coordinates(4,:)
            Av1 = CrossProduct(vector1,vector2)
            Av2 = CrossProduct(vector3,vector4)
            AreaAprox = abs(Norm(Av1))/2 + abs(Norm(Av2))/2
        end if
    end function Area
    ! 3. Sorting
    subroutine Sort(vector,n)
        implicit none
        integer, intent(in)                 :: n
        integer, allocatable, intent(inout) :: vector(:)
        ! internal variables
        integer                             :: i, j, temp
        integer, allocatable                :: dumb(:)
        ! Ordenado de valores
        do i = 1, n-1
            do j = 1, n-i
                if (vector(j) > vector(j+1)) then
                    temp = vector(j)
                    vector(j) = vector(j+1)
                    vector(j+1) = temp
                endif
            end do
        end do
        ! eliminar 0's y repetidos
        do i = 2, n, 1
            if (Vector(i-1).eq.vector(i)) then 
            vector(i-1) = 0
            end if
        end do
        vector = pack(vector,vector.gt.0)
    end subroutine Sort
    ! ----------------- SUBROUTINES FOR STRUCTURE INFORMATION ----------------- 
    ! 1. Input Analysis Type
    subroutine SetAnalysisType(Self,AnalysisType)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: AnalysisType
        Self%AnalysisType = AnalysisType
        if ((Self%AnalysisType.eq.'PlaneStress').or.(Self%AnalysisType.eq.'PlaneStrain')) then
            Self%DimAnalysis = 2
        elseif (Self%AnalysisType.eq.'SolidIso') then 
            Self%DimAnalysis = 3
        else
            stop "ERROR, Setting AnalysisType"
        end if
    end subroutine SetAnalysisType
    ! 2. Input Element Type
    subroutine SetElementType(Self,ElementType)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: ElementType
        Self%ElementType = ElementType
        if ( Self%ElementType.eq.'tria3' ) then; Self%Npe = 3 ; end if
        if ( Self%ElementType.eq.'tria6' ) then; Self%Npe = 6 ; end if
        if ( Self%ElementType.eq.'quad4' ) then; Self%Npe = 4 ; end if
        if ( Self%ElementType.eq.'quad8' ) then; Self%Npe = 8 ; end if
        if ( Self%ElementType.eq.'tetra4' ) then; Self%Npe = 4 ; end if
        if ( Self%ElementType.eq.'tetra10' ) then; Self%Npe = 10 ; end if
        if ( Self%ElementType.eq.'hexa8' ) then; Self%Npe = 8 ; end if
        if ( Self%ElementType.eq.'hexa20' ) then; Self%Npe = 20 ; end if
    end subroutine SetElementType
    ! 3. Input Thickness (only for 3D cases)
    subroutine SetThickness(Self,Thickness)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        double precision, intent(in)                                :: Thickness
        Self%Thickness = Thickness
    end subroutine SetThickness
    ! 4. Input Young Modulus
    subroutine SetYoungModulus(Self,YoungModulus)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        double precision, intent(in)                                :: YoungModulus
        Self%YoungModulus = YoungModulus
    end subroutine SetYoungModulus
    ! 5. Input Poisson Modulus
    subroutine SetPoissonModulus(Self,PoissonModulus)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        double precision, intent(in)                                :: PoissonModulus
        Self%PoissonModulus = PoissonModulus
    end subroutine SetPoissonModulus
    ! 6. Input Gauss Aproximation
    subroutine SetGaussAprox(Self,Gauss)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self 
        integer, intent(in)                                         :: Gauss
        double precision                                            :: LB, UB, Coef1, Coef2
        if ((Gauss.le.0).or.(Gauss.gt.5)) stop "ERROR, GuassQuadrature greater than 5 or equal to zero"
        Self%QuadGauss = Gauss
        select case (Gauss)
            case (1)
                allocate(Self%GaussPoint(1))
                allocate(Self%GaussWeights(1))
                Self%GaussPoint = [0.00000000d0]
                Self%GaussWeights = [2.00000000d0]
            case (2)
                allocate(Self%GaussPoint(2))
                allocate(Self%GaussWeights(2))
                Self%GaussPoint = [-0.57735026d0,0.57735026d0]
                Self%GaussWeights = [1.00000000d0,1.00000000d0]
            case (3)
                allocate(Self%GaussPoint(3))
                allocate(Self%GaussWeights(3))
                Self%GaussPoint = [-0.77459666d0,0.00000000d0,0.77459666d0]
                Self%GaussWeights = [0.55555555d0,0.88888888d0,0.55555555d0]
            case (4)
                allocate(Self%GaussPoint(4))
                allocate(Self%GaussWeights(4))
                Self%GaussPoint = [-0.86113631d0,-0.33998104d0,0.33998104d0,0.86113631d0]
                Self%GaussWeights = [0.34785484d0,0.65214515d0,0.65214515d0,0.34785484d0]
            case (5)
                allocate(Self%GaussPoint(5))
                allocate(Self%GaussWeights(5))
                Self%GaussPoint = [-0.90617984d0,-0.53846931d0,0.00000000d0,0.53846931d0,0.90617984d0]
                Self%GaussWeights = [0.23692688d0,0.47862867d0,0.56888888d0,0.47862867d0,0.23692688d0]
        end select
        LB = 0.0d0;  UB = 1.0d0
        Coef1 = (UB - LB)/2.0d0
        Coef2 = (UB + LB)/2.0d0
        ! changing coordinates and points due to boundary changes
        if ((Self%ElementType.eq.'tria3').or.(Self%ElementType.eq.'tria6').or. &
            (Self%ElementType.eq.'tetra4').or.(Self%ElementType.eq.'tetra10')) then
            Self%GaussPoint = Self%GaussPoint*Coef1 + Coef2
            Self%GaussWeights = Self%GaussWeights*Coef1
        end if
    end subroutine SetGaussAprox

    ! Input Hammer Aproximation
    subroutine SetHammerAprox(Self,Gauss)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self 
        integer, intent(in)                                         :: Gauss
        double precision                                            :: LB, UB, Coef1, Coef2
        if (Gauss.le.0.or.Gauss.gt.5) stop "ERROR, GuassQuadrature (greater than 5 or <= to zero)" 
        Self%QuadGauss = Gauss
        select case (Gauss)
            case (1)
                allocate(Self%GaussPoint(1))
                allocate(Self%GaussWeights(1))
                Self%GaussPoint = [0.00000000d0]
                Self%GaussWeights = [2.00000000d0]
            case (2)
                allocate(Self%GaussPoint(2))
                allocate(Self%GaussWeights(2))
                Self%GaussPoint = [-0.57735026d0,0.57735026d0]
                Self%GaussWeights = [1.00000000d0,1.00000000d0]
            case (3)
                allocate(Self%GaussPoint(3))
                allocate(Self%GaussWeights(3))
                Self%GaussPoint = [-0.77459666d0,0.00000000d0,0.77459666d0]
                Self%GaussWeights = [0.55555555d0,0.88888888d0,0.55555555d0]
            case (4)
                allocate(Self%GaussPoint(4))
                allocate(Self%GaussWeights(4))
                Self%GaussPoint = [-0.86113631d0,-0.33998104d0,0.33998104d0,0.86113631d0]
                Self%GaussWeights = [0.34785484d0,0.65214515d0,0.65214515d0,0.34785484d0]
            case (5)
                allocate(Self%GaussPoint(5))
                allocate(Self%GaussWeights(5))
                Self%GaussPoint = [-0.90617984d0,-0.53846931d0,0.00000000d0,0.53846931d0,0.90617984d0]
                Self%GaussWeights = [0.23692688d0,0.47862867d0,0.56888888d0,0.47862867d0,0.23692688d0]
        end select
    end subroutine SetHammerAprox

    ! 7. Input Coordinates
    subroutine SetCoordinates(Self,Path)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        call ReadingfileDP(Path,Self%N,Self%DimAnalysis,Self%Coordinates)
        write(unit=*, fmt=*) '- Coordinates'
        Self%MatCoordinates = Self%Coordinates
        !call FilePrinting(Self%Coordinates,'DataResults/.InternalData/CheckingCoordinatesLecture.txt')
    end subroutine SetCoordinates
    ! 8. Input Connectivity
    subroutine SetConnectivity(Self,Path)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        integer                                                     :: i,j,k
        call ReadingfileInteger(Path,Self%Ne,Self%Npe,Self%ConnectivityN)
        allocate(Self%ConnectivityD(Self%Ne,Self%Npe*Self%DimAnalysis))
        do i = 1, Self%Ne, 1
            do j = 1, Self%Npe, 1
                do k = Self%DimAnalysis-1, 0, -1
                    Self%ConnectivityD(i,j*Self%DimAnalysis - k) = Self%ConnectivityN(i,j)*Self%DimAnalysis - k
                end do
            end do
        end do
        write(unit=*, fmt=*) '- Connectivity'
        !call FilePrinting(Self%ConnectivityN,'DataResults/.InternalData/CheckingConnectivityNLecture.txt')
        !call FilePrinting(Self%ConnectivityD,'DataResults/.InternalData/CheckingConnectivityDLecture.txt')
    end subroutine SetConnectivity
    ! 9. Input Boundary Conditions
    subroutine SetBondaryConditions(Self,Path)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        integer, dimension(:,:), allocatable                        :: BoundaryC
        integer                                                     :: i,j,k
        call ReadingfileInteger(Path,Self%NBc,(Self%DimAnalysis+1),BoundaryC)
        i = (Self%N)*(Self%DimAnalysis) - count(BoundaryC(:,2:).eq.1)
        j = count(BoundaryC(:,2:).eq.1)
        k = 1
        allocate(Self%BaseFreeD(i),Self%FreeD(i))
        allocate(Self%BaseFixedD(j),Self%FixedD(j))
        Self%BaseFreeD = 0
        Self%BaseFixedD = 0
        ! constrained?
        do i = 1, Self%NBc, 1
            do j = 1, Self%DimAnalysis, 1
                if (BoundaryC(i,j+1).eq.1) then
                    Self%BaseFixedD(k) = BoundaryC(i,1)*Self%DimAnalysis - (Self%DimAnalysis-j)
                    k = k + 1
                else
                    cycle
                end if
            end do
        end do
        ! free!
        k = 1
        do i = 1, Self%N*Self%DimAnalysis, 1
            if (any(Self%BaseFixedD.eq.i)) then
                cycle
            else
                Self%BaseFreeD(k) = i
                k = k + 1
            end if
        end do
        Self%FixedD = Self%BaseFixedD
        Self%FreeD = Self%BaseFreeD 
        write(unit=*, fmt=*) '- Boundary conditions/Constrains'
        !call FilePrinting(Self%BaseFreeD,'V','DataResults/.InternalData/BaseFreeDof.txt')
        !call FilePrinting(Self%BaseFixedD,'V','DataResults/.InternalData/BaseFixedDof.txt')
    end subroutine SetBondaryConditions
    ! 10. Input Point Loads
    subroutine SetPointLoads(Self,Path)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        integer                                                     :: i,j,node,dim
        double precision, dimension(:,:), allocatable               :: PointLoads
        dim = Self%DimAnalysis
        call ReadingfileDP(Path,Self%Npl,(dim+1),PointLoads)
        allocate(Self%FGlobal_PL(Self%N*dim)) ; Self%FGlobal_PL = 0.0d0
        ! assembly load vector
        do i = 1, Self%Npl, 1
            node = int(PointLoads(i,1))
            do j = 1, dim, 1
                Self%FGlobal_PL(node*dim-(dim-j)) = Self%FGlobal_PL(node*dim-(j-dim)) + PointLoads(i,j+1)
            end do
        end do
        write(unit=*, fmt=*) '- Point Loads'
        !call FilePrinting(Self%FGlobal_PL,'V','DataResults/.InternalData/CheckingPointLoadsLecture.txt')
    end subroutine SetPointLoads
    ! 11. Input Distributed Loads
    subroutine SetDistributedLoads(Self,Path)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        integer                                                     :: n,i,j,k,Npf,dim,node
        double precision                                            :: FaceArea, length
        double precision, dimension(:), allocatable                 :: Vector, ForcePerNode
        double precision, dimension(:,:), allocatable               :: DistributedLoads
        double precision, dimension(:,:), allocatable               :: Coordinates
        ! calculation and assembly of load vectors
        dim = Self%DimAnalysis
        allocate(Self%FGlobal_DL(Self%N*dim)); Self%FGlobal_DL = 0.0d0;
        allocate(ForcePerNode(dim));              ForcePerNode = 0.0d0;
        allocate(Vector(dim));                          vector = 0.0d0;
        ! nodes per face
        if (Self%ElementType.eq.'tria3'.or.Self%ElementType.eq.'quad4') then; Npf = 2; end if
        if (Self%ElementType.eq.'tria6'.or.Self%ElementType.eq.'quad8') then; Npf = 3; end if
        if (Self%ElementType.eq.'tetra4')                               then; Npf = 3; n=3; allocate(Coordinates(3,3)); end if
        if (Self%ElementType.eq.'tetra10')                              then; Npf = 6; n=3; allocate(Coordinates(3,3)); end if
        if (Self%ElementType.eq.'hexa8')                                then; Npf = 4; n=4; allocate(Coordinates(4,3)); end if
        if (Self%ElementType.eq.'hexa20')                               then; Npf = 8; n=4; allocate(Coordinates(4,3)); end if
        ! reading file
        call ReadingfileDP(Path,Self%Ndl,(Npf+dim),DistributedLoads)
        select case (dim)
            case (2)
                do i = 1, Self%Ndl, 1
                    Vector = Self%Coordinates(DistributedLoads(i,2),1:dim) - Self%Coordinates(DistributedLoads(i,1),1:dim)
                    length = Norm(Vector)
                    FaceArea = length*Self%Thickness
                    ForcePerNode = (FaceArea/Npf)*(DistributedLoads(i,(Npf+1):))
                    do j = 1, Npf, 1
                        node = int(DistributedLoads(i,j))
                        do k = 1, dim, 1
                            Self%FGlobal_DL(node*dim-(dim-k)) = ForcePerNode(k)
                        end do
                    end do
                end do
            case (3)
                do i = 1, Self%Ndl, 1
                    Coordinates = Self%Coordinates(int(DistributedLoads(i,1:n)),:)
                    FaceArea = Area(Coordinates,Self%ElementType)
                    ForcePerNode = (FaceArea/Npf)*DistributedLoads(i,(Npf+1):)
                    do j = 1, Npf, 1
                        node = int(DistributedLoads(i,j))
                        do k = 1, dim, 1
                            Self%FGlobal_DL(node*dim-(dim-k)) = ForcePerNode(k)
                        end do
                    end do
                end do
        end select
        write(unit=*, fmt=*) '- Distributed Loads'
        !call FilePrinting(Self%FGlobal_PL,'V','DataResults/.InternalData/CheckingDistributedLoadsLecture.txt')
    end subroutine SetDistributedLoads
    ! 12. Read Files(Connectivity,Coordinates,BoundaryConditions,etc...)
    subroutine ReadFiles(Self)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        call SetCoordinates(Self,'input/Coordinates.txt')
        call SetConnectivity(Self,'input/Connectivity.txt')
        call SetBondaryConditions(Self,'input/Constrains.txt')
        call SetPointLoads(Self,'input/PointLoads.txt')
        call SetDistributedLoads(Self,'input/DistributedLoads.txt')
    end subroutine ReadFiles
    ! 13. Input Load Incremental
    subroutine SetLoadIncremental(Self,Incremental)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        integer                                                     :: Incremental
        Self%Incremental = Incremental
    end subroutine SetLoadIncremental
    ! 14. Input Convergence Tolerance
    subroutine SetConvergenceTolerance(Self,Tolerance)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        double precision                                            :: Tolerance
        Self%Tolerance = Tolerance
    end subroutine SetConvergenceTolerance
    ! 15. Input Max Iteration
    subroutine SetMaxFEAIteration(Self,MaxFEAIteration)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        integer                                                     :: MaxFEAIteration
        Self%MaxFEAIteration = MaxFEAIteration
    end subroutine SetMaxFEAIteration
    ! 16. Kinematic formulation: 'TL' (Total Lagrangian) or 'UL' (Updated Lagrangian).
    !     Independent of SetLoadIncremental -- see the note in the type definition.
    subroutine SetFormulation(Self,Formulation)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        character(len=*), intent(in)                                    :: Formulation
        if ((Formulation.ne.'TL').and.(Formulation.ne.'UL')) then
            stop "ERROR, in SetFormulation: only 'TL' or 'UL' are accepted"
        end if
        Self%Formulation = Formulation
    end subroutine SetFormulation

    ! ----------------------------------------------------------------- !
    !                    MULTI-MATERIAL PROPERTIES                       !
    ! ----------------------------------------------------------------- !
    ! 17. Number of candidate materials
    subroutine SetNMaterial(Self,NMaterial)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        integer, intent(in)                                             :: NMaterial
        Self%NMaterial = NMaterial
    end subroutine SetNMaterial
    ! 18. Young's modulus of each material, e.g. SetMaterialProperties(M,[1000.0d0,5000.0d0])
    subroutine SetMaterialProperties(Self,EMaterial)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        double precision, dimension(:), intent(in)                      :: EMaterial
        if (allocated(Self%EMaterial)) deallocate(Self%EMaterial)
        allocate(Self%EMaterial(size(EMaterial)))
        Self%EMaterial = EMaterial
        ! keep the single-material field consistent, some legacy routines read it
        Self%YoungModulus = EMaterial(1)
    end subroutine SetMaterialProperties
    ! 19. Poisson's ratio of each material, SAME ORDER as SetMaterialProperties
    subroutine SetPoissonModulusMaterial(Self,PoissonMaterial)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        double precision, dimension(:), intent(in)                      :: PoissonMaterial
        if (allocated(Self%PoissonMaterial)) deallocate(Self%PoissonMaterial)
        allocate(Self%PoissonMaterial(size(PoissonMaterial)))
        Self%PoissonMaterial = PoissonMaterial
        Self%PoissonVoid     = PoissonMaterial(1)
        Self%PoissonModulus  = PoissonMaterial(1)
    end subroutine SetPoissonModulusMaterial
    ! 20. Young's modulus of the void phase (E_void in Eq. 4 of Zheng et al. 2024)
    subroutine SetVoidModulus(Self,EVoid)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        double precision, intent(in)                                    :: EVoid
        Self%EVoid = EVoid
    end subroutine SetVoidModulus

    ! 21. Constitutive tensor of every phase, built ONCE before the optimization loop.
    !     DMaterial(i,:,:) for i = 1..NMaterial are the real materials, and
    !     DMaterial(NMaterial+1,:,:) is the void phase. Reuses ElasticityTensor by
    !     temporarily swapping (YoungModulus, PoissonModulus) -- the same trick the
    !     linear multi-material code uses in BuildMaterialK0Basis.
    subroutine BuildMaterialDTensors(Self)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        integer                                                         :: mat,ncomp
        double precision                                                :: EBackup,VBackup
        double precision, dimension(:,:), allocatable                   :: D
        if (.not.allocated(Self%EMaterial))       stop "ERROR: call SetMaterialProperties before BuildMaterialDTensors"
        if (.not.allocated(Self%PoissonMaterial)) stop "ERROR: call SetPoissonModulusMaterial before BuildMaterialDTensors"
        if (size(Self%EMaterial).ne.Self%NMaterial) stop "ERROR: size(EMaterial) /= NMaterial"
        if (size(Self%PoissonMaterial).ne.Self%NMaterial) stop "ERROR: size(PoissonMaterial) /= NMaterial"
        if (Self%DimAnalysis.eq.2) then; ncomp = 3; else; ncomp = 6; end if
        if (allocated(Self%DMaterial)) deallocate(Self%DMaterial)
        allocate(Self%DMaterial(Self%NMaterial+1,ncomp,ncomp))
        EBackup = Self%YoungModulus
        VBackup = Self%PoissonModulus
        do mat = 1, Self%NMaterial+1, 1
            if (mat.le.Self%NMaterial) then
                Self%YoungModulus   = Self%EMaterial(mat)
                Self%PoissonModulus = Self%PoissonMaterial(mat)
            else
                Self%YoungModulus   = Self%EVoid
                Self%PoissonModulus = Self%PoissonVoid
            end if
            call ElasticityTensor(Self,D)
            Self%DMaterial(mat,:,:) = D
            deallocate(D)
        end do
        Self%YoungModulus   = EBackup
        Self%PoissonModulus = VBackup
    end subroutine BuildMaterialDTensors

    ! 22. ---------- THE CORE OF THE MULTI-MATERIAL / NONLINEAR FUSION ----------
    !     Interpolated constitutive tensor of one element:
    !
    !         D_e = D_void + sum_i  Psi_i(e)^n * ( D_i - D_void )
    !
    !     with Psi_i the mapping-based interpolation weights (Eq. 4 of Zheng et al.
    !     2024) and n the penalization exponent.
    !
    !     WHY THIS IS ENOUGH. In the Total Lagrangian / St. Venant-Kirchhoff setting
    !     used by this code both the tangent stiffness and the internal force vector
    !     are LINEAR in D:
    !            K_T   = INT B^T D B dV  +  INT B_G^T H(S) B_G dV ,  S = D:E_GL
    !            f_int = INT B^T S dV    =  INT B^T D E_GL dV
    !     so interpolating D is exactly equivalent to interpolating the strain energy
    !         Psi_e(F) = sum_i w_i * Psi_i(F)
    !     which is what a general hyperelastic multi-material formulation would
    !     require. Here it comes for free because the StVK energy (1/2 E:D:E) is
    !     itself linear in D. Consequence: the whole Newton-Raphson machinery,
    !     the sparse assembly and the MA86 solver are reused UNCHANGED.
    !
    !     PsiRow: Psi(el,1:NMaterial) for the element being integrated.
    subroutine InterpolatedD(Self,PsiRow,PenalFactor,D)
        implicit none
        class(FEM_MMNL_Base), intent(in)                                :: Self
        double precision, dimension(:), intent(in)                      :: PsiRow
        double precision, intent(in)                                    :: PenalFactor
        double precision, dimension(:,:), allocatable, intent(inout)    :: D
        integer                                                         :: mat,ncomp
        ncomp = size(Self%DMaterial,2)
        if (.not.allocated(D)) allocate(D(ncomp,ncomp))
        D = Self%DMaterial(Self%NMaterial+1,:,:)                        ! void base
        do mat = 1, Self%NMaterial, 1
            D = D + (PsiRow(mat)**PenalFactor) * &
                    (Self%DMaterial(mat,:,:) - Self%DMaterial(Self%NMaterial+1,:,:))
        end do
    end subroutine InterpolatedD
    ! ------------ FINITE ELEMENT ANALYSISS FUNCTIONS AND SUBROUTINES ------------
    ! 1. Derivative of the shape functions
    subroutine DiffFormFunction(Self,DiffFunction,e,n,z)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                             :: Self
        double precision, intent(inout)                             :: e, n
        double precision, intent(inout), optional                   :: z
        double precision, dimension(:,:), allocatable, intent(out)  :: DiffFunction
        ! dN1/de     dN2/de     dN3/de      ....     dNn/de
        ! dN1/dn     dN2/dn     dN3/dn      ....     dNn/dn
        ! dN1/dz     dN2/dz     dN3/dz      ....     dNn/dz (caso 3D)
        if (Self%ElementType.eq.'tria3') then
            allocate(DiffFunction(2,3))
            !  line 1
            DiffFunction(1,1) = 1.0d0
            DiffFunction(1,2) = 0.0d0
            DiffFunction(1,3) = - 1.0d0
            !  line 2
            DiffFunction(2,1) = 0.0d0
            DiffFunction(2,2) = 1.0d0
            DiffFunction(2,3) = - 1.0d0
        elseif (Self%ElementType.eq.'tria6') then
            allocate(DiffFunction(2,6))
            !  line 1
            DiffFunction(1,1) = 4.0d0*e + 4.0d0*n - 3.0d0 
            DiffFunction(1,2) = 4.0d0*e - 1.0d0 
            DiffFunction(1,3) = 0.0d0 
            DiffFunction(1,4) = 4.0d0 - 4.0d0*n - 8.0d0*e 
            DiffFunction(1,5) = 4.0d0*n 
            DiffFunction(1,6) = - 4.0d0*n 
            !  line 2
            DiffFunction(2,1) = 4.0d0*e + 4.0d0*n - 3.0d0 
            DiffFunction(2,2) = 0.0d0 
            DiffFunction(2,3) = 4.0d0*n - 1.0d0 
            DiffFunction(2,4) = - 4.0d0*e 
            DiffFunction(2,5) = 4.0d0*e 
            DiffFunction(2,6) = 4.0d0 - 8.0d0*n - 4.0d0*e 
        elseif (Self%ElementType.eq.'quad4') then
            allocate(DiffFunction(2,4))
            !  line 1
            DiffFunction(1,1) = n/4.0d0 - 1.0d0/4.0d0 
            DiffFunction(1,2) = 1.0d0/4.0d0 - n/4.0d0 
            DiffFunction(1,3) = n/4.0d0 + 1.0d0/4.0d0 
            DiffFunction(1,4) = - n/4.0d0 - 1.0d0/4.0d0 
            !  line 2
            DiffFunction(2,1) = e/4.0d0 - 1.0d0/4.0d0 
            DiffFunction(2,2) = - e/4.0d0 - 1.0d0/4.0d0 
            DiffFunction(2,3) = e/4.0d0 + 1.0d0/4.0d0 
            DiffFunction(2,4) = 1.0d0/4.0d0 - e/4.0d0 
        elseif (Self%ElementType.eq.'quad8') then
            allocate(DiffFunction(2,8))
            !  line 1
            DiffFunction(1,1) = - (e/4.0d0 - 1.0d0/4.0d0)*(n - 1.0d0) - ((n - 1.0d0)*(e + n + 1.0d0))/4.0d0 
            DiffFunction(1,2) = ((n - 1.0d0)*(n - e + 1.0d0))/4.0d0 - (e/4.0d0 + 1.0d0/4.0d0)*(n - 1.0d0) 
            DiffFunction(1,3) = (e/4.0d0 + 1.0d0/4.0d0)*(n + 1.0d0) + ((n + 1.0d0)*(e + n - 1.0d0))/4.0d0 
            DiffFunction(1,4) = (e/4.0d0 - 1.0d0/4.0d0)*(n + 1.0d0) + ((n + 1.0d0)*(e - n + 1.0d0))/4.0d0 
            DiffFunction(1,5) = e*(n - 1.0d0) 
            DiffFunction(1,6) = 1.0d0/2.0d0 - n**2.0d0/2.0d0 
            DiffFunction(1,7) = - e*(n + 1.0d0) 
            DiffFunction(1,8) = n**2.0d0/2.0d0 - 1.0d0/2.0d0 
            !  line 2
            DiffFunction(2,1) = - (e/4.0d0 - 1.0d0/4.0d0)*(n - 1.0d0) - (e/4.0d0 - 1.0d0/4.0d0)*(e + n + 1.0d0) 
            DiffFunction(2,2) = (e/4.0d0 + 1.0d0/4.0d0)*(n - e + 1.0d0) + (e/4.0d0 + 1.0d0/4.0d0)*(n - 1.0d0) 
            DiffFunction(2,3) = (e/4.0d0 + 1.0d0/4.0d0)*(n + 1.0d0) + (e/4.0d0 + 1.0d0/4.0d0)*(e + n - 1.0d0) 
            DiffFunction(2,4) = (e/4.0d0 - 1.0d0/4.0d0)*(e - n + 1.0d0) - (e/4.0d0 - 1.0d0/4.0d0)*(n + 1.0d0) 
            DiffFunction(2,5) = e**2.0d0/2.0d0 - 1.0d0/2.0d0 
            DiffFunction(2,6) = - 2.0d0*n*(e/2.0d0 + 1.0d0/2.0d0) 
            DiffFunction(2,7) = 1.0d0/2.0d0 - e**2.0d0/2.0d0 
            DiffFunction(2,8) = 2.0d0*n*(e/2.0d0 - 1.0d0/2.0d0) 
        elseif (Self%ElementType.eq.'tetra4') then
            allocate(DiffFunction(3,4))
            !  line 1
            DiffFunction(1,1) = - 1.0d0 
            DiffFunction(1,2) = 1.0d0
            DiffFunction(1,3) = 0.0d0
            DiffFunction(1,4) = 0.0d0
            !  line 2
            DiffFunction(2,1) = - 1.0d0
            DiffFunction(2,2) = 0.0d0
            DiffFunction(2,3) = 1.0d0
            DiffFunction(2,4) = 0.0d0
            !  line 3
            DiffFunction(3,1) = - 1.0d0
            DiffFunction(3,2) = 0.0d0
            DiffFunction(3,3) = 0.0d0
            DiffFunction(3,4) = 1.0d0
        elseif (Self%ElementType.eq.'tetra10') then
            allocate(DiffFunction(3,10))
            !  line 1
            DiffFunction(1,1) = 4.0d0*e + 4.0d0*n + 4.0d0*z - 3.0d0 
            DiffFunction(1,2) = 4.0d0*e - 1.0d0 
            DiffFunction(1,3) = 0.0d0 
            DiffFunction(1,4) = 0.0d0 
            DiffFunction(1,5) = 4.0d0 - 4.0d0*n - 4.0d0*z - 8.0d0*e 
            DiffFunction(1,6) = 4.0d0*n 
            DiffFunction(1,7) = - 4.0d0*n 
            DiffFunction(1,8) = 4.0d0*z 
            DiffFunction(1,9) = 0 
            DiffFunction(1,10) = - 4.0d0*z 
            !  line 2
            DiffFunction(2,1) = 4.0d0*e + 4.0d0*n + 4.0d0*z - 3.0d0 
            DiffFunction(2,2) = 0.0d0 
            DiffFunction(2,3) = 4.0d0*n - 1.0d0 
            DiffFunction(2,4) = 0.0d0 
            DiffFunction(2,5) = - 4.0d0*e 
            DiffFunction(2,6) = 4.0d0*e 
            DiffFunction(2,7) = 4.0d0 - 8.0d0*n - 4.0d0*z - 4.0d0*e 
            DiffFunction(2,8) = 0.0d0 
            DiffFunction(2,9) = 4.0d0*z 
            DiffFunction(2,10) = - 4.0d0*z 
            !  line 3
            DiffFunction(3,1) = 4.0d0*e + 4.0d0*n + 4.0d0*z - 3.0d0 
            DiffFunction(3,2) = 0.0d0 
            DiffFunction(3,3) = 0.0d0 
            DiffFunction(3,4) = 4.0d0*z - 1.0d0 
            DiffFunction(3,5) = - 4.0d0*e 
            DiffFunction(3,6) = 0.0d0 
            DiffFunction(3,7) = - 4.0d0*n 
            DiffFunction(3,8) = 4.0d0*e 
            DiffFunction(3,9) = 4.0d0*n 
            DiffFunction(3,10) = 4.0d0 - 4.0d0*n - 8.0d0*z - 4.0d0*e 
        elseif (Self%ElementType.eq.'hexa8') then
            allocate(DiffFunction(3,8))
            !  line 1
            DiffFunction(1,1) = - ((n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,2) = ((n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,3) = - ((n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,4) = ((n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,5) = ((n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,6) = - ((n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,7) = ((n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,8) = - ((n + 1.0d0)*(z + 1.0d0))/8.0d0 
            !  line 2
            DiffFunction(2,1) = - (e/8.0d0 - 1.0d0/8.0d0)*(z - 1.0d0) 
            DiffFunction(2,2) = (e/8.0d0 + 1.0d0/8.0d0)*(z - 1.0d0) 
            DiffFunction(2,3) = - (e/8.0d0 + 1.0d0/8.0d0)*(z - 1.0d0) 
            DiffFunction(2,4) = (e/8.0d0 - 1.0d0/8.0d0)*(z - 1.0d0) 
            DiffFunction(2,5) = (e/8.0d0 - 1.0d0/8.0d0)*(z + 1.0d0) 
            DiffFunction(2,6) = - (e/8.0d0 + 1.0d0/8.0d0)*(z + 1.0d0) 
            DiffFunction(2,7) = (e/8.0d0 + 1.0d0/8.0d0)*(z + 1.0d0) 
            DiffFunction(2,8) = - (e/8.0d0 - 1.0d0/8.0d0)*(z + 1.0d0) 
            !  line 3
            DiffFunction(3,1) = - (e/8.0d0 - 1.0d0/8.0d0)*(n - 1.0d0) 
            DiffFunction(3,2) = (e/8.0d0 + 1.0d0/8.0d0)*(n - 1.0d0) 
            DiffFunction(3,3) = - (e/8.0d0 + 1.0d0/8.0d0)*(n + 1.0d0) 
            DiffFunction(3,4) = (e/8.0d0 - 1.0d0/8.0d0)*(n + 1.0d0) 
            DiffFunction(3,5) = (e/8.0d0 - 1.0d0/8.0d0)*(n - 1.0d0) 
            DiffFunction(3,6) = - (e/8.0d0 + 1.0d0/8.0d0)*(n - 1.0d0) 
            DiffFunction(3,7) = (e/8.0d0 + 1.0d0/8.0d0)*(n + 1.0d0) 
            DiffFunction(3,8) = - (e/8.0d0 - 1.0d0/8.0d0)*(n + 1.0d0) 
        elseif (Self%ElementType.eq.'hexa20') then
            allocate(DiffFunction(3,20))
            !  line 1
            DiffFunction(1,1) = - (e*n*z*(n - 1.0d0)*(z - 1.0d0))/8.0d0 - (n*z*(e - 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,2) = (e*n*z*(n - 1.0d0)*(z - 1.0d0))/8.0d0 + (n*z*(e + 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,3) = - (e*n*z*(n + 1.0d0)*(z - 1.0d0))/8.0d0 - (n*z*(e + 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,4) = (e*n*z*(n + 1.0d0)*(z - 1.0d0))/8.0d0 + (n*z*(e - 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,5) = (e*n*z*(n - 1.0d0)*(z + 1.0d0))/8.0d0 + (n*z*(e - 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,6) = - (e*n*z*(n - 1.0d0)*(z + 1.0d0))/8.0d0 - (n*z*(e + 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,7) = (e*n*z*(n + 1.0d0)*(z + 1.0d0))/8.0d0 + (n*z*(e + 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,8) = - (e*n*z*(n + 1.0d0)*(z + 1.0d0))/8.0d0 - (n*z*(e - 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,9) = - (e*(n - 1.0d0)*(z - 1.0d0))/2.0d0 
            DiffFunction(1,10) = (e*(n + 1.0d0)*(z - 1.0d0))/2.0d0 
            DiffFunction(1,11) = - (e*(n + 1.0d0)*(z + 1.0d0))/2.0d0 
            DiffFunction(1,12) = (e*(n - 1.0d0)*(z + 1.0d0))/2.0d0 
            DiffFunction(1,13) = - ((n**2.0d0 - 1.0d0)*(z - 1.0d0))/4.0d0 
            DiffFunction(1,14) = ((n**2.0d0 - 1.0d0)*(z - 1.0d0))/4.0d0 
            DiffFunction(1,15) = - ((n**2.0d0 - 1.0d0)*(z + 1.0d0))/4.0d0 
            DiffFunction(1,16) = ((n**2.0d0 - 1.0d0)*(z + 1.0d0))/4.0d0 
            DiffFunction(1,17) = - ((z**2.0d0 - 1.0d0)*(n - 1.0d0))/4.0d0 
            DiffFunction(1,18) = ((z**2.0d0 - 1.0d0)*(n - 1.0d0))/4.0d0 
            DiffFunction(1,19) = - ((z**2.0d0 - 1.0d0)*(n + 1.0d0))/4.0d0 
            DiffFunction(1,20) = ((z**2.0d0 - 1.0d0)*(n + 1.0d0))/4.0d0 
            !  line 2
            DiffFunction(2,1) = - (e*n*z*(e - 1.0d0)*(z - 1.0d0))/8.0d0 - (e*z*(e - 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(2,2) = (e*n*z*(e + 1.0d0)*(z - 1.0d0))/8.0d0 + (e*z*(e + 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(2,3) = - (e*n*z*(e + 1.0d0)*(z - 1.0d0))/8.0d0 - (e*z*(e + 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(2,4) = (e*n*z*(e - 1.0d0)*(z - 1.0d0))/8.0d0 + (e*z*(e - 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(2,5) = (e*n*z*(e - 1.0d0)*(z + 1.0d0))/8.0d0 + (e*z*(e - 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(2,6) = - (e*n*z*(e + 1.0d0)*(z + 1.0d0))/8.0d0 - (e*z*(e + 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(2,7) = (e*n*z*(e + 1.0d0)*(z + 1.0d0))/8.0d0 + (e*z*(e + 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(2,8) = - (e*n*z*(e - 1.0d0)*(z + 1.0d0))/8.0d0 - (e*z*(e - 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(2,9) = - (e**2.0d0/4.0d0 - 1.0d0/4.0d0)*(z - 1.0d0) 
            DiffFunction(2,10) = (e**2.0d0/4.0d0 - 1.0d0/4.0d0)*(z - 1.0d0) 
            DiffFunction(2,11) = - (e**2.0d0/4.0d0 - 1.0d0/4.0d0)*(z + 1.0d0) 
            DiffFunction(2,12) = (e**2.0d0/4.0d0 - 1.0d0/4.0d0)*(z + 1.0d0) 
            DiffFunction(2,13) = - 2.0d0*n*(e/4.0d0 - 1.0d0/4.0d0)*(z - 1.0d0) 
            DiffFunction(2,14) = 2.0d0*n*(e/4.0d0 + 1.0d0/4.0d0)*(z - 1.0d0) 
            DiffFunction(2,15) = - 2.0d0*n*(e/4.0d0 + 1.0d0/4.0d0)*(z + 1.0d0) 
            DiffFunction(2,16) = 2.0d0*n*(e/4.0d0 - 1.0d0/4.0d0)*(z + 1.0d0) 
            DiffFunction(2,17) = - (e/4.0d0 - 1.0d0/4.0d0)*(z**2.0d0 - 1.0d0) 
            DiffFunction(2,18) = (e/4.0d0 + 1.0d0/4.0d0)*(z**2.0d0 - 1.0d0) 
            DiffFunction(2,19) = - (e/4.0d0 + 1.0d0/4.0d0)*(z**2.0d0 - 1.0d0) 
            DiffFunction(2,20) = (e/4.0d0 - 1.0d0/4.0d0)*(z**2.0d0 - 1.0d0) 
            !  line 3
            DiffFunction(3,1) = - (e*n*z*(e - 1.0d0)*(n - 1.0d0))/8.0d0 - (e*n*(e - 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(3,2) = (e*n*z*(e + 1.0d0)*(n - 1.0d0))/8.0d0 + (e*n*(e + 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(3,3) = - (e*n*z*(e + 1.0d0)*(n + 1.0d0))/8.0d0 - (e*n*(e + 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(3,4) = (e*n*z*(e - 1.0d0)*(n + 1.0d0))/8.0d0 + (e*n*(e - 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(3,5) = (e*n*z*(e - 1.0d0)*(n - 1.0d0))/8.0d0 + (e*n*(e - 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(3,6) = - (e*n*z*(e + 1.0d0)*(n - 1.0d0))/8.0d0 - (e*n*(e + 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(3,7) = (e*n*z*(e + 1.0d0)*(n + 1.0d0))/8.0d0 + (e*n*(e + 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(3,8) = - (e*n*z*(e - 1.0d0)*(n + 1.0d0))/8.0d0 - (e*n*(e - 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(3,9) = - (e**2.0d0/4.0d0 - 1.0d0/4.0d0)*(n - 1.0d0) 
            DiffFunction(3,10) = (e**2.0d0/4.0d0 - 1.0d0/4.0d0)*(n + 1.0d0) 
            DiffFunction(3,11) = - (e**2.0d0/4.0d0 - 1.0d0/4.0d0)*(n + 1.0d0) 
            DiffFunction(3,12) = (e**2.0d0/4.0d0 - 1.0d0/4.0d0)*(n - 1.0d0) 
            DiffFunction(3,13) = - (e/4.0d0 - 1.0d0/4.0d0)*(n**2.0d0 - 1.0d0) 
            DiffFunction(3,14) = (e/4.0d0 + 1.0d0/4.0d0)*(n**2.0d0 - 1.0d0) 
            DiffFunction(3,15) = - (e/4.0d0 + 1.0d0/4.0d0)*(n**2.0d0 - 1.0d0) 
            DiffFunction(3,16) = (e/4.0d0 - 1.0d0/4.0d0)*(n**2.0d0 - 1.0d0) 
            DiffFunction(3,17) = - 2.0d0*z*(e/4.0d0 - 1.0d0/4.0d0)*(n - 1.0d0) 
            DiffFunction(3,18) = 2.0d0*z*(e/4.0d0 + 1.0d0/4.0d0)*(n - 1.0d0) 
            DiffFunction(3,19) = - 2.0d0*z*(e/4.0d0 + 1.0d0/4.0d0)*(n + 1.0d0) 
            DiffFunction(3,20) = 2.0d0*z*(e/4.0d0 - 1.0d0/4.0d0)*(n + 1.0d0) 
        else
            stop "ERROR, in DiffFormFunction, problem with ElementType"
        end if
    end subroutine DiffFormFunction
    ! 2. Elascitity Tensor
    subroutine ElasticityTensor(Self,ETensor)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                                   :: Self    
        double precision, dimension(:,:), allocatable, intent(inout)      :: ETensor
        double precision                                                  :: E,V,Constant
        E = Self%YoungModulus
        V = Self%PoissonModulus    
        if (Self%AnalysisType.eq.'PlaneStress') then
            allocate(ETensor(3,3))
            Constant = E/(1.0d0 - V**2.0d0)
            ETensor(1,:) = [1.0d0,V,0.0d0]
            ETensor(2,:) = [V,1.0d0,0.0d0]
            ETensor(3,:) = [0.0d0,0.0d0,(1.0d0-V)/2.0d0]
            ETensor = Constant*ETensor
        elseif (Self%AnalysisType.eq.'PlaneStrain') then
            allocate(ETensor(3,3))
            Constant = E/((1.0d0 + V)*(1.0d0 - 2.0d0*V))
            ETensor(1,:) = [1.0d0-V,V,0.0d0]
            ETensor(2,:) = [V,1.0d0-V,0.0d0]
            ETensor(3,:) = [0.0d0,0.0d0,(1.0d0-2.0d0*V)/2.0d0]
            ETensor = Constant*ETensor
        elseif (Self%AnalysisType.eq.'SolidIso') then
            allocate(ETensor(6,6))
            Constant = E*(1.0d0-v)/((1.0d0+V)*(1.0d0-2.0d0*V))
            ETensor(1,:) = [1.0d0,V/(1.0d0-v),V/(1.0d0-v),0.0d0,0.0d0,0.0d0]
            ETensor(2,:) = [V/(1.0d0-v),1.0d0,V/(1.0d0-v),0.0d0,0.0d0,0.0d0]
            ETensor(3,:) = [V/(1.0d0-v),V/(1.0d0-v),1.0d0,0.0d0,0.0d0,0.0d0]
            ETensor(4,:) = [0.0d0,0.0d0,0.0d0,(1.0d0-2.0d0*V)/(2.0d0*(1.0d0-v)),0.0d0,0.0d0]
            ETensor(5,:) = [0.0d0,0.0d0,0.0d0,0.0d0,(1.0d0-2.0d0*V)/(2.0d0*(1.0d0-v)),0.0d0]
            ETensor(6,:) = [0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,(1.0d0-2.0d0*V)/(2.0d0*(1.0d0-v))]
            ETensor = Constant*ETensor
        end if
    end subroutine ElasticityTensor
    ! 3. Node-Element interaction
    subroutine PreAssemblyRoutine(self)
        implicit none
        class(FEM_MMNL_Base), intent(inout)                           :: Self
        integer                                                     :: i,j,k,k_local,i1,IndexRow,IndexCol
        integer, dimension(:), allocatable                          :: InPosition
        logical, dimension(:), allocatable                          :: InLogical
        integer, dimension(:,:), allocatable                        :: Node_Interaction,Elem_Interaction
        ! -------------------------------------------------------------------------------------- !
        ! note: In this part, a preliminary assembly of the stiffness matrix (in sparse form)    !
        !       is made. The idea is to locate the position of each element of the matrix in     !
        !       the sparse vector so that the assembly of the global matrix in each iteration    !
        !       can be faster.                                                                   !
        ! -------------------------------------------------------------------------------------- !
        ! get de element and DoF incidence for each DoF (only the free degres)
        j = size(Self%FreeD)
        ! maximum of 50 elements per node
        allocate(Elem_Interaction(j,50))
        Elem_Interaction = 0
        ! considering the quad20 20*20 
        allocate(Node_Interaction(j,600))
        Node_Interaction = 0

        !$OMP PARALLEL DO PRIVATE(InLogical,InPosition,j,k) SHARED(Self,Elem_Interaction,Node_Interaction)
        do i = 1, size(Self%FreeD), 1
            ! element interaction
            InLogical = any(Self%ConnectivityD.eq.Self%FreeD(i),2)
            InPosition = pack([(k,k=1,Self%Ne)],InLogical)
            j = count(InLogical)
            Elem_Interaction(i,1:j) = InPosition 
            ! DoF interaction
            InPosition = reshape(Self%ConnectivityD(Elem_Interaction(i,1:j),:),[j*(Self%Npe)*Self%DimAnalysis])
            ! Removing fixed Dof
            do k = 1, size(InPosition), 1
                if(any(InPosition(k).eq.Self%FixedD)) InPosition(k) = 0
            end do
            ! Sorting ang eliminating 0s
            call sort(InPosition,j*Self%Npe*Self%DimAnalysis)
            InLogical = InPosition.ge.Self%FreeD(i)
            InPosition = pack(InPosition,InLogical)
            j = size(InPosition)
            Node_Interaction(i,1:j) = InPosition
        end do
        !$OMP END PARALLEL DO

        ! this print the interaction of the Free DoF with the others
        !call FilePrinting(Elem_Interaction,'DataResults/.InternalData/Elem_Interaction.txt')
        !call FilePrinting(Node_Interaction,'DataResults/.InternalData/Node_Interaction.txt')

        ! pre-assembly
        i = count(Node_Interaction.ne.0)
        allocate(Self%index_i_KGlobal(i)) ; Self%index_i_KGlobal = 0
        allocate(Self%index_j_KGlobal(i)) ; Self%index_j_KGlobal = 0
        allocate(Self%Location_KGlobal(Self%Ne,Self%Npe*Self%DimAnalysis,Self%Npe*Self%DimAnalysis))
        Self%Location_KGlobal = 0
        k = 1
        
        !do i = 1, size(Self%FreeD), 1                           ! Col
        !    do j = 1, count(Node_Interaction(i,:).ne.0), 1      ! row
        !        Self%index_i_KGlobal(k) = findloc(Self%FreeD,Node_Interaction(i,j),1)   ! row-indx
        !        Self%index_j_KGlobal(k) = i                                             ! col-indx
        !        do i1 = 1, count(Elem_Interaction(i,:).ne.0), 1
        !            IndexRow = findloc(Self%ConnectivityD(Elem_Interaction(i,i1),:),Node_Interaction(i,j),1)
        !            IndexCol = findloc(Self%ConnectivityD(Elem_Interaction(i,i1),:),Node_Interaction(i,1),1)
        !            if (IndexCol.eq.0.or.IndexRow.eq.0) cycle
        !            Self%Location_KGlobal(Elem_Interaction(i,i1),IndexRow,IndexCol) = k
        !        end do
        !        k = k + 1
        !    end do
        !end do
        do i = 1, size(Self%FreeD), 1                               ! Col
            do j = 1, count(Node_Interaction(i,:).ne.0), 1          ! Row
                k_local = k
                k = k + 1
                Self%index_i_KGlobal(k_local) = findloc(Self%FreeD,Node_Interaction(i,j),1) ! Row-indx
                Self%index_j_KGlobal(k_local) = i                                           ! Col-indx
                do i1 = 1, count(Elem_Interaction(i,:).ne.0), 1
                    IndexRow = findloc(Self%ConnectivityD(Elem_Interaction(i,i1),:),Node_Interaction(i,j),1)
                    IndexCol = findloc(Self%ConnectivityD(Elem_Interaction(i,i1),:),Node_Interaction(i,1),1)
                    if ((IndexCol.eq.0).or.(IndexRow.eq.0)) cycle
                    Self%Location_KGlobal(Elem_Interaction(i,i1),IndexRow,IndexCol) = k_local
                end do
                !k = k + 1
            end do
        end do

    end subroutine PreAssemblyRoutine

end module Base_FEA_MMNL_Module