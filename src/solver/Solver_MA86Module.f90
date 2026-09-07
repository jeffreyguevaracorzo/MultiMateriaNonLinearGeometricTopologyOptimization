module Solver_MA86Module
    implicit none
    ! Lapack, BLAS, OpenMP and Metis are needed in some routines/functions
    Interface SparseSystemMA86Solver
        module procedure RealSparseSystemSolverMA86_BVector
        module procedure RealSparseSystemSolverMA86_BMatrix
        module procedure DPSparseSystemSolverMA86_BMatrix
        module procedure DPSparseSystemSolverMA86_BVector
    end interface SparseSystemMA86Solver

    contains
    ! ----------------- ALL FUNCTIONS AND SUBROUTINES -----------------
    ! Sparse solver for a Non-definide system
    ! 2.1 Real system (B Vector)
    function RealSparseSystemSolverMA86_BVector(RowVectorA,ColVectorA,ValueVectorA,ValueVectorB) result(X)
        use hsl_ma86_single
        use hsl_mc68_single
        use hsl_mc69_single
        implicit none
        type(mc68_control)                                              :: control68
        type(mc68_info)                                                 :: info68
        type(ma86_keep)                                                 :: keep
        type(ma86_control)                                              :: control
        type(ma86_info)                                                 :: info
        integer                                                         :: n, ne, lmap, flag
        integer, dimension(:), allocatable                              :: ptr, row, order, map
        integer, dimension(:), allocatable, intent(in)                  :: ColVectorA, RowVectorA
        real, dimension(:), allocatable, intent(in)                     :: ValueVectorA, ValueVectorB
        real, dimension(:), allocatable                                 :: val,x
        n = size(ValueVectorB)
        ne = size(ValueVectorA)
        X = ValueVectorB
        ! Convert to HSL standard format
        allocate(ptr(n+1))
        call mc69_coord_convert(HSL_MATRIX_REAL_SYM_INDEF, n, n, ne, RowVectorA, ColVectorA, &
            ptr, row, flag, val_in=ValueVectorA, val_out=val, lmap=lmap, map=map)
        call stop_on_bad_flag("mc69_coord_convert", flag)
        ! Call mc68 to find a fill reducing ordering (1=AMD)
        allocate(order(n))
        call mc68_order(1, n, ptr, row, order, control68, info68)
        call stop_on_bad_flag("mc68_order", info68%flag)
        ! Analyse
        call ma86_analyse(n, ptr, row, order, keep, control, info)
        call stop_on_bad_flag("analyse", info%flag)
        ! Factor
        call ma86_factor(n, ptr, row, val, order, keep, control, info)
        call stop_on_bad_flag("factor", info%flag)
        ! Solve
        call ma86_solve(x, order, keep, control, info)
        call stop_on_bad_flag("solve", info%flag)
        ! Finalize
        call ma86_finalise(keep, control)
    end function RealSparseSystemSolverMA86_BVector
    ! 2.2 Real system (B Matrix)
    function RealSparseSystemSolverMA86_BMatrix(RowVectorA,ColVectorA,ValueVectorA,ValueVectorB) result(X)
        use hsl_ma86_single
        use hsl_mc68_single
        use hsl_mc69_single
        implicit none
        type(mc68_control)                                              :: control68
        type(mc68_info)                                                 :: info68
        type(ma86_keep)                                                 :: keep
        type(ma86_control)                                              :: control
        type(ma86_info)                                                 :: info
        integer                                                         :: n, ne, nrhs, lmap, flag
        integer, dimension(:), allocatable                              :: ptr, row, order, map
        integer, dimension(:), allocatable, intent(in)                  :: ColVectorA, RowVectorA
        real, dimension(:), allocatable, intent(in)                     :: ValueVectorA
        real, dimension(:,:), allocatable, intent(in)                   :: ValueVectorB
        real, dimension(:), allocatable                                 :: val
        real, dimension(:,:), allocatable                               :: X
        n = size(ValueVectorB,1)
        ne = size(ValueVectorA)
        nrhs = size(ValueVectorB,2)
        X = ValueVectorB
        ! Convert to HSL standard format
        allocate(ptr(n+1))
        call mc69_coord_convert(HSL_MATRIX_REAL_SYM_INDEF, n, n, ne, RowVectorA, ColVectorA, &
            ptr, row, flag, val_in=ValueVectorA, val_out=val, lmap=lmap, map=map)
        call stop_on_bad_flag("mc69_coord_convert", flag)
        ! Call mc68 to find a fill reducing ordering (1=AMD)
        allocate(order(n))
        call mc68_order(1, n, ptr, row, order, control68, info68)
        call stop_on_bad_flag("mc68_order", info68%flag)
        ! Analyse
        call ma86_analyse(n, ptr, row, order, keep, control, info)
        call stop_on_bad_flag("analyse", info%flag)
        ! Factor and solve
        call ma86_factor_solve(n, ptr, row, val, order, keep, control, info, nrhs, n, X)
        call stop_on_bad_flag("factor", info%flag)
        ! Finalize
        call ma86_finalise(keep, control)
    end function RealSparseSystemSolverMA86_BMatrix
    ! 2.3 DP system (B Vector)
    function DPSparseSystemSolverMA86_BVector(RowVectorA,ColVectorA,ValueVectorA,ValueVectorB) result(X)
        use hsl_ma86_double
        use hsl_mc68_double
        use hsl_mc69_double
        implicit none
        type(mc68_control)                                              :: control68
        type(mc68_info)                                                 :: info68
        type(ma86_keep)                                                 :: keep
        type(ma86_control)                                              :: control
        type(ma86_info)                                                 :: info
        integer                                                         :: n, ne, lmap, flag
        integer, dimension(:), allocatable                              :: ptr, row, order, map
        integer, dimension(:), allocatable                              :: ColVectorA, RowVectorA
        double precision, dimension(:), allocatable, intent(in)         :: ValueVectorA, ValueVectorB
        double precision, dimension(:), allocatable                     :: val,X
        n = size(ValueVectorB)
        ne = size(ValueVectorA)
        X = ValueVectorB
        ! Convert to HSL standard format
        allocate(ptr(n+1))
        call mc69_coord_convert(HSL_MATRIX_REAL_SYM_INDEF, n, n, ne, RowVectorA, ColVectorA, &
            ptr, row, flag, val_in=ValueVectorA, val_out=val, lmap=lmap, map=map)
        call stop_on_bad_flag("mc69_coord_convert", flag)
        ! Call mc68 to find a fill reducing ordering (1=AMD)
        allocate(order(n))
        call mc68_order(1, n, ptr, row, order, control68, info68)
        call stop_on_bad_flag("mc68_order", info68%flag)
        ! Analyse
        call ma86_analyse(n, ptr, row, order, keep, control, info)
        call stop_on_bad_flag("analyse", info%flag)
        ! Factor
        call ma86_factor(n, ptr, row, val, order, keep, control, info)
        call stop_on_bad_flag("factor", info%flag)
        ! Solve
        call ma86_solve(X, order, keep, control, info)
        call stop_on_bad_flag("solve", info%flag)
        ! Finalize
        call ma86_finalise(keep, control)
    end function DPSparseSystemSolverMA86_BVector
    ! 2.4 DP system (B matrix)
    function DPSparseSystemSolverMA86_BMatrix(RowVectorA,ColVectorA,ValueVectorA,ValueVectorB) result(X)
        use hsl_ma86_double
        use hsl_mc68_double
        use hsl_mc69_double
        implicit none
        type(mc68_control)                                              :: control68
        type(mc68_info)                                                 :: info68
        type(ma86_keep)                                                 :: keep
        type(ma86_control)                                              :: control
        type(ma86_info)                                                 :: info
        integer                                                         :: n, ne, nrhs, lmap, flag
        integer, dimension(:), allocatable                              :: ptr, row, order, map
        integer, dimension(:), allocatable, intent(in)                  :: ColVectorA, RowVectorA
        double precision, dimension(:), allocatable, intent(in)         :: ValueVectorA
        double precision, dimension(:,:), allocatable, intent(in)       :: ValueVectorB
        double precision, dimension(:), allocatable                     :: val
        double precision, dimension(:,:), allocatable                   :: X
        n = size(ValueVectorB,1)
        ne = size(ValueVectorA)
        nrhs = size(ValueVectorB,2)
        X = ValueVectorB
        ! Convert to HSL standard format
        allocate(ptr(n+1))
        call mc69_coord_convert(HSL_MATRIX_REAL_SYM_INDEF, n, n, ne, RowVectorA, ColVectorA, &
            ptr, row, flag, val_in=ValueVectorA, val_out=val, lmap=lmap, map=map)
        call stop_on_bad_flag("mc69_coord_convert", flag)
        ! Call mc68 to find a fill reducing ordering (1=AMD)
        allocate(order(n))
        call mc68_order(1, n, ptr, row, order, control68, info68)
        call stop_on_bad_flag("mc68_order", info68%flag)
        ! Analyse
        call ma86_analyse(n, ptr, row, order, keep, control, info)
        call stop_on_bad_flag("analyse", info%flag)
        ! Factor and solve
        call ma86_factor_solve(n, ptr, row, val, order, keep, control, info, nrhs, n, X)
        call stop_on_bad_flag("factor", info%flag)
        ! Finalize
        call ma86_finalise(keep, control)
    end function DPSparseSystemSolverMA86_BMatrix
    ! Warring signal
    subroutine stop_on_bad_flag(context, flag)
        character(len=*), intent(in)                                 :: context
        integer, intent(in)                                          :: flag
        if(flag.eq.0) return
        write(*,*) "Failure during ", context, " with flag = ", flag
        stop
    end subroutine stop_on_bad_flag

end module Solver_MA86Module