module MMNL_Paraview_Module
    ! Extends the single-material Paraview post-processing (EnSight Gold format)
    ! to multi-material topology optimization results (MMNLTop type).
    !
    ! NOTE: this module is self-contained on purpose. It does NOT "use
    ! Paraview_Module" because that module's GetParaviewFiles subroutine
    ! depends on Optimization_module (your single-material module), which
    ! would drag an unrelated/unnecessary dependency in here. Instead, the
    ! two writer routines it needs (geometry + scalar-per-element file, both
    ! double precision) are reproduced verbatim below from Paraview_Module.f90.
    implicit none
    contains

    ! ---- verbatim copy of GenerateGeometryParaviewFileDP (Paraview_Module.f90) ----
    subroutine GenerateGeometryParaviewFileDP_MM(Path,Element,DimAnalysis,Connect,Coord)
        implicit none
        character(len=*), intent(in)                               :: Path
        character(len=*), intent(in)                               :: Element
        integer, intent(in)                                        :: DimAnalysis
        integer, dimension(:,:), allocatable, intent(in)           :: Connect
        double precision, dimension(:,:), allocatable, intent(in)  :: Coord
        integer                                                    :: i,ios,iounit
        open(newunit=iounit, file=Path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name_geometryfile"
            write(unit=iounit,fmt='(a)') 'This is the 1st description line of the EnSight Gold geometry example'
            write(unit=iounit,fmt='(a)') 'This is the 2st description line of the EnSight Gold geometry example'
            write(unit=iounit,fmt='(a)') 'node id given'
            write(unit=iounit,fmt='(a)') 'element id given'
            write(unit=iounit,fmt='(a)') 'extents'
            write(unit=iounit,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=iounit,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=iounit,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=iounit,fmt='(a)') 'part'
            write(unit=iounit,fmt='(a)') '1'
            if (DimAnalysis.eq.2) then
                write(unit=iounit,fmt='(a)') '2D uns-elements (description line for part 1)'
            elseif (DimAnalysis.eq.3) then
                write(unit=iounit,fmt='(a)') '3D uns-elements (description line for part 1)'
            else
                stop "ERROR DimAnalysis generating paraview geometry file"
            end if
            write(unit=iounit,fmt='(a)') 'coordinates'
            write(unit=iounit,fmt='(a)') ''
            write(unit=iounit,fmt=*) size(Coord,1)
            do i = 1, size(Coord,1), 1
                write(unit=iounit,fmt=*) i
            end do
            write(unit=iounit,fmt='(a)') ''
            do i = 1, size(Coord,1), 1
                write(unit=iounit,fmt=*) Coord(i,1)
            end do
            write(unit=iounit,fmt='(a)') ''
            do i = 1, size(Coord,1), 1
                write(unit=iounit,fmt=*) Coord(i,2)
            end do
            write(unit=iounit,fmt='(a)') ''
            do i = 1, size(Coord,1), 1
                if (DimAnalysis.eq.2) then
                    write(unit=iounit,fmt=*) 0.0d0
                elseif (DimAnalysis.eq.3) then
                    write(unit=iounit,fmt=*) Coord(i,3)
                else
                    stop "ERROR DimAnalysis generating paraview geometry file"
                end if
            end do
            write(unit=iounit,fmt='(a)') ''
            write(unit=iounit,fmt='(a)') Element
            write(unit=iounit,fmt=*) size(Connect,1)
            do i = 1, size(Connect,1), 1
                write(unit=iounit,fmt=*) i
            end do
            write(unit=iounit,fmt='(a)') ''
            do i = 1, size(Connect,1), 1
                write(unit=iounit,fmt=*) Connect(i,:)
            end do
        close(iounit)
    end subroutine GenerateGeometryParaviewFileDP_MM

    ! ---- verbatim copy of GenerateEscalarParaviewFileDP (Paraview_Module.f90) ----
    subroutine GenerateEscalarParaviewFileDP_MM(Path,Element,Result)
        implicit none
        character(len=*), intent(in)                               :: Path
        character(len=*), intent(in)                               :: Element
        double precision, dimension(:), allocatable, intent(in)    :: Result
        integer                                                    :: i,ios,iounit
        open(newunit=iounit, file=Path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name_Coordfile"
            write(unit=iounit,fmt='(a)') 'Scalar File'
            write(unit=iounit,fmt='(a)') 'part'
            write(unit=iounit,fmt='(a)') '1'
            write(unit=iounit,fmt='(a)') Element
            do i = 1, size(Result), 1
                write(unit=iounit,fmt=*) Result(i)
            end do
        close(iounit)
    end subroutine GenerateEscalarParaviewFileDP_MM

    subroutine MMNLParaviewPostProcessing(Self,Path,ResultName)
        use MMNL_Optimization_Module
        implicit none
        class(MMNLTop), intent(inout)                                :: Self
        character(len=*), intent(in)                               :: Path        ! e.g. 'output/paraview'
        character(len=*), intent(in)                               :: ResultName  ! e.g. 'MultiMaterialResult'
        ! internal variables
        double precision, dimension(:), allocatable                :: ResultDP
        character(len=200)                                         :: Path1,CaseName,ScalarFileName
        character(len=20)                                          :: MatLabel
        integer                                                    :: i,ios,iounit,NKeep
        logical, dimension(:), allocatable                         :: KeepMask
        integer, dimension(:,:), allocatable                       :: ConnectivityFiltered

        ! ---------------------------------------------------------------
        ! 0. Make sure the output folder exists (e.g. output/paraview)
        ! ---------------------------------------------------------------
        call system('mkdir -p ' // trim(Path))

        ! ---------------------------------------------------------------
        ! 1. Case file (the one you open in ParaView): geometry + N+1 scalars
        !    (MaterialIndex, plus one projected density field per material)
        ! ---------------------------------------------------------------
        CaseName = trim(Path) // '/' // trim(ResultName) // '.case'
        open(newunit=iounit, file=CaseName, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name_casefile (MMNLParaviewPostProcessing)"
            write(unit=iounit,fmt='(a)') 'FORMAT'
            write(unit=iounit,fmt='(a)') 'type: ensight gold'
            write(unit=iounit,fmt='(a)') ''
            write(unit=iounit,fmt='(a)') 'GEOMETRY'
            write(unit=iounit,fmt='(a)') 'model:     ' // trim(ResultName) // '.geom'
            write(unit=iounit,fmt='(a)') ''
            write(unit=iounit,fmt='(a)') 'VARIABLE'
            ! -- material index (0 = void, i = material i, from FinalTopologyMulti) --
            ScalarFileName = trim(ResultName) // '_MaterialIndex' // '.esca'
            write(unit=iounit,fmt='(a)') 'scalar per element:   MaterialIndex ' // trim(ScalarFileName)
            ! -- one projected-density field per material --
            do i = 1, Self%NMaterial, 1
                write(MatLabel,'(A,I0)') 'Density_Material', i
                ScalarFileName = trim(ResultName) // '_' // trim(MatLabel) // '.esca'
                write(unit=iounit,fmt='(a)') 'scalar per element:   ' // trim(MatLabel) // ' ' // trim(ScalarFileName)
            end do
        close(iounit)

        ! ---------------------------------------------------------------
        ! 2. Geometry file. Elements whose max density is below the
        !    user-defined PostProcesingFilter (SetPostProcesingFilterTO)
        !    are EXCLUDED entirely from the mesh written to .geom -- not
        !    just relabeled -- so very-low-density "gray" leftovers do not
        !    clutter the ParaView view. Uses the same max(xProj) >= filter
        !    criterion as FinalTopologyMulti, so MaterialIndex and the
        !    geometry stay consistent with each other.
        ! ---------------------------------------------------------------
        if (.not.allocated(Self%MaterialIndex)) then
            write(*,*) 'WARNING: Self%MaterialIndex not allocated -> calling FinalTopologyMulti(Self) now.'
            call FinalTopologyMulti(Self)
        end if
        KeepMask = Self%MaterialIndex.gt.0
        NKeep = count(KeepMask)
        if (NKeep.eq.0) then
            write(*,*) 'WARNING: PostProcesingFilter = ', Self%PostProcesingFilter, &
                       ' removed ALL elements. Nothing to write to ParaView.'
            return
        end if
        allocate(ConnectivityFiltered(NKeep,Self%Npe))
        ConnectivityFiltered = Self%ConnectivityN(pack([(i,i=1,Self%Ne)],KeepMask),:)

        Path1 = trim(Path) // '/' // trim(ResultName) // '.geom'
        call GenerateGeometryParaviewFileDP_MM(Path1,Self%ElementType,Self%DimAnalysis,ConnectivityFiltered,Self%Coordinates)

        ! ---------------------------------------------------------------
        ! 3. Scalar files (filtered the same way, so they line up 1-to-1
        !    with ConnectivityFiltered / the .geom element order)
        ! ---------------------------------------------------------------
        ! -- material index (integer -> double, so it plots as a clean step field) --
        allocate(ResultDP(NKeep))
        ResultDP = dble(pack(Self%MaterialIndex,KeepMask))
        Path1 = trim(Path) // '/' // trim(ResultName) // '_MaterialIndex' // '.esca'
        call GenerateEscalarParaviewFileDP_MM(Path1,Self%ElementType,ResultDP)
        deallocate(ResultDP)

        ! -- one .esca per material, using the projected density x~~_i (Self%xProj) --
        do i = 1, Self%NMaterial, 1
            allocate(ResultDP(NKeep))
            ResultDP = pack(Self%xProj(:,i),KeepMask)
            write(MatLabel,'(A,I0)') 'Density_Material', i
            Path1 = trim(Path) // '/' // trim(ResultName) // '_' // trim(MatLabel) // '.esca'
            call GenerateEscalarParaviewFileDP_MM(Path1,Self%ElementType,ResultDP)
            deallocate(ResultDP)
        end do
    end subroutine MMNLParaviewPostProcessing

end module MMNL_Paraview_Module