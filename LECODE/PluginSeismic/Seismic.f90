module Seismic

    use file_data
    use mpidata
    use strata_data
    use sediment_data

    implicit none

    public

    integer :: SeisPts, LinPts

    real( tkind ), dimension( :,: ), allocatable :: seismic_data
    real( tkind ), dimension( :,: ), allocatable :: envelop

contains

    ! ============================================================================
    !> Subroutine define_seismic_points
    !! Find the seismic points.
    !<
    ! ============================================================================
    subroutine define_seismic_points

        logical :: record, xline

        integer :: k, gid, maxPts, p, n, ks, linePts, m

        integer, dimension( nproc ) :: sdisp, nums

        real( tkind ) :: minh, maxh, minv, maxv, miv, mih, mav, mah, mgz, poro

        real( tkind ),dimension( : ),allocatable :: hor, vert, dat, dat2, low, high, lineX
        real( tkind ),dimension( : ),allocatable :: ghor, gvert, gdat, gdat2, glow, ghigh, glineX

        if( iam == 0 ) print*,'Create seismic line'

        if( seismicX_min /= seismicX_max )then
            xline = .true.
            linePts = int( abs( seismicX_max - seismicX_min ) / strat_dx ) + 1
        else
            xline = .false.
            linePts = int( abs( seismicY_max - seismicY_min ) / strat_dx ) + 1
        endif
        maxPts = G_strat%total_layer * linePts
        if( allocated( low ) ) deallocate( low )
        if( allocated( high ) ) deallocate( high )
        if( allocated( lineX ) ) deallocate( lineX )
        allocate( low( linePts ), high( linePts ), lineX( linePts ) )

        minv = 1.e9_8
        maxv = -1.e9_8
        minh = 1.e9_8
        maxh = -1.e9_8

        if( allocated( hor ) ) deallocate( hor )
        if( allocated( vert ) ) deallocate( vert )
        if( allocated( dat ) ) deallocate( dat )
        if( allocated( dat2 ) ) deallocate( dat2 )
        allocate( hor( maxPts ), vert( maxPts ), dat( maxPts ), dat2( maxPts ) )

        p = 0
        m = 0
        do k = 1, L_strat%nodes_nb

            if( locv_coord( k, 1 ) >= seismicX_min .and. locv_coord( k, 1 ) <= seismicX_max )then
                if( locv_coord( k, 2 ) >= seismicY_min .and. locv_coord( k, 2 ) <= seismicY_max )then
                    gid = locv_gid( k )
                    record = .true.
                    if( gv_halo( gid ) > 1 )then
                        record = .false.
                        if( gv_ShareID( gid, 1 ) == iam )then
                            if( gv_ShareID( gid, 2 ) > iam ) record = .true.
                        elseif( gv_ShareID( gid, 2 ) == iam )then
                            if( gv_ShareID( gid, 1 ) > iam ) record = .true.
                        endif
                    endif

                    if( record )then
                        do n = 1 + InitDep, locv_NbLay( k )
                            if( strat_thick( k, n ) > 0.0_8 )then
                                p = p + 1
                                if( xline )then
                                    hor( p ) = locv_coord( k, 1 )
                                else
                                    hor( p ) = locv_coord( k, 2 )
                                endif
                                minh = min( hor( p ), minh )
                                maxh = max( hor( p ), maxh )
                                vert( p ) = strat_zelev( k, n )
                                minv = min( vert( p ), minv )
                                maxv = max( vert( p ), maxv )
                                mgz = 0.0_8
                                poro = 0.0_8
                                do ks = 1, totgrn
                                    mgz = real( mgz + sediment( ks )%density *  strat_sedh( k, n, ks ) )
                                    poro = poro + strat_porosity( k, n, ks ) *  strat_sedh( k, n, ks )
                                enddo
                                dat( p ) = real( mgz / strat_thick( k, n ) )
                                dat2( p ) = poro / strat_thick( k, n )
                            endif
                        enddo
                        m = m + 1
                        low( m ) = strat_zelev( k, 1+InitDep )
                        high( m ) = strat_zelev( k,  locv_NbLay( k ) )
                        lineX( m ) = locv_coord( k, 2 )
                        if( xline ) lineX( m ) = locv_coord( k, 1 )
                    endif

                endif
            endif

        enddo
        nums = 0
        call mpi_allgather( p, 1, int_type, nums, 1, int_type, SPModel_comm_world, ierr )
        call mpi_allreduce( p, SeisPts, 1, int_type, sum_type, SPModel_comm_world, ierr )
        sdisp( 1 ) = 0
        do k = 1, nproc - 1
            sdisp( k + 1 ) = sdisp( k ) + nums( k )
        enddo
        call mpi_allreduce( minv, miv, 1, dbl_type, min_type, SPModel_comm_world, ierr )
        call mpi_allreduce( minh, mih, 1, dbl_type, min_type, SPModel_comm_world, ierr )
        call mpi_allreduce( maxv, mav, 1, dbl_type, max_type, SPModel_comm_world, ierr )
        call mpi_allreduce( maxh, mah, 1, dbl_type, max_type, SPModel_comm_world, ierr )

        if( allocated( ghor ) ) deallocate( ghor )
        if( allocated( gvert ) ) deallocate( gvert )
        if( allocated( gdat ) ) deallocate( gdat )
        if( allocated( gdat2 ) ) deallocate( gdat2 )
        allocate( ghor( SeisPts ), gvert( SeisPts ), gdat( SeisPts ), gdat2( SeisPts ) )

        ! Merge all points to first processor
        call mpi_gatherv( hor, p, dbl_type, ghor, nums, sdisp, dbl_type, 0, SPModel_comm_world, ierr )
        call mpi_gatherv( vert, p, dbl_type, gvert, nums, sdisp, dbl_type, 0, SPModel_comm_world, ierr )
        call mpi_gatherv( dat, p, dbl_type, gdat, nums, sdisp, dbl_type, 0, SPModel_comm_world, ierr )
        call mpi_gatherv( dat2, p, dbl_type, gdat2, nums, sdisp, dbl_type, 0, SPModel_comm_world, ierr )

        nums = 0
        call mpi_allgather( m, 1, int_type, nums, 1, int_type, SPModel_comm_world, ierr )
        call mpi_allreduce( m, LinPts, 1, int_type, sum_type, SPModel_comm_world, ierr )
        sdisp( 1 ) = 0
        do k = 1, nproc - 1
            sdisp( k + 1 ) = sdisp( k ) + nums( k )
        enddo
        if( allocated( glow ) ) deallocate( glow )
        if( allocated( ghigh ) ) deallocate( ghigh )
        if( allocated( glineX ) ) deallocate( glineX )
        allocate( glow( LinPts ), ghigh( LinPts ), glineX( LinPts ) )

        ! Merge all points to first processor
        call mpi_gatherv( low, m, dbl_type, glow, nums, sdisp, dbl_type, 0, SPModel_comm_world, ierr )
        call mpi_gatherv( high, m, dbl_type, ghigh, nums, sdisp, dbl_type, 0, SPModel_comm_world, ierr )
        call mpi_gatherv( lineX, m, dbl_type, glineX, nums, sdisp, dbl_type, 0, SPModel_comm_world, ierr )

        if( iam == 0 )then
            if( allocated( seismic_data ) ) deallocate( seismic_data )
            allocate( seismic_data( SeisPts, 4  ) )
            seismic_data( 1:SeisPts, 1 ) = ghor
            seismic_data( 1:SeisPts, 2 ) = gvert
            seismic_data( 1:SeisPts, 3 ) = gdat
            seismic_data( 1:SeisPts, 4 ) = gdat2
            call create_seismic_line
            deallocate( seismic_data )
            print*,'Seismic line extension'
            print*,'Nb of points: ',SeisPts
            print*,'Horizontal min/max: ',mih,mah
            print*,'Vertical min/max: ',miv,mav
            if( allocated( envelop ) ) deallocate( envelop )
            allocate( envelop( LinPts, 3  ) )
            envelop( 1:LinPts, 1 ) = glineX
            envelop( 1:LinPts, 2 ) = glow
            envelop( 1:LinPts, 3 ) = ghigh
            call create_seismic_envelop
            deallocate( envelop )
        endif

        deallocate( hor, vert, dat, low, high, lineX )
        deallocate( ghor, gvert, gdat, glow, ghigh, glineX )

        return

    end subroutine define_seismic_points
    ! ============================================================================
    !> Subroutine create_seismic_line
    !! Create a seismic line.
    !<
    ! ============================================================================
    subroutine create_seismic_line

        character(len=128) :: file

        integer :: iunit, ios, k

        ! Create the top surface
        iunit = 20

        fseismic = 'seismic.csv'
        file = ''
        file = fseismic
        call addpath2( file )
        fseismic = file

        open(iunit,file=fseismic,status="replace",action="write",iostat=ios)
        rewind(iunit)

        do k = 1, SeisPts
            write(iunit,*) seismic_data( k, 1 ),',',seismic_data( k, 2 ),',',seismic_data( k, 3 ),&
                ',',seismic_data( k, 4 )
        enddo

        close( iunit )

        return

    end subroutine create_seismic_line
    ! ============================================================================
    !> Subroutine create_seismic_envelop
    !! Create envelop mask.
    !<
    ! ============================================================================
    subroutine create_seismic_envelop

        character(len=128) :: file

        integer :: iunit, ios, k

        ! Create the top surface
        iunit = 20

        fseismic = 'envelop.csv'
        file = ''
        file = fseismic
        call addpath2( file )
        fseismic = file

        open(iunit,file=fseismic,status="replace",action="write",iostat=ios)
        rewind(iunit)

        do k = 1, LinPts
            write(iunit,*) envelop( k, 1 ),',',envelop( k, 2 ),',',envelop( k, 3 )
        enddo

        close( iunit )

        return

    end subroutine create_seismic_envelop
    ! ============================================================================

end module Seismic
