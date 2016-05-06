! ============================================================================
! Name        : Rain_init.f90
! Author      : Tristan Salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file Rain_init.f90
!!
!! Rain_init computes the inflow into the system according to rainfall.
!<
! ============================================================================
module rain_initial

    use file_data
    use mpidata
    use flux_data
    use time_data
    use fwalker_data
    use forces_data
    use error_data
    use strata_data
    use param_data
    use strata_data
    use mod_sort
    use TIN_data
    use strata_ini
    use TIN_function

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine read_rain_file
    !! Subroutine read_rain_file reads the rainfall file.
    !<
    ! ============================================================================
    subroutine read_rain_file

        ! Parameters declaration
        logical :: found
        integer :: k, n, i, maxrainarea
        integer :: iu, ios

        character(len=128) :: text

        ! In case the rain is defined as a global map we allocate the rainfall values for each
        ! stratal mesh points
        if( .not. rain_region )then

            do n = 1, rain_event

                ! Open rainfall map file
                iu = 49 + iam
                inquire( file=frainmap( n ), exist=found )
                if( found) then
                    open( iu, file=frainmap( n ), status="old", action="read", iostat=ios )
                    rewind( iu )
                else
                    attempt = RAIN_FILE
                endif

                ! Ensure that all processes have managed to open the file.
                call completion

                maxrain = G_strat%nodes_nb
                rain_areaNb( n ) = 1

                do k = 1, G_strat%nodes_nb
                    read(iu,*) i, rain_h( n , k )
                    if( rain_h( n , k ) < 0.0_8 .or. rain_h( n , k ) > 1000.0_8 )then
                        print*,'Problem in rainfall map file',n,'at vertex',k,'rain value equal',rain_h( n , k )
                        stop
                    endif
                enddo

                ! Close file
                close(iu)
            enddo

            ! Exit subroutine
            return

        endif

        ! Open rainfall declaration file when rain is defined by regions
        iu = 49 + iam
        inquire( file=frain, exist=found )
        if( found) then
            open( iu, file=frain, status="old", action="read", iostat=ios )
            rewind( iu )
        else
            attempt = RAIN_FILE
        endif

        ! Ensure that all processes have managed to open the file.
        call completion

        ! Get the total number of rainfall events
        k = 0
        do
            read(iu,'(a128)',iostat=ios,end=30) text
            if( text(1:1) == '#')then
                k = k + 1
            endif
        enddo
30  continue
        rewind(iu)
        rain_event = k

        ! Clip the maximum number of rainfall points
        if( G_strat%nodes_nb < maxrain ) maxrain = G_strat%nodes_nb

        ! Allocate rain arrays
        allocate( rain_areaNb( rain_event ) )
        allocate( rain_tend( rain_event ) )
        allocate( rain_tstart( rain_event ) )

        ! Read the rain region file and register events time information
        k = 0
        do
            if( k == 0 .or. k == rain_event )read(iu,'(a128)',iostat=ios,end=50) text
            if( text(1:1) == '*' ) goto 50
            if( text(1:1) == '#' )then
                k = k + 1
                i = len(text)
                read(text(2:i),*,iostat=ios) rain_tstart( k ), rain_tend( k )
                if( ios /= 0 )then
                    attempt = RAIN_CYC
                endif
                n = 0
                do
                    read(iu,'(a128)',iostat=ios,end=40) text
                    if( text(1:1) == '#' ) goto 20
                    if( text(1:1) == '*' ) goto 50
                    n = n + 1
                    rain_areaNb( k ) = n
                enddo
            endif
20          continue
        enddo
40      continue
50      continue

        ! Ensure that all processes have managed to read properly the file.
        call completion
        rewind(iu)

        ! Clip maximum rainfall points by region
        maxrainarea = 1
        do k = 1, rain_event
            n = rain_areaNb( k )
            maxrainarea = max( maxrainarea, n )
        enddo
        if( maxrainarea >  G_strat%nodes_nb )  maxrainarea = G_strat%nodes_nb

        ! Allocate min, max arrays and rainfall values
        allocate( rain_xmin( rain_event, maxrainarea ) )
        allocate( rain_xmax( rain_event, maxrainarea ) )
        allocate( rain_ymin( rain_event, maxrainarea ) )
        allocate( rain_ymax( rain_event, maxrainarea ) )
        allocate( rain_h( rain_event, maxrainarea ) )
        allocate( rain_min( rain_event, maxrainarea ) )
        allocate( rain_max( rain_event, maxrainarea ) )

        ! Read the rain region file and register rainfall parameters information
        k = 0
        do
            read(iu,'(a128)',iostat=ios,end=60) text
            if( text(1:1) == '#' )then
                k = k + 1
                do n = 1, rain_areaNb( k )
                    read(iu,*,iostat=ios,end=60) rain_xmin( k, n ), rain_xmax( k, n ), rain_ymin( k , n ), &
                        rain_ymax( k , n ), rain_h( k , n ), rain_min( k , n ), rain_max( k , n )
                    if( rain_min( k , n ) > rain_max( k , n ) .or. rain_max( k , n ) < 0 .or. &
                        rain_max( k , n ) > 1.0e6_8 )then
                        print*,'Problem in rainfall region event',k,'region',n,'rain max elevation value of',rain_max( k , n )
                        stop
                    endif
                enddo
            endif
        enddo
60      continue

        ! Close rainfall region file
        close(iu)

        return

    end subroutine read_rain_file
    ! ============================================================================
    !> Subroutine register_rain_arrays
    !! Subroutine register_rain_arrays allocate the rainfall values according to the
    !! simulation time.
    !<
    ! ============================================================================
    subroutine register_rain_arrays

        ! Variables declaration
        integer :: k, k1, k2, know
        integer :: temp_src_nb, temp_raincount, temp_rainnb
        integer,dimension( : ),allocatable :: rainPt,  rainRegion, lpickRegion, lpickPt

        ! Rain event index for the actual simulation time
        know = 0
        rain_event_nb: do k = 1, rain_event
            if( rain_tstart( k ) <= tnow .and. rain_tend( k ) > tnow )then
                know = k
                exit rain_event_nb
            endif
        enddo rain_event_nb

        ! In case there is no rain event for the considered period exit this function
        if( know == 0 .or. know > rain_event ) return

        ! Find nodes within the region
        allocate( rainPt( L_strat%nodes_nb ), rainRegion( L_strat%nodes_nb ) )
        temp_rainnb = 0
        temp_raincount = 0

        ! In case rain is defined by regions
        if( rain_region )then
            rainRegion = -1
            ! For each points in the local stratal grid
            do k = 1, L_strat%nodes_nb
                k1 = locv_gid( k )

                ! Loop over the rainfall region
                rain_region_nb: do k2 = 1, rain_areaNb( know )

                    ! Elevation needs to be above sea level
                    if( elev_record( k1 ) > gsea%actual_sea )then

                        ! Clip elevation based on rainfall min and max elevation
                        if( elev_record( k1 ) >= rain_min( know , k2 ) .and. elev_record( k1 ) <= rain_max( know , k2 ) )then

                            ! Clip position based on rain region extent
                            if( gv_coord( k1, 1 ) <= rain_xmax( know , k2 ) .and. gv_coord( k1, 1 ) >= rain_xmin( know , k2 ) )then
                                if( gv_coord( k1, 2 ) <= rain_ymax( know , k2 ) .and. &
                                    gv_coord( k1, 2 ) >= rain_ymin( know , k2 ) )then

                                    ! Increment number of rain points
                                    temp_raincount = temp_raincount + 1

                                    ! Check flow accumulation value against rain accumulation parameters
                                    ! and prevent rainfall point to be initially positioned on simulation border
                                    if( facc( k1 ) <= rain_facc .and. gv_coord( k1, 1 ) > strat_xo .and. &
                                        gv_coord( k1, 2 ) > strat_yo .and. gv_coord( k1, 1 ) < strat_xo + strat_dx &
                                        * ( strat_X - 1 ) .and. gv_coord( k1, 2 ) < strat_yo + strat_dx * ( strat_Y - 1 ) )then
                                        temp_rainnb = temp_rainnb + 1
                                        rainPt( temp_rainnb ) = k1
                                        rainRegion( temp_rainnb ) = k2
                                        exit rain_region_nb
                                    endif
                                endif
                            endif
                        endif
                    endif
                enddo rain_region_nb
            enddo

        ! In case rain is defined using a rainfall map
        else

            rainRegion = 1

            ! For each points in the local stratal grid
            do k = 1, L_strat%nodes_nb
                k1 = locv_gid( k )

                ! Elevation needs to be above sea level and rainfall value positive
                if( elev_record( k1 ) >= gsea%actual_sea .and. rain_h( know , k1 ) > 0.0_8 )then

                    ! Increment number of rain points
                    temp_raincount = temp_raincount + 1

                    ! Check flow accumulation value against rain accumulation parameters
                    ! and prevent rainfall point to be initially positioned on simulation border
                    if( facc( k1 ) <= rain_facc .and. gv_coord( k1, 1 ) > strat_xo .and. &
                        gv_coord( k1, 2 ) > strat_yo .and. gv_coord( k1, 1 ) < strat_xo + strat_dx * &
                        ( strat_X - 1 ) .and. gv_coord( k1, 2 ) < strat_yo + strat_dx * ( strat_Y - 1 ) )then
                        temp_rainnb = temp_rainnb + 1
                        rainPt( temp_rainnb ) = k1
                    endif
                endif
            enddo
        endif

        ! Transfer values globally
        call mpi_allreduce( temp_rainnb, rain_src, 1, int_type, sum_type, SPModel_comm_world, ierr )
        call mpi_allreduce( temp_raincount, rain_total, 1, int_type, sum_type, SPModel_comm_world, ierr )

        ! In case there is no rainfall sources exit function
        if( rain_src == 0 )then
            deallocate(  rainPt, rainRegion )
            return
        endif

        ! Allocate rainfall arrays for flow walkers definition
        if( allocated( frain_node ) )then
            deallocate( frain_node, frain_xposition, frain_yposition, frain_qfl )
        endif
        allocate( frain_node( temp_rainnb  ) )
        allocate( frain_xposition( temp_rainnb  ) )
        allocate( frain_yposition( temp_rainnb  ) )
        allocate( frain_qfl( temp_rainnb  ) )
        frain_node = -1

        ! Declare arrays to pick up randomly some points within the selected group
        if( allocated( lpickPt ) ) deallocate( lpickPt )
        allocate( lpickPt( temp_rainnb ), lpickRegion( temp_rainnb ) )

        ! Check if the number of points where possible rain could be initialise is
        ! greater than the maximum number of rain FW
        if( rain_src >= rain_nb )then
            k = rain_src - rain_nb
            temp_src_nb = temp_rainnb - int( k * real( temp_rainnb ) / real( rain_src ) + 1 )
            if( temp_src_nb < 1 ) temp_src_nb = 0
            call random_pickup( temp_rainnb, rainPt, rainRegion, lpickPt, lpickRegion, temp_src_nb )
        ! There is sufficient points just pick some of them
        else
            call random_pickup( temp_rainnb, rainPt, rainRegion, lpickPt, lpickRegion, temp_rainnb )
            temp_src_nb = temp_rainnb
        endif

        ! Create local array of rain flow walkers based on random pick up
        temp_rainnb = 0
        do k = 1, temp_src_nb
            k1 = lpickPt( k )
            k2 = lpickRegion( k )
            if( k1 > 0 )then
                temp_rainnb = temp_rainnb + 1
                ! Position
                frain_xposition( temp_rainnb ) = gv_coord( k1, 1 )
                frain_yposition( temp_rainnb ) = gv_coord( k1, 2 )
                frain_node( temp_rainnb ) = k1
                ! If rainfall is declared by region
                if( rain_region )then
                    frain_qfl( temp_rainnb ) = dble( rain_h( know , k2 ) )
                ! Otherwise this is a rainfall map
                else
                    frain_qfl( temp_rainnb ) = dble( rain_h( know , k1 ) )
                endif
            endif
        enddo

        ! Deallocate some arrays
        deallocate(  rainPt, rainRegion, lpickPt, lpickRegion )

        ! Allocate new flow walkers
        rain_src = temp_rainnb
        rain_total = temp_raincount

        ! Create the flow walker based on the rainfall distribution
        call define_rain_flow_walker( know )

        return

    end subroutine register_rain_arrays
    ! ============================================================================
    !> Subroutine define_rain_flow_walker
    !! Subroutine define_rain_flow_walker calculates the inflow by allocating new flow
    !! walkers that flow into the system.
    !<
    ! ============================================================================
    subroutine define_rain_flow_walker( rain_cycle )

        ! Variables declaration
        integer :: reloc, relocT
        integer :: i, k, kn, ks, id, nsourceactive
        integer :: rain_cycle, i1, m, seed, idS
        integer,dimension( : ), allocatable :: fce, fcS

        real( tkind ) :: randomVal, qflw, flowtime, fv
        real( tkind ) :: vel, minbas, distX, distY, dist, newwidth

        sourceactive = 0

        ! Localise flow walkers in the TIN grid using a kdtree search
        allocate( fce( rain_src ), fcS( rain_src ) )
        call find_TINfaces_kdtree( rain_src, frain_xposition, frain_yposition, mxnghb, fce )

        ! Check the faces validity
        do ks = 1, rain_src
            if( fce( ks ) < 0 )then
                print*,'Something went wrong when looking for flow walkers faces in the TIN.'
                stop
            endif
        enddo

        ! Find the closest stratal node for each flow walker
        fcS = -1
        do ks = 1, rain_src
            call closest_stratal_node( frain_xposition( ks ), frain_yposition( ks ), fcS( ks ), idS )
        enddo

        ! Determine the duration of rainfall
        flowtime = rain_int
        if( rain_tend( rain_cycle ) - rain_tstart( rain_cycle ) < rain_int )&
            flowtime = ( rain_tend( rain_cycle ) - rain_tstart( rain_cycle ) )

        ! Check that the number of flow walkers is not
        ! greater than the maximum allowed
        reloc = 0
        relocT = 0
        k = num_fw + rain_src
        if( k > maxfw ) reloc = 1
        call mpi_allreduce( reloc, relocT, 1, int_type, max_type, SPModel_comm_world, ierr )

        ! If this is the case we need to reallocate some space for flow walkers
        if( relocT == 1 ) call reallocate_flow_walker_space( k )

        ! Register rain as active flow walkers
        k = num_fw
        do i = 1, rain_src
            sourceactive = 1
            k = k + 1
            if( k > maxfw )then

                attempt = IN_MAXFW
                exit

            else

                ! Main parameters defining flow walker
                fa_refid( k ) = -1
                fa_faceID( k , 1:numold ) = 0
                fa_xpos( k ) = real( frain_xposition( i ) )
                fa_ypos( k ) = real( frain_yposition( i ) )

                ! Rain water flux
                qflw = frain_qfl( i ) * strat_dx**2 * rain_total / rain_src
                qflw = qflw * flowtime

                ! Flow walker volume, height and discharge
                fa_volume( k ) = real( qflw )
                fa_h( k ) =  rain_elev
                fa_density( k ) = fluid_density
                fa_discharge( k ) = frain_qfl( i ) * strat_dx**2 &
                    * rain_total / ( rain_src * secyear )

                ! Initially the rain flow walker doesn't transport any sediment
                fa_sedcharge( k , : ) = 0.0_8
                fa_count( k ) = 0
                fa_type( k ) = tprain
                fa_pgrad( k ) = 0.0_8

                ! Flow walker velocity
                vel = transport%fvmin * 10.0_8
                call random_seed( size=seed )
                call random_number( harvest=randomVal )
                vel = 2.0_8 * vel * ( randomVal - 0.5_8 )

                ! Find the minimal elevation in the neighborhood
                id = frain_node( i )
                minbas =  elev_record( id )
                m = 0
                do kn = 1, 8
                    i1 = gv_ngbID( id, kn )
                    if( i1 > 0 )then
                        if( minbas >  elev_record( i1 ) )then
                            minbas = elev_record( i1 )
                            m = i1
                        endif
                    endif
                enddo

                ! Get the distance between the rain point and the minimum elevation
                distX = 1.0_8
                distY = 1.0_8
                dist = 1.0_8
                if( i1 > 0 .and. m > 0 )then
                    distX = gv_coord( m, 1 ) - gv_coord( id, 1 )
                    distY = gv_coord( m, 2 ) - gv_coord( id, 2 )
                    dist = sqrt( distX**2 + distY**2 )
                endif

                ! Update velocity vector
                fa_xvel( k ) = dble( vel * distX / dist )
                fa_yvel( k ) = dble( vel * distY / dist )

                ! Calculate the Euclidean norm of the FW's current velocity.
                fv = sqrt( fa_xvel( k )**2 + fa_yvel( k )**2 )

                ! Calculate FW's width
                newwidth = dble( fa_discharge( k ) / ( fv * fa_h( k ) ) )
                newwidth = min( newwidth, transport%wdtmax )
                newwidth = max( newwidth, transport%wdtmin )
                fa_width( k ) = newwidth
                if( fce( i ) > G_tin%faces_nb ) print*,'Something went wrong when allocating rain face.'
                if( fce( i ) > 0 ) fa_inbound( k ) = 0
                if( fce( i ) <= 0 ) fa_inbound( k ) = 1
                fa_faceID( k , 1 ) = fce( i )
                if( fce( i ) <= 0 ) print*,'Something went wrong when declaring rain flow walkers.'
                if( fce( i ) <= 0 ) stop
            endif
        enddo

        ! Ensure that all processes have managed to declare properly the rain flow walkers.
        call completion

        ! Update number of sources
        num_fw = k
        nsourceactive = sourceactive

        ! Deallocate some local arrays
        if( allocated( fce ) ) deallocate( fce, fcS )

        ! Transfer values globally
        call mpi_allreduce( num_fw, tot_fw, 1, int_type, sum_type, SPModel_comm_world, ierr )
        call mpi_allreduce( nsourceactive, sourceactive, 1, int_type, max_type, SPModel_comm_world, ierr )

        return

    end subroutine define_rain_flow_walker
    ! ============================================================================

end module rain_initial
