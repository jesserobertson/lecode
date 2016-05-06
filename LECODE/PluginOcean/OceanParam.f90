module Oceanplug

    use hdf5
    use file_data
    use mpidata
    use strata_data
    use TIN_data
    use sediment_data
    use forces_data
    use param_data
    use flux_data
    use time_data
    use ocean_data
    use mod_diffuse

    implicit none

    real(tkind), dimension( : ),allocatable :: points

contains

    ! ============================================================================
    !> Subroutine create_oceanplug_surface
    !! Create the ocean surface for the Ocean plugin.
    !<
    ! ============================================================================
    subroutine create_oceanplug_surface

        integer :: iunit, ios, n

        ! Open a file tunnel
        iunit = 20
        open(iunit,file=foceanxyz,status="replace",action="write",iostat=ios)
        rewind(iunit)

        ! Write topographic file relative to sea level
        do n = 1, G_strat%nodes_nb
            write(iunit, *) gv_coord( n, 1 ), gv_coord( n, 2 ), elev_record( n ) - gsea%actual_sea
        enddo
        close( iunit )

        ! Update/Create the poseidon file
        open(iunit,file=poseidon,status="replace",action="write",iostat=ios)
        rewind(iunit)

        write(iunit,'(a1)')'O'
        write(iunit,'(a1)')' '

        close( iunit )

        ! Allocate total load transport locally
        if( .not. allocated( soulsby_total_load ) ) allocate( soulsby_total_load( L_strat%nodes_nb, totgrn ) )

        if( .not. allocated( total_load ) ) allocate( total_load( DL_nodesNb, totgrn ) )
        if( .not. allocated( transportX ) ) allocate( transportX( DL_nodesNb, totgrn ) )
        if( .not. allocated( transportY ) ) allocate( transportY( DL_nodesNb, totgrn ) )

        return

    end subroutine create_oceanplug_surface
    ! ============================================================================
    !> Subroutine ocean_plugin_wait_function
    !! Wait for Ocean to compute the waves & currents fields.
    !<
    ! ============================================================================
    subroutine ocean_plugin_wait_function

        integer :: iunit, ios, scn, scn1, sgp, err

        character(len=1) :: charac
        character(len=128) :: stg

        charac = 'O'

        if( iam == 0 )then
112     continue
            do while( charac /= 'L' )
                iunit = 12
                ! Read the maestro file
                open(iunit,file=poseidon,status="old",action="read",iostat=ios)
                rewind(iunit)
                ! Read character
                read(iunit,'(a1)',err=112,end=112,iostat=err) charac
                if( err /= 0 )  charac = 'O'
                close( iunit )
                call Sleep( 1 )
           enddo
       endif
       call mpi_barrier( SPModel_comm_world, ierr )

       ! Find the appropriate climatic forces
       do scn = 1, forecast_nb
            if( tnow >= hindcast( scn )%tstart .and. tnow < hindcast( scn )%tend )then
                scn1 = scn
            endif
       enddo

       ! Define ocean arrays
       if( allocated( ocean_current ) ) deallocate( ocean_current )
       if( allocated( ocean_wave ) ) deallocate( ocean_wave )
       allocate( ocean_current( hindcast( scn1 )%cnb, G_strat%nodes_nb, 2 ) )
       allocate( ocean_wave( hindcast( scn1 )%cnb, G_strat%nodes_nb, 2 ) )
       active_hindcast = scn1

       do sgp = 1, hindcast( scn1 )%cnb

            foceanout = 'ocean_scn'
            call noblnk( foceanout )
            call noblnk( foceanout )
            call append_nb( foceanout,sgp )
            stg = '.h5'
            call append_str( foceanout,stg )
            call noblnk(foceanout)
            call addpath4( foceanout )

            ! Read waves & currents fields
            call read_waves_currents_values( sgp )

       enddo

       return

    end subroutine ocean_plugin_wait_function
    ! ============================================================================
    !> Subroutine read_waves_currents_values
    !! Reads waves & currents fields.
    !<
    ! ============================================================================
    subroutine read_waves_currents_values( sgp )

        logical :: compression, simple

        character(len=128) :: text

        integer :: hdferr, filter, filter_all, rank, n, sgp, p

        integer(hid_t) :: file_id, d_spc
        integer(hid_t) :: dset_id, dtype_id
        integer(hsize_t),dimension( 2 ) :: dims, maxdims

        real( tkind ) :: water_depth

        ! Initialize predefined datatypes
        call h5open_f( hdferr )
        call h5zfilter_avail_f( h5z_filter_deflate_f, compression, hdferr )
        call h5zget_filter_info_f( h5z_filter_deflate_f, filter, hdferr )

        filter_all = ior( h5z_filter_encode_enabled_f, h5z_filter_decode_enabled_f )
        if( filter_all .ne. filter ) compression = .false.

        ! Open the file collectively.
        call h5fopen_f( foceanout, h5f_acc_rdonly_f, file_id, hdferr )

        text = ''
        text = "/WavesCirculation"
        call h5dopen_f( file_id, trim(text), dset_id, hdferr )
        call h5dget_type_f( dset_id, dtype_id, hdferr )
        call h5dget_space_f( dset_id, d_spc, hdferr )
        call h5sis_simple_f( d_spc, simple, hdferr )
        call h5sget_simple_extent_ndims_f( d_spc, rank, hdferr )
        call h5sget_simple_extent_dims_f( d_spc, dims, maxdims, hdferr )
        if( .not. allocated( points ) ) allocate( points( dims( 1 ) * dims( 2 ) ) )
        call h5dread_f( dset_id, h5t_native_double, points, dims, hdferr )

        n = 0
        do p = 1, G_strat%nodes_nb
            if( current_on .and. .not. wave_on )then
                ocean_current( sgp, p, 1 ) = points( n + 1 )
                ocean_current( sgp, p, 2 ) = points( n + 2 )
                ocean_wave( sgp, p, 1:2 ) = 0.0_8
                n = n + 2
            elseif( wave_on .and. .not. current_on )then
                ocean_current( sgp, p, 1:2 ) = 0.0_8
                ocean_wave( sgp, p, 1 ) = points( n + 1 )
                ocean_wave( sgp, p, 2 ) = points( n + 2 )
                n = n + 2
            elseif( wave_on .and. current_on )then
                ocean_current( sgp, p, 1 ) = points( n + 1 )
                ocean_current( sgp, p, 2 ) = points( n + 2 )
                ocean_wave( sgp, p, 1 ) = points( n + 3 )
                ocean_wave( sgp, p, 2 ) = points( n + 4 )
                n = n + 4
            endif
            water_depth = gsea%actual_sea - elev_record( p )
            if( water_depth < tor .or. water_depth > max_WC_action_depth )then
                ocean_current( sgp, p, 1 ) = 0.0_8
                ocean_current( sgp, p, 2 ) = 0.0_8
                ocean_wave( sgp, p, 1 ) = 0.0_8
                ocean_wave( sgp, p, 2 ) = 0.0_8
            endif

        enddo

        ! Close the dataset
        call h5sclose_f( d_spc, hdferr )
        call h5dclose_f( dset_id, hdferr )
        ! Close the file.
        call h5fclose_f( file_id, hdferr )
        ! Close interface
        call h5close_f( hdferr )

        ! Update border
        do p = 1, G_strat%nodes_nb

            ! West
            if( gv_coord( p, 1 ) == strat_xo )then
                ocean_current( sgp, p, 1 ) = ocean_current( sgp, p+1, 1 )
                ocean_current( sgp, p, 2 ) = ocean_current( sgp, p+1, 2 )
                ocean_wave( sgp, p, 1 ) = ocean_wave( sgp, p+1, 1 )
                ocean_wave( sgp, p, 2 ) = ocean_wave( sgp, p+1, 2 )
            endif

            ! East
            if( gv_coord( p, 1 ) == gv_coord( G_strat%nodes_nb, 1 ) )then
                ocean_current( sgp, p, 1 ) = ocean_current( sgp, p-1, 1 )
                ocean_current( sgp, p, 2 ) = ocean_current( sgp, p-1, 2 )
                ocean_wave( sgp, p, 1 ) = ocean_wave( sgp, p-1, 1 )
                ocean_wave( sgp, p, 2 ) = ocean_wave( sgp, p-1, 2 )
            endif

            ! North
            if( gv_coord( p, 2 ) == gv_coord( G_strat%nodes_nb, 2 )  )then
                ocean_current( sgp, p, 1 ) = ocean_current( sgp, p-strat_X-1, 1 )
                ocean_current( sgp, p, 2 ) = ocean_current( sgp, p-strat_X-1, 2 )
                ocean_wave( sgp, p, 1 ) = ocean_wave( sgp, p-strat_X-1, 1 )
                ocean_wave( sgp, p, 2 ) = ocean_wave( sgp, p-strat_X-1, 2 )
            endif

            ! South
            if( gv_coord( p, 2 ) == strat_yo )then
                ocean_current( sgp, p, 1 ) = ocean_current( sgp, p+strat_X+1, 1 )
                ocean_current( sgp, p, 2 ) = ocean_current( sgp, p+strat_X+1, 2 )
                ocean_wave( sgp, p, 1 ) = ocean_wave( sgp, p+strat_X+1, 1 )
                ocean_wave( sgp, p, 2 ) = ocean_wave( sgp, p+strat_X+1, 2 )
            endif

        enddo

        ! Update corners
        ! SW
        p = 1
        ocean_current( sgp, p, 1 ) = ocean_current( sgp, p+strat_X+1, 1 )
        ocean_current( sgp, p, 2 ) = ocean_current( sgp, p+strat_X+1, 2 )
        ocean_wave( sgp, p, 1 ) = ocean_wave( sgp, p+strat_X+1, 1 )
        ocean_wave( sgp, p, 2 ) = ocean_wave( sgp, p+strat_X+1, 2 )
        ! SE
        p = strat_X
        ocean_current( sgp, p, 1 ) = ocean_current( sgp, p+strat_X-1, 1 )
        ocean_current( sgp, p, 2 ) = ocean_current( sgp, p+strat_X-1, 2 )
        ocean_wave( sgp, p, 1 ) = ocean_wave( sgp, p+strat_X-1, 1 )
        ocean_wave( sgp, p, 2 ) = ocean_wave( sgp, p+strat_X-1, 2 )
        ! NW
        p = strat_X * ( strat_Y - 1 ) + 1
        ocean_current( sgp, p, 1 ) = ocean_current( sgp, p-strat_X+1, 1 )
        ocean_current( sgp, p, 2 ) = ocean_current( sgp, p-strat_X+1, 2 )
        ocean_wave( sgp, p, 1 ) = ocean_wave( sgp, p-strat_X+1, 1 )
        ocean_wave( sgp, p, 2 ) = ocean_wave( sgp, p-strat_X+1, 2 )
        ! NE
        p = strat_X * strat_Y
        ocean_current( sgp, p, 1 ) = ocean_current( sgp, p-strat_X-1, 1 )
        ocean_current( sgp, p, 2 ) = ocean_current( sgp, p-strat_X-1, 2 )
        ocean_wave( sgp, p, 1 ) = ocean_wave( sgp, p-strat_X-1, 1 )
        ocean_wave( sgp, p, 2 ) = ocean_wave( sgp, p-strat_X-1, 2 )

        return

    end subroutine read_waves_currents_values
    ! ============================================================================

end module Oceanplug
