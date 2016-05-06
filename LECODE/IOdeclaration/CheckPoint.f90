! ============================================================================
! Name        : CheckPoint.f90
! Author      : Tristan Salles
! Created on: Jan 15, 2013
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file CheckPoint.f90
!
! Description : CheckPoint is used to write the information usefull for continuing a SPModel
! experiment.
!
!<
! ============================================================================
module check_point

    use hdf5
    use file_data
    use mpidata
    use time_data
    use FoX_wxml
    use strata_data
    use sediment_data

    implicit none

    public

    integer, dimension( : ), allocatable :: checkgid

contains

    ! ============================================================================
    !> Subroutine checkpt_layers
    !! Subroutine checkpt_layers gather sediment layers information to the processors.
    !<
    ! ============================================================================
    subroutine checkpt_layers

        integer :: i, k, nid

        allocate( rec_NbLay( L_strat%nodes_nb ) )
        allocate( rec_topelev( L_strat%nodes_nb ) )
        allocate( rec_coord( L_strat%nodes_nb, 3 ) )

        allocate( rec_slayID( L_strat%nodes_nb, L_strat%layersID ) )
        allocate( rec_shardness( L_strat%nodes_nb, L_strat%layersID ) )
        allocate( rec_szelev( L_strat%nodes_nb, L_strat%layersID ) )
        allocate( rec_sporosity( L_strat%nodes_nb, L_strat%layersID, totgrn ) )
        allocate( rec_sthick( L_strat%nodes_nb, L_strat%layersID ) )
        allocate( rec_sfac( L_strat%nodes_nb, L_strat%layersID ) )
        allocate( rec_ssedh( L_strat%nodes_nb, L_strat%layersID, totgrn ) )

        allocate( checkgid( L_strat%nodes_nb ) )
        rec_NbLay( : ) = -1
        do i = 1, L_strat%nodes_nb
            nid = locv_gid( i )
            if( nid > 0 .and. rec_NbLay( i ) == -1 )then
                checkgid( i ) = nid
                rec_NbLay( i ) = locv_NbLay( i )
                rec_coord( i, 1:2 ) = gv_coord( nid, 1:2 )
                rec_coord( i, 3 ) = locv_coord( i, 3 )
                rec_topelev( i ) = elev_record( nid )
                do k = 1, rec_NbLay( i )
                    rec_ssedh( i , k , 1:totgrn ) = 0.0_8
                    rec_slayID( i , k ) = strat_layID( i, k )
                    rec_shardness( i , k ) = strat_hardness(  i ,  k  )
                    rec_szelev( i, k ) = strat_zelev(  i ,  k  )
                    rec_sporosity( i , k, 1:totgrn ) = strat_porosity(  i ,  k, 1:totgrn  )
                    rec_sthick( i , k ) = strat_thick(  i ,  k  )
                    rec_sfac( i , k ) = strat_fac(  i ,  k  )
                    rec_ssedh( i , k , 1:totgrn ) = strat_sedh(  i ,  k ,  1:totgrn  )
                enddo
            endif
        enddo

        return

    end subroutine checkpt_layers
    ! ============================================================================
    !> Subroutine checkpointer
    !!
    !! Subroutine checkpointer records layers space evolution through time.
    !! \param iter
    !<
    ! ============================================================================
    subroutine checkpointer( iter )

        integer, intent( in ) :: iter

        logical :: compression

        integer :: k, p, q, n, rank
        integer :: hdferr
        real(tkind), dimension( : ),allocatable :: nodes, lays

        integer(hid_t) :: file_id, plist_id
        integer(hid_t) :: filespace, dset_id
        integer(hsize_t),dimension(2) :: dims

        character(len=128) :: text, checkf, stg

        checkf = ''
        checkf = 'checkpoint'
        call noblnk( checkf )
        text = '.'
        call append_str( checkf,text )
        call append_zero( checkf,iter )
        stg = '.p'
        call append_str( checkf,stg )
        call append_nb( checkf,iam )
        stg = '.h5'
        call append_str( checkf,stg )
        call noblnk(checkf)
        call addpath2( checkf )

        ! Gather all information from the processors
        call checkpt_layers

        ! Initialize predefined datatypes
        call h5open_f( hdferr )
        call h5zfilter_avail_f( h5z_filter_deflate_f, compression, hdferr )

        ! Setup file access property list for MPI-IO access.
        call h5pcreate_f( h5p_file_access_f, plist_id, hdferr )
        ! Create the file collectively.
        call h5fcreate_f( checkf, h5f_acc_excl_f, file_id, hdferr )

        ! ========================
        ! Stratal coordinates
        ! ========================
        dims( 1 ) = 6
        dims( 2 ) = L_strat%nodes_nb
        rank = 2
        ! Create dataspace and opens it for access
        call h5screate_simple_f( rank, dims, filespace, hdferr )
        allocate( nodes( 6 * dims( 2 ) ) )
        n = 0
        do p = 1, L_strat%nodes_nb
            nodes( n + 1 ) = real( checkgid( p ) )
            nodes( n + 2 ) = rec_coord( p , 1 )
            nodes( n + 3 ) = rec_coord( p , 2 )
            nodes( n + 4 ) = rec_coord( p , 3 )
            nodes( n + 5 ) = rec_topelev( p )
            nodes( n + 6 ) = rec_NbLay( p )
            n = n + 6
        enddo
        text = ''
        text = "/DataVertices"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, hdferr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, hdferr )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, nodes, dims, hdferr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, hdferr )
            call h5pset_deflate_f( plist_id, 9, hdferr )
            call h5pset_chunk_f( plist_id, rank, dims, hdferr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, hdferr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, nodes, dims, hdferr )
            call h5pclose_f( plist_id, hdferr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, hdferr )
        call h5sclose_f( filespace, hdferr )

        ! ========================
        ! Deposition space parameters
        ! ========================
        dims( 1 ) = 5 + 2*totgrn
        dims( 2 ) = L_strat%nodes_nb * L_strat%layersID
        rank = 2
        ! Create dataspace and opens it for access
        call h5screate_simple_f( rank, dims, filespace, hdferr )
        allocate( lays( dims( 1 ) * dims( 2 ) ) )
        n = 0
        do p = 1, L_strat%nodes_nb
            do k = 1, rec_NbLay( p )
                n = n + 1
                lays( n ) = rec_slayID( p , k )
                n = n + 1
                lays( n ) = rec_shardness( p , k )
                n = n + 1
                lays( n ) = rec_szelev( p, k )
                n = n + 1
                lays( n ) = rec_sthick( p , k )
                n = n + 1
                lays( n ) = rec_sfac( p , k )
                do q = 1, totgrn
                    n = n + 1
                    lays( n ) = rec_ssedh( p , k , q )
                enddo
                do q = 1, totgrn
                    n = n + 1
                    lays( n ) = rec_sporosity( p , k , q )
                enddo
            enddo
            if( rec_NbLay( p ) + 1 <=  L_strat%layersID )then
                do k = rec_NbLay( p ) + 1, L_strat%layersID
                    n = n + 1
                    lays( n ) = 0.0_8
                    n = n + 1
                    lays( n ) = 0.0_8
                    n = n + 1
                    lays( n ) = 0.0_8
                    n = n + 1
                    lays( n ) = 0.0_8
                    n = n + 1
                    lays( n ) = 0.0_8
                    do q = 1, 2*totgrn
                        n = n + 1
                        lays( n ) = 0.0_8
                    enddo
                enddo
            endif
        enddo
        text = ''
        text = "/DepositsLayer"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, hdferr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, hdferr )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, lays, dims, hdferr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, hdferr )
            call h5pset_deflate_f( plist_id, 9, hdferr )
            call h5pset_chunk_f( plist_id, rank, dims, hdferr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, hdferr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, lays, dims, hdferr )
            call h5pclose_f( plist_id, hdferr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, hdferr )
        call h5sclose_f( filespace, hdferr )

        ! Close the file.
        call h5fclose_f( file_id, hdferr )
        ! Close interface
        call h5close_f( hdferr )

        deallocate( lays, nodes, checkgid )

        return

    end subroutine checkpointer
    ! ============================================================================
    !> Subroutine define_restarted_simulation
    !! Subroutine define_restarted_simulation reads and creates a new simulation based on previous run.
    !<
    ! ============================================================================
    subroutine define_restarted_simulation

        logical :: simple, compression

        integer :: k, p, q, n, rank, m, nbnodes
        integer :: hdferr, filter, filter_all, maxlay
        integer, dimension( : ),allocatable :: lgid

        real(tkind), dimension( : ),allocatable :: nodes, lays
        real(tkind) :: h

        integer(hid_t) :: file_id, d_spc
        integer(hid_t) :: dset_id, dtype_id
        integer(hsize_t),dimension( 2 ) :: dims, maxdims

        character(len=128) :: text, stg

        do m = 1, restart_proc

            ! Previous files for restart
            frestart = restartfolder
            call noblnk( frestart )
            stg = '/runfiles/checkpoint.'
            call append_str( frestart, stg )
            call append_zero( frestart, restart_iter )
            stg = '.p'
            call append_str( frestart,stg )
            call append_nb( frestart,m-1 )
            stg = '.h5'
            call append_str( frestart, stg )
            call noblnk( frestart )

            ! Initialize predefined datatypes
            call h5open_f( hdferr )

            call h5zfilter_avail_f( h5z_filter_deflate_f, compression, hdferr )
            call h5zget_filter_info_f( h5z_filter_deflate_f, filter, hdferr )
            filter_all = ior( h5z_filter_encode_enabled_f, h5z_filter_decode_enabled_f )
            if( filter_all .ne. filter ) compression = .false.

            ! Open the file collectively.
            call h5fopen_f( frestart, h5f_acc_rdonly_f, file_id, hdferr )

            !-----------------------------------------------
            ! Layer coordinates
            text = ''
            text = "/DataVertices"
            call h5dopen_f( file_id, trim(text), dset_id, hdferr )
            call h5dget_type_f( dset_id, dtype_id, hdferr )
            call h5dget_space_f( dset_id, d_spc, hdferr )
            call h5sis_simple_f( d_spc, simple, hdferr )
            call h5sget_simple_extent_ndims_f( d_spc, rank, hdferr )
            call h5sget_simple_extent_dims_f( d_spc, dims, maxdims, hdferr )
            allocate( lgid( dims( 2 ) ) )
            allocate( nodes( dims( 1 ) * dims( 2 ) ) )
            call h5dread_f( dset_id, h5t_native_double, nodes, dims, hdferr )

            n = 0
            nbnodes = int( dims( 2 ), 4 )
            do k = 1, nbnodes
                lgid( k ) =  int( nodes( n + 1 ) )
                rec_coord( int( nodes( n + 1 ) ), 1 ) = nodes( n + 2 )
                rec_coord( int( nodes( n + 1 ) ), 2 ) = nodes( n + 3 )
                rec_coord( int( nodes( n + 1 ) ), 3 ) = nodes( n + 4 )
                rec_topelev( int( nodes( n + 1 ) ) ) = nodes( n + 5 )
                rec_NbLay( int( nodes( n + 1 ) ) ) = int( nodes( n + 6 ) )
                n = n + 6
            enddo
            maxlay = L_strat%layersID
            if( .not. allocated( rec_slayID ) ) allocate( rec_slayID( G_strat%nodes_nb, maxlay ) )
            if( .not. allocated( rec_shardness ) ) allocate( rec_shardness( G_strat%nodes_nb, maxlay ) )
            if( .not. allocated( rec_szelev ) ) allocate( rec_szelev( G_strat%nodes_nb, maxlay ) )
            if( .not. allocated( rec_sporosity ) ) allocate( rec_sporosity( G_strat%nodes_nb, maxlay, totgrn ) )
            if( .not. allocated( rec_sthick ) ) allocate( rec_sthick( G_strat%nodes_nb, maxlay ) )
            if( .not. allocated( rec_sfac ) ) allocate( rec_sfac( G_strat%nodes_nb, maxlay ) )
            if( .not. allocated( rec_ssedh ) ) allocate( rec_ssedh( G_strat%nodes_nb, maxlay, totgrn ) )

            ! Close the dataset
            call h5sclose_f( d_spc, hdferr )
            call h5dclose_f( dset_id, hdferr )

            !-----------------------------------------------
            ! Deposition space parameters
            text = ''
            text = "/DepositsLayer"
            call h5dopen_f( file_id, trim(text), dset_id, hdferr )
            call h5dget_type_f( dset_id, dtype_id, hdferr )
            call h5dget_space_f( dset_id, d_spc, hdferr )
            call h5sis_simple_f( d_spc, simple, hdferr )
            call h5sget_simple_extent_ndims_f( d_spc, rank, hdferr )
            call h5sget_simple_extent_dims_f( d_spc, dims, maxdims, hdferr )
            allocate( lays( dims( 1 ) * dims( 2 ) ) )
            call h5dread_f( dset_id, h5t_native_double, lays, dims, hdferr )
            n = 0
            do k = 1, nbnodes
                do p = 1, rec_NbLay( lgid( k ) )
                    n = n + 1
                    rec_slayID( lgid( k ), p ) = int( lays( n ) )
                    n = n + 1
                    rec_shardness( lgid( k ), p ) = lays( n )
                    n = n + 1
                    rec_szelev( lgid( k ), p ) = lays( n )
                    n = n + 1
                    rec_sthick( lgid( k ), p ) = lays( n )
                    n = n + 1
                    rec_sfac( lgid( k ), p ) = lays( n )
                    h = 0.0_8
                    do q = 1, totgrn
                        n = n + 1
                        rec_ssedh( lgid( k ), p, q ) = lays( n )
                        h = dble( h + rec_ssedh( lgid( k ), p, q ) )
                    enddo
                    do q = 1, totgrn
                        n = n + 1
                        rec_sporosity( lgid( k ), p, q ) = lays( n )
                    enddo
                    if( h /= rec_sthick( lgid( k ), p ) .and. h /= 0.0_8 )then
                        do q = 1, totgrn
                            rec_ssedh( lgid( k ), p, q ) = rec_ssedh( lgid( k ), p, q ) * &
                                rec_sthick( lgid( k ), p ) / h
                        enddo
                    endif
                    h = 0.0_8
                    do q = 1, totgrn
                        h = dble( h + rec_ssedh( lgid( k ), p, q ) )
                    enddo
                    rec_sthick( lgid( k ), p ) = h
                enddo
                if( rec_NbLay( lgid( k ) ) + 1 <=  L_strat%layersID )then
                    do p = rec_NbLay( lgid( k ) ) + 1, L_strat%layersID
                        n = n + 5 + 2*totgrn
                    enddo
                endif
            enddo

            ! Close the dataset
            call h5sclose_f( d_spc, hdferr )
            call h5dclose_f( dset_id, hdferr )
            ! Close the file.
            call h5fclose_f( file_id, hdferr )
            ! Close interface
            call h5close_f( hdferr )
            deallocate( nodes, lays, lgid )

        enddo

        do k = 1, G_strat%nodes_nb
            gv_coord( k, 1:3 ) = rec_coord( k, 1:3 )
            elev_record( k ) = rec_topelev( k )
        enddo

        return

    end subroutine define_restarted_simulation
    ! ============================================================================

end module check_point
