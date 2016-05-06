! ============================================================================
! Name        : StrataO.f90
! Author      : Tristan Salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file StrataO.f90
!
! Description : StrataO is used to generate Hdf5 output of the stratal mesh during SPModel simulation.
!
!<
! ============================================================================
module strata_out

    use hdf5
    use SL_out
    use TIN_out
    use mpidata
    use file_data
    use time_data
    use FoX_wxml
    use FWLK_out
    use strata_data
    use sediment_data

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine xdmf_output
    !! Subroutine xdmf_output controls SP Model output frequency.
    !! \param iter
    !<
    ! ============================================================================
    subroutine xdmf_output( iter )

        integer :: iter

        ! TIN
        if( iam == 0 )then
            call TIN_hdf5( iter )
            call TIN_xmf( iter )
            call TIN_series( iter )
            if( gsea%output )then
                call SL_xmf( iter )
                call SL_series( iter )
            endif
        endif
        ! Stratigraphy
        call Strata_hdf5( iter, gporo%compaction )
        call Strata_xmf( iter, gporo%compaction )
        if( iam == 0 ) call Strata_series( iter )
        ! flow walkers
        call FWLK_hdf5( iter )
        call FW_xmf( iter )
        if( iam == 0 ) call FW_series( iter )

        call mpi_barrier( SPModel_comm_world,ierr )

    end subroutine xdmf_output
    ! ============================================================================
    !> Subroutine Strata_hdf5
    !! Subroutine Strata_hdf5 uses Hdf5 to output SP Model stratal evolution.
    !! It records both connectivity and vertices.
    !<
    ! ============================================================================
    subroutine Strata_hdf5( iter, compac )

        logical :: compac, compression

        ! Parameters Declaration
        integer :: totnodes, totelems, id, id0, i, rank, iter, k, pid, idt, ks, gid, c, nnb, idt2
        integer, dimension( : ),allocatable :: layid, nlay, connect

        real(tkind) :: mgz
        real(tkind), dimension( : ),allocatable :: nodes, thick, poro, mz, fac, dis, Ethick, Emz
        real(tkind), dimension( :,: ),allocatable :: prop, Eprop

        character(len=128) :: text, stg, file

        integer(hid_t) :: file_id, plist_id
        integer(hid_t) :: filespace, dset_id
        integer(hsize_t),dimension(2) :: dims

        file = ''
        file = 'SMesh'
        call noblnk( file )
        stg = '.'
        call append_str( file,stg )
        call append_zero( file,iter )
        stg = '.p'
        call append_str( file,stg )
        call append_nb( file,iam )
        stg = '.h5'
        call append_str( file,stg )
        call addpath1( file )
        totnodes = L_strat%nodes_nb * ( L_strat%layersID  )
        totelems = L_strat%faces_nb * ( L_strat%layersID - 1 )
        allocate( fac(totnodes) )
        allocate( dis(totnodes) )
        allocate( mz(totnodes) )
        allocate( nlay(totnodes) )
        allocate( thick(totnodes) )
        allocate( poro(totnodes) )
        allocate( prop(totnodes,totgrn) )
        allocate( nodes(3*totnodes) )
        allocate( connect(8*totelems) )
        allocate( layid(L_strat%nodes_nb) )
        allocate( Emz(totelems) )
        allocate( Ethick(totelems) )
        allocate( Eprop(totelems,totgrn) )

        check_int = 0

        ! Create nodes arrays
        id = 1
        id0 = 1
        idt = 1
        idt2 = 1
        layid = 0
        mz = 0.0_8
        thick = 0.0_8
        poro = 0.0_8
        prop = 0.0_8
        Eprop = 0.0_8
        Ethick = 0.0_8

        do i = 1, L_strat%nodes_nb
            gid = locv_gid( i )
            check_int( id0 ) = gid
            check_int( id0 + 1 ) = locv_NbLay( i )
            nodes( id ) = real( locv_coord( i, 1 ) )
            nodes( id + 1 ) = real( locv_coord( i, 2 ) )
            nodes( id + 2 ) = real( locv_coord( i, 3 ) )
            thick( idt ) = 0.0_8
            poro( idt ) = 0.0_8
            prop( idt,1:totgrn ) = 0.0_8
            fac( idt ) = 0.0_8
            mz( idt ) = 0.0_8
            dis( idt ) = real( locv_vrate( i ) )
            nlay( idt ) = 1
            idt = idt + 1
            idt2 = idt2 + totgrn
            id0 = id0 + 2
            id = id + 3
        enddo

        do k = 1, L_strat%layersID - 1
            do i = 1, L_strat%nodes_nb
                gid = locv_gid( i )
                nodes( id ) = real( locv_coord( i, 1 ) )
                nodes( id + 1 ) = real( locv_coord( i, 2 ) )
                if( layid( i ) < locv_NbLay( i ) )then
                    if( strat_layID( i , layid( i ) + 1 ) == k )then
                        layid( i ) = layid( i ) + 1
                    endif
                endif
                if( strat_layID( i , layid( i ) ) /= k )then
                    pid = 3 * ( L_strat%nodes_nb * ( k - 1 ) + ( i  ) )
                    if( udw ) pid = pid - 1
                    nodes( id + 2 ) = nodes( pid )
                    dis( idt ) = real( locv_vrate( i ) )
                    mgz = 0.0_8
                    poro( idt ) = 0.0_8
                    if( strat_thick( i, layid( i ) ) > 0.0_8 )then
                        do ks = 1, totgrn
                            prop( idt, ks ) = real( strat_sedh( i, layid( i ), ks ) )
                            mgz = real( mgz + sediment( ks )%diameter *  strat_sedh( i, layid( i ), ks ) )
                            poro( idt ) =  poro( idt ) + real( strat_porosity( i, layid( i ), ks ) * &
                                prop( idt, ks )/strat_thick( i, layid( i ) )  )
                        enddo
                        mz( idt ) = real( mgz / strat_thick( i, layid( i ) ) )
                    else
                        mz( idt ) = real( mgz )
                    endif
                    fac( idt ) = real( strat_fac( i, layid( i ) ) )
                    thick( idt ) =  0.0_8
                    mz( idt ) = mz( idt ) * 1000.0_8
                else
                    nodes( id + 2 ) = real( strat_zelev( i, layid( i ) ) )
                    thick( idt ) =  real( strat_thick( i, layid( i ) ) )
                    fac( idt ) = real( strat_fac( i, layid( i ) ) )
                    dis( idt ) = real( locv_vrate( i ) )
                    mgz = 0.0_8
                    poro( idt ) = 0.0_8
                    if( strat_thick( i, layid( i ) ) > 0.0_8 )then
                        do ks = 1, totgrn
                            prop( idt, ks ) = real( strat_sedh( i, layid( i ), ks ) )
                            mgz = real( mgz + sediment( ks )%diameter * strat_sedh( i, layid( i ), ks ) )
                            poro( idt ) =  poro( idt ) + real( strat_porosity( i, layid( i ), ks ) * prop( idt, ks ) )
                        enddo
                        mz( idt ) = real( mgz / strat_thick( i, layid( i ) ) )
                        poro( idt ) = real( poro( idt ) / strat_thick( i, layid( i ) )  )
                    else
                        mz( idt ) = real( mgz )
                    endif
                    mz( idt ) = mz( idt ) * 1000.0_8
                endif
                id = id + 3
                nlay( idt ) = k + 1
                idt = idt + 1
                idt2 = idt2 + totgrn
            enddo
        enddo
        ! Create cell data
        id = 1
        idt = 1
        do i = 1, L_strat%faces_nb
            connect( id : id + 3 ) = locf_points( i, 1:4 )
            id = id + 4
            connect( id : id + 3 ) = locf_points( i, 1:4 ) + L_strat%nodes_nb
            Emz( idt ) = 0.0_8
            Eprop( idt, : ) = 0.0_8
            Ethick( idt ) = 0.0_8
            do c = 1, 4
                Ethick( idt ) = real( Ethick( idt ) + thick( connect( id - 1 + c ) ) )
                Emz( idt ) = real( Emz( idt ) + mz( connect( id - 1 + c ) ) )
                do ks = 1, totgrn
                    Eprop( idt, ks ) = real( Eprop( idt, ks ) + prop( connect( id - 1 + c ), ks ) )
                enddo
            enddo
            Emz( idt ) = Emz( idt ) * 0.25_8
            Ethick( idt ) = Ethick( idt ) * 0.25_8
            if( Ethick( idt ) > 0.0_8 )then
                do ks = 1, totgrn
                    Eprop( idt, ks ) = Eprop( idt, ks ) * 0.25_8
                    Eprop( idt, ks ) = Eprop( idt, ks ) / Ethick( idt )
                    Eprop( idt, ks ) = Eprop( idt, ks ) * 100.0_8
                    if( Eprop( idt, ks ) > 100.0_8 ) Eprop( idt, ks ) = 100.0_8
                enddo
            else
                Eprop( idt, : ) = 0.0_8
            endif
            idt = idt + 1
            id = id + 4
        enddo
        do k = 1, L_strat%layersID - 2
            do i = 1, L_strat%faces_nb
                Ethick( idt ) = 0.0_8
                connect( id : id + 3 ) = locf_points( i, 1:4 ) + k * L_strat%nodes_nb
                id = id + 4
                connect( id : id + 3 ) = locf_points( i, 1:4 ) + ( k + 1 ) * L_strat%nodes_nb
                nnb = 0
                Emz( idt ) = 0.0_8
                Eprop( idt, : ) = 0.0_8
                do c = 1, 4
                    Ethick( idt ) = real( Ethick( idt ) + thick( connect( id - 1 + c ) ) )
                    if( thick( connect( id - 1 + c ) ) > 0.0_8 )then
                        if( mz( connect( id - 1 + c ) ) > 0.0_8 ) nnb = nnb + 1
                        Emz( idt ) = real( Emz( idt ) + mz( connect( id - 1 + c ) ) )
                        do ks = 1, totgrn
                            Eprop( idt, ks ) = real( Eprop( idt, ks ) + prop( connect( id - 1 + c ), ks ) )
                        enddo
                    endif
                enddo
                Ethick( idt ) = Ethick( idt ) * 0.25_8
                if( nnb == 0 ) Emz( idt ) = 0.0_8
                if( nnb > 0 ) Emz( idt ) = Emz( idt ) / real( nnb )
                if( Ethick( idt ) > 0.0_8 )then
                    do ks = 1, totgrn
                        Eprop( idt, ks ) = Eprop( idt, ks ) * 0.25_8
                        Eprop( idt, ks ) = Eprop( idt, ks ) / Ethick( idt )
                        Eprop( idt, ks ) = Eprop( idt, ks ) * 100.0_8
                        if( Eprop( idt, ks ) > 100.0_8 ) Eprop( idt, ks ) = 100.0_8
                    enddo
                else
                    Eprop( idt, : ) = 0.0_8
                endif
                id = id + 4
                idt = idt + 1
            enddo
        enddo

        ! Initialize predefined datatypes
        call h5open_f( ierr )
        call h5zfilter_avail_f( h5z_filter_deflate_f, compression, ierr )

        ! Setup file access property list for MPI-IO access.
        call h5pcreate_f( h5p_file_access_f, plist_id, ierr )
        ! Create the file collectively.
        call h5fcreate_f( file, h5f_acc_trunc_f, file_id, ierr, access_prp = plist_id )

        ! ========================
        ! The Coordinates - vertices
        ! ========================
        dims( 1 ) = 3
        dims( 2 ) = totnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/vertices"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, nodes, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, nodes, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Thickness attribute
        ! ========================
        dims( 1 ) = 1
        dims( 2 ) = totnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/thick"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the dataset
            call h5dwrite_f( dset_id, h5t_native_double, thick, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the dataset
            call h5dwrite_f( dset_id, h5t_native_double, thick, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        if( compac )then
            ! ========================
            ! Porosity attribute
            ! ========================
            dims( 1 ) = 1
            dims( 2 ) = totnodes
            rank = 2
            call h5screate_simple_f( rank, dims, filespace, ierr )
            text = ''
            text = "/porosity"
            if( .not. compression )then
                ! Create the dataset with default properties
                call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
                ! Create property list for collective dataset write
                call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
                ! Write the dataset
                call h5dwrite_f( dset_id, h5t_native_double, poro, dims, ierr, &
                    file_space_id = filespace, xfer_prp = plist_id )
            else
                ! Create property list for collective dataset write
                call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
                call h5pset_deflate_f( plist_id, 9, ierr )
                call h5pset_chunk_f( plist_id, rank, dims, ierr )
                ! Create the dataset with default properties
                call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
                ! Write the dataset
                call h5dwrite_f( dset_id, h5t_native_double, poro, dims, ierr )
                call h5pclose_f( plist_id, ierr )
            endif
            ! Close the dataset
            call h5dclose_f( dset_id, ierr )
            call h5sclose_f( filespace, ierr )
        endif

        ! ========================
        ! Layer number type attribute
        ! ========================
        dims( 1 ) = 1
        dims( 2 ) = totnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/layerID"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_integer, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the dataset
            call h5dwrite_f( dset_id, h5t_native_integer, nlay, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_integer , filespace, dset_id, ierr, plist_id )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_integer, nlay, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Displacement attribute
        ! ========================
        dims( 1 ) = 1
        dims( 2 ) = totnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/vdisp"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, dis, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, dis, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Mean grainsize attribute
        ! ========================
        dims( 1 ) = 1
        dims( 2 ) = totnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/mgz"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, mz, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, mz, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Facies attribute
        ! ========================
        dims( 1 ) = 1
        dims( 2 ) = totnodes
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/fac"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, fac, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, fac, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Element thickness
        ! ========================
        dims( 1 ) = 1
        dims( 2 ) = totelems
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/ethick"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, Ethick, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, Ethick, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Element mean grainsize
        ! ========================
        dims( 1 ) = 1
        dims( 2 ) = totelems
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        text = ''
        text = "/emgz"
        if( .not. compression )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, Emz, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, Emz, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Element proportion
        ! ========================
        do ks = 1, totgrn
            dims( 1 ) = 1
            dims( 2 ) = totelems
            rank = 2
            call h5screate_simple_f( rank, dims, filespace, ierr )
            text = ''
            text = "/eperc"
            call append_nb( text,ks )
            if( .not. compression )then
                ! Create the dataset with default properties
                call h5dcreate_f( file_id, trim(text), h5t_native_double, filespace, dset_id, ierr )
                ! Create property list for collective dataset write
                call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
                ! Write the datase
                call h5dwrite_f( dset_id, h5t_native_double, Eprop( 1:totelems, ks ), dims, ierr, &
                    file_space_id = filespace, xfer_prp = plist_id )
            else
                ! Create property list for collective dataset write
                call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
                call h5pset_deflate_f( plist_id, 9, ierr )
                call h5pset_chunk_f( plist_id, rank, dims, ierr )
                ! Create the dataset with default properties
                call h5dcreate_f( file_id, trim(text), h5t_native_double , filespace, dset_id, ierr, plist_id )
                ! Write the datase
                call h5dwrite_f( dset_id, h5t_native_double, Eprop( 1:totelems, ks ), dims, ierr )
                call h5pclose_f( plist_id, ierr )
            endif
            ! Close the dataset
            call h5dclose_f( dset_id, ierr )
            call h5sclose_f( filespace, ierr )
        enddo

        ! Close the file.
        call h5fclose_f( file_id, ierr )
        ! Close interface
        call h5close_f( ierr )

        deallocate( nodes, layid, thick, nlay, poro, fac, mz, prop, dis, Emz, Eprop, Ethick, connect )

        return

    end subroutine Strata_hdf5
    ! ============================================================================
    !> Subroutine Strata_xmf
    !! Subroutine Strata_xmf generates the Strata XdmF file for the considered time simulation.
    !! It builds the file which calls all the hdf5 created.
    !<
    ! ============================================================================
    subroutine Strata_xmf( iter, compac )

        ! Parameters Declaration
        type(xmlf_t) :: xf

        logical :: compac
        integer :: iter, totnodes, k, ks, nX, nY, nbZ
        integer,dimension( nproc ) :: nbX, nbY

        character(len=128) :: str, filename, file, h5file, stg
        character(len=128) :: filename1, filename2, filename3, filename4, filename5
        character(len=128) :: filename6, filename7, filename8, filename9, filename10, filename11, txt

        totnodes = L_strat%nodes_nb * ( L_strat%layersID )

        ! Define decomposition for AEDyC SP component
        nX = int( ( locv_coord( L_strat%nodes_nb, 1 ) - locv_coord( 1, 1 ) ) / strat_dx + 1 )
        nY = int( ( locv_coord( L_strat%nodes_nb, 2 ) - locv_coord( 1, 2 ) ) / strat_dx + 1 )
        nbZ = L_strat%layersID
        call mpi_gather( nX,1,int_type, nbX,1,int_type,0,SPModel_comm_world,ierr)
        call mpi_gather( nY,1,int_type, nbY,1,int_type,0,SPModel_comm_world,ierr)

        if( iam == 0 )then
            file = ''
            file = 'Mesh'
            call noblnk( file )
            str = '.'
            call append_str( file,str )
            call append_zero( file,iter )
            str = '.xmf'
            call append_str( file,str )
            call addpath1( file )
            call xml_OpenFile( file, xf )
            ! Header
            call xml_AddDOCTYPE(xf, "Xdmf", "Xdmf.dtd")
            call xml_DeclareNamespace(xf, "http://www.w3.org/2001/XInclude", "xi")
            call xml_NewElement(xf, "Xdmf" )
            call xml_AddAttribute(xf, "Version", "2.0" )
            call xml_NewElement( xf, "Domain" )
            call xml_NewElement( xf, "Grid" )
            call xml_AddAttribute( xf, "GridType", "Collection" )
            call xml_AddAttribute( xf, "CollectionType", "Spatial" )
            call xml_NewElement( xf, "Time" )
            call xml_AddAttribute( xf, "Type", "Single" )
            call xml_AddAttribute( xf, "Value", next_display-time_display )
            call xml_EndElement( xf, "Time" )
            ! Loop over processors
            do k = 1, nproc
                h5file = ''
                h5file = 'SMesh'
                call noblnk( h5file )
                stg = '.'
                call append_str( h5file,stg )
                call append_zero( h5file,iter )
                stg = '.p'
                call append_str( h5file,stg )
                call append_nb( h5file,k-1 )
                stg = '.h5'
                call append_str( h5file,stg )
                filename = h5file
                call noblnk( filename )
                filename1 = filename
                filename2 = filename
                filename3 = filename
                filename4 = filename
                filename5 = filename
                filename6 = filename
                filename7 = filename
                filename8 = filename
                filename9 = filename
                filename10 = filename
                str = ':/vertices'
                call append_str( filename1,str )
                str = ':/thick'
                call append_str( filename2,str )
                str = ':/layerID'
                call append_str( filename3,str )
                str = ':/porosity'
                call append_str( filename4,str )
                str = ':/vdisp'
                call append_str( filename5,str )
                str = ':/fac'
                call append_str( filename6,str )
                str = ':/mgz'
                call append_str( filename7,str )
                str = ':/ethick'
                call append_str( filename8,str )
                str = ':/emgz'
                call append_str( filename9,str )
                str = ':/eperc'
                call append_str( filename10,str )
                ! Block begin
                call xml_NewElement( xf, "Grid" )
                str = 'SBlock.'
                call append_zero( str,iter )
                stg = '.p'
                call append_str( str,stg )
                call append_nb( str,k - 1 )
                call xml_AddAttribute( xf, "Name", trim( str ) )
                call xml_NewElement( xf, "Topology" )
                call xml_AddAttribute( xf, "TopologyType", "3DSMesh" )
                str = ' '
                call append_nb2( str, nbZ )
                call append_nb2( str, nbY( k ) )
                call append_nb2( str, nbX( k ) )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_EndElement( xf, "Topology" )
                ! Geometry
                call xml_NewElement( xf, "Geometry" )
                call xml_AddAttribute( xf, "Type", "XYZ" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                str = ' '
                call append_nb2( str, nbX( k ) )
                call append_nb2( str, nbY( k ) )
                call append_nb2( str, nbZ )
                call append_nb2( str, 3 )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename1 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Geometry" )
                ! Layer thickness
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Node" )
                call xml_AddAttribute( xf, "Name", "Layer Thickness" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                str = ' '
                call append_nb2( str, nbX( k ) )
                call append_nb2( str, nbY( k ) )
                call append_nb2( str, nbZ )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename2 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                ! Layer number
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Node" )
                call xml_AddAttribute( xf, "Name", "Layer Index" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Int" )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename3 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                ! Porosity
                if( compac )then
                    call xml_NewElement( xf, "Attribute" )
                    call xml_AddAttribute( xf, "Type", "Scalar" )
                    call xml_AddAttribute( xf, "Center", "Node" )
                    call xml_AddAttribute( xf, "Name", "Porosity" )
                    call xml_NewElement( xf, "DataItem" )
                    call xml_AddAttribute( xf, "Format", "HDF" )
                    call xml_AddAttribute( xf, "NumberType", "Float" )
                    call xml_AddAttribute( xf, "Precision", "4" )
                    call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                    call xml_AddCharacters( xf, trim( filename4 ) )
                    call xml_EndElement( xf, "DataItem" )
                    call xml_EndElement( xf, "Attribute" )
                endif
                ! Vertical displacement Rate
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Node" )
                call xml_AddAttribute( xf, "Name", "Vertical Displacement Rate" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename5 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                ! Facies
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Node" )
                call xml_AddAttribute( xf, "Name", "Relative Elevation" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename6 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                ! Mean grain size
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Node" )
                call xml_AddAttribute( xf, "Name", "Mean Grainsize" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename7 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                ! Cell thickness
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Cell" )
                call xml_AddAttribute( xf, "Name", "Cell Thickness" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                str = ' '
                call append_nb2( str, nbX( k ) - 1 )
                call append_nb2( str, nbY( k ) - 1 )
                call append_nb2( str, nbZ  - 1 )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename8 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                ! Cell mean grainsize
                call xml_NewElement( xf, "Attribute" )
                call xml_AddAttribute( xf, "Type", "Scalar" )
                call xml_AddAttribute( xf, "Center", "Cell" )
                call xml_AddAttribute( xf, "Name", "Cell Mean Grain Size" )
                call xml_NewElement( xf, "DataItem" )
                call xml_AddAttribute( xf, "Format", "HDF" )
                call xml_AddAttribute( xf, "NumberType", "Float" )
                call xml_AddAttribute( xf, "Precision", "4" )
                call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                call xml_AddCharacters( xf, trim( filename9 ) )
                call xml_EndElement( xf, "DataItem" )
                call xml_EndElement( xf, "Attribute" )
                ! Cell proportion
                do ks = 1, totgrn
                    filename11 = filename10
                    call append_nb( filename11, ks )
                    txt = 'Prop_'
                    call append_str( txt,material_name( ks ) )
                    call xml_NewElement( xf, "Attribute" )
                    call xml_AddAttribute( xf, "Type", "Scalar" )
                    call xml_AddAttribute( xf, "Center", "Cell" )
                    call xml_AddAttribute( xf, "Name", trim( txt ) )
                    call xml_NewElement( xf, "DataItem" )
                    call xml_AddAttribute( xf, "Format", "HDF" )
                    call xml_AddAttribute( xf, "NumberType", "Float" )
                    call xml_AddAttribute( xf, "Precision", "8" )
                    call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                    call xml_AddCharacters( xf, trim( filename11 ) )
                    call xml_EndElement( xf, "DataItem" )
                    call xml_EndElement( xf, "Attribute" )
                enddo
                call xml_EndElement( xf, "Grid" )
            enddo
            ! Footer
            call xml_EndElement( xf, "Grid" )
            call xml_EndElement( xf, "Domain" )
            call xml_EndElement( xf, "Xdmf" )
            call xml_Close( xf )
        endif

        return

    end subroutine Strata_xmf
    ! ============================================================================
    !> Subroutine Strata_series
    !! Subroutine Strata_series generates the XmL file for TIN surface.
    !<
    ! ============================================================================
    subroutine Strata_series( iter )

        ! Parameters Declaration
        type(xmlf_t) :: xf

        integer :: i, iter, it0
        character(len=128) :: filename, str, fname

        filename = 'Stratal_series.xdmf'
        call addpath1(filename)

        ! Header
        call xml_OpenFile( filename, xf )
        call xml_AddDOCTYPE(xf, "Xdmf", "Xdmf.dtd")
        call xml_DeclareNamespace(xf, "http://www.w3.org/2001/XInclude", "xi")
        call xml_NewElement(xf, "Xdmf" )
        call xml_AddAttribute(xf, "Version", "2.0" )
        call xml_NewElement( xf, "Domain" )
        call xml_NewElement( xf, "Grid" )
        call xml_AddAttribute( xf, "GridType", "Collection" )
        call xml_AddAttribute( xf, "CollectionType", "Temporal" )
        it0 = 1
        if( restart ) it0 = restart_iter + 1
        ! Loop over time step
        do i = it0, iter + 1
            ! Grid name
            fname = ''
            fname = 'Mesh'
            call noblnk( fname )
            str = '.'
            call append_str( fname,str )
            call append_zero( fname,i - 1)
            str = '.xmf'
            call append_str( fname,str )
            call xml_NewElement( xf, "xi:include" )
            call xml_AddAttribute( xf, "href", trim( fname ) )
            call xml_AddAttribute( xf, "xpointer", "xpointer(//Xdmf/Domain/Grid)" )
            call xml_EndElement( xf, "xi:include" )
        enddo
        ! Footer
        call xml_EndElement( xf, "Grid" )
        call xml_EndElement( xf, "Domain" )
        call xml_EndElement( xf, "Xdmf" )
        call xml_Close( xf )

        return

    end subroutine Strata_series
  ! ============================================================================

end module strata_out
