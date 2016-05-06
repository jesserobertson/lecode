! ============================================================================
! Name        : Strata_init.f90
! Author      : tristan salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file Strata_init.f90
!!
!! Strata_init is called at the beginning of the simulation to
!! initialize some parameters.
!!
!<
! ============================================================================
module strata_ini

    use file_data
    use mpidata
    use time_data
    use fwalker_data
    use error_data
    use strata_data
    use forces_data
    use mod_diffuse
    use sediment_data
    use kdtree2_module
    use kdtree2_precision_module

    implicit none

    public

    ! Stratal sediment thickness
    real( tkind ),dimension(:,:),allocatable :: init_hard
    ! Stratal sediment thickness
    real( tkind ),dimension(:,:,:),allocatable :: init_sedh

contains

    ! ============================================================================
    !> Subroutine generate_global_stratigraphy()
    !! Generate the initial stratigraphic grid.
    !<
    ! ============================================================================
    subroutine generate_global_stratigraphy

        ! Parameters Declaration
        logical :: found
        integer :: n, iunit, ios, p, k, nb, id, l

        real( tkind ) :: tolerance

        ! Check file exist
        iunit = 40 + iam
        inquire(FILE=fstrata, EXIST=found)
        if(found)then
            open(iunit,file=fstrata,status="old",action="read",iostat=ios)
            rewind(iunit)
        else
            attempt = STRAT_FILE
        endif

        ! Ensure that all processes have successfully read the node file.
        call completion

        ! Check the first line of the node file
        nb = strat_X * strat_Y
        G_strat%nodes_nb = nb

        ! Allocate stratal grid arrays
        allocate( gv_SharenID( G_strat%nodes_nb, nproc ), gv_orientation( nb )  )
        allocate( gv_halo( nb ), gv_ShareID( nb, 2 ), gv_ngbID( nb, 8 ), elev_record( nb ) )
        allocate( gv_vdisp( nb ), gv_coord( nb, 3 ) )

        ! Allocate the nodes values
        G_strat%faces_nb = ( strat_X - 1 ) * ( strat_Y -1 )
        do n = 1, G_strat%nodes_nb
            read(iunit, *) k, gv_coord( n, 1 ), gv_coord( n, 2 ), gv_coord( n, 3 )
            gv_halo( n ) = 0
            gv_ShareID( n, 1:2) = -1
        enddo

        ! Close file
        close( iunit )

        if( allocated( gf_pid ) ) deallocate( gf_pid )
        if( allocated( gf_points ) ) deallocate( gf_points )
        if( allocated( gf_XYcoord ) ) deallocate( gf_XYcoord )
        allocate( gf_pid( G_strat%faces_nb ) )
        allocate( gf_points( G_strat%faces_nb, 4 ) )
        allocate( gf_XYcoord( G_strat%faces_nb, 2 ) )

        id = 1
        k = 1
        tolerance = 1.e-4_8
        ! Allocate the faces values
        do n = 1, G_strat%faces_nb
            ! Points
            gf_points( n, 1 ) = id
            gf_points( n, 2 ) = id + 1
            gf_points( n, 3 ) = strat_X + id + 1
            gf_points( n, 4 ) = strat_X + id
            ! Coordinates
            gf_XYcoord( n, 1 ) = real( gv_coord( id, 1 ) + strat_dx * 0.50_8 )
            gf_XYcoord( n, 2 ) = real( gv_coord( id, 2 ) + strat_dx * 0.50_8 )
            id = id + 1
            if( id == k * strat_X )then
                id = k * strat_X + 1
                k = k + 1
            endif
        enddo

        ! Allocate the initial deposit if it exist and if this is not a
        ! restart simulation
        if( InitDep > 0 .and. .not. restart )then

            if( allocated( init_hard ) ) deallocate( init_hard )
            if( allocated( init_sedh ) ) deallocate( init_sedh )
            allocate( init_hard( InitDep, G_strat%nodes_nb ) )
            allocate( init_sedh( InitDep, G_strat%nodes_nb, max_grn ) )

            ! Read through the initial deposit layers
            do n = 1 , InitDep

                ! Open deposit file
                inquire(FILE=fdep( n ), EXIST=found)
                if(found)then
                    open(iunit,file=fdep( n ),status="old",action="read",iostat=ios)
                    rewind(iunit)
                else
                    attempt = DEPO_FILE
                endif

                ! Ensure that all processes have successfully open the file.
                call completion

                ! Read the sediment thicknesses
                do p = 1, G_strat%nodes_nb
                    read(iunit, *) l, init_sedh( n, p, 1:totgrn ), init_hard( n, p )
                enddo

                close( iunit )
            enddo

        ! Otherwise if there is not initial deposit define a ghost layer at the bottom
        elseif( InitDep == 0 .and. .not. restart )then

            InitDep = 1
            if( allocated( init_hard ) ) deallocate( init_hard )
            if( allocated( init_sedh ) ) deallocate( init_sedh )
            allocate( init_hard( InitDep, G_strat%nodes_nb ) )
            allocate( init_sedh( InitDep, G_strat%nodes_nb, max_grn ) )

            do p = 1, G_strat%nodes_nb
                do l = 1, totgrn
                    init_sedh( 1, p, l ) = 0.0_8
                enddo
                init_hard( 1, p ) = 1.0_8
            enddo
        ! If this is a restart simulation check numer of layer already existing
        else

            L_strat%layersID = InitDep + restart_iter + 1
            allocate( rec_NbLay( G_strat%nodes_nb ) )
            allocate( rec_topelev( G_strat%nodes_nb ) )
            allocate( rec_coord( G_strat%nodes_nb, 3 ) )

        endif

        allocate( pick_refine( G_strat%nodes_nb ) )
        pick_refine = 0

        return

    end subroutine generate_global_stratigraphy
    ! ============================================================================
    !> Subroutine generate_local_stratigraphy()
    !! Allocate the stratigraphy to each processors depending on the band partitioning.
    !<
    ! ============================================================================
    subroutine generate_local_stratigraphy

        ! Parameters Declaration
        integer :: nb_lay, n, id, p, GID, m, k, ll, l, r
        integer :: nID, p1, strat_width, pass, totnb

        integer, dimension( nproc ) :: valnb
        integer, dimension( : ), allocatable ::  lid
        integer, dimension( : ), allocatable ::  lidT
        integer,dimension( : ), allocatable :: h_gid

        real( tkind ) :: zh

        ! Allocate arrays
        allocate( uplift( G_strat%nodes_nb ) )
        allocate( lid( G_strat%nodes_nb ), lidT( G_strat%nodes_nb*nproc ) )

        ! Define the total number of sedimentary layers for the simulation
        if( .not. restart )then
            nb_lay = InitDep + 1
            L_strat%layersID = nb_lay
            G_strat%total_layer = nb_lay + int ( ( time_end - time_start ) / time_display )
        else
            G_strat%total_layer = L_strat%layersID + int ( ( time_end - time_start ) / time_display )
        endif

        ! Find the number of faces lying in each partition
        id = 0
        do n = 1, G_strat%faces_nb
            if( gf_pid( n ) == iam )then
                id = id + 1
                L_strat%faces_nb = id
            endif
        enddo

        ! Allocate faces parameters arrays
        if( allocated( locf_gid ) ) deallocate( locf_gid )
        if( allocated( locf_points ) ) deallocate( locf_points )
        if( allocated( locf_XYcoord ) ) deallocate( locf_XYcoord )
        allocate( locf_gid( L_strat%faces_nb ) )
        allocate( locf_points( L_strat%faces_nb, 4 ) )
        allocate( locf_XYcoord( L_strat%faces_nb, 2 ) )

        ! Compute faces center coordinates
        id = 0
        do n = 1, G_strat%faces_nb
            if( gf_pid( n ) == iam )then
                id = id + 1
                locf_gid( id ) = n
                L_strat%faces_nb = id
                locf_XYcoord( id, 1:2 ) = real( gf_XYcoord( n, 1:2 ) )
            endif
        enddo

        ! Define processor id for each node based on faces
        do n = 1, G_strat%faces_nb
            do k = 1, 4
                nID = gf_points( n, k )
                if( nID > 0 .and. nID <= G_strat%nodes_nb )then
                   if( gv_halo( niD ) > 0 )then
                      do m = 1, gv_halo( niD )
                         if( gv_ShareID( nID, m ) == gf_pid( n ) )goto 20
                      enddo
                   endif
                   gv_halo( niD ) =  gv_halo( niD ) + 1
                   if( gv_halo( niD ) > 2 )print*,'Something went wrong when allocating processor ID to nodes'
                   gv_ShareID( nID, gv_halo( niD ) ) = gf_pid( n )
                else
                   print*,'Something went wrong when looking for the points ID of stratal faces'
                   stop
                endif
20          continue
            enddo
        enddo

        ! Define local nodes number for each processors
        id = 0
        do n = 1, G_strat%nodes_nb
            do k = 1, gv_halo( n )
                if( gv_ShareID( n, k ) == iam )then
                    id = id + 1
                    L_strat%nodes_nb = id
                endif
            enddo
        enddo

        ! Now that we have the number of points on each processor
        ! define the local points arrays
        if( allocated( locv_NbLay ) ) deallocate( locv_NbLay )
        if( allocated( locv_gid ) ) deallocate( locv_gid )
        if( allocated( locv_vdisp ) ) deallocate( locv_vdisp )
        if( allocated( locv_vrate ) ) deallocate( locv_vrate )
        if( allocated( locv_coord ) ) deallocate( locv_coord )
        allocate( locv_NbLay( L_strat%nodes_nb ) )
        allocate( locv_gid( L_strat%nodes_nb ) )
        allocate( locv_vdisp( L_strat%nodes_nb ) )
        allocate( locv_vrate( L_strat%nodes_nb ) )
        allocate( locv_coord( L_strat%nodes_nb, 3 ) )

        ! Define the stratigraphic arrays
        if( allocated( strat_layID ) ) deallocate( strat_layID )
        if( allocated( strat_hardness ) ) deallocate( strat_hardness )
        if( allocated( strat_zelev ) ) deallocate( strat_zelev )
        if( allocated( strat_porosity ) ) deallocate( strat_porosity )
        if( allocated( strat_thick ) ) deallocate( strat_thick )
        if( allocated( strat_fac ) ) deallocate( strat_fac )
        if( allocated( strat_sedh ) ) deallocate( strat_sedh )
        allocate( strat_layID( L_strat%nodes_nb, G_strat%total_layer ) )
        allocate( strat_hardness( L_strat%nodes_nb, G_strat%total_layer ) )
        allocate( strat_thick( L_strat%nodes_nb, G_strat%total_layer ) )
        allocate( strat_zelev( L_strat%nodes_nb, G_strat%total_layer ) )
        allocate( strat_porosity( L_strat%nodes_nb, G_strat%total_layer, max_grn ) )
        allocate( strat_fac( L_strat%nodes_nb, G_strat%total_layer ) )
        allocate( strat_sedh( L_strat%nodes_nb, G_strat%total_layer, max_grn ) )
        allocate( check_int( L_strat%nodes_nb*2 ) )

        ! Define the offset
        disp( 1 ) = 0
        do k = 1, nproc - 1
            disp( k + 1 ) = disp( k ) + G_strat%nodes_nb
        enddo

        ! For global points find local id
        valnb = G_strat%nodes_nb
        id = 0
        do n = 1, G_strat%nodes_nb
            lid( n ) = 0
            do k = 1, gv_halo( n )
                if( gv_ShareID( n, k ) == iam )then
                    id = id + 1
                    lid( n ) = id
                    locv_gid( id ) = n
                    L_strat%nodes_nb = id
                    locv_NbLay( id ) = 1
                    locv_vdisp( id ) = 0.0_8
                    locv_vrate( id ) = 0.0_8
                    locv_coord( id, 1:3 ) = real( gv_coord( n, 1:3 ) )
                    strat_layID(  id , : ) = 0
                endif
            enddo
        enddo

        ! Gather information to all processors
        lidT = -1
        call mpi_allgatherv(lid,G_strat%nodes_nb,int_type,lidT,valnb,disp,int_type,SPModel_comm_world,ierr)

        ! Find shared processors containing the node
        do n = 1, G_strat%nodes_nb
            do k = 1, nproc
                id = n + ( k - 1 )*G_strat%nodes_nb
                gv_SharenID( n, k ) =  lidT( id )
            enddo
        enddo

        ! Compute faces points ID
        p = 1
        p1 = 1
        strat_width = strat_X
        if( strat_X > strat_Y ) strat_width = strat_col( iam + 1 )
        do k = 1, L_strat%faces_nb
            locf_points( k, 1 ) = p
            locf_points( k, 2 ) = p + 1
            locf_points( k, 3 ) = p + strat_width
            locf_points( k, 4 ) = p + strat_width + 1
            if( p1 < strat_width - 1 )then
                p1 = p1 + 1
                p = p + 1
            else
                p1 = 1
                p = p + 2
            endif
        enddo

        ! Create the basement
        do n = 1, L_strat%nodes_nb
            strat_hardness(  n ,  1  ) = 1.0_8
            strat_layID( n, 1 ) = 1
            strat_thick(  n ,  1  ) = 0.0_8
            strat_porosity(  n ,  1, 1:max_grn  ) = 0.0_8
            strat_zelev(  n ,  1  ) = real( locv_coord( n, 3 ) )
            strat_fac(  n ,  1  ) = strat_zelev(  n ,  1  ) - gsea%actual_sea
            strat_sedh( n, 1, 1:max_grn ) = 0.0_8
        enddo

        ! Create the initial deposits
        if( InitDep > 0 .and. .not. restart )then
            locv_NbLay = 1
            do p = 2, nb_lay

                do n = 1, L_strat%nodes_nb
                    GID = locv_gid( n )
                    zh = 0.0_8

                    do k = 1, totgrn
                        zh = zh + init_sedh( p-1, GID, k )
                    enddo

                    if( zh > 0.0_8 )then
                        ll = locv_NbLay( n ) + 1
                        strat_sedh(  n ,  ll ,  :  ) = 0.0_8

                        do k = 1, totgrn
                            strat_sedh(  n ,  ll ,  k  ) =  real( init_sedh( p-1, GID, k ) )
                        enddo

                        strat_hardness(  n ,  ll  ) =  init_hard( p-1, GID )
                        if( strat_hardness(  n ,  ll  ) == 0.0_8 ) strat_hardness(  n ,  ll  ) = 1.0_8
                        strat_thick(  n ,  ll  ) = real( zh )
                        strat_porosity(  n ,  ll, 1:totgrn  ) = 0.0_8
                        strat_layID( n, ll ) = p
                        locv_NbLay( n ) = ll
                        strat_zelev(  n ,  ll  ) = real( strat_zelev(  n ,  ll - 1  ) + zh )
                        strat_fac(  n ,  ll  ) = real( strat_zelev(  n ,  ll  ) - gsea%actual_sea )
                    endif
                enddo
            enddo
            ! Deallocate local arrays
            deallocate( init_sedh, init_hard )

        ! In case this is a restart from a previous simulation
        ! copy values from recorded arrays
        else

            do n = 1, L_strat%nodes_nb
                GID = locv_gid( n )
                locv_NbLay( n ) = rec_NbLay( GID )
                do p = 1, rec_NbLay( GID )
                    strat_layID( n, p ) =  rec_slayID( GID , p )
                    strat_hardness(  n ,  p  ) =  rec_shardness( GID , p )
                    strat_zelev(  n ,  p  ) =  real( rec_szelev( GID, p ) )
                    strat_porosity(  n ,  p ,  1:totgrn  ) =  real( rec_sporosity( GID , p,  1:totgrn ) )
                    strat_thick(  n ,  p  ) =  real( rec_sthick( GID , p ) )
                    strat_fac(  n ,  p  ) =  rec_sfac( GID , p )
                    strat_sedh(  n ,  p ,  1:totgrn  ) =  real( rec_ssedh( GID , p , 1:totgrn ) )
                enddo
            enddo

            ! Deallocate record arrays
            deallocate( rec_coord, rec_topelev, rec_NbLay, rec_ssedh )
            deallocate( rec_slayID, rec_shardness, rec_szelev )
            deallocate( rec_sporosity, rec_sthick, rec_sfac )

        endif

        ! Broadcast top elevation to global grid
        allocate( proc_elev( G_strat%nodes_nb ) )
        if( .not. restart )then
            proc_elev = -1.0e6_8
            do n = 1, L_strat%nodes_nb
                GID = locv_gid( n )
                ll = locv_NbLay( n )
                proc_elev( GID ) = strat_zelev(  n ,  ll  )
            enddo
            call mpi_allreduce( proc_elev, elev_record, G_strat%nodes_nb, dbl_type, &
                max_type, SPModel_comm_world, ierr )
        endif

        ! Allocate neighbourhood
        gv_ngbID = -1
        do n = 1, G_strat%nodes_nb
            ! Corners
            if( n == 1 )then
                gv_ngbID( n, 1 ) = n + strat_X
                gv_ngbID( n, 2 ) = n + strat_X + 1
                gv_ngbID( n, 3 ) = n + 1
            elseif( n == strat_X )then
                gv_ngbID( n, 1 ) = n + strat_X
                gv_ngbID( n, 7 ) = n - 1
                gv_ngbID( n, 8 ) = n + strat_X - 1
            elseif( n == G_strat%nodes_nb - strat_X + 1 )then
                gv_ngbID( n, 3 ) = n + 1
                gv_ngbID( n, 4 ) = n - strat_X + 1
                gv_ngbID( n, 5 ) = n - strat_X
            elseif( n == G_strat%nodes_nb )then
                gv_ngbID( n, 5 ) = n - strat_X
                gv_ngbID( n, 6 ) = n - strat_X - 1
                gv_ngbID( n, 7 ) = n - 1
            ! Borders
            elseif( gv_coord( n, 2 ) == strat_yo )then
                gv_ngbID( n, 1 ) = n + strat_X
                gv_ngbID( n, 2 ) = n + strat_X + 1
                gv_ngbID( n, 3 ) = n + 1
                gv_ngbID( n, 7 ) = n - 1
                gv_ngbID( n, 8 ) = n + strat_X - 1
            elseif( gv_coord( n, 1 ) == strat_xo )then
                gv_ngbID( n, 1 ) = n + strat_X
                gv_ngbID( n, 2 ) = n + strat_X + 1
                gv_ngbID( n, 3 ) = n + 1
                gv_ngbID( n, 4 ) = n - strat_X + 1
                gv_ngbID( n, 5 ) = n - strat_X
            elseif( gv_coord( n, 1 ) == strat_xo + ( strat_X - 1 ) * strat_dx )then
                gv_ngbID( n, 1 ) = n + strat_X
                gv_ngbID( n, 5 ) = n - strat_X
                gv_ngbID( n, 6 ) = n - strat_X - 1
                gv_ngbID( n, 7 ) = n - 1
                gv_ngbID( n, 8 ) = n + strat_X - 1
            elseif( gv_coord( n, 2 ) == strat_yo + ( strat_Y - 1 ) * strat_dx )then
                gv_ngbID( n, 3 ) = n + 1
                gv_ngbID( n, 4 ) = n - strat_X + 1
                gv_ngbID( n, 5 ) = n - strat_X
                gv_ngbID( n, 6 ) = n - strat_X - 1
                gv_ngbID( n, 7 ) = n - 1
            else
                gv_ngbID( n, 1 ) = n + strat_X
                gv_ngbID( n, 2 ) = n + strat_X + 1
                gv_ngbID( n, 3 ) = n + 1
                gv_ngbID( n, 4 ) = n - strat_X + 1
                gv_ngbID( n, 5 ) = n - strat_X
                gv_ngbID( n, 6 ) = n - strat_X - 1
                gv_ngbID( n, 7 ) = n - 1
                gv_ngbID( n, 8 ) = n + strat_X - 1
            endif
        enddo

        ! Get the border of the halo points
        m = 0
        do n = 1, G_strat%nodes_nb
            if( gv_halo( n ) > 1 ) m = m + 1
        enddo

        m = 0
        do n = 1, G_strat%nodes_nb
            if( gv_halo( n ) == 1 )then
                pass = 0
                do k = 1, 8
                    id = gv_ngbID( n, k )
                    if( id > 0 )then
                        if( gv_halo( id ) > 1 .and. pass == 0 )then
                            if( gv_ShareID( n, 1 ) == iam )then
                                pass = 1
                                m = m + 1
                            endif
                        endif
                    endif
                enddo
            endif
        enddo

        ! Define halo border arrays
        haloborder_nb = m
        if( allocated( halo_proc ) ) deallocate( halo_proc )
        if( allocated( halo_lid ) ) deallocate( halo_lid )
        if( allocated( halo_gid ) ) deallocate( halo_gid, halo1_elev, halo2_elev )
        allocate( halo_proc( haloborder_nb ) )
        if( iam == 0 .or. iam == nproc - 1 )then
            totnb = haloborder_nb*2
        else
            totnb = haloborder_nb
        endif
        allocate( h_gid( totnb ), halo_gid( totnb, nproc ) )
        allocate( halo_lid( totnb ) )
        allocate( halo_send( totnb ) )

        ! Define global and local id of points in the halo borders
        m = 0
        l = 0
        r = int( totnb * 0.5 )
        halo_send = -1
        halo_proc = -1
        halo_lid = -1
        h_gid = -1
        do n = 1, G_strat%nodes_nb
            if( gv_halo( n ) == 1 )then
                pass = 0
                do k = 1, 8
                    id = gv_ngbID( n, k )
                    if( id > 0 )then
                        if( gv_halo( id ) > 1 .and. pass == 0 )then
                            if( gv_ShareID( n, 1 ) == iam )then
                                m = m + 1
                                pass = 1
                                halo_proc( m ) = gv_ShareID( n, 1 )
                                if( gv_ShareID( id, 1 ) == iam )then
                                    if( gv_ShareID( id, 2 ) > iam )then
                                        r = r + 1
                                        halo_send( r ) = gv_ShareID( id, 2 )
                                        halo_lid( r ) = gv_SharenID( n, halo_proc( m ) + 1 )
                                        h_gid( r ) = n
                                    else
                                        l = l + 1
                                        halo_send( l ) = gv_ShareID( id, 2 )
                                        halo_lid( l ) = gv_SharenID( n, halo_proc( m ) + 1 )
                                        h_gid( l ) = n
                                    endif
                                elseif( gv_ShareID( id, 2 ) == iam )then
                                    if( gv_ShareID( id, 1 ) > iam )then
                                        r = r + 1
                                        halo_send( r ) = gv_ShareID( id, 1 )
                                        halo_lid( r ) = gv_SharenID( n, halo_proc( m ) + 1 )
                                        h_gid( r ) = n
                                    else
                                        l = l + 1
                                        halo_send( l ) = gv_ShareID( id, 1 )
                                        halo_lid( l ) = gv_SharenID( n, halo_proc( m ) + 1 )
                                        h_gid( l ) = n
                                    endif
                                else
                                    print*,'Something went wrong when looking for stratal halo values'
                                    halo_send( m ) = -10
                                endif

                            endif
                        endif
                    endif
                enddo
            endif
        enddo

        ! Gather halo global ids from each processors
        call mpi_allgather( h_gid, totnb, int_type, halo_gid, totnb, int_type, SPModel_comm_world, ierr )

        haloproc_nb = 0
        if( iam == 0 .or. iam == nproc - 1 ) haloproc_nb = haloborder_nb
        if( iam /= 0 .and. iam /= nproc - 1 ) haloproc_nb = int( haloborder_nb*0.5 )

        ! Ids of the points that will be send to other processors
        allocate( halo_rcv2( haloproc_nb ), halo_rcv1( haloproc_nb ) )
        halo_rcv2 = -1
        halo_rcv1 = -1
        halo_ngb = -1
        if( iam > 0 )then
            halo_ngb( 1 ) = iam - 1
            halo_rcv2( 1:haloproc_nb ) = halo_gid( haloproc_nb+1:totnb, iam )
        endif
        if( iam < nproc - 1 )then
            halo_ngb( 2 ) = iam + 1
            halo_rcv1( 1:haloproc_nb ) = halo_gid( 1:haloproc_nb, iam+2 )
        endif

        allocate( halo1_elev( haloproc_nb ) )
        allocate( halo2_elev( haloproc_nb ) )
        deallocate( lid, lidT, h_gid )
        haloborder_nb = totnb

        return

    end subroutine generate_local_stratigraphy
    ! ============================================================================
    !> Subroutine create_diffusion_grid
    !! Subroutine create_diffusion_grid create the grid for diffusion.
    !<
    ! ============================================================================
    subroutine create_diffusion_grid

        logical :: inX

        integer :: k, i, p, kn, x1, y1
        integer, dimension( strat_col( iam + 1 ), strat_Y ) :: pid
        real( tkind ), dimension( strat_col( iam + 1 ), strat_Y, 3 ) :: xyz

        ! Find the neighbors processor ID
        diff_cols = strat_col( iam + 1 ) + 2
        if( iam == 0 )then
            diff_neighbors( 1 ) = -1
        else
            diff_neighbors( 1 ) = iam - 1
        endif
        if( iam == nproc - 1 .or. nproc == 1 )then
            diff_neighbors( 2 ) = -1
        else
            diff_neighbors( 2 ) = iam + 1
        endif

        ! Find the number of nodes on the diffusion grid
        if( strat_Y >= strat_X )then
            DL_nodesNb = diff_cols * ( strat_X + 2 )
            inX = .false.
        else
            DL_nodesNb = diff_cols * ( strat_Y + 2 )
            inX = .true.
        endif

        ! Allocate diffusion arrays
        allocate( diff_ID( DL_nodesNb ) )
        allocate( diff_ref( DL_nodesNb ) )
        allocate( diff_coord( DL_nodesNb, 3 ) )
        allocate( diff_ngb( DL_nodesNb, 8 ) )
        diff_ngb = -1
        diff_ID = -1
        diff_ref = 0
        diff_coord = -1.0e5_8

        ! Create the diffusion grid coordinates
        p = 1

        ! In case we split grid along Y
        if( .not. inX )then
            k = strat_X + 4
            do i = 1, L_strat%nodes_nb
                diff_coord( k ,1:3 ) = locv_coord( i, 1:3 )
                diff_ID( k ) = i
                diff_ref( k ) = 0
                ! West
                if( p == 1 )then
                    diff_coord( k-1 ,1 ) = locv_coord( i, 1 ) - strat_dx
                    diff_coord( k-1 ,2 ) = locv_coord( i, 2 )
                    diff_ref( k-1 ) = -3
                    if( bounds( 3 ) == 2 .or. bounds( 3 ) == 0 ) diff_ID( k-1 ) = 0
                    if( bounds( 3 ) == 1 ) diff_coord( k-1, 3 ) = 1.0e5_8
                    if( bounds( 3 ) == 0 ) diff_ref( k-1 ) = -7
                endif
                diff_ngb( k, 1 ) = k - strat_X - 2 - 1
                diff_ngb( k, 2 ) = k - strat_X - 2
                diff_ngb( k, 3 ) = k - strat_X - 2 + 1
                diff_ngb( k, 4 ) = k + 1
                diff_ngb( k, 5 ) = k + strat_X + 2 + 1
                diff_ngb( k, 6 ) = k + strat_X + 2
                diff_ngb( k, 7 ) = k + strat_X + 2 - 1
                diff_ngb( k, 8 ) = k - 1
                p = p + 1
                if( p == strat_X + 1 )then
                    ! East
                    diff_coord( k+1 ,1 ) = locv_coord( i, 1 ) + strat_dx
                    diff_coord( k+1 ,2 ) = locv_coord( i, 2 )
                    diff_ref( k+1 ) = -4
                    if( bounds( 4 ) == 2 .or. bounds( 4 ) == 0 ) diff_ID( k+1 ) = 0
                    if( bounds( 4 ) == 1 ) diff_coord( k+1, 3 ) = 1.0e5_8
                    if( bounds( 4 ) == 0 ) diff_ref( k+1 ) = -8
                    p = 1
                    k = k + 3
                else
                    k = k + 1
                endif
            enddo

            ! South
            diff_ref( 1:strat_X+2 ) = -2
            if( bounds( 2 ) == 2 .or. bounds( 2 ) == 0 ) diff_ID( 1:strat_X+2 ) = 0
            if( bounds( 2 ) == 0 ) diff_ref( 1:strat_X ) = -5
            if( iam == 0 )then
                if( bounds( 2 ) == 1 ) diff_coord( 1:strat_X+2 ,3 ) = 1.0e5_8
            else
                ! West
                if( bounds( 3 ) == 1 ) diff_coord( 1,3 ) = 1.0e5_8
                if( bounds( 3 ) == 2 .or. bounds( 3 ) == 0 ) diff_ID( 1 ) = 0
                if( bounds( 3 ) == 2 ) diff_ref( 1 ) = -3
                if( bounds( 3 ) == 0 ) diff_ref( 1 ) = -7
                ! East
                if( bounds( 4 ) == 1 ) diff_coord( strat_X+2,3 ) = 1.0e5_8
                if( bounds( 4 ) == 2 .or. bounds( 1 ) == 0 ) diff_ID( strat_X+2 ) = 0
                if( bounds( 4 ) == 2 ) diff_ref( strat_X+2 ) = -4
                if( bounds( 4 ) == 0 ) diff_ref( strat_X+2 ) = -8
            endif
            i = 1
            diff_coord( i, 1 ) = locv_coord( i, 1 ) - strat_dx
            diff_coord( i, 2 ) = locv_coord( i, 2 ) - strat_dx
            do k = 2, strat_X + 1
                diff_coord( k, 1 ) = locv_coord( i, 1 )
                diff_coord( k, 2 ) = locv_coord( i, 2 ) - strat_dx
                i = i + 1
            enddo
            diff_coord( strat_X + 2 ,1 ) = locv_coord( strat_X, 1 ) + strat_dx
            diff_coord( strat_X + 2 ,2 ) = locv_coord( strat_X, 2 ) - strat_dx

            ! North
            i = L_strat%nodes_nb - strat_X + 1
            kn = DL_nodesNb - ( strat_X + 2 ) + 1
            diff_ref( kn:DL_nodesNb ) = -1
            if( bounds( 1 ) == 2 .or. bounds( 1 ) == 0 ) diff_ID( kn:DL_nodesNb ) = 0
            if( bounds( 1 ) == 0 ) diff_ref( kn:DL_nodesNb ) = -6
            if( iam == nproc - 1 )then
                if( bounds( 1 ) == 1 ) diff_coord( kn:DL_nodesNb,3 ) = 1.0e5_8
            else
                ! West
                if( bounds( 3 ) == 1 ) diff_coord( kn,3 ) = 1.0e5_8
                if( bounds( 3 ) == 2 .or. bounds( 3 ) == 0 ) diff_ID( kn ) = 0
                if( bounds( 3 ) == 2 ) diff_ref( kn ) = -3
                if( bounds( 3 ) == 0 ) diff_ref( kn ) = -7
                ! East
                if( bounds( 4 ) == 1 ) diff_coord( DL_nodesNb,3 ) = 1.0e5_8
                if( bounds( 4 ) == 2 .or. bounds( 1 ) == 0 ) diff_ID( DL_nodesNb ) = 0
                if( bounds( 4 ) == 2 ) diff_ref( DL_nodesNb ) = -4
                if( bounds( 4 ) == 0 ) diff_ref( DL_nodesNb ) = -8
            endif
            diff_coord( kn, 1 ) = locv_coord( i, 1 ) - strat_dx
            diff_coord( kn, 2 ) = locv_coord( i, 2 ) + strat_dx
            do k = kn+1, DL_nodesNb - 1
                diff_coord( k, 1 ) = locv_coord( i, 1 )
                diff_coord( k, 2 ) = locv_coord( i, 2 ) + strat_dx
                i = i + 1
            enddo
            diff_coord( DL_nodesNb, 1 ) = locv_coord( L_strat%nodes_nb, 1 ) + strat_dx
            diff_coord( DL_nodesNb, 2 ) = locv_coord( L_strat%nodes_nb, 2 ) + strat_dx

        ! Otherwise we split grid along X
        else
            x1 = 1
            y1 = 1
            pid = 0
            do i = 1, L_strat%nodes_nb
                xyz( x1, y1, 1:3 ) = locv_coord( i, 1:3 )
                pid( x1, y1 ) = i
                x1 = x1 + 1
                if( x1 > strat_col( iam + 1 ) )then
                    x1 = 1
                    y1 = y1 + 1
                endif
            enddo
            k = strat_Y + 3
            do x1 = 1, strat_col( iam + 1 )
                ! South
                diff_coord( k,1 ) = xyz( x1, 1, 1 )
                diff_coord( k,2 ) = xyz( x1, 1, 2 ) - strat_dx
                diff_ref( k ) = -2
                diff_ngb( k, 1 ) = k - strat_Y - 2 - 1
                diff_ngb( k, 2 ) = k - strat_Y - 2
                diff_ngb( k, 3 ) = k - strat_Y - 2 + 1
                diff_ngb( k, 4 ) = k + 1
                diff_ngb( k, 5 ) = k + strat_Y + 2 + 1
                diff_ngb( k, 6 ) = k + strat_Y + 2
                diff_ngb( k, 7 ) = k + strat_Y + 2 - 1
                diff_ngb( k, 8 ) = k - 1
                if( bounds( 2 ) == 2 .or. bounds( 2 ) == 0 ) diff_ID( k ) = 0
                if( bounds( 2 ) == 1 ) diff_coord( k, 3 ) = 1.0e5_8
                if( bounds( 2 ) == 0 ) diff_ref( k ) = -5
                k = k + 1
                do y1 = 1, strat_Y
                    diff_coord( k, 1:3 ) = xyz( x1, y1, 1:3 )
                    diff_ID( k ) = pid( x1, y1 )
                    diff_ngb( k, 1 ) = k - strat_Y - 2 - 1
                    diff_ngb( k, 2 ) = k - strat_Y - 2
                    diff_ngb( k, 3 ) = k - strat_Y - 2 + 1
                    diff_ngb( k, 4 ) = k + 1
                    diff_ngb( k, 5 ) = k + strat_Y + 2 + 1
                    diff_ngb( k, 6 ) = k + strat_Y + 2
                    diff_ngb( k, 7 ) = k + strat_Y + 2 - 1
                    diff_ngb( k, 8 ) = k - 1
                    diff_ref( k ) = 0
                    k = k + 1
                enddo
                ! North
                diff_coord( k,1 ) = xyz( x1, strat_Y, 1 )
                diff_coord( k,2 ) = xyz( x1, strat_Y, 2 ) + strat_dx
                diff_ngb( k, 1 ) = k - strat_Y - 2 - 1
                diff_ngb( k, 2 ) = k - strat_Y - 2
                diff_ngb( k, 3 ) = k - strat_Y - 2 + 1
                diff_ngb( k, 4 ) = k + 1
                diff_ngb( k, 5 ) = k + strat_Y + 2 + 1
                diff_ngb( k, 6 ) = k + strat_Y + 2
                diff_ngb( k, 7 ) = k + strat_Y + 2 - 1
                diff_ngb( k, 8 ) = k - 1
                diff_ref( k ) = -1
                if( bounds( 1 ) == 2 .or. bounds( 1 ) == 0 ) diff_ID( k ) = 0
                if( bounds( 1 ) == 1 ) diff_coord( k, 3 ) = 1.0e5_8
                if( bounds( 1 ) == 0 ) diff_ref( k ) = -6
                k = k + 1
            enddo
            ! West
            diff_ref( 1:strat_Y+2 ) = -3
            if( bounds( 3 ) == 2 .or. bounds( 3 ) == 0 ) diff_ID( 1:strat_Y+2 ) = 0
            if( bounds( 3 ) == 0 ) diff_ref( 1:strat_Y+2 ) = -7
            if( iam == 0 )then
                if( bounds( 3 ) == 1 ) diff_coord( 1:strat_Y+2 ,3 ) = 1.0e5_8
            else
                ! South
                if( bounds( 2 ) == 1 ) diff_coord( 1,3 ) = 1.0e5_8
                if( bounds( 2 ) == 2 .or. bounds( 2 ) == 0 ) diff_ID( 1 ) = 0
                if( bounds( 2 ) == 2 ) diff_ref( 1 ) = -2
                if( bounds( 2 ) == 0 ) diff_ref( 1 ) = -5
                ! North
                if( bounds( 1 ) == 1 ) diff_coord( strat_Y+2,3 ) = 1.0e5_8
                if( bounds( 1 ) == 2 .or. bounds( 1 ) == 0 ) diff_ID( strat_Y+2 ) = 0
                if( bounds( 1 ) == 2 ) diff_ref( strat_Y+2 ) = -1
                if( bounds( 1 ) == 0 ) diff_ref( strat_Y+2 ) = -6
            endif
            diff_coord( 1, 1 ) = xyz( 1, 1, 1 ) - strat_dx
            diff_coord( 1, 2 ) = xyz( 1, 1, 2 ) - strat_dx
            do k = 2, strat_Y + 1
                diff_coord( k, 1 ) = diff_coord( 1, 1 )
                diff_coord( k, 2 ) = diff_coord( k-1, 2 ) + strat_dx
            enddo
            diff_coord( strat_Y+2,1 ) = diff_coord( 1, 1 )
            diff_coord( strat_Y+2,2 ) = xyz( 1, strat_Y, 2 ) + strat_dx
            ! East
            i = L_strat%nodes_nb - ( strat_Y + 2 ) + 1
            kn = DL_nodesNb -  ( strat_Y + 2 ) + 1
            diff_ref( kn:DL_nodesNb ) = -4
            if( bounds( 4 ) == 2 .or. bounds( 4 ) == 0 ) diff_ID( kn:DL_nodesNb ) = 0
            if( bounds( 4 ) == 0 ) diff_ref( kn:DL_nodesNb ) = -8
            if( iam == nproc - 1 )then
                if( bounds( 4 ) == 1 ) diff_coord( kn:DL_nodesNb,3 ) = 1.0e5_8
            else
                ! South
                if( bounds( 2 ) == 1 ) diff_coord( kn,3 ) = 1.0e5_8
                if( bounds( 2 ) == 2 .or. bounds( 2 ) == 0 ) diff_ID( kn ) = 0
                if( bounds( 2 ) == 2 ) diff_ref( kn ) = -2
                if( bounds( 2 ) == 0 ) diff_ref( kn ) = -5
                ! North
                if( bounds( 1 ) == 1 ) diff_coord( DL_nodesNb,3 ) = 1.0e5_8
                if( bounds( 1 ) == 2 .or. bounds( 1 ) == 0 ) diff_ID( DL_nodesNb ) = 0
                if( bounds( 1 ) == 2 ) diff_ref( DL_nodesNb ) = -1
                if( bounds( 1 ) == 0 ) diff_ref( DL_nodesNb ) = -6
            endif
            diff_coord( kn, 1 ) = xyz( strat_col( iam + 1 ), 1, 1 ) + strat_dx
            diff_coord( kn, 2 ) = xyz( strat_col( iam + 1 ), 1, 2 ) - strat_dx
            do k = kn+1, DL_nodesNb - 1
                diff_coord( k, 1 ) = diff_coord( kn, 1 )
                diff_coord( k, 2 ) = diff_coord( k-1, 2 ) + strat_dx
            enddo
            diff_coord( DL_nodesNb, 1 ) = diff_coord( kn, 1 )
            diff_coord( DL_nodesNb, 2 ) = xyz( strat_col( iam + 1 ), strat_Y, 2 ) + strat_dx
        endif

        ! Define additional diffusion arrays
        allocate( top_sedh( L_strat%nodes_nb, totgrn ) )
        allocate( top_sedprev( L_strat%nodes_nb, totgrn ) )
        top_sedh = 0.0_8
        top_sedprev = 0.0_8

        allocate( depo( L_strat%nodes_nb, totgrn ) )
        allocate( dstart( L_strat%nodes_nb, totgrn ) )
        allocate( topmax( L_strat%nodes_nb ) )
        allocate( cdif( L_strat%nodes_nb ) )
        allocate( difo( DL_nodesNb ) )
        allocate( difp( DL_nodesNb ) )

        return

    end subroutine create_diffusion_grid
    ! ============================================================================
    !> Subroutine closest_stratal_node
    !! Subroutine closest_stratal_node is used to determine the closest stratal node.
    !<
    ! ============================================================================
    subroutine closest_stratal_node( xpos, ypos, fcS, snode )

        integer :: Xid, Yid, snode, fcS, id, n
        integer, dimension( 4 ) :: N_id

        real( tkind ) :: xpos, ypos, mindist, dist

        Xid = int( ( xpos - strat_xo ) / strat_dx ) + 1
        Yid = int( ( ypos - strat_yo ) / strat_dx ) + 1
        if( Yid == 1 .and. Xid < strat_X )then
            fcS = Xid
        elseif( Yid == 1 .and. Xid == strat_X )then
            fcS = Xid - 1
        elseif( Xid == 1 .and. Yid < strat_Y )then
            fcS = ( Yid - 1 )*( strat_X - 1 ) + 1
        elseif( Yid == strat_Y .and. Xid < strat_X )then
            fcS = ( Yid - 2 )*( strat_X - 1 ) + Xid
        elseif( Yid == strat_Y .and. Xid == strat_X )then
            fcS = ( Yid - 1 )*( Xid - 1 )
        elseif( Yid < strat_Y .and. Xid == strat_X )then
            fcS = Yid * ( Xid - 1 )
        else
            fcS = ( strat_X - 1 ) * ( Yid - 1 ) + Xid
        endif

        ! Check picked face
        id = gf_points( fcS, 1 )
        if( xpos - gv_coord( id, 1 ) > strat_dx .or. &
            ypos - gv_coord( id, 2 ) > strat_dx ) then
            print*,'Problem the flow walker is not in the proper stratal face'
            print*,'id',Xid,strat_X,Yid,strat_Y,strat_dx
            print*,'strat',gv_coord( id, 1:2 )
            print*,'fw',xpos,ypos, fcS
            stop
        endif

        ! Find the closest point
        snode = -1
        mindist = 2.0_8 * strat_dx
        N_id = gf_points( fcS, 1:4 )
        do n = 1, 4
            dist = sqrt( ( xpos - gv_coord( N_id( n ), 1 ) )**2 + &
                ( ypos - gv_coord( N_id( n ), 2 ) )**2 )
            if( dist < mindist .and. N_id( n ) > 0 )then
                mindist = dist
                snode = N_id( n )
            endif
        enddo
        if( snode <= 0 )then
           print*,'Problem locating point ID',xpos,ypos
           print*,N_id
           stop
        endif

        return

    end subroutine closest_stratal_node
    ! ============================================================================
    !> Subroutine find_points_rectangular_shape
    !! Subroutine find_points_rectangular_shape given a rectangle shape check whether or
    !! not a point is inside the polygon.
    !! \param x0, y0, x, y, n, l, m
    !<
    ! ============================================================================
    subroutine find_points_rectangular_shape( xa, ya, xb, yb, n, l, m )

        integer, intent( in ) :: n
        real( tkind ), intent( in ) :: xa, ya, xb( n ), yb( n )

        integer, intent( out ) :: l, m

        integer :: i, n0

        real( tkind )  :: x0, y0, x(n), y(n)
        real( tkind ) :: angle, eps, pi, pi2, sum, theta, theta1, thetai, tol, u, v

        eps = epsilon(1.0_8)
        x0 = xa
        y0 = ya
        x(:) = xb(:)
        y(:) = yb(:)
        n0 = n
        if (x(1) == x(n) .and. y(1) == y(n)) n0 = n - 1
        pi = atan2( 0.0_8, -1.0_8 )
        pi2 = 2.0_8*pi
        tol = 4.0_8*eps*pi
        l = -1
        m = 0
        u = x(1) - x0
        v = y(1) - y0
        if (u == 0.0_8 .and. v == 0.0_8) go to 20
        if (n0 < 2) return
        theta1 = atan2(v, u)
        sum = 0.0
        theta = theta1
        do i = 2, n0
            u = x(i) - x0
            v = y(i) - y0
            if (u == 0.0_8 .and. v == 0.0_8) go to 20
            thetai = atan2(v, u)
            angle = abs(thetai - theta)
            if (abs(angle - pi) < tol) go to 20
            if (angle > pi) angle = angle - pi2
            if (theta > thetai) angle = -angle
            sum = sum + angle
            theta = thetai
        enddo
        angle = abs(theta1 - theta)
        if (abs(angle - pi) < tol) go to 20
        if (angle > pi) angle = angle - pi2
        if (theta > theta1) angle = -angle
        sum = sum + angle
        m = int( abs(sum)/pi2 + 0.2_8 )
        if(m == 0) return
        l = 1
        if (sum < 0.0_8) m = -m

        return

20      l = 0

        return

    end subroutine find_points_rectangular_shape
    ! ============================================================================
    !> Function compute_distance()
    !! Function distance is used to return the distance between two nodes of a considered face
    !! \param node1, node2
    !<
    ! ============================================================================
    function compute_distance( node1, node2  ) result( dist )

        real( tkind ), intent( in ) :: node1( 3 ), node2( 3 )
        real( tkind ) :: dist

        dist = ( node1( 1 ) - node2( 1 ) )**2
        dist = dist + ( node1( 2 )  - node2( 2 ) )**2
        dist = sqrt( dble( dist ) )

    end function compute_distance
    ! ============================================================================

end module strata_ini

