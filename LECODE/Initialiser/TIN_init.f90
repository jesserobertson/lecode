! ============================================================================
! Name        : TIN_init.f90
! Author      : tristan salles
! Created on: Aug 16, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file TIN_init.f90
! ============================================================================
!
! Description : TIN_function is used to defined several functions needed for the TIN grid.
!
!<
! ============================================================================
module TIN_function

    use file_data
    use mpidata
    use TIN_data
    use error_data
    use strata_data
    use param_data
    use kdtree2_module
    use kdtree2_precision_module

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine create_adaptive_TIN_grid
    !! Subroutine create_adaptive_TIN_grid create a TIN using the refinement values
    !! defined by the user.
    !<
    ! ============================================================================
    subroutine create_adaptive_TIN_grid

        ! Parameters Declaration
        integer :: iunit, ios, n, nbpt, id, val, i, p, iseed, gfc, N_id( 4 )

        real( tkind ) :: estp, wstp, nstp, sstp, ran, ran1
        real( tkind ) :: xlow, ylow, xhigh, yhigh, x, y, x1, x2, y1, y2, dh
        real( tkind ), dimension( 2 ) :: f1, f2, f3, f4, res

        ! Get the vertex from the stratal mesh
        nbpt = 0
        nstp = strat_xo
        sstp = strat_xo
        wstp = strat_yo + strat_dx
        estp = strat_yo + strat_dx
        xlow = strat_xo + grd_LR
        ylow = strat_yo + grd_LR
        xhigh = strat_xo + grd_HR
        yhigh = strat_yo + grd_HR
        pick_refine = 0

        do n = 1, G_strat%nodes_nb

            ! Border & low resolution nodes
            ! South side
            if( gv_coord( n, 2 ) == strat_yo )then
                if( gv_coord( n, 1 ) == sstp )then
                    pick_refine( n ) = 1
                    nbpt = nbpt + 1
                    sstp = sstp + strat_dx
                endif
               ! North side
            elseif( gv_coord( n, 2 ) == strat_yo + strat_dx * ( strat_Y- 1 ) )then
                if( gv_coord( n, 1 ) == nstp )then
                    pick_refine( n ) = 1
                    nbpt = nbpt + 1
                    nstp = nstp + strat_dx
                endif
               ! West side
            elseif( gv_coord( n, 1 ) == strat_xo )then
                if( gv_coord( n, 2 ) == wstp )then
                    pick_refine( n ) = 1
                    nbpt = nbpt + 1
                    wstp = wstp + strat_dx
                endif
               ! East side
            elseif( gv_coord( n, 1 ) == strat_xo + strat_dx * ( strat_X- 1 ) )then
                if( gv_coord( n, 2 ) == estp )then
                    pick_refine( n ) = 1
                    nbpt = nbpt + 1
                    estp = estp + strat_dx
                endif
               ! Low resolution
            elseif( gv_coord( n, 1 ) == xlow .and. gv_coord( n, 2 ) == ylow )then
                pick_refine( n ) = 2
                nbpt = nbpt + 1
                if( xlow + grd_LR < strat_xo + strat_dx * ( strat_X - 1 ) )then
                    xlow = xlow + grd_LR
                    if( strat_xo + strat_dx * ( strat_X - 1 ) - xlow < strat_dx )then
                        xlow = strat_xo + grd_LR
                        ylow = ylow + grd_LR
                        if( strat_yo + strat_dx * ( strat_Y - 1 ) - ylow < strat_dx )&
                            ylow = strat_yo + strat_dx * strat_Y
                    endif
                else
                    xlow = strat_xo + grd_LR
                    ylow = ylow + grd_LR
                    if( strat_yo + strat_dx * ( strat_Y - 1 ) - ylow < strat_dx )&
                        ylow = strat_yo + strat_dx * strat_Y
                endif
            endif

            ! High resolution & gradient criteria
            if( grd_HR < grd_LR )then
                if( gv_coord( n, 1 ) == xhigh .and. gv_coord( n, 2 ) == yhigh )then
                    if( facc( n ) >= accumulation_refine .and. pick_refine( n ) == 0  )then
                        pick_refine( n ) = 3
                        nbpt = nbpt + 1
                    endif
                    if( xhigh + grd_HR < strat_xo + strat_dx * ( strat_X - 1 ) )then
                        xhigh = xhigh + grd_HR
                        if( strat_xo + strat_dx * ( strat_X - 1 ) - xhigh < strat_dx )then
                            xhigh = strat_xo + grd_HR
                            yhigh = yhigh + grd_HR
                            if( strat_yo + strat_dx * ( strat_Y - 1 ) - yhigh < strat_dx )&
                                yhigh = strat_yo + strat_dx * strat_Y
                        endif
                    else
                        xhigh = strat_xo + grd_HR
                        yhigh = yhigh + grd_HR
                        if( strat_yo + strat_dx * ( strat_Y - 1 ) - yhigh < strat_dx )&
                            yhigh = strat_yo + strat_dx * strat_Y
                    endif
                endif
            endif
        enddo

        ! Create the TIN refined file
        if( iam == 0 )then
            fnode = 'TIN.node'
            call addpath2( fnode )
            iunit = 20
            open(iunit,file=fnode,status="replace",action="write",iostat=ios)
            rewind(iunit)
        endif
        just_nodes = nbpt
        G_tin%nodes_nb = nbpt + ( strat_X + strat_Y + 2 ) * 2
        if( iam == 0 )write(iunit, '(I10,1X,3(I2,1X))') G_tin%nodes_nb, 2, 1, 0

        if( allocated( flow_accu ) ) deallocate( flow_accu )
        allocate( flow_accu( just_nodes ) )
        id = 0

        if( .not. allocated( Tin_Gproc ) ) allocate( Tin_Gproc( nproc ) )
        Tin_Gproc = 0
        Tif_Gproc = 0
        tinv_ngb = 0
        tinv_pid = -1
        do n = 1, G_strat%nodes_nb
            if( pick_refine( n ) >= 1 )then
                id = id + 1
                if( pick_refine( n ) == 1 )then
                    x = real( gv_coord( n, 1 ) )
                    y = real( gv_coord( n, 2 ) )
                    tinv_coord( id, 1:2 ) = gv_coord( n, 1:2 )
                    tinv_coord( id, 3 ) = wdem( n )
                    flow_accu( id ) = facc( n )
                    tinv_boundary( id ) = 0
                    tinv_pid( id ) = gv_ShareID( n, 1 )
                    Tin_Gproc( tinv_pid( id ) + 1 ) = Tin_Gproc( tinv_pid( id ) + 1 ) + 1
                elseif( pick_refine( n ) == 2 )then
                    call random_seed( size= iseed )
                    call random_number( harvest= ran )
                    x = real( gv_coord( n, 1 ) + grd_LR * 0.1_8 * ( ran * 2.0_8 - 1.0_8 ) )
                    call random_seed( size= iseed )
                    call random_number( harvest= ran1 )
                    y = real( gv_coord( n, 2 ) + grd_LR * 0.1_8 * ( ran1 * 2.0_8 - 1.0_8 ) )
                    tinv_coord( id, 1 ) = x
                    tinv_coord( id, 2 ) = y
                    ! Get global face id containing the point
                    call check_point_TIN_face_containing_point( x, y, gfc )
                    tinv_pid( id ) = gf_pid( gfc )
                    Tin_Gproc( tinv_pid( id ) + 1 ) = Tin_Gproc( tinv_pid( id ) + 1 ) + 1
                    N_id = gf_points( gfc, 1:4 )
                    x1 = gv_coord( N_id( 1 ), 1 )
                    x2 = gv_coord( N_id( 3 ), 1 )
                    y1 = gv_coord( N_id( 1 ), 2 )
                    y2 = gv_coord( N_id( 3 ), 2 )
                    f1( 1 ) = wdem( N_id( 1 ) )
                    f2( 1 ) = wdem( N_id( 2 ) )
                    f3( 1 ) = wdem( N_id( 4 ) )
                    f4( 1 ) = wdem( N_id( 3 ) )
                    f1( 2 ) = facc( N_id( 1 ) )
                    f2( 2 ) = facc( N_id( 2 ) )
                    f3( 2 ) = facc( N_id( 4 ) )
                    f4( 2 ) = facc( N_id( 3 ) )
                    call bilinear_interp( x, y, x1, x2, y1, y2, f1, f2, f3, f4, res )
                    tinv_coord( id, 3 ) = res( 1 )
                    flow_accu( id ) = res( 2 )
                    tinv_boundary( id ) = 0
                else
                    call random_seed( size= iseed )
                    call random_number( harvest= ran )
                    x = real( gv_coord( n, 1 ) + grd_HR * 0.2_8 * ( ran * 2.0_8 - 1.0_8 ) )
                    call random_seed( size= iseed )
                    call random_number( harvest= ran1 )
                    y = real( gv_coord( n, 2 ) + grd_HR * 0.2_8 * ( ran1 * 2.0_8 - 1.0_8 ) )
                    tinv_coord( id, 1 ) = x
                    tinv_coord( id, 2 ) = y
                    ! Get global face id containing the point
                    call check_point_TIN_face_containing_point( x, y, gfc )
                    tinv_pid( id ) = gf_pid( gfc )
                    Tin_Gproc( tinv_pid( id ) + 1 ) = Tin_Gproc( tinv_pid( id ) + 1 ) + 1
                    N_id = gf_points( gfc, 1:4 )
                    x1 = gv_coord( N_id( 1 ), 1 )
                    x2 = gv_coord( N_id( 3 ), 1 )
                    y1 = gv_coord( N_id( 1 ), 2 )
                    y2 = gv_coord( N_id( 3 ), 2 )
                    f1( 1 ) = wdem( N_id( 1 ) )
                    f2( 1 ) = wdem( N_id( 2 ) )
                    f3( 1 ) = wdem( N_id( 4 ) )
                    f4( 1 ) = wdem( N_id( 3 ) )
                    f1( 2 ) = facc( N_id( 1 ) )
                    f2( 2 ) = facc( N_id( 2 ) )
                    f3( 2 ) = facc( N_id( 4 ) )
                    f4( 2 ) = facc( N_id( 3 ) )
                    call bilinear_interp( x, y, x1, x2, y1, y2, f1, f2, f3, f4, res )
                    tinv_coord( id, 3 ) = res( 1 )
                    flow_accu( id ) = res( 2 )
                    tinv_boundary( id ) = 0

                endif
                if( iam == 0 ) write(iunit,'(I10,1X,2(F16.3,1X),F4.1,1X,I2)')id,x,y,0.0,0
            endif
        enddo

        if( id /= just_nodes )&
            print*,'Something went wrong: values mismatch in building TIN grid'

        ! Then create the ghost coordinates
        val =  just_nodes + 1
        ! 1- the south ghosts
        y = strat_yo - strat_dx
        p = 1
        do i = 1, strat_X + 2
            x = strat_xo + ( i - 2 )*strat_dx
            if( iam == 0 ) write(iunit, '(I10,1X,2(F16.3,1X),F4.1,1X,I2)') val, x, y, 0.0, 1
            tinv_coord( val, 1 ) = x
            tinv_coord( val, 2 ) = y
            tinv_boundary( val ) = 1
            if( bounds( 2 ) == 2 ) tinv_boundary( val ) = 2
            if( bounds( 2 ) == 1 ) tinv_boundary( val ) = 3
            if( bounds( 2 ) == 3 ) tinv_boundary( val ) = 4
            if( tinv_boundary( val ) == 1 ) tinv_coord( val,  3 ) = wdem( p )
            if( tinv_boundary( val ) == 2 )then
                dh = wdem( p + strat_X ) - wdem( p )
                tinv_coord( val,  3 ) = wdem( p ) - dh
            endif
            if( tinv_boundary( val ) == 3 ) tinv_coord( val,  3 ) = 1.0e5_8
            if( tinv_boundary( val ) == 4 ) tinv_coord( val,  3 ) = -1.0e5_8
            tinv_pid( val ) = gv_ShareID( p, 1 )
            Tin_Gproc( tinv_pid( val ) + 1 ) = Tin_Gproc( tinv_pid( val ) + 1 ) + 1
            if( i > 2 .and. i < strat_X + 1 ) p = p + 1
            val = val + 1
        enddo

        ! 2- the north ghosts
        y = strat_yo + strat_Y * strat_dx
        p = G_strat%nodes_nb - strat_X + 1
        do i = 1, strat_X + 2
            x = strat_xo + ( i - 2 )*strat_dx
            if( iam == 0 ) write(iunit, '(I10,1X,2(F16.3,1X),F4.1,1X,I2)') val, x, y, 0.0, 1
            tinv_coord( val, 1 ) = x
            tinv_coord( val, 2 ) = y
            tinv_boundary( val ) = 1
            if( bounds( 1 ) == 2 ) tinv_boundary( val ) = 2
            if( bounds( 1 ) == 1 ) tinv_boundary( val ) = 3
            if( bounds( 1 ) == 3 ) tinv_boundary( val ) = 4
            if( tinv_boundary( val ) == 1 ) tinv_coord( val,  3 ) = wdem( p )
            if( tinv_boundary( val ) == 2 )then
                dh = wdem( p - strat_X ) - wdem( p )
                tinv_coord( val,  3 ) = wdem( p ) - dh
            endif
            if( tinv_boundary( val ) == 3 ) tinv_coord( val,  3 ) = 1.0e5_8
            if( tinv_boundary( val ) == 4 ) tinv_coord( val,  3 ) = -1.0e5_8
            tinv_pid( val ) = gv_ShareID( p, 1 )
            Tin_Gproc( tinv_pid( val ) + 1 ) = Tin_Gproc( tinv_pid( val ) + 1 ) + 1
            if( i > 2 .and. i < strat_X + 1 ) p = p + 1
            val = val + 1
        enddo

        ! 3- the west ghosts
        x = strat_xo - strat_dx
        p = 1
        do i = 1, strat_Y
            y = strat_yo + ( i - 1 )*strat_dx
            if( iam == 0 ) write(iunit, '(I10,1X,2(F16.3,1X),F4.1,1X,I2)') val, x, y, 0.0, 1
            tinv_coord( val, 1 ) = x
            tinv_coord( val, 2 ) = y
            tinv_boundary( val ) = 1
            if( bounds( 3 ) == 2 ) tinv_boundary( val ) = 2
            if( bounds( 3 ) == 1 ) tinv_boundary( val ) = 3
            if( bounds( 3 ) == 3 ) tinv_boundary( val ) = 4
            if( tinv_boundary( val ) == 2 )then
                dh = wdem( p + 1 ) - wdem( p )
                tinv_coord( val,  3 ) = wdem( p ) - dh
            endif
            if( tinv_boundary( val ) == 1 ) tinv_coord( val,  3 ) = wdem( p )
            if( tinv_boundary( val ) == 3 ) tinv_coord( val,  3 ) = 1.0e5_8
            if( tinv_boundary( val ) == 4 ) tinv_coord( val,  3 ) = -1.0e5_8
            tinv_pid( val ) = gv_ShareID( p, 1 )
            Tin_Gproc( tinv_pid( val ) + 1 ) = Tin_Gproc( tinv_pid( val ) + 1 ) + 1
            p = p + strat_X
            val = val + 1
        enddo

        ! 4- the east ghosts
        x = strat_xo + strat_X * strat_dx
        p = strat_X
        do i = 1, strat_Y
            y = strat_yo + ( i - 1 )*strat_dx
            if( iam == 0 ) write(iunit, '(I10,1X,2(F16.3,1X),F4.1,1X,I2)') val, x, y, 0.0, 1
            tinv_coord( val, 1 ) = x
            tinv_coord( val, 2 ) = y
            tinv_boundary( val ) = 1
            if( bounds( 4 ) == 2 ) tinv_boundary( val ) = 2
            if( bounds( 4 ) == 1 ) tinv_boundary( val ) = 3
            if( bounds( 4 ) == 3 ) tinv_boundary( val ) = 4
            if( tinv_boundary( val ) == 2 )then
                dh = wdem( p - 1 ) - wdem( p )
                tinv_coord( val,  3 ) = wdem( p ) - dh
            endif
            if( tinv_boundary( val ) == 1 ) tinv_coord( val,  3 ) = wdem( p )
            if( tinv_boundary( val ) == 3 ) tinv_coord( val,  3 ) = 1.0e5_8
            if( tinv_boundary( val ) == 4 ) tinv_coord( val,  3 ) = -1.0e5_8
            tinv_pid( val ) = gv_ShareID( p, 1 )
            Tin_Gproc( tinv_pid( val ) + 1 ) = Tin_Gproc( tinv_pid( val ) + 1 ) + 1
            p = p + strat_X
            val = val + 1
        enddo

        if( iam == 0 ) write(iunit,*)' '
        if( iam == 0 ) close( iunit )

        ftin = 'GTIN'

        if( tinv_coord( 1, 2 ) /= tinv_coord( 2, 2 ) )&
            attempt = TIN_NOD
        if( tinv_coord( 2, 1 ) - tinv_coord( 1, 1 ) /= strat_dx )&
            attempt = TIN_SPAC

        ! Ensure that all processes have successfully completed subroutine.
        call completion

        return

    end subroutine create_adaptive_TIN_grid
    ! ============================================================================
    !> Subroutine search_tin_face_neighbors_kdtree
    !! Subroutine search_tin_face_neighbors_kdtree determines the faces neighbors
    !! for the TIN grid.
    !<
    ! ============================================================================
    subroutine search_tin_face_neighbors_kdtree( maxneigh )

        ! Parameters Declaration
        logical :: inside

        integer :: k, kn, fid, n1, n2, nb, maxneigh
        integer, dimension( 3 ) :: nids, nids2
        integer, dimension( mxnghb ) :: nid

        real( tkind ),dimension( :, : ), allocatable :: Fdata

        type(kdtree2), pointer :: Ftree
        type(kdtree2_result), dimension( : ), allocatable :: FRslt

        allocate( Fdata( 2,  G_tin%faces_nb ) )
        allocate( FRslt( maxneigh ) )

        ! Create the kd-tree
        do k = 1, G_tin%faces_nb
            Fdata( 1, k ) =  tinf_centroid( k, 1 )
            Fdata( 2, k ) =  tinf_centroid( k, 2 )
        enddo
        Ftree => kdtree2_create( Fdata, sort = .true., rearrange = .true. )

        ! For each TIN faces find the neighboring faces IDs
        nid = -1
        do k = 1,  G_tin%faces_nb
            nids = tinf_nids( k, 1:3 )
            nb = 0
            call kdtree2_n_nearest_around_point(Ftree, idxin=k, nn=maxneigh, correltime=1, results=FRslt)
            do kn = 1, maxneigh
                fid = FRslt( kn )%idx
                nids2 = tinf_nids( fid, 1:3 )
                inside = .false.

                same_points_face: do n1 = 1, 3
                    do n2 = 1, 3
                        if( nids2( n1 ) == nids( n2 ) )then
                            inside = .true.
                            exit same_points_face
                        endif
                    enddo
                enddo same_points_face

                if( inside )then
                    nb = nb + 1
                    if( nb > mxnghb )&
                        print*,'Something went worng the number of TIN face neighbors is greater than the maximum allowed.'
                    nid( nb ) = fid
                    if( nid( nb ) < 1 .or. nid( nb ) > G_tin%faces_nb )then
                        print*,'Something went worng when looking for face neighbors'
                        stop
                    endif
                endif
            enddo
            tinf_nghb( k ) = nb
            fids( k, 1:nb ) = nid( 1:nb )
        enddo

        ! Destroy kd-tree
        call kdtree2_destroy(Ftree)
        deallocate( Fdata, FRslt )

        return

    end subroutine search_tin_face_neighbors_kdtree
    ! ============================================================================
    !> Subroutine global_find_TINfaces_kdtree
    !! Subroutine global_find_TINfaces_kdtree determines the indexes of the faces containing a list of points.
    !<
    ! ============================================================================
    subroutine global_find_TINfaces_kdtree( nb, listx, listy, maxneigh, fce )

        ! Parameters Declaration
        logical :: inside

        integer :: k, kn, fid, nb, maxneigh, ffid
        integer, dimension( nb ) :: lfce, fce

        real( tkind ) :: maxdist, pt( 2 )
        real( tkind ),dimension( :, : ), allocatable :: Fdata

        real( tkind ),dimension( nb ) :: listx, listy

        type(kdtree2), pointer :: Ftree
        type(kdtree2_result), dimension( maxneigh ) :: FRslt

        allocate( Fdata( 2, Tif_Gproc ) )

        ! Build the kd-tree
        do kn = 1, Tif_Gproc
            k = LF_tin( kn )
            Fdata( 1, kn ) =  tinf_centroid( k, 1 )
            Fdata( 2, kn ) =  tinf_centroid( k, 2 )
        enddo
        Ftree => kdtree2_create( Fdata, sort = .true., rearrange = .true. )

        fce = -1
        lfce = -1
        maxdist = (G_tin%maximal_dist + 1.0_8)**2

        ! For each point check the nearest face id
        do k = 1, nb
            inside = .false.
            pt( 1 ) = listx( k )
            pt( 2 ) = listy( k )
            call kdtree2_n_nearest(Ftree, pt, nn=maxneigh, results=FRslt)
                lp: do kn = 1, maxneigh
                    if( FRslt( kn )%idx > 0 .and. FRslt( kn )%idx <= Tif_Gproc )then
                        fid = FRslt( kn )%idx
                        ffid = LF_tin( fid )
                        if( ffid < 1 .or. ffid > G_tin%faces_nb )then
                            print*,'Something went wrong when locating face for flow walker'
                            stop
                        endif
                        ! Check is the point is in the picked face number
                        call check_point_TIN_face( ffid, listx( k ), listy( k ), inside )
                        if( inside )then
                            lfce( k ) = ffid
                            exit lp
                        endif
                    endif
                enddo lp
        enddo

        ! Destroy kd-tree
        call kdtree2_destroy(Ftree)
        deallocate( Fdata )

        call mpi_allreduce( lfce, fce, nb, int_type, max_type, SPModel_comm_world, ierr )

        return

    end subroutine global_find_TINfaces_kdtree
    ! ============================================================================
    !> Subroutine find_TINfaces_kdtree
    !! Subroutine find_TINfaces_kdtree determines the indexes of the faces containing a list of points.
    !<
    ! ============================================================================
    subroutine find_TINfaces_kdtree( nb, listx, listy, maxneigh, fce )

        ! Parameters Declaration
        logical :: inside

        integer :: k, kn, fid, nb, maxneigh, ffid
        integer, dimension( nb ) :: fce

        real( tkind ) :: pt( 2 )
        real( tkind ),dimension( :, : ), allocatable :: Fdata

        real( tkind ),dimension( nb ) :: listx, listy

        type(kdtree2), pointer :: Ftree
        type(kdtree2_result), dimension( maxneigh ) :: FRslt

        allocate( Fdata( 2, Tif_Gproc ) )

        ! Create the kd-tree
        do kn = 1, Tif_Gproc
            k = LF_tin( kn )
            Fdata( 1, kn ) =  tinf_centroid( k, 1 )
            Fdata( 2, kn ) =  tinf_centroid( k, 2 )
        enddo
        Ftree => kdtree2_create( Fdata, sort = .true., rearrange = .true. )

        ! For each points find the nearest TIN face
        fce = -1
        do k = 1, nb
            inside = .false.
            pt( 1 ) = listx( k )
            pt( 2 ) = listy( k )
            call kdtree2_n_nearest(Ftree, pt, nn=maxneigh, results=FRslt)

            locate_point_in_face: do kn = 1, maxneigh
                if( FRslt( kn )%idx > 0 .and. FRslt( kn )%idx <= Tif_Gproc )then
                    fid = FRslt( kn )%idx
                    ffid = LF_tin( fid )
                    if( ffid < 1 .or. ffid > G_tin%faces_nb )then
                        print*,'Something went wrong when locating face for flow walker.'
                        stop
                    endif
                    ! Is the point within the considered face
                    call check_point_TIN_face( ffid, listx( k ), listy( k ), inside )

                    if( inside )then
                        fce( k ) = ffid
                        exit locate_point_in_face
                    endif
                endif
            enddo locate_point_in_face

        enddo

        ! Destroy kd-tree
        call kdtree2_destroy(Ftree)
        deallocate( Fdata )

        return

    end subroutine find_TINfaces_kdtree
    ! ============================================================================
    !> Subroutine check_point_TIN_face()
    !!
    !! Subroutine check_point_TIN_face finds the face containing the flow walker.
    !! \param facenb, x0, y0, inside
    !<
    ! ============================================================================
    subroutine check_point_TIN_face( facenb, x0, y0, inside )

        logical :: inside

        integer :: facenb, k, m, nid( 3 )
        real( tkind ) :: x( 3 ), y( 3 ), x0, y0

        inside = .false.

        nid = tinf_nids( facenb, 1:3 )

        ! Get face points
        x( 1 ) = tinv_coord( nid( 1 ),  1 )
        y( 1 ) = tinv_coord( nid( 1 ),  2 )
        x( 2 ) = tinv_coord( nid( 2 ),  1 )
        y( 2 ) = tinv_coord( nid( 2 ),  2 )
        x( 3 ) = tinv_coord( nid( 3 ),  1 )
        y( 3 ) = tinv_coord( nid( 3 ),  2 )

        call is_point_in_triangle( x0, y0, x, y, k, m )

        if( k >= 0 )inside = .true.

        return

    end subroutine check_point_TIN_face
    ! ============================================================================
    !> Subroutine is_point_in_triangle
    !! Subroutine is_point_in_triangle given a triangle check whether or not the point is
    !! inside the polygon.
    !! \param x0, y0, x, y, l, m
    !<
    ! ============================================================================
    subroutine is_point_in_triangle( xa, ya, xb, yb, l, m )

        integer, intent( out ) :: l, m

        real( tkind ), intent( in ) :: xa, ya, xb(3), yb(3)
        real( tkind ) :: det0, det1, det2

        l = -1
        m = 0

        det0 = ( xb( 2 ) - xb( 1 ) )*( ya - yb( 1 ) ) - ( yb( 2 ) - yb( 1 ) )*( xa - xb( 1 ) )
        det1 = ( xb( 3 ) - xb( 2 ) )*( ya - yb( 2 ) ) - ( yb( 3 ) - yb( 2 ) )*( xa - xb( 2 ) )
        det2 = ( xb( 1 ) - xb( 3 ) )*( ya - yb( 3 ) ) - ( yb( 1 ) - yb( 3 ) )*( xa - xb( 3 ) )

        if( det0 >= 0 .and. det1 >= 0 .and. det2 >= 0 )then
            l = 1
        elseif( det0 <= 0 .and. det1 <= 0 .and. det2 <= 0 )then
            l = 1
        endif

        return

    end subroutine is_point_in_triangle
    ! ============================================================================
    !> Subroutine check_point_TIN_face_containing_point
    !! Subroutine check_point_TIN_face_containing_point is used to determine the face ID containing the point.
    !<
    ! ============================================================================
    subroutine check_point_TIN_face_containing_point( x, y, fcS )

        integer :: Xid, Yid, fcS, id
        real( tkind ) :: x, y

        ! Define the propable row and column number
        Xid = int( ( x - strat_xo ) / strat_dx ) + 1
        Yid = int( ( y - strat_yo ) / strat_dx ) + 1

        ! Define the face id accordingly
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
            fcS = ( strat_X - 1 )* ( Yid - 1 ) + Xid
        endif

        ! Check that picked face is containing the point
        id = gf_points( fcS, 1 )
        if( x - gv_coord( id, 1 ) > strat_dx .or. y - gv_coord( id, 2 ) > strat_dx ) then
            print*,'Something went wrong when looking for the stratal face containing the TIN node'
            stop
        endif

        return

    end subroutine check_point_TIN_face_containing_point
    ! ============================================================================
    !> Subroutine bilinear_interp
    !! Subroutine bilinear_interp performs bilinear interpolation.
    !<
    ! ============================================================================
    subroutine bilinear_interp( x, y, x1, x2, y1, y2, f1, f2, f3, f4, res )

        real( tkind ) :: x, y, x1, x2, y1, y2
        real( tkind ), dimension( 2 ) :: f1, f2, f3, f4, res

        res( 1 ) = -1.e6_8
        res( 2 ) = -100.0_8
        res( 1 ) = f1( 1 ) * ( x2 - x )*( y2 - y ) + f2( 1 ) * ( x - x1 )*( y2 - y ) + &
            f3( 1 ) * ( x2 - x )*( y - y1 ) + f4( 1 ) * ( x - x1 )*( y - y1 )

        res( 2 ) = f1( 2 ) * ( x2 - x )*( y2 - y ) + f2( 2 ) * ( x - x1 )*( y2 - y ) + &
            f3( 2 ) * ( x2 - x )*( y - y1 ) + f4( 2 ) * ( x - x1 )*( y - y1 )

        res( 1 ) = res( 1 ) / ( strat_dx * strat_dx )
        res( 2 ) = res( 2 ) / ( strat_dx * strat_dx )

        return

    end subroutine bilinear_interp
    ! ============================================================================

end module TIN_function

