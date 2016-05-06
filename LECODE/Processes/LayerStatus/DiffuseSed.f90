! ============================================================================
! Name        : DiffuseSed.f90
! Author      : tristan salles
! Created on: March 11, 2013
! Copyright (C) 2013 CSIRO
! ============================================================================
!> \file DiffuseSed.f90
!!
!! DiffuseSed simulates sediment repartition over topography according to slope.
!!
!! The file contains the following subroutine:
!! \arg \c DiffuseSed
!! \arg \c redistribute_sediment_by_diffusion
!<
! ============================================================================
module mod_seddiffusion

    use flux_data
    use mpidata
    use strata_data
    use param_data
    use mod_diffuse
    use TIN_function
    use mod_direction
    use sediment_data
    use diffuse_fct
    use TIN_surface

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine perform_top_sediment_layer_diffusion
    !! Subroutine perform_top_sediment_layer_diffusion mimics slope failure by geometric redistribution of sediment
    !! that does not meet a slope criterion which is function of grain size.
    !<
    ! ============================================================================
    subroutine perform_top_sediment_layer_diffusion( dotin )

        integer :: dotin

        dstart( :, : ) = 0.0_8

        ! Update grid and face values
        call build_current_layer_diffusion_grid

        ! Distribute sediment according to slope criteria
        call redistribute_sediment_by_diffusion

        ! Redistribute sediment into active cells
        call finalise_top_sediment_layer_diffusion

        ! Update the TIN grid elevation
        if( dotin == 1 ) call generate_TIN_grid

        ! Find the global pathways
        if( dotin == 1 ) call find_global_pathways_G8_Paik

        return

    end subroutine perform_top_sediment_layer_diffusion
    ! ============================================================================
    !> Subroutine redistribute_sediment_by_diffusion
    !!
    !! Subroutine redistribute_sediment_by_diffusion redistributes sediment geometrically using slope criteria.
    !! Slope criteria are functions of grain size and sediment is diffused to downstream nodes.
    !!
    !! Diffusion of sediment is subdivided into diffusion cycles.
    !! In each cycle all excess material is diffused iteratively until everything is deposited.
    !!
    !! Diffusion only takes place to connected nodes.
    !! All material that has been deposited at unconnected nodes is not checked for slope.
    !!
    !! List of variables:
    !! \arg \c idown - size of kernel
    !! \arg \c depo - redistributed sediment
    !! \arg \c difo - sediment to be deposited at current node
    !! \arg \c difp - sediment to be deposited to downstream nodes
    !! \arg \c dstart - fraction of sediment used when starting diffusion
    !! \arg \c cdif - fraction of sediment to be deposited at current node
    !! \arg \c movout - sediment moved out of simulation area
    !! \arg \c spac - accomodation space at downstream nodes
    !! \arg \c icrow - row position change of downstream nodes in kernel
    !! \arg \c iccol - col position change of downstream nodes in kernel
    !!
    !<
    ! ============================================================================
    subroutine redistribute_sediment_by_diffusion

        logical :: inX

        integer :: k, kd, ke, ks, vp, iter, locid, ndown, vd, vn
        integer :: request, request2, request3, request4, vl, vr, p, colW

        integer, dimension( 8 ) :: down
        integer, dimension( mpi_status_size ) :: status1, status2
        integer, dimension( mpi_status_size ) :: status3, status4

        real( tkind ) :: difmax, difsed, fc
        real( tkind ) :: spacsum, cc

        real( tkind ), dimension( 8 ) :: spac

        ! Allocate diffusion mpi passing arrays
        if( .not. allocated( tdifp ) )then
            if( strat_X > strat_Y )then
                allocate( tdifp( strat_Y ) )
                allocate( bdifp( strat_Y ) )
                allocate( difps1( strat_Y ) )
                allocate( difps2( strat_Y ) )
            else
                allocate( tdifp( strat_X ) )
                allocate( bdifp( strat_X ) )
                allocate( difps1( strat_X ) )
                allocate( difps2( strat_X ) )
            endif
        endif

        inX = .false.
        colW = strat_X
        if( strat_X > strat_Y ) colW = strat_Y
        if( strat_X > strat_Y ) inX = .true.

        ! Initialize deposition arrays
        depo( :, : ) = 0.0_8

        ! Deposit sediment in diff_nb increments
        do kd = 1, diff_nb

            ! Deposit grain sizes from maximum to minimum sediment grain weight
            do ke = 1, totgrn
                ks = dep_order( ke )

                ! Initialise sediments
                difp = 0.0_8
                difo = 0.0_8
                difps1 = 0.0_8
                difps2 = 0.0_8
                difsed = 0.0_8

                ! Update diffusion based on initial sediment class thickness
                do vn = 1 , DL_nodesNb
                    vp = diff_ID( vn )
                    if( vp > 0 ) difo( vn ) = dstart( vp, ks )
                enddo

                ! Move sediment to downstream nodes as long as there is sediment in the buffer.
                difmax = 1.0e6_8
                iter = 0
                sediment_diffusion: do

                    if( difmax <= diff_res .or. iter >= max_it_cyc ) exit sediment_diffusion
                    difsed = 0.0_8
                    difmax = 0.0_8
                    iter = iter + 1
                    topmax = toplimit

                    ! Determine the fraction of buffer that will be deposited at current node
                    cdif = 0.0_8

                    ! Get maximum elevation at current nodes to ensure sediment stability
                    call find_maximum_elevation( ks )

                    bh = 0.0_8
                    th = 0.0_8

                    ! Update partition halo & border elevation based on new elevation
                    if( diff_neighbors( 1 ) >= 0 )then
                        vl = 2 * ( colW+2 ) + 2
                        vr = 3 * ( colW+2 ) - 1
                        th( 1:colW ) = diff_coord( vl:vr, 3 )
                        call mpi_isend( th, colW, dbl_type, diff_neighbors( 1 ), 11, &
                            SPModel_comm_world, request, ierr )
                        call mpi_request_free( request, ierr )
                        call mpi_irecv( bh, colW, dbl_type, diff_neighbors( 1 ), 13, &
                            SPModel_comm_world,  request2, ierr )
                        call mpi_wait( request2, status2, ierr )
                        vl = 2
                        vr = colW + 2 - 1
                        diff_coord( vl:vr, 3 ) = bh( 1:colW )
                    endif

                    if( diff_neighbors( 2 ) >= 0 )then
                        vl = DL_nodesNb - 3 * ( colW+2 ) + 2
                        vr = DL_nodesNb - 2 * ( colW+2 ) - 1
                        bh( 1:colW ) = diff_coord( vl:vr, 3 )
                        call mpi_isend( bh, colW, dbl_type, diff_neighbors( 2 ), 13, &
                            SPModel_comm_world, request2, ierr )
                        call mpi_request_free( request2, ierr )
                        call mpi_irecv( th, colW, dbl_type, diff_neighbors( 2 ), 11, &
                            SPModel_comm_world,  request, ierr )
                        call mpi_wait( request, status1, ierr )
                        vl = DL_nodesNb - ( colW + 2 ) + 2
                        vr = DL_nodesNb - 1
                        diff_coord( vl:vr, 3 ) = th( 1:colW )
                    endif
                    
                    ! Calculate individual and cumulative accomodation space for each downstream node
                    do vd = 1, DL_nodesNb
                        vp = diff_ID( vd )
                        if( vp > 0 )then
                            if( difo( vd ) > 0.0_8 )then
                                if( cdif( vp ) < 1.0_8 )then
                                    spacsum = 0.0_8
                                    ndown = 0
                                    down(:) = 0
                                    spac(:) = 0.0_8

                                    call get_available_space_remaining( vd, spacsum, spac, down, ndown )
                                    cc = difo( vd ) * ( 1.0_8 - cdif( vp ) )
                                    difsed = difsed + cc

                                    ! If enough downstream space is available weight sediment by
                                    ! available space and diffuse sediment to downstream nodes
                                    if( spacsum >= cc )then

                                        ! Loop over neighboring nodes
                                        do k = 1, 8
                                            locid = diff_ngb( vd, k )
                                            if( diff_ID( locid ) > 0 )then
                                                difp( locid ) = difp( locid ) + cc * spac( k ) / spacsum

                                                ! First line in stratal grid
                                                if( ( inX .and. locv_coord( diff_ID( locid ), 1 ) &
                                                    == locv_coord( 1, 1 ) ) .or. ( .not. inX .and. &
                                                    locv_coord( diff_ID( locid ), 2 ) == locv_coord( 1, 2 ) ) )then
                                                    if( vd > 2 * ( colW + 2 ) + 1 .and. vd < 3 * ( colW + 2 ) )then
                                                        p = locid - ( colW + 2 + 1 )
                                                        difps1( p ) = difps1( p ) + cc * spac( k ) / spacsum
                                                    endif
                                                endif

                                                ! Last line in stratal grid
                                                if( ( inX .and. locv_coord( diff_ID( locid ), 1 ) &
                                                    == locv_coord( L_strat%nodes_nb, 1 ) ) .or. ( .not. inX .and. &
                                                    locv_coord( diff_ID( locid ), 2 ) == locv_coord( L_strat%nodes_nb, 2 ) ) )then
                                                    ! Third line from top in diffusion grid
                                                    if( vd > DL_nodesNb - 3 * ( colW + 2 ) + 1 .and. &
                                                        vd < DL_nodesNb - 2 * ( colW + 2 ) )then
                                                        p = locid - ( DL_nodesNb - 2 * ( colW + 2 ) + 1 )
                                                        difps2( p ) = difps2( p ) +  cc * spac( k ) / spacsum
                                                    endif
                                                endif
                                            endif
                                        enddo

                                    ! If not enough downstream space is available weight sediment
                                    ! by available space and distribute rest equally among
                                    ! connected downstream nodes or if no downstream nodes are
                                    ! connected assume that enough downstream space is available
                                    ! and weight sediment by available space
                                    else
                                        if( ndown > 0 )then
                                            fc = ( cc - spacsum )/dble( ndown )

                                            ! Loop over neighboring nodes
                                            do k = 1, 8
                                                locid = diff_ngb( vd ,k )
                                                if( diff_ID( locid ) > 0 )then
                                                    if( down( k ) == 1 )then
                                                        difp( locid ) = difp( locid ) + spac( k ) + fc

                                                        ! First line in stratal grid
                                                        if( ( inX .and. locv_coord( diff_ID( locid ), 1 ) &
                                                            == locv_coord( 1, 1 ) ) .or. ( .not. inX .and. &
                                                            locv_coord( diff_ID( locid ), 2 ) == locv_coord( 1, 2 ) ) )then
                                                            ! Third line in diffusion grid
                                                            if( vd > 2 * ( colW + 2 ) + 1 .and. vd < 3 * ( colW + 2 ) )then
                                                                p = locid - ( colW + 2 + 1 )
                                                                difps1( p ) = difps1( p ) + spac( k ) + fc
                                                            endif
                                                        endif

                                                        ! Last line in stratal grid
                                                        if( ( inX .and. locv_coord( diff_ID( locid ), 1 ) &
                                                            == locv_coord( L_strat%nodes_nb, 1 ) ) .or. ( .not. inX .and. &
                                                            locv_coord( diff_ID( locid ), 2 ) == locv_coord( L_strat%nodes_nb, 2 ))&
                                                            )then
                                                            ! Third line from top in diffusion grid
                                                            if( vd > DL_nodesNb - 3 * ( colW + 2 ) + 1 .and. &
                                                                vd < DL_nodesNb - 2 * ( colW + 2 ) )then
                                                                p = locid - ( DL_nodesNb - 2 * ( colW + 2 ) + 1 )
                                                                difps2( p ) = difps2( p ) + spac( k ) + fc
                                                            endif
                                                        endif
                                                    else
                                                        difp( locid ) = difp( locid ) + spac( k )

                                                        ! First line in stratal grid
                                                        if( ( inX .and. locv_coord( diff_ID( locid ), 1 ) &
                                                            == locv_coord( 1, 1 ) ) .or. ( .not. inX .and. &
                                                            locv_coord( diff_ID( locid ), 2 ) == locv_coord( 1, 2 ) ) )then
                                                            ! Third line in diffusion grid
                                                            if( vd > 2 * ( colW + 2 ) + 1 .and. vd < 3 * ( colW + 2 ) )then
                                                                p = locid - ( colW + 2 + 1 )
                                                                difps1( p ) = difps1( p )  + spac( k )
                                                            endif
                                                        endif

                                                        ! Last line in stratal grid
                                                        if( ( inX .and. locv_coord( diff_ID( locid ), 1 ) &
                                                            == locv_coord( L_strat%nodes_nb, 1 ) ) .or. ( .not. inX .and. &
                                                            locv_coord( diff_ID( locid ), 2 ) == locv_coord( L_strat%nodes_nb, 2 ))&
                                                            )then
                                                            ! Third line from top in diffusion grid
                                                            if( vd > DL_nodesNb - 3 * ( colW + 2 ) + 1 .and. &
                                                                vd < DL_nodesNb - 2 * ( colW + 2 ) )then
                                                                p = locid - ( DL_nodesNb - 2 * ( colW + 2 ) + 1 )
                                                                difps2( p ) = difps2( p ) + spac( k )
                                                            endif
                                                        endif
                                                    endif
                                                endif
                                            enddo
                                        ! If there are no downtream nodes, just smear it all around
                                        ! and hope the next iteration looks after it
                                        else
                                            fc = cc / 8.0_8

                                            ! Loop over neighboring nodes
                                            do k = 1, 8
                                                locid = diff_ngb( vd ,k )
                                                if( diff_ID( locid ) > 0 )then
                                                    difp( locid ) = difp( locid ) + fc

                                                    ! First line in stratal grid
                                                    if( ( inX .and. locv_coord( diff_ID( locid ), 1 ) &
                                                        == locv_coord( 1, 1 ) ) .or. ( .not. inX .and. &
                                                        locv_coord( diff_ID( locid ), 2 ) == locv_coord( 1, 2 ) ) )then
                                                        ! Third line in diffusion grid
                                                        if( vd > 2 * ( colW + 2 ) + 1 .and. vd < 3 * ( colW + 2 ) )then
                                                            p = locid - ( colW + 2 + 1 )
                                                            difps1( p ) = difps1( p ) + fc
                                                        endif
                                                    endif

                                                    ! Last line in stratal grid
                                                    if( ( inX .and. locv_coord( diff_ID( locid ), 1 ) &
                                                        == locv_coord( L_strat%nodes_nb, 1 ) ) .or. ( .not. inX .and. &
                                                        locv_coord( diff_ID( locid ), 2 ) == locv_coord( L_strat%nodes_nb, 2 ))&
                                                        )then
                                                        ! Third line from top in diffusion grid
                                                        if( vd > DL_nodesNb - 3 * ( colW + 2 ) + 1 .and. &
                                                            vd < DL_nodesNb - 2 * ( colW + 2 ) )then
                                                            p = locid - ( DL_nodesNb - 2 * ( colW + 2 ) + 1 )
                                                            difps2( p ) = difps2( p ) + fc
                                                        endif
                                                    endif
                                                endif
                                            enddo
                                        endif
                                    endif
                                endif
                            endif
                        endif
                    enddo

                    ! Update border thickness value to be diffuse
                    bdifp = 0.0_8
                    tdifp = 0.0_8
                    if( diff_neighbors( 1 ) >= 0 )then
                        call mpi_isend( difps1, colW, dbl_type, diff_neighbors( 1 ), &
                            22, SPModel_comm_world, request4, ierr )
                        call mpi_request_free( request4, ierr )
                        call mpi_irecv( bdifp, colW, dbl_type, diff_neighbors( 1 ), 23, &
                            SPModel_comm_world,  request3, ierr )
                        call mpi_wait( request3, status3, ierr )
                    endif

                    if( diff_neighbors( 2 ) >= 0 )then
                        call mpi_isend( difps2, colW, dbl_type, &
                            diff_neighbors( 2 ), 23, SPModel_comm_world, request3, ierr )
                        call mpi_request_free( request3, ierr )
                        call mpi_irecv( tdifp, colW, dbl_type, diff_neighbors( 2 ), 22, &
                            SPModel_comm_world,  request4, ierr )
                        call mpi_wait( request4, status4, ierr )
                    endif

                    do k = 0, colW - 1
                        if( diff_neighbors( 1 ) >= 0 )then
                            vl = ( colW+2 ) + 2
                            difp( vl + k ) = difp( vl + k ) + bdifp( k + 1 )
                        endif
                        if( diff_neighbors( 2 ) >= 0 )then
                            vl = DL_nodesNb - 2 * ( colW+2 ) + 2
                            difp( vl + k ) = difp( vl + k ) + tdifp( k + 1 )
                        endif
                    enddo

                    ! Store sediment still to be diffused in difo for next iteration
                    difo = difp
                    difp = 0.0_8
                    difps1 = 0.0_8
                    difps2 = 0.0_8

                    ! Determine max. diffusion residual
                    call mpi_allreduce( difsed, difmax, 1, dbl_type, sum_type, SPModel_comm_world, ierr )

                enddo sediment_diffusion
            enddo

        enddo

        return

    end subroutine redistribute_sediment_by_diffusion
    ! ============================================================================


end module mod_seddiffusion
