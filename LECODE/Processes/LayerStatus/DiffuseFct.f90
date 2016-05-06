! ============================================================================
! Name        : DiffuseFct.f90
! Author      : tristan salles
! Created on: March 11, 2013
! Copyright (C) 2013 CSIRO
! ============================================================================
!> \file DiffuseFct.f90
!!
!! DiffuseFct set of functions used for slope diffusion.
!!
!! The file contains the following subroutine:
!! \arg \c build_current_layer_diffusion_grid
!! \arg \c sedimentredistribution
!! \arg \c maxHtopo
!! \arg \c depositedfraction
!! \arg \c get_available_space_remaining
!<
! ============================================================================
module diffuse_fct

    use mpidata
    use mod_diffuse
    use strata_data
    use param_data
    use forces_data
    use TIN_function
    use sediment_data
    use time_data
    use mod_compaction

    implicit none
    public

    real( tkind ), dimension(:), allocatable :: bh, th
    real( tkind ), dimension( : ), allocatable :: tdifp, bdifp
    real( tkind ), dimension( : ), allocatable :: difps1, difps2

contains

    ! ============================================================================
    !> Subroutine build_current_layer_diffusion_grid
    !! Subroutine build_current_layer_diffusion_grid determines nodes with deposition.
    !<
    ! ============================================================================
    subroutine build_current_layer_diffusion_grid

        integer :: i, ks, lid, lid1, lid2, lb, ld, vr, vl
        integer :: req, req2, colW

        integer, dimension( mpi_status_size ) :: stat1, stat2

        real( tkind ) :: dh, ddc, xmax, ymax

        if( .not. allocated( bh ) )then
            if( strat_X > strat_Y )then
                allocate( bh( strat_Y ) )
                allocate( th( strat_Y ) )
            else
                allocate( bh( strat_X ) )
                allocate( th( strat_X ) )
            endif
        endif
        colW = strat_X
        if( strat_X > strat_Y ) colW = strat_Y

        ! Update local diffusion grid top elevation
        do i = 1 , DL_nodesNb
            lid = diff_ID( i )
            if( lid > 0 )then
                lb = locv_NbLay( lid )
                ld = strat_layID( lid, lb )
                diff_coord( i, 3 ) = strat_zelev( lid, lb )
                dh = 0.0_8
                if(  L_strat%layersID == ld .and. strat_thick( lid, lb ) > 0.0_8 )then
                    dh = strat_thick( lid, lb )

                    do ks = 1, totgrn
                        dstart( lid, ks ) = 0.0_8
                        ddc = strat_sedh( lid, lb, ks ) - top_sedprev( lid, ks )
                        if( ddc > 0.0_8 )then
                            dstart( lid, ks ) = ddc / dble( diff_nb )
                            top_sedh( lid, ks ) = top_sedprev( lid ,ks )
                        else
                            top_sedh( lid, ks ) = strat_sedh( lid, lb, ks )
                            top_sedprev( lid ,ks ) = top_sedh( lid, ks )
                        endif
                    enddo
                    strat_thick( lid, lb ) = 0.0_8
                    do ks = 1, totgrn
                        strat_thick( lid, lb ) = strat_thick( lid, lb ) + top_sedh( lid, ks )
                    enddo
                    diff_coord( i, 3 ) = diff_coord( i, 3 ) - ( dh - strat_thick(  lid ,  lb  ) )
                else
                    top_sedprev( lid, 1:totgrn ) = 0.0_8
                    top_sedh( lid, 1:totgrn ) = 0.0_8
                endif
            endif
        enddo

        ! Update halo elevation
        th = 0.0_8
        bh = 0.0_8
        if( diff_neighbors( 1 ) >= 0 )then
            vl = 2 * ( colW+2 ) + 2
            vr = 3 * ( colW+2 ) - 1
            call mpi_isend( diff_coord( vl:vr, 3 ), &
                colW, dbl_type, diff_neighbors( 1 ), 121, &
                SPModel_comm_world, req, ierr )
            call mpi_request_free( req, ierr )
            vl = 2
            vr = colW + 2 - 1
            call mpi_irecv( bh, colW, dbl_type, diff_neighbors( 1 ), 131, &
                SPModel_comm_world,  req2, ierr )
            call mpi_wait( req2, stat2, ierr )
            diff_coord( vl:vr, 3 ) = bh( 1:colW )
        endif

        if( diff_neighbors( 2 ) >= 0 )then
            vl = DL_nodesNb - 3 * ( colW+2 ) + 2
            vr = DL_nodesNb - 2 * ( colW+2 ) - 1
            call mpi_isend( diff_coord( vl:vr, 3 ), colW, dbl_type, &
                diff_neighbors( 2 ), 131, SPModel_comm_world, req2, ierr )
            call mpi_request_free( req2, ierr )
            vl = DL_nodesNb - ( colW+2 ) + 2
            vr = DL_nodesNb - 1
            call mpi_irecv( th, colW, dbl_type, diff_neighbors( 2 ), 121, &
                SPModel_comm_world,  req, ierr )
            call mpi_wait( req, stat1, ierr )
            diff_coord( vl:vr, 3 ) = th( 1:colW )
        endif

        ! Update borders
        xmax = gv_coord( G_strat%nodes_nb, 1 )
        ymax = gv_coord( G_strat%nodes_nb, 2 )
        do i = 1 , DL_nodesNb
            if( diff_coord( i, 1 ) < strat_xo .or. diff_coord( i, 2 ) < strat_yo .or. &
                diff_coord( i, 1 ) > xmax .or. diff_coord( i, 2 ) > ymax )then
                lid = diff_ID( i )
                if( lid == -1 .and. diff_coord( i, 3 ) < 7.e4_8 )then
                    diff_coord( i, 3 ) = -1.0e8_8
                elseif( lid == 0 )then
                    lid = diff_ref( i )

                    ! South slope
                    if( lid == -2 )then
                        if( colW == strat_X )then
                            lid1 = i + strat_X + 2
                            lid2 = lid1 + strat_X + 2
                        else
                            lid1 = i + 1
                            lid2 = lid1 + 1
                        endif
                        dh = diff_coord( lid2, 3 ) - diff_coord( lid1, 3 )
                        diff_coord( i, 3 ) = diff_coord( lid1, 3 ) - dh

                    ! South flat
                    elseif( lid == -5 )then
                        lid1 = i + 1
                        if( colW == strat_X ) lid1 = i + strat_X + 2
                        diff_coord( i, 3 ) = diff_coord( lid1, 3 )

                    ! North slope
                    elseif( lid == -1 )then
                        if( colW == strat_X )then
                            lid1 = i - strat_X - 2
                            lid2 = lid1 - strat_X - 2
                        else
                            lid1 = i - 1
                            lid2 = lid1 - 1
                        endif
                        dh = diff_coord( lid2, 3 ) - diff_coord( lid1, 3 )
                        diff_coord( i ,3 ) = diff_coord( lid1, 3 ) - dh

                    ! North flat
                    elseif( lid == -6 )then
                        lid1 = i - 1
                        if( colW == strat_X ) lid1 = i - strat_X - 2
                        diff_coord( i ,3 ) = diff_coord( lid1, 3 )

                    ! West slope
                    elseif( lid == -3 )then
                        if( colW == strat_X )then
                            lid1 = i + 1
                            lid2 = lid1 + 1
                        else
                            lid1 = i + strat_Y + 2
                            lid2 = lid1 + strat_Y + 2
                        endif
                        dh = diff_coord( lid2, 3 ) - diff_coord( lid1, 3 )
                        diff_coord( i ,3 ) = diff_coord( lid1, 3 ) - dh

                    ! West flat
                    elseif( lid == -7 )then
                        lid1 = i + strat_Y + 2
                        if( colW == strat_X ) lid1 = i + 1
                        diff_coord( i ,3 ) = diff_coord( lid1, 3 )

                    ! East slope
                    elseif( lid == -4 )then
                        if( colW == strat_X )then
                            lid1 = i - 1
                            lid2 = lid1 - 1
                        else
                            lid1 = i - strat_Y - 2
                            lid2 = lid1 - strat_Y - 2
                        endif
                        dh = diff_coord( lid2, 3 ) - diff_coord( lid1, 3 )
                        diff_coord( i ,3 ) = diff_coord( lid1, 3 ) - dh

                    ! East flat
                    elseif( lid == -8 )then
                        lid1 = i - strat_Y - 2
                        if( colW == strat_X ) lid1 = i - 1
                        diff_coord( i ,3 ) = diff_coord( lid1, 3 )
                    endif
                endif

            endif
        enddo

        return

    end subroutine build_current_layer_diffusion_grid
    ! ============================================================================
    !> Subroutine find_maximum_elevation
    !! Subroutine find_maximum_elevation determines what the maximum topographic height
    !!  that each grid node can support whilst being stable and deposit the corresponding elevation.
    !! \param ks
    !<
    ! ============================================================================
    subroutine find_maximum_elevation( ks )

        integer :: k, locid, ks, i, lid

        real( tkind ) :: dist, elev, topnew, toph( DL_nodesNb )

        ! Get the maximum topographic elevation for each node
        do i = 1, DL_nodesNb
            lid = diff_ID( i )
            toph( i ) = diff_coord( i, 3 )
            if( lid > 0 )then
                if( difo( i ) > 0.0_8 )then

                    ! Loop over neighboring nodes
                    do k = 1, 8
                        locid = diff_ngb( i, k )

                        if( locid > 0 .and. locid <= DL_nodesNb )then

                            elev = diff_coord( locid, 3 )
                            dist = ( diff_coord( locid, 2 ) - diff_coord( i, 2 ) )**2
                            dist = dist + ( diff_coord( locid ,1 ) - diff_coord( i, 1 ) )**2
                            dist = sqrt( dist )

                            ! In case node is below sea level
                            if( elev < gsea%actual_sea )then
                                topnew = elev + sediment( ks )%slp_marine * dist
                                if( topnew > gsea%actual_sea + tor .and. sediment( ks )%slp_marine > 0.0_8 )&
                                    topnew = gsea%actual_sea + sediment( ks )%slp_aerial * &
                                    ( dist + ( elev - gsea%actual_sea ) / sediment( ks )%slp_marine )

                            ! Otherwise node is above water
                            else
                                topnew = elev + sediment( ks )%slp_aerial * dist
                            endif

                            topmax( lid ) = min( topmax( lid ), topnew )
                        endif
                    enddo
                endif
            endif
        enddo

        ! Fraction of buffer that will be deposited at current node
        do i = 1, DL_nodesNb
            lid = diff_ID( i )
            if( lid > 0 )then
                if( difo( i ) > 0.0_8 )then
                    if( topmax( lid ) - diff_coord( i, 3 ) > 0.0_8 )then
                        topnew = diff_coord( i, 3 ) + difo( i )
                        if( topnew <= topmax( lid ) )then
                            cdif( lid ) = 1.0_8
                            toph( i ) = topnew
                        else
                            cdif( lid ) = ( topmax( lid ) - diff_coord( i, 3 ) )
                            cdif( lid ) = cdif( lid ) / ( topnew - diff_coord( i, 3 ) )
                            toph( i ) = topmax( lid )
                        endif
                    endif
                    depo( lid, ks ) = depo( lid, ks ) +  difo( i ) * cdif( lid )
                endif
            endif
        enddo

        ! Update local diffusion grid top elevation
        diff_coord( 1:DL_nodesNb,3 ) = toph( 1:DL_nodesNb )

        return

    end subroutine find_maximum_elevation
    ! ============================================================================
    !> Subroutine get_available_space_remaining
    !! Subroutine get_available_space_remaining determines the space available for diffusion around a specific node.
    !! \param vp, spacsum, nconc, spac, down, ndown
    !<
    ! ============================================================================
    subroutine get_available_space_remaining( vd, spacsum, spac, down, ndown )

        integer :: k, locid, ndown, vd
        integer, dimension( 8 ) :: down

        real( tkind ) :: sp, dist, spacsum, slop
        real( tkind ), dimension( 8 ) :: spac

        ! Loop over neighboring cell
        do k = 1, 8
            locid = diff_ngb( vd, k )
            sp = diff_coord( vd, 3 ) - diff_coord( locid, 3 )
            dist = ( diff_coord( locid, 2 ) - diff_coord( vd, 2 ) )**2
            dist = dist + ( diff_coord( locid, 1 ) - diff_coord( vd, 1 ) )**2
            dist = sqrt( dist )

            ! If there is some space
            if( sp > 0 )then
                slop = sp / dist
                if( slop >= minimum_slp )then
                    ndown = ndown + 1
                    down( k ) = 1
                    spac( k ) = sp
                    spacsum = spacsum + sp
                endif
            endif

        enddo

        return

    end subroutine get_available_space_remaining
    ! ============================================================================
    !> Subroutine finalise_top_sediment_layer_diffusion
    !! Subroutine finalise_top_sediment_layer_diffusion updates top sediment layers information as well as the global stratal
    !! on the strata grid and diffusion grid.
    !<
    ! ============================================================================
    subroutine finalise_top_sediment_layer_diffusion

        integer :: k, ks, laynb, toplay, gid

        real( tkind ) :: th, nh

        ! Update the stratigraphic layer locally
        proc_elev = -1.0e6_8

        do k = 1, L_strat%nodes_nb
            gid = locv_gid( k )
            laynb = locv_NbLay( k )
            toplay = strat_layID(  k ,  laynb  )
            th = 0.0_8

            ! Get the updated thickness
            do ks = 1, totgrn
                top_sedh( k ,ks ) = top_sedprev( k ,ks ) + depo( k, ks )
                top_sedprev( k ,ks ) = top_sedh( k ,ks )
                th =  th + top_sedh( k ,ks )
            enddo

            ! In case there is something to deposit
            if( th > tor )then

                ! If the top layer ID corresponds to the current layer ID
                if( L_strat%layersID == toplay )then
                    th = 0.0_8
                    ! Update layer sediment class thicknesses
                    do ks = 1, totgrn
                        th = th + top_sedh( k, ks )
                        strat_sedh( k, laynb, ks ) = top_sedh( k, ks )
                        if( top_sedh( k, ks ) > 0.0_8 .and. gporo%compaction ) &
                            strat_porosity(  k, laynb, ks ) = porosity( ks, 1 )
                    enddo
                    ! Update layer parameters
                    strat_thick(  k, laynb ) = th
                    nh = strat_zelev(  k , laynb - 1  ) + strat_thick(  k, laynb )
                    strat_zelev(  k, laynb ) = nh
                    proc_elev( gid ) = nh
                    strat_fac( k, laynb ) = nh - gsea%actual_sea
                ! Otherwise create a new layer
                else
                    locv_NbLay( k ) = locv_NbLay( k ) + 1
                    laynb = locv_NbLay( k )
                    if( locv_NbLay( k ) > G_strat%total_layer )print*,'Something went wrong went adding layer after diffusion'
                    if( locv_NbLay( k ) > G_strat%total_layer ) stop
                    strat_layID( k, laynb ) = L_strat%layersID
                    strat_hardness(  k, laynb ) = 1.0_8
                    strat_porosity(  k, laynb, 1:totgrn ) = 0.0_8
                    toplay = locv_NbLay( k )
                    ! Update layer sediment class thicknesses
                    th = 0.0_8
                    do ks = 1, totgrn
                        th = th + top_sedh( k, ks )
                        strat_sedh( k, laynb, ks ) = top_sedh( k, ks )
                        if( top_sedh( k, ks ) > 0.0_8 .and. gporo%compaction ) &
                            strat_porosity(  k, laynb, ks ) = porosity( ks, 1 )
                    enddo
                    ! Update layer parameters
                    strat_thick( k, laynb ) = th
                    nh = strat_zelev( k, laynb - 1 ) + strat_thick( k, laynb )
                    strat_zelev( k, laynb ) = nh
                    proc_elev( gid ) = nh
                    strat_fac( k, laynb ) = nh - gsea%actual_sea
                endif

            ! In case the layer is empty
            else
                ! If a layer was there delete it
                if( L_strat%layersID == toplay )then
                    strat_thick(  k, laynb ) = 0.0_8
                    strat_porosity(  k, laynb, 1:totgrn ) = 0.0_8
                    strat_layID( k, laynb ) = -1
                    locv_NbLay( k ) = locv_NbLay( k ) - 1
                    proc_elev( gid ) = strat_zelev(  k, locv_NbLay( k ) )
                else
                    proc_elev( gid ) = strat_zelev(  k, laynb )
                endif
            endif
        enddo

        ! Update elevation
        call mpi_allreduce( proc_elev, elev_record , G_strat%nodes_nb, &
            dbl_type, max_type, SPModel_comm_world, ierr )

        ! Update diffusion grid local elevation
        do k = 1, DL_nodesNb
            gid = diff_ID( k )
            if( gid > 0 )then
                th = 0.0_8
                do ks = 1, totgrn
                    th =  th + top_sedh( gid, ks )
                enddo
                diff_coord( k, 3 ) = diff_coord( k, 3 ) + th
            endif
        enddo

        return

    end subroutine finalise_top_sediment_layer_diffusion
    ! ============================================================================

end module diffuse_fct

