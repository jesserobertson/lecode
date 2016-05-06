module ocean_transport

    use mpidata
    use strata_data
    use forces_data
    use param_data
    use flux_data
    use time_data
    use ocean_data
    use mod_diffuse
    use strata_update
    use sediment_data

    implicit none

    ! Ocean time step for morphological changes
    real( tkind ) :: ocean_time

    ! Number of ocean morphological update steps
    integer, parameter :: ocean_step = 20

    ! Bed roughness
    real( tkind ), parameter :: zo = 0.006_8

    ! von Karman constant
    real( tkind ), parameter :: vonKarman = 0.41_8

    ! Kinematic viscosity
    real( tkind ), parameter :: kinematic_viscosity = 2.e-6_8

    ! Halo borders
    real( tkind ), dimension(:,:), allocatable :: bq, tq

contains

    ! ============================================================================
    !> Subroutine cmpt_ocean_transport
    !! Main function to compute sediment transport by waves and currents.
    !<
    ! ============================================================================
    subroutine cmpt_ocean_transport

        integer :: sgp, nbclass, n, ks, kk, gid, layer_nb, ostp

        real( tkind ) :: water_depth, epsk, time_perc
        real( tkind ), dimension( totgrn ) :: prop, ph, pe

        nbclass = hindcast( active_hindcast )%cnb

        ! Proceed to morphological changes in several steps
        do ostp = 1, ocean_step

            do sgp = 1, nbclass

                ! Incremental time required to simulate morphological evolution
                time_perc = hindcast( active_hindcast )%subgroup( sgp )%perc
                ocean_time = secyear * time_perc * ocean_time_step / ( transport%morpho )
                ocean_time = ocean_time / ocean_step

                ! Find local component of sediment transport using Soulsby-van Rijn formula
                do n = 1, L_strat%nodes_nb

                    ! Get stratal node global ID
                    gid = locv_gid( n )

                    ! Get the water depth
                    water_depth = gsea%actual_sea - elev_record( gid )

                    soulsby_total_load( n, 1:totgrn ) = 0.0_8

                    if( water_depth > 0.01_8 .and. ( ocean_current( sgp, gid, 1 ) > tor .or. &
                        ocean_wave( sgp, gid, 1 ) > tor ) )then

                        ! Update active layer
                        call build_active_nodelayer( ocean_morph_th, n )

                        ! Find top layer number
                        layer_nb = locv_NbLay( n )

                        ! Get the fraction of each grain present in the top layer
                        prop( 1:totgrn ) = 0.0_8
                        if( strat_thick( n, layer_nb ) > tor )then

                            do ks = 1, totgrn
                                prop( ks ) = strat_sedh(  n,  layer_nb ,  ks  ) / strat_thick( n, layer_nb )
                            enddo

                            ! Exposed / hidden probabilities in the top layer
                            do kk = 1, totgrn
                                ph( kk ) = 0.0_8
                                pe( kk ) = 0.0_8
                                ! Get hidden probabilities
                                do ks = 1, totgrn
                                    ph( kk ) = ph( kk ) + prop( ks ) * sediment( ks )%diameter / &
                                        ( sediment( ks )%diameter + sediment( kk )%diameter )
                                enddo
                                if( ph( kk ) < tor )then
                                    print*,'Something went wrong computing Soulsby hidden probabilities.'
                                    print*,ph( kk ), prop(1:totgrn),strat_sedh(  n,  layer_nb ,  1:totgrn  )
                                    stop
                                endif

                                ! Deduce exposed probabilities
                                pe( kk ) =  1.0_8 - ph( kk )
                                if( pe( kk ) < 0.0_8 )then
                                    ph( kk ) = 1.0_8
                                    pe( kk ) = 0.0_8
                                endif

                            enddo

                            ! Equilibrium total-load transport induced by waves and current (m2/s)
                            do ks = 1, totgrn
                                ! Fractional bed composition
                                epsk = ( pe( ks ) / ph( ks ) ) !**(-0.75_8)
                                if( prop( ks ) > 0.0_8 )then
                                    ! Soulsby van Rijn total load transport
                                    call soulsby_vanRijn_total_load_transport( sgp, gid, ks, prop( ks ), &
                                        epsk, strat_sedh(  n,  layer_nb ,  ks  ), soulsby_total_load( n, ks ) )
                                endif
                            enddo

                        endif
                    endif
                enddo

                ! Get total load vector components for diffusion grid nodes
                call update_total_load_vector_component( sgp )

                ! Get yearly variation in sea bed elevation
                call continuity_equation_sediment_transport

            enddo

        enddo

        ! Update local stratigraphic elevation
        proc_elev = -1.0e6_8
        do n = 1, L_strat%nodes_nb
            kk = locv_gid( n )
            layer_nb = locv_NbLay( n )
            proc_elev( kk ) = strat_zelev(  n ,  layer_nb  )
        enddo

        ! Update elevation
        call mpi_allreduce( proc_elev, elev_record , G_strat%nodes_nb, &
            dbl_type, max_type, SPModel_comm_world, ierr )

        return

    end subroutine cmpt_ocean_transport
    ! ============================================================================
    !> Subroutine soulsby_vanRijn_total_load_transport
    !!
    !! Soulsby (1997) developed an analytical expression similar to the formulation with wave-
    !! current interaction proposed by van Rijn (1993). This Soulsby approach evaluates load
    !! and suspended sediment transport.
    !!
    !<
    ! ============================================================================
    subroutine soulsby_vanRijn_total_load_transport( sgp, gid, ks, pks, epsk, hks, qt )

        integer :: sgp, gid, ks, p

        real( tkind ) :: water_depth, relative_density, d_star, hks, zp, max_eroh
        real( tkind ) :: Cd, Asb, Ass, epsk, ucrk, qt, pks, comb_vel, zmin

        ! Get the water depth
        water_depth = gsea%actual_sea - elev_record( gid )

        relative_density = sediment( ks )%density / sea_density

        ! Dimensionless grain size
        d_star = sediment( ks )%diameter * ( ( relative_density - 1 ) * gravity / &
                kinematic_viscosity**2.0_8 )**( 1/3 )

        ! Current drag coefficient
        Cd = ( vonKarman / log( water_depth / zo ) )**2.0_8

        ! Bed load related coefficient
        Asb = 0.005_8 * water_depth * ( sediment( ks )%diameter / water_depth )**1.2_8 / &
            ( ( relative_density - 1 ) * gravity * sediment( ks )%diameter )**1.2_8

        ! Suspended load related coefficient
        Ass = 0.012_8 * sediment( ks )%diameter * d_star**(-0.6_8 ) / &
            ( ( relative_density - 1 ) * gravity * sediment( ks )%diameter )**1.2_8

        ! Critical velocity of incipient motion
        if( sediment( ks )%diameter <= 0.0005_8 )then
            ucrk = 0.19_8 * sediment( ks )%diameter**0.1_8 * log( 4.0_8 * water_depth / &
                sediment( ks )%diameter )
        else
            ucrk = 8.5_8 * sediment( ks )%diameter**0.6_8 * log( 4.0_8 * water_depth / &
                sediment( ks )%diameter )
        endif

        ! Combined wave current velocity
        comb_vel = sqrt( ocean_current( sgp, gid, 1 )**2.0_8 + 0.018 * ocean_wave( sgp, gid, 1 )**2.0_8 / Cd )

        ! Is the combined wave current velocity above threshold of motion
        if( comb_vel > epsk * ucrk )then
            ! Soulsby-van Rijn formula (m2/s)
            qt = ( Asb + Ass ) * ocean_current( sgp, gid, 1 ) * pks * &
                ( comb_vel - epsk * ucrk )**2.4_8
            qt = max( qt, 0.0_8 )

            ! Check sediment availbility
            if( qt * ocean_time / strat_dx > hks ) qt = hks * strat_dx / ocean_time

            ! Check surrounding elevation to prevent formation of hole
            zmin = 1.e8_8
            do p = 1, 8
                if( gv_ngbID( gid, p ) <= 0 )then
                    zp = elev_record( gid )
                else
                    zp = elev_record( gv_ngbID( gid, p ) )
                endif
                zmin = min( zp, zmin )
            enddo
            max_eroh = elev_record( gid ) - zmin

            if( max_eroh <= 0.0_8 )then
                qt = 0.0_8
            else
                ! Morphological sediment thickness limiter
                if( qt * ocean_time / strat_dx > max_eroh ) &
                    qt = max_eroh * strat_dx / ocean_time
                ! Morphological sediment thickness limiter
                if( qt * ocean_time / strat_dx > ocean_morph_th ) &
                     qt = ocean_morph_th * strat_dx / ocean_time
            endif
        else
            qt = 0.0_8
        endif

        return

    end subroutine soulsby_vanRijn_total_load_transport
    ! ============================================================================
    !> Subroutine update_total_load_vector_component
    !!
    !! Get the Soulsby-van Rijn total load transport vector over the diffusion grid to manage boundary
    !! conditions as well as halo cells.
    !!
    !<
    ! ============================================================================
    subroutine update_total_load_vector_component( sgp )

        integer :: i, lid, vr, vl, lid1, sgp, nx, ny, locId
        integer :: req, req2, colW, gid, ks, layer_nb

        integer, dimension( mpi_status_size ) :: stat1, stat2

        real( tkind ) :: porosity, prop

        ! Get the size of the vectors to send to neighbors
        if( .not. allocated( bq ) )then
            if( strat_X > strat_Y )then
                allocate( bq( strat_Y, totgrn ) )
                allocate( tq( strat_Y, totgrn ) )
            else
                allocate( bq( strat_X, totgrn ) )
                allocate( tq( strat_X, totgrn ) )
            endif
        endif
        colW = strat_X
        if( strat_X > strat_Y ) colW = strat_Y

        ! Update total load on each local diffusion grid
        do i = 1 , DL_nodesNb
            lid = diff_ID( i )
            if( lid > 0 ) total_load( i, 1:totgrn ) =  soulsby_total_load( lid, 1:totgrn )
        enddo

        ! Update total load halo nodes
        tq = 0.0_8
        bq = 0.0_8
        if( diff_neighbors( 1 ) >= 0 )then
            vl = 2 * ( colW+2 ) + 2
            vr = 3 * ( colW+2 ) - 1
            call mpi_isend( total_load( vl:vr, 1:totgrn ), colW*totgrn, dbl_type, diff_neighbors( 1 ), &
                121, SPModel_comm_world, req, ierr )
            call mpi_request_free( req, ierr )
            vl = 2
            vr = colW + 2 - 1
            call mpi_irecv( bq, colW*totgrn, dbl_type, diff_neighbors( 1 ), 131, &
                SPModel_comm_world,  req2, ierr )
            call mpi_wait( req2, stat2, ierr )
            total_load( vl:vr, 1:totgrn ) = bq( 1:colW, 1:totgrn )
        endif

        if( diff_neighbors( 2 ) >= 0 )then
            vl = DL_nodesNb - 3 * ( colW+2 ) + 2
            vr = DL_nodesNb - 2 * ( colW+2 ) - 1
            call mpi_isend( total_load( vl:vr, 1:totgrn ), colW*totgrn, dbl_type, diff_neighbors( 2 ), &
                131, SPModel_comm_world, req2, ierr )
            call mpi_request_free( req2, ierr )
            vl = DL_nodesNb - ( colW+2 ) + 2
            vr = DL_nodesNb - 1
            call mpi_irecv( tq, colW*totgrn, dbl_type, diff_neighbors( 2 ), 121, &
                SPModel_comm_world,  req, ierr )
            call mpi_wait( req, stat1, ierr )
            total_load( vl:vr, 1:totgrn ) = tq( 1:colW, 1:totgrn )
        endif

        ! Compute total load transport vector components
        do i = 1 , DL_nodesNb
            lid = diff_ID( i )
            transportX( i, 1:totgrn ) = 0.0_8
            transportY( i, 1:totgrn ) = 0.0_8
            ! If nodes is within the processor partitionned area
            if( lid > 0 )then
                gid = locv_gid( lid )
                do ks = 1, totgrn
                    ! The sediment flux is considered approximately aligned with the current, although wave effects can cause
                    ! deviations of the sediment transport direction by up to 15 degree relative to the current direction (Soulsby, 1997).
                    ! Thus, the sediment transport vector, in two coordinates, can be written as:
                    if( ocean_current( sgp, gid, 1 ) > 0.0_8 )then
                        transportX( i, ks ) = total_load( i, ks ) * cos( ocean_current( sgp, gid, 2 ) )
                        transportY( i, ks ) = total_load( i, ks ) * sin( ocean_current( sgp, gid, 2 ) )
                    ! In case there is no current use wave direction
                    elseif( ocean_wave( sgp, gid, 1 ) > 0.0_8 )then
                        transportX( i, ks ) = total_load( i, ks ) * cos( ocean_wave( sgp, gid, 2 ) )
                        transportY( i, ks ) = total_load( i, ks ) * sin( ocean_wave( sgp, gid, 2 ) )
                    endif
                enddo

            ! If we are on a halo border
            elseif( diff_coord( i, 1 ) >= strat_xo .and. diff_coord( i, 2 ) >= strat_yo .and. &
                diff_coord( i, 1 ) <= gv_coord( G_strat%nodes_nb, 1 ) .and. &
                diff_coord( i, 2 ) <= gv_coord( G_strat%nodes_nb, 2 ) )then

                ! Find the global ID
                nx = int( ( diff_coord( i, 1 ) - strat_xo ) / strat_dx ) + 1
                ny = int( ( diff_coord( i, 2 ) - strat_yo ) / strat_dx )
                gid = ny * strat_X + nx

                do ks = 1, totgrn
                    ! The sediment flux is considered approximately aligned with the current, although wave effects can cause
                    ! deviations of the sediment transport direction by up to 15 degree relative to the current direction (Soulsby, 1997).
                    ! Thus, the sediment transport vector, in two coordinates, can be written as:
                    if( ocean_current( sgp, gid, 1 ) > 0.0_8 )then
                        transportX( i, ks ) = total_load( i, ks ) * cos( ocean_current( sgp, gid, 2 ) )
                        transportY( i, ks ) = total_load( i, ks ) * sin( ocean_current( sgp, gid, 2 ) )
                    ! In case there is no current use wave direction
                    elseif( ocean_wave( sgp, gid, 1 ) > 0.0_8 )then
                        transportX( i, ks ) = total_load( i, ks ) * cos( ocean_wave( sgp, gid, 2 ) )
                        transportY( i, ks ) = total_load( i, ks ) * sin( ocean_wave( sgp, gid, 2 ) )
                    endif
                enddo

            endif
        enddo

        ! Check sediment availability
        do  i = 1, DL_nodesNb
            locId = diff_ID( i )
            if( locId > 0 )then
                ! Get top layer number
                layer_nb = locv_NbLay( locId )
                ! Get top layer porosity
                porosity = 0.0_8
                prop = 0.0_8
                if( strat_thick(  locId,  layer_nb  ) > 0.0_8 )then
                    do ks = 1, totgrn
                        prop = strat_sedh(  locId,  layer_nb,  ks  ) / strat_thick(  locId,  layer_nb  )
                        porosity = porosity + strat_porosity( locId, layer_nb, ks ) * prop
                    enddo
                    ! Ensure sufficient sediment exist
                    do ks = 1, totgrn
                        if( strat_sedh(  locId,  layer_nb,  ks  ) > 0.0_8 )then
                            if( ocean_time * ( abs( transportX( i, ks ) ) + abs( transportY( i, ks ) ) ) / ( strat_dx * &
                                ( 1 - porosity ) ) > strat_sedh(  locId,  layer_nb,  ks  ) )then

                                prop = strat_sedh(  locId,  layer_nb,  ks  ) / ( ocean_time * ( abs( transportX( i, ks ) ) + &
                                    abs( transportY( i, ks ) ) ) / ( strat_dx * ( 1 - porosity ) ) )
                                transportX( i, ks ) = prop * transportX( i, ks )
                                transportY( i, ks ) = prop * transportY( i, ks )
                            endif
                        else
                            transportX( i, ks ) = 0.0_8
                            transportY( i, ks ) = 0.0_8
                        endif
                    enddo
                else
                    total_load( i, 1:totgrn ) = 0.0_8
                endif
            endif
        enddo

        ! Update total load vector components border values
        do i = 1 , DL_nodesNb
            if( diff_coord( i, 1 ) < strat_xo .or. diff_coord( i, 2 ) < strat_yo .or. &
                diff_coord( i, 1 ) > gv_coord( G_strat%nodes_nb, 1 ) .or. &
                diff_coord( i, 2 ) > gv_coord( G_strat%nodes_nb, 2 ) )then

                lid = diff_ID( i )
                if( lid == -1 .and. diff_coord( i, 3 ) < 7.e4_8 )then
                    transportX( i, 1:totgrn ) = 0.0_8
                    transportY( i, 1:totgrn ) = 0.0_8
                elseif( lid == 0 )then
                    lid = diff_ref( i )
                    ! South
                    if( lid == -2 .or. lid == -5 )then
                        lid1 = i + 1
                        if( colW == strat_X ) lid1 = i + strat_X + 2
                        transportX( i, 1:totgrn ) = transportX( lid1, 1:totgrn )
                        transportY( i, 1:totgrn ) = transportY( lid1, 1:totgrn )
                    ! North
                    elseif( lid == -1 .or. lid == -6 )then
                        lid1 = i - 1
                        if( colW == strat_X ) lid1 = i - strat_X - 2
                        transportX( i, 1:totgrn ) = transportX( lid1, 1:totgrn )
                        transportY( i, 1:totgrn ) = transportY( lid1, 1:totgrn )
                    ! West
                    elseif( lid == -3 .or. lid == -7 )then
                        lid1 = i + strat_Y + 2
                        if( colW == strat_X ) lid1 = i + 1
                        transportX( i, 1:totgrn ) = transportX( lid1, 1:totgrn )
                        transportY( i, 1:totgrn ) = transportY( lid1, 1:totgrn )
                    ! East
                    elseif( lid == -4 .or. lid == -8 )then
                        lid1 = i - strat_Y - 2
                        if( colW == strat_X ) lid1 = i - 1
                        transportX( i, 1:totgrn ) = transportX( lid1, 1:totgrn )
                        transportY( i, 1:totgrn ) = transportY( lid1, 1:totgrn )
                    endif
                endif

            endif

        enddo

        return

    end subroutine update_total_load_vector_component
    ! ============================================================================
    !> Subroutine continuity_equation_sediment_transport
    !!
    !! The rate of bed level changes is determined from the equation for conservation of mass,
    !! for bed material, using the instantaneous sediment transport rates. The values of the
    !! transport rates are calculated at each point of the computational domain at each
    !! morphological time step.
    !!
    !<
    ! ============================================================================
    subroutine continuity_equation_sediment_transport

        integer :: n, locId, n1, n2, layer_nb, ks, colW

        real( tkind ) :: delta_qx, delta_qy, porosity, prop

        real( tkind ), dimension( L_strat%nodes_nb, totgrn ) :: delta_z

        colW = strat_X
        if( strat_X > strat_Y ) colW = strat_Y

        ! Yearly variation of seabed morphology
        delta_z = 0.0_8

        do n = 1, DL_nodesNb

            locId = diff_ID( n )

            if( locId > 0 )then

                ! Get top layer number
                layer_nb = locv_NbLay( locId )

                ! Get top layer porosity
                porosity = 0.0_8
                prop = 0.0_8
                if( strat_thick(  locId,  layer_nb  ) > 0.0_8 )then
                    do ks = 1, totgrn
                        prop = strat_sedh(  locId,  layer_nb,  ks  ) / strat_thick(  locId,  layer_nb  )
                        porosity = porosity + strat_porosity( locId, layer_nb, ks ) * prop
                    enddo
                endif
                delta_z( locId, : ) = 0.0_8
                do ks = 1, totgrn

                    ! Central difference scheme for sediment transport along X
                    n1 = n + 1
                    n2 = n - 1
                    delta_qx = -abs( transportX( n, ks ) )
                    if( transportX( n1, ks ) < 0.0_8 )then
                        delta_qx = delta_qx + abs(transportX( n1, ks ))
                    endif
                    if( transportX( n2, ks ) > 0.0_8 )then
                        delta_qx = delta_qx + abs(transportX( n2, ks ))
                    endif
                    delta_qx = delta_qx / strat_dx

                    ! Central difference scheme for sediment transport along Y
                    n1 = n + ( colW + 2 )
                    n2 = n - ( colW + 2 )

                    delta_qy = -abs( transportY( n, ks ) )
                    if( transportY( n1, ks ) < 0.0_8 )then
                        delta_qy = delta_qy + abs(transportY( n1, ks ))
                    endif

                    if( transportY( n2, ks ) > 0.0_8 )then
                        delta_qy = delta_qy + abs(transportY( n2, ks ))
                    endif

                    delta_qy = delta_qy / strat_dx

                    ! Get the change of morphology during time step
                    delta_z( locId, ks ) = delta_z( locId, ks ) + ocean_time * ( delta_qx + delta_qy ) / ( 1 - porosity )
                    if( -delta_z( locId, ks ) > strat_sedh(  locId,  layer_nb,  ks  ) )then
                        if( -delta_z( locId, ks ) - strat_sedh(  locId,  layer_nb,  ks  ) > 0.001_8 )then
                            print*,'Something went wrong adjusting the morphological changes.'
                            print*, -delta_z( locId, ks ) - strat_sedh(  locId,  layer_nb,  ks  )
                            stop
                        else
                            delta_z( locId, ks ) = -strat_sedh(  locId,  layer_nb,  ks  )
                        endif
                    endif
                enddo
            endif
        enddo

        ! Found the elevation changes, adjust sedimentary layers accordingly
        do n = 1, L_strat%nodes_nb
            do ks = 1, totgrn
                ! In case of erosion check sediment availabilty and update top layer
                if( delta_z( n, ks ) < -tor )&
                    call ocean_erosion_layer( n, ks, -delta_z( n, ks ) )
            enddo

            do ks = 1, totgrn
                ! In case of deposition update top layer
                if( delta_z( n, ks ) > tor )&
                    call ocean_deposit_layer( n, ks, delta_z( n, ks ) )
            enddo

        enddo

        do n = 1, L_strat%nodes_nb
            ! Check the modified sedimentary layer
            call check_ocean_sedimentary_layer( n )
        enddo

        return

    end subroutine continuity_equation_sediment_transport
    ! ============================================================================
    !> Subroutine ocean_deposit_layer
    !! Subroutine ocean_deposit_layer gets new sediment class composition due to deposition.
    !! \param kn, ks, dcomp
    !<
    ! ============================================================================
    subroutine ocean_deposit_layer( locID, ks, dcomp  )

        integer :: ks, lb, locID, layID

        real( tkind ) :: dcomp

        ! Find local ID and layer number
        lb = locv_NbLay( locID )
        layID = strat_layID(  locID , lb )

        ! Create a new layer
        if( layID < L_strat%layersID )then
            lb = lb + 1
            strat_sedh( locID, lb, 1:max_grn ) = 0.0_8
            locv_NbLay( locID ) = lb
            strat_porosity(  locID,  lb, 1:max_grn ) = 0.0_8
            strat_layID(  locID ,  lb  ) = L_strat%layersID
            strat_thick(  locID ,  lb  ) = dcomp
            strat_zelev(  locID ,  lb  ) = dcomp + strat_zelev(  locID ,  lb - 1  )
            strat_hardness(  locID ,  lb  ) = 1.0_8
            strat_sedh(  locID ,  lb ,  ks  ) = dcomp
            if( gporo%compaction ) strat_porosity(  locID,  lb, ks ) = porosity( ks, 1 )
        ! Add on an existing layer
        else
            strat_thick(  locID ,  lb  ) = dcomp + strat_thick(  locID ,  lb )
            strat_zelev(  locID ,  lb  ) = dcomp + strat_zelev(  locID ,  lb )
            strat_sedh(  locID ,  lb ,  ks  ) = dcomp + strat_sedh(  locID ,  lb ,  ks  )
            if( gporo%compaction ) strat_porosity(  locID,  lb, ks ) = porosity( ks, 1 )
        endif

        return

    end subroutine ocean_deposit_layer
    ! ============================================================================
    !> Subroutine ocean_erosion_layer
    !! Subroutine ocean_erosion_layer gets new sediment class composition due to erosion.
    !! \param kn, ks, delta_z
    !<
    ! ============================================================================
    subroutine ocean_erosion_layer( locID, ks, ecomp  )

        integer :: ks, lb, locID, layID

        real( tkind ) :: ecomp, th_sed

        ! Find local ID and layer number
        lb = locv_NbLay( locID )
        layID = strat_layID(  locID, lb )

        ! Take into account the layer hardness
        ecomp = ecomp / strat_hardness( locID, lb )

        ! If erosion greater than sediment available in the layer
        ! take it all
        if( strat_sedh( locID , lb, ks ) < ecomp )&
            ecomp = strat_sedh( locID, lb, ks )

        ! Adjust sediment layer
        th_sed = strat_sedh( locID, lb, ks )
        strat_sedh( locID, lb, ks ) = -ecomp + strat_sedh( locID, lb, ks )
        if( strat_sedh( locID, lb, ks ) < tor )then
            strat_sedh( locID, lb, ks ) = 0.0_8
            ecomp = th_sed
        endif

        strat_thick( locID, lb ) = -ecomp + strat_thick( locID, lb )
        strat_zelev( locID, lb ) = -ecomp + strat_zelev( locID, lb )

        ! If not sufficient sediment in the top layer just delete this sediment layer
        if( strat_thick( locID, lb ) < tor .and. lb >= 2 )then
            strat_thick( locID, lb ) = 0.0_8
            strat_zelev( locID, lb ) = strat_zelev(  locID ,  lb - 1  )
            strat_sedh( locID, lb, 1:max_grn ) = 0.0_8
            strat_porosity(  locID,  lb, 1:max_grn ) = 0.0_8
            if( locv_NbLay( locID ) > 1 ) locv_NbLay( locID ) = lb - 1
        endif

        return

    end subroutine ocean_erosion_layer
    ! ============================================================================
    !> Subroutine check_ocean_sedimentary_layer
    !! Subroutine check_ocean_sedimentary_layer check consistency in the modified sedimentary layers.
    !! \param kn
    !<
    ! ============================================================================
    subroutine check_ocean_sedimentary_layer( kn  )

        integer :: kn, ks, lb0, lb, locID, layID

        real( tkind ) :: h, layh


        ! Find local ID and layer number
        locID = kn
        lb0 = locv_NbLay( locID )
        layID = strat_layID(  locID,  lb0  )

        if( locv_NbLay( locID ) == 1 ) return

        lb = lb0
        if( strat_layID(  locID,  lb0 - 1  ) > 1 ) lb = lb0 - 1

10  continue

        ! In case a layer still exist
        if( strat_layID(  locID,  lb  ) > 1 )then

            ! Get layer thickness
            layh = strat_thick(  locID,  lb  )

            ! If not enough sediment just delete the layer
            if( layh <= tor )then
                strat_zelev(  locID ,  lb  ) = strat_zelev(  locID ,  lb - 1  )
                strat_thick(  locID ,  lb  ) = 0.0_8
                strat_sedh( locID, lb, 1:max_grn ) = 0.0_8
                strat_porosity(  locID,  lb, 1:max_grn ) = 0.0_8
                if( locv_NbLay( locID ) > 1 ) locv_NbLay( locID ) = lb - 1
                if( lb < lb0 )then
                    lb = lb0
                    goto 10
                endif
            endif

            ! Get total thickness of sediment layer based on sediment classes thicknesses
            h = 0.0_8
            do ks = 1, totgrn
                if( strat_sedh(  locID ,  lb ,  ks  ) < tor ) strat_sedh(  locID ,  lb ,  ks  )  = 0.0_8
                h = h +  strat_sedh(  locID ,  lb ,  ks  )
            enddo

            ! In case there is a mismatch between calculated thickness
            ! and recorded one make it match
            if( h /= layh .and. h > tor )then
                strat_thick(  locID ,  lb  ) = h
                strat_zelev(  locID ,  lb  ) = strat_zelev(  locID ,  lb - 1  ) + h
            elseif( h <= tor )then
                strat_zelev(  locID ,  lb  ) = strat_zelev(  locID ,  lb - 1  )
                strat_thick(  locID ,  lb  ) = 0.0_8
                strat_sedh( locID, lb, 1:max_grn ) = 0.0_8
                strat_porosity( locID, lb, 1:max_grn ) = 0.0_8
                if( locv_NbLay( locID ) > 1 ) locv_NbLay( locID ) = lb - 1
            endif

        endif

        ! Check we are at the top layer otherwise keep checking layers.
        if( lb < lb0 )then
            lb = lb0
            goto 10
        endif

        return

    end subroutine check_ocean_sedimentary_layer
    ! ============================================================================

end module ocean_transport
