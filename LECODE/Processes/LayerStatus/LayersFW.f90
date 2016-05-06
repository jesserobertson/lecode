! ============================================================================
! Name        : LayersFW.f90
! Author      : tristan salles
! Created on: Aug 17, 2012
! Copyright (C) 2012 CSIRO
!============================================================================
!> \file LayersFW.f90
!!
!! File LayersFW determines the new concentration and the deposit layer caracteristics.
!!
!<
! ============================================================================
module mod_layersFW

    use file_data
    use mpidata
    use flux_data
    use time_data
    use mod_sort
    use error_data
    use fwalker_data
    use quicksort
    use strata_ini
    use strata_data
    use param_data
    use forces_data
    use TIN_function
    use mod_diffuse
    use strata_update
    use sediment_data
    use mod_functionsfw
    use mod_compaction

    implicit none

    integer, dimension( : ), allocatable :: fw_gid

contains

    ! ============================================================================
    !> Subroutine apply_erosion_deposition_rules
    !! Subroutine apply_erosion_deposition_rules change the sediment load of the flow
    !! walkers by eroding or depositing to the grid.
    !!
    !! If sediment concentration divided by sediment transportability is greater than
    !! flow's transport capacity, deposition occurs.  If sediment concentration
    !! divided by sediment transportability is less than transport capacity, erosion
    !! occurs provided that critical shear stress for material at water-sediment
    !! interface is exceeded.
    !<
    ! ============================================================================
    subroutine apply_erosion_deposition_rules

        logical :: conc_transfert

        integer :: k, kn, ks, layID, ko, p, gid

        integer, dimension( num_fwed ) :: forder
        real( tkind ), dimension( num_fwed ) :: fvs

        real( tkind ) :: fall, fv, alpha, thick, ethick, dthick, prop, porosity
        real( tkind ) :: factor, hardness, layh, eth, maxdep, dtt, compe
        real( tkind ), dimension( totgrn ) :: laycomp, conceq, conc, erodep
        real( tkind ), dimension( totgrn ) :: ecomp, dcomp, edconc

        real( tkind ) :: szt

        ! Allocate erosion deposition array
        if( allocated( concED ) ) deallocate( concED )
        allocate( concED( num_fwed, totgrn ) )

        ! Arrange flow walkers by increasing velocity
        fvs = 0.0_8
        forder = -1
        do k = 1, num_fwed
            if( ed_bound( k ) /= 1 )&
                fvs( k ) = ed_vel( k )
            forder( k ) = k
        enddo

        if( num_fw > 0 ) call quick_sort( fvs, forder )
        if( flagx )then
            top_sedprev = 0.0_8
            top_sedh = 0.0_8
            flagx = .false.
        endif
        
        eromax = -1.0e6_8

        ! Compute erosion or deposition
        do ko = 1, num_fwed

            k = forder( ko )
            concED( k, 1:totgrn ) = 0.0_8

            if( ed_bound( k ) /= 1 )then

                fv = ed_vel( k )
                conc_transfert = .true.

                do p = 1, ed_AInb( k )

                    gid = ed_AIpts( k, p )

                    if( gid > 0 )then
                        if( gv_SharenID( gid, iam + 1 ) > 0 )then

                            kn = gv_SharenID( gid, iam + 1 )

                            ! If the FW belongs to another processor
                            if( ed_procid( k ) /= iam )then
                                ! If the node sits at the border between 2 processors
                                if( gv_SharenID( gid, ed_procid( k ) + 1 ) > 0 )then
                                    ! Set a no transfert flag of FW concentration changes
                                    conc_transfert = .false.
                                endif
                            endif

                            ! Force deposit of flow walker below sediment load threshold
                            if( ed_bound( k ) /= 0 )then
                                dcomp = 0.0_8
                                eth = 0.0_8
                                do ks = 1, totgrn
                                    dcomp( ks ) = ed_sed( k , ks ) / ( ed_AInb( k ) * &
                                        strat_dx**2.0_8 )
                                    eth = dcomp( ks ) + eth
                                    call add_deposit_layer( kn, ks, dcomp( ks ) )
                                    ed_sed( k , ks ) = 0.0_8
                                enddo
                                goto 10
                            endif

                            ! Update active layer
                            call build_active_nodelayer( transport%erosion_limit, kn )

                            ! Get the thickness of the top layer and the composition information
                            call get_layer_composition( kn, layh, laycomp, porosity, hardness, layID )

                            ! Flow concentration
                            conc = 0.0_8

                            do ks = 1, totgrn
                                conc( ks ) = ed_sed( k , ks ) / ed_vol( k )
                            enddo

                            ! Equilibrium concentration
                            conceq = 0.0_8
                            call get_sediment_equilibrium_concentration( kn, k, conceq )

                            ! alpha: Adaptation length coefficient ( Armanini ) alpha is set to 1.0,
                            ! values range from 0.25 to 5.0 according to Wu (2008)
                            alpha = 1.0_8
                            do ks = 1, totgrn
                                fall = sediment( ks )%vfall
                                erodep( ks ) = ( conc( ks ) - conceq( ks ) ) * fall
                                if( conceq( ks ) == 0.0_8 .and. ed_type( k ) >= tprain &
                                    .and. ed_type( k ) < tpturb .and. fv < transport%fvmin )then
                                    alpha = 2.5_8
                                    dtt = rain_int * secyear
                                    erodep( ks ) = erodep( ks ) * alpha * dtt / ed_h( k )
                                    erodep( ks ) = erodep( ks ) / ( 1 - porosity )
                                elseif( ed_type( k ) >= tpplum1 .and. erodep( ks ) < 0.0_8  )then
                                    erodep( ks ) = 0.0_8
                                else
                                    dtt = flow_int * secyear
                                    erodep( ks ) = erodep( ks ) * alpha * dtt / ed_h( k )
                                    erodep( ks ) = erodep( ks ) / ( 1 - porosity )
                                endif
                            enddo
                            szt = elev_record( kn )

                            ! Get new sedimentary composition due to erosion and sedimentation
                            thick = 0.0_8
                            ethick = 0.0_8
                            dthick = 0.0_8
                            ecomp = 0.0_8
                            dcomp = 0.0_8

                            ! Maximum thickness allowed to be eroded
                            if( eromax( kn ) == -1.0e6_8 )then
                                call get_maximum_erosion_thickness( kn, eromax( kn ) )
                                eromax( kn ) = -eromax( kn )
                            endif

                            do ks = 1, totgrn
                                ecomp( ks ) = 0.0_8
                                dcomp( ks ) = 0.0_8
                                if( erodep( ks ) < 0.0_8 )then
                                    ecomp( ks ) = erodep( ks ) !/  * ed_vol( k ) / &
                                     !   ( ed_AInb( k ) * strat_dx**2.0_8 * hardness )
                                    if( ecomp( ks ) > -tor ) ecomp( ks ) = 0.0_8
                                    thick = thick + ecomp( ks )
                                    ethick = ethick - ecomp( ks )
                                elseif( erodep( ks ) > 0.0_8 )then
                                    dcomp( ks ) = erodep( ks ) !* ed_vol( k ) / &
                                     !   ( ed_AInb( k ) * strat_dx**2.0_8  )
                                     maxdep = ed_sed( k , ks ) / ( ed_AInb( k ) * &
                                        strat_dx**2.0_8 )
                                    if( maxdep < dcomp( ks ) ) dcomp( ks ) = maxdep
                                    if( dcomp( ks ) < tor ) dcomp( ks ) = 0.0_8
                                    thick = thick + dcomp( ks )
                                    dthick = dthick + dcomp( ks )
                                endif
                            enddo
                            
                            ! Limited erosion based on neighbors elevation
                            if( -ethick < eromax( kn ) )then
                                eth = 0.0_8
                                do ks = 1, totgrn
                                    if( erodep( ks ) < 0.0_8 )then
                                        prop = -ecomp( ks ) / ethick
                                        ecomp( ks ) = eromax( kn ) * prop
                                        eth = eth - ecomp( ks )
                                    endif
                                enddo
                                eromax( kn ) = 0.0_8
                            else
                                eromax( kn ) = eromax( kn ) + ethick
                            endif

                            ! Limited deposition based on source elevation
                            factor = transport%fac_limit
                            if( ed_type( k ) >= tpplum1 )then
                                if( ed_type( k ) == tpplum1 )then
                                    factor = 10.0_8 *  transport%fac_limit
                                else
                                    factor = 50.0_8 *  transport%fac_limit
                                endif
                            endif
                            if( fv < transport%fvmin ) factor = 10.0_8 * transport%fac_limit

!                            if( dthick > factor * ed_h( k ) )then
!                                do ks = 1, totgrn
!                                    if( erodep( ks ) > 0.0_8 )then
!                                        prop = dcomp( ks ) / dthick
!                                        dcomp( ks ) = factor * ed_h( k ) * prop
!                                        if( dcomp( ks ) < tor ) dcomp( ks ) = 0.0_8
!                                    endif
!                                enddo
!                            endif

                            ! Update the sedimentary layer due to erosion
                            edconc = 0.0_8
                            do ks = 1, totgrn
                                if( -ecomp( ks ) > tor )then

                                    compe = 0.0_8
                                    call get_eroded_layer_composition( kn, ks, -ecomp( ks ), compe )
                                    ecomp( ks ) = -compe

                                    ! Volume accumulated within the flow
                                    edconc( ks ) = -ecomp( ks ) * ( 1 - porosity ) * strat_dx**2.0_8
                                    if( conc_transfert ) concED( k, ks ) = concED( k, ks ) + edconc( ks )
                                endif
                            enddo
                            call update_top_layer_parameters( kn )

                            ! Update the sedimentary layer due to deposition
                            do ks = 1, totgrn
                                if( dcomp( ks ) > tor )then
                                    call add_deposit_layer( kn, ks, dcomp( ks ) )

                                    ! Volume deposited from the flow
                                    edconc( ks ) = - dcomp( ks ) * ( 1 - porosity ) * strat_dx**2.0_8
                                    if( conc_transfert ) concED( k, ks ) = concED( k, ks ) + edconc( ks )
                                endif
                            enddo

                            ! Check the modified sedimentary layer
                            call check_modified_sedimentary_layer( kn )

                            ! Update top porosity values
                            if( gporo%compaction ) call update_top_layer_porosity( kn )

                        endif

                    endif

                enddo

            endif

10      continue

        enddo

        ! Merge FW sedimentary composition
        call update_flow_walkers_concentration_and_topography

        return

    end subroutine apply_erosion_deposition_rules
    ! ============================================================================
    !> Subroutine combine_rainwalker
    !! Subroutine combine_rainwalker combines rain flow walkers when they are within the same face.
    !!
    !<
    ! ============================================================================
    subroutine combine_rainwalkers

        integer :: k, p, ks

        real( tkind ) :: v1, v2, v3, A1, A2, A3, v

        ! For each flow walkers
        do k = 1, num_fw - 1

            ! If this is a rain walker which could be combined
            if( fa_inbound( k ) == 0 .and. fa_type( k ) >= tprain .and. &
                fa_type( k ) < tprain + combine .and. fw_gid( k ) > 0 )then

                ! Check other flow walkers
                do p = k + 1, num_fw

                    ! If rain walker can be merged
                    if( fa_inbound( p ) == 0 .and. fa_type( p ) >= tprain .and. &
                        fa_type( p ) < tprain + combine .and. fw_gid( k ) == fw_gid( p )  )then

                        fa_inbound( p ) = 1

                        ! Wet areas
                        A1 = fa_width( k ) * fa_h( k )
                        A2 = fa_width( p ) * fa_h( p )

                        ! Calculate the Euclidean norm of the FW's current velocity.
                        v1 = sqrt( fa_xvel( k )**2 + fa_yvel( k )**2 )
                        v2 = sqrt( fa_xvel( p )**2 + fa_yvel( p )**2 )

                        ! Get new width and height
                        fa_width( k ) = max( fa_width( k ), fa_width( p ) )
                        fa_h( k ) = max( fa_h( k ), fa_h( p ) )

                        ! Get new wet area
                        A3 = fa_width( k ) * fa_h( k )
                        v3 = ( A1 * v1 + A2 * v2 ) / A3

                        fa_discharge( k ) =  fa_discharge( k ) + fa_discharge( p )
                        fa_xvel( k ) =  0.5_8 * ( fa_xvel( k ) + fa_xvel( p ) )
                        fa_yvel( k ) =  0.5_8 * ( fa_yvel( k ) + fa_yvel( p ) )

                        ! Calculate the Euclidean norm of the FW's current velocity.
                        v = sqrt( fa_xvel( k )**2 + fa_yvel( k )**2 )

                        if( v /= 0.0_8 )then
                            fa_xvel( k ) = fa_xvel( k ) * v3 / v
                            fa_yvel( k ) = fa_yvel( k ) * v3 / v
                        else
                            fa_xvel( k ) = max( fa_xvel( k ), fa_xvel( p ) )
                            fa_yvel( k ) = max( fa_yvel( k ), fa_yvel( p ) )
                        endif
                        ks = p
                        if(  fa_volume( k ) >  fa_volume( p ) ) ks = k
                        fa_xpos( k ) = fa_xpos( ks )
                        fa_ypos( k ) = fa_ypos( ks )
                        fa_zpos( k ) = fa_zpos( ks )
                        fa_xslp( k ) =  fa_xslp( ks )
                        fa_yslp( k ) =  fa_yslp( ks )

                        fa_density( k ) = 0.0_8
                        do ks = 1, totgrn
                            fa_sedcharge( k , ks ) = fa_sedcharge( k , ks ) + fa_sedcharge( p , ks )
                            fa_density( k ) = fa_density( k ) + fa_sedcharge( k , ks ) * sediment( ks )%density
                        enddo
                        fa_volume( k ) =  fa_volume( k ) + fa_volume( p )
                        fa_discharge( k ) =  fa_discharge( k ) + fa_discharge( p )
                        fa_density( k ) = fa_density( k ) / fa_volume( k ) + fluid_density

                        ! Delete merged rain walker
                        fa_h( p ) = 0.0_8
                        fa_xvel( p ) = 0.0_8
                        fa_yvel( p ) = 0.0_8
                        fa_density( p ) = 0.0_8
                        fa_sedcharge( p , 1:totgrn ) = 0.0_8
                        if( fa_volume( k ) >= rainstreamQ ) fa_type( k ) =  tprain + combine

                    endif
                enddo
            endif
        enddo

        return

    end subroutine combine_rainwalkers
    ! ============================================================================
    !> Subroutine update_stratigraphic_layers()
    !! Allocate the stratigraphy to each processors depending on the partitioning.
    !<
    ! ============================================================================
    subroutine update_stratigraphic_layers

        ! Parameters Declaration
        integer :: nb_lay, n, p, GID, m, k

        real( tkind ) :: zh

        ! Check stratigraphic layers
        do n = 1, L_strat%nodes_nb
            do p = locv_NbLay( n ), 2, -1
                if( strat_thick(  n ,  p  ) == 0.0_8 )then
                    do m = p, locv_NbLay( n ) - 1
                        strat_layID( n, m ) = strat_layID( n, m+1 )
                        strat_hardness( n, m ) = strat_hardness( n, m+1 )
                        strat_zelev( n, m ) = strat_zelev( n, m+1 )
                        strat_thick( n, m ) = strat_thick( n, m+1 )
                        strat_fac( n, m ) = strat_fac( n, m+1 )
                        strat_porosity( n, m, 1:max_grn ) = strat_porosity( n, m+1, 1:max_grn )
                        strat_sedh( n, m, 1:max_grn ) = strat_sedh( n, m+1, 1:max_grn )
                    enddo
                    strat_layID( n , locv_NbLay( n ) ) = -1
                    locv_NbLay( n ) = locv_NbLay( n ) - 1
                endif
            enddo

            do p = 2, locv_NbLay( n )
                zh = 0.0_8
                do k = 1, totgrn
                    zh = zh + strat_sedh(  n ,  p ,  k  )
                enddo
                strat_thick(  n ,  p  ) = zh
                strat_zelev(  n ,  p  ) = strat_zelev(  n ,  p - 1  ) + zh
            enddo
        enddo

        ! Broadcast top elevation to global grid
        proc_elev = -1.0e6_8
        do n = 1, L_strat%nodes_nb
            GID = locv_gid( n )
            nb_lay = locv_NbLay( n )
            proc_elev( GID ) = strat_zelev(  n ,  nb_lay  )
        enddo
        call mpi_allreduce( proc_elev, elev_record , G_strat%nodes_nb, &
            dbl_type, max_type, SPModel_comm_world, ierr )

        return

    end subroutine update_stratigraphic_layers
    ! ============================================================================
    
end module mod_layersFW
