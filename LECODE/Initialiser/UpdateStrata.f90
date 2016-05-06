! ============================================================================
! Name        : UpdateStrata.f90
! Author      : tristan salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file UpdateStrata.f90
!!
!! UpdateStrata updates stratigraphic mesh through time..
!!
!<
! ============================================================================
module strata_update

    use file_data
    use mpidata
    use time_data
    use strata_ini
    use strata_data
    use strata_parser
    use sediment_data

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine build_active_nodelayer()
    !!
    !! The thickness of the top sedimentary layer has a minimum thickness equivalent to za.
    !! If the top layer is thicker than za, no action is required. If the top layer is less than za thick,
    !! then the top layer thickness is increased by entraining sediment mass from deeper layers until
    !! the top layer thickness equals za. If sediment from deeper than the second layer is mixed into
    !! the top layer, the bottom layer is split to enforce a constant number of layers and conservation
    !! of sediment mass.
    !! Allocate the stratigraphy to each processors depending on the partitioning.
    !<
    ! ============================================================================
    subroutine build_active_nodelayer( hlimit, n )

        ! Parameters Declaration
        integer :: n, p, ks, nb_lay, lb, layid

        real( tkind ) :: th, cum_h, hardness, prop, h_grab, hlimit, th2
        real( tkind ), dimension( totgrn ) ::  sed_comp, bot_lay_sed

        ! Get top layer number
        nb_lay = locv_NbLay( n )

        ! Check top sediment layer thickness
        if( strat_thick(  n ,  nb_lay  ) < hlimit .and. nb_lay > 2 )then

            ! Loop through the layers from top to bottom to find
            ! the layers that needs to be activated
            th = 0.0_8
            cum_h = 0.0_8
            hardness = 0.0_8
            sed_comp = 0.0_8
            h_grab = 0.0_8
            layid = 1

            recursive_thickness: do p = locv_NbLay( n ), 2, -1

                th = th + strat_thick( n, p )

                ! In case all the sedimentary layer is taken to the active layer
                if( th <= hlimit )then
                    do ks = 1, totgrn
                        sed_comp( ks ) = sed_comp( ks ) + strat_sedh( n, p, ks )
                        strat_sedh(  n, p, ks ) = 0.0_8
                        strat_porosity( n, p, ks ) = 0.0_8
                    enddo
                    hardness = max( hardness, strat_hardness( n, p ) )
                    cum_h = cum_h + strat_thick( n, p )
                    layid = max( strat_layID(  n,  p  ), layid )
                    strat_thick(  n, p ) = 0
                    strat_hardness( n, p ) = 0
                    nb_lay = p
                endif

                if( th > hlimit )then

                    h_grab = hlimit - cum_h
                    bot_lay_sed = 0.0_8

                    ! Get proportion of sediment present in the layer
                    ! And adjust layer thickness
                    th2 = 0.0_8
                    do ks = 1, totgrn
                        prop = strat_sedh( n, p, ks ) / strat_thick( n, p )
                        sed_comp( ks ) = sed_comp( ks ) + h_grab * prop
                        bot_lay_sed( ks ) = strat_sedh( n, p, ks ) - h_grab * prop
                        strat_sedh(  n, p, ks ) = strat_sedh( n, p, ks ) - h_grab * prop
                        th2 = th2 + strat_sedh(  n, p, ks )
                    enddo
                    strat_thick(  n, p ) = th2
                    hardness = min( hardness, strat_hardness( n, p ) )
                    strat_zelev(  n, p ) = strat_zelev(  n, p-1 ) + th2
                    exit recursive_thickness
                endif

            enddo recursive_thickness

            ! Now move the combined layers to form the active one
            locv_NbLay( n ) =  nb_lay
            lb = locv_NbLay( n )

            if( lb >  L_strat%layersID )then
                print*,"Something went wrong when updating the active layer index."
                print*,iam,n,lb,L_strat%layersID
                stop
            endif

            strat_layID(  n,  lb  ) = layid
            strat_sedh( n, lb, 1:max_grn ) = 0.0_8
            strat_porosity( n, lb, 1:max_grn ) = 0.0_8
            th = 0.0_8

            do ks = 1, totgrn
                th = th + sed_comp( ks )
                strat_sedh(  n, lb, ks ) = sed_comp( ks )
                if( gporo%compaction ) strat_porosity(  n,  lb, ks  ) = porosity( ks, 1 )
            enddo

            strat_hardness( n, lb ) = hardness
            strat_thick( n, lb ) = th
            strat_zelev( n, lb ) = strat_zelev(  n,  lb - 1  ) + th

            if( abs( th - hlimit ) > 0.0001_8 .and. lb > 2 )then
                print*,"Something went wrong when updating the active layer thickness."
                print*,iam,n,th,hlimit,lb,th-hlimit
                stop
            endif

        endif

        return

    end subroutine build_active_nodelayer
    ! ============================================================================

end module strata_update
