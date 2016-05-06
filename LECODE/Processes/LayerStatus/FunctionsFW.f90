! ============================================================================
! Name        : FunctionsFW.f90
! Author      : tristan salles
! Created on: Aug 17, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file FunctionsFW.f90
!!
!! File FunctionsFW determines the new concentration and the deposit layer caracteristics.
!!
!<
! ============================================================================
module mod_functionsfw

    use file_data
    use mpidata
    use flux_data
    use time_data
    use mod_sort
    use error_data
    use fwalker_data
    use strata_data
    use param_data
    use forces_data
    use TIN_function
    use sediment_data

    implicit none

    type flowConc
        integer :: id = -1
        real(tkind), dimension(max_grn) :: conc
    end type flowConc

    integer :: mpi_concfw

    real( tkind ), dimension( :,: ), allocatable :: concED

    public

contains

    ! ============================================================================
    !> Subroutine build_mpi_flowconc_type
    !! Subroutine build_mpi_flowconc_type defines a mpi data type for passing FW
    !! concentration between processors for erosion deposition.
    !<
    ! ============================================================================
    subroutine build_mpi_flowconc_type

        integer, parameter :: count = 2
        integer, dimension( count ) :: array_of_types, array_of_block_lengths
        integer( kind = mpi_address_kind ), dimension( count ) :: array_of_displacements, array_of_address

        type( flowConc ) :: sample

        integer :: k

        array_of_types=(/int_type,dbl_type/)

        array_of_block_lengths = (/1,max_grn/)

        ! Set up the displacements and base types for the FW variables.
        call mpi_get_address(sample%id, array_of_address( 1 ), ierr)
        call mpi_get_address(sample%conc, array_of_address( 2 ), ierr)

        do k = 1, 2
            array_of_displacements( k ) = array_of_address( k ) - array_of_address( 1 )
        enddo

        ! Create and commit the type for MPI.
        call mpi_type_create_struct(count, &
            array_of_block_lengths, &
            array_of_displacements, &
            array_of_types, &
            mpi_concfw, ierr)
        call mpi_type_commit(mpi_concfw, ierr)

        return

    end subroutine build_mpi_flowconc_type
    ! ============================================================================
    !> Subroutine get_sediment_equilibrium_concentration
    !! Subroutine get_sediment_equilibrium_concentration computes equilibrium sediment
    !! concentration of each size class of sediment load.
    !<
    ! ============================================================================
    subroutine get_sediment_equilibrium_concentration( nid, fw, conceq )

        integer :: fw, nid, lb, locID, ks, kk, p, sed_nb

        real( tkind ) :: ratio, sumprop, d50, fv, mid

        real( tkind ), dimension( totgrn ) :: sorted_dia, perc_sed, cum_perc
        real( tkind ), dimension( totgrn ) :: sflux, bflux, conceq, prop, ph, pe

        ! Find the top layer
        locID = nid
        lb = locv_NbLay( locID )
        if( strat_thick(  locID ,  lb  ) == 0.0_8 .and. strat_layID(  locID ,  lb  ) > 1 )then
            lb = lb - 1
            locv_NbLay( locID ) = lb
        endif
        conceq = 0.0_8

        ! If the layer exist
        if( strat_thick(  locID ,  lb  ) > 0.0_8 )then

            ! Get the fraction of each grain present in the top layer
            sumprop = 0.0_8
            do ks = 1, totgrn
                prop( ks ) = strat_sedh(  locID ,  lb ,  ks  )
                prop( ks ) = prop( ks ) /  strat_thick(  locID ,  lb  )
                sumprop = sumprop + prop( ks )
            enddo

            ! Ensure the sediment proportions are suming to 1.0
            ! minor numerical precision
            if( abs( sumprop - 1.0_8 ) > tor )then
                ratio = 1.0_8 / sumprop
                do ks = 1, totgrn
                    prop( ks ) = prop( ks ) * ratio
                enddo
            endif

            ! Get the particles size percentatge and cumulative distribution
            p = 0
            perc_sed = 0.0_8
            cum_perc = 0.0_8
            do ks = totgrn, 1 , -1
                ! We order the sediment by increase diameter
                p = p + 1
                sorted_dia( p ) = sediment( ks )%diameter
                perc_sed( p ) = prop( ks ) * 100.0_8
                if( perc_sed( p ) > 0.0_8 )  sed_nb = sed_nb + 1
                if( p == 1 )then
                    cum_perc( p ) = perc_sed( p )
                else
                    cum_perc( p ) =  cum_perc( p - 1 ) + perc_sed( p )
                endif
            enddo

            ! Find d10, d50, d90 grain size for this particular distribution
            if( sed_nb > 1 )then
                mid = ( sed_nb + 1. ) * 0.5_8
                do ks = 1, totgrn - 1
                    if( mid >= ks .and. mid < ks + 1 )then
                        d50 =  sorted_dia( ks - 1 ) + ( 50.0_8 - cum_perc( ks ) ) * ( sorted_dia( ks + 1 ) - &
                            sorted_dia( ks ) ) / cum_perc( ks + 1 )
                    endif
                enddo
            else
                find_sed: do ks = 1, totgrn
                    if( perc_sed( ks ) > 0.0_8 )then
                        d50 = sorted_dia( ks )
                        exit find_sed
                    endif
                enddo find_sed
            endif
            if( totgrn == 1 ) d50 = sediment( 1 )%diameter

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
                    print*,'Something went wrong computing hidden probabilities.'
                    stop
                endif
                ! Deduce exposed probabilities
                pe( kk ) =  1.0_8 - ph( kk )
                if( pe( kk ) < 0.0_8 )then
                    ph( kk ) = 1.0_8
                    pe( kk ) = 0.0_8
                endif
            enddo

            ! Suspension mode
            sflux = 0.0_8
            if( ed_type( fw ) < tprain .or. ed_type( fw ) >= tprain + combine ) &
                call suspension_sediment_transport( fw, pe, ph, prop, sflux )

            ! Bedload mode
            bflux = 0.0_8
            if( ed_type( fw ) < tpplum1 ) &
                call bedload_sediment_transport( fw, pe, ph, prop, d50, bflux )

            ! Equilibrium sediment concentration
            fv = ed_vel( fw )
            if( fv * ed_h( fw ) > 0.0_8 )then
                do ks = 1, totgrn
                    conceq( ks ) = ( bflux( ks ) + sflux( ks ) ) / ( fv * ed_h( fw ) )
                enddo
            else
                conceq = 0.0_8
            endif
        endif

        return

    end subroutine get_sediment_equilibrium_concentration
    ! ============================================================================
    !> Subroutine bedload_sediment_transport
    !! Subroutine bedload_sediment_transport determines bedload sediment transport.
    !! Wu et al. (2000) proposed a formula to calculate the fractional bed load transport capacity.
    !<
    ! ============================================================================
    subroutine bedload_sediment_transport( fw, pe, ph, prop, d50, bflux )

        integer :: fw, ks

        real( tkind ) :: d50, ntilde
        real( tkind ) :: theta, m, fv, Sf, tau_b, szt, xnn

        real( tkind ), dimension( totgrn ) :: bflux, prop, ph, pe, tau_c, psi

        theta = 0.03_8
        m = -0.6_8

        ! Set Manning's coefficient
        szt = ed_zpos( fw )
        xnn = manning_fct( gsea%actual_sea - szt, ed_dens( fw ) )

        ntilde = d50**(1.0_8/6.0_8) / 20.0_8

        ! Bottom shear stress
        fv = ed_vel( fw )
        Sf = xnn**2 * fv**2 / ( ed_h( fw )**( 4.0_8 / 3.0_8 ) )
        tau_b = ed_dens( fw ) * gravity * ed_h( fw ) * Sf

        ! Criterion for incipient motion and bedload rate
        psi( : ) = 0.0_8
        tau_c( : ) = 0.0_8
        do ks = 1, totgrn
            tau_c( ks ) = theta * sediment( ks )%diameter * ( pe( ks ) / ph( ks ) )**m
            tau_c( ks ) = tau_c( ks ) * ( sediment( ks )%density - ed_dens( fw ) ) * gravity
            if( tau_c( ks ) <= 0.0_8 )then
                bflux( ks ) = 0.0_8
            else
                if( ( ntilde / xnn )**1.5_8 * tau_b / tau_c( ks ) - 1.0_8 < 0.0_8 )then
                    psi( ks ) = 0.0_8
                    bflux( ks ) = 0.0_8
                else
                    psi( ks ) = 0.0053_8 * (  ( ntilde / xnn )**1.5_8 * tau_b / tau_c( ks ) - 1.0_8 )**2.2_8
                    bflux( ks ) = psi( ks ) * prop( ks ) * sqrt( sediment( ks )%diameter**3 * gravity * &
                        ( sediment( ks )%density  / ed_dens( fw ) - 1.0_8 ) )
                endif
            endif
        enddo

        return

    end subroutine bedload_sediment_transport
    ! ============================================================================
    !> Subroutine suspension_sediment_transport
    !! Subroutine suspension_sediment_transport determines suspension sediment transport capacity.
    !! Wu et al (2000) proposed a formula to calculate the fractional suspended load transport capacity.
    !<
    ! ============================================================================
    subroutine suspension_sediment_transport( fw, pe, ph, prop, sflux )

        integer :: fw, ks

        real( tkind ) :: theta, m, fv, Sf, tau_b, szt, xnn, fall

        real( tkind ), dimension( totgrn ) :: sflux, prop, ph, pe, tau_c, psi

        theta = 0.03_8
        m = -0.6_8

        ! Set Manning's coefficient
        szt = ed_zpos( fw )
        xnn = manning_fct( gsea%actual_sea - szt, ed_dens( fw ) )

        ! Bottom shear stress
        fv = ed_vel( fw )
        Sf = xnn**2 * fv**2 / ( ed_h( fw )**( 4.0_8 / 3.0_8 ) )
        tau_b = ed_dens( fw ) * gravity * ed_h( fw ) * Sf

        ! Criterion for incipient motion and suspended load rate
        psi( : ) = 0.0_8
        tau_c( : ) = 0.0_8
        do ks = 1, totgrn
            if( sediment( ks )%diameter > 0.001_8 * 0.0625_8 )then
                tau_c( ks ) = theta * sediment( ks )%diameter
                tau_c( ks ) = tau_c( ks ) * ( pe( ks ) / ph( ks ) )**m
                tau_c( ks ) = tau_c( ks ) * ( sediment( ks )%density - ed_dens( fw ) ) * gravity
                if(  tau_c( ks ) > 0.0_8 )then
                    if(  tau_b / tau_c( ks ) - 1 < 0.0_8 )then
                        psi( ks ) = 0.0_8
                        sflux( ks ) = 0.0_8
                    else
                        fall = sediment( ks )%vfall
                        if( fv < transport%fvmin .or. ed_type( fw ) >= tpplum1 ) fall = 1.0e3_8 * fall
                        psi( ks ) = 2.62e-5_8 * ( ( tau_b / tau_c( ks ) - 1 ) * fv / fall )**1.74_8
                        sflux( ks ) = psi( ks ) * prop( ks ) * sqrt( sediment( ks )%diameter**3 * gravity * &
                            ( sediment( ks )%density  / ed_dens( fw ) - 1.0_8 ) )
                    endif
                else
                    sflux( ks ) = 0.0_8
                endif
            else
                sflux( ks ) = 0.0_8
            endif
        enddo

        return

    end subroutine suspension_sediment_transport
    ! ============================================================================
    !> Subroutine get_layer_composition
    !! Subroutine get_layer_composition gets top layer thickness and composition information.
    !! \param kn, layh, laycomp, porosity, hardness, layID
    !<
    ! ============================================================================
    subroutine get_layer_composition( kn, layh, laycomp, porosity, hardness, layID  )

        integer :: kn, ks, lb, locID, layID

        real( tkind ) :: layh, h, porosity, hardness
        real( tkind ), dimension( totgrn ) :: laycomp

        ! Find local ID and layer number
        locID = kn
        lb = locv_NbLay( locID )
        layID = strat_layID(  locID ,  lb  )

        ! Initialise information
        layh = 0.0_8
        laycomp = 0.0_8
        porosity = 0.0_8
        hardness = 1.0_8

        if( locv_NbLay( locID ) == 1 ) return

        ! If sedimentary layer thickness exists
        if( strat_thick(  locID ,  lb  ) > 0.0_8 )then
            layh = strat_thick(  locID ,  lb  )
            h = 0.0_8
            do ks = 1, totgrn
                laycomp( ks ) = strat_sedh(  locID ,  lb ,  ks  )
                h = h + laycomp( ks )
                porosity = porosity + strat_porosity(  locID ,  lb, ks  ) * laycomp( ks )
            enddo
            if( abs( h - layh ) > tor ) strat_thick(  locID ,  lb  ) = h
            if( abs( h - layh ) > tor ) strat_zelev(  locID ,  lb  ) = &
                strat_zelev(  locID ,  lb - 1  ) + h
            if( abs( h - layh ) > tor ) layh = h
            porosity = porosity / strat_thick(  locID ,  lb  )
            hardness = strat_hardness(  locID ,  lb  )

        ! Otherwise takes the underlying layer and update
        ! sedimentary layer stack
        elseif( strat_layID(  locID ,  lb  ) > 1. .and. &
            strat_thick(  locID ,  lb - 1  ) > 0.0_8 )then
            lb = lb - 1
            locv_NbLay( locID ) = lb
            layh = strat_thick(  locID ,  lb  )
            h = 0.0_8
            do ks = 1, totgrn
                laycomp( ks ) = strat_sedh(  locID ,  lb ,  ks  )
                h = h + laycomp( ks )
                porosity = porosity + strat_porosity(  locID ,  lb, ks  ) * laycomp( ks )
            enddo
            if( abs( h - layh ) > tor ) strat_thick(  locID ,  lb  ) = h
            if( abs( h - layh ) > tor ) strat_zelev(  locID ,  lb  ) = &
                strat_zelev(  locID ,  lb - 1  ) + h
            if( abs( h - layh ) > tor ) layh = h
            porosity = porosity / strat_thick(  locID ,  lb  )
            hardness = strat_hardness(  locID ,  lb  )

        endif

        ! Return layer ID
        layID = strat_layID(  locID ,  lb  )

        return

    end subroutine get_layer_composition
    ! ============================================================================
    !> Subroutine update_top_layer_parameters
    !! Subroutine update_top_layer_parameters checks if the top layer still remain.
    !! \param kn
    !<
    ! ============================================================================
    subroutine update_top_layer_parameters( kn  )

        integer :: kn, ks, lb, locID

        real( tkind ) :: layh, h
        real( tkind ), dimension( totgrn ) :: laycomp

        ! Find local ID and layer number
        locID = kn
        lb = locv_NbLay( locID )

        ! Initialise information
        layh = 0.0_8
        h = 0.0_8
        laycomp = 0.0_8

        if( locv_NbLay( locID ) == 1 ) return

        ! Sedimentary layer thickness
        layh = strat_thick(  locID ,  lb  )

        ! If not sufficient sediment in the layer just delete the layer
        if( layh <= tor )then
            strat_zelev(  locID ,  lb  ) = strat_zelev(  locID ,  lb - 1  )
            strat_thick(  locID ,  lb  ) = 0.0_8
            strat_sedh( locID, lb, 1:max_grn ) = 0.0_8
            strat_porosity(  locID ,  lb, 1:max_grn  ) = 0.0_8
            if( locv_NbLay( locID ) > 1 ) locv_NbLay( locID ) = lb - 1
        endif

        ! Get the sedimentary layer composition
        do ks = 1, totgrn
            laycomp( ks ) = strat_sedh(  locID ,  lb ,  ks  )
            if( laycomp( ks ) < tor ) laycomp( ks ) = 0.0_8
            strat_porosity(  locID ,  lb, ks  ) = 0.0_8
            h = h + laycomp( ks )
        enddo

        ! In case there is a mismatch between calculated thickness
        ! and recorded one make it match
        if( layh /= h )then
            if( h <= tor )then
                strat_zelev(  locID ,  lb  ) = strat_zelev(  locID ,  lb - 1  )
                strat_thick(  locID ,  lb  ) = 0.0_8
                strat_sedh( locID, lb, 1:max_grn ) = 0.0_8
                strat_porosity(  locID ,  lb, 1:max_grn  ) = 0.0_8
                if( locv_NbLay( locID ) > 1 ) locv_NbLay( locID ) = lb - 1
            else
                strat_thick(  locID ,  lb  ) = h
                strat_zelev(  locID ,  lb  ) = h + strat_zelev(  locID ,  lb - 1  )
                do ks = 1, totgrn
                    strat_sedh(  locID ,  lb ,  ks  ) = laycomp( ks )
                    if( laycomp( ks ) > 0.0_8 .and. gporo%compaction .and. &
                        strat_porosity(  locID ,  lb, ks  ) == 0.0_8 ) &
                        strat_porosity(  locID ,  lb, ks  ) = porosity( ks, 1 )
                enddo
            endif
        endif

        return

    end subroutine update_top_layer_parameters
    ! ============================================================================
    !> Subroutine get_eroded_layer_composition
    !! Subroutine get_eroded_layer_composition gets new sediment class composition due to erosion.
    !! \param kn, ks, ecomp
    !<
    ! ============================================================================
    subroutine get_eroded_layer_composition( kn, ks, comp, newC  )

        integer :: kn, ks, lb, locID

        real( tkind ) :: comp
        real( tkind ), intent( out ) :: newC

        ! Find local ID and layer number
        locID = kn
        lb = locv_NbLay( locID )
        newC = 0

        ! Take into account the layer hardness
        comp = comp / strat_hardness(  locID ,  lb  )

        ! If erosion greater than sediment available in the layer
        ! take it all
        if( strat_sedh(  locID ,  lb ,  ks  ) < comp )&
            comp = strat_sedh(  locID ,  lb ,  ks  )

        ! Adjust sediment layer
        strat_sedh(  locID ,  lb ,  ks  ) = -comp + strat_sedh(  locID ,  lb ,  ks  )
        strat_thick(  locID ,  lb  ) = -comp + strat_thick(  locID ,  lb  )
        strat_zelev(  locID ,  lb  ) = -comp + strat_zelev(  locID ,  lb  )
        if( strat_sedh(  locID ,  lb ,  ks  ) > 0.0_8 .and. gporo%compaction .and. &
            strat_porosity( locID ,  lb, ks  ) == 0.0_8 ) &
            strat_porosity(  locID ,  lb, ks  ) = porosity( ks, 1 )

        ! If not sufficient sediment in the layer just delete the sediment layer
        if( strat_thick(  locID ,  lb  ) < tor )then
            strat_thick(  locID ,  lb  ) = 0.0_8
            strat_zelev(  locID ,  lb  ) = strat_zelev(  locID ,  lb - 1  )
            strat_sedh( locID, lb, 1:max_grn ) = 0.0_8
            strat_porosity(  locID ,  lb, 1:max_grn  ) = 0.0_8
        endif
        newC = comp

        return

    end subroutine get_eroded_layer_composition
    ! ============================================================================
    !> Subroutine add_deposit_layer
    !! Subroutine add_deposit_layer gets new sediment class composition due to deposition.
    !! \param kn, ks, dcomp
    !<
    ! ============================================================================
    subroutine add_deposit_layer( kn, ks, dcomp  )

        integer :: kn, ks, lb, locID, layID

        real( tkind ) :: dcomp

        ! Find local ID and layer number
        locID = kn
        lb = locv_NbLay( locID )
        layID = strat_layID(  locID ,  lb  )

        ! Create a new layer
        if( layID < L_strat%layersID )then
            lb = lb + 1
            strat_sedh( locID, lb, 1:max_grn ) = 0.0_8
            strat_porosity(  locID ,  lb, 1:max_grn  ) = 0.0_8
            locv_NbLay( locID ) = lb
            strat_layID(  locID ,  lb  ) = L_strat%layersID
            strat_thick(  locID ,  lb  ) = dcomp
            strat_zelev(  locID ,  lb  ) = dcomp + strat_zelev(  locID ,  lb - 1  )
            strat_hardness(  locID ,  lb  ) = 1.0_8
            strat_sedh(  locID ,  lb ,  ks  ) = dcomp
            if( gporo%compaction ) strat_porosity(  locID,  lb, ks  ) = porosity( ks, 1 )
        ! Add on an existing layer
        else
            strat_thick(  locID ,  lb  ) = dcomp + strat_thick(  locID ,  lb )
            strat_zelev(  locID ,  lb  ) = dcomp + strat_zelev(  locID ,  lb )
            strat_sedh(  locID ,  lb ,  ks  ) = dcomp + strat_sedh(  locID ,  lb ,  ks  )
            if( strat_sedh(  locID ,  lb ,  ks  ) > 0.0_8 .and. gporo%compaction .and. &
                strat_porosity( locID ,  lb, ks  ) == 0.0_8 ) &
                strat_porosity(  locID,  lb, ks  ) = porosity( ks, 1 )
        endif

        return

    end subroutine add_deposit_layer
    ! ============================================================================
    !> Subroutine check_modified_sedimentary_layer
    !! Subroutine check_modified_sedimentary_layer check consistency in the modified sedimentary layers.
    !! \param kn
    !<
    ! ============================================================================
    subroutine check_modified_sedimentary_layer( kn  )

        integer :: kn, ks, lb0, lb, locID, layID

        real( tkind ) :: h, layh


        ! Find local ID and layer number
        locID = kn
        lb0 = locv_NbLay( locID )
        layID = strat_layID(  locID ,  lb0  )

        if( locv_NbLay( locID ) == 1 ) return

        lb = lb0
        if( strat_layID(  locID ,  lb0 - 1  ) > 1 ) lb = lb0 - 1

10  continue

        ! In case a layer still exist
        if( strat_layID(  locID ,  lb  ) > 1 )then

            ! Get layer thickness
            layh = strat_thick(  locID ,  lb  )

            ! If not enough sediment just delete the layer
            if( layh <= tor )then
                strat_zelev(  locID ,  lb  ) = strat_zelev(  locID ,  lb - 1  )
                strat_thick(  locID ,  lb  ) = 0.0_8
                strat_sedh( locID, lb, 1:max_grn ) = 0.0_8
                strat_porosity(  locID ,  lb, 1:max_grn  ) = 0.0_8
                if( locv_NbLay( locID ) > 1 ) locv_NbLay( locID ) = lb - 1
                if( lb < lb0 )then
                    lb = lb0
                    goto 10
                endif
            endif

            ! Get total thickness of sediment layer based on sediment classes thicknesses
            h = 0.0_8
            do ks = 1, totgrn
                if( strat_sedh(  locID ,  lb ,  ks  ) < tor )then
                    strat_sedh(  locID ,  lb ,  ks  )  = 0.0_8
                    strat_porosity(  locID ,  lb, ks  )  = 0.0_8
                endif
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
                strat_porosity(  locID ,  lb, 1:max_grn  ) = 0.0_8
                if( locv_NbLay( locID ) > 1 ) locv_NbLay( locID ) = lb - 1
            endif

        endif

        ! Check we are at the top layer otherwise keep checking layers.
        if( lb < lb0 )then
            lb = lb0
            goto 10
        endif
        
        return

    end subroutine check_modified_sedimentary_layer
    ! ============================================================================
    !> Subroutine get_maximum_erosion_thickness
    !! Subroutine get_maximum_erosion_thickness get the maximum erosion thickness
    !! \param nid, erom
    !<
    ! ============================================================================
    subroutine get_maximum_erosion_thickness( nid, erom )

        integer :: nid, ke, kk, gid

        real( tkind ) :: szt, minz, erom

        gid = locv_gid( nid )
        erom = 0.0_8

        ! Current layer elevation
        szt = elev_record( gid )
        minz = szt

        ! Find neighborhood minimum elevation
        do ke = 1, 8
            kk = gv_ngbID( gid, ke )
            if(  kk > 0 ) minz = min( minz, elev_record( kk ) )
        enddo

        ! Erosion maximum to prevent formation of holes
        erom = szt - minz
        if( erom < 0.0_8 ) erom = 0.0_8

        if( szt > gsea%actual_sea .and. minz < gsea%actual_sea ) erom = szt - gsea%actual_sea
        erom = min( erom, transport%erosion_limit )

        return

    end subroutine get_maximum_erosion_thickness
    ! ============================================================================
    !> Subroutine update_flow_walkers_concentration_and_topography
    !! Subroutine update_flow_walkers_concentration_and_topography update the flow
    !! walkers variables.
    !<
    ! ============================================================================
    subroutine update_flow_walkers_concentration_and_topography

        integer :: k, ks, nblay, l, r, req, req2

        integer, dimension( mpi_status_size ) :: stat1, stat2

        real( tkind ) :: sndr( haloproc_nb ), sndl( haloproc_nb )

        ! Update flow walkers parameters
        call send_back_flow_walkers_concentration

        ! Update local stratigraphic elevation
        do k = 1, L_strat%nodes_nb
            ks = locv_gid( k )
            nblay = locv_NbLay( k )
            elev_record( ks ) = strat_zelev(  k ,  nblay  )
        enddo

        ! Broadcast halo top elevation to global grid
        l = 0
        r = 0
        sndr = -1.0e6_8
        sndl = -1.0e6_8
        do k = 1, haloborder_nb
            if( halo_send( k ) == iam - 1 )then
                ks = halo_lid( k )
                if( ks > 0 )then
                    l = l + 1
                    nblay = locv_NbLay( ks )
                    sndr( l ) = strat_zelev(  ks ,  nblay  )
                    sndr( l ) = strat_zelev(  ks ,  nblay  )
                 endif
            elseif( halo_send( k ) == iam + 1 )then
                ks = halo_lid( k )
                if( ks > 0 )then
                    r = r + 1
                    nblay = locv_NbLay( ks )
                    sndl( r ) = strat_zelev(  ks ,  nblay  )
                 endif
            endif
        enddo

        if( halo_ngb( 1 ) >= 0 )then
            call mpi_isend( sndr, haloproc_nb, dbl_type, halo_ngb( 1 ), 2121, &
                SPModel_comm_world, req, ierr )
            call mpi_request_free( req, ierr )
            call mpi_recv( halo1_elev, haloproc_nb, dbl_type, halo_ngb( 1 ), 2131, &
                SPModel_comm_world,  stat2, ierr )
        endif

        if( halo_ngb( 2 ) >= 0 )then
            call mpi_isend( sndl, haloproc_nb, dbl_type,halo_ngb( 2 ), 2131, &
                SPModel_comm_world, req2, ierr )
            call mpi_request_free( req2, ierr )
            call mpi_recv( halo2_elev, haloproc_nb, dbl_type, halo_ngb( 2 ), 2121, &
                SPModel_comm_world,  stat1, ierr )
        endif

        do k = 1, haloproc_nb
            ks = halo_rcv1( k )
            if( ks > 0 ) elev_record( ks ) = halo2_elev( k )
            ks = halo_rcv2( k )
            if( ks > 0 ) elev_record( ks ) = halo1_elev( k )
        enddo

        return

    end subroutine update_flow_walkers_concentration_and_topography
    ! ============================================================================
    !> Subroutine send_back_flow_walkers_concentration
    !! Subroutine send_back_flow_walkers_concentration is used to transfer FW concentration
    !! back to where it belongs.
    !<
    ! ============================================================================
    subroutine send_back_flow_walkers_concentration

        integer :: k, kr, fid, i
        integer :: outgoingCount, incomingCount

        integer, dimension( nproc ) :: sendCounts, sendOffsets, sendIdx
        integer, dimension( nproc ) :: recvCounts, recvOffsets
        integer, dimension( num_fwed ) :: sendTargets

        type( flowConc ), dimension(:), allocatable :: outgoingFWs
        type( flowConc ), dimension(:), allocatable :: incomingFWs
        type( flowConc ) :: tempFW

        real( tkind ) :: slsed
        real( tkind ), dimension( num_fw, totgrn ) :: newconc

        newconc = 0.0_8

        ! Run over all flow walkers, and check the location of the influenced nodes on the
        ! stratigraphy. Then work out which processor they belong on, keeping track of how
        ! many will be going to that processor in total.
        sendCounts = 0
        sendTargets = -1
        do k = 1, num_fwed
            if( ed_bound( k ) /= 1 )then
                ! If the FW belongs somewhere else.
                if( ed_procid( k ) /= iam )then
                    sendTargets( k ) = ed_procid( k )
                    sendCounts( ed_procid( k ) + 1 )  = sendCounts( ed_procid( k ) + 1 ) + 1
                else
                    refCheckLocal: do kr = 1, num_fw
                        if( ed_refid( k ) == fa_refid( kr ) )then
                            do i = 1, totgrn
                                newconc( kr, i ) = newconc( kr, i ) + concED( k, i )
                            enddo
                            exit refCheckLocal
                        endif
                    enddo refCheckLocal
                endif
            endif
        enddo

        ! Count how many particles I'm sending, and allocate my outgoing FW buffer.
        outgoingCount = 0
        do k = 1, nproc
            outgoingCount = outgoingCount + sendCounts( k )
        enddo
        allocate( outgoingFWs( outgoingCount ) )

        ! Make sure that every processor knows how many FWs it will be receiving from
        ! every other processor.
        recvCounts = 0
        call mpi_alltoall( sendCounts, 1, int_type, recvCounts, &
            1, int_type, SPModel_comm_world, ierr )

        ! Count how many new FWs I'll be handled, calculate offsets,
        ! and allocate my incoming FW buffer.
        incomingCount = 0
        recvOffsets = 0
        do k = 1, nproc
            ! Set the offset for data flowing from this process.
            recvOffsets( k ) = incomingCount

            ! Incorporate the incoming particles into our count.
            incomingCount = incomingCount + recvCounts( k )

            ! Sanity-check.
            if( k == iam + 1 .and. recvCounts( k ) /= 0 )then
                print *, "Something went wrong: trying to receive particles from myself..."
                stop
            endif
        enddo
        allocate( incomingFWs( incomingCount ) )

        !  Set up our outgoing offset array.
        sendOffsets = 0
        outgoingCount = 0
        do k = 1, nproc
            sendOffsets( k ) = outgoingCount
            outgoingCount = outgoingCount + sendCounts( k )
        enddo

        ! Marshall the FWs into the output array.
        sendIdx = sendOffsets
        do k = 1, num_fwed

            if( sendTargets( k ) >= 0 )then

                ! Move our offset along by one.
                sendIdx( sendTargets( k ) + 1 ) = sendIdx( sendTargets( k ) + 1 ) + 1

                ! Populate the FW datatype.
                tempFW%id = ed_refid( k )
                tempFW%conc = 0.0_8
                tempFW%conc( 1:totgrn ) =  concED( k, 1:totgrn )

                ! Copy the FW into the marshalling area.
                outgoingFWs( sendIdx(sendTargets( k ) + 1 ) ) = tempFW
            endif
        enddo

        ! Exchange FWs.
        call mpi_alltoallv( outgoingFWs, sendCounts, sendOffsets, mpi_concfw, &
            incomingFWs, recvCounts, recvOffsets, mpi_concfw, SPModel_comm_world, ierr )

        ! Unpack our received particles into our FW array.
        do k = 1, incomingCount
            fid = incomingFWs( k )%id
            refCheck: do kr = 1, num_fw
                if( fid == fa_refid( kr ) )then
                    do i = 1, totgrn
                        newconc( kr, i ) = newconc( kr, i ) + incomingFWs( k )%conc( i )
                    enddo
                    exit refCheck
                endif
            enddo refCheck
        enddo

        !  Change volume of sediments transfered during erosion and deposition
        do k = 1, num_fw

            if( fa_inbound( k ) /= 0 )then
                fa_sedcharge( k, 1:totgrn ) = 0.0_8
                fa_inbound( k ) = 1
            else
                slsed = 0.0_8
                fa_density( k ) = 0.0_8
                do i = 1, totgrn
                    fa_volume( k ) = fa_volume( k ) + newconc( k, i )
                enddo

                do i = 1, totgrn
                    fa_sedcharge( k , i ) = fa_sedcharge( k , i ) + newconc( k, i )
                    if( fa_sedcharge( k , i ) < 0.0_8 .and. fa_sedcharge( k , i ) < -1.e-4 )&
                        print*,'Warning: flow walker sediment flux was negative and set to zero.',fa_sedcharge( k , i )
                    if( fa_sedcharge( k , i ) < tor ) fa_sedcharge( k , i ) = 0.0_8
                    slsed = slsed + fa_sedcharge( k , i ) * sediment( i )%density
                enddo

                fa_density( k ) = slsed / fa_volume( k ) + fluid_density
                if( fa_type( k ) < tpplum1 )then
                    if( fa_density( k ) >= sea_density )then
                        if( fa_type( k ) < tpturb .and. fa_zpos( k ) <= gsea%actual_sea )&
                            fa_type( k ) = tpturb
                    endif
                endif

            endif

        enddo

        return

    end subroutine send_back_flow_walkers_concentration
    ! ============================================================================

end module mod_functionsfw

