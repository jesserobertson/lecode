! ============================================================================
! Name        : Porosity_init.f90
! Author      : tristan salles
! Created on: Aug 15, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
! ============================================================================
!> Module compaction
!!
!! Porosity is used to read porosity look-up table and allocate porosity values over the
!! mesh vertices. It is also used to modify mesh coordinates based on the induced compaction.
!! The porosity is computed using a user defined lookup table as proposed in Sedsim.
!!
!<
! ============================================================================
module mod_compaction

    use file_data
    use mpidata
    use flux_data
    use error_data
    use strata_data
    use forces_data
    use sediment_data

    implicit none

    private

    public :: cmpt_compaction, porosity_init, assign_porosity_table, update_top_layer_porosity

contains

    ! ============================================================================
    !> Subroutine cmpt_compaction()
    !! calculates new porosities for each sediment cell caused by
    !! compaction due to an additional deposited layer.
    !! This subroutine calculates the lithostatic pressure of the overlying column
    !! and calls porosity_function() which calculates the porosity as a function of grain
    !! size distribution and effective pressure.
    !! Then it moves the topographic elevation by the compactional subsidence
    !! caused by compaction of all layers below the surface.
    !<
    ! ============================================================================
    subroutine cmpt_compaction

        integer :: kn, kl, nbl
        integer :: laynb, gid
      
        real( tkind ) :: Plith

        ! Loop through surface faces and calculate porosity for the new top deposit layer.
        call compute_top_layer_porosity

        ! Update porosities of underlying column due to compaction of newly deposited sediment
        do kn = 1, L_strat%nodes_nb
            nbl = locv_NbLay( kn )
            laynb = strat_layID( kn, nbl )
            if( laynb == L_strat%layersID .and. strat_thick(  kn,  nbl  ) > 0.0_8 )then
                ! Lithostatic pressure
                Plith = 0.0_8
                ! Get porosity from top to bottom layers
                do kl = nbl - 1, 2, -1
                    call compute_underlying_porosity( kl, kn, Plith, nbl )
                enddo
            endif
        enddo

        ! Broadcast top elevation to global grid
        proc_elev = -1.0e6_8
        do kn = 1, L_strat%nodes_nb
            nbl = locv_NbLay( kn )
            gid = locv_gid( kn )
            proc_elev( gid ) = strat_zelev(  kn ,  nbl  )
        enddo
        call mpi_allreduce( proc_elev, elev_record, G_strat%nodes_nb, dbl_type, &
            max_type, SPModel_comm_world, ierr )

        return

    end subroutine cmpt_compaction
    ! ============================================================================
    !> Subroutine update_top_layer_porosity()
    !! update top layer porosities.
    !<
    ! ============================================================================
    subroutine  update_top_layer_porosity( k )

        integer :: k, laynb, ks, nbl

        real( tkind ) :: in_poro

        nbl = locv_NbLay( k )
        laynb = strat_layID( k, nbl )

        if( laynb > 1 )then

            ! Get porosity value and update sediment class and overall layer thicknesses
            strat_thick(  k ,  nbl  ) = 0.0_8
            do ks = 1, totgrn
                in_poro =  strat_porosity(  k, nbl, ks  )
                strat_porosity(  k, nbl, ks  ) = porosity_function( ks, 0.0_8, in_poro )
                ! In case the porosity was not set for the considered layer virtually increase
                ! the sediment thickness to take into account the void space based on the
                ! newly calculated porosity value
                if( in_poro == 0.0_8 )then
                    strat_sedh(  k ,  nbl ,  ks  ) = strat_sedh(  k ,  nbl ,  ks  ) * &
                        ( 1 + strat_porosity(  k , nbl, ks  ) )
                ! In case the porosity has decreased then decrease the sediment thickness
                ! accordingly to take into account the compaction of the layer.
                elseif( in_poro > strat_porosity(  k, nbl, ks  ) )then
                    strat_sedh(  k ,  nbl ,  ks  ) = strat_sedh(  k ,  nbl ,  ks  ) * &
                        ( 1 + strat_porosity( k, nbl, ks ) - in_poro )
                endif
                ! Update the layer thickness accordingly
                strat_thick(  k ,  nbl  ) = strat_thick(  k ,  nbl  ) + strat_sedh(  k ,  nbl ,  ks  )
            enddo

            ! Update top layer elevation
            strat_zelev(  k ,  nbl  ) = strat_zelev(  k,  nbl - 1 ) + &
                strat_thick(  k ,  nbl  )
        endif

        return

    end subroutine  update_top_layer_porosity
    ! ============================================================================
    !> Subroutine compute_top_layer_porosity()
    !! computes new layer porosities.
    !<
    ! ============================================================================
    subroutine  compute_top_layer_porosity

        integer :: k, laynb, ks, nbl

        real( tkind ) :: in_poro

        do k = 1, L_strat%nodes_nb
            nbl = locv_NbLay( k )
            laynb = strat_layID( k, nbl )

            if( L_strat%layersID == laynb )then

                ! Get porosity value and update sediment class and overall layer thicknesses
                strat_thick(  k ,  nbl  ) = 0.0_8
                do ks = 1, totgrn
                    in_poro =  strat_porosity(  k, nbl, ks  )
                    strat_porosity(  k, nbl, ks  ) = porosity_function( ks, 0.0_8, in_poro )
                    ! In case the porosity was not set for the considered layer virtually increase
                    ! the sediment thickness to take into account the void space based on the
                    ! newly calculated porosity value
                    if( in_poro == 0.0_8 )then
                        strat_sedh(  k ,  nbl ,  ks  ) = strat_sedh(  k ,  nbl ,  ks  ) * &
                            ( 1 + strat_porosity(  k ,  nbl, ks  ) )
                    ! In case the porosity has decreased then decrease the sediment thickness
                    ! accordingly to take into account the compaction of the layer.
                    elseif( in_poro > strat_porosity(  k, nbl, ks  ) )then
                        strat_sedh(  k ,  nbl ,  ks  ) = strat_sedh(  k ,  nbl ,  ks  ) * &
                            ( 1 + strat_porosity( k, nbl, ks ) - in_poro )
                    endif
                    ! Update the layer thickness accordingly
                    strat_thick(  k ,  nbl  ) = strat_thick(  k ,  nbl  ) + strat_sedh(  k ,  nbl ,  ks  )
                enddo

                ! Update top layer elevation
                strat_zelev(  k ,  nbl  ) = strat_zelev(  k,  nbl - 1 ) + &
                    strat_thick(  k ,  nbl  )
            endif

        enddo

        return

    end subroutine  compute_top_layer_porosity
    ! ============================================================================
    !> Subroutine compute_underlying_porosity()
    !! calculates porosities and compaction values for vertices which coincided with top ones.
    !! \param kn, kl, nid, Plith, mass, nbl
    !<
    ! ============================================================================
    subroutine  compute_underlying_porosity( kl, nid, Plith, nbl )

        integer :: ks, kl, kll, nid, nbl

        real( tkind ) :: laysub, thick, in_poro
        real( tkind ) :: mass, Plith

        mass = 0.0_8
        do ks = 1, totgrn
            ! Compute upper layer sediment mass
            mass = mass + strat_sedh(  nid,  kl + 1,  ks  ) * &
                sediment( ks )%density * ( 1 - strat_porosity(  nid,  kl + 1, ks  ) )
            ! Compute mass due to water in sediment pore space
            mass = mass + strat_sedh(  nid,  kl + 1,  ks  ) * &
                sea_density * strat_porosity(  nid,  kl + 1, ks  )
        enddo

        ! If cell has sediment, increment cumulative lithostatic pressure
        if( mass > tor )then

            ! Compute recursively lithostatic pressure for current sediment layer depth
            Plith = Plith + 1.0e-6_8 * mass * gravity

            strat_thick(  nid ,  kl  ) = 0.0_8
            laysub=0.0_8

            do ks = 1, totgrn
                in_poro =  strat_porosity(  nid, kl, ks  )
                thick = strat_sedh(  nid,  nbl ,  ks  )
                strat_porosity(  nid, kl, ks  ) = porosity_function( ks, Plith, in_poro )
                ! In case the porosity was not set for the considered layer virtually increase
                ! the sediment thickness to take into account the void space based on the
                ! newly calculated porosity value
                if( in_poro == 0.0_8 )then
                    strat_sedh(  nid ,  nbl ,  ks  ) = strat_sedh(  nid,  nbl ,  ks  ) * &
                        ( 1 + strat_porosity(  nid,  nbl, ks  ) )
                ! In case the porosity has decreased then decrease the sediment thickness
                ! accordingly to take into account the compaction of the layer.
                elseif( in_poro > strat_porosity(  nid, nbl, ks  ) )then
                    strat_sedh( nid,  nbl ,  ks  ) = strat_sedh(  nid,  nbl ,  ks  ) * &
                        ( 1 + strat_porosity( nid, nbl, ks ) - in_poro )
                endif
                laysub = laysub + thick - strat_sedh(  nid,  nbl ,  ks  )
                strat_thick(  nid ,  kl  ) = strat_thick(  nid ,  kl  ) + strat_sedh(  nid,  nbl ,  ks  )
                strat_zelev(  nid ,  kl  ) = strat_zelev(  nid ,  kl - 1  ) + strat_thick(  nid ,  kl  )
            enddo

            ! Recursively update layer elevations due to subsidence
            if( laysub > tor )then
                do kll = kl + 1, nbl
                    strat_zelev(  nid ,  kll  ) = strat_zelev(  nid ,  kll  ) - laysub
                enddo
            endif

        endif

        return

    end subroutine  compute_underlying_porosity
    ! ============================================================================
    !> Subroutine porosity_init()
    !! Calculates initial deposit porosities for each sediment cell.
    !<
    ! ============================================================================
    subroutine porosity_init

        integer :: k, kl, ks, nid, nbl
        integer :: laynb

        real( tkind ) :: in_poro, Plith, mass

        ! Loop through surface faces and calculate porosity for the initial deposit layers.
        do k = 1, L_strat%nodes_nb
            nbl = locv_NbLay( k )
            nid = k
            Plith = 0.0_8
            do kl = nbl, 2, -1
                mass = 0.0_8
                laynb = strat_layID(  k ,  kl  )
                if( L_strat%layersID <= laynb .and. laynb > 1  )then
                    if( L_strat%layersID < laynb )then
                        do ks = 1, totgrn
                            mass = mass + strat_sedh(  nid ,  kl + 1 ,  ks  ) * &
                                sediment( ks )%density * ( 1 - strat_porosity(  nid ,  kl + 1, ks  ) )
                            ! Mass due to water in pore space
                            mass = mass + strat_sedh(  nid ,  kl + 1 ,  ks  ) * &
                                sea_density * strat_porosity(  nid ,  kl + 1, ks  )
                        enddo
                    endif
                    Plith = Plith + 1.0e-6_8 * mass * gravity
                    do ks = 1, totgrn
                        in_poro =  strat_porosity(  k, nbl, ks  )
                        strat_porosity(  nid ,  kl, ks  ) = porosity_function( ks, Plith, in_poro )
                    enddo
                endif
            enddo
        enddo

        call mpi_barrier(SPModel_comm_world,ierr)

        return

    end subroutine porosity_init
    ! ============================================================================
    !> Subroutine porosity_function()
    !! Computes a variable porosity based upon effective subsurface pressure.
    !! \param dz_poro, prate, in_poro, out_poro
    !<
    ! ============================================================================
    function porosity_function( grn, LithPressure, previous_porosity ) result( new_porosity )

        integer :: grn, k

        real( tkind ) :: LithPressure, previous_porosity, new_porosity, interpolate_porosity

        new_porosity = previous_porosity

        loop_over_pressure: do k = 2, gporo%ePnb
            if( effPressure( k - 1 ) <= LithPressure .and. &
                effPressure( k ) > LithPressure ) exit loop_over_pressure
        enddo loop_over_pressure

        ! Calculate interpolated porosity based on effective subsurface pressure.
        interpolate_porosity = ( porosity( grn, k ) - porosity( grn, k - 1 ) ) / ( effPressure( k ) - effPressure( k - 1 ) )
        interpolate_porosity = interpolate_porosity * ( LithPressure - porosity( grn, k - 1 ) )
        interpolate_porosity = interpolate_porosity + porosity( grn, k - 1 )

        if( interpolate_porosity < previous_porosity .or. previous_porosity == 0.0_8 )&
            new_porosity = interpolate_porosity

    end function porosity_function
    ! ============================================================================
    !> Subroutine assign_porosity_table()
    !! Assigns the porosity look-up table values from the input file.
    !<
    ! ============================================================================
    subroutine assign_porosity_table

        logical :: found
        integer :: i, maxl, l, d
        integer :: iu, ios
        character( len=200 ) :: line

        ! Open porosity look-up table file
        iu = 80 + iam
        inquire( file=fporosity, exist=found )
        if(found)then
            open(iu,file=fporosity,status="old",action="read",iostat=ios)
            rewind(iu)
        else
            write(*,*)'Warning: the input file for porosity looku-up table cannot be found'
            write(*,*)' The code is looking for the following name:',trim(fporosity)
            stop
        endif

        ! Number of lines in the porosity file
        maxl = 1 + totgrn
        l = 0
        d = 0

        ! Read parameters
        do while ( l < maxl )
            read(iu,'(a200)') line
            i = len_trim( line )
            if( line(i:i) == char(13) ) i = i-1
            if( line(1:2) /= '* ')then
                ! Read effective pressure
                if( l == 0 )then
                    read( line(2:i),* ) effPressure( 1:gporo%ePnb )
                ! Read porosity values
                elseif( l < maxl - 1)then
                    d = d + 1
                    read( line(2:i),* ) porosity( d,1:gporo%ePnb )
                endif
                l = l + 1
            endif
        enddo

        ! Close file
        close( iu )

        if( d /= totgrn ) print*,'Something went wrong when reading porosity look-up table',d,totgrn

        return

    end subroutine assign_porosity_table
  ! ============================================================================

end module mod_compaction
