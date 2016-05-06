! ============================================================================
! Name        : ExtForces_init.f90
! Author      : Tristan Salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file ExtForces_init.f90
!!
!! ExtForces_init is used to initialise external forces (sea-level/displacement)
!!
!<
! ============================================================================
module mod_displacements

    use file_data
    use mpidata
    use time_data
    use error_data
    use strata_data
    use forces_data
    use param_data

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine assign_vertical_displacements_rate()
    !! Performs a simple vertical displacement fields over the nodes of the stratal mesh.
    !! It generates the vertical_displacements field which contains the
    !! displacements in the Z direction for the simulated period.
    !! \param period
    !<
    ! ============================================================================
    subroutine assign_vertical_displacements_rate( period )

        ! Parameters Declaration
        logical :: found

        integer,intent( in ) :: period

        integer :: k, gid, i, kn
        integer :: iu, ios

        if( .not. allocated( gdisp_fill ) ) allocate( gdisp_fill( gdisp%event ) )

        ! First displacement period allocated required arrays
        if( period == 1 )then
            if( gdisp_time( 1, 1 ) > time_start ) vdisp%event = 1

            ! Find the vertical displacement event number
            do k = 1, gdisp%event
                vdisp%event = vdisp%event + 1
                if( gdisp_time( k, 1 ) >= gdisp_time( k, 2 ) ) attempt = DISP_TIME1
                if( k < gdisp%event )then
                    if( gdisp_time( k, 2 ) > gdisp_time( k + 1, 1 ) ) attempt = DISP_TIME2
                    if( gdisp_time( k, 2 ) /= gdisp_time( k + 1, 1 ) ) vdisp%event = vdisp%event + 1
                endif
            enddo

            ! Ensure that all processes have successfully completed subroutine.
            call completion

            ! If last displacement period stops before the end of the simulation
            ! add a new event during the last period
            if( gdisp_time( gdisp%event, 2 ) < time_end ) vdisp%event = vdisp%event + 1

            ! In case the number of vertical displacement is different from the number of
            ! of tectonic events defined in the XmL input file defines the new event parameters
            if( vdisp%event /= gdisp%event )then

                allocate( vdisp_time( vdisp%event, 2 ), vdisp_fill( vdisp%event ) )
                kn = 1
                if( gdisp_time( 1, 1 ) > time_start )then
                    vdisp_time( kn, 1 ) = time_start
                    vdisp_time( kn, 2 ) = gdisp_time( 1, 1 )
                    vdisp_fill( kn ) = 0
                    kn = kn + 1
                endif

                do k = 1, gdisp%event
                    vdisp_time( kn, 1 ) = gdisp_time( k, 1 )
                    vdisp_time( kn, 2 ) = gdisp_time( k, 2 )
                    vdisp_fill( kn ) = k
                    kn = kn + 1
                    if( k < gdisp%event )then
                        if( gdisp_time( k, 2 ) /= gdisp_time( k + 1, 1 ) )then
                                vdisp_time( kn, 1 ) = gdisp_time( k, 2 )
                                vdisp_time( kn, 2 ) = gdisp_time( k + 1, 1 )
                                vdisp_fill( kn ) = 0
                                kn = kn + 1
                        endif
                    endif
                enddo

                if( gdisp_time( gdisp%event, 2 ) < time_end )then
                    vdisp_time( kn, 2 ) = time_end
                    vdisp_time( kn, 1 ) = gdisp_time( gdisp%event, 2 )
                    vdisp_fill( kn ) = 0
                endif
                deallocate( gdisp_time, gdisp_fill )
                allocate( gdisp_time( vdisp%event, 2 ), gdisp_fill( vdisp%event ) )

                gdisp%event =  vdisp%event
                do k = 1,  vdisp%event
                    gdisp_time( k, 1 ) = vdisp_time( k, 1 )
                    gdisp_time( k, 2 ) = vdisp_time( k, 2 )
                    gdisp_fill( k ) = vdisp_fill( k )
                enddo
                deallocate( vdisp_time, vdisp_fill )

            ! Otherwise define the displacement event number directly
            else
                do k = 1, gdisp%event
                    gdisp_fill( k ) = k
                enddo
            endif

        endif

        ! Find actual event number
        gdisp%actual = 0
        do k = 1, gdisp%event
            if( tnow < gdisp_time( k, 2 ) .and. tnow >= gdisp_time( k, 1 ) ) gdisp%actual = k
        enddo

        ! Open displacements field file for the considered event number
        kn = gdisp%actual
        if( gdisp_fill( kn ) > 0 )then
            inquire( file=fdisp( gdisp_fill( kn )), exist=found )
            if( .not. found ) attempt = DISP_FILE
            call completion
            ! Allocate displacement values
            iu = 80 + iam
            open(iu,file=fdisp( gdisp_fill( kn )),status="old",action="read",iostat=ios)
            rewind( iu )
            ! Assign coordinate array
            uplift = 0.0_8
            do i = 1, G_strat%nodes_nb
                read( iu, * ) k, gv_vdisp( k )
            enddo
            close( iu )
        ! If this event doesn't correspond to any defined one initialise displacement values to 0
        else
            gv_vdisp = 0.0_8
        endif
        
        ! Allocate local processor displacement rate (m/s) values
        do k = 1, L_strat%nodes_nb
            gid = locv_gid( k )
            locv_vdisp( k ) = real( gv_vdisp( gid ) )
        enddo
        
        return

    end subroutine assign_vertical_displacements_rate
    ! ============================================================================
    !> Subroutine update_local_displacements()
    !! Update displacement on the mesh node according to the displacement fields.
    !<
    ! ============================================================================
    subroutine update_local_displacements

        ! Parameters Declaration
        integer :: k
        integer :: period

        ! Take displacement event period number
        period = gdisp%actual

        ! Assign array
        do k =1, L_strat%nodes_nb
            ! Compute displacement in metres
            locv_vrate( k ) = dble( locv_vdisp( k ) / &
                ( (gdisp_time( period, 2 ) - gdisp_time( period, 1 ) ) * secyear ) )
        enddo

        return

    end subroutine update_local_displacements
    ! ============================================================================
    !> Subroutine apply_displacements()
    !! Apply displacement on the mesh node according to the displacement fields.
    !! This will read the previous SPModel mesh and create a new one based on
    !! displacements values and previous node and elements connectivity.
    !<
    ! ============================================================================
    subroutine apply_displacements

        ! Parameters Declaration
        integer :: k, n, GID, nb_lay

        ! Get the elapsed time since previous displacement call
        gdisp%disptime = ( tnow - gdisp%lastdisp )* secyear

        ! Apply the displacement to the deposit layers
        do k = 1, L_strat%nodes_nb
            ! First update the local stratal grid basement
            locv_coord( k, 3 ) = real( locv_coord( k, 3 ) + &
                locv_vrate( k ) * gdisp%disptime )
            ! For each stratal point performs the displacement
            do n = 1, locv_NbLay( k )
                strat_zelev(  k ,  n  ) = real( strat_zelev(  k ,  n  ) + &
                    locv_vrate( k ) * gdisp%disptime )
            enddo
        enddo

        ! Broadcast top elevation to global grid
        proc_elev = -1.0e6_8
        elev_record = -1.0e6_8
        do n = 1, L_strat%nodes_nb
            nb_lay = locv_NbLay( n )
            GID = locv_gid( n )
            proc_elev( GID ) = strat_zelev(  n ,  nb_lay  )
        enddo

        call mpi_allreduce( proc_elev, elev_record, G_strat%nodes_nb, dbl_type, &
            max_type, SPModel_comm_world, ierr )

        return

    end subroutine apply_displacements
    ! ============================================================================
    !> Subroutine load_displacements()
    !! load_displacements is used to load a new displacement field based on user defined
    !! displacement data.
    !! \param period
    !<
    ! ============================================================================
    subroutine load_displacements( period )

        ! Parameters Declaration
        integer,intent( in ) :: period

        ! Get new displacement fields
        call assign_vertical_displacements_rate( period )

        ! Update displacement locally
        call update_local_displacements

        ! Set displacement flag to false as we've just updated it
        new_disp = .false.

        return

    end subroutine load_displacements
   ! ============================================================================

end module mod_displacements
! ============================================================================
module mod_sealevel

    use file_data
    use mpidata
    use time_data
    use error_data
    use param_data
    use forces_data

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine read_sealevel_file()
    !! Reads formatted data of sea level fluctuations.
    !<
    ! ============================================================================
    subroutine read_sealevel_file

        ! Parameters Declaration
        logical :: found
        integer :: i, k, i2
        integer :: iu, ios

        character(len=128) :: text

        ! Find and open the sea level file
        iu = 45 + iam
        inquire( file=fsea, exist=found )
        if(found)then
            open(iu,file=fsea,status="old",action="read",iostat=ios)
            rewind(iu)
        else
            print*,'Cannot find sea level file.'
            attempt = SEA_FILE
        endif

        ! Find the number of events
        k = 0
        do
            read(iu,'(a128)',iostat=ios,end=30) text
            k = k + 1
        enddo
30      continue
        rewind(iu)

        ! Ensure that all processes have successfully read the file.
        call completion

        ! Define the number of sea fluctuation events & allocate array
        gsea%event = k
        allocate( sealvl( 2, gsea%event ) )

        ! Sea-level fluctuations
        do k = 1, gsea%event
            read(iu, '(a128)' ,iostat=ios, end=60 )text
            i2=len_trim(text)
            if( text( i2:i2 ) == char(13) )then
                i2=i2-1
            endif
            read(text(1:i2),*,iostat=ios)( sealvl( i, k ), i = 1, 2 )
            if( ios /= 0 ) attempt = SEA_EOF
            if( k > 1 )then
                if( sealvl( 1, k ) <= sealvl( 1, k - 1 ) )&
                    attempt = SEA_ORDER
            endif
        enddo

        ! Ensure that sea level files has been properly ordered.
        call completion

60      close(iu)

        gsea%last_sea = 0.0_8

        ! Update the sea level value
        call eustatism

        return

    end subroutine read_sealevel_file
    ! ============================================================================
    !> Subroutine eustatism()
    !! Performs the change of water elevation according to sea-level fluctuations.
    !<
    ! ============================================================================
    subroutine eustatism

        ! Parameters Declaration
        integer :: i
        type( sl_par ) :: sl_param

        if( .not. gsea%sealevel ) return

        gsea%last_sea = gsea%actual_sea

        i = 1
        ! Find time for the considered event
        do
            if( i > gsea%event )exit
            if( tnow < sealvl( 1, i ) )exit
            sl_param%time1 = sealvl( 1, i )
            sl_param%sea1 = sealvl( 2, i )
            i = i + 1
        enddo

        ! Linear interpolation of sea level elevation based on current time
        if( i <= gsea%event )then
            sl_param%time2 = sealvl( 1, i )
            sl_param%sea2 = sealvl( 2, i )
            sl_param%tsim = tnow
            call sealvl_interpolation( sl_param )
            gsea%actual_sea = sl_param%sl_int
        ! If event greater than max defined allocate with last value
        else
            gsea%actual_sea = sl_param%sea1
            sl_param%time2 = sealvl( 1, i - 1 )

        endif

        ! Declare the next sea level fluctuation update
        if( ( sl_param%time2 < next_coarse_dt ) .and. ( sl_param%time2 - tnow ) > time_tolerance )&
            next_coarse_dt = sl_param%time2 + tor

        return

    end subroutine eustatism
    ! ============================================================================
    !> Subroutine sealvl_interpolation()
    !! Performs the sea-level interpolation function.
    !!
    !! \arg \c use precision_data
    !!
    !! \param sl_param
    !<
    ! ============================================================================
    subroutine sealvl_interpolation( sl_param )

        ! Parameters Declaration
        type( sl_par ) :: sl_param

        sl_param%sl_int = ( sl_param%tsim - sl_param%time1 ) / ( sl_param%time2 - sl_param%time1 )
        sl_param%sl_int = sl_param%sl_int * ( sl_param%sea2 - sl_param%sea1 )
        sl_param%sl_int = sl_param%sl_int + sl_param%sea1

        return

    end subroutine sealvl_interpolation
    ! ============================================================================

end module mod_sealevel
! ============================================================================
module mod_tempsal

    use file_data
    use mpidata
    use time_data
    use error_data
    use param_data
    use forces_data

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine read_temperature_salinity_file()
    !! Reads formatted data of ocean temperature and salinity fluctuations.
    !<
    ! ============================================================================
    subroutine read_temperature_salinity_file

        ! Parameters Declaration
        logical :: found
        integer :: i, k, i2
        integer :: iu, ios

        character(len=128) :: text

        ! Find and open the temperature and salinity file
        iu = 45 + iam
        inquire( file=ftempsal, exist=found )
        if(found)then
            open(iu,file=ftempsal,status="old",action="read",iostat=ios)
            rewind(iu)
        else
            print*,'Cannot find temperature and salinity file.'
            attempt = SEA_FILE
        endif

        ! Find the number of events
        k = 0
        do
            read(iu,'(a128)',iostat=ios,end=30) text
            k = k + 1
        enddo
30      continue
        rewind(iu)

        ! Ensure that all processes have successfully read the file.
        call completion

        ! Define the number of temperature and salinity fluctuation events & allocate array
        gtempsal%event = k
        allocate( tempsal( 3, gtempsal%event ) )

        ! Temperature and salinity fluctuations
        do k = 1, gtempsal%event
            read(iu, '(a128)' ,iostat=ios, end=60 )text
            i2=len_trim(text)
            if( text( i2:i2 ) == char(13) )then
                i2=i2-1
            endif
            read(text(1:i2),*,iostat=ios)( tempsal( i, k ), i = 1, 3 )
            if( ios /= 0 ) attempt = SEA_EOF
            if( k > 1 )then
                if( tempsal( 1, k ) <= tempsal( 1, k - 1 ) )&
                    attempt = SEA_ORDER
            endif
        enddo

        ! Ensure that temperature and salinity file has been properly ordered.
        call completion

60      close(iu)

        gtempsal%actual_tempsal = 0.0_8
        gtempsal%last_tempsal = gtempsal%actual_tempsal

        ! Update the temperature and salinity values
!        call get_temperature_salinity

        return

    end subroutine read_temperature_salinity_file
    ! ============================================================================
    !> Subroutine get_temperature_salinity()
    !! Performs the change of water temperature and salinity.
    !<
    ! ============================================================================
    subroutine get_temperature_salinity

        ! Parameters Declaration
        integer :: i
        type( st_par ) :: ts_param

        if( .not. gtempsal%tpsl ) return

        gtempsal%last_tempsal = gtempsal%actual_tempsal

        i = 1
        ! Find time for the considered event
        do
            if( i > gtempsal%event )exit
            if( tnow < tempsal( 1, i ) )exit
            ts_param%time1 = tempsal( 1, i )
            ts_param%ts1 = tempsal( 2:3, i )
            i = i + 1
        enddo

        ! Linear interpolation of sea level elevation based on current time
        if( i <= gtempsal%event )then
            ts_param%time2 = tempsal( 1, i )
            ts_param%ts2 = tempsal( 2:3, i )
            ts_param%tsim = tnow
            call temperature_salinity_interpolation( ts_param )
            gtempsal%actual_tempsal = ts_param%ts_int
        ! If event greater than max defined allocate with last value
        else
            gtempsal%actual_tempsal = ts_param%ts1
            ts_param%time2 = tempsal( 1, i - 1 )

        endif

        ! Declare the next sea level fluctuation update
        if( ( ts_param%time2 < next_coarse_dt ) .and. ( ts_param%time2 - tnow ) > time_tolerance )&
            next_coarse_dt = ts_param%time2 + tor

        return

    end subroutine get_temperature_salinity
    ! ============================================================================
    !> Subroutine temperature_salinity_interpolation()
    !! Performs the temperature and salinity interpolation function.
    !!
    !! \arg \c use precision_data
    !!
    !! \param ts_param
    !<
    ! ============================================================================
    subroutine temperature_salinity_interpolation( ts_param )

        ! Parameters Declaration
        type( st_par ) :: ts_param

        ! Temperature
        ts_param%ts_int( 1 ) = ( ts_param%tsim - ts_param%time1 ) / ( ts_param%time2 - ts_param%time1 )
        ts_param%ts_int( 1 ) = ts_param%ts_int( 1 ) * ( ts_param%ts2( 1 ) - ts_param%ts1( 1 ) )
        ts_param%ts_int( 1 ) = ts_param%ts_int( 1 ) + ts_param%ts1( 1 )

        ! Salinity
        ts_param%ts_int( 2 ) = ( ts_param%tsim - ts_param%time1 ) / ( ts_param%time2 - ts_param%time1 )
        ts_param%ts_int( 2 ) = ts_param%ts_int( 2 ) * ( ts_param%ts2( 2 ) - ts_param%ts1( 2 ) )
        ts_param%ts_int( 2 ) = ts_param%ts_int( 2 ) + ts_param%ts1( 2 )

        return

    end subroutine temperature_salinity_interpolation
    ! ============================================================================

end module mod_tempsal
! ============================================================================
module mod_hemipelagic

    use file_data
    use mpidata
    use time_data
    use error_data
    use param_data
    use forces_data
    use sediment_data
    use strata_data
    use flux_data

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine read_hemilpelagic_file()
    !! Reads formatted data of hemipelagic sedimentation rates.
    !<
    ! ============================================================================
    subroutine read_hemilpelagic_file

        ! Parameters Declaration
        logical :: found
        integer :: i, k, i2
        integer :: iu, ios

        character(len=128) :: text

        ! Find and open the sea level file
        iu = 45 + iam
        inquire( file=fhemi, exist=found )
        if(found)then
            open(iu,file=fhemi,status="old",action="read",iostat=ios)
            rewind(iu)
        else
            print*,'Cannot find hemipelagic file.'
            attempt = SEA_FILE
        endif

        ! Find the number of hemipelagic sedimentation rate
        k = 0
        do
            read(iu,'(a128)',iostat=ios,end=30) text
            k = k + 1
        enddo
30      continue
        rewind(iu)

        ! Ensure that hemipelagic file has been found.
        call completion

        ! Allocate hemipelagic array
        contourhemi_nb = k
        if( allocated( hemipelagic ) ) deallocate( hemipelagic )
        allocate( hemipelagic( 2, contourhemi_nb ) )

        ! Read hemipelagic rate for each contour
        do k = 1, contourhemi_nb
            read(iu, '(a128)' ,iostat=ios, end=60 )text
            i2=len_trim(text)
            if( text( i2:i2 ) == char(13) )then
                i2=i2-1
            endif
            read(text(1:i2),*,iostat=ios)( hemipelagic( i, k ), i = 1, 2 )
            if( ios /= 0 )then
                print*,'hemipelagic file reading problem'
                stop
            endif
            if( k > 1 )then
                if( hemipelagic( 1, k-1 ) <= hemipelagic( 1, k ) )then
                    print*,'hemipelagic should be sorted in decreasing contour value'
                    stop
                endif
            endif
        enddo

60      close(iu)

        return

    end subroutine read_hemilpelagic_file
    ! ============================================================================
    !> Subroutine hemipelagites()
    !! Performs the sedimentation based on hemipelagic rate.
    !<
    ! ============================================================================
    subroutine hemipelagites

        integer :: k, i, gid, lb, layID

        real( tkind ) :: hemi_rate1, hemi_contour1, hemi_rate2, hemi_contour2
        real( tkind ) :: ztop, hemi_contour, hemi_rate, hemi_time, hemi_dep


        proc_elev = -1.0e6_8

        do k = 1, L_strat%nodes_nb

            gid = locv_gid( k )

            ! Sea bed elevation
            ztop = elev_record( gid )
            proc_elev( gid ) = elev_record( gid )

            ! As long as we are under water
            if( ztop < gsea%actual_sea )then
                i = 1
                ! Find the corresponding hemipelagic rate for the
                ! considered bathymetry.
                pick_hemipelagic_contour: do
                    if( i > contourhemi_nb )exit
                    if( ztop - gsea%actual_sea > hemipelagic( 1, i ) ) &
                        exit pick_hemipelagic_contour
                    hemi_contour1 = hemipelagic( 1, i )
                    hemi_rate1 = hemipelagic( 2, i )
                    i = i + 1
                enddo pick_hemipelagic_contour

                ! Update hemipelagic rate based on contour number
                hemi_rate = 0.0_8
                if( i <= contourhemi_nb )then
                    hemi_contour2 = hemipelagic( 1, i )
                    hemi_rate2 = hemipelagic( 2, i )
                    hemi_contour = ztop - gsea%actual_sea
                    ! Interpolate hemipelagic rate based on distance between current
                    ! depth and contours
                    call hemi_interpolation( hemi_contour, hemi_rate1, hemi_contour1, &
                        hemi_rate2, hemi_contour2, hemi_rate )
                else
                    hemi_rate = hemipelagic( 2, contourhemi_nb )
                endif

                ! Elapsed time since last call to hemipelagite function
                hemi_time = tnow - last_hemi

                ! Get the deposit thickness from hemipelagic rate ( m/yr )
                hemi_dep = hemi_rate * hemi_time

                ! In case there is sufficient deposit, add a new layer
                if( hemi_dep > 0.001_8 )then

                    lb = locv_NbLay( k )
                    layID = strat_layID(  k, lb )

                    ! Create a new layer in case a current layer doesn't exist
                    if( layID < L_strat%layersID )then
                        lb = lb + 1
                        strat_sedh( k, lb, 1:max_grn ) = 0.0_8
                        locv_NbLay( k ) = lb
                        strat_layID(  k ,  lb  ) = L_strat%layersID
                        strat_thick(  k ,  lb  ) = hemi_dep
                        strat_zelev(  k ,  lb  ) = hemi_dep + &
                            strat_zelev(  k ,  lb - 1  )
                        strat_porosity(  k ,  lb, 1:max_grn  ) = 0.0_8
                        strat_hardness(  k ,  lb  ) = 1.0_8
                        strat_sedh(  k ,  lb ,  hemi_mat  ) = hemi_dep
                        if( gporo%compaction ) &
                            strat_porosity(  k,  lb, hemi_mat ) = porosity( hemi_mat, 1 )

                    ! Otherwise add on an existing layer
                    else
                        strat_thick(  k ,  lb  ) = hemi_dep + &
                            strat_thick(  k ,  lb )
                        strat_zelev(  k ,  lb  ) = hemi_dep + &
                            strat_zelev(  k ,  lb )
                        strat_sedh(  k ,  lb ,  hemi_mat  ) = hemi_dep + &
                            strat_sedh(  k ,  lb ,  hemi_mat  )
                        if( gporo%compaction ) &
                            strat_porosity(  k,  lb, hemi_mat ) = porosity( hemi_mat, 1 )
                    endif

                    ! Update top elevation
                    proc_elev( gid ) = strat_zelev(  k ,  lb  )
                endif

            endif
        enddo

        ! Update global elevation
        call mpi_allreduce( proc_elev, elev_record , G_strat%nodes_nb, &
            dbl_type, max_type, SPModel_comm_world, ierr )

        ! Update last call to hemipelagite function to actual time
        last_hemi = tnow

        return

    end subroutine hemipelagites
    ! ============================================================================
    !> Subroutine hemi_interpolation()
    !! Performs the hemipelagic sedimentation rate interpolation function.
    !!
    !! \arg \c use precision_data
    !!
    !! \param hemi_contour, hemi_rate1, hemi_contour1, hemi_rate2, hemi_contour2, hemi_rate
    !<
    ! ============================================================================
    subroutine hemi_interpolation( hemi_contour, hemi_rate1, hemi_contour1, hemi_rate2, &
        hemi_contour2, hemi_rate )

        real( tkind ) :: hemi_contour, hemi_rate1, hemi_contour1, hemi_rate2, hemi_contour2, hemi_rate

        hemi_rate = ( hemi_contour - hemi_contour1 ) / ( hemi_contour2 - hemi_contour1 )
        hemi_rate = hemi_rate * ( hemi_rate2 - hemi_rate1 )
        hemi_rate = hemi_rate + hemi_rate1

        return

    end subroutine hemi_interpolation
    ! ============================================================================

end module mod_hemipelagic
! ============================================================================

