! ============================================================================
! Name        : FillDepressionDEM.f90
! Author      : Tristan Salles
! Created on: May 14, 2013
! Copyright (C) 2013 CSIRO
! ============================================================================
!> \file  FillDepressionDEM.f90
!!
!! Description : Sequential DEM preprocessing algorithm
!!
!! Reference:
!! Planchon, O., Darboux, F., 2001. A fast, simple, and versatile algorithm to fill the depressions
!! of digital elevation models. Catena 46(2003), 159-176.
!!
!<
! ============================================================================
module mod_filldepression

    use file_data
    use mpidata
    use strata_data
    use TIN_data
    use sediment_data
    use forces_data
    use param_data
    use flux_data
    use time_data

    implicit none

    public

    integer, dimension( 8 ) :: direct8x = (/ 1, 1, 1, 0, -1, -1, -1, 0 /)
    integer, dimension( 8 ) :: direct8y = (/ -1, 0, 1, 1, 1, 0, -1, -1 /)


    real( tkind ), dimension( : ), allocatable :: slp, curv

contains

    ! ============================================================================
    !> Subroutine planchon_dem_fill_algorithm
    !! No sink algorithm.
    !<
    ! ============================================================================
    subroutine planchon_dem_fill_algorithm

        logical :: flag
        integer :: p, k, k2, i, j, row, col

        flag = .true.

        ! In case we need to fill holes within the DEM
        if( step_fill > 0.0_8 )then

            ! Update DEM borders elevation values
            ! We add a ghost border on each face of the simulation area.
            p = 1
            do k = 1, G_strat%nodes_nb
                ! On East West border fix DEM elevation
                if( p == 1 .or. p == strat_X )then
                    wdem( k ) = elev_record( k )
                    if( p == 1 ) p = p + 1
                    if( p == strat_X ) p = 1
                ! On South border
                elseif( k <= strat_X )then
                    wdem( k ) = elev_record( k )
                    p = p + 1
                ! On the North border
                elseif( k >= G_strat%nodes_nb - strat_X )then
                    wdem( k ) = elev_record( k )
                    p = p + 1
                ! Else give high values
                else
                    wdem( k ) = 100000.0_8
                    p = p + 1
                endif

            enddo

            ! Now find the sinks and fill them usng Planchon's method
            do while( flag )
                flag = .false.

                do j = 2, strat_Y - 1
                    do i = 2, strat_X - 1
                        k = ( j - 1 ) * strat_X + i
                        ! In case DEM elevation is greater than top stratal grid elevation
                        if( wdem( k ) > elev_record( k ) )then
                            ! Look at the neighbors
                            do p = 1, 8

                                col = i + direct8x( p )
                                row = j + direct8y( p )
                                k2 = ( row - 1 ) * strat_X + col
                                ! In case stratal grid elevation greater than neighbors DEM elevation
                                if( elev_record( k ) >= wdem( k2 ) + step_fill )then
                                    wdem( k ) = elev_record( k )
                                    flag = .true.
                                ! Otherwise it is a sink and we perform sink filling
                                else
                                    if( wdem( k ) > wdem( k2 ) + step_fill )then
                                        wdem( k ) = wdem( k2 ) + step_fill
                                        flag = .true.
                                    endif
                                endif

                            enddo

                        endif

                    enddo
                enddo
            enddo

        ! In case the sink filling algorithm is not used just update the DEM elevation to
        ! the top stratal elevation
        else
            wdem = elev_record

        endif

        return

    end subroutine planchon_dem_fill_algorithm
    ! ============================================================================
    !> Subroutine create_geographical_grid
    !! Creates SRTM ASCII files used for flow accumulation calculation.
    !<
    ! ============================================================================
    subroutine create_geographical_grid

        ! Parameters Declaration
        integer :: iunit, ios, k, p
        character(len=128) :: filegeo, stg

        ! Open a file tunnel
        iunit = 10
        filegeo = fdem
        open(iunit,file=filegeo,status="replace",action="write",iostat=ios)
        ! Comply to SRTM ASCII format
        stg = 'ncols'
        write(10,101) trim( stg), strat_X + 2
        stg = 'nrows'
        write(10,101) trim( stg), strat_Y + 2
        stg = 'xllcorner'
        write(10,102) trim( stg), strat_xo - strat_dx
        stg = 'yllcorner'
        write(10,102) trim( stg), strat_yo - strat_dx
        stg = 'cellsize'
        write(10,103) trim( stg), strat_dx
        stg = 'NODATA_value -999999.99'
        write(10,104) trim( stg)

101     format( a5,1x,i3 )
102     format( a9,1x,f12.3 )
103     format( a8,1x,f12.3 )
104     format( a23 )

        ! Write DEM elevation

        ! First line
        write(10,106,advance='no') wdem( 1 )
        do k = 1, strat_X-1
            write(10,106,advance='no') wdem( k )
        enddo
        write(10,107) wdem( k+1 ),wdem( k+1 )

        p = 1
        do k = 1, G_strat%nodes_nb
            if( p < strat_X )then
                if( p == 1 )then
                    write(10,105,advance='no') wdem( k ), wdem( k )
                else
                    write(10,106,advance='no') wdem( k )
                endif
                p = p + 1
            else
                write(10,107) wdem( k ), wdem( k )
                p = 1
            endif
        enddo

        ! Last line
        p = G_strat%nodes_nb - strat_X + 1
        write(10,106,advance='no') wdem( p )
        do k = p, G_strat%nodes_nb - 1
            write(10,106,advance='no') wdem( k )
        enddo
        write(10,107) wdem( G_strat%nodes_nb ),wdem( G_strat%nodes_nb )

105     format( f12.3,1x,f12.3,1x )
106     format( f12.3,1x )
107     format( f12.3,1x,f12.3 )

        close( iunit )

        ! Flow accumulation computation is performed in a CPP file
        call flowacc( fdemcpp, faccucpp )

        return

    end subroutine create_geographical_grid
    ! ============================================================================
    !> Subroutine smooth_function
    !! Creates a smoothing version of the flow elevation.
    !<
    ! ============================================================================
    subroutine smooth_function( input, output )

        integer :: n, m
        real( tkind ), dimension( strat_X, strat_Y ) ::  output, input

        n = strat_X
        m = strat_Y

        output( 2:n-1, 2:m-1 ) = &
            ( input( 1:n-2, 1:m-2 ) + input( 1:n-2, 2:m-1 ) + input( 1:n-2, 3:m ) &
            + input( 2:n-1, 1:m-2 ) + input( 2:n-1, 2:m-1 ) + input( 2:n-1, 3:m ) &
            + input( 3:n, 1:m-2 ) + input( 3:n, 2:m-1 ) + input( 3:n, 3:m ) ) / 9.0_8


        output( (/1,n/), : ) = input( (/1,n/), : )
        output( 2:n-1, (/1,m/) ) = input( 2:n-1, (/1,m/) )

        return

    end subroutine smooth_function
    ! ============================================================================
    !> Subroutine read_flow_accumulaton
    !! Read flow accumulation values.
    !<
    ! ============================================================================
    subroutine read_flow_accumulaton

        ! Parameters Declaration
        integer :: iunit, j, p, ios, k, n
        character(len=128) :: filegeo, stg

        real( tkind ), dimension( strat_X + 2 ) ::  facb
        real( tkind ), dimension( strat_X, strat_Y ) ::  flac, sflac

        iunit = 10
        filegeo = faccu

        ! Check CPP function is finished before starting to read the
        ! flow accumulation file
        do while( p < strat_Y + 8 )
            call count_lines_file( filegeo, p )
        enddo

        ! Flow accumulation function is done, open the created file
        open(iunit,file=filegeo,action="read",iostat=ios)

        ! Read SRTM format
        read(10,*) stg
        read(10,*) stg
        read(10,*) stg
        read(10,*) stg
        read(10,*) stg
        read(10,*) stg

        ! Read flow accumulation values from file.
        ! Be aware there is an additional ghost border in the flow accumulation file
        p = 1
        n = 1
        read( 10,* ) facb( 1 : strat_X + 2 )
        do j = 1, strat_Y
            read( 10,* )facb( 1 : strat_X + 2 )
            facc( n : n + strat_X - 1 ) = facb( 3 : strat_X + 2 )
            flac( 1:strat_X, j ) = facc( n : n + strat_X - 1 )
            n = strat_X + n
        enddo

        ! Done reading close file
        close( iunit )

        ! Smooth the flow accumulation based on local variation
        call smooth_function( flac, sflac  )

        ! Update the flow accumulation values for the global stratal grid
        p = 1
        do j = 1, strat_Y
            do k = 1, strat_X
                facc( p ) = sflac( k, j )
                p = p + 1
            enddo
        enddo

        return

    end subroutine read_flow_accumulaton
    ! ============================================================================
    !> Subroutine surface_parameters_Zevenbergen_Thorne
    !! Zevenbergen and Thorne (1987) proposed a procedure to output surface parameters such as slope angle,
    !! aspect angle, profile curvature and plan curvature through calculating coefficients of the quadratic surface (polynomial)
    !! fitted to the nine neighborhood points. Here we just output the slope and plan curvature.
    !<
    ! ============================================================================
    subroutine surface_parameters_Zevenbergen_Thorne

        integer :: gid, nid( 8 ), p, k

        real( tkind ) :: d, e, f, g, h, zk, zp( 8 ), zmin

        ! Check that arrays for curvature and slope have been initialised
        if( .not. allocated( curv ) ) allocate( curv( L_strat%nodes_nb ) )
        if( .not. allocated( slp ) ) allocate( slp( L_strat%nodes_nb ) )

        do k = 1, L_strat%nodes_nb
            gid = locv_gid( k )
            nid = gv_ngbID( gid, 1:8 )

            ! Top elevation values
            zk = elev_record( gid )

            ! Find the minimum value in the neighborhood
            zmin = 1.e8_8
            do p = 1, 8
                if( nid( p ) <= 0 )then
                    zp( p ) = zk
                else
                    zp( p ) = elev_record( nid( p ) )
                endif
                zmin = min( zp( p ), zmin )
            enddo

            ! Define local parameters
            d = ( 0.5_8 * ( zp( 7 ) + zp( 3 ) ) - zk ) / ( strat_dx**2.0_8 )
            e = ( 0.5_8 * ( zp( 5 ) + zp( 1 ) ) - zk ) / ( strat_dx**2.0_8 )
            f = ( -zp( 6 ) + zp( 4 ) + zp( 8 ) - zp( 2 ) ) / ( 4.0_8 * strat_dx**2.0_8 )
            g = ( zp( 3 ) - zp( 7 ) ) / ( strat_dx * 2.0_8 )
            h = ( zp( 5 ) - zp( 1 ) ) / ( strat_dx * 2.0_8 )

            ! Compute slope angle
            if( g * g + h * h > 0.0_8 )then
                slp( k ) = sqrt( g * g + h * h )
            else
                slp( k ) = 0.0_8
            endif
            if( zk > zmin ) slp( k ) = ( zk - zmin ) / strat_dx

            ! Compute plan curvature
            if( g * g + h * h /= 0.0_8  )then
                curv( k ) = -2.0_8 * ( d * h * h + e * g * g - f * g * h ) / ( g * g + h * h )
            else
                curv( k ) = 0.0_8
            endif

        enddo

        return

    end subroutine surface_parameters_Zevenbergen_Thorne
    ! ============================================================================
    !> Subroutine cmpt_massmovement
    !! Compute mass movement.
    !<
    ! ============================================================================
    subroutine cmpt_massmovement

        integer :: k, p, nb, newlaynb, ks, layID, gid, nid( 8 ), bordernode

        real( tkind ) :: sua, limit, th_ero, zk
        real( tkind ) :: layh, prop, th, temp_h, minelev

        real( tkind ), dimension( totgrn ) :: temp_layh

        ! Get the curvature and slope to find possible areas prone to mass movements
        call surface_parameters_Zevenbergen_Thorne

        do k = 1, L_strat%nodes_nb

            gid = locv_gid( k )

            ! Top elevation
            zk = elev_record( gid )

            ! Based on flow accumulation compute
            ! surface of the upslope area
            sua =  facc( gid ) * ( strat_dx )**2.0_8

            ! From value for stepness and concavity index compute threshold limit
            limit = 1.e6_8
            if( zk >= gsea%actual_sea )then
                limit = mass_mov( 1 )%steepness( 2 ) * sua **( -mass_mov( 1 )%concavity( 2 ) )
                if( sua < mass_mov( 1 )%upslope * 1000.0_8**2_8  ) &
                    limit = mass_mov( 1 )%steepness( 1 ) * sua **( -mass_mov( 1 )%concavity( 1 ) )
            else
                limit = mass_mov( 2 )%steepness( 2 ) * sua **( -mass_mov( 2 )%concavity( 2 ) )
                if( sua < mass_mov( 2 )%upslope *  1000.0_8**2_8 ) &
                    limit = mass_mov( 2 )%steepness( 1 ) * sua **( -mass_mov( 2 )%concavity( 1 ) )
            endif

            ! Check if we are on a border
            bordernode = 0
            ! North border
            if( gv_coord( gid,  2 ) == strat_yo + ( strat_Y -  1 ) * strat_dx )then
                bordernode = 1
            ! South border
            elseif( gv_coord( gid,  2 ) == strat_yo )then
                bordernode = 1
            endif

            ! West border
            if( gv_coord( gid,  1 ) == strat_xo  )then
                bordernode = 1
            ! East border
            elseif( gv_coord( gid,  1 ) == strat_xo + ( strat_X -  1 ) * strat_dx )then
                bordernode = 1
            endif

            ! Mass movement occurs if the gradient is above the threshold for motion
            ! and the curvature is lower than the one defined in the input file.
            if( ( slp( k ) >= limit .and. curv( k ) < mass_curv .and. facc( gid ) > mass_acc ) .or. &
                ( bordernode == 1 .and. slp( k ) >= limit .and. curv( k ) < mass_curv ) )then

                ! Compute mass movement thickness to reach equilibrium slope
                th_ero = strat_dx * ( slp( k ) - limit )

                ! Find the minimum elevation in the neighborhood
                nid = gv_ngbID( gid, 1:8 )
                minelev = 1.e6_8
                do p = 1, 8
                    if( nid( p ) > 0 )then
                        minelev = min( minelev,  elev_record( nid( p ) ) )
                    endif
                enddo

                ! Adjust erosion thickness to prevent formation of holes
                if( th_ero > zk - minelev ) th_ero = zk - minelev
                if( th_ero < 0.0_8 ) goto 10

                ! Get the required thickness from the sedimentary layer
                ! Initialise top layer
                nb = locv_NbLay( k )
                newlaynb = locv_NbLay( k )
                layID = strat_layID( k, nb )

                ! If we are on the basement don't do anything
                if( nb == 1 ) goto 10

                ! If top layer ID is equal to current layer ID don't do anything
                ! the diffusion function will move sediment anyway
                if( layID == L_strat%layersID .and. th_ero <= &
                    strat_thick(  k ,  nb  ) ) goto 10

                ! Initialise top layer sedimentary composition
                temp_h = 0.0_8
                temp_layh( 1:totgrn ) = 0.0_8

                ! If top layer ID is equal to current layer ID take the top layer thickness
                ! and find the erosion thickness which remains and which should be merge
                ! to top layer
                if( layID == L_strat%layersID )then
                    ! Update remaining erosion thickness
                    th_ero = th_ero - strat_thick(  k ,  nb  )
                    ! Copy temporary top layer composition
                    temp_layh( 1:totgrn ) = strat_sedh(  k ,  nb ,  1:totgrn  )
                    temp_h = strat_thick(  k ,  nb  )
                    nb = nb - 1

                ! If top layer is not active we will create one as the diffusion process
                ! only occurs on top layer
                else
                    newlaynb = newlaynb + 1
                    strat_sedh(  k ,  nb + 1,  1:totgrn  ) = 0.0_8
                    strat_thick( k, nb + 1 ) = 0.0_8
                    strat_layID( k, nb + 1 ) = layID
                    strat_porosity(  k ,  nb + 1, 1:totgrn  ) = 0.0_8
                    strat_hardness(  k ,  nb + 1  ) = 1.0_8
                    strat_zelev(  k ,  nb + 1  ) = &
                        strat_zelev(  k ,  nb  )
                endif

                ! Loop over the underlying sedimentary stack
                move_sediment_to_top: do p =  nb, 1, -1
                    if( p == 1 ) exit move_sediment_to_top
                    layh = strat_thick(  k ,  p  )

                    ! All the layer is activated and everything has to be moved to the top layer
                    if( th_ero >= layh )then
                        newlaynb = newlaynb - 1
                        do ks = 1, totgrn
                            temp_layh( ks ) = temp_layh( ks ) + strat_sedh(  k ,  p ,  ks  )
                            temp_h = temp_h + strat_sedh(  k ,  p ,  ks  )
                        enddo
                        th_ero = th_ero - strat_thick(  k ,  p  )
                        strat_thick(  k ,  p  ) = 0.0_8
                        strat_sedh(  k ,  p ,  1:totgrn  ) = 0.0_8
                        strat_porosity(  k ,  p, 1:totgrn  ) = 0.0_8
                        strat_zelev(  k ,  p  ) = strat_zelev(  k ,  p - 1  )

                    ! Just a portion of the layer is activated and needs to be moved to the top layer
                    else
                        th = 0.0_8
                        if( strat_thick(  k ,  p  ) == 0.0_8 )then
                            print*,iam,strat_thick( k, p ), k,p,nb,layh
                            print*,'Something went wrong during mass movement computation.'
                            stop
                        endif

                        ! Proportion of sediment to move to top layer
                        prop = th_ero / strat_thick(  k ,  p  )

                        do ks = 1, totgrn
                            temp_layh( ks ) = temp_layh( ks ) + prop * strat_sedh(  k ,  p ,  ks  )
                            temp_h = temp_h + prop * strat_sedh(  k ,  p ,  ks  )
                            strat_thick(  k ,  p  ) = strat_thick(  k ,  p  ) - prop * strat_sedh(  k ,  p ,  ks  )
                            strat_zelev(  k ,  p  ) = strat_zelev(  k ,  p  ) - prop * strat_sedh(  k ,  p ,  ks  )
                            strat_sedh(  k ,  p ,  ks  ) = ( 1 - prop ) * strat_sedh(  k ,  p ,  ks  )
                            th = th + strat_sedh(  k ,  p ,  ks  )
                        enddo
                        th_ero = 0.0_8

                        exit move_sediment_to_top
                    endif
                enddo move_sediment_to_top

                ! If some underlying layers have been activated update sedimentary stack
                if( temp_h > 0.0_8 )then
                    nb = newlaynb
                    locv_NbLay( k ) = nb
                    strat_layID( k, nb ) = L_strat%layersID
                    strat_thick(  k ,  nb  ) = temp_h
                    th = 0.0_8
                    do ks = 1, totgrn
                        strat_sedh(  k ,  nb ,  ks  ) = temp_layh( ks )
                        th = th  + temp_layh( ks )
                    enddo
                    strat_zelev(  k ,  nb  ) = th + strat_zelev(  k ,  nb - 1  )

                    if( abs( th - temp_h ) > tor )then
                        print*,'Something went wrong when updating layers during mass movement computation.'
                        stop
                    endif
                endif
            endif

10      continue

        enddo

        return

    end subroutine cmpt_massmovement
    ! ============================================================================



end module mod_filldepression
