module RMSim

    use file_data
    use mpidata
    use forces_data
    use strata_data
    use sediment_data

    implicit none

    public


contains

    ! ============================================================================
    !> Subroutine create_geocellular_model
    !! Create the RMS eclipse simulator file.
    !<
    ! ============================================================================
    subroutine create_geocellular_model

       character(len=128) :: file

        integer :: iunit, ios, k, nb1, nb2, i, j, nid, ks, p, m, ks2
        integer :: rms_nx, rms_ny, rms_nblay, ptid( 2 )

        real :: poro, prop, sumprop

        integer, dimension( :,: ), allocatable :: rms_nc
        integer, dimension( :,: ), allocatable :: layid
        integer, dimension( :,:,: ), allocatable :: rms_id

        real, dimension( :,: ), allocatable :: rms_h, rms_poro
        real, dimension( :,:,: ), allocatable :: rms_prop

        ! Create the file
        iunit = 20

        frms = 'geocellular.grdecl'
        file = ''
        file = frms
        call addpath2( file )
        frms = file

        open(iunit,file=frms,status="replace",action="write",iostat=ios)
        rewind(iunit)

        ! SPECGRID defines the number of cells in X, Y, Z directions
        write(iunit,101)'SPECGRID'
101     format( a8 )

        ! First find the closest points coordinates based on user defined extension
        nb1 = int( ( RMSX_min - strat_xo ) / strat_dx )
        nb2 = int( ( RMSX_max - strat_xo ) / strat_dx )
        if( nb1 < strat_X .and. nb2 < strat_X )then
            RMSX_min = strat_xo + strat_dx * nb1
            RMSX_max = strat_xo + strat_dx * nb2
        else
            print*,'Problem in geocellular X extension declaration', nb1, nb2
            stop
        endif
        nb1 = int( ( RMSY_min - strat_yo ) / strat_dx )
        nb2 = int( ( RMSY_max - strat_yo ) / strat_dx )
        if( nb1 < strat_Y .and. nb2 < strat_Y )then
            RMSY_min = strat_yo + strat_dx * nb1
            RMSY_max = strat_yo + strat_dx * nb2
        else
            print*,'Problem in geocellular Y extension declaration', nb1, nb2
            stop
        endif

        ! Define grdecl grid specifications
        rms_nx = int( ( RMSX_max - RMSX_min ) / strat_dx )
        rms_ny = int( ( RMSY_max - RMSY_min ) / strat_dx )
        rms_nblay = L_strat%layersID - 1 - InitDep
        write(iunit,119) rms_nx, rms_ny, rms_nblay, 1, ' F /'
119     format( i8,i8,i8,i3,a4 )
        write(iunit,*)' '

        ! COORD defines the X, Y coordinates of the corner points defining the geo-cells
        ! we use an infinte X-Y line along Z ccordinates
        write(iunit,120)'COORD'
120     format( a5 )
        ! Allocate local id of the points within the stratigraphic mesh
        if( allocated( rms_nc ) ) deallocate( rms_nc )
        allocate( rms_nc( rms_nx + 1, rms_ny + 1 ) )

        rms_nc = -1
        do k = 1, L_strat%nodes_nb
            if( locv_coord( k, 1 ) >= RMSX_min .and.   locv_coord( k, 1 ) <= RMSX_max )then
                if( locv_coord( k, 2 ) >= RMSY_min .and.   locv_coord( k, 2 ) <= RMSY_max )then
                    nb1 = int( ( locv_coord( k, 1 ) - RMSX_min ) / strat_dx ) + 1
                    nb2 = int( ( locv_coord( k, 2 ) - RMSY_min ) / strat_dx ) + 1
                    rms_nc( nb1, nb2 ) = locv_gid( k )
                    write(iunit,*) real( locv_coord( k, 1:2 ) ),0.0, real( locv_coord( k, 1:2 ) ),0.0
                endif
            endif
        enddo
        write(iunit,*) '/ '
        write(iunit,*)' '

        ! Allocate local id of the points for the geo-cellular faces
        if( allocated( rms_id ) ) deallocate( rms_id )
        allocate( rms_id( rms_nx, rms_ny, 4 ) )

        rms_id = -1
        nb2 = 1
        nb1 = 1
        do k = 1, L_strat%faces_nb
            if( locf_XYcoord( k, 2 ) > RMSY_min .and.   locf_XYcoord( k, 2 ) < RMSY_max )then
                if( locf_XYcoord( k, 1 ) > RMSX_min .and.   locf_XYcoord( k, 1 ) < RMSX_max )then
                    rms_id( nb1, nb2, 1:4 ) = locf_points( k, 1:4 )
                    nb1 = nb1 + 1
                    if( locf_XYcoord( k, 1 ) + strat_dx > RMSX_max )then
                        nb1 = 1
                        nb2 = nb2 + 1
                    endif
                endif
            endif
        enddo

        ! Allocate sedimentary layer ids for each points of the geo-cellular model
        if( allocated( layid ) ) deallocate( layid )
        allocate( layid(rms_nx+1, rms_ny+1) )
        layid = 0

        ! Allocate sedimentary layer thickness for each points of the geo-cellular model
        if( allocated( rms_h ) ) deallocate( rms_h )
        allocate( rms_h(L_strat%nodes_nb,L_strat%layersID - 1) )
        rms_h = 1.0e6

        ! Allocate sedimentary layer porosity for each points of the geo-cellular model
        if( allocated( rms_poro ) ) deallocate( rms_poro )
        allocate( rms_poro(L_strat%nodes_nb,L_strat%layersID - 1) )
        rms_poro = 0.0

        ! Allocate sedimentary layer proportions for each points of the geo-cellular model
        if( allocated( rms_prop ) ) deallocate( rms_prop )
        allocate( rms_prop(L_strat%nodes_nb,L_strat%layersID - 1, totgrn) )
        rms_prop = 0.0

        ! Declare values for allocated fields
        do k = 1, L_strat%layersID - 1

            do i = 1, rms_nx + 1
                do j = 1, rms_ny + 1
                    nid = rms_nc( i, j )
                    if( layid( i, j ) < locv_NbLay( nid ) )then
                        if( strat_layID( nid , layid( i, j ) + 1 ) == k )then
                            layid( i, j ) = layid( i, j ) + 1
                        endif
                    endif
                    rms_poro( nid, k ) =  0.0
                    if( strat_layID( nid , layid( i, j ) ) /= k )then
                        rms_h( nid, k ) = rms_h( nid, k - 1 )
                        rms_prop( nid, k, 1:totgrn ) =  0.0
                        rms_poro( nid, k ) = 0.0
                    else
                        rms_h( nid, k ) = real( gsea%actual_sea - strat_zelev( nid, layid( i, j ) ) )
                        if( strat_thick( nid, layid( i, j ) ) == 0.0 )then
                            rms_prop( nid, k, 1:totgrn ) = 0.0
                            rms_poro( nid, k ) = 0.0
                        else
                            do ks = 1, totgrn
                                rms_prop( nid, k, ks ) =  real( strat_sedh( nid, layid( i, j ), ks ) )
                                rms_poro( nid, k ) =  rms_poro( nid, k ) + real( strat_porosity( nid, layid( i, j ), ks ) * &
                                    rms_prop( nid, k, ks )  / strat_thick( nid, layid( i, j ) ) )
                            enddo
                        endif
                    endif
                enddo
             enddo

        enddo

        ! ZCORN defines layers top and bottom faces elevations
        write(iunit,121) 'ZCORN'
121     format( a5 )
        do k = L_strat%layersID - 1, InitDep + 1, -1

            ! Top elevation
            do j = 1, rms_ny
                ! faces south
                do i = 1, rms_nx
                    ptid( 1:2 ) = rms_id( i, j, 1:2 )
                    if( i < rms_nx )then
                        write( iunit, 105, advance='no' )rms_h( ptid( 1 ), k ), rms_h( ptid( 2 ), k )
                   else
                        write( iunit, 107 )rms_h( ptid( 1 ), k ), rms_h( ptid( 2 ), k )
                   endif
                enddo
                ! faces north
                do i = 1, rms_nx
                    ptid( 1:2 ) = rms_id( i, j, 3:4 )
                    if( i < rms_nx )then
                        write( iunit, 105, advance='no' )rms_h( ptid( 1 ), k ), rms_h( ptid( 2 ), k )
                   else
                        write( iunit, 107 )rms_h( ptid( 1 ), k ), rms_h( ptid( 2 ), k )
                   endif
                enddo
            enddo

            ! Bottom elevation
            do j = 1, rms_ny
                ! faces south
                do i = 1, rms_nx
                    ptid( 1:2 ) = rms_id( i, j, 1:2 )
                    if( i < rms_nx )then
                        write( iunit, 105, advance='no' )rms_h( ptid( 1 ), k-1 ), rms_h( ptid( 2 ), k-1 )
                   else
                        write( iunit, 107 )rms_h( ptid( 1 ), k-1 ), rms_h( ptid( 2 ), k-1 )
                   endif
                enddo
                ! faces north
                do i = 1, rms_nx
                    ptid( 1:2 ) = rms_id( i, j, 3:4 )
                    if( i < rms_nx )then
                        write( iunit, 105, advance='no' )rms_h( ptid( 1 ), k-1 ), rms_h( ptid( 2 ), k-1 )
                   else
                        write( iunit, 107 )rms_h( ptid( 1 ), k-1 ), rms_h( ptid( 2 ), k-1 )
                   endif
                enddo

            enddo

        enddo
        write(iunit,*) '/ '

105     format( f12.3,1x,f12.3,1x )
107     format( f12.3,1x,f12.3 )


        ! ACTNUM values
        write(iunit,*)' '
        write(iunit,122) 'ACTNUM'
122     format( a6 )
        do k = L_strat%layersID - 1, InitDep + 1, -1

            do j = 1, rms_ny
                do i = 1, rms_nx
                    p = 0
                    m = 0
                    do ks = 1, 4
                        if( rms_poro( rms_id( i, j, ks ), k ) > 0.0 ) p = p + 1
                    enddo
                    if( p > 0 ) m = 1

                    if( i < rms_nx )then
                        write( iunit, 112, advance='no' )m
                   else
                        write( iunit, 113 )m
                   endif
                enddo
            enddo
        enddo
        write(iunit,*) '/ '


        ! POROSITY values
        write(iunit,*)' '
        write(iunit,123) 'POROSITY'
123     format( a8 )

        do k = L_strat%layersID - 1, InitDep + 1, -1

            do j = 1, rms_ny
                do i = 1, rms_nx
                    poro = 0.0
                    p = 0
                    do ks = 1, 4
                        poro = poro + rms_poro( rms_id( i, j, ks ), k )
                        if( rms_poro( rms_id( i, j, ks ), k ) > 0.0 ) p = p + 1
                    enddo
                    if( p > 0 ) poro = poro / p
                    if( i < rms_nx )then
                        write( iunit, 106, advance='no' )poro
                   else
                        write( iunit, 108 )poro
                   endif
                enddo
            enddo
        enddo
        write(iunit,*) '/ '

        ! PROPORTION values
        do ks = 1, totgrn

            write(iunit,*)' '
            if( ks < 10 ) write(iunit,110) 'PROP',ks
            if( ks > 9 ) write(iunit,111) 'PROP',ks

            do k = L_strat%layersID - 1, InitDep + 1, -1

                do j = 1, rms_ny
                    do i = 1, rms_nx
                        prop = 0.0
                        ! Get total thickness
                        sumprop = 0
                        do ks2 = 1, totgrn
                            do m = 1, 4
                                sumprop = sumprop + rms_prop( rms_id( i, j, m ), k, ks2 )
                            enddo
                        enddo
                        do m = 1, 4
                            prop = prop + rms_prop( rms_id( i, j, m ), k, ks )
                        enddo
                        if( sumprop > 0.0 )then
                            prop = prop / sumprop
                        endif
                        if( i < rms_nx )then
                            write( iunit, 106, advance='no' )prop
                       else
                            write( iunit, 108 )prop
                       endif
                    enddo
                enddo
            enddo
            write(iunit,*) '/ '
        enddo

106     format( f12.3,1x )
108     format( f12.3 )
110     format( a4,i1 )
111     format( a4,i2 )
112     format( i1,1x )
113     format( i1 )

        close( iunit )

        deallocate( rms_id, rms_h, layid, rms_nc, rms_poro, rms_prop )

        return

    end subroutine create_geocellular_model
    ! ============================================================================

end module RMSim
