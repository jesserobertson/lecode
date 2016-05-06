module UDWplug2

    use file_data
    use mpidata
    use strata_data
    use sediment_data

    implicit none

    public

contains
    subroutine luke

    integer :: p, k, i, ks
    integer, dimension( : ),allocatable :: layid
    real, dimension(:), allocatable :: x

    allocate( layid(L_strat%nodes_nb) )

    ! For each variables
    allocate( x( L_strat%nodes_nb * L_strat%layersID ) )
    y
    z
    allocate( sedh( L_strat%nodes_nb * L_strat%layersID, totgrn ) )

    p = 0
    layid = 0

    ! For each stratigraphic layer
    do k = 1, L_strat%layersID - 1
        ! For each processor
        do i = 1, L_strat%nodes_nb
            p = p + 1
            gid = locv_gid( i )
            x( p ) = real( locv_coord( i, 1 ) )
            y( p ) = real( locv_coord( i, 2 ) )
            ! Number of layers for the considered point
            ! The local layer nb of the point
            if( layid( i ) < locv_NbLay( i ) )then
                if( strat_layID( i , layid( i ) + 1 ) == k )&
                    layid( i ) = layid( i ) + 1
            endif

            ! If the layer is absent take the elebation of the bottom point
            if( strat_layID( i , layid( i ) ) /= k )then
                z( p ) = z( p - L_strat%nodes_nb )
            ! Else if present just records the elevation for the layer
            else
                z( p ) = real( strat_zelev( i, layid( i ) ) )
            endif
            sedh( p, 1:totgrn ) = real( strat_sedh( i, layid( i ), 1:totgrn ) )
        enddo

    enddo

    end subroutine luke
end module UDWplug2
