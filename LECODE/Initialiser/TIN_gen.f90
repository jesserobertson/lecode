! ============================================================================
! Name        : TIN_gen.f90
! Author      : tristan salles
! Created on: Aug 16, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file TIN_gen.f90
!
! Description : TIN_gen is used to create the visualisation triangular surface used in SPModel simulation.
!
!<
! ============================================================================
module TIN_surface

    use file_data
    use TIN_data
    use error_data
    use TIN_function
    use strata_data
    use mod_filldepression

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine generate_TIN_surface
    !! Subroutine generate_TIN_surface generates the TIN surface.
    !<
    ! ============================================================================
    subroutine generate_TIN_surface

        ! Parameters Declaration
        character(len=128) :: fnde

        ! Call Triangle to built the TIN based on cloud points
        if( iam == 0 )then
            fnde = ''
            fnde = fnode
            call noblnk( fnde )
            fnde( len( fnde ) : len( fnde ) ) = CHAR(0)
            call trianglegen( 0, fnde )
        endif
        call mpi_barrier(SPModel_comm_world,ierr)

        ! Allocate the TIN grid
        call build_TIN_grid_allocation

        return

    end subroutine generate_TIN_surface
    ! ============================================================================
    !> Subroutine build_TIN_grid_allocation
    !! Subroutine build_TIN_grid_allocation allocates TIN parameters: nodes coordinates,
    !! faces and grid.
    !<
    ! ============================================================================
    subroutine build_TIN_grid_allocation

        ! Parameters Declaration
        logical :: found

        integer :: iunit, ios, n, kn, ks, kw, ke, nid( 3 ), fid, id( 3 ), totid, p1
        integer :: km, kp, nb, dimension, skip, skip1, skip2, k, maxneigh, nb1

        character(len=128) :: tElem

        ! Name of Triangle output files
        iunit = 444
        tElem = 'TIN.1.ele'
        call addpath2( tElem )

        if( allocated( LN_tin ) )deallocate( LN_tin )
        allocate( LN_tin( Tin_Gproc( iam + 1 ) ) )

        LN_tin = 0
        k = 0
        do n = 1, G_tin%nodes_nb
            if( tinv_pid( n ) == iam )then
                k = k + 1
                LN_tin( k ) = n
            endif
        enddo

        ! Read elements
        inquire(FILE=tElem, EXIST=found)
        if(found)then
            open(iunit,file=tElem,status="old",action="read",iostat=ios)
            rewind(iunit)
        else
            attempt = TIN_FILE
        endif

        ! Ensure that all processes have successfully read the element file.
        call completion
        read(iunit, *) G_tin%faces_nb, dimension, skip

        just_faces = 0

        G_tin%minimal_dist = 1000.0_8 * strat_dx
        G_tin%maximal_dist = 0.0_8

        ! Read the TIN file face values
        do n = 1, G_tin%faces_nb
            read(iunit, *,end=123) nb, skip, skip1, skip2
            tinf_nids( n, 1 ) = skip
            tinf_nids( n, 2 ) = skip1
            tinf_nids( n, 3 ) = skip2
            tinf_pid( n, 1 ) =  tinv_pid( skip )
            tinf_pid( n, 2 ) =  tinv_pid( skip1 )
            tinf_pid( n, 3 ) =  tinv_pid( skip2 )
            if( tinf_pid( n, 1 ) == tinf_pid( n, 2 ) ) tinf_pid( n, 2 ) = -1
            if( tinf_pid( n, 1 ) == tinf_pid( n, 3 ) ) tinf_pid( n, 3 ) = -1
            if( tinf_pid( n, 2 ) == tinf_pid( n, 3 ) ) tinf_pid( n, 3 ) = -1
            tinf_boundary( n ) = 0
            tinf_nghb( n ) = 0
            tinf_centroid( n, 1:3 ) = cmp_centroid( n )
            tinf_length( n, 1:2 ) = cmp_lenght( n )
            G_tin%minimal_dist = min( G_tin%minimal_dist, tinf_length( n, 1 ) )
            G_tin%maximal_dist = max( G_tin%maximal_dist, tinf_length( n, 2 ) )
            totid = 0
            do k = 1, 3
                kn = tinf_nids( n, k )
                if( tinf_pid( n, k ) == iam ) Tif_Gproc = Tif_Gproc + 1
                if( tinv_boundary( kn )  >= 1 ) totid = totid + 1
            enddo
            if( totid >= 1 ) tinf_boundary( n ) = totid
            If( tinf_boundary( n ) == 0 ) just_faces = just_faces + 1
        enddo
123     continue

        ! Close element file
        close( iunit )

        if( allocated( LF_tin ) )deallocate( LF_tin )
        allocate( LF_tin( Tif_Gproc ) )
        LF_tin = 0
        k = 0
        do n = 1, G_tin%faces_nb
            do kp = 1, 3
                if( tinf_pid( n, kp ) == iam )then
                    k = k + 1
                    LF_tin( k ) = n
                endif
            enddo
        enddo

        ! Built face neighbors
        maxneigh = mxnghb
        fids = -1
        tinf_nghb = 0
        if( Tif_Gproc < maxneigh ) maxneigh =  Tif_Gproc - 1
        call search_tin_face_neighbors_kdtree( maxneigh )

        ! Built TIN node neighbors
        do n = 1, G_tin%faces_nb
            nb1 = tinf_nghb( n )
            nid(1:3) = tinf_nids( n, 1:3)
            do k = 1, 3
                if( tinv_ngb( nid( k ) ) /= 0 ) goto 21
                if( k == 1 )then
                    ngbID( nid( k ), 1 ) = nid( 2 )
                    ngbID( nid( k ), 2 ) = nid( 3 )
                elseif( k == 2 )then
                    ngbID( nid( k ), 1 ) = nid( 1 )
                    ngbID( nid( k ), 2 ) = nid( 3 )
                else
                    ngbID( nid( k ), 1 ) = nid( 1 )
                    ngbID( nid( k ), 2 ) = nid( 2 )
                endif
                kp = 2

                do p1 = 1, nb1
                    fid = fids( n, p1 )
                    found = .false.
                    if( fid < 1 .or. fid > G_tin%faces_nb )then
                        print*,'Something went wrong when looking for TIN faces.'
                        stop
                    endif
                    id = tinf_nids( fid, 1:3 )
                    do km = 1, 3
                        if( id( km ) == nid( k ) ) found = .true.
                    enddo

                    if( found )then
                        do km = 1, 3
                            if( id( km ) /= nid( k ) ) kp = kp + 1
                            if( id( km ) /= nid( k ) ) ngbID( nid( k ), kp ) = id( km )
                        enddo
                        ke = kp
                        do km = 1, kp - 1
                            do ks = km + 1, kp
                                if( ngbID( nid( k ), ks ) == ngbID( nid( k ), km ) )then
                                    ngbID( nid( k ), ks ) = -1
                                    ke = ke - 1
                                endif
                            enddo
                        enddo
                        do km = 1, kp - 1
                            if( ngbID( nid( k ), km ) == -1 )then
                                do kw = km + 1, kp
                                    ngbID( nid( k ), kw - 1 ) = ngbID( nid( k ), kw )
                                enddo
                            endif
                        enddo
                        if( ke < kp )then
                            do km = ke + 1, kp
                                ngbID( nid( k ), km ) = -1
                            enddo
                        endif
                        kp = ke
                    endif
                enddo
                tinv_ngb( nid( k ) ) = kp
21              continue
            enddo
        enddo

        return

    end subroutine build_TIN_grid_allocation
    ! ============================================================================
    !> Subroutine generate_TIN_grid
    !! Subroutine generate_TIN_grid update a TIN using the refinement values defined by the user.
    !<
    ! ============================================================================
    subroutine generate_TIN_grid

        if( iam == 0 )then
            ! Planchon algorithm
            call planchon_dem_fill_algorithm
            ! Create DEM geographical grid for flow accumulation calculation
            call create_geographical_grid
            ! Read computed flow accumulation
            call read_flow_accumulaton
        endif

        ! Broadcast flow accumulation and DEM elevation
        call mpi_bcast( wdem,G_strat%nodes_nb,dbl_type,0,SPModel_comm_world,ierr )
        call mpi_bcast( facc,G_strat%nodes_nb,dbl_type,0,SPModel_comm_world,ierr )

        ! Define the grid based on flow accumulation
        call create_adaptive_TIN_grid

        ! Generate the triangular grid
        call generate_TIN_surface

        return

    end subroutine generate_TIN_grid
    ! ============================================================================

end module TIN_surface

