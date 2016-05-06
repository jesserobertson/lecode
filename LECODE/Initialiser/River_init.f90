! ============================================================================
! Name        : River_init.f90
! Author      : tristan salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
!
! ============================================================================
!> \file River_init.f90
!!
!! River_init computes the inflow into the system according to sources and slopes.
!!
!<
! ============================================================================
module mod_inflow

    use file_data
    use mpidata
    use flux_data
    use time_data
    use mod_sort
    use error_data
    use fwalker_data
    use strata_data
    use strata_ini
    use param_data
    use forces_data
    use TIN_function
    use sediment_data

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine define_river_fws
    !! Subroutine define_river_fws introduces rivers by allocating new flow walkers into
    !! the simulation area.
    !<
    ! ============================================================================
    subroutine define_river_fws

        ! Variables declaration
        logical :: marine

        integer :: i, j, n, k, i0, ks, p, newsrc, km, seed, nids( 3 )
        integer :: srcproc( nproc ), idsrc_proc( nproc, 2 ), reloc, relocT, idS

        real( tkind ) :: slsed, randomVal, ranges, fv, newwidth

        integer,dimension(:),allocatable :: srcid, fce, fcS, rcv
        real( tkind ),dimension( : ),allocatable :: listfwx, listfwy

        ! Record currrent number of flow walkers
        k = num_fw

        ! Define local parameters
        newsrc = 0
        srcproc = 0
        idsrc_proc = 0
        sourceactive = 0

        ! Find the number of river sources that needs to be define for the
        ! current time step
        do i = 1, num_src
            if( tnow >= fws_tstrt( i ) .and. tnow < fws_tend( i ) - time_tolerance )then
                newsrc = newsrc + 1
            endif
        enddo

        ! In case there is no rivers exit the subroutine
        if( newsrc == 0 ) return

        ! Allocate local arrays
        allocate( srcid( newsrc ), fce( newsrc ), fcS( newsrc ) )
        allocate( listfwx( newsrc ), listfwy( newsrc ), rcv( newsrc ) )

        ! In first processor define the position of the new rivers
        if( iam == 0 )then
            j = 1

            ! Loop over the river sources
            do i = 1, num_src

                ! Find the ones active for the current time step
                if( tnow >= fws_tstrt( i ) .and. tnow < fws_tend( i ) - time_tolerance )then

                    srcid( j ) = i
                    ! Move randomly the river position based on defined range
                    call random_seed(size=seed)
                    call random_number(harvest=randomVal)
                    ranges = 0.0_8
                    listfwx( j ) = real( fws_xposition( i ) + fws_xrange( i ) * ( 2.0_8 * randomVal - 1.0_8 ) )
                    listfwy( j ) = real( fws_yposition( i ) + fws_yrange( i ) * ( 2.0_8 * randomVal - 1.0_8 ) )
                    if( fws_xrange( i ) == 0 .and. fws_yrange( i ) == 0 )then
                        ranges = 0.001_8
                        listfwx( j ) = real( fws_xposition( i ) + ranges * ( 2.0_8 * randomVal - 1.0_8 ) )
                        listfwy( j ) = real( fws_yposition( i ) + ranges * ( 2.0_8 * randomVal - 1.0_8 ) )
                    endif
                    j = j + 1
                endif

            enddo
        endif

        ! Broadcast river position values to all processors
        call mpi_bcast( srcid, newsrc, int_type, 0, SPModel_comm_world,ierr )
        call mpi_bcast( listfwx( 1:newsrc ), newsrc, dbl_type, 0, SPModel_comm_world,ierr )
        call mpi_bcast( listfwy( 1:newsrc ), newsrc, dbl_type, 0, SPModel_comm_world,ierr )

        ! Find the TIN face containing the river flow walker
        call global_find_TINfaces_kdtree( newsrc, listfwx, listfwy, 5*mxnghb, fce )

        do ks = 1, newsrc
            if( fce( ks ) < 0 ) print*,'Something went wrong when locating river flow walkers.'
            if( fce( ks ) < 0 ) stop
        enddo

        ! Find the position of the flow walker within the stratal grid
        srcproc = 0
        rcv = -1
        do ks = 1, newsrc
            call closest_stratal_node( listfwx( ks ), listfwy( ks ), fcS( ks ), idS )
            lpp:do n = 1, 2
                i0 = gv_ShareID( idS, n )
                if( i0 >= 0 )then
                    srcproc( i0 + 1 ) = srcproc( i0 + 1 ) + 1
                    if( iam == i0 ) rcv( srcproc( i0 + 1 ) ) = ks
                    exit lpp
                endif
            enddo lpp
        enddo

        ! Check that the number of flow walkers is not
        ! greater than the maximum allowed
        reloc = 0
        relocT = 0
        k = num_fw + srcproc( iam + 1 )
        if( k > maxfw ) reloc = 1
        call mpi_allreduce( reloc, relocT, 1, int_type, max_type, SPModel_comm_world, ierr )
        if( relocT == 1 ) call reallocate_flow_walker_space( k )

        ! Register river streams as active flow walkers
        k = num_fw
        km = 0

        define_rivers_fw: do i = 1, srcproc( iam + 1 )
            i0 = rcv( i )

            ! Check that the rivers are active for the current time step
            if( tnow >= fws_tstrt( srcid( i0 ) ) .and. tnow < fws_tend( srcid( i0 ) ) - time_tolerance )then

                sourceactive = 1

                do j =  1, fws_num( srcid( i0 ) )
                    nids = tinf_nids( fce( i0 ), 1:3 )
                    marine = .false.
                    ! In case the source has been defined as a river and is sitting below sea level
                    ! do not create this flow walker
                    do p = 1, 3
                        if( gsea%actual_sea > tinv_coord( nids( p ), 3 ) .and. fws_type( srcid( i0 ) ) == 0 ) &
                            marine = .true.
                    enddo
                    if( .not. marine )then
                        k = k + 1
                        km = km + 1
                        if( k > maxfw )then
                            attempt = IN_MAXFW
                            exit define_rivers_fw
                        else

                            fa_refid( k ) = -1
                            fa_faceID( k , : ) = 0

                            ! Main parameters defining flow walker
                            fa_xpos( k ) = listfwx( i0 )
                            fa_ypos( k ) = listfwy( i0 )

                            ! Add a random component to the velocities
                            call random_seed(size=seed)
                            call random_number(harvest=randomVal)
                            fa_xvel( k ) = dble( fws_xvel( srcid( i0 ) ) * ( 1.0_8 + 0.5_8 * randomVal ) )
                            call random_seed(size=seed)
                            call random_number(harvest=randomVal)
                            fa_yvel( k ) = dble( fws_yvel( srcid( i0 ) ) * ( 1.0_8 + 0.5_8 * randomVal ) )
                            fv = sqrt( fa_xvel( k )**2 + fa_yvel( k )**2 )
                            fa_volume( k ) = dble( fws_volume( srcid( i0 ) ) )
                            fa_h( k ) = dble( fws_srch( srcid( i0 ) ) )

                            ! Update width
                            newwidth = dble( fws_qfl( srcid( i0 ) ) / ( fv * fws_srch( srcid( i0 ) ) ) )
                            newwidth = min( newwidth, transport%wdtmax )
                            newwidth = max( newwidth, transport%wdtmin )
                            fa_width( k ) = newwidth

                            ! Define discharge
                            fa_discharge( k ) = dble( fws_qfl( srcid( i0 ) ) )

                            ! Define type of flow
                            if( fws_type( srcid( i0 ) ) == 1 )then
                                fa_type( k ) = tpturb
                            else
                                fa_type( k ) = 0
                            endif

                            ! Define sediment properties and density
                            slsed = 0.0_8
                            fa_sedcharge( k , : ) = 0.0_8
                            do ks = 1, silgrn
                                fa_sedcharge( k , ks ) = dble( fws_sedcharge( srcid( i0 ), ks ) )
                                slsed = slsed + fa_sedcharge( k , ks ) * sediment( ks )%density
                            enddo
                            fa_density( k ) = dble( slsed / fa_volume( k ) + fluid_density )

                            fa_count( k ) = 0
                            fa_pgrad( k ) = 0.0_8
                            fa_faceID( k , 1 ) = fce( i )
                            fa_inbound( k ) = 0

                        endif
                    endif
                enddo
            endif
        enddo define_rivers_fw

        ! Ensure that all processes have managed to declare properly the river flow walkers.
        call completion

        ! Deallocate local arrays
        deallocate( srcid, fce, listfwx, listfwy, fcS, rcv )

        ! Update number of flow walkers
        num_fw = k

        ! Transfer value globally
        call mpi_allreduce( num_fw, tot_fw, 1, int_type, sum_type, SPModel_comm_world, ierr )

        return

    end subroutine define_river_fws
    ! ============================================================================

end module mod_inflow
