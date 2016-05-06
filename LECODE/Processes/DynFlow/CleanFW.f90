! ============================================================================
! Name        : CleanFW.f90
! Author      : tristan salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file CleanFW.f90
!!
!! CleanFW checks for flow walkers that are out of bounds and replaces memory
!! used to store properties of these flow walkers by flow walkers that are
!! still within the area, thereby freeing memory for additional flow walkers.
!!
!<
! ============================================================================
module mod_cleanfw

    use mpidata
    use mod_sort
    use fwalker_data
    use param_data
    use TIN_function

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine clean_flow_walkers
    !! Subroutine clean_flow_walkers checks for flow walkers that are out of bounds and replaces
    !! memory used to store properties of these flow walkers by flow walkers that are
    !! still within the area, thereby freeing memory for additional flow walkers.
    !<
    ! ============================================================================
    subroutine clean_flow_walkers

        implicit none

        integer :: k, ks
        logical :: lastvalid

        ! Adjust flow walker mass balance
        k = 1
        do while( k <= num_fw )
            lastvalid = .false.
            do
                if( lastvalid ) exit
                ! As long as we have some flow walkers
                if( num_fw > 0 ) then
                    ! Is the flow wlaker out of the simulation
                    if( fa_inbound( num_fw ) == 1 )then
                        num_fw = num_fw - 1
                    else
                        lastvalid = .true.
                    endif
                else
                    lastvalid = .true.
                endif
            enddo

            if( num_fw > k )then

                ! Move flow walker to the free position
                if( fa_inbound( k ) == 1 )then

                    ! Copy active flow walker in the free position
                    fa_refid( k ) = fa_refid( num_fw )
                    fa_nbPtAI( k ) = fa_nbPtAI( num_fw )
                    fa_inbound( k ) = fa_inbound( num_fw )
                    fa_type( k ) = fa_type( num_fw )
                    fa_count( k ) = fa_count( num_fw )
                    fa_xpos( k ) = fa_xpos( num_fw )
                    fa_ypos( k ) = fa_ypos( num_fw )
                    fa_zpos( k ) = fa_zpos( num_fw )
                    fa_xslp( k ) = fa_xslp( num_fw )
                    fa_yslp( k ) = fa_yslp( num_fw )
                    fa_h( k ) = fa_h( num_fw )
                    fa_xvel( k ) = fa_xvel( num_fw )
                    fa_yvel( k ) = fa_yvel( num_fw )
                    fa_xvelO( k ) = fa_xvelO( num_fw )
                    fa_yvelO( k ) = fa_yvelO( num_fw )
                    fa_volume( k ) = fa_volume( num_fw )
                    fa_discharge( k ) = fa_discharge( num_fw )
                    fa_density( k ) = fa_density( num_fw )
                    fa_pgrad( k ) = fa_pgrad( num_fw )
                    fa_Cfric( k ) = fa_Cfric( num_fw )
                    fa_accx( k ) = fa_accx( num_fw )
                    fa_accy( k ) = fa_accy( num_fw )
                    fa_pzfe( k ) = fa_pzfe( num_fw )
                    fa_pcoeff( k ) = fa_pcoeff( num_fw )
                    fa_pdist( k ) = fa_pdist( num_fw )
                    fa_pvelx( k ) = fa_pvelx( num_fw )
                    fa_pvely( k ) = fa_pvely( num_fw )
                    do ks = 1, numold
                        fa_faceID( k, ks ) = fa_faceID( num_fw, ks )
                    enddo
                    if( fa_faceID( k, 1 ) > G_tin%faces_nb )then
                        print*,"Something went wrong in cleaning face ID:",&
                            k,fa_faceID( k , 1 ),G_tin%faces_nb
                        stop
                    endif
                    do ks = 1, maxPtsFW
                        fa_ptsAIfw( k, ks ) = fa_ptsAIfw( num_fw, ks )
                    enddo
                    do ks = 1, totgrn
                        fa_sedcharge( k, ks ) = fa_sedcharge( num_fw, ks )
                    enddo
                    ! Update number of flow walker in the current processor
                    num_fw = num_fw - 1
                endif
            endif
            k = k + 1
        enddo

        return

    end subroutine clean_flow_walkers
    ! ============================================================================
    !> Subroutine update_flow_walker_TINposition
    !! Subroutine update_flow_walker_TINposition reset flow walkers values after a display.
    !<
    ! ============================================================================
    subroutine update_flow_walker_TINposition

        integer :: ks
        integer,dimension( num_fw ):: fce
        integer, dimension( nproc ) :: numel

        do ks = 1, num_fw
            fa_faceID( ks , : ) = 0
            fa_ptsAIfw( ks, : ) = 0
        enddo

        numel = 0
        call mpi_allgather(num_fw,1,int_type,numel,1,int_type,SPModel_comm_world,ierr)

        call find_TINfaces_kdtree( tot_fw, fa_xpos, fa_ypos, mxnghb, fce )
        do ks = 1, num_fw
            if( fce( ks  ) > 0 ) fa_inbound( ks ) = 0
            if( fce( ks ) < 0 ) fa_inbound( ks ) = 1
            fa_faceID( ks , 1 ) = fce( ks )
            if( fce( ks ) < 0 ) print*,'Something went wrong when updating FW position in TIN grid'
            if( fce( ks ) < 0 ) stop
            if( fce( ks ) > G_tin%faces_nb ) &
                print*,'Something went wrong when updating FW position in TIN grid...'
            if( fce( ks ) > G_tin%faces_nb ) stop
        enddo

        return

    end subroutine update_flow_walker_TINposition
    !===================================================

end module mod_cleanfw
