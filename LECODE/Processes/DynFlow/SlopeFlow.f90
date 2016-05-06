! ============================================================================
! Name        : SlopeFlow.f90
! Author      : tristan salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file SlopeFlow.f90
!!
!! slope_flow computes the slope for each flow walker.
!!
!! The file contains the following subroutine:
!! \arg \c travel_distance_vector
!! \arg \c compute_slope_along_travel_direction
!<
! ============================================================================
module mod_slopefw

    use TIN_out
    use mpidata
    use mod_sort
    use strata_ini
    use fwalker_data
    use param_data
    use forces_data
    use kdtree2_module
    use kdtree2_precision_module

    implicit none

    public

    private :: travel_distance_vector

contains

    ! ============================================================================
    !> Subroutine travel_distance_vector
    !! Subroutine travel_distance_vector determines which of the 3 surrounding face points
    !! is in the direction of movement.
    !! \param fw_id, fc, pt_id
    !<
    ! ============================================================================
    subroutine travel_distance_vector( fw_id, fc, pt_id )

        logical :: intersect, in

        type(kdtree2), pointer :: tree
        type(kdtree2_result), dimension( 1 ) :: results

        integer :: k, p1, p2, p3, fw_id, fc, nnb, pt_id( 3 ), rate, fid( 3 )

        real( tkind ) :: grid_val( 2, 4 ), newpt( 2 ), pt1( 2 ), pt2( 2 ), pt3( 2 ), fv, multi
        real( tkind ) :: s(3), lside, ptintersect( 2 ), d1, d2, d0, ratio1, ratio2

        ! Get face's points coordinates
        fid = tinf_nids( fc, 1:3 )
        grid_val( 1, 1 ) = real( tinv_coord( fid( 1 ),  1 ) )
        grid_val( 2, 1 ) = real( tinv_coord( fid( 1 ),  2 ) )
        grid_val( 1, 2 ) = real( tinv_coord( fid( 2 ),  1 ) )
        grid_val( 2, 2 ) = real( tinv_coord( fid( 2 ),  2 ) )
        grid_val( 1, 3 ) = real( tinv_coord( fid( 3 ),  1 ) )
        grid_val( 2, 3 ) = real( tinv_coord( fid( 3 ),  2 ) )

        ! Add the flow walker coordinates
        grid_val( 1, 4 ) = real( fa_xpos( fw_id ) )
        grid_val( 2, 4 ) = real( fa_ypos( fw_id ) )

        ptintersect( : ) = -9999999999.99_8
        nnb = 1
        rate = 1

        ! If flow walker velocity not null
        fv = sqrt( ( fa_xvel( fw_id ) )**2 + ( fa_yvel( fw_id ) )**2 )

        if( fv > 0.0_8 )then

            ! Get the maximum side
            lside = 0.0_8
            s( 1 ) = sqrt( ( grid_val( 1, 1 ) - grid_val( 1, 2 ) )**2 + ( grid_val( 2, 1 ) - grid_val( 2, 2 ) )**2 )
            lside = max( lside, s( 1 ) )
            s( 2 ) = sqrt( ( grid_val( 1, 3 ) - grid_val( 1, 1 ) )**2 + ( grid_val( 2, 3 ) - grid_val( 2, 1 ) )**2 )
            lside = max( lside, s( 2 ) )
            s( 3 ) = sqrt( ( grid_val( 1, 2 ) - grid_val( 1, 3 ) )**2 + ( grid_val( 2, 2 ) - grid_val( 2, 3 ) )**2 )
            lside = max( lside, s( 3 ) )
            lside = lside * 1000.0_8
            multi = dble( 100.0_8 * lside / fv )

20          continue

            ! Based on flow walker velocity and position defined
            ! point in the direction of travel
            newpt( 1 ) = grid_val( 1, 4 ) + ( multi * fa_xvel( fw_id ) ) * rate
            newpt( 2 ) = grid_val( 2, 4 ) + ( multi * fa_yvel( fw_id ) ) * rate
            pt3( 1 ) = grid_val( 1, 4 )
            pt3( 2 ) = grid_val( 2, 4 )
            in = .false.

            intersectionSegments: do k = 1, 3
                ! Allocate points
                p1= k
                pt1( 1 ) = grid_val( 1, k ) * rate
                pt1( 2 ) = grid_val( 2, k ) * rate
                if( k == 3 )then
                    p2 = k - 2
                    pt2( 1 ) = grid_val( 1, k - 2 ) * rate
                    pt2( 2 ) = grid_val( 2, k - 2 ) * rate
                else
                    p2 = k + 1
                    pt2( 1 ) = grid_val( 1, k + 1 ) * rate
                    pt2( 2 ) = grid_val( 2, k + 1 ) * rate
                endif
                if( k == 1 ) p3 = 3
                if( k == 2 ) p3 = 1
                if( k == 3 ) p3 = 2
                ! Check if segments intersect
                call find_intersection_between_segments( pt1, pt2, pt3, newpt, intersect, ptintersect )

                if( intersect ) then
                    in = .true.
                    d0 = sqrt( ( pt2( 1 ) - pt1( 1 ) )**2 + ( pt2( 2 ) - pt1( 2 ) )**2 )
                    d1 = sqrt( ( ptintersect( 1 ) - pt1( 1 ) )**2 + ( ptintersect( 2 ) - pt1( 2 ) )**2 )
                    d2 = sqrt( ( ptintersect( 1 ) - pt2( 1 ) )**2 + ( ptintersect( 2 ) - pt2( 2 ) )**2 )
                    ratio1 = d1 / d0
                    ratio2 = d2 / d0
                    if( ratio1 <= ratio2 )then
                        pt_id( 1 ) = fid( p1 )
                        pt_id( 2 ) = fid( p2 )
                        pt_id( 3 ) = fid( p3 )
                    else
                        pt_id( 1 ) = fid( p2 )
                        pt_id( 2 ) = fid( p1 )
                        pt_id( 3 ) = fid( p3 )
                    endif
                    exit intersectionSegments
                endif
            enddo intersectionSegments

            if( .not. in )then
                if( rate == 1 )then
                    rate = 100000
                    goto 20
                elseif( rate == 100000 )then
                    goto 21
                endif
            else
                goto 22
            endif
        endif

21      continue

        ! Create a kd-tree
        tree => kdtree2_create( grid_val, sort = .true., rearrange = .true. )
        call kdtree2_n_nearest_around_point(tree, idxin=4, nn=nnb, correltime=1, results=results)
        call kdtree2_destroy(tree)
        if( results( 1 )%idx == 1 )then
            pt_id( 1 ) = fid( 1 )
            pt_id( 2 ) = fid( 2 )
            pt_id( 3 ) = fid( 3 )
        elseif( results( 1 )%idx == 2 )then
            pt_id( 1 ) = fid( 2 )
            pt_id( 2 ) = fid( 1 )
            pt_id( 3 ) = fid( 3 )
        else
            pt_id( 1 ) = fid( 3 )
            pt_id( 2 ) = fid( 2 )
            pt_id( 3 ) = fid( 1 )
        endif

22      continue

        return

   end subroutine travel_distance_vector
   ! ============================================================================
   !> Subroutine compute_slope_along_travel_direction
   !! Subroutine compute_slope_along_travel_direction calculates the angle and direction.
   !! \param fw_nb, face_nb, xsb, ysb
   !<
   ! ============================================================================
   subroutine compute_slope_along_travel_direction( fw_nb, face_nb, xsb, ysb )

       integer :: fw_nb, face_nb, loc_pt( 3 ), i, kpt, kpt2, k
       real( tkind ) :: xsb, ysb, xsa, ysa, delts, delts2, maxtravel,diff, deltx, delty

       xsb = fa_xslp( fw_nb )
       ysb = fa_yslp( fw_nb )

       ! Look at the direction of travel and determines which of the 3 surrounding face points
       ! is in the direction of movement
       call travel_distance_vector( fw_nb, face_nb, loc_pt )

       kpt = -1
       kpt2 = -1
       maxtravel = strat_dx
       i = loc_pt( 1 )

       if( i > 0 )then
           kpt = i
           ! Compute maximal travel distance
           do k = 1, tinv_ngb( kpt )
               kpt2 = ngbID( kpt, k )
               delts = compute_distance( tinv_coord(kpt, 1:3 ), tinv_coord(kpt2, 1:3 ) )
               maxtravel = min( maxtravel, delts )
           enddo
           maxtravel = max( strat_dx, maxtravel )
       endif

       ! Choose the second pointing direction to compute slope angle
       if( G8_pathways( kpt ) <= 0 )then
           xsb = 0.0_8
           ysb = 0.0_8
           return
       endif

       if( G8_pathways( kpt ) > -9 )then
           if( G8_pathways( kpt ) > 0 )then
               kpt2 = G8_pathways( kpt )
               if( G8_pathways( kpt2 ) > 0 )then
                   if( kpt2 == loc_pt( 1 ) .or. kpt2 == loc_pt( 2 ) .or. kpt2 == loc_pt( 3 ) )then
                       if( G8_pathways( kpt2 ) > 0 ) kpt2 = G8_pathways( kpt2 )
                   endif
               else
                    kpt2 = kpt
               endif
               xsa = xsb
               ysa = ysb

               ! On a plain keep local slope as default
               if( fa_zpos( fw_nb ) == tinv_coord( kpt,  3 ) .and. &
                   tinv_coord( kpt2,  3 ) == fa_zpos( fw_nb ) ) return

               ! At the interface with ocean use first pointing direction
               if( fa_zpos( fw_nb ) >=  gsea%actual_sea .and. &
                   tinv_coord( kpt,  3 ) < gsea%actual_sea )then
                   deltx = dble( fa_xpos( fw_nb ) - tinv_coord( kpt,  1 ) )
                   delty = dble( fa_ypos( fw_nb ) - tinv_coord( kpt,  2 ) )
                   delts2 = dble( deltx**2 + delty**2 )
                   if( delts2 /= 0.0_8 )then
                       diff = int(( fa_zpos( fw_nb ) + fa_h( fw_nb ) )*1.0e4_8 - &
                           gsea%actual_sea *1.0e4_8 ) * 1.0e-4_8
                       xsa = dble( deltx * diff / delts2 )
                       ysa = dble( delty * diff / delts2 )
                   endif
                   xsb = dble( xsa )
                   ysb = dble( ysa )

               ! Near the ocean interface on the second direction
               elseif( fa_zpos( fw_nb ) >=  gsea%actual_sea .and. &
                   tinv_coord( kpt2,  3 ) < gsea%actual_sea )then
                   deltx = dble( fa_xpos( fw_nb ) - tinv_coord( kpt2,  1 ) )
                   delty = dble( fa_ypos( fw_nb ) - tinv_coord( kpt2,  2 ) )
                   delts2 = dble( deltx**2 + delty**2 )
                   if( delts2 /= 0.0_8 )then
                       diff = int(( fa_zpos( fw_nb ) + fa_h( fw_nb ) )*1.0e4_8 - &
                           gsea%actual_sea * 1.0e4_8 ) * 1.0e-4_8
                       xsa = dble( deltx * diff / delts2 )
                       ysa = dble( delty * diff / delts2 )
                   endif
                   xsb = dble( xsa )
                   ysb = dble( ysa )

               else
                   deltx = dble( fa_xpos( fw_nb ) - tinv_coord( kpt2,  1 ) )
                   delty = dble( fa_ypos( fw_nb ) - tinv_coord( kpt2,  2 ) )
                   delts2 = dble( deltx**2 + delty**2 )
                   if( delts2 /= 0.0_8 )then
                       diff = int( ( fa_zpos( fw_nb ) + fa_h( fw_nb ) )*1.0e4_8 - &
                           tinv_coord( kpt2,  3 ) *1.0e4_8 ) * 1.0e-4_8
                       xsa = dble( deltx * diff / delts2 )
                       ysa = dble( delty * diff / delts2 )
                   endif
                   xsb = dble( xsa )
                   ysb = dble( ysa )
               endif
           endif

       endif

       return

   end subroutine compute_slope_along_travel_direction
   ! ============================================================================

end module mod_slopefw
