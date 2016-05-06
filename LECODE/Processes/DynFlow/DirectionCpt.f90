! ============================================================================
! Name        : DirectionCpt.f90
! Author      : tristan salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file DirectionCpt.f90
!!
!! D8_method utilises an adaptation of the the Global D8 method proposed by Paik Kyungrock
!! to determine the flow walkers flowpaths  (Journal of Geophysical Reseach vol 113, 2008).
!<
! ============================================================================
module mod_direction

   use mpidata
   use TIN_data
   use time_data
   use mod_sort
   use strata_ini
   use error_data
   use fwalker_data
   use strata_data
   use param_data
   use forces_data

   implicit none

   public

contains

   ! ============================================================================
   !> Subroutine find_global_pathways_G8_Paik
   !! Subroutine find_global_pathways_G8_Paik utilises an adaptation of the Global D8
   !! method proposed by Paik Kyungrock to determine the flow walkers flowpaths  (JGR vol 113, 2008).
   !<
   ! ============================================================================
   subroutine find_global_pathways_G8_Paik

      integer :: i, k, it, i0, i1, j
      integer :: ilast
      integer, dimension( : ), allocatable :: dir

      real( tkind ) :: maxel
      real( tkind ) :: dist
      real( tkind ) :: elast, gradmax, grad( mxnghb )

      allocate( dir( G_tin%nodes_nb ) )

      ! Allocate direction path
      if( allocated( G8_pathways ) ) deallocate( G8_pathways )
      allocate( G8_pathways( G_tin%nodes_nb ) )

      ! Mark boundaries as termination points
      dir = -10
      do

         ! Set the initial value to very deep
         maxel = -1.0e8_8

         ! Scan through and find the highest remaining point on the surface
         do j = 1, Tin_Gproc( iam + 1 )
            i = LN_tin( j )
            if( tinv_coord( i,  3 ) >= maxel &
            .and. dir( i ) == -10 )then
               it = i
               maxel = real( tinv_coord( i,  3 ) )
            endif
         enddo

         ! Maxel should now contain the highest remaining point or
         ! if no points are left, then it should contain the initial value
         if(maxel == -1.0e8_8) exit
         ilast = it
         elast = maxel
         i1 = -1
         loop: do
            ! From the last position, locate the path of steepest decent
            gradmax = -1.0e16_8
            grad = -1.0e16_8
            it = ilast
            do k = 1, tinv_ngb( it ) 
               i0 = ngbID( it, k )
               if( i0 < 1 .or. i0 > G_tin%nodes_nb )print*,'Something went wrong in Paik algorithm.'
               grad( k ) = elast - tinv_coord( i0, 3 )
               dist = compute_distance( tinv_coord( ilast, 1:3 ), tinv_coord( i0, 1:3 ) )
               if( dist > 0 )then
                  grad( k ) = dble( grad( k ) / dist )
               else
                  grad( k ) = 0.0_8
               endif
               if( gradmax <= grad( k ) )then
                  if( i1 > 0 .and. i0 > 0 )then
                     gradmax = grad( k )
                     i1 = i0
                  elseif( i1 < 0 .and. i0 < 0 )then
                     gradmax = grad( k )
                     i1 = i0
                  elseif( i1 < 0 .and. i0 > 0 )then
                     gradmax = grad( k )
                     i1 = i0
                  elseif( i1 > 0 .and. i0 < 0 )then
                     if( gradmax == grad( k ) )then
                        gradmax = grad( k )
                        i1 = i1
                     else
                        gradmax = grad( k )
                        i1 = i0
                     endif
                  endif
               endif
            enddo
            if( gradmax <= 0.0_8 )then
               ! We are in a hollow, or on a flat plane there is nowhere to go.
               ! Set a default value of -9
               dir( ilast ) = -9
               exit loop
            else
               dir( ilast ) = i1
               if( i1 < 0 )then
                  dir( ilast ) = -9
                  exit loop
               endif
               if( dir( i1 ) /= -10 ) &
               exit loop
               if( tinv_pid( i1 ) /= iam )&
               exit loop
               ilast = i1
               elast = tinv_coord( i1, 3 )
            endif
         enddo loop
      enddo

      call mpi_allreduce( dir, G8_pathways, G_tin%nodes_nb, int_type, &
      max_type, SPModel_comm_world, ierr )

      deallocate( dir )

      return

   end subroutine find_global_pathways_G8_Paik
   ! ============================================================================
   !> Subroutine find_slope_within_TINface
   !! Subroutine find_slope_within_TINface computes the average slope for the considered
   !! flow walker based on the local plane for the cloud of points defined by the flow
   !! walker position as well as the 3 neighbors of the closest point of the flow walker.
   !! \param facenb, fwnb, sx, sy
   !<
   ! ============================================================================
   subroutine find_slope_within_TINface( facenb, fwnb, sx, sy )

      integer :: facenb, fwnb, fid( 3 )
      real( tkind ) :: x, y, x1, y1, z1, x2, y2, z2, x3, y3, z3
      real( tkind ) :: A, B, C, D, sx, sy

      ! Get flow walker point
      x = real( fa_xpos( fwnb ) )
      y = real( fa_ypos( fwnb ) )

      ! Get face points
      fid = tinf_nids( facenb, 1:3 )
      x1 = real( tinv_coord( fid( 1 ),  1 ) )
      y1 = real( tinv_coord( fid( 1 ),  2 ) )
      z1 = real( tinv_coord( fid( 1 ),  3 ) )
      x2 = real( tinv_coord( fid( 2 ),  1 ) )
      y2 = real( tinv_coord( fid( 2 ),  2 ) )
      z2 = real( tinv_coord( fid( 2 ),  3 ) )
      x3 = real( tinv_coord( fid( 3 ),  1 ) )
      y3 = real( tinv_coord( fid( 3 ),  2 ) )
      z3 = real( tinv_coord( fid( 3 ),  3 ) )

      if( fa_zpos( fwnb ) >  gsea%actual_sea .and. fa_type( fwnb ) < tpturb )then
         z1 = max( z1, gsea%actual_sea )
         z2 = max( z2, gsea%actual_sea )
         z3 = max( z3, gsea%actual_sea )
      endif

      A = y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2)
      B = z1*(x2 - x3) + z2*(x3 - x1) + z3*(x1 - x2)
      C = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)
      D = x1*(y2*z3 - y3*z2) + x2*(y3*z1 - y1*z3) + x3*(y1*z2 - y2*z1)

      ! Plane equation: A*x + B*y + C*z = D
      if( C /= 0.0_8 )then
          sx = -A/C
          sy = -B/C
          sx = dble( int( sx*1.0e5_8 ) * 1.0e-5_8 )
          sy = dble( int( sy*1.0e5_8 ) * 1.0e-5_8 )
      else
          sx = 0.0_8
          sy = 0.0_8
      endif

      return

   end subroutine find_slope_within_TINface
   ! ============================================================================
   !> Subroutine find_elevation_within_TINface
   !! Subroutine find_elevation_within_TINface finds the z-coordinate of the flow walker
   !! inside a triangle based on the equation of the plane using the 3 points
   !! which defined a triangular faced of the deposit grid.
   !! \param facenb, fwnb, zfe
   !<
   ! ============================================================================
   subroutine find_elevation_within_TINface( facenb, fwnb, z )

      integer :: facenb, fwnb, fid( 3 )
      real( tkind ) :: x, y, z, x1, y1, z1, x2, y2, z2, x3, y3, z3
      real( tkind ) :: A, B, C, D

      ! Get flow walker point
      x = real( fa_xpos( fwnb ) )
      y = real( fa_ypos( fwnb ) )

      ! Get face points
      fid = tinf_nids( facenb, 1:3 )
      x1 = real( tinv_coord( fid( 1 ),  1 ) )
      y1 = real( tinv_coord( fid( 1 ),  2 ) )
      z1 = real( tinv_coord( fid( 1 ),  3 ) )
      x2 = real( tinv_coord( fid( 2 ),  1 ) )
      y2 = real( tinv_coord( fid( 2 ),  2 ) )
      z2 = real( tinv_coord( fid( 2 ),  3 ) )
      x3 = real( tinv_coord( fid( 3 ),  1 ) )
      y3 = real( tinv_coord( fid( 3 ),  2 ) )
      z3 = real( tinv_coord( fid( 3 ),  3 ) )

      if( x == x1 .and. y == y1 )then
         z = z1
         return
      elseif( x == x2 .and. y == y2 )then
         z = z2
         return
      elseif(  x == x3 .and. y == y3 )then
         z = z3
         return
      endif
      A = y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2)
      B = z1*(x2 - x3) + z2*(x3 - x1) + z3*(x1 - x2)
      C = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)
      D = x1*(y2*z3 - y3*z2) + x2*(y3*z1 - y1*z3) + x3*(y1*z2 - y2*z1)

      ! Plane equation: A*x + B*y + C*z = D
      z = D - A*x - B*y
      if( C /= 0.0_8 )then
         z = dble( z / C )
      else
         attempt = CBASEVAL
         print*,'face nb',facenb,G_tin%faces_nb
         print*,'flow element position',x,y
         print*,'TIN face pt1',fid( 1 ),tinv_coord( fid( 1 ),1:3)
         print*,'TIN face pt2',fid( 2 ),tinv_coord( fid( 2 ),1:3)
         print*,'TIN face pt3',fid( 3 ),tinv_coord( fid( 3 ),1:3)
      endif

      return

   end subroutine find_elevation_within_TINface
   ! ============================================================================
   !> Subroutine record_flow_walker_faceID
   !! Subroutine record_flow_walker_faceID finds the flow walker direction.
   !! \param fc, fw_id
   !<
   ! ============================================================================
   subroutine record_flow_walker_faceID( fc, fw_id )

      integer :: k, fw_id, fc

      ! Move face ID from one increment to free the first position
      do k = numold-1, 1, -1
        fa_faceID( fw_id , k + 1 ) = fa_faceID( fw_id , k )
      enddo

      ! Add the new face ID at first position
      fa_faceID( fw_id , 1 ) = fc
      if( fa_faceID( fw_id , 1 ) > G_tin%faces_nb )then
        print*,'Something went wrong when recording FW face ID.'
        stop
      endif

      return

   end subroutine record_flow_walker_faceID
  ! ============================================================================

end module mod_direction
