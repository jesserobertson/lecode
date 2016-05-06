! ============================================================================
! Name        : Sort_algo.f90
! Author      : tristan salles
! Created on: Aug 17, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file Sort_algo.f90
!!
!<
! ============================================================================
module mod_sort

   use fwalker_data
   use TIN_function

   implicit none

   public

   integer, parameter :: Isw = 10

   type Limits
      integer :: Ileft, Iright
   end type Limits

contains

   ! ============================================================================
   !> subroutine find_intersection_between_segments()
   !!  return false if no solution exists.
   !! \param pt1, pt2, pt3, pt4, solution, intersect
   !<
   ! ============================================================================
   subroutine find_intersection_between_segments( pt1, pt2, pt3, pt4, solution, intersect )

      ! Variables declaration
      logical :: solution
      real( tkind ) :: pt1( 2 ), pt2( 2 ), pt3( 2 ), pt4( 2 )
      real( tkind ) :: p13( 2 ), p43( 2 ), p21( 2 ), intersect( 2 )
      real( tkind ) :: num_a, num_b , denom, mua, mub
      integer :: sol, headingpt

      solution = .false.
      headingpt = 0
      intersect( : ) = -99999999.99_8
      p21( 1 ) = pt2( 1 ) - pt1( 1 )
      p21( 2 ) = pt2( 2 ) - pt1( 2 )
      p43( 1 ) = pt4( 1 ) - pt3( 1 )
      p43( 2 ) = pt4( 2 ) - pt3( 2 )
      p13( 1 ) = pt1( 1 ) - pt3( 1 )
      p13( 2 ) = pt1( 2 ) - pt3( 2 )
      denom = p43( 2 ) * p21( 1 ) - p43( 1 ) * p21( 2 )
      num_a = p43( 1 ) * p13( 2 ) - p43( 2 ) * p13( 1 )
      num_b = p21( 1 ) * p13( 2 ) - p21( 2 ) * p13( 1 )
      if ( denom == 0.0_8 )then
         ! Coincident lines
         if( num_a == 0.0_8 .and. num_b == 0.0_8 )then
            sol = 1
         else
            ! Parallel lines
            sol = 2
         endif
      else
         mua = num_a / denom
         mub = num_b / denom
         ! Intersecting
         if( mua >= 0.0_8 .and. mua <= 1.0_8 .and. mub >= 0.0_8 .and. mub <= 1.0_8 )then
            intersect( 1 ) = pt1( 1 ) + mua * ( pt2( 1 ) - pt1( 1 ) )
            intersect( 2 ) = pt1( 2 ) + mua * ( pt2( 2 ) - pt1( 2 ) )
            sol = 3
            ! Not intersecting
         else
            sol = 4
         endif
      endif
      if( sol == 3 ) solution = .true.

      return

   end subroutine find_intersection_between_segments
   ! ============================================================================
   !> subroutine random_pickup
   !! subroutine random_pickup picks up random indexes from a list of indexes.
   !! \param nb, list_in, height_in, list_out, height_out, maxnb
   !<
   ! ============================================================================
   subroutine random_pickup( nb, list_in, height_in, list_out, height_out, maxnb  )

      ! Variables declaration
      integer :: nb, i, picked, k, maxnb

      integer,dimension( nb ):: list_inter, height_inter
      integer,dimension( nb ),intent( in ) :: list_in, height_in
      integer,dimension( maxnb ),intent( out ) :: list_out, height_out

      list_out = -1
      list_inter = -1
      height_out = -1
      height_inter = -1
      do i = 1, nb
         picked = pickone( nb )
         list_inter( i ) = list_in( picked )
         height_inter( i ) = height_in( picked )
      enddo
      k = 0
      loop: do i = 1, nb
         picked = pickone( nb )
         k = k + 1
         list_out( k ) = list_inter( picked )
         height_out( k ) = height_inter( picked )
         if( k == maxnb ) exit loop
      enddo loop

      return

   end subroutine random_pickup
   ! ============================================================================
   !> Function pickone()
   !! Function pickone is used to return a random integer from a list
   !! \param nb
   !<
   ! ============================================================================
   function pickone( nb ) result( randint )

      ! Variables declaration
      integer, intent( in ) :: nb
      integer :: randint, seed

      real( tkind ) :: rand

      call random_seed(size=seed)
      call random_number( harvest=rand )
      randint = 1 + int( real( nb ) * rand )
      if( randint > nb ) randint = nb

   end function pickone
   ! ============================================================================
   !> subroutine reallocate_flow_walker_space
   !! subroutine reallocate_flow_walker_space change the maximum number of flow walker based on the simulation requirements.
   !! \param fwmax
   !<
   ! ============================================================================
   subroutine reallocate_flow_walker_space( fwmax )

       ! parameters declaration
       integer :: fwmax, maxT, oldmaxfw, ks
       ! FW global reference
       integer,dimension(:),allocatable :: fa_refid1
       ! Number of point in the area of influence for ero/dep rules
       integer,dimension(:),allocatable :: fa_nbPtAI1
       ! FW inside the simulation domain(0), if outside (1), if forced deposition (2)
       integer,dimension(:),allocatable :: fa_inbound1
       ! FW flow type
       integer,dimension(:),allocatable :: fa_type1
       ! FW TIN step counter
       integer,dimension(:),allocatable :: fa_count1
       ! FW x position
       real( tkind ),dimension(:),allocatable :: fa_xpos1
       ! FW y position
       real( tkind ),dimension(:),allocatable :: fa_ypos1
       ! FW z position
       real( tkind ),dimension(:),allocatable :: fa_zpos1
       ! X slope experienced by FW
       real( tkind ),dimension(:),allocatable :: fa_xslp1
       ! Y slope experienced by FW
       real( tkind ),dimension(:),allocatable :: fa_yslp1
       ! FW height
       real( tkind ),dimension(:),allocatable :: fa_h1
       ! Stream width
       real( tkind ),dimension(:),allocatable :: fa_width1
       ! FW x velocity
       real( tkind ),dimension(:),allocatable :: fa_xvel1
       ! FW y velocity
       real( tkind ),dimension(:),allocatable :: fa_yvel1
       ! FW x velocity solution for the ODE
       real( tkind ),dimension(:),allocatable :: fa_xvelO1
       ! FW y velocity solution for the ODE
       real( tkind ),dimension(:),allocatable :: fa_yvelO1
       ! FW flow volume
       real( tkind ),dimension(:),allocatable :: fa_volume1
       ! FW flow discharge (m3/s)
       real( tkind ),dimension(:),allocatable :: fa_discharge1
       ! FW stream density
       real( tkind ),dimension(:),allocatable :: fa_density1
       ! FW previous gradient value
       real( tkind ),dimension(:),allocatable :: fa_pgrad1
       ! C1 coefficient of bottom friction
       real( tkind ),dimension(:),allocatable :: fa_Cfric1
       ! Flow walker acceleration along X direction
       real( tkind ),dimension(:),allocatable :: fa_accx1
       ! Flow walker acceleration along Y direction
       real( tkind ),dimension(:),allocatable :: fa_accy1
       ! Zone of flow establishment for the plume
       real( tkind ),dimension(:),allocatable :: fa_pzfe1
       ! Plume coefficient of decceleration
       real( tkind ),dimension(:),allocatable :: fa_pcoeff1
       ! Plume trajectory distance from the mouth
       real( tkind ),dimension(:),allocatable :: fa_pdist1
       ! Plume initial x velocity at the mouth
       real( tkind ),dimension(:),allocatable :: fa_pvelx1
       ! Plume initial y velocity at the mouth
       real( tkind ),dimension(:),allocatable :: fa_pvely1
       ! FW all previous face ID on the TIN grid
       integer,dimension(:,:),allocatable :: fa_faceID1
       ! FW GIDs point in the area of influence for ero/dep rules
       integer,dimension(:,:),allocatable :: fa_ptsAIfw1
       ! FW sediment discharge (m3/s)
       real( tkind ),dimension(:,:),allocatable :: fa_sedcharge1

       print*,'flow walkers reallocation',iam,maxfw,fwmax
       call mpi_barrier( SPModel_comm_world, ierr )

      oldmaxfw = maxfw
      do while( maxfw < fwmax )
         maxfw = maxfw * 2
      enddo

      if( maxfw > oldmaxfw )then
          call mpi_allreduce( maxfw, maxT, 1, int_type, max_type, SPModel_comm_world, ierr )
          maxfw = maxT
          allocate( fa_refid1( num_fw ) )
          allocate( fa_nbPtAI1( num_fw ) )
          allocate( fa_inbound1( num_fw ) )
          allocate( fa_type1( num_fw ) )
          allocate( fa_count1( num_fw ) )
          allocate( fa_xpos1( num_fw ) )
          allocate( fa_ypos1( num_fw ) )
          allocate( fa_zpos1( num_fw ) )
          allocate( fa_xslp1( num_fw ) )
          allocate( fa_yslp1( num_fw ) )
          allocate( fa_width1( num_fw ) )
          allocate( fa_h1( num_fw ) )
          allocate( fa_xvel1( num_fw ) )
          allocate( fa_yvel1( num_fw ) )
          allocate( fa_xvelO1( num_fw ) )
          allocate( fa_yvelO1( num_fw ) )
          allocate( fa_volume1( num_fw ) )
          allocate( fa_discharge1( num_fw ) )
          allocate( fa_density1( num_fw ) )
          allocate( fa_pgrad1( num_fw ) )
          allocate( fa_Cfric1( num_fw ) )
          allocate( fa_accx1( num_fw ) )
          allocate( fa_accy1( num_fw ) )
          allocate( fa_pzfe1( num_fw ) )
          allocate( fa_pcoeff1( num_fw ) )
          allocate( fa_pdist1( num_fw ) )
          allocate( fa_pvelx1( num_fw ) )
          allocate( fa_pvely1( num_fw ) )
          allocate( fa_faceID1( num_fw, numold ) )
          allocate( fa_ptsAIfw1( num_fw, maxPtsFW ) )
          allocate( fa_sedcharge1( num_fw, totgrn ) )

          fa_refid1 = fa_refid
          fa_nbPtAI1 = fa_nbPtAI
          fa_inbound1 = fa_inbound
          fa_type1 = fa_type
          fa_width1 = fa_width
          fa_count1 = fa_count
          fa_xpos1 = fa_xpos
          fa_ypos1 = fa_ypos
          fa_zpos1 = fa_zpos
          fa_xslp1 = fa_xslp
          fa_yslp1 = fa_yslp
          fa_h1 = fa_h
          fa_xvel1 = fa_xvel
          fa_yvel1 = fa_yvel
          fa_xvelO1 = fa_xvelO
          fa_yvelO1 = fa_yvelO
          fa_volume1 = fa_volume
          fa_discharge1 = fa_discharge
          fa_density1 = fa_density
          fa_pgrad1 = fa_pgrad
          fa_Cfric1 = fa_Cfric
          fa_accx1 = fa_accx
          fa_accy1 = fa_accy
          fa_pzfe1 = fa_pzfe
          fa_pcoeff1 = fa_pcoeff
          fa_pdist1 = fa_pdist
          fa_pvelx1 = fa_pvelx
          fa_pvely1 = fa_pvely
          fa_faceID1 = fa_faceID
          fa_ptsAIfw1 = fa_ptsAIfw
          fa_sedcharge1 = fa_sedcharge
          deallocate( fa_refid, fa_nbPtAI, fa_inbound, fa_type, fa_count, fa_width )
          deallocate( fa_xpos, fa_ypos, fa_zpos, fa_xslp, fa_yslp, fa_h, fa_xvel )
          deallocate( fa_yvel, fa_xvelO, fa_yvelO, fa_volume, fa_discharge, fa_density )
          deallocate( fa_pgrad, fa_Cfric, fa_accx, fa_accy, fa_pzfe, fa_pcoeff, fa_pdist )
          deallocate( fa_pvelx, fa_pvely, fa_faceID, fa_ptsAIfw, fa_sedcharge )

          allocate( fa_refid( maxfw ) )
          allocate( fa_nbPtAI( maxfw ) )
          allocate( fa_width( maxfw ) )
          allocate( fa_inbound( maxfw ) )
          allocate( fa_type( maxfw ) )
          allocate( fa_count( maxfw ) )
          allocate( fa_xpos( maxfw ) )
          allocate( fa_ypos( maxfw ) )
          allocate( fa_zpos( maxfw ) )
          allocate( fa_xslp( maxfw ) )
          allocate( fa_yslp( maxfw ) )
          allocate( fa_h( maxfw ) )
          allocate( fa_xvel( maxfw ) )
          allocate( fa_yvel( maxfw ) )
          allocate( fa_xvelO( maxfw ) )
          allocate( fa_yvelO( maxfw ) )
          allocate( fa_volume( maxfw ) )
          allocate( fa_discharge( maxfw ) )
          allocate( fa_density( maxfw ) )
          allocate( fa_pgrad( maxfw ) )
          allocate( fa_Cfric( maxfw ) )
          allocate( fa_accx( maxfw ) )
          allocate( fa_accy( maxfw ) )
          allocate( fa_pzfe( maxfw ) )
          allocate( fa_pcoeff( maxfw ) )
          allocate( fa_pdist( maxfw ) )
          allocate( fa_pvelx( maxfw ) )
          allocate( fa_pvely( maxfw ) )
          allocate( fa_faceID( maxfw, numold ) )
          allocate( fa_ptsAIfw( maxfw, maxPtsFW ) )
          allocate( fa_sedcharge( maxfw, totgrn ) )

          ! Reallocate erosion deposition arrays
          if( allocated( ed_refid ) )then
            deallocate( ed_refid, ed_procid, ed_bound, ed_type )
            deallocate( ed_AIpts, ed_AInb, ed_zpos, ed_h, ed_vel )
            deallocate( ed_vol, ed_dens, ed_sed )
          endif
          allocate( ed_refid( maxfw ) )
          allocate( ed_procid( maxfw ) )
          allocate( ed_bound( maxfw ) )
          allocate( ed_type( maxfw ) )
          allocate( ed_AIpts( maxfw, maxPtsFW ) )
          allocate( ed_AInb( maxfw ) )
          allocate( ed_zpos( maxfw ) )
          allocate( ed_h( maxfw ) )
          allocate( ed_vel( maxfw ) )
          allocate( ed_vol( maxfw ) )
          allocate( ed_dens( maxfw ) )
          allocate( ed_sed( maxfw, 10 ) )

          fa_inbound = 1
          fa_faceID = 0
          fa_sedcharge = 0.0_8

          fa_refid( 1:num_fw ) = fa_refid1( 1:num_fw )
          fa_nbPtAI( 1:num_fw ) = fa_nbPtAI1( 1:num_fw )
          fa_width( 1:num_fw ) = fa_width1( 1:num_fw )
          fa_inbound( 1:num_fw ) = fa_inbound1( 1:num_fw )
          fa_type( 1:num_fw ) = fa_type1( 1:num_fw )
          fa_count( 1:num_fw ) = fa_count1( 1:num_fw )
          fa_xpos( 1:num_fw ) = fa_xpos1( 1:num_fw )
          fa_ypos( 1:num_fw ) = fa_ypos1( 1:num_fw )
          fa_zpos( 1:num_fw ) = fa_zpos1( 1:num_fw )
          fa_xslp( 1:num_fw ) = fa_xslp1( 1:num_fw )
          fa_yslp( 1:num_fw ) = fa_yslp1( 1:num_fw )
          fa_h( 1:num_fw ) = fa_h1( 1:num_fw )
          fa_xvel( 1:num_fw ) = fa_xvel1( 1:num_fw )
          fa_yvel( 1:num_fw ) = fa_yvel1( 1:num_fw )
          fa_xvelO( 1:num_fw ) = fa_xvelO1( 1:num_fw )
          fa_yvelO( 1:num_fw ) = fa_yvelO1( 1:num_fw )
          fa_volume( 1:num_fw ) = fa_volume1( 1:num_fw )
          fa_discharge( 1:num_fw ) = fa_discharge1( 1:num_fw )
          fa_density( 1:num_fw ) = fa_density1( 1:num_fw )
          fa_pgrad( 1:num_fw ) = fa_pgrad1( 1:num_fw )
          fa_Cfric( 1:num_fw ) = fa_Cfric1( 1:num_fw )
          fa_accx( 1:num_fw ) = fa_accx1( 1:num_fw )
          fa_accy( 1:num_fw ) = fa_accy1( 1:num_fw )
          fa_pzfe( 1:num_fw ) = fa_pzfe1( 1:num_fw )
          fa_pcoeff( 1:num_fw ) = fa_pcoeff1( 1:num_fw )
          fa_pdist( 1:num_fw ) = fa_pdist1( 1:num_fw )
          fa_pvelx( 1:num_fw ) = fa_pvelx1( 1:num_fw )
          fa_pvely( 1:num_fw ) = fa_pvely1( 1:num_fw )
          do ks = 1, numold
            fa_faceID( 1:num_fw, ks ) = fa_faceID1( 1:num_fw, ks )
          enddo
          do ks = 1, maxPtsFW
            fa_ptsAIfw( 1:num_fw, ks ) = fa_ptsAIfw1( 1:num_fw, ks )
          enddo
          do ks = 1, totgrn
            fa_sedcharge( 1:num_fw, ks ) = fa_sedcharge1( 1:num_fw, ks )
          enddo
          deallocate( fa_refid1, fa_nbPtAI1, fa_inbound1, fa_type1, fa_count1, fa_width1 )
          deallocate( fa_xpos1, fa_ypos1, fa_zpos1, fa_xslp1, fa_yslp1, fa_h1, fa_xvel1 )
          deallocate( fa_yvel1, fa_xvelO1, fa_yvelO1, fa_volume1, fa_discharge1, fa_density1 )
          deallocate( fa_pgrad1, fa_Cfric1, fa_accx1, fa_accy1, fa_pzfe1, fa_pcoeff1, fa_pdist1 )
          deallocate( fa_pvelx1, fa_pvely1, fa_faceID1, fa_ptsAIfw1, fa_sedcharge1 )
      endif

      return

   end subroutine reallocate_flow_walker_space
   ! ============================================================================

end module mod_sort
! ============================================================================
! Module quicksort arranges array elements from smallest to largest
!
! grabbed from A millers web site http://users.bigpond.net.au/amiller/
! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.
! pjr added module declaration
! mvr modified integer array to intent inout - may now be any integer
!     array that gets sorted along with associated real array
! ============================================================================
module quicksort


   use precision_data

   implicit none

   save

   private

   public quick_sort

contains

  ! ============================================================================
  recursive subroutine quick_sort(list, order)

    real( tkind ), dimension(:), intent( inout ) :: list
    integer, dimension(:), intent( inout ) :: order

    call quick_sort_1(1, size(list))

    contains

        ! ============================================================================
        recursive subroutine quick_sort_1(left_end, right_end)

            integer, intent( in ) :: left_end, right_end

            ! Local variables
            integer :: i, j, itemp
            real( tkind ) :: reference, temp
            integer, parameter :: max_simple_sort_size = 6

            if (right_end < left_end + max_simple_sort_size) then
                ! Use interchange sort for small lists
                call interchange_sort(left_end, right_end)
            else
                ! Use partition ("quick") sort
                reference = list((left_end + right_end)/2)
                i = left_end - 1; j = right_end + 1

                do
                    ! Scan list from left end until element >= reference is found
                    do
                        i = i + 1
                        if (list(i) >= reference) exit
                    end do
                    ! Scan list from right end until element <= reference is found
                    do
                        j = j - 1
                        if (list(j) <= reference) exit
                    end do

                    if (i < j) then
                        ! Swap two out-of-order elements
                        temp = list(i); list(i) = list(j); list(j) = temp
                        itemp = order(i); order(i) = order(j); order(j) = itemp
                    else if (i == j) then
                        i = i + 1
                        exit
                    else
                        exit
                    end if
            end do

            if (left_end < j) call quick_sort_1(left_end, j)
            if (i < right_end) call quick_sort_1(i, right_end)
        end if

      end subroutine quick_sort_1
      ! ============================================================================
      subroutine interchange_sort(left_end, right_end)

        integer, intent( in ) :: left_end, right_end

        !  Local variables
        integer :: i, j, itemp
        real( tkind ) :: temp

        do i = left_end, right_end - 1
            do j = i+1, right_end
                if (list(i) > list(j)) then
                    temp = list(i); list(i) = list(j); list(j) = temp
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                end if
            end do
         end do

     end subroutine interchange_sort
     ! ============================================================================

  end subroutine quick_sort
  ! ============================================================================

end module quicksort


