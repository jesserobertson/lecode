! ============================================================================
! Name        : ClassTIN.f90
! Author      : tristan salles
! Created on: Aug 15, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file ClassTIN.f90
!
! Description : ClassTIN is used to encompass several parameters and types used in TIN definition,
! declaration and evolution.
!
!<
! ============================================================================
module TIN_data

   use precision_data
   use strata_data

   implicit none

   public

   ! Maximum number of neighbors for a given node
   integer, parameter :: mxnghb = 25
   ! Number of nodes on the initial unstructured grid
   integer :: just_nodes
   ! Number of faces on the initial unstructured grid
   integer :: just_faces
   ! Number of TIN faces per partition
   integer :: Tif_Gproc
   ! Number of TIN nodes per partition
   integer, dimension( : ), allocatable :: Tin_Gproc
   ! Global TIN nodes ID in the local partition
   integer, dimension( : ), allocatable :: LN_tin
   ! Global TIN faces ID in the local partition
   integer, dimension( : ), allocatable :: LF_tin

   ! Low & high resolution values for the TIN grid
   real( tkind ) :: grd_LR, grd_HR
   ! Flow accumulation threshold for refinement
   real( tkind ) :: accumulation_refine
   ! Filling step for depression
   real( tkind ) :: step_fill

   ! Flowpath
   integer,dimension(:),allocatable :: G8_pathways
   ! Flow accumulation
   real( tkind ),dimension(:),allocatable :: flow_accu

   !> TIN node type
   ! TIN node coordinates
   real( tkind ),dimension( :,: ),allocatable :: tinv_coord
   ! Boundary of the considered vertex
   integer,dimension( : ),allocatable :: tinv_boundary
   ! Number of connected nodes
   integer,dimension( : ),allocatable :: tinv_ngb
   ! Processor ID containing the node
   integer,dimension( : ),allocatable :: tinv_pid
   ! Number of connected nodes
   integer, dimension( :,: ), allocatable :: ngbID

   !> TIN face type
   ! Is the face on the border (0 if not else >0 if yes)
   integer,dimension( : ),allocatable :: tinf_boundary
   ! Number of neighbors faces
   integer,dimension( : ),allocatable :: tinf_nghb
   ! Processor ID containing the face
   integer,dimension( :,: ),allocatable :: tinf_pid
   ! IDs of the nodes making the face
   integer,dimension( :,: ),allocatable :: tinf_nids
   ! Centroid of the face
   real( tkind ),dimension( :,: ),allocatable :: tinf_centroid
   ! Minimum and maximum lengths of the edges
   real( tkind ),dimension( :,: ),allocatable :: tinf_length
   ! IDs of the neighboring faces
   integer, dimension( :,: ), allocatable :: fids

   !> TIN grid type
   type TIN_grid
      ! Total number of nodes in the TIN
      integer :: nodes_nb
      ! Total number of faces in the TIN
      integer :: faces_nb
      ! Minimal distance of the TIN grid edges
      real( tkind ) :: minimal_dist
      ! Maximal distance of the TIN grid edges
      real( tkind ) :: maximal_dist
   end type TIN_grid
   type( TIN_grid ) :: G_tin

contains

   ! ============================================================================
   !> Function cmp_centroid()
   !! Function cmp_centroid is used to compute the centroid of the face.
   !! \param fce
   !<
   ! ============================================================================
   function cmp_centroid( fce ) result( centroid )

      integer :: fce, k
      integer, dimension( 3 ) :: nids

      real( tkind ), dimension( 3 ) :: centroid

      nids = tinf_nids( fce , 1:3 )
      centroid( 1: 3 ) = 0.0_8
      do k = 1, 3
         centroid( 1 ) = centroid( 1 ) + tinv_coord( nids( k ),  1 )
         centroid( 2 ) = centroid( 2 ) + tinv_coord( nids( k ),  2 )
         centroid( 3 ) = centroid( 3 ) + tinv_coord( nids( k ),  3 )
      enddo
      do k = 1 , 3
         centroid( k ) = centroid( k ) / 3
      enddo

   end function cmp_centroid
   ! ============================================================================
   !> Function cmp_lenght()
   !! Function cmp_lenght is used to compute the min and max  lenght of the edges of a  face.
   !! \param fce
   !<
   ! ============================================================================
   function cmp_lenght( fce ) result( lenght )

      integer :: fce
      integer, dimension( 3 ) :: nids

      real( tkind ) :: mxdist, mndist, l1, l2, l3
      real( tkind ), dimension( 2 ) :: lenght

      nids = tinf_nids( fce, 1:3 )
      mxdist = 0.0_8
      mndist = 1000.0_8 * strat_dx
      l1 = sqrt( ( tinv_coord( nids( 2 ),  1 ) - tinv_coord( nids( 1 ),  1 ) )**2 + &
      ( tinv_coord( nids( 2 ),  2 ) - tinv_coord( nids( 1 ),  2 ) )**2 )
      mndist = min( mndist, l1 )
      mxdist = max( mxdist, l1 )
      l2 = sqrt( ( tinv_coord( nids( 3 ),  1 ) - tinv_coord( nids( 1 ),  1 ) )**2 + &
      ( tinv_coord( nids( 3 ),  2 ) - tinv_coord( nids( 1 ),  2 ) )**2 )
      mndist = min( mndist, l2 )
      mxdist = max( mxdist, l2 )
      l3 = sqrt( ( tinv_coord( nids( 2 ),  1 ) - tinv_coord( nids( 3 ),  1 ) )**2 + &
      ( tinv_coord( nids( 2 ),  2 ) - tinv_coord( nids( 3 ),  2 ) )**2 )
      mndist = min( mndist, l3 )
      mxdist = max( mxdist, l3 )
      lenght( 1 ) = mndist
      lenght( 2 ) = mxdist

   end function cmp_lenght
  ! ============================================================================

end module TIN_data
! ============================================================================

