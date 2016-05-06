! ============================================================================
! Name        : TINI.f90
! Author      : tristan salles
! Created on: Aug 16, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file TINI.f90
!
! Description : TINI is used to gather the information within the SPModel XmL input file to build
! the TIN used in SPModel.
!
!<
! ============================================================================
module TIN_parser

   use file_data
   use FoX_sax
   use TIN_data
   use strata_data
   use error_data
   use precision_data
   use FoX_common

   implicit none

   public

   ! Refinement Parameters
   logical, save :: refinesection = .false.
   logical, save :: in_RMesh = .false.
   logical, save :: in_HR = .false.
   logical, save :: in_LR = .false.
   logical, save :: in_SR = .false.
   logical, save :: in_SG = .false.
   logical, save :: in_SF = .false.
   character(len=128), save :: lr,hr,sr,sg,sf
   ! Boundary Parameters
   logical, save :: in_Boundary = .false.
   logical, save :: in_North = .false.
   logical, save :: in_South = .false.
   logical, save :: in_West = .false.
   logical, save :: in_East = .false.
   character(len=128), save :: bnorth, bsouth, beast, bwest

contains

   ! ============================================================================
   subroutine startDocument_handler

   end subroutine startDocument_handler
   ! ============================================================================
   subroutine startElement_handler(namespaceURI, localname, name, atts)

      character(len=*), intent(in) :: namespaceURI
      character(len=*), intent(in) :: localname
      character(len=*), intent(in) :: name
      type(dictionary_t), intent(in) :: atts

      ! Refine Element
      if (name=='TIN_refine') in_RMesh = .true.
      if (name=='TIN_refine') refinesection = .true.
      if(in_RMesh) call SrmeshElement_handler(name)
      ! Boundary Element
      if (name=='Boundary') in_Boundary = .true.
      if(in_Boundary) call SboundaryElement_handler(name)

   end subroutine startElement_handler
   ! ============================================================================
   subroutine endElement_handler(namespaceURI, localname, name)

      character(len=*), intent(in) :: namespaceURI
      character(len=*), intent(in) :: localname
      character(len=*), intent(in) :: name

      ! Refine Element
      call ErmeshElement_handler(name)
      ! Boundary Element
      call EboundaryElement_handler(name)

   end subroutine endElement_handler
   ! ============================================================================
   subroutine characters_handler(chars)

      character(len=*), intent(in) :: chars

      ! Refine Element
      if(in_RMesh) call rmesh_characters_handler(chars)
      ! Get Boundary Element
      if(in_Boundary) call boundary_characters_handler(chars)

   end subroutine characters_handler
   ! ============================================================================
   subroutine SrmeshElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='HighRes') in_HR = .true.
      if (name=='LowRes') in_LR = .true.
      if (name=='FlowAcc') in_SG = .true.
      if (name=='StepFill') in_SF = .true.

   end subroutine SrmeshElement_handler
   ! ============================================================================
   subroutine SboundaryElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='North') in_North = .true.
      if (name=='South') in_South = .true.
      if (name=='East') in_East = .true.
      if (name=='West') in_West = .true.

   end subroutine SboundaryElement_handler
   ! ============================================================================
   subroutine ErmeshElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='TIN_refine') in_RMesh = .false.
      if (name=='HighRes') in_HR = .false.
      if (name=='LowRes') in_LR = .false.
      if (name=='FlowAcc') in_SG = .false.
      if (name=='StepFill') in_SF = .false.

   end subroutine ErmeshElement_handler
   ! ============================================================================
   subroutine EboundaryElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='North') in_North = .false.
      if (name=='South') in_South = .false.
      if (name=='East') in_East = .false.
      if (name=='West') in_West = .false.

   end subroutine EboundaryElement_handler
   ! ============================================================================
   subroutine rmesh_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if( in_HR)then
         hr = chars
         call rts(hr,grd_HR)
      elseif (in_LR) then
         lr = chars
         call rts(lr,grd_LR)
      elseif (in_SG) then
         sg = chars
         call rts(sg,accumulation_refine)
      elseif (in_SF) then
         sf = chars
         call rts(sf,step_fill)
      endif

   end subroutine rmesh_characters_handler
   ! ============================================================================
   subroutine boundary_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_North) then
         bnorth = chars
         call rts(bnorth,bounds(1))
      elseif (in_South) then
         bsouth = chars
         call rts(bsouth,bounds(2))
      elseif (in_East) then
         beast = chars
         call rts(beast,bounds(4))
      elseif (in_West) then
         bwest = chars
         call rts(bwest,bounds(3))
      endif

   end subroutine boundary_characters_handler
   ! ============================================================================
   subroutine endDocument_handler

   end subroutine endDocument_handler
   ! ============================================================================
   !> Subroutine xml_TIN()
   !! reads input that defines the SP Model experiment (XmL)
   !<
   ! ============================================================================
   subroutine xml_TIN

      type(xml_t) :: xf
      integer :: ios

      ! Open file
      call open_xml_file(xf, finput , ios)
      if (ios/=0) then
         print*,'---------------------'
         print*, 'Error opening input file for parsing XmL'
         print*,'---------------------'
      endif
      ! Parse file at first and find data
      call parse(xf, &
      startDocument_handler = startDocument_handler, &
      startElement_handler = startElement_handler, &
      endElement_handler = endElement_handler, &
      characters_handler = characters_handler, &
      endDocument_handler = endDocument_handler &
      )
      ! Close file
      call close_xml_t(xf)

   end subroutine xml_TIN
   ! ============================================================================
   !> Subroutine broadcast_tin()
   !! Broadcast initial parameters which define the TIN surface.
   !<
   ! ============================================================================
   subroutine broadcast_tin

      ! Refine parameters
      call mpi_bcast( grd_LR,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( grd_HR,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( accumulation_refine,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( step_fill,1,dbl_type,0,SPModel_comm_world,ierr )
      ! Broadcast boundaries parameters
      call mpi_bcast( bounds,4,int_type,0,SPModel_comm_world,ierr )

      return

   end subroutine broadcast_tin
  ! ============================================================================

end module TIN_parser
