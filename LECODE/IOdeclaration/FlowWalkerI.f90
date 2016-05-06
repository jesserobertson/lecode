! ============================================================================
! Name        : FlowWalkerI.f90
! Author      : Tristan Salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file FlowWalkerI.f90
!
! Description : FlowWalkerI is used to gather the information within the SPModel XmL input file to build
! the external flow walker variables used in SPModel.
!
!<
! ============================================================================
module FWLK_parser

   use file_data
   use mpidata
   use FoX_sax
   use flux_data
   use fwalker_data
   use FoX_common

   implicit none

   public

   ! Density Parameters
   logical, save :: densitysection = .false.
   logical, save :: in_fDensity = .false.
   logical, save :: in_fwdensity = .false.
   logical, save :: in_seadensity = .false.
   character(len=128), save :: seadensity, fwdensity
   ! Plume Parameters
   logical, save :: in_Plume = .false.
   logical, save :: in_oceanVx = .false.
   logical, save :: in_oceanVy = .false.
   logical, save :: in_riverMouth = .false.
   logical, save :: in_plumeLenght = .false.
   character(len=128), save :: oceanVx, oceanVy, riverMouth, plumeLenght
   ! Manning Parameters
   logical, save :: in_Manning = .false.
   logical, save :: in_Open = .false.
   logical, save :: in_Hyper = .false.
   logical, save :: in_Hypo = .false.
   character(len=128), save :: mopen, mhyper, mhypo
   ! Transport Parameters
   logical, save :: transportsection = .false.
   logical, save :: in_Transport = .false.
   logical, save :: in_fwmxh = .false.
   logical, save :: in_fwmnh = .false.
   logical, save :: in_fwmxw = .false.
   logical, save :: in_fwmnw = .false.
   logical, save :: in_fwmnv = .false.
   logical, save :: in_fwmxs = .false.
   logical, save :: in_sedload = .false.
   logical, save :: in_erolim = .false.
   logical, save :: in_hdlim = .false.
   character(len=128), save :: fwmxh, fwmnh, fwmxw, fwmnw, fwmnv
   character(len=128), save :: sedload, erolim, fwmxs, hdlim

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

      ! Density Element
      if (name=='WaterDensities') in_fDensity = .true.
      if(in_fDensity) densitysection = .true.
      if(in_fDensity) call SdensityElement_handler(name)
      ! Plume Element
      if (name=='Plume') in_Plume = .true.
      if(in_Plume) call SplumeElement_handler(name)
      ! Manning Element
      if (name=='Manning') in_Manning = .true.
      if(in_Manning) call SmanningElement_handler(name)
      ! Transport Element
      if (name=='SedimentTransportParameters') in_Transport = .true.
      if(in_Transport) transportsection = .true.
      if(in_Transport) call StransportElement_handler(name)

   end subroutine startElement_handler
   ! ============================================================================
   subroutine endElement_handler(namespaceURI, localname, name)

      character(len=*), intent(in) :: namespaceURI
      character(len=*), intent(in) :: localname
      character(len=*), intent(in) :: name

      ! Density Element
      call EdensityElement_handler(name)
      ! Plume Element
      call EplumeElement_handler(name)
      ! Manning Element
      call EmanningElement_handler(name)
      ! Transport Element
      call EtransportElement_handler(name)

   end subroutine endElement_handler
   ! ============================================================================
   subroutine characters_handler(chars)

      character(len=*), intent(in) :: chars

      ! Get Density Element
      if(in_fDensity) call density_characters_handler(chars)
      ! Get Plume Element
      if(in_Plume) call plume_characters_handler(chars)
      ! Get Manning Element
      if(in_Manning) call manning_characters_handler(chars)
      ! Get Transport Element
      if(in_Transport) call transport_characters_handler(chars)

   end subroutine characters_handler
   ! ============================================================================
   subroutine SdensityElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='enteringFluid') in_fwdensity = .true.
      if (name=='seaDensity') in_seadensity = .true.

   end subroutine SdensityElement_handler
   ! ============================================================================
   subroutine SplumeElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='oceanVx') in_oceanVx = .true.
      if (name=='oceanVy') in_oceanVy = .true.
      if (name=='riverMouth') in_riverMouth = .true.
      if (name=='plumeLenght') in_plumeLenght = .true.

   end subroutine SplumeElement_handler
   ! ============================================================================
   subroutine SmanningElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='openchannelFlow') in_Open = .true.
      if (name=='hyperpycnalFlow') in_Hyper = .true.
      if (name=='hypopycnalFlow') in_Hypo = .true.

   end subroutine SmanningElement_handler
   ! ============================================================================
   subroutine StransportElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='streamMinDepth') in_fwmnh = .true.
      if (name=='streamMaxDepth') in_fwmxh = .true.
      if (name=='streamMinWidth') in_fwmnw = .true.
      if (name=='streamMaxWidth') in_fwmxw = .true.
      if (name=='maxFWStepPerCell') in_fwmxs = .true.
      if (name=='flowWalkerMinVelocity') in_fwmnv = .true.
      if (name=='minSedimentLoad') in_sedload = .true.
      if (name=='flowWalkerErosionLimit') in_erolim = .true.
      if (name=='flowWalkerDepositionLimit') in_hdlim = .true.

   end subroutine StransportElement_handler
   ! ============================================================================
   subroutine EdensityElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='enteringFluid') in_fwdensity = .false.
      if (name=='seaDensity') in_seadensity = .false.

   end subroutine EdensityElement_handler
   ! ============================================================================
   subroutine EplumeElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='oceanVx') in_oceanVx = .false.
      if (name=='oceanVy') in_oceanVy = .false.
      if (name=='riverMouth') in_riverMouth = .false.
      if (name=='plumeLenght') in_plumeLenght = .false.

   end subroutine EplumeElement_handler
   ! ============================================================================
   subroutine EmanningElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='openchannelFlow') in_Open = .false.
      if (name=='hyperpycnalFlow') in_Hyper = .false.
      if (name=='hypopycnalFlow') in_Hypo = .false.

   end subroutine EmanningElement_handler
   ! ============================================================================
   subroutine EtransportElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='streamMinDepth') in_fwmnh = .false.
      if (name=='streamMaxDepth') in_fwmxh = .false.
      if (name=='streamMinWidth') in_fwmnw = .false.
      if (name=='streamMaxWidth') in_fwmxw = .false.
      if (name=='maxFWStepPerCell') in_fwmxs = .false.
      if (name=='flowWalkerMinVelocity') in_fwmnv = .false.
      if (name=='minSedimentLoad') in_sedload = .false.
      if (name=='flowWalkerErosionLimit') in_erolim = .false.
      if (name=='flowWalkerDepositionLimit') in_hdlim = .false.

   end subroutine EtransportElement_handler
   ! ============================================================================
   subroutine density_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_fwdensity) then
         fwdensity = chars
         call rts(fwdensity,fluid_density)
      elseif (in_seadensity) then
         seadensity = chars
         call rts(seadensity,sea_density)
      endif

   end subroutine density_characters_handler
   ! ============================================================================
   subroutine plume_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_oceanVx) then
         plume_cmpt = .true.
         oceanVx = chars
         call rts(oceanVx,ocean_flow(1))
      elseif (in_oceanVy) then
         oceanVy = chars
         call rts(oceanVy,ocean_flow(2))
      elseif (in_riverMouth) then
         riverMouth = chars
         call rts(riverMouth,river_mouth)
      elseif (in_plumeLenght) then
         plumeLenght = chars
         call rts(plumeLenght,plume_lenght)
      endif

   end subroutine plume_characters_handler
   ! ============================================================================
   subroutine manning_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_Open) then
         mopen = chars
         call rts(mopen,manning_open)
      elseif (in_Hyper) then
         mhyper = chars
         call rts(mhyper,manning_hyper)
      elseif (in_Hypo) then
         mhypo = chars
         call rts(mhypo,manning_hypo)
      endif

   end subroutine manning_characters_handler
   ! ============================================================================
   subroutine transport_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_fwmxh) then
         fwmxh = chars
         call rts(fwmxh,transport%dzimax)
      elseif (in_fwmnh) then
         fwmnh = chars
         call rts(fwmnh,transport%dzimin)
      elseif (in_fwmxw) then
         fwmxw = chars
         call rts(fwmxw,transport%wdtmax)
      elseif (in_fwmnw) then
         fwmnw = chars
         call rts(fwmnw,transport%wdtmin)
      elseif (in_fwmxs) then
         fwmxs = chars
         call rts(fwmxs,transport%stepmax)
      elseif (in_fwmnv) then
         fwmnv = chars
         call rts(fwmnv,transport%fvmin)
      elseif (in_sedload) then
         sedload = chars
         call rts(sedload,transport%slsmin)
      elseif (in_erolim) then
         erolim = chars
         call rts(erolim,transport%erosion_limit)
      elseif (in_hdlim) then
         hdlim = chars
         call rts(hdlim,transport%fac_limit)
      endif

   end subroutine transport_characters_handler
   ! ============================================================================
   subroutine endDocument_handler

   end subroutine endDocument_handler
   ! ============================================================================
   !> Subroutine xml_FWLK()
   !! reads input that defines the SP Model experiment (XmL)
   !<
   ! ============================================================================
   subroutine xml_FWLK

      type(xml_t) :: xf
      integer :: ios

      ocean_flow = 0.0_8
      transport%fac_limit = 1.0_8
      fluid_density = 1015.0_8
      sea_density = 1027.0_8
      manning_open = 0.025_8
      manning_hyper = 0.01_8
      manning_hypo = 0.08_8
      transport%dzimin = 0.1_8
      transport%dzimax = 10.0_8
      transport%wdtmin = 1.0_8
      transport%wdtmax = 500.0_8
      transport%fvmin = 0.00001_8
      transport%morpho = 5.0e10_8
      transport%slsmin = 1.0e-12_8
      transport%erosion_limit = 0.25_8
      transport%stepmax = 50

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

   end subroutine xml_FWLK
   ! ============================================================================
   !> Subroutine broadcast_fwlk()
   !! Broadcast initial parameters which define the flow walker.
   !<
   ! ============================================================================
   subroutine broadcast_fwlk

      transport%sed_dt = 1.0e4_8

      ! Broadcast densities parameters
      call mpi_bcast( fluid_density,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( sea_density,1,dbl_type,0,SPModel_comm_world,ierr )
      ! Broadcast manning coefficient parameters
      call mpi_bcast( manning_open,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( manning_hyper,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( manning_hypo,1,dbl_type,0,SPModel_comm_world,ierr )
      ! Broadcast sediment transport parameters
      call mpi_bcast( transport%dzimin,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( transport%dzimax,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( transport%wdtmin,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( transport%wdtmax,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( transport%fvmin,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( transport%slsmin,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( transport%erosion_limit,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( transport%stepmax,1,int_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( transport%fac_limit,1,dbl_type,0,SPModel_comm_world,ierr )
      transport%slsmin = transport%slsmin * 1.e9_8
      ! Broadcast plume parameters
      call mpi_bcast( plume_cmpt,1,lgc_type,0,SPModel_comm_world,ierr )
      if( plume_cmpt )then
          call mpi_bcast( ocean_flow,2,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( river_mouth,1,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( plume_lenght,1,dbl_type,0,SPModel_comm_world,ierr )
      else
          ocean_flow = 0.0_8
          river_mouth = 0.0_8
          plume_lenght = 0.0_8
      endif

      return

   end subroutine broadcast_fwlk
  ! ============================================================================

end module FWLK_parser
