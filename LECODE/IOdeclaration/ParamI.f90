! ============================================================================
! Name        : ParamI.f90
! Author      : Tristan Salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file ParamI.f90
!
! Description : ParamI is used to gather the information within the SPModel XmL input file to build
! the simulation parameters used in SPModel.
!
!<
! ============================================================================
module simulation_parser

   use file_data
   use FoX_sax
   use mpidata
   use time_data
   use FoX_common

   implicit none

   public

   ! Time Parameters
   logical, save :: timesection = .false.
   logical, save :: in_Time = .false.
   logical, save :: in_startTime = .false.
   logical, save :: in_endTime = .false.
   logical, save :: in_samplingInterval = .false.
   logical, save :: in_flowInterval = .false.
   logical, save :: in_rainInterval = .false.
   logical, save :: in_slpInterval = .false.
   logical, save :: in_displayTime = .false.
   logical, save :: in_outputTime = .false.
   character(len=128), save :: startTime,endTime,slopeInterval
   character(len=128), save :: displayTime,samplingInterval,flowInterval,rainInterval,outputTime

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

      ! Time Element
      if (name=='Time') in_Time = .true.
      if(in_Time) timesection = .true.
      if(in_Time) call StimeElement_handler(name)

   end subroutine startElement_handler

   !============================================================================
   subroutine endElement_handler(namespaceURI, localname, name)

      character(len=*), intent(in) :: namespaceURI
      character(len=*), intent(in) :: localname
      character(len=*), intent(in) :: name

      ! Time Element
      call EtimeElement_handler(name)

   end subroutine endElement_handler
   ! ============================================================================
   subroutine characters_handler(chars)

      character(len=*), intent(in) :: chars

      ! Get Time Parameters
      if (in_Time) call time_characters_handler(chars)

   end subroutine characters_handler
   ! ============================================================================
   subroutine StimeElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='startTime') in_startTime = .true.
      if (name=='endTime') in_endTime = .true.
      if (name=='layerTime') in_displayTime = .true.
      if (name=='outputTime') in_outputTime = .true.
      if (name=='riverInterval') in_flowInterval = .true.
      if (name=='rainInterval') in_rainInterval = .true.
      if (name=='samplingInterval') in_samplingInterval = .true.
      if (name=='SlopeInterval') in_slpInterval = .true.

   end subroutine StimeElement_handler
   ! ============================================================================
   subroutine EtimeElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='Time') in_Time = .false.
      if (name=='startTime') in_startTime = .false.
      if (name=='endTime') in_endTime = .false.
      if (name=='layerTime') in_displayTime = .false.
      if (name=='outputTime') in_outputTime = .false.
      if (name=='riverInterval') in_flowInterval = .false.
      if (name=='rainInterval') in_rainInterval = .false.
      if (name=='samplingInterval') in_samplingInterval = .false.
      if (name=='SlopeInterval') in_slpInterval = .false.

   end subroutine EtimeElement_handler
   ! ============================================================================
   subroutine time_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_startTime) then
         startTime = chars
         call rts(startTime,time_start)
      elseif(in_endTime)then
         endTime = chars
         call rts(endTime,time_end)
      elseif(in_outputTime)then
         outputTime = chars
         call rts(outputTime,time_output)
      elseif(in_displayTime)then
         displayTime = chars
         call rts(displayTime,time_display)
      elseif(in_flowInterval)then
         flowInterval = chars
         call rts(flowInterval,flow_int)
      elseif(in_rainInterval)then
         rainInterval = chars
         call rts(rainInterval,rain_int)
      elseif(in_samplingInterval)then
         samplingInterval = chars
         call rts(samplingInterval,sampling_int)
      elseif(in_slpInterval)then
         slopeInterval = chars
         call rts(slopeInterval,slope_int)
      endif

   end subroutine time_characters_handler
   ! ============================================================================
   subroutine endDocument_handler

   end subroutine endDocument_handler
   ! ============================================================================
   !> Subroutine xml_Sim()
   !! reads input that defines the SP Model experiment (XmL)
   !<
   ! ============================================================================
   subroutine xml_Sim

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

   end subroutine xml_Sim
   ! ============================================================================
   !> Subroutine broadcast_sim()
   !! Broadcast initial parameters which define the simulation parameters.
   !<
   ! ============================================================================
   subroutine broadcast_sim

      ! Broadcast time parameters
      call mpi_bcast( timesection,1,lgc_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( time_start,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( time_end,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( time_display,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( time_output,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( flow_int,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( sampling_int,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( slope_int,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( rain_int,1,dbl_type,0,SPModel_comm_world,ierr )
      tnow = time_start

      ! Check events time steps
      if( flow_int > time_display ) flow_int = time_display
      if( time_output < time_display ) time_output = time_display
      if( sampling_int > time_display ) sampling_int = time_display
      if( slope_int > time_display ) slope_int = time_display
      if( rain_int > time_display ) rain_int = time_display

      return

   end subroutine broadcast_sim
  ! ============================================================================

end module simulation_parser
