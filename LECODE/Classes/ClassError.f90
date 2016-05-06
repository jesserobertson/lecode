! ============================================================================
! Name        : ClassError.f90
! Author      : tristan salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file ClassError.f90
!
! Description :  It encapsulates the error module used by SPModel during compilation.
!
!<
! ============================================================================
!> Module error_data gives error signal for common SPModel problems.
!<
module error_data

   use mpidata
   use precision_data

   implicit none

   public

   ! Error when reading name of experiment
   integer, parameter :: ARG_FILE = 1
   ! Error when reading name of displacements file
   integer, parameter :: DISP_FILE = 2
   ! Error in reading the format of the displacements file
   integer, parameter :: DISP_TIME1 = 3
   ! Error when reading the time period of the displacements field
   integer, parameter :: DISP_TIME2 = 4
   ! Error when reading name of deposition file
   integer, parameter :: DEPO_FILE = 14
   ! Error declaration Stratal grid
   integer, parameter :: STRAT_DEC = 15
   ! Error looking for Stratal grid
   integer, parameter :: STRAT_FILE = 16
   ! Error when reading name of sea level file
   integer, parameter :: SEA_EOF = 17
   ! Error reading sea level data
   integer, parameter :: SEA_ORDER = 18
   ! Error looking for sea level file
   integer, parameter :: SEA_FILE = 19
   ! Error section time
   integer, parameter :: SEC_TIME = 20
   ! Error when reading sources number
   integer, parameter :: SRCNBS = 22
   ! Error when reading sources percentages XML file
   integer, parameter :: SRCPER = 23
   ! Error when reading sources time XML file
   integer, parameter :: SRCTIME = 24
   ! Error when reading sources x poistion XML file
   integer, parameter :: SRCXPOS = 25
   ! Error when reading sources y poistion XML file
   integer, parameter :: SRCYPOS = 26
   ! Error when defining the simulation start and end time
   integer, parameter :: TIME_DEF = 27
   ! Error when defining the simulation display time step
   integer, parameter :: TIME_DISP = 28
   ! Error when defining the simulation sampling interval
   integer, parameter :: TIME_SPI = 29
   ! Error when looking for rainfall file
   integer, parameter :: RAIN_FILE = 30
   ! Error when looking for rainfall cyclicity
   integer, parameter :: RAIN_CYC = 31
   ! Error when looking for rainfall file
   integer, parameter :: IN_MAXFW = 32
   ! Error when computing plane equatuion
   integer, parameter :: CBASEVAL = 33
   ! Error computing next time step in CKRK algo
   integer, parameter :: DTNNEG = 34
   ! Error computing  intersecting line
   integer, parameter :: INTERSECTING = 35
   ! Error interpolation TIN grid
   integer, parameter :: TIN_INTERP = 5
   ! Error resolution TIN grid
   integer, parameter :: TIN_RESO1 = 6
   ! Error resolution TIN grid
   integer, parameter :: TIN_RESO2 = 7
   ! Error SW corner TIN grid
   integer, parameter :: TIN_CORNER = 8
   ! Error declaration TIN grid
   integer, parameter :: TIN_DEC = 9
   ! Error spacing resoltion of the TIN grid
   integer, parameter :: TIN_SPAC = 10
   ! Error nodes number of the TIN grid
   integer, parameter :: TIN_NOD = 11
   ! Error when calling triangle
   integer, parameter :: TIN_TRIC = 12
   ! Error looking for TIN grid file
   integer, parameter :: TIN_FILE = 13
   ! Error looking for TIN elements
   integer, parameter :: TIN_INPUT = 21

contains

   ! ============================================================================
   !> Subroutine completion
   !! Subroutine completion synchronizes multiple processors.
   !! It also keeps track of incurred errors in each processor.
   !<
   ! ============================================================================
   subroutine completion

      integer :: every, contrib

      every = 0
      contrib = attempt

      if(1 == nproc)then
         every = attempt
      else
         call mpi_allreduce(contrib, every, 1, int_type, max_type, SPModel_comm_world,ierr)
      endif
      if( every > 0 )then
         attempt = every
         if( iam == 0 .and. attempt > 0 ) write(*,*)'An error occurred with code = ',attempt
         if( iam == 0 ) call error_signal
         call mpi_finalize( ierr )
         stop
      endif

   end subroutine completion
   ! ============================================================================
   !> Subroutine error_signal
   !! Subroutine error_signal gives information about the error in SPModel.
   !<
   ! ============================================================================
   subroutine error_signal

      print*,'-----------------------------------------------'
      print*,'Explanation of the error :'
      if( attempt == ARG_FILE )&
      print*,'--> Error reading argument file.'
      if( attempt == SEC_TIME )&
      print*,'--> Error Time section missing in the input file.'
      if( attempt == SEA_FILE )&
      print*,'--> Error sea level file does not exist.'
      if( attempt == SEA_ORDER )&
      print*,'--> Error sea level file file not ordered in increasing time.'
      if( attempt == SEA_EOF )&
      print*,'--> Error during allocation of sea-level variables.'
      if( attempt == STRAT_FILE )&
      print*,'--> Error StrataGrid file does not exist.'
      if( attempt == STRAT_DEC )&
      print*,'--> Error StrataGrid nb of nodes not equal discretisation number for X and Y.'
      if( attempt == DEPO_FILE )&
      print*,'--> Error deposit file does not exist.'
      if( attempt == DISP_TIME1 )&
      print*,'--> Error reading displacement times.'
      if( attempt == DISP_TIME2 )&
      print*,'--> Error reading displacement times.'
      if( attempt == DISP_FILE )&
      print*,'--> Error input file for vertical displacement fields cannot be found.'
      if( attempt == SRCNBS )&
      print*,'--> Error number of sources differs from the allocated total sources number.'
      if( attempt == SRCPER )&
      print*,'--> Error sediment ratios dont add to 100% for a source.'
      if( attempt == SRCTIME )&
      print*,'--> Error source start time > end time.'
      if( attempt == SRCXPOS )&
      print*,'--> Error X position of source out of range.'
      if( attempt == SRCYPOS )&
      print*,'--> Error Y position of source out of range.'
      if( attempt == TIME_DEF )&
      print*,'--> Error start time > end time.'
      if( attempt == TIME_DISP )&
      print*,'--> Error in display and/or output time definition.'
      if( attempt == TIME_SPI )&
      print*,'--> Error in sampling interval definition.'
      if( attempt == RAIN_FILE )&
      print*,'--> Error when looking for rainfall file.'
      if( attempt == RAIN_CYC )&
      print*,'--> Error when reading cyclicity in rainfall file.'
      if( attempt == IN_MAXFW )&
      print*,'--> Error to many flow walkers to allocate.'
      if( attempt == CBASEVAL )&
      print*,'--> Error computing plane equation.'
      if( attempt == DTNNEG )&
      print*,'--> Error computing next time step.'
      if( attempt == INTERSECTING )&
      print*,'--> Error computing intersecting line.'
      if( attempt == TIN_INPUT )&
      print*,'--> Error <TIN_surface> or <TIN_refine> element missing in the XmL.'
      if( attempt == TIN_FILE )&
      print*,'--> Error TIN_Grid file does not exist.'
      if( attempt == TIN_TRIC )&
      print*,'--> Error loading the file created by Triangle.'
      if( attempt == TIN_NOD )&
      print*,'--> Error TIN nodes not ordered in increasing X values first.'
      if( attempt == TIN_SPAC )&
      print*,'--> Error TIN nodes spacing on the border differs from the value given in the input file.'
      if( attempt == TIN_DEC )&
      print*,'--> Error the TIN_Grid declaration.'
      if( attempt == TIN_CORNER )&
      print*,'--> Error TIN SW corner do not match with stratigraphic mesh.'
      if( attempt == TIN_RESO1 )&
      print*,'--> Error resolution along borders not multiple of stratal ones.'
      if( attempt == TIN_RESO2 )&
      print*,'--> Error high/low resolution space not multiple of stratal spacing.'
      if( attempt == TIN_INTERP )&
      print*,'--> Error not enough point for interpolation.'
      print*,'-----------------------------------------------'

      return

   end subroutine error_signal
  ! ============================================================================

end module error_data
! ============================================================================

