! ============================================================================
! Name        : ClassGeneral.f90
! Author      : tristan salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file ClassGeneral.f90
!
! Description :  It encapsulates some general parameters used by SPModel during compilation.
!
!<
! ============================================================================
!> Module precision_data gives the data precision type (single or double)
!<
module precision_data

   ! INTEL COMPILER
!   use ifport

   ! Single or double precision
   integer, parameter :: singlep = kind(0.0)
   integer, parameter :: doublep = kind(0.0d0)

   ! Dont't forget to change all dbl_type in mpi_real
   ! and all h5t_native_double in h5t_native_real
   ! and in kdtree2.f90 change sp to dp at line 22

!   integer, parameter :: tkind = singlep
   integer, parameter :: tkind = doublep
   integer, parameter :: sizesp = 4
   integer, parameter :: sizedp = 8

   ! MPI integer type communicator
   integer :: int_type
   ! MPI double type communicator
   integer :: dbl_type
   ! MPI logical type communicator
   integer :: lgc_type
   ! MPI max type communicator
   integer :: max_type
   ! MPI min type communicator
   integer :: min_type
   ! MPI sum type communicator
   integer :: sum_type
   ! SPModel communicator
   integer :: SPModel_comm_world
   ! New communicator
   integer :: new_comm

contains

   ! ============================================================================
   !> Subroutine term_command.
   !!
   !<
   ! ============================================================================
   subroutine term_command( cmds )

        logical( 4 ) :: result
        character(len=128) :: cmds

        result = .false.

        ! INTEL FORTRAN COMPILER
!        result = systemqq( cmds )
        ! GNU FORTRAN COMPILER
        call system( cmds )

        return

   end subroutine term_command
   ! ============================================================================
   !> Subroutine count_lines_file.
   !! Counts the number of line in a given file.
   !<
   ! ============================================================================
   subroutine count_lines_file( filetoread, nlines )

        integer :: lunit, stat, nlines
        character(len=128) :: filetoread

        lunit = 12098
        nlines = 0
        open( lunit, file = filetoread )
        lp: do
            read( lunit,*, iostat = stat )
            if( stat /= 0 ) exit lp
            nlines = nlines + 1
        enddo lp

        close (lunit)

        return

   end subroutine count_lines_file
   ! ============================================================================

end module precision_data
! ============================================================================
!> Module mpi encapsulates mpi data as well as SPModel module
!<
module mpidata

#include "mpif.h"

   ! Error code
   integer :: ierr
   ! Processor ID
   integer :: iam
   ! Number of processors
   integer :: nproc
   ! Track error messages through simulation
   integer :: attempt

   ! Declaration of created mpi carthesian grid
   integer :: dim_size(0:1), coords(0:1)
   integer, parameter :: ndims = 2
   logical, parameter :: reorder = .true.
   logical, dimension(2), parameter :: periods(0:1) =  .false.

   ! Nb of bytes from the beginning of the structure
   integer, dimension(:), allocatable :: disp

end module mpidata
! ============================================================================
!> Module time_data encapsulates characteristics of the simulation time
!<
module time_data

   use precision_data

   ! Check if outputs should be written
   logical :: output_sim
   ! Check if check pointing is needed
   logical :: checkpointing
   ! Check if horizontal deformation is needed
   logical :: deformlayer
   ! Check if restart is needed
   logical :: restart
   ! Check if underworld is needed
   logical :: udw
   ! Check if underworld coupling is needed
   integer :: udw_coupling
   ! Check pointing frequency
   integer :: checkfreq
   ! Restart iteration
   integer :: restart_iter
   ! Number of processor from previous run
   integer :: restart_proc
   ! Start and end time of the simulation
   real( tkind ) :: time_start, time_end
   ! Layer time step
   real( tkind ) :: time_display
   ! Output time step
   real( tkind ) :: time_output
   ! Next time step
   real( tkind ) :: time_next
   ! Rain sampling interval
   real( tkind ) :: rain_int
   ! Flow sampling interval
   real( tkind ) :: flow_int
   ! Sampling interval
   real( tkind ) :: sampling_int
   ! Diffusion slope interval
   real( tkind ) :: slope_int
   ! Current simulation time
   real( tkind ) :: tnow
   ! Simulation time step
   real( tkind ) :: dt
   ! Simulation next time step
   real( tkind ) :: dtnext
   ! Simulation coarse time step
   real( tkind ) :: coarse_dt
   ! Simulation next coarse time step
   real( tkind ) :: next_coarse_dt
   ! Simulation fine time step
   real( tkind ) :: fine_dt
   ! Next displacement time
   real( tkind ) :: next_displacement
   ! Next ocean circulation time
   real( tkind ) :: next_oceancirc
   ! Next organic/carbonate time
   real( tkind ) :: next_carborg
   ! Next slope diffusion time
   real( tkind ) :: next_slpdiff
   ! Next layer time
   real( tkind ) :: next_display
   ! Next rain time
   real( tkind ) :: next_rain
   ! Next output time
   real( tkind ) :: next_output
   ! Next hemipelagites time
   real( tkind ) :: last_hemi
   ! Carbonates/organics time step
   real( tkind ) :: time_carborg

end module time_data
! ============================================================================
!> Module param_data encapsulates parameters definition
!<
module param_data

   use precision_data

   ! Embeds the previously defined flag values
   integer :: nextevent, oldevent
   ! Undefined flag
   integer, parameter :: SPM_UNDEFINED = 0
   ! Quit simulation flag
   integer, parameter :: SPM_QUIT = 1
   ! Displacement event flag
   integer, parameter :: SPM_DISPLACEMENT = 2
   ! Display event flag
   integer, parameter :: SPM_DISPLAY = 3
   ! Stepping flag
   integer, parameter :: SPM_STEP = 4
   ! Initialisation flag
   integer, parameter :: SPM_INIT = 5
   ! Flow walker flag
   integer, parameter :: SPM_FLUID = 6
   ! Source flow walkers flag
   integer, parameter :: SPM_INFLOW = 7
   ! Sea level fluctuation time step flag
   integer, parameter :: SPM_SEALEVEL = 8
   ! Slope diffusion time step flag
   integer, parameter :: SPM_SLPDIFFUSION = 9
   ! Rain flow walkers flag
   integer, parameter :: SPM_RAINFALL = 10
   ! Ocean circulation flag
   integer, parameter :: SPM_OCEAN = 11
   ! Carbonate/organic evolution flag
   integer, parameter :: SPM_CARBORG = 12

   ! Zero check accuracy
   real( tkind ), parameter :: tor = 1.0e-8
   ! Zero time check accuracy
   real( tkind ) :: time_tolerance
   ! Number of seconds per year
   real( tkind ), parameter :: secyear = 31536000.0_8
   ! Inverse of number of seconds per year
   real( tkind ), parameter :: secyrinv = 1.0_8 / 31536000.0_8
   ! Pi number
   real( tkind ), parameter :: pi = 4.0_8 * atan( 1.0_8 )

   ! Type of deposition processes classification
   ! Alluvial deposits
   integer, parameter :: tpriv = 2
   ! Rain deposits
   integer, parameter :: tprain = 100
   ! Turbidites
   integer, parameter :: tpturb = 300
   ! Plume deposits ZFE
   integer, parameter :: tpplum1 = 400
   ! Plume deposits ZEF
   integer, parameter :: tpplum2 = 401

   ! Step counter and flags
   logical :: new_disp, fexistx, flagx
   integer :: stepout, stepin, rcd, rcdin, step, stepc


end module param_data
! ============================================================================

