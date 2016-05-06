! ============================================================================
! Name        : ClassForces.f90
! Author      : tristan salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
! ============================================================================
!> \file ClassForces.f90
!
! Description :  It encapsulates the external forces module used by SPModel during compilation.
!
!<
! ============================================================================
!> Module extforces_data gives the sea-level fluctuation and displacements data
!<
module forces_data

   use precision_data

   implicit none

   public

   !> Sea level fluctuations type
   type sea_fluc
      ! Output sea level
      logical  :: output
      ! Flag for sea level fluctuations
      logical :: sealevel
      ! Number of sea fluctuation events
      integer  :: event
      ! Old sea level value
      real( tkind ) :: last_sea
      ! Actual sea level value
      real( tkind ) :: actual_sea
   end type sea_fluc
   type( sea_fluc ) :: gsea
   ! Sea-level elevation data
   real( tkind ),dimension(:,:),allocatable :: sealvl

   !> Temperature and salinity fluctuations type
   type temperature_salinity_fluc
      ! Flag for sea level fluctuations
      logical :: tpsl
      ! Old temperature and salinity values
      real( tkind ), dimension( 2 ) :: last_tempsal
      ! Number of temperature and salinity fluctuation events
      integer  :: event
      ! Actual sea level value
      real( tkind ), dimension( 2 ) :: actual_tempsal
   end type temperature_salinity_fluc
   type( temperature_salinity_fluc ) :: gtempsal

   ! Sea temperature salinity data
   real( tkind ),dimension(:,:),allocatable :: tempsal

   !> Porosity type
   type porosity_lt
      ! Flag for sea level fluctuations
      logical :: compaction
      ! Effective pressure number
      integer :: ePnb
   end type porosity_lt
   type( porosity_lt ) :: gporo
   ! Effective pressure
   real( tkind ),dimension(:),allocatable :: effPressure
   ! Porosity table values
   real( tkind ),dimension(:,:),allocatable :: porosity

   !> Vertical displacement type
   type vert_disp
      ! Flag for displacements
      logical :: vdisp
      ! Actual displacement event
      integer :: actual
      ! Displacements event numbers
      integer  :: event
      ! Time from last displacement call
      real( tkind ) :: lastdisp
      ! Time elapsed between current and previous displacement time
      real( tkind ) :: disptime
   end type vert_disp
   type( vert_disp ) :: gdisp, vdisp
   ! Displacement file unit to get a continuous displacement field
   integer( tkind ),dimension(:),allocatable :: gdisp_fill, vdisp_fill
   ! Displacement time values in the input file
   real( tkind ),dimension(:,:),allocatable :: gdisp_time, vdisp_time

   !> Sea-level parameters
   type sl_par
      ! Lower sea-level elevation
      real( tkind )  :: sea1
      ! Upper sea-level elevation
      real( tkind ) :: sea2
      ! Lower sea-level time
      real( tkind ) :: time1
      ! Upper sea-level time
      real( tkind ) :: time2
      ! Simulation current time
      real( tkind ) :: tsim
      ! Interpolated sea-level elevation
      real( tkind ) :: sl_int
   end type sl_par

   !> Temperature and salinity parameters
   type st_par
      ! Lower temperature and salinity
      real( tkind ), dimension( 2 )  :: ts1
      ! Upper temperature and salinity
      real( tkind ), dimension( 2 ) :: ts2
      ! Lower temperature and salinity time
      real( tkind ) :: time1
      ! Upper temperature and salinity time
      real( tkind ) :: time2
      ! Simulation current time
      real( tkind ) :: tsim
      ! Interpolated temperature and salinity values
      real( tkind ), dimension( 2 ) :: ts_int
   end type st_par

end module forces_data
! ============================================================================
!> Module ocean_data set up parameters for ocean coupling
!<
module ocean_data

    use precision_data

    implicit none

    public

    ! Ocean plugin flag
    logical :: ocean_plug

    ! Ocean forces flag
    logical :: current_on
    logical :: wave_on

    ! Forecast parameters
    integer :: forecast_nb

    ! Active hindcast value
    integer :: active_hindcast

    ! Ocean module morphodynamic thickness
    real( tkind ) :: ocean_morph_th

    ! Ocean module morphodynamic time step
    real( tkind ) :: ocean_time_step

    ! Number of hindcast scenarios
    type hindcast_param
        ! Percentage of wind/wave subgroup.
        real :: perc
    end type hindcast_param

    ! Maximum depth of wave/current action
    real( tkind ) :: max_WC_action_depth

    type hindcast_def
        ! Class number
        integer :: cnb
        ! Time start
        real :: tstart
        ! Time end
        real :: tend
        ! Subgroud parameters
        type( hindcast_param ), dimension( : ), allocatable :: subgroup
    end type hindcast_def
    type( hindcast_def ), dimension( : ), allocatable :: hindcast

    ! Ocean currents and waves arrays
    real( tkind ),dimension(:,:,:),allocatable :: ocean_current, ocean_wave

    ! Ocean currents and waves total load components
    real( tkind ),dimension(:,:),allocatable :: soulsby_total_load, total_load
    real( tkind ),dimension(:,:),allocatable :: transportX, transportY



end module ocean_data
! ============================================================================
