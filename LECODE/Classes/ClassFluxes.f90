! ============================================================================
! Name        : ClassFluxes.f90
! Author      : tristan salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file ClassFluxes.f90
!
! Description :  It encapsulates the influx used by SPModel during compilation.
!
!<
! ============================================================================
!> Module flux_data gives the influx data
!<
module flux_data

   use param_data
   use precision_data

   implicit none

   public

   logical :: geomorphology, rain_region, hemi_flag
   ! Flow walkers density
   real( tkind ) :: fluid_density
   ! Sea water density
   real( tkind ) :: sea_density
   ! Gravitational constant
   real( tkind ),parameter :: gravity = 9.81_8
   ! Manning coefficient for open channel flow
   real( tkind ) :: manning_open
   ! Manning coefficient for hyperpycnal flow
   real( tkind ) :: manning_hyper
   ! Manning coefficient for hypopycnal flow
   real( tkind ) :: manning_hypo
   ! Number of river sources
   integer :: num_src
   ! Maximum elevation of river src
   real( tkind ) :: max_elev_src
   ! Number of rain sources for the considered time
   integer :: rain_src
   ! Maximum number of rain sources for the considered time
   integer :: rain_total
   ! Maximum number of flow walker due to rain released each time
   integer :: rain_nb
   ! Number of rain events during the simulation
   integer :: rain_event
   ! Initial elevation of rain
   real( tkind ) :: rain_elev
   ! Rain flow accumulation threshold
   integer :: rain_facc
   ! Rain max combination value
   integer :: combine
   ! Rain max flow rate before being considered as a river stream
   real( tkind ) :: rainstreamQ

   !> Hemipelagic parameters
   integer :: contourhemi_nb
   real( tkind ),dimension(:,:),allocatable :: hemipelagic

   !> River source parameter type
   ! Souce splitting number
   integer,dimension(:),allocatable :: fws_num
   ! Source flow type
   integer,dimension(:),allocatable :: fws_type
   ! Source start time
   real( tkind ),dimension(:),allocatable :: fws_tstrt
   ! Source end time
   real( tkind ),dimension(:),allocatable :: fws_tend
   ! Source seasonal percentage
   real( tkind ),dimension(:),allocatable :: fws_perc
   ! Source x position
   real( tkind ),dimension(:),allocatable :: fws_xposition
   ! Source y position
   real( tkind ),dimension(:),allocatable :: fws_yposition
   ! Source x range
   real( tkind ),dimension(:),allocatable :: fws_xrange
   ! Source y range
   real( tkind ),dimension(:),allocatable :: fws_yrange
   ! Source height
   real( tkind ),dimension(:),allocatable :: fws_srch
   ! Source x velocity
   real( tkind ),dimension(:),allocatable :: fws_xvel
   ! Source y velocity
   real( tkind ),dimension(:),allocatable :: fws_yvel
   ! Source flow discharge
   real( tkind ),dimension(:),allocatable :: fws_qfl
   ! Source flow volume
   real( tkind ),dimension(:),allocatable :: fws_volume
   ! Source sediment concentration
   real( tkind ),dimension(:),allocatable :: fws_sedconc
   ! Source sediment percentage
   real( tkind ),dimension(:,:),allocatable :: fws_sedperc
   ! Source sediment discharge (m3/s)
   real( tkind ),dimension(:,:),allocatable :: fws_sedcharge

   !> Rain parameter type
   ! Source referenced node
   integer,dimension(:),allocatable :: frain_node
   ! Source x position
   real( tkind ),dimension(:),allocatable :: frain_xposition
   ! Source y position
   real( tkind ),dimension(:),allocatable :: frain_yposition
   ! Source flow discharge
   real( tkind ),dimension(:),allocatable :: frain_qfl

  ! Rainfall data type

  ! Number of rainfall region
  integer,dimension(:),allocatable :: rain_areaNb
  ! Starting time of the considered rainfall period
  real( tkind ),dimension(:),allocatable :: rain_tstart
  ! Ending time of the considered rainfall period
  real( tkind ),dimension(:),allocatable :: rain_tend
  ! Minimal X coordinates of the defined rainfall region
  real( tkind ),dimension(:,:),allocatable :: rain_xmin
  ! Maximal X coordinates of the defined rainfall region
  real( tkind ),dimension(:,:),allocatable :: rain_xmax
  ! Minimal Y coordinates of the defined rainfall region
  real( tkind ),dimension(:,:),allocatable :: rain_ymin
  ! Maximal Y coordinates of the defined rainfall region
  real( tkind ),dimension(:,:),allocatable :: rain_ymax
  ! Rate of rainfall in the considered region ( m/y )
  real( tkind ),dimension(:,:),allocatable :: rain_h
  ! Minimal rainfall altitude in the considered region ( m )
  real( tkind ),dimension(:,:),allocatable :: rain_min
  ! Maximal rainfall altitude in the considered region ( m )
  real( tkind ),dimension(:,:),allocatable :: rain_max

   !> Transport parameters type declaration
   type :: transport_parameters
      ! Sedimentation time step factor
      real( tkind ) :: sed_dt
      ! Sedimentation limited thickness factor
      real( tkind ) :: fac_limit
      ! Minimum depth of flow walkers
      real( tkind ) :: dzimin
      ! Maximum depth of flow walkers
      real( tkind ) :: dzimax
      ! Minimum width of flow walkers
      real( tkind ) :: wdtmin
      ! Maximum width of flow walkers
      real( tkind ) :: wdtmax
      ! Maximum number of steps for a FW
      integer :: stepmax
      ! Morphological factor
      real( tkind ) :: morpho
      ! Minimum velocity of flow walkers
      real( tkind ) :: fvmin
      ! Minimum ratio of sediment load
      real( tkind ) :: slsmin
      ! Erosion limitation parameter
      real( tkind ) :: erosion_limit
   end type transport_parameters
   type(transport_parameters) :: transport

end module flux_data
! ============================================================================
!> Module fwalker_data encapsulates characteristics of the flow walkers.
!<
module fwalker_data

   use flux_data
   use sediment_data
   use precision_data

   implicit none

   logical :: plume_cmpt

   ! Maximum number of flow walkers on the simulation area
   integer :: maxfw = 10000
   ! Maximum number of flow walkers for rain on the simulation area
   integer :: maxrain = 5000
   ! Maximum number of flow walkers recorded for output
   integer :: maxfwinfo = 50000
   ! Number of grid nodes recording flow walkers properties
   integer, parameter :: numold = 500
   ! Number of flow walkers per processor
   integer :: num_fw
   ! Number of flow walkers per processor for ero/dep
   integer :: num_fwed
   ! Number of flow walkers per processor which will be visualized
   integer :: num_fw_o
   ! Number of flow walkers per timestep
   integer :: tot_fw
   ! Number of active sources
   integer :: sourceactive
   ! Max number of stratigraphic nodes possibily activated for ero/dep
   ! for each FW
   integer :: maxPtsFW
   ! Maximum lenght of a considered plume
   real( tkind ) :: plume_lenght
   ! Averaged size of a river mouth
   real( tkind ) :: river_mouth
   ! Ocean flow current
   real( tkind ), dimension( 2 ) :: ocean_flow
   ! Erosion maximum per time step
   real( tkind ), dimension( : ), allocatable :: eromax

   !> Stream classification type
   ! Number of stream classes number
   integer :: stcl_nb
   type :: stream_class
      ! Stream slope
      real( tkind ) :: slope
      ! Depth channel
      real( tkind ) :: cdepth
      ! Width:Depth ratio
      real( tkind ) :: wd_ratio
  end type stream_class
  type( stream_class ),dimension(:),allocatable :: streamclass

  !> flow walkers parameter type
  ! FW global reference
  integer,dimension(:),allocatable :: fa_refid
  ! Number of point in the area of influence for ero/dep rules
  integer,dimension(:),allocatable :: fa_nbPtAI
  ! FW inside the simulation domain(0), if outside (1), if forced deposition (2)
  integer,dimension(:),allocatable :: fa_inbound
  ! FW flow type
  integer,dimension(:),allocatable :: fa_type
  ! FW TIN step counter
  integer,dimension(:),allocatable :: fa_count
  ! FW x position
  real( tkind ),dimension(:),allocatable :: fa_xpos
  ! FW y position
  real( tkind ),dimension(:),allocatable :: fa_ypos
  ! FW z position
  real( tkind ),dimension(:),allocatable :: fa_zpos
  ! X slope experienced by FW
  real( tkind ),dimension(:),allocatable :: fa_xslp
  ! Y slope experienced by FW
  real( tkind ),dimension(:),allocatable :: fa_yslp
  ! FW height
  real( tkind ),dimension(:),allocatable :: fa_h
  ! Stream width
  real( tkind ),dimension(:),allocatable :: fa_width
  ! FW x velocity
  real( tkind ),dimension(:),allocatable :: fa_xvel
  ! FW y velocity
  real( tkind ),dimension(:),allocatable :: fa_yvel
  ! FW x velocity solution for the ODE
  real( tkind ),dimension(:),allocatable :: fa_xvelO
  ! FW y velocity solution for the ODE
  real( tkind ),dimension(:),allocatable :: fa_yvelO
  ! FW flow volume
  real( tkind ),dimension(:),allocatable :: fa_volume
  ! FW flow discharge (m3/s)
  real( tkind ),dimension(:),allocatable :: fa_discharge
  ! FW stream density
  real( tkind ),dimension(:),allocatable :: fa_density
  ! FW previous gradient value
  real( tkind ),dimension(:),allocatable :: fa_pgrad
  ! C1 coefficient of bottom friction
  real( tkind ),dimension(:),allocatable :: fa_Cfric
  ! Flow walker acceleration along X direction
  real( tkind ),dimension(:),allocatable :: fa_accx
  ! Flow walker acceleration along Y direction
  real( tkind ),dimension(:),allocatable :: fa_accy
  ! Zone of flow establishment for the plume
  real( tkind ),dimension(:),allocatable :: fa_pzfe
  ! Plume coefficient of decceleration
  real( tkind ),dimension(:),allocatable :: fa_pcoeff
  ! Plume trajectory distance from the mouth
  real( tkind ),dimension(:),allocatable :: fa_pdist
  ! Plume initial x velocity at the mouth
  real( tkind ),dimension(:),allocatable :: fa_pvelx
  ! Plume initial y velocity at the mouth
  real( tkind ),dimension(:),allocatable :: fa_pvely
  ! FW all previous face ID on the TIN grid
  integer,dimension(:,:),allocatable :: fa_faceID
  ! FW GIDs point in the area of influence for ero/dep rules
  integer,dimension(:,:),allocatable :: fa_ptsAIfw
  ! FW sediment discharge (m3/s)
  real( tkind ),dimension(:,:),allocatable :: fa_sedcharge

  !> Flow walker information type
  ! FW flow type
  integer,dimension(:),allocatable :: recfw_type
  ! FW x position
  real( tkind ),dimension(:),allocatable :: recfw_xpos
  ! FW y position
  real( tkind ),dimension(:),allocatable :: recfw_ypos
  ! FW z position
  real( tkind ),dimension(:),allocatable :: recfw_zpos
  ! FW x velocity
  real( tkind ),dimension(:),allocatable :: recfw_xvel
  ! FW y velocity
  real( tkind ),dimension(:),allocatable :: recfw_yvel
  ! FW source height
  real( tkind ),dimension(:),allocatable :: recfw_srch
  ! FW source width
  real( tkind ),dimension(:),allocatable :: recfw_width
  ! Flow discharge
  real( tkind ),dimension(:),allocatable :: recfw_volume
  ! FW sediment charge
  real( tkind ),dimension(:),allocatable :: recfw_sedcharge

   real( tkind ), parameter :: maxdepth=0.50_8
   real( tkind ), parameter :: maxdepinv=2.0_8

   ! FW local node id
   integer,dimension(:),allocatable :: ed_refid
   ! FW local processor id
   integer,dimension(:),allocatable :: ed_procid
   ! FW inside the simulation domain(0), if outside (1)
   integer,dimension(:),allocatable :: ed_bound
   ! FW flow type
   integer,dimension(:),allocatable :: ed_type
   ! Number of points in the area of influence of the FW
   integer,dimension(:),allocatable :: ed_AInb
   ! GIDs of points in the area of influence of the FW
   integer,dimension(:,:),allocatable :: ed_AIpts
   ! FW z position
   real( tkind ),dimension(:),allocatable :: ed_zpos
   ! FW height
   real( tkind ),dimension(:),allocatable :: ed_h
   ! FW velocity
   real( tkind ),dimension(:),allocatable :: ed_vel
   ! FW volume
   real( tkind ),dimension(:),allocatable :: ed_vol
   ! FW stream density
   real( tkind ),dimension(:),allocatable :: ed_dens
   ! FW sediment discharge in stream
   real( tkind ),dimension(:,:),allocatable :: ed_sed

contains

   ! ============================================================================
   !> Function manning_fct
   !! Function manning sets manning's coefficient.
   !! When flow walker's location is in the sea and the density of fresh water
   !! including suspended sediment is less than density of sea water, hypopycnal
   !! flow occurs, and manning's coefficient for hypopycnal flow is applied.
   !! when flow walker is located in the sea and the density of fresh water,
   !! including suspended sediment, is greater than density of sea water, then
   !! hyperpycnal flow occurs, and manning's coefficient for hyperpycnal flow is
   !! employed.
   !! \param depth, den
   !<
   ! ============================================================================
   function manning_fct(depth, den) result( manning )

      real( tkind ) :: manning, depth, den

      if( depth > 0.0_8 )then
         if( den > sea_density )then
            if( depth > maxdepth )then
               manning = manning_hyper
            else
               manning = manning_open - ( manning_open - manning_hyper )* depth * maxdepinv
            endif
         elseif( den <= sea_density )then
            manning = manning_hypo
         endif
      else
         manning = manning_open
      endif

      return

   end function manning_fct
   ! ============================================================================

end module fwalker_data
