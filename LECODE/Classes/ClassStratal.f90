! ============================================================================
! Name        : ClassStratal.f90
! Author      : tristan salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file ClassStratal.f90
!
! Description :  It encapsulates the sediment and stratal module used by SPModel during compilation.
!
!<
! ============================================================================
!> Module sediment_data encapsulates characteristics of the sediments used
!<
module sediment_data

    use precision_data

    integer, parameter :: max_grn = 10

    ! Number of siliciclastics
    integer :: silgrn
    ! Number of carbonates
    integer :: carbgrn
    ! Number of organics
    integer :: orggrn
    ! Number of total grains
    integer :: totgrn
    ! Hemipelagic material number
    integer :: hemi_mat
    ! Minimum slope value (dz/dx)
    real( tkind ) :: minimum_slp
    ! Order of grain sizes deposition
    integer,dimension(:),allocatable :: dep_order
    ! Material name
    character(len=128), dimension( : ), allocatable :: material_name

    ! Carbonate/organic membership functions number
    integer :: membership_fct_nb
    ! Carbonate/organic fuzzy functions number
    integer :: fuzzy_fct_nb
    ! Carbonate/organic fuzzy rules number
    integer :: fuzzy_rl_nb

    ! Carbonates/organics membership function declaration
    type membership_fct
        ! Membership function name
        character(len=128) :: name
        ! Defined membership variable ID
        integer :: variableID
        ! Sediment class ID
        integer :: sedclassID
        ! Number of points defining the membership function
        integer :: nb_points
        ! Coordinates of each points
        real( tkind ), dimension( 20, 2 ) :: xypoints
    end type membership_fct
    type(membership_fct),dimension(:),allocatable :: membership

    ! Carbonates/organics fuzzy function declaration
    type fuzzy_fct
        ! Fuzzy function name
        character(len=128) :: name
        ! Defined fuzzy variable ID
        integer :: variableID
        ! Number of points defining the fuzzy function
        integer :: nb_points
        ! Coordinates of each points
        real( tkind ), dimension( 20, 2 ) :: xypoints
    end type fuzzy_fct
    type(fuzzy_fct),dimension(:),allocatable :: fuzzyfct

    ! Carbonates/organics fuzzy rule declaration
    type fuzzy_rule
        ! Sediment class ID
        integer :: sedclassID
        ! Fuzzy function name ID
        integer :: fuzzyID
        ! Number of membership function to apply
        integer :: membershipNb
        !  ID of applied membership function
        integer, dimension( 20 ) :: membershipID
    end type fuzzy_rule
    type(fuzzy_rule),dimension(:),allocatable :: fuzzyfy

    ! Combined distribution declaration
    type fuzzy_combined
        ! Number of points defining the fuzzy output distribution
        integer :: nb_points
        ! Coordinates of each points
        real( tkind ), dimension( 40, 2 ) :: xypoints
    end type fuzzy_combined
    type(fuzzy_combined),dimension(:,:),allocatable :: fuzzycomb

    !> Sediment parameters type
    type sediment_parameters
        ! Transport mode suspension (1) or bedload (0)
        integer :: transport
        ! Diameter in millimetres
        real( tkind ) :: diameter
        ! Density in kg.m-3
        real( tkind ) :: density
        ! Fall velocity
        real( tkind ) :: vfall
        ! Marine stability slope
        real( tkind ) :: slp_marine
        ! Aerial stability slope
        real( tkind ) :: slp_aerial
    end type sediment_parameters
    type(sediment_parameters),dimension(:),allocatable :: sediment

    ! Carbonates/organics evolution flags
    integer :: mat_nearest_flag, carb_shore_flag, carb_river_flag, carb_ocean_flag
    integer :: carb_sediment_flag, carb_slope_flag, carb_burial_flag, salinity_flag, temperature_flag
    integer :: carb_valley_flag, carb_depth_flag, carb_age_flag, carb_exposed_flag

    ! Sum of sediment deposited directly since the last call (m)
    real( tkind ), dimension( : ), allocatable :: carb_sedimentation

end module sediment_data
! ============================================================================
!> Module strata_data encapsulates characteristics of the stratigraphic layers
!<
module strata_data

    use precision_data
    use kdtree2_module
    use kdtree2_precision_module

    implicit none

    public

    logical :: seismic_plug, rms_plug

    ! Initial deposit layer number
    integer :: InitDep
    ! Number of points on the initial cartesian grid
    integer :: strat_X, strat_Y
    ! Xo & Yo coordinates of the SW corner of the initial TIN
    real( tkind ) :: strat_xo, strat_yo
    ! Grid spacing on the borders
    real( tkind ) :: strat_dx
    ! Fine grid spacing for horizontal deformation
    real( tkind ) :: fine_dx
    ! Grid partitioning values
    integer,dimension(:),allocatable :: zolt_part
    ! Output integer
    integer, dimension( : ),allocatable :: check_int
    ! Record elevation
    real( tkind ),dimension(:),allocatable :: elev_record
    ! Record elevation per processor
    real( tkind ),dimension(:),allocatable :: proc_elev
    ! Record halo border elevation per processor
    real( tkind ),dimension(:),allocatable :: halo1_elev
    ! Record halo border elevation per processor
    real( tkind ),dimension(:),allocatable :: halo2_elev
    ! Record uplift
    real( tkind ),dimension(:),allocatable :: uplift

    ! Diffusion grid partitioning values
    integer :: GlobalNodes
    integer,dimension(:),allocatable :: NodesLocal
    integer,dimension(:),allocatable :: DispLocal

    ! Nodes to be used for creating the TIN based on stratal grid
    integer, dimension(:),allocatable :: pick_refine
    ! Boundaries of the simulated model
    integer,dimension( 4 ) :: bounds

    ! Seismic lines
    real( tkind ) :: seismicX_min, seismicX_max
    real( tkind ) :: seismicY_min, seismicY_max

    ! RMS geocellular extension
    real( tkind ) :: RMSX_min, RMSX_max
    real( tkind ) :: RMSY_min, RMSY_max

    !--------------------------------------------------------------
    ! Local stratigraphic grid declaration

    !> Deposit definition

    ! Number of columns on each partition
    integer, dimension( : ), allocatable :: strat_col
    ! Index of deposited layer
    integer,dimension(:,:),allocatable :: strat_layID
    ! Hardness of the deposit
    real( tkind ),dimension(:,:),allocatable :: strat_hardness
    ! Stratal elevation
    real( tkind ),dimension(:,:),allocatable :: strat_zelev
    ! Stratal thickness
    real( tkind ),dimension(:,:),allocatable :: strat_thick
    ! Stratal facies
    real( tkind ),dimension(:,:),allocatable :: strat_fac
    ! Stratal sediment thickness
    real( tkind ),dimension(:,:,:),allocatable :: strat_sedh
    ! Stratal porosity
    real( tkind ),dimension(:,:,:),allocatable :: strat_porosity

    !> Local node type
    ! Number of layers deposited on the considered vertex
    integer,dimension(:),allocatable :: locv_NbLay
    ! Global ID of the considered vertex
    integer,dimension(:),allocatable :: locv_gid
    ! Displacements for the considered node
    real( tkind ),dimension(:),allocatable :: locv_vdisp
    ! Rate of displacements for the considered node
    real( tkind ),dimension(:),allocatable :: locv_vrate
    ! XYZ position of the base of the considered vertex stratigraphic layer
    real( tkind ),dimension(:,:),allocatable :: locv_coord

    !> Local stratigraphic face
    ! Global ID of the considered face
    integer,dimension(:),allocatable :: locf_gid
    ! Stratal face neighbors
    integer,dimension(:,:),allocatable :: locf_points
    ! Coordinate of the centre of the 2d face
    real( tkind ),dimension(:,:),allocatable :: locf_XYcoord

    !> Local stratigraphic grid
    type loc_grid
        ! Current layer index
        integer :: layersID
        ! Total number of nodes locally
        integer :: nodes_nb
        ! Total number of faces locally
        integer :: faces_nb
    end type loc_grid
    type( loc_grid ) :: L_strat

    !--------------------------------------------------------------
    ! Global stratigraphic grid declaration

    !> Global stratigraphic grid
    type strat_grid
        ! Total number of layers for the entire simulation
        integer :: total_layer
        ! Total number of nodes in the stratigraphic grid
        integer :: nodes_nb
        ! Total number of faces in the stratigraphic grid
        integer :: faces_nb
    end type strat_grid
    type( strat_grid ) :: G_strat

    !> Global node type
    ! Number of processors sharing the node (max = 2)
    integer,dimension( : ),allocatable :: gv_halo
    ! Processor ID sharing the node
    integer,dimension( :,: ),allocatable :: gv_ShareID
    ! Local node ID on each processor sharing the node
    integer,dimension( :,: ),allocatable :: gv_SharenID
    ! Neighbors ID
    integer,dimension( :,: ),allocatable :: gv_ngbID
    ! Displacements for the considered node
    real( tkind ),dimension( : ),allocatable :: gv_vdisp
    ! Node XYZ position
    real( tkind ),dimension( :,: ),allocatable :: gv_coord
    ! Coastline orientation in radian
    real( tkind ),dimension( : ),allocatable :: gv_orientation

    !> Global stratigraphic face
    ! Processor ID containing the face
    integer,dimension( : ),allocatable :: gf_pid
    ! Stratal face nodes
    integer,dimension( :,: ),allocatable :: gf_points
    ! Coordinate of the centre of the 2d face
    real( tkind ),dimension( :,: ),allocatable :: gf_XYcoord

    !> Global halo border grid points
    ! Halo neighbor proc id
    integer :: halo_ngb( 2 )
    ! Number of halo border points
    integer :: haloborder_nb
    ! Number of halo points on processor
    integer :: haloproc_nb
    ! Processor ID of the halo node
    integer,dimension( : ), allocatable :: halo_proc
    ! Local ID of the halo node
    integer,dimension( : ), allocatable :: halo_lid
    ! Global ID of the halo node
    integer,dimension( :,: ), allocatable :: halo_gid
    ! Processor where we need to send the elevation
    integer,dimension( : ), allocatable :: halo_send
    ! Global ID of the halo node send from left to the processor
    integer,dimension( : ), allocatable :: halo_rcv2
    ! Global ID of the halo node send from right to the processor
    integer,dimension( : ), allocatable :: halo_rcv1

     !--------------------------------------------------------------
     ! Hot start stratigraphic grid declaration

    !> Stratigraphy type
    ! Index of deposited layer
    integer,dimension( :,: ), allocatable :: rec_slayID
    ! Hardness of the deposit
    real( tkind ),dimension( :,: ), allocatable :: rec_shardness
    ! Stratal elevation
    real( tkind ),dimension( :,: ), allocatable :: rec_szelev
    ! Stratal porosity
    real( tkind ),dimension( :,:,: ), allocatable :: rec_sporosity
    ! Stratal thickness
    real( tkind ),dimension( :,: ), allocatable :: rec_sthick
    ! Stratal facies
    real( tkind ),dimension( :,: ), allocatable :: rec_sfac
    ! Stratal sediment thickness
    real( tkind ), dimension( :,:,: ), allocatable :: rec_ssedh

    ! Number of layers deposited on the considered vertex
    integer,dimension( : ), allocatable :: rec_NbLay
    ! Top elevation of the stratigraphic stack
    real( tkind ),dimension( : ), allocatable :: rec_topelev
    ! XYZ position of the base of the considered vertex stratigraphic layer
    real( tkind ),dimension( :,: ), allocatable :: rec_coord
    ! Flow accumulation and no sink topographic dem
    real( tkind ), dimension( : ), allocatable :: wdem, facc

    ! Mass movement parameters
    ! Mass movement flag
    logical :: massmove
    ! Curvature threshold for mass movement initiation
    real( tkind ) :: mass_curv
    ! Flow accumulation threshold for mass movement initiation
    integer :: mass_acc
    type massmv
        ! Surface upslope area ( m2 )
        real( tkind ) :: upslope
        ! Steepness index
        real( tkind ), dimension( 2 ) :: steepness
        ! Concavity index
        real( tkind ), dimension( 2 ) :: concavity
    end type massmv
    type(massmv),dimension( 2 ) :: mass_mov

end module strata_data
! ============================================================================
module mod_diffuse

    use precision_data

    implicit none

    integer, parameter :: diff_nb = 2
    integer, parameter :: max_it_cyc = 10000

    real( tkind ), parameter :: diff_res = 1.e-4_8
    real( tkind ), parameter :: toplimit = 1.e10_8


    ! Number of point in the local grid
    integer :: DL_nodesNb
    ! Columns number in the current local grid
    integer :: diff_cols
    ! Processor ID for bottom and top neighbors respectively
    integer, dimension( 2 ) :: diff_neighbors
    ! Local ID in the stratigraphic grid
    integer, dimension( : ), allocatable :: diff_ID
    ! Define reference node within the local stratal grid
    integer, dimension( : ), allocatable :: diff_ref
    ! Vertex information for the diffusion grid
    real( tkind ), dimension( :,: ), allocatable :: diff_coord
    ! Neighbors ID information for the diffusion grid
    integer, dimension( :,: ), allocatable :: diff_ngb
    ! Stratal sediment composition
    real( tkind ), dimension(:,:),allocatable :: top_sedh
    ! Previous stratal sediment composition
    real( tkind ), dimension(:,:),allocatable :: top_sedprev

    real( tkind ), dimension(:,:), allocatable :: depo
    real( tkind ), dimension(:,:), allocatable :: dstart
    real( tkind ), dimension( : ), allocatable :: topmax
    real( tkind ), dimension( : ), allocatable :: cdif
    real( tkind ), dimension( : ), allocatable :: difo
    real( tkind ), dimension( : ), allocatable :: difp


end module mod_diffuse
! ============================================================================
