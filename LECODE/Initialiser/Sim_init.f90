! ============================================================================
! Name        : Simu_init.f90
! Author      : Tristan Salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file  Simu_init.f90
!!
!! Description : initialises the parameters used in the experiment based on input file
!!
!<
! ============================================================================
module mod_simulation

    use file_data
    use mpidata
    use time_data
    use TIN_out
    use ocean_io
    use ocean_data
    use error_data
    use strata_out
    use strata_data
    use TIN_data
    use forces_data
    use param_data
    use precision_data
    use mod_compaction
    use mod_displacements

    use strata_ini
    use params_init
    use check_point
    use rain_initial
    use TIN_surface
    use TIN_function

    use FWLK_parser
    use TIN_parser
    use Flux_parser
    use strata_parser
    use Forces_parser
    use carborg_parser
    use simulation_parser
    use mod_filldepression

    use mod_hemipelagic

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine spm_initial
    !! Subroutine spm_initial initializes the SPModel model.
    !<
    ! ============================================================================
    subroutine spm_initial

        integer :: iunit

        ! Initialise some flags
        rcdin = 0
        step = 0
        stepin = 0
        stepc = 0
        stepout = 0
        attempt = 0
        fexistx = .false.
        allocate( disp( nproc ) )

        ! Read simulation parameters from XmL
        call simulation_initialisation

        ! Define grids
        call simulation_grids_declaration

        ! Initialisation of processes variables
        nextevent = SPM_INIT
        tnow = time_start
        if( hemi_flag ) last_hemi = tnow

        ! Create diffusion grid
        call create_diffusion_grid

        ! Initialise porosity
        if( gporo%compaction .and. .not. restart ) call porosity_init

        L_strat%layersID = L_strat%layersID + 1
        rcd = 0
        if( restart ) stepout = restart_iter
        if( restart ) stepc = restart_iter

        ! Define simulation time step for each processes
        if( tnow == time_start .and. next_display == 0.0_8 )then
            time_next = flow_int + tnow
            next_display = next_display + tnow
            next_output = next_output + tnow
            next_coarse_dt = coarse_dt + tnow
            next_slpdiff = slope_int + tnow
            next_rain = time_end + 1.e20_8
            if( geomorphology ) next_rain = rain_int + tnow
            if( gdisp%event  > 0 )then
                next_displacement = next_displacement + tnow
            else
                next_displacement = time_end + 1.e20_8
            endif
        endif

        ! Create simulation initial output
        tnow = next_display
        next_display = time_display + next_display
        call xdmf_output( stepout )

        next_output = next_output + time_output
        stepout = stepout + 1
        num_fw_o = 0
        rcd = 0
        rcdin = 0
        stepin = stepin + 1
        stepc = stepc + 1
        fexistx = .false.
        flagx = .true.

        ! Delete flatten XmL file from the local directory
        if( iam == 0 )then
            iunit = 423
            open( iunit, file=finput, form='unformatted' )
            close( iunit, status='delete' )
        endif

        return

    end subroutine spm_initial
    ! ============================================================================
    !> Subroutine simulation_initialisation
    !! Allocates all parameters that defines an experiment
    !<
    ! ============================================================================
    subroutine simulation_initialisation

        ! Parameters Declaration
        logical :: found
        integer :: l1, l2, i, j, spare1, spare2

        real( tkind ) :: ratcheck

        character(len=128) :: command, stg, fildir

        ! Default output directory
        outdir=''
        outdir(1:7)='outputs'
        outputs=''
        runfiles=''
        outputs(1:7)='outputs'
        runfiles(1:8)='runfiles'
        InitDep = 0
        gdisp%event = 0
        udw_coupling = 0
        gsea%sealevel = .false.
        gsea%output = .true.
        gporo%compaction = .false.
        new_disp = .false.
        geomorphology = .false.
        hemi_flag = .false.
        rain_region = .true.
        restart = .false.
        udw = .false.
        deformlayer = .false.
        plume_cmpt = .false.

        ! Default densities
        fluid_density = 1015.0_8
        sea_density = 1027.0_8

        ! Default Manning coefficient
        manning_hypo = 0.070_8
        manning_hyper = 0.020_8
        manning_open = 0.025_8
        transport%erosion_limit = 100.0_8

        ! Default rain parameters
        rain_elev = 0.0_8
        rain_nb = 1
        fine_dx = 0.0_8
        num_src = 0
        bounds(1:4) = 0

        ! Parsing input file
        rain_int = 1.e6
        flow_int = 1.e6
        if( iam == 0 )then
            call xml_Strata
            call xml_TIN
            call xml_Forces
            call xml_Sim
            call xml_FWLK
            call xml_Flux
        endif

        ! Broadcast values globally
        call broadcast_strata
        call broadcast_sim
        call broadcast_tin
        call broadcast_forces
        call broadcast_fwlk
        call broadcast_flux

        ! Time section
        if( timesection )then
            if( sampling_int <= 0 ) attempt = TIME_SPI
            if( time_output <= 0 ) attempt = TIME_DISP
            if( time_display <= 0 ) attempt = TIME_DISP
            if( time_start >= time_end ) attempt = TIME_DEF
            if( sampling_int > time_display ) sampling_int = time_display
        else
            attempt = SEC_TIME
        endif

        ! Ensure that all processes have successfully read the initial parameters.
        call completion

        ! Carbonates/organics parameters
        if( silgrn < totgrn )then
            if( iam == 0 ) call xml_CarbOrgs
            call broadcast_carborgs
        endif
        if( sampling_int < slope_int ) sampling_int = slope_int

        ! Hemipelagic section
        hemi_mat = 0
        if( hemi_flag )then
            do i = 1, totgrn
                stg = material_name( i )
                if( stg( 1:4 ) == 'Hemi' .or. stg( 1:4 ) == 'hemi' ) hemi_mat = i
            enddo
            if( hemi_mat == 0 )then
                 if( iam == 0 ) &
                    print*,'WARNING: To use hemipelagic sedimentation it is mandatory to named Hemi one of the materials'
                 stop
            endif
            ! Read hemipelagic file
            call read_hemilpelagic_file
        endif

        ! Sources section
        do i = 1,num_src
            ratcheck = 0.0_8
            do j = 1, silgrn
                ratcheck=ratcheck+fws_sedperc( i , j )
            enddo
            if( ratcheck /= 100.0_8 ) attempt = SRCPER
            if( fws_tstrt( i ) > fws_tend( i ) ) attempt = SRCTIME
            spare1 = int((fws_xposition( i ) - strat_xo)/strat_dx + 1 )
            spare2 = int((fws_yposition( i ) - strat_yo)/strat_dx + 1 )
            if( spare1 < 1 .or. spare1 > strat_X ) attempt = SRCXPOS
            if( spare2 < 1 .or. spare2 > strat_Y ) attempt = SRCYPOS
        enddo

        ! Ensure the proper declaration of rivers.
        call completion

        ! Create the output directories
        if( iam == 0 )then
            call noblnk(outdir)
            fildir = outdir
            stg = '/runfiles/TIN.node'
            call noblnk(stg)
            call append_str( fildir,stg )
            call noblnk(fildir)
            inquire(file=fildir,exist=found)
            if( .not. found )then
                command=' '
                command(1:6)='mkdir '
                l1 = len_trim(outdir)
                command(7:l1+7)=outdir
                call term_command( command )
                command(l1+7:l1+7)='/'
                stg=''
                stg(1:l1+7)=command(1:l1+7)
                stg(l1+8:l1+15)=outputs(1:7)
                call term_command( stg )
                stg=''
                stg(1:l1+7)=command(1:l1+7)
                stg(l1+8:l1+16)=runfiles(1:8)
                call term_command( stg )
                outdir1(1:l1)=outdir(1:l1)
                outdir1(l1+1:l1+1)='/'
                outdir1(l1+2:l1+9)=outputs(1:7)
                outdir2(1:l1)=outdir(1:l1)
                outdir2(l1+1:l1+1)='/'
                outdir2(l1+2:l1+10)=runfiles(1:8)
            else
                command=' '
                command(1:6) = 'rm -r '
                l1 = len_trim(outdir)
                command(7:l1+7)=outdir
                call term_command( command )
                command=' '
                command(1:6)='mkdir '
                l1 = len_trim(outdir)
                command(7:l1+7)=outdir
                call term_command( command )
                command(l1+7:l1+7)='/'
                stg=''
                stg(1:l1+7)=command(1:l1+7)
                stg(l1+8:l1+15)=outputs(1:7)
                call term_command( stg )
                stg=''
                stg(1:l1+7)=command(1:l1+7)
                stg(l1+8:l1+16)=runfiles(1:8)
                call term_command( stg )
                outdir1(1:l1)=outdir(1:l1)
                outdir1(l1+1:l1+1)='/'
                outdir1(l1+2:l1+9)=outputs(1:7)
                outdir2(1:l1)=outdir(1:l1)
                outdir2(l1+1:l1+1)='/'
                outdir2(l1+2:l1+10)=runfiles(1:8)
            endif

            ! Underworld Lecode shared folder
            if( udw_plug )then
                call noblnk(outdir3)
                fildir = outdir3
                stg = '/topsurface.vtk'
                call noblnk(stg)
                call append_str( fildir,stg )
                call noblnk(fildir)
                inquire(file=fildir,exist=found)
                if( .not. found )then
                    command=' '
                    command(1:6)='mkdir '
                    l1 = len_trim(outdir3)
                    command(7:l1+7)=outdir3
                    call term_command( command )
                else
                    command=' '
                    command(1:6) = 'rm -r '
                    l1 = len_trim(outdir3)
                    command(7:l1+7)=outdir3
                    call term_command( command )
                    command=' '
                    command(1:6)='mkdir '
                    l1 = len_trim(outdir3)
                    command(7:l1+7)=outdir3
                    call term_command( command )
                endif

            endif

            ! Copy XML file in the output folder
            command = ''
            command(1:3) = 'cp '
            l1 = len_trim(finput)
            command(4:l1+4) = finput(1:l1)
            command(l1+4:l1+4)= ' '
            l2 = len_trim(outdir)
            command(l1+5:l1+5+l2) = outdir(1:l2)
            command(l1+5+l2:l1+5+l2) = '/'
            call term_command( command )

        endif

        ! Broadcast file name
        call mpi_bcast( outdir,128,mpi_character,0,SPModel_comm_world,ierr )
        call mpi_bcast( outdir1,128,mpi_character,0,SPModel_comm_world,ierr )
        call mpi_bcast( outdir2,128,mpi_character,0,SPModel_comm_world,ierr )
        call mpi_bcast( outdir3,128,mpi_character,0,SPModel_comm_world,ierr )

        ! Allocate porosity table
        if( gporo%compaction )then
            allocate( effPressure( gporo%ePnb ) )
            allocate( porosity( totgrn, gporo%ePnb ) )
            call assign_porosity_table
        endif

        ! Get the Ocean parameters if the coupling exists
        if( iam == 0 ) call xml_Ocean
        call broadcast_ocean

        return

    end subroutine simulation_initialisation
    ! ============================================================================
    !> Subroutine simulation_grids_declaration
    !! Allocates all mesh and grid that define an experiment
    !<
    ! ============================================================================
    subroutine simulation_grids_declaration

        ! Generate the stratigraphic grid globally
        call generate_global_stratigraphy

        ! In case of restart use previous parameters
        if( restart ) call define_restarted_simulation

        ! Define band partitioning
        allocate( strat_col( nproc ) )
        call define_band_partitioning

        ! Define the stratigraphic architecture locally
        call generate_local_stratigraphy

        ! Flow accumulation files declaration
        ! NOTE: this needs to be pass to a CPP function
        fdem = 'demfill.asc'
        call addpath2( fdem )
        fdemcpp = ''
        fdemcpp = fdem
        call noblnk( fdemcpp )
        fdemcpp( len( fdemcpp ) : len( fdemcpp ) ) = CHAR(0)
        faccu = 'FA.asc'
        call addpath2( faccu )
        faccucpp = ''
        faccucpp = faccu
        call noblnk( faccucpp )
        faccucpp( len( faccucpp ) : len( faccucpp ) ) = CHAR(0)

        ! Arrays declaration for flow accumulation computation
        if( allocated( wdem ) ) deallocate( wdem )
        allocate( wdem( G_strat%nodes_nb ) )
        if( allocated( facc ) ) deallocate( facc )
        allocate( facc( G_strat%nodes_nb ) )

        ! Allocate TIN grid arrays
        if( allocated( tinv_pid ) ) deallocate( tinv_pid )
        if( allocated( ngbID ) ) deallocate( ngbID )
        if( allocated( tinv_ngb ) ) deallocate( tinv_ngb )
        if( allocated( tinv_coord ) ) deallocate( tinv_coord )
        if( allocated( tinv_boundary ) ) deallocate( tinv_boundary )

        allocate( ngbID( G_strat%nodes_nb, mxnghb ) )
        allocate( tinv_ngb( G_strat%nodes_nb ) )
        allocate( tinv_coord( G_strat%nodes_nb, 3 ) )
        allocate( tinv_boundary( G_strat%nodes_nb ) )
        allocate( tinv_pid( G_strat%nodes_nb ) )

        if( allocated( tinf_boundary ) ) deallocate( tinf_boundary )
        if( allocated( tinf_nghb ) ) deallocate( tinf_nghb )
        if( allocated( tinf_pid ) ) deallocate( tinf_pid )
        if( allocated( tinf_nids ) ) deallocate( tinf_nids )
        if( allocated( tinf_centroid ) ) deallocate( tinf_centroid )
        if( allocated( tinf_length ) ) deallocate( tinf_length )

        allocate( tinf_boundary( 2*G_strat%faces_nb ) )
        allocate( tinf_nghb( 2*G_strat%faces_nb ) )
        allocate( tinf_pid( 2*G_strat%faces_nb, 3 ) )
        allocate( tinf_nids( 2*G_strat%faces_nb, 3 ) )
        allocate( tinf_centroid( 2*G_strat%faces_nb, 3 ) )
        allocate( tinf_length( 2*G_strat%faces_nb, 2 ) )
        allocate( fids( 2*G_strat%faces_nb, mxnghb ) )

        ! Allocate erosion array
        allocate( eromax( G_strat%nodes_nb )  )

        ! Create the initial TIN
        call generate_TIN_grid

        ! Define additional parameters
        call spm_parameters_initialise

        return

    end subroutine simulation_grids_declaration
    ! ============================================================================
    !> Subroutine define_band_partitioning
    !! Subroutine define_band_partitioning used to send the partitioned grid values.
    !<
    ! ============================================================================
    subroutine define_band_partitioning

        integer :: i, avecol, extras, cols, k

        real( tkind ) :: maxP( nproc ), ptmax

        ! Check higher cell number in X or Y direction
        if( strat_Y >= strat_X )then
            avecol = int( ( strat_Y - 1 )/ nproc )
            extras = mod( ( strat_Y - 1 ), nproc )
        else
            avecol = int( ( strat_X - 1 )/ nproc )
            extras = mod( ( strat_X - 1 ), nproc )
        endif

        if( iam == 0 )then

            ! Get the number of columns for each processor
            do k = 0, nproc - 1
                cols = avecol
                if( k < extras ) cols = avecol + 1
                if( k == 0 )then
                    if( strat_Y >= strat_X )then
                        maxP( k + 1 ) = strat_yo + cols * strat_dx
                    else
                        maxP( k + 1 ) = strat_xo + cols * strat_dx
                    endif
                else
                    maxP( k + 1 ) = maxP( k ) + cols * strat_dx
                endif
                strat_col( k+1 ) = cols + 1
            enddo

            ! Loop over the faces to find the processor they belong to.
            do i = 1, G_strat%faces_nb
                partition_face: do k = 1, nproc
                    ptmax = gf_XYcoord( i, 1 )
                    if( strat_Y >= strat_X ) ptmax = gf_XYcoord( i, 2 )
                    if( ptmax < maxP( k ) )then
                        gf_pid( i ) = k - 1
                        exit partition_face
                    endif
                enddo partition_face
            enddo
        endif

        ! Broadcast values globally
        call mpi_bcast(strat_col, nproc, int_type, 0, SPModel_comm_world, ierr)
        call mpi_bcast(gf_pid, G_strat%faces_nb, int_type, 0, SPModel_comm_world, ierr)

        return

    end subroutine define_band_partitioning
    ! ============================================================================

end module mod_simulation
