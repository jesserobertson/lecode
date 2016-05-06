! ============================================================================
! Name        : FunctionIdle.f90
! Author      : tristan salles
! Created on: Aug 17, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file idle_function.f90
!!
!! idle_function is called when wathever main subroutines stopped. It
!! then determines which subroutine will be call next.
!!
!<
! ============================================================================
module mod_idle

    use flux_data
    use mpidata
    use time_data
    use fwalker_data
    use strata_data
    use mod_recordfw
    use mod_cleanfw
    use mod_inflow
    use rain_initial
    use mod_inflow
    use strata_out
    use strata_ini
    use param_data
    use forces_data
    use mod_slopefw
    use TIN_function
    use TIN_surface
    use mod_sealevel
    use mod_direction
    use mod_layersFW
    use mod_rungekutta
    use mod_compaction
    use mod_displacements
    use check_point
    use sediment_data
    use mod_seddiffusion
    use mod_filldepression
    use mod_advancefw
    use UDWplug1
    use Seismic
    use RMSim
    use ocean_data
    use Oceanplug
    use ocean_transport
    use mod_hemipelagic
    use kdtree2_module
    use mod_edtransfer
    use mod_buildcarbs
    use kdtree2_precision_module

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine spm_run
    !! Subroutine spm_run computes surface processes and sediment transport.
    !<
    ! ============================================================================
    subroutine spm_run

        integer :: endflag, n

        ! Update flow walkers
        if( tnow == time_start )then

            ! Build FWs mpi type
            call build_mpi_flowconc_type
            call build_mpi_flowwalker_type
            call build_mpi_flowtransfer_type

            ! Find river flow walkers
            if( num_src > 0 ) call define_river_fws

            ! Find rain flow walkers
            if( geomorphology ) call register_rain_arrays

            ! Define FW's global IDs
            call define_flow_walkers_global_ids

            ! Find global paths
            call find_global_pathways_G8_Paik
            oldevent = SPM_DISPLAY

            ! Find next process
            call control_process_events

            ! Carbonates/organics declaration of initial elevation
            if( silgrn < totgrn )then
                if( allocated( carb_sedimentation ) ) deallocate( carb_sedimentation )
                allocate( carb_sedimentation( L_strat%nodes_nb ) )
                do n = 1, L_strat%nodes_nb
                    carb_sedimentation( n ) = elev_record( locv_gid( n ) )
                enddo
                next_carborg = tnow
            else
                next_carborg = time_end + 1.e20_8
            endif
        endif

        ! Allocate wave current circulation fields for Ocean plugin
        if( tnow == time_start .and. ocean_plug )then
            call create_oceanplug_surface
            call ocean_plugin_wait_function
            next_oceancirc = tnow
        else
            next_oceancirc = time_end + 1.e20_8
        endif

        ! Allocate displacement fields for UDW plugin
        if( tnow == time_start .and. udw_plug )then
            gdisp%disptime = 0.0_8
            call create_vtk_topsurface( 0 )
            call udw_plugin_wait_function
            new_disp = .true.
            ! Assign displacement rate
            call assign_vertical_displacements_rate( 1 )
            ! Update displacement
            call update_local_displacements
        endif

        endflag = 0
        ! Find event to perform
        event_control: do

            ! If we've reach the end of the simulation
            if( nextevent == SPM_QUIT )then
                ! If seismic line needs to be output
                if( seismic_plug ) call define_seismic_points

                ! If RMS plugin is turned on create the geocellular model
                if( rms_plug .and. nproc == 1 ) call create_geocellular_model
                if( rms_plug .and. nproc > 1 .and. iam == 0 )&
                    print*,'Warning: The RMS plugin only works on serial version of LECODE.'
                endflag = 1

                ! Exit the event control loop
                exit event_control

            ! If we've to write a new display
            elseif( nextevent == SPM_DISPLAY .and. tnow > time_start )then

                if( stepin > 0 )then

                    ! Get compaction layers
                    if( gporo%compaction ) call cmpt_compaction

                    ! Compute mass movement
                    if( massmove ) call cmpt_massmovement

                    ! Deposit hemipelagites
                    if( hemi_flag ) call hemipelagites

                    ! Compute diffusion
                    call perform_top_sediment_layer_diffusion( 1 )

                    ! If check pointing on, records sedimentary layers
                    if( checkpointing .and. mod( stepc, checkfreq ) == 0 )then
                        call checkpointer( stepc )
                        ! Deallocate check point arrays
                        deallocate( rec_coord, rec_topelev, rec_NbLay, rec_ssedh )
                        deallocate( rec_slayID, rec_shardness, rec_szelev )
                        deallocate( rec_sporosity, rec_sthick, rec_sfac )

                    ! If clast time step, records sedimentary layers
                    elseif( tnow == time_end )then
                        call checkpointer( stepc )
                        ! Deallocate check point arrays
                        deallocate( rec_coord, rec_topelev, rec_NbLay, rec_ssedh )
                        deallocate( rec_slayID, rec_shardness, rec_szelev )
                        deallocate( rec_sporosity, rec_sthick, rec_sfac )
                    endif

                    ! Increment layers ID
                    L_strat%layersID = L_strat%layersID + 1

                    ! Check stratigraphic layers
                    if( tnow <= time_end ) call update_stratigraphic_layers

                endif

                ! Check if an output is required
                if( stepout == 0 .or. output_sim .or. tnow == time_end )then
                    call xdmf_output( stepout )
                    stepout = stepout + 1
                endif

                ! Update global parameters
                num_fw_o = 0
                rcd = 0
                rcdin = 0
                stepin = stepin + 1
                stepc = stepc + 1

                ! Deallocate record arrays
                if( fexistx )then
                    deallocate( recfw_type, recfw_xpos, recfw_ypos, recfw_zpos )
                    deallocate( recfw_xvel, recfw_yvel, recfw_srch, recfw_width )
                    deallocate( recfw_volume, recfw_sedcharge )
                endif

                ! If some flow walkers are still active update their TIN's positions
                if( tot_fw > 0 ) call update_flow_walker_TINposition
                fexistx = .false.
                flagx = .true.

            ! In case of new displacements
            elseif( nextevent == SPM_DISPLACEMENT )then

                ! Apply displacements
                call apply_displacements

                ! Update displacement rates if required
                if( new_disp .and. tnow < time_end )then
                    if( udw_plug .and. tnow > time_start )then
                        call create_vtk_topsurface( 0 )
                        call udw_plugin_wait_function
                    endif
                    call load_displacements( 0 )
                endif
                gdisp%lastdisp = tnow
                rcd = 0

            ! In case of rain event
            elseif( nextevent == SPM_RAINFALL )then

                ! Get new rain flow walkers
                if( geomorphology ) call register_rain_arrays

                ! Check flow walkers
                call clean_flow_walkers

                ! Define their IDs
                call define_flow_walkers_global_ids
                rcd = 0

            ! In case of inflow walkers
            elseif( nextevent == SPM_INFLOW )then

                ! Get new river flow walkers
                if( num_src > 0 ) call define_river_fws

                ! Check flow walkers
                call clean_flow_walkers

                ! Define their IDs
                call define_flow_walkers_global_ids
                rcd = 0

            ! Update sea level fluctuations
            elseif( nextevent == SPM_SEALEVEL )then

                call eustatism
                rcd = 0

            ! Move flow walkers over the simulation area
            elseif( nextevent == SPM_OCEAN )then

                if( tnow > time_start .and. abs( time_end - tnow ) > tor )then
                    ! Output SPM relative elevation
                    call create_oceanplug_surface

                    ! Wait for Ocean plugin to finish
                    call ocean_plugin_wait_function
                endif

                ! Move sediments according to wave & currents induced circulation
                if( abs( time_end - tnow ) > tor )then

                    ! Update morphological changes
                    call cmpt_ocean_transport

                    ! Compute diffusion
                    call perform_top_sediment_layer_diffusion( 0 )

                endif

            ! Move flow walkers over the simulation area
            elseif( nextevent == SPM_STEP )then

                ! If there are active flow walkers
                if( tot_fw > 0 )then
                    fexistx = .true.

                    ! Update the slope vectors of the local flow walkers.
                    call calculate_local_slope

                    ! Update the acceleration vectors of the local flow walkers.
                    call calculate_acceleration

                    ! Apply Cash-Karp Runge Kutta method, and update timestep parameters.
                    ! This also updates the FWs' acceleration and velocity components.
                    call advance_runge_kutta

                    ! Update current time
                    tnow = tnow + fine_dt

                    ! Clean up flow walkers
                    call clean_flow_walkers

                    ! Update number of FW
                    call mpi_allreduce( num_fw, tot_fw, 1, int_type, sum_type, SPModel_comm_world, ierr )

                    ! Exchange the flow walkers.
                    if( tot_fw >  0 ) call exchange_flow_walkers

                    ! Apply deposition/erosion rules
                    if( tot_fw >  0 ) call apply_erosion_deposition_rules

                    ! Record flow walkers positions and velocities
                    if( rcd == 0 )then
                        call record_flow_walkers( rcdin )
                        rcd = rcd + 1
                    else
                        rcd = rcd + 1
                        if( rcd == 50 ) rcd = 0
                    endif
                    rcdin = 1

                ! If there is no more flow walkers
                else

                    ! Update time step
                    dt = shortest_interval( ) * secyear
                    tnow = tnow + shortest_interval( )
                    rcd = 0
                endif

            ! Perform carbonate/organic evolution
            elseif( nextevent == SPM_CARBORG )then

                ! Update carbonate/organic based on fuzzy logic
                call compute_carborgs_evolution

                ! Perform diffusion
                call perform_top_sediment_layer_diffusion( 0 )

            ! If sediments need to be diffused
            elseif( nextevent == SPM_SLPDIFFUSION )then

                ! Check for mass movements
                if( massmove ) call cmpt_massmovement

                ! Add hemipelagites deposits
                if( hemi_flag ) call hemipelagites

                ! Perform diffusion
                call perform_top_sediment_layer_diffusion( 1 )

            endif

            ! Get the new process event
            oldevent = nextevent
            call control_process_events

            call mpi_barrier( SPModel_comm_world,ierr )
            if(iam == 0 .and. nextevent /= SPM_STEP )print*,'-----------------------------------------------'
            if(iam == 0 .and. nextevent /= SPM_STEP )write(*,101)'Current time : ',tnow, &
                ' next calling event nb : ',nextevent
            if(iam == 0 .and. nextevent == SPM_STEP .and. step >= 250 )&
                write(*,102)'Nb of flow walkers : ',tot_fw
            if( step == 250 ) step = 0
            step = step + 1

        enddo event_control

101     format( a15, f16.3, a25, i3 )
102     format( a21, i6 )

        return

    end subroutine spm_run
    ! ============================================================================
    !> Subroutine define_flow_walkers_global_ids
    !! Subroutine define_flow_walkers_global_ids controls the global ID of the flow walkers
    !! after rain and river initialisation.
    !<
    ! ============================================================================
    subroutine define_flow_walkers_global_ids

        integer :: k, pid, kn, p, maxref, refT, locid
        integer, dimension( nproc ) :: numel

        maxref = 0
        do k = 1, num_fw
            maxref = max( maxref, fa_refid( k ) )
        enddo
        call mpi_allreduce( maxref, refT, 1, int_type, max_type, SPModel_comm_world, ierr )

        ! Gather all flow walkers to the processors
        numel = 0
        call mpi_allgather( num_fw, 1, int_type, numel, 1, int_type, SPModel_comm_world, ierr )
        disp( 1 ) = 0
        do k = 1, nproc - 1
            disp( k + 1 ) = disp( k ) + numel( k )
        enddo
        call mpi_allreduce( num_fw, tot_fw, 1, int_type, sum_type, SPModel_comm_world, ierr )

        p = 0
        do k = 1, tot_fw
            pid = -1
            ! Find the sending processor
            do while ( pid < 0 )
                if( k > disp( p + 1 ) ) kn = p + 1
                if( k <= disp( p + 1 ) ) pid = p - 1
                if( kn > nproc - 1 ) pid = nproc - 1
                if( kn > nproc - 1 ) kn = nproc - 1
                p = kn
            enddo
            if( p > 0 ) p = p - 1
            if( pid  > nproc .and. iam == 0 ) print*,'Something went wrong allocating processors ID for FW global ref computation'

            ! Check if the considered flow walker needs to be send
            locid = k - disp( pid + 1 )
            if( iam == pid )then
                if( fa_refid( locid ) == -1 ) fa_refid( locid ) = refT + k
            endif
        enddo

        return

    end subroutine define_flow_walkers_global_ids
    ! ============================================================================
    !> Subroutine calculate_local_slope
    !! Subroutine calculate_local_slope runs through all flow walkers in the process
    !! and determines their slope vector x/y components.
    !<
    ! ============================================================================
    subroutine calculate_local_slope

        integer :: i, k

        real( tkind ) :: xsa, ysa

        ! Walk over all processor-local flow walkers.
        do k = 1, num_fw

            ! Initialise temporary slope component variables.
            xsa = 0.0_8
            ysa = 0.0_8

            ! Only calculate slope if the FW is inside the simulation boundary.
            if( fa_inbound( k ) == 0 )then

                ! Get the last TIN face occupied by the FW.
                i = fa_faceID( k , 1 )
                if( i > G_tin%faces_nb )then
                    print*,iam,k
                    print*,"Something went wrong, face number is outside range:",i,k,fa_faceID( k , 1 ),G_tin%faces_nb
                    stop
                endif

                ! Calculate the z-position of the FW.
                call find_elevation_within_TINface( i, k, fa_zpos( k ) )

                ! If the three vertices of the face have the same z-coordinate
                ! as the FW, the face is horizontal and there is no slope.
                if( tinv_coord( tinf_nids( i, 1 ),  3 ) == fa_zpos( k ) .and. &
                    tinv_coord( tinf_nids( i, 2 ), 3 ) == fa_zpos( k ) .and. &
                    tinv_coord( tinf_nids( i, 3 ), 3 ) == fa_zpos( k ) )then
                    xsa = 0.0_8
                    ysa = 0.0_8

                ! If not, but if the FW is above sea level or denser than sea water,
                ! calculate slope as an average gradient of the FW point and
                ! surrounding points.
                elseif( fa_zpos( k ) > gsea%actual_sea .or. fa_density( k ) >= sea_density ) then
                    call find_slope_within_TINface( i, k, xsa, ysa )

                ! Otherwise, assume the FW feels no slope.
                else
                    xsa = 0.0_8
                    ysa = 0.0_8
                endif
            endif

           ! Store the temp slope components.
            fa_xslp( k ) = dble( xsa )
            fa_yslp( k ) = dble( ysa )
        enddo

        ! Ensure that all processes have successfully completed subroutine.
        call completion

        return

    end subroutine calculate_local_slope
    ! ============================================================================
    !> Subroutine control_process_events
    !! Subroutine control_process_events controls SP Model execution by generating events for main
    !! processing loop.
    !<
    ! ============================================================================
    subroutine control_process_events

        integer    :: tsea, period

        real( tkind ) :: elapsedtime
        real( tkind ), save  :: previous_time

        elapsedtime = tnow - previous_time
        previous_time = tnow
        nextevent = SPM_UNDEFINED
        period = gdisp%actual

        ! Next ocean circulation event
        if( tnow >= next_oceancirc .and. ocean_plug )then
            tnow = next_oceancirc
            next_oceancirc = next_oceancirc + ocean_time_step
            nextevent = SPM_OCEAN

        ! Next carbonate/organic event
        elseif( tnow >= next_carborg .and. totgrn - silgrn > 0 )then
            tnow = next_carborg
            next_carborg = next_carborg + time_carborg
            nextevent = SPM_CARBORG

        ! Next displacement event
        elseif( tnow >= next_displacement .and. next_display > time_start )then
            tnow = next_displacement
            next_displacement = sampling_int + next_displacement

            if(gdisp_time( period , 2 ) < next_displacement ) &
                next_displacement = gdisp_time( period , 2 )

            if( gdisp_time( period , 2 ) <= tnow .and. tnow < time_end ) then
                new_disp = .true.
                next_displacement = next_display
                if( gdisp_time( period + 1, 2 ) > time_start .and. gdisp_time( period  + 1, 2 ) < next_displacement ) &
                    next_displacement = gdisp_time( period  + 1, 2 )
            endif

            if( tnow >= time_end ) next_displacement = sampling_int + time_end
            nextevent = SPM_DISPLACEMENT

        ! Next diffusion event
        elseif( tnow >= next_slpdiff )then
            tnow = next_slpdiff
            next_slpdiff = slope_int + next_slpdiff
            nextevent = SPM_SLPDIFFUSION

        ! Next sea level event
        elseif( tnow >= next_coarse_dt )then
            tsea = int( ( tnow - time_start + tor )/coarse_dt ) + 1
            next_coarse_dt = tsea * coarse_dt + time_start
            nextevent = SPM_SEALEVEL

        ! Next display event
        elseif( tnow >= next_display )then
            tnow = next_display
            next_display = time_display + next_display
            output_sim = .false.
            if( tnow >= next_output )then
                output_sim = .true.
                next_output = time_output + next_output
            endif
            nextevent = SPM_DISPLAY

        ! Next rain event
        elseif( tnow >= next_rain .and. tnow < time_end )then
            tnow = next_rain
            next_rain = next_rain + rain_int
            nextevent = SPM_RAINFALL

        ! Next inflow event
        elseif( tnow >= time_next .and. tnow < time_end )then
            time_next = time_next + flow_int
            nextevent = SPM_INFLOW

        ! Final event
        elseif( tnow >= time_end )then
            nextevent = SPM_QUIT

        ! Flow walker evolution event
        else
            nextevent = SPM_STEP
        endif

        return

    end subroutine control_process_events
    ! ============================================================================

end module mod_idle
