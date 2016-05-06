! ============================================================================
! Name        : AdvanceFW.f90
! Author      : tristan salles
! Created on: March 20, 2014
! Copyright (C) 2014 CSIRO
! ============================================================================
!> \file AdvanceFW.f90
!!
!! AdvanceFW advances flow walkers through time.
!!
!<
! ============================================================================
module mod_advancefw

    use flux_data
    use strata_ini
    use mpidata
    use time_data
    use fwalker_data
    use strata_data
    use param_data
    use mod_sealevel
    use TIN_function
    use TIN_surface
    use mod_layersFW
    use mod_rungekutta
    use mod_direction
    use mod_slopefw
    use mod_edtransfer
    use sediment_data

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine build_mpi_flowtransfer_type
    !! Subroutine build_mpi_flowtransfer_type defines a mpi data type for passing FW between
    !! processors for erosion deposition.
    !<
    ! ============================================================================
    subroutine build_mpi_flowtransfer_type

        integer, parameter :: count = 12
        integer, dimension( count ) :: array_of_types, array_of_block_lengths
        integer( kind = mpi_address_kind ), dimension( count ) :: array_of_displacements, array_of_address

        type( flowTransfer ) :: sample

        integer :: k

        array_of_types=(/int_type,int_type,int_type,&
            int_type,dbl_type,dbl_type,dbl_type,dbl_type,&
            dbl_type,int_type,int_type,dbl_type/)

        array_of_block_lengths = (/1,1,1,1,1,1,1,1,1,1,20,max_grn/)

        ! Set up the displacements and base types for the FW variables.
        call mpi_get_address(sample%procid, array_of_address( 1 ), ierr)
        call mpi_get_address(sample%refid, array_of_address( 2 ), ierr)
        call mpi_get_address(sample%bound, array_of_address( 3 ), ierr)
        call mpi_get_address(sample%type, array_of_address( 4 ), ierr)
        call mpi_get_address(sample%zpos, array_of_address( 5 ), ierr)
        call mpi_get_address(sample%height, array_of_address( 6 ), ierr)
        call mpi_get_address(sample%volume, array_of_address( 7 ), ierr)
        call mpi_get_address(sample%density, array_of_address( 8 ), ierr)
        call mpi_get_address(sample%velocity, array_of_address( 9 ), ierr)
        call mpi_get_address(sample%AInb, array_of_address( 10 ), ierr)
        call mpi_get_address(sample%AIpts, array_of_address( 11 ), ierr)
        call mpi_get_address(sample%sedcharge, array_of_address( 12 ), ierr)

        do k = 1, 12
            array_of_displacements( k ) = array_of_address( k ) - array_of_address( 1 )
        enddo

        ! Create and commit the type for MPI.
        call mpi_type_create_struct(count, &
            array_of_block_lengths, &
            array_of_displacements, &
            array_of_types, &
            mpi_flowtransfer, ierr)
        call mpi_type_commit(mpi_flowtransfer, ierr)

        ! Allocate erosion deposition arrays
        allocate( ed_refid( maxfw ) )
        allocate( ed_procid( maxfw ) )
        allocate( ed_bound( maxfw ) )
        allocate( ed_type( maxfw ) )
        allocate( ed_AIpts( maxfw, maxPtsFW ) )
        allocate( ed_AInb( maxfw ) )
        allocate( ed_zpos( maxfw ) )
        allocate( ed_h( maxfw ) )
        allocate( ed_vel( maxfw ) )
        allocate( ed_vol( maxfw ) )
        allocate( ed_dens( maxfw ) )
        allocate( ed_sed( maxfw, max_grn ) )

        return

    end subroutine build_mpi_flowtransfer_type
    ! ============================================================================
    !> Subroutine build_mpi_flowwalker_type
    !! Subroutine build_mpi_flowwalker_type defines a mpi data type for passing FW between
    !! processors.
    !<
    ! ============================================================================
    subroutine build_mpi_flowwalker_type

        integer, parameter :: count = 25
        integer, dimension( count ) :: array_of_types, array_of_block_lengths
        integer( kind = mpi_address_kind ), dimension( count ) :: array_of_displacements, array_of_address

        type( flowWalker ) :: sample

        integer :: k

        array_of_types=(/int_type,int_type,int_type,int_type,&
            dbl_type,dbl_type,dbl_type,dbl_type,dbl_type,&
            dbl_type,dbl_type,dbl_type,dbl_type,dbl_type,&
            dbl_type,dbl_type,dbl_type,dbl_type,dbl_type,&
            dbl_type,dbl_type,int_type,int_type,dbl_type,&
            int_type/)

        array_of_block_lengths = (/1,1,1,1,1,1,1,1,1,1,1,1,1,1,&
            1,1,1,1,1,1,1,1,20,max_grn,numold/)

        ! Set up the displacements and base types for the FW variables.
        call mpi_get_address(sample%refid, array_of_address( 1 ), ierr)
        call mpi_get_address(sample%inbound, array_of_address( 2 ), ierr)
        call mpi_get_address(sample%count, array_of_address( 3 ), ierr)
        call mpi_get_address(sample%type, array_of_address( 4 ), ierr)
        call mpi_get_address(sample%xpos, array_of_address( 5 ), ierr)
        call mpi_get_address(sample%ypos, array_of_address( 6 ), ierr)
        call mpi_get_address(sample%zpos, array_of_address( 7 ), ierr)
        call mpi_get_address(sample%xslp, array_of_address( 8 ), ierr)
        call mpi_get_address(sample%yslp, array_of_address( 9 ), ierr)
        call mpi_get_address(sample%pgrad, array_of_address( 10 ), ierr)
        call mpi_get_address(sample%height, array_of_address( 11 ), ierr)
        call mpi_get_address(sample%width, array_of_address( 12 ), ierr)
        call mpi_get_address(sample%xvel, array_of_address( 13 ), ierr)
        call mpi_get_address(sample%yvel, array_of_address( 14 ), ierr)
        call mpi_get_address(sample%volume, array_of_address( 15 ), ierr)
        call mpi_get_address(sample%density, array_of_address( 16 ), ierr)
        call mpi_get_address(sample%pzfe, array_of_address( 17 ), ierr)
        call mpi_get_address(sample%pcoeff, array_of_address( 18 ), ierr)
        call mpi_get_address(sample%pdist, array_of_address( 19 ), ierr)
        call mpi_get_address(sample%pvelx, array_of_address( 20 ), ierr)
        call mpi_get_address(sample%pvely, array_of_address( 21 ), ierr)
        call mpi_get_address(sample%AInb, array_of_address( 22 ), ierr)
        call mpi_get_address(sample%AIpts, array_of_address( 23 ), ierr)
        call mpi_get_address(sample%sedcharge, array_of_address( 24 ), ierr)
        call mpi_get_address(sample%faceid, array_of_address( 25 ), ierr)

        do k = 1, 25
            array_of_displacements( k ) = array_of_address( k ) - array_of_address( 1 )
        enddo

        ! Create and commit the type for MPI.
        call mpi_type_create_struct(count, &
            array_of_block_lengths, &
            array_of_displacements, &
            array_of_types, &
            mpi_flowwalker, ierr)
        call mpi_type_commit(mpi_flowwalker, ierr)

        return

    end subroutine build_mpi_flowwalker_type
    ! ============================================================================
    !> Subroutine calculate_acceleration
    !! Subroutine calculate_acceleration considers gravity, Manning coefficient,
    !! and plume status to calculate the x/y components of each local flow walker's
    !! acceleration vector.
    !<
    ! ============================================================================
    subroutine calculate_acceleration

        integer :: k, fce
        real( tkind ) :: xsb, ysb, xnn, fv2, ga

        ! Loop over all flow walkers.
        do k = 1, num_fw

            ! Only update acceleration if FW is within simulation boundaries.
            if( fa_inbound( k ) == 0 )then

                ! Get the last TIN face occupied by the FW.
                fce = fa_faceID( k , 1 )

                ! Increment the FW step counter.
                fa_count( k ) = fa_count( k ) + 1

                ! Calculate the Manning coefficient of the FW.
                xnn = manning_fct( gsea%actual_sea - fa_zpos( k ), &
                    fa_density( k ) )

                ! Calculate the effect of gravity, scaled by the ratio of FW density
                ! to seawater density.
                ga = gravity
                if( fa_zpos( k ) <=  gsea%actual_sea )then
                    ga = ga * ( fa_density( k ) - sea_density ) / fa_density( k )
                endif

                ! Calculate the FW's coefficient of friction.
                fa_Cfric( k ) = dble( gravity * xnn**2 / fa_h( k )**(4/3) )

                ! If the FW represents a flow above sea-level, or if the
                ! FW has density greater than or equal to sea-water density,
                ! set the slope as derived from previous location(s) of the FW.
                if( ( fa_zpos( k ) > gsea%actual_sea .and. fa_type( k ) < tpplum1 ) &
                    .or. fa_density( k ) >= sea_density ) then
                    call compute_slope_along_travel_direction( k, fce, xsb, ysb )
                    fa_xslp( k ) = xsb
                    fa_yslp( k ) = ysb

                ! Otherwise, if we are considering plumes, convert the FW into a plume element
                ! using the Syvitski plume model.
                elseif( plume_cmpt )then

                     ! Calculate the Euclidean norm of the FW's current velocity.
                     fv2 = sqrt( fa_xvel( k )**2 + fa_yvel( k )**2 )

                    ! Calculate plume parameters.
                    if( fa_type( k ) < tpplum1 .and. fv2 > 0.0_8 )then
                        fa_type( k ) = tpplum1
!                        fa_pzfe( k ) = 7.0_8 * sqrt( pi * 0.5_8 * ( river_mouth )**2 )
                        ! Based on Syvitski Plume model the following approximation is used
                        fa_pzfe( k ) = 5.0_8 * river_mouth
                        fa_pcoeff( k ) =  -log( 0.001_8 /  fv2 ) / ( ( plume_lenght - fa_pzfe( k ) ) )
                        fa_pdist( k ) = 0.0_8
                        fa_pvelx( k ) = fa_xvel( k )
                        fa_pvely( k ) = fa_yvel( k )
                    endif
                    xsb = 0.0_8
                    ysb = 0.0_8
                ! Otherwise, assume zero slope.
                else
                    xsb = 0.0_8
                    ysb = 0.0_8
                    fa_xslp( k ) = xsb
                    fa_yslp( k ) = ysb
                endif

                ! Set the acceleration components according to slope and gravity.
                fa_accx( k ) = dble( - ga * xsb )
                fa_accy( k ) = dble( - ga * ysb )
            endif

        enddo

        ! Ensure that all processes have successfully completed subroutine.
        call completion

        return

    end subroutine calculate_acceleration
    ! ============================================================================
    !> Subroutine advance_runge_kutta
    !! Subroutine advance_runge_kutta repeatedly applies a Runge-Kutta solver,
    !! and adapts the timestep width according to observed error. The subroutine
    !! also handles timestep evolution for the plume case.
    !<
    ! ============================================================================
    subroutine advance_runge_kutta

        integer :: k, faceLocation, nodeLocation

        real(tkind) :: large_dt, dtplume

        ! Initialise delta-t values.
        large_dt = shortest_interval( ) * secyear / transport%sed_dt
        if (tnow == time_start) dtnext = large_dt
        dt = min( 5.0e2_8, dtnext )
        if (dt < tor) dt = tor

        ! Run adaptive Runge-Kutta solver and get the next timestep width.
        ! TODO: this uses an all-reduce at each Runge-Kutta call, which
        ! is going to be painful for scaling. Can we do without it?
        dtnext = runge_kutta()

        ! Handle the timestep evolution in case we're dealing with (potential) plumes.
        if( plume_cmpt )then
            dtplume = plume_step()
            dtnext = min(dtplume, dtnext)
        endif

        ! Update the flow walkers' acc/vel components.
        call update_velocity_and_acceleration(large_dt)

        ! Update various parameters for FWs, including height,
        ! width, sediment load, etc.
        call advance_flow_walkers

        ! Update fine timestep.
        fine_dt = dble( dt * transport%sed_dt / secyear )

        ! Ensure that all processes have successfully completed subroutine.
        call completion

        ! If we need to combine the various rain flowwalkers, then set up the
        ! fw_gid array and go nuts.
        if( combine == 2 )then

            ! Set up fw_gid
            allocate( fw_gid( num_fw ) )
            fw_gid = -1
            do k = 1, num_fw
                if( fa_inbound( k ) /= 1 )then
                    call closest_stratal_node( fa_xpos( k ), fa_ypos( k ), faceLocation, nodeLocation )
                    fw_gid( k ) = nodeLocation
                endif
            enddo

            ! Combine the rainwalkers.
            call combine_rainwalkers
            deallocate( fw_gid )

        endif

        ! Update all processes about the number of flow walkers.
        call mpi_allreduce( num_fw, tot_fw, 1, int_type, sum_type, SPModel_comm_world, ierr )

        ! Find the nodes in the area of influence of the FW
        if( num_fw > 0 ) call points_area_of_influence

        return

    end subroutine advance_runge_kutta
    ! ============================================================================
    !> Subroutine update_velocity_and_acceleration
    !! Subroutine update_velocity_and_acceleration updates flow walker acceleration and velocity.
    !! \param large_dt
    !<
    ! ============================================================================
    subroutine update_velocity_and_acceleration( large_dt )

        integer :: k, tclas

        real( tkind ) :: sc1, sc2, large_dt, frac, xpos, ypos
        real( tkind ) :: fwacc, grdnew, ududx, dudt !, dir

        do k = 1, num_fw

            ! In case the FW is a gravity current
            if ( fa_type( k ) >= tpturb .and. fa_type( k ) < tpplum1 )then

                fwacc = sqrt( fa_xvelO( k )**2 + fa_yvelO( k )**2 )
                fwacc = fwacc - sqrt( fa_xvel( k )**2 + fa_yvel( k )**2 )
                fwacc = fwacc / dt
                grdnew = sqrt( fa_accx( k )**2 + fa_accy( k )**2 )

                ! approximation of u*du/dx =grdnew - gradold
                ududx = grdnew - fa_pgrad( k )

                ! approximation of dudt= DU/dt-grdnew
                dudt = fwacc - ududx

                if( fwacc > 0.0_8 )then
                    tclas = 0
                else
                     ! Accumulative flow class (must be 1 or 2)
                    if( ududx > 0.0_8 )then

                        if( abs( dudt / ududx ) < 1.0_8 )then
                            ! ududdx is strongly +ve, must be type 1
                            tclas = 1
                        else
                            ! ududx is not very strong, class as 2
                            tclas = 2
                        endif
                    ! Depletive flow class (must be 2, 3, 4, or 5)
                    else
                        ! Waxing flow type (must be 4 or 5)
                        if( dudt > 0.0_8 )then
                            if( abs( ududx / dudt ) < 1.0_8 )then
                                ! dudt is strongly +ve, must be type 5
                                tclas = 5
                            else
                                ! dydt isn't very strong, class as 4
                                tclas = 4
                            endif

                        elseif( dudt == 0.0_8 )then
                            tclas = 4

                        else
                            ! ok, both terms are negative , must be 2, 3 or 4
                            if( abs( ududx / dudt ) < 0.25_8 )then
                                tclas = 2
                            elseif( abs( ududx / dudt ) > 4.0_8 )then
                                tclas = 4
                            else
                                tclas = 3
                            endif
                        endif
                    endif
                endif
                fa_type( k ) = tpturb + tclas
                fa_pgrad( k ) = grdnew

            ! In case the FW is a plume
            elseif( fa_type( k ) >= tpplum1 )then

                xpos = fa_xpos( k )
                ypos = fa_ypos( k )

                ! Zone of Flow Establishement (ZFE)
                if( fa_pzfe( k ) > fa_pdist( k ) )then
                    fa_xvelO( k ) = fa_xvel( k )
                    fa_yvelO( k ) = fa_yvel( k )
                    fa_xpos( k ) = real( xpos + fa_xvelO( k ) * dt )
                    fa_ypos( k ) = real( ypos + fa_yvelO( k ) * dt )
                    fa_pdist( k ) = fa_pdist( k ) + sqrt( ( fa_xpos( k ) - xpos )**2 + &
                        ( fa_ypos( k ) - ypos )**2 )
                    fa_type( k ) = tpplum1

                ! Zone of Established Flow (ZEF)
                else
                    frac = ( fa_pdist( k ) - fa_pzfe( k ) ) / ( plume_lenght )
                    if( frac > 1.0_8 ) frac = 1.0_8
                    fa_xpos( k ) = real( xpos + fa_xvel( k ) * dt + &
                        frac * ocean_flow( 1 ) * dt )
                    fa_ypos( k ) = real( ypos + fa_yvel( k ) * dt + &
                        frac * ocean_flow( 2 ) * dt )
                    fa_pdist( k ) = fa_pdist( k ) + sqrt( ( fa_xpos( k ) - xpos )**2 + &
                        ( fa_ypos( k ) - ypos )**2 )
                    ! Get plume centerline direction
!                    call angle_alpha( xpos, ypos, fa_xpos( k ), fa_ypos( k ), dir )
                    fa_xvelO( k ) = fa_pvelx( k ) * exp( - fa_pcoeff( k ) * &
                        ( fa_pdist( k ) - fa_pzfe( k ) ) )
                    fa_yvelO( k ) = fa_pvely( k ) * exp( - fa_pcoeff( k ) * &
                        ( fa_pdist( k ) - fa_pzfe( k ) ) )
                    fa_type( k ) = tpplum2
                endif
            endif
        enddo

        if( dt > large_dt )then
            sc1 = ( dt - large_dt ) / dt
            sc2 = large_dt / dt
            dt = large_dt
            do k = 1, num_fw
                fa_xvelO( k ) = dble( fa_xvelO( k ) * sc2 + fa_xvel( k ) * sc1 )
                fa_yvelO( k ) = dble( fa_yvelO( k ) * sc2 + fa_yvel( k ) * sc1 )
            enddo
        endif

        return

    end subroutine update_velocity_and_acceleration
    ! ============================================================================
    !> Subroutine advance_flow_walkers
    !! Subroutine advance_flow_walkers updates flow walker position and discharge.
    !<
    ! ============================================================================
    subroutine advance_flow_walkers

        integer :: i, k, face, nFaces
        logical :: inside
        real(tkind) :: dist, zposition, tempWidth
        real(tkind) :: widthDepthRatio, fv1, fv2, grad, sedLoad, area

        do k = 1, num_fw

            ! If flow walker is in the simulation space, update flow walker positions
            ! and calculate their (potentially) new TIN faces.
            if( fa_inbound( k ) == 0 )then

                ! Using new delta-t and updated velocity, find the new position of the flow walker
                if( fa_type( k ) < tpplum1 )then
                    fa_xpos( k ) = real( fa_xpos( k ) + fa_xvelO( k ) * dt )
                    fa_ypos( k ) = real( fa_ypos( k ) + fa_yvelO( k ) * dt )
                endif
                face = fa_faceID( k , 1 )
                if( face > G_tin%faces_nb ) &
                    print*,"Something went wrong, face number is outside range...",face,G_tin%faces_nb

                ! Tag FW as outside simulation, until proved otherwise.
                fa_inbound( k ) = 1

                ! Find the face number within the TIN grid
                ! First check if we are still in the same face number
                face = fa_faceID( k, 1)
                call check_point_TIN_face( face, fa_xpos( k ), fa_ypos( k ), inside )
                if( inside )then
                    fa_inbound( k ) = 0

                ! If we have moved to another face, look in all the TIN faces neighbouring the prior face.
                else
                    neighbourSearch: do i = 1, tinf_nghb( face )
                        if( fids(face, i) > 0 )then
                            face = fids(face, i)
                            call check_point_TIN_face( face, fa_xpos( k ), fa_ypos( k ), inside )
                            if( inside ) then
                                fa_inbound( k ) = 0
                                exit neighbourSearch
                            endif
                         endif
                    enddo neighbourSearch


                    ! In case we still cannot find the face in the vicinity of previous one
                    ! do a search all over the TIN grid
                    if( fa_inbound( k ) == 1 )then
                        completeSearch: do i = 1, G_tin%faces_nb
                            face = i
                            call check_point_TIN_face( face, fa_xpos( k ), fa_ypos( k ), inside )
                            if( inside ) then
                                fa_inbound( k ) = 0
                                exit completeSearch
                            endif
                        enddo completeSearch
                    endif

                    ! If we are still outside the sim, something's gone seriously wrong.
                    if( fa_inbound( k ) == 1 )then
                        print *, 'Something went wrong locating new TIN face for FW.'
                        stop
                    endif

                endif

                ! Are we on the border of the mesh?
                if (tinf_boundary(face) > 0) then
                    fa_inbound( k ) = 1
                endif

            endif

            ! Based on new TIN face surface equation compute the slope
            if( fa_inbound( k ) == 0 .and. fa_type( k ) < tpplum1 )then

                ! Get z-position
                call find_elevation_within_TINface( face, k, zposition )

                ! If the three vertices of the face have the same z-coordinate
                ! as the FW, the face is horizontal and there is no slope.
                if( tinv_coord( tinf_nids( face, 1 ), 3 ) == zposition .and. &
                    tinv_coord( tinf_nids( face, 2 ), 3 ) == zposition .and. &
                    tinv_coord( tinf_nids( face, 3 ), 3 ) == zposition )then
                    fa_xslp( k ) = 0.0_8
                    fa_yslp( k ) = 0.0_8

                ! If not, but if the FW is above sea level or denser than sea water,
                ! calculate slope as an average gradient of the FW point and
                ! surrounding points.
                elseif( fa_zpos( k ) > gsea%actual_sea .or. fa_density( k ) >= sea_density ) then
                    call find_slope_within_TINface(  face, k, fa_xslp( k ), fa_yslp( k ) )

                ! Otherwise, assume the FW feels no slope.
                else
                    fa_xslp( k ) = 0.0_8
                    fa_yslp( k ) = 0.0_8

                endif

                ! Calculate the Euclidean norm of the FW's current and previous velocity.
                fv1 = sqrt( fa_xvel( k )**2 + fa_yvel( k )**2 )
                fv2 = sqrt( fa_xvelO( k )**2 + fa_yvelO( k )**2 )

                ! FW's travelling distance during time step delta-t
                dist = sqrt( fa_xvelO( k )**2 + fa_yvelO( k )**2 ) * dt

                ! Compute flow walker height and width based on stream classification
                if( dist > 0.0_8 )then

                    ! Compute slope between current and previous FW's position
                    grad = ( fa_zpos( k ) - zposition ) / dist

                    ! Based on flow discarge conservation A1v1 = A2v2 which gives the
                    ! new flow area normal to the velocity vector
                    area = fa_h( k ) *  fa_width( k ) * fv1 / fv2

                    ! If the FW is flowing upslope due to its acceleration make the assumption
                    ! that the flow width is conserved
                    if( grad < 0.0_8 )then
                        ! Conserve flow width
                        tempWidth = fa_width( k )
                        tempWidth = min(tempWidth, transport%wdtmax)
                        tempWidth = max(tempWidth, transport%wdtmin)
                        ! Update flow height
                        fa_h( k ) = area / tempWidth

                    ! Otherwise update FW height and width based on stream classification
                    else
                        widthDepthRatio = streamclass_fct( grad )
                        tempWidth = sqrt( area / widthDepthRatio ) * widthDepthRatio
                        tempWidth = min(tempWidth, transport%wdtmax)
                        tempWidth = max(tempWidth, transport%wdtmin)
                        ! Update flow height
                        fa_h( k ) = tempWidth / widthDepthRatio
                    endif

                    ! Clip height to within specified boundaries.
                    fa_h( k ) = min(transport%dzimax, fa_h( k ))
                    fa_h( k ) = max(transport%dzimin, fa_h( k ))

                else
                    tempWidth = fa_width( k )
                endif

                ! Consider entrainment for turbidity currents.
                if( fa_type( k ) >= tpturb ) call water_entrainment( k, fv2 )

                ! Set discharge and width.
                fa_discharge( k ) = fv2 * fa_h( k ) * tempWidth
                fa_width( k ) = tempWidth

            ! Consider plume FWs.
            elseif( fa_type( k ) >= tpplum1 )then

                ! Calculate the Euclidean norm of the FW's current and previous velocity.
                fv2 = sqrt( fa_xvelO( k )**2 + fa_yvelO( k )**2 )
                fa_discharge( k ) = fv2 * fa_h( k ) * fa_width( k )
                zposition = gsea%actual_sea

            else
                zposition = 1.e6_8
            endif

            ! Update velocity components and z-position, and identify
            ! state changes in flow walkers.
            if( fa_inbound( k ) == 0 )then

                fa_xvel( k ) = dble( fa_xvelO( k ) )
                fa_yvel( k ) = dble( fa_yvelO( k ) )
                fa_zpos( k ) = dble( zposition )

                ! Update TIN traversal history.
                if (face > 0 .and. dist > 0.0_8) then
                    call record_flow_walker_faceID(face, k)
                end if

                ! Count how many times the FW has been sitting in the same TIN face.
                nFaces = 0
                sameFaceCounter: do i = 2, numold
                    if( face == fa_faceID(k, i) ) then
                        nFaces = nFaces + 1
                    end if
                    if (fa_faceID(k, i) == 0) then
                        exit sameFaceCounter
                    end if
                end do sameFaceCounter

                ! Detect transition to ocean-based turbulent flow.
                if (fa_zpos( k ) < gsea%actual_sea &
                    .and. fa_type( k ) < tpturb &
                    .and. fa_density( k ) >= sea_density) then
                    fa_type( k ) = tpturb
                end if

                ! In the absence of plume consideration, detect FWs that are entering
                ! the sea and will be deposited rapidly.
                if (fa_zpos( k ) < gsea%actual_sea &
                    .and. .not. plume_cmpt &
                    .and. fa_type( k ) < tpturb) then
                    fa_inbound( k ) = 2
                end if

                ! If the FW is a plume element above sea level, zero its velocity components
                ! and set it to forced deposit.
                if (fa_zpos( k ) > gsea%actual_sea &
                    .and. fa_type( k ) >= tpplum1) then
                    fa_xvel( k ) = 0
                    fa_yvel( k ) = 0
                    fa_inbound( k ) = 2
                end if

                ! Calculate the sediment load of the FW.
                sedLoad = 0.0_8
                do i = 1, totgrn
                    sedLoad = sedLoad + fa_sedcharge(k, i) * sediment(i)%density
                end do

                ! Adjust sediment load by FW type.
                if (fa_type( k ) >= tprain .and. fa_type( k ) < tpturb) then
                    sedLoad = sedLoad / rain_int
                else
                    sedLoad = sedLoad / flow_int
                end if

                ! Check if FW is forced to deposit due to low sediment weight and high step count.
                if (fa_count( k ) >= 100 .and. sedLoad < transport%slsmin) then
                    fa_inbound( k ) = 2
                end if

                ! Check if rain FW is forced to deposit due to high step count.
                if (fa_count( k ) >= 5e5 .and. fa_type( k ) >= tprain .and. fa_type( k ) < tpturb) then
                    fa_inbound( k ) = 2
                end if

                ! Check if plume FW is forced to deposit due to exceeding distance from mouth.
                if (fa_type( k ) >= tpplum1 .and. fa_pdist( k ) >= plume_lenght * 0.99_8) then
                    fa_inbound( k ) = 2
                end if

                ! Check if plume FW is forced to deposit due to low movement.
                if (fa_type( k ) >= tpplum1 .and. fv2 < transport%fvmin) then
                    fa_inbound( k ) = 2
                end if

                ! Check if FW is forced to deposit due to being stationary for too long.
                if (nFaces > transport%stepmax) then
                    fa_inbound( k ) = 2
                end if
            end if

        enddo

        return

    end subroutine advance_flow_walkers
    ! ============================================================================
    !> Function shortest_interval
    !! Function shortest_interval determines the shortest interval time
    !<
    ! ============================================================================
    function shortest_interval( ) result( eval )

        real( tkind ) :: eval

        eval = 1.e16_8
        if( gdisp%event  > 0 ) eval = min( eval, next_displacement - tnow )
        eval = min( eval, next_display - tnow )
        eval = min( eval, time_next - tnow )
        eval = min( eval, next_coarse_dt - tnow )
        eval = min( eval, next_oceancirc - tnow )
        eval = min( eval, next_carborg - tnow )
        eval = min( eval, next_rain - tnow )
        eval = min( eval, next_slpdiff - tnow )
        if( eval < time_tolerance ) eval = time_tolerance

    end function shortest_interval
    ! ============================================================================
    !> Function streamclass_fct()
    !! Function streamclass_fct is used to compute the width/depth ratio based on slope.
    !! \param slp
    !<
    ! ============================================================================
    function streamclass_fct( slp ) result( Psi )

        integer :: k
        real( tkind ) :: slp, wd
        real( tkind ) :: Psi

        if( slp < streamclass( 1 )%slope )then
            wd = ( streamclass( 2 )%wd_ratio - streamclass( 1 )%wd_ratio ) / &
                ( streamclass( 2 )%slope - streamclass( 1 )%slope )
            Psi = ( streamclass( 1 )%wd_ratio + wd * ( slp - streamclass( 1 )%slope ) )
        elseif( slp > streamclass( stcl_nb )%slope )then
            Psi = streamclass( stcl_nb )%wd_ratio
        else
            lp: do k = 1, stcl_nb - 1
                if( slp >= streamclass( k )%slope .and. slp < streamclass( k + 1 )%slope )then
                    wd = ( streamclass( k + 1 )%wd_ratio - streamclass( k )%wd_ratio ) / &
                        ( streamclass( k + 1 )%slope - streamclass( k )%slope )
                    Psi = ( streamclass( k )%wd_ratio + wd * ( slp - streamclass( k )%slope ) )
                    exit lp
                endif
            enddo lp
        endif


    end function streamclass_fct
    ! ============================================================================
    !> Subroutine water_entrainment
    !! Subroutine water_entrainment computes water entrainment for turbidity currents.
    !! \param k, fv2
    !<
    ! ============================================================================
    subroutine water_entrainment( k, fv2 )

        integer :: k

        real( tkind ) :: fv2, Ri, ew

        ! Water entrainment for turbidity currents
        if( fa_type( k ) < tpplum1 .and. fa_type( k ) == tpturb + 1 .and. fv2 > 0.0_8 )then
            Ri = ( fa_density( k ) - sea_density ) * gravity / fa_density( k )
            Ri = Ri * fa_h( k ) / fv2**2
            ew = 1.0_8 + 718.0_8 * Ri**2.4_8
            ew = sqrt( ew )
            ew = 0.075_8 / ew
            ew = ew * dt * sqrt( fv2 )

            fa_volume( k ) = fa_volume( k ) + ew * strat_dx**2.0_8

        endif

        return

    end subroutine water_entrainment
    ! ============================================================================
    !> Subroutine angle_alpha
    !! Subroutine angle_alpha computes the angle between previous and current plume point positions.
    !!
    !<
    ! ============================================================================
    subroutine angle_alpha( c1x, c1y, c3x, c3y, alpha )

        real( tkind ) :: alpha, a, b, c, div
        real( tkind ) :: c1x, c1y, c3x, c3y
        real( tkind ), dimension( 2 ) :: c1, c2, c3

        c1( 1 ) = c1x
        c1( 2 ) = c1y
        c2( 1 ) = c1( 1 ) + 3.0_8 * strat_dx
        c2( 2 ) = c1( 2 )
        c3( 1 ) = c3x
        c3( 2 ) = c3y

        a = sqrt( ( c2( 1 ) - c3( 1 ) )**2 + ( c2( 2 ) - c3( 2 ) ) **2 )
        b = sqrt( ( c1( 1 ) - c3( 1 ) )**2 + ( c1( 2 ) - c3( 2 ) ) **2 )
        c = sqrt( ( c2( 1 ) - c1( 1 ) )**2 + ( c2( 2 ) - c1( 2 ) ) **2 )

        div = ( b**2 + c**2 - a**2 ) / ( 2.0_8 * b * c )

        if( div > 1.0_8 ) div = 1.0_8
        if( div < -1.0_8 ) div = -1.0_8

        alpha = acos( div )

        return

    end subroutine angle_alpha
    ! ============================================================================

end module mod_advancefw
