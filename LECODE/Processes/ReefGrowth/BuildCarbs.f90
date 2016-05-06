! ============================================================================
! Name        : BuildCarbs.f90
! Author      : tristan salles
! Created on: April 05, 2014
! Copyright (C) 2014 CSIRO
! ============================================================================
!> \file BuildCarbs.f90
!!
!! File BuildCarbs determines reefs evolution. The approach is based on fuzzy logic,
!! inspired by Ulf Nordlund code FUZZIM and is similar to Chris Dyt implementation in
!! Sedsim.
!!
!<
! ============================================================================
module mod_buildcarbs

    use mpidata
    use file_data
    use mod_sort
    use flux_data
    use time_data
    use error_data
    use ocean_data
    use strata_data
    use param_data
    use forces_data
    use mod_diffuse
    use mod_tempsal
    use sediment_data
    use mod_functioncarbs

    implicit none

    logical, dimension( 4 ) :: action_flag, action_subflag

    public

contains

    ! ============================================================================
    !> Subroutine compute_carborgs_evolution
    !! This subroutine determines carbs/orgs evolution. The approach is based on fuzzy logic,
    !! inspired by Ulf Nordlund code FUZZIM and is similar to Chris Dyt implementation in
    !! Sedsim.
    !<
    ! ============================================================================
    subroutine compute_carborgs_evolution

        integer :: n, k, ks, act, laynb, kl, dmat( totgrn )

        real( tkind ) :: max_vel_ocean, sumdep
        real( tkind ) :: carborg( 4, totgrn ), carborg_sed( totgrn )
        real( tkind ) :: carborg_subsurface( 4, totgrn, L_strat%layersID )

        proc_elev = -1.0e6_8

        if( flagx )then
            top_sedprev = 0.0_8
            top_sedh = 0.0_8
            flagx = .false.
        endif

        dmat = 0

        if( tnow == time_start ) call allocate_step_function

        ! Update the current temperature and salinity
        if( salinity_flag == 1 .or. temperature_flag == 1 ) call get_temperature_salinity

        ! Get the distance to nearest shoreline and rivers
        call get_distance_to_nearest_shore_and_rivers

        ! Get the distance to nearest materials
        if( mat_nearest_flag == 1 )then
            do n = 1, membership_fct_nb
                if( membership( n )%sedclassID > 0 )then
                    dmat( membership( n )%sedclassID ) = 1
                     call get_distance_to_nearest_materials( membership( n )%sedclassID )
                endif
            enddo
        endif

        ! Get slope from neighborhood nodes
        call get_averaged_slope_from_neighbours

        ! Get burial depth of each sedimentary layers
        call get_sedimentary_layer_burial_depth

        do n = 1, L_strat%nodes_nb

            ! Reset action flag
            action_flag( 1:4 ) = .false.
            action_subflag( 1:4 ) = .false.

            ! Get top elevation
            proc_elev( locv_gid( n ) ) = strat_zelev(  n ,  locv_NbLay( n )  )

            ! Find carbonate/organic paremeter values based on local declaration
            if( .not. allocated( carborg_parameter ) ) allocate( carborg_parameter( 13+totgrn ) )
            carborg_parameter = 0.0_8

            ! Water depth
            if( carb_depth_flag == 1 ) carborg_parameter( 1 ) = gsea%actual_sea - elev_record( locv_gid( n ) )
            ! Ocean circulation flag
            if( carb_ocean_flag == 1 )then
                max_vel_ocean = 0.0_8
                do k = 1, hindcast( active_hindcast )%cnb
                    max_vel_ocean = max( max_vel_ocean, ocean_current( k, locv_gid( n ), 1 ) )
                    max_vel_ocean = max( max_vel_ocean, ocean_wave( k, locv_gid( n ), 1 ) )
                enddo
                carborg_parameter( 2 ) = max_vel_ocean
            endif
            ! Shore distance
            if( carb_shore_flag == 1 ) carborg_parameter( 3 ) = distance_shore( n )
            ! River distance
            if( carb_river_flag == 1 ) carborg_parameter( 4 ) = distance_river( n )
            ! Material distance
            do k = 1, totgrn
                carborg_parameter( 4+k ) = -1.0_8
                if( mat_nearest_flag == 1 .and. dmat( k ) == 1 ) carborg_parameter( 4+k ) = distance_mats( n, k )
            enddo
            ! Sedimentation
            if( carb_sediment_flag == 1 ) call get_carborg_sedimentation( n, carborg_parameter( 5+totgrn ) )
            ! Temperature
            if( temperature_flag == 1 ) carborg_parameter( 6+totgrn ) = gtempsal%actual_tempsal( 1 )
            ! Salinity
            if( salinity_flag == 1 ) carborg_parameter( 7+totgrn ) = gtempsal%actual_tempsal( 2 )
            ! Valley
            if( carb_valley_flag == 1 ) carborg_parameter( 8+totgrn ) = nslp_carbs( n )
            ! Gradient
            if( carb_slope_flag == 1 ) carborg_parameter( 9+totgrn ) = aslp_carbs( n )
            ! Exposed time
            if( carb_exposed_flag == 1 ) call get_carborg_exposed_time( n, carborg_parameter( 10+totgrn ) )
            ! Sub-surface parameters
            if( carb_burial_flag == 1 .or. carb_age_flag == 1 )then
                do k = 1, locv_NbLay( n )
                    if( carb_burial_flag == 1 ) carborg_parameter( 11+totgrn ) = burial_carbs( n, k )
                    if( carb_age_flag == 1 ) call get_carborg_buriedlayer_time( n, k, carborg_parameter( 12+totgrn ) )
                    ! For each sedimentary layer compute fuzzy function for subsurface effect
                    call compute_fuzzy_logic( carborg, 1 )
                    do ks = 1, totgrn
                        do act = 1, 4
                            carborg_subsurface( act, ks, k ) = carborg( act, ks )
                            if( carborg_subsurface( act, ks, k ) > 0.0_8 ) action_subflag( act ) = .true.
                        enddo
                    enddo
                enddo
            endif
            ! Reset burial to zero for surface effects
            if( carb_burial_flag == 1 ) carborg_parameter( 11+totgrn ) = 0.0_8
            ! When on the surface, age is set to exposed time
            if( carb_age_flag == 1 ) carborg_parameter( 12+totgrn ) = carborg_parameter( 10+totgrn )
            ! All carbonate parameters have been initialised work on surface effects
            call compute_fuzzy_logic( carborg, 0 )
      
            ! Growth
            act = 1
            if( action_flag( act ) )then
                sumdep = 0.0_8
                do ks = 1, totgrn
                    if( carborg( act, ks ) > tor )then
                        carborg_sed( ks ) = carborg( act, ks )
                    else
                        carborg_sed( ks ) = 0.0_8
                    endif
                    sumdep = sumdep + carborg_sed( ks )
                enddo
               if( sumdep > tor ) call add_carborg_deposit( n, carborg_sed, sumdep )
            endif

            ! Erosion
            act = 2
            if( action_flag( act ) )then
                sumdep = 0
                do ks = 1, totgrn
                    if( carborg( act, ks ) > tor )then
                        carborg_sed( ks ) = carborg( act, ks )
                    else
                        carborg_sed( ks ) = 0.0_8
                    endif
                    sumdep = sumdep + carborg_sed( ks )
                enddo
                if( sumdep > tor ) call erode_carborg_deposit( n, carborg_sed )
            endif

            ! Porosity
            act = 3
            if( action_subflag( act ) )then
                ! First work on the sub-surface
                laynb = locv_NbLay( n )
                if( gporo%compaction .and. laynb > 1 )then
                    do kl = 1, laynb
                        do ks = 1, totgrn
                            ! For each considered sediment check if any fuzzy rule has been applied
                            if( strat_sedh( n, kl, ks ) > 0.0_8 .and. carborg_subsurface( act, ks, kl ) > 0.0_8 )then
                                strat_porosity( n,  kl, ks  ) = min( carborg_subsurface( act, ks, kl ), &
                                    strat_porosity( n,  kl, ks  ) )
                            endif
                        enddo
                    enddo
                endif
            endif

            ! Then check top sedimentary layer
            if( action_flag( act ) )then
                if( gporo%compaction )then
                    do ks = 1, totgrn
                        ! For each considered sediment check if any fuzzy rule has been applied
                        if( strat_sedh( n, laynb, ks ) > 0.0_8 .and. carborg( act, ks ) > 0.0_8 )then
                            if( strat_porosity( n,  laynb, ks  ) == 0.0_8 )then
                                strat_porosity( n,  laynb, ks  ) = carborg( act, ks )
                            else
                                strat_porosity( n,  laynb, ks  ) = min( carborg( act, ks ), &
                                    strat_porosity( n,  laynb, ks  ) )
                            endif
                        endif
                    enddo
                endif
            endif

            ! Hardness
            act = 4
            if( action_subflag( act ) )then
                ! First work on the sub-surface
                laynb = locv_NbLay( n )
                if( laynb > 1 )then
                    do kl = 1, laynb
                        do ks = 1, totgrn
                            ! For each considered sediment check if any fuzzy rule has been applied
                            if( strat_sedh( n, kl, ks ) > 0.0_8 .and. carborg_subsurface( act, ks, kl ) > 0.0_8 )then
                                strat_hardness(  n,  kl  ) = max( strat_hardness(  n,  kl  ), &
                                    carborg_subsurface( act, ks, kl ) )
                            endif
                        enddo
                    enddo
                endif
            endif

            ! Then check top sedimentary layer
            if( action_flag( act ) )then
                do ks = 1, totgrn
                    ! For each considered sediment check if any fuzzy rule has been applied
                    if( strat_sedh( n, laynb, ks ) > 0.0_8 .and. carborg( act, ks ) > 0.0_8 )then
                        strat_hardness(  n,  laynb  ) = max( strat_hardness(  n,  laynb  ), &
                                    carborg_subsurface( act, ks, laynb ) )
                    endif
                enddo
            endif

        enddo

        ! Update global elevation
        call mpi_allreduce( proc_elev, elev_record , G_strat%nodes_nb, &
            dbl_type, max_type, SPModel_comm_world, ierr )

        return

    end subroutine compute_carborgs_evolution
    ! ============================================================================
    !> Subroutine compute_fuzzy_logic
    !! Subroutine compute_fuzzy_logic compute carbonates growth using fuzzy logic.
    !!
    !<
    ! ============================================================================
    subroutine compute_fuzzy_logic( carborg, subsurface )

        integer :: mfn, nrl, sedID, act, subsurface, mfnID, ks

        real( tkind ) :: temp
        real( tkind ) :: val( membership_fct_nb )
        real( tkind ) :: carborg( 4, totgrn )

        ! 1st step: Fuzzyfication
        ! Return a fuzzy value ([0,1]) for each active membership function
        do mfn = 1, membership_fct_nb
            act = membership( mfn )%variableID
            if( act > 12 .and. subsurface == 0 ) goto 11
            if( act < 13 .and. subsurface == 1 ) goto 11
            if( act < 5 )then
                temp = carborg_parameter( act )
            elseif( act == 5 )then
                ks = membership( mfn )%sedclassID
                temp = carborg_parameter( 4 + ks )
            else
                temp = carborg_parameter( act + totgrn - 1 )
            endif
            val( mfn ) = fuzzify_membership_function( mfn, temp )
11          continue
        enddo

        ! Reset carbonates/organics local values
        carborg = 0.0_8
        do ks = 1, totgrn
            do act = 1, 4
                fuzzycomb( act, ks )%nb_points = 0
            enddo
        enddo

        ! Generate the value for each fuzzy rule based on membership and fuzzy functions.
        do nrl = 1, fuzzy_rl_nb

            ! Get the fuzzy function control type associated with the considered rule
            act = fuzzyfct( fuzzyfy( nrl )%fuzzyID )%variableID

            ! In case we are looking at sub-surface rules
            ! do not use surface fuzzy control
            if( act < 3 .and. subsurface == 1 ) goto 10

            ! In case we are looking at surface rules
            ! do not use subsurface fuzzy control
            if( act > 2 .and. subsurface == 0 ) goto 10

            ! 2nd step: Fuzzy combination using Zadeh T-norm operator with
            ! intersection only.
            ! Find the rule strength from the list of membership functions
            temp = 1.0_8
            do mfn = 1, fuzzyfy( nrl )%membershipNb
                mfnID = fuzzyfy( nrl )%membershipID( mfn )
                temp = min( temp, val( mfnID ) )
            enddo

            ! Get the carbonate/organic type defined for this fuzzy rule
            sedID = fuzzyfy( nrl )%sedclassID

            ! 3rd step: Combined fuzzy output distribution for each sediment type
            if( temp > 0 ) call fuzzy_distribution( act, sedID, fuzzyfy( nrl )%fuzzyID, temp )

10          continue

        enddo

        ! 4th step: Apply the defuzzification using the Centre of Area fuzzy technique
        ! which returns a crisp value for each type of sediment
        do ks = 1, totgrn
            do act = 1, 4
                carborg( act, ks ) = defuzzify_function( act, ks )
                ! For deposition/erosion type fuzzy rule convert the crisp value to metres
                if(  act <= 2 ) carborg( act, ks ) = carborg( act, ks ) * time_carborg
                if( carborg( act, ks ) > 0 ) action_flag( act ) = .true.
            enddo
        enddo

        return

    end subroutine compute_fuzzy_logic
    ! ============================================================================
    !> Function fuzzify_membership_function
    !! Function fuzzify_membership_function computes the solution of the membership
    !! function for the considered value. This corresponds to the fuzzyfication step in
    !! fuzzy logic systems.
    !!
    !<
    ! ============================================================================
    function fuzzify_membership_function( mfn, temp )  result( rslt )

        integer :: mfn, m

        real( tkind ) :: rslt, temp

        ! Scan through to find the appropriate bracket values for the computed variable
        if( temp < membership( mfn )%xypoints( 1, 1 ) )then
            rslt = membership( mfn )%xypoints( 1, 2 )
        elseif( temp > membership( mfn )%xypoints( membership( mfn )%nb_points, 1 ) )then
            rslt = membership( mfn )%xypoints( membership( mfn )%nb_points, 2 )
        else
            mbsfct: do m = 1, membership( mfn )%nb_points
                if( temp <= membership( mfn )%xypoints( m, 1 ) ) exit mbsfct
            enddo mbsfct

            ! Solve the rule by linear interpolation
            if( m > 1 .and. m < membership( mfn )%nb_points )then
                if( membership( mfn )%xypoints( m, 1 ) == membership( mfn )%xypoints( m-1, 1 ) )then
                    rslt = membership( mfn )%xypoints( m, 1 )
                else
                    rslt = membership( mfn )%xypoints( m, 2 ) - membership( mfn )%xypoints( m-1, 2 )
                    rslt = rslt * ( temp - membership( mfn )%xypoints( m-1, 1 )  )
                    rslt = rslt / ( membership( mfn )%xypoints( m, 1 ) - membership( mfn )%xypoints( m-1, 1 ) )
                    rslt = rslt + membership( mfn )%xypoints( m-1, 2 )
                endif
            else
                rslt = membership( mfn )%xypoints( m, 2 )
            endif
        endif

        if( rslt > 1.0 .or. rslt < 0.0 )then
            print*,'Something went wrong when fuzzifying the membership function.', mfn, temp, m, rslt, membership( mfn )%nb_points
            stop
        endif

    end function fuzzify_membership_function
    ! ============================================================================
    !> Function defuzzify_function
    !! Function defuzzify_function computes the solution of the fuzzy logic
    !! rule based on output distribution and centre of area technique.
    !!
    !<
    ! ============================================================================
    function defuzzify_function( act, sed )  result( rslt )

        integer :: m, act, sed

        real( tkind ) :: rslt, sum_denom, sum_num

        rslt = 0.0_8
        sum_num = 0.0_8
        sum_denom = 0.0_8

        ! Apply the centre of area technique
        if( fuzzycomb( act, sed )%nb_points > 0 )then
            do m = 1, fuzzycomb( act, sed )%nb_points
                sum_num = sum_num + fuzzycomb( act, sed )%xypoints( m, 1 ) * &
                    fuzzycomb( act, sed )%xypoints( m, 2 )
                sum_denom = sum_denom + fuzzycomb( act, sed )%xypoints( m, 2 )
            enddo
        endif

        if( sum_denom /= 0.0_8 ) rslt = sum_num / sum_denom

    end function defuzzify_function
    ! ============================================================================
    !> Subroutine fuzzy_distribution
    !! Subroutine fuzzy_distribution computes the solution of the fuzzy logic
    !! rules for each type of sediments to form the fuzzy output distribution .
    !!
    !<
    ! ============================================================================
    subroutine fuzzy_distribution( act, sed, fuzzid, temp )

        logical :: intersect

        integer :: sed, fuzzid, interpts
        integer :: m, nb, cnb, act, p, n, newp, pt, ptt, ptt2

        real( tkind ) :: temp, rsl, rslt( 2 ), nx( 40 )
        real( tkind ), dimension( 40, 2 ) :: temp_xy, clean_xy, intersect_xy
        real( tkind ), dimension( 40, 2 ) :: xy_new, new_xy

        ! Create a temporary set of points containing the clipped region
        ! for the considered fuzzy function ID
        nb = 1
        temp_xy( nb, 1 ) = fuzzyfct( fuzzid )%xypoints( 1, 1 )
        temp_xy( nb, 2 ) = fuzzyfct( fuzzid )%xypoints( 1, 2 )

        do m = 2, fuzzyfct( fuzzid )%nb_points
            if( ( fuzzyfct( fuzzid )%xypoints( m, 2 ) > temp .and. fuzzyfct( fuzzid )%xypoints( m-1, 2 ) < temp ) .or. &
                ( fuzzyfct( fuzzid )%xypoints( m, 2 ) < temp .and. fuzzyfct( fuzzid )%xypoints( m-1, 2 ) > temp ) )then
                 ! Find the x coordinate for the intersection
                 rsl = fuzzyfct( fuzzid )%xypoints( m, 1 ) - fuzzyfct( fuzzid )%xypoints( m-1, 1 )
                 rsl = rsl * ( temp - fuzzyfct( fuzzid )%xypoints( m-1, 2 )  )
                 rsl = rsl / ( fuzzyfct( fuzzid )%xypoints( m, 2 ) - fuzzyfct( fuzzid )%xypoints( m-1, 2 ) )
                 rsl = rsl + fuzzyfct( fuzzid )%xypoints( m-1, 1 )
                 ! Add new point to output fuzzy distribution
                 nb = nb + 1
                 temp_xy( nb, 1 ) = rsl
                 temp_xy( nb, 2 ) = temp
            elseif( fuzzyfct( fuzzid )%xypoints( m, 2 ) >= temp )then
                 ! Add new point to output fuzzy distribution
                 nb = nb + 1
                 temp_xy( nb, 1 ) = fuzzyfct( fuzzid )%xypoints( m, 1 )
                 temp_xy( nb, 2 ) = temp
            elseif( fuzzyfct( fuzzid )%xypoints( m, 2 ) < temp )then
                 ! Add new point to output fuzzy distribution
                 nb = nb + 1
                 temp_xy( nb, 1 ) = fuzzyfct( fuzzid )%xypoints( m, 1 )
                 temp_xy( nb, 2 ) = fuzzyfct( fuzzid )%xypoints( m, 2 )
            endif
        enddo

        ! Clean the points
        call clean_sorted_pointset( temp_xy, nb, clean_xy, cnb )

        ! Combine new points with existing ones
        if( fuzzycomb( act, sed )%nb_points == 0 )then

            fuzzycomb( act, sed )%nb_points = cnb
            do m = 1, cnb
                fuzzycomb( act, sed )%xypoints( m, 1:2 ) = clean_xy( m, 1:2 )
            enddo

        else

            ! Compute intersection points between fuzzy cliped distributions
            interpts = 0
            do p = 1, fuzzycomb( act, sed )%nb_points - 1
                do n = 1, cnb - 1
                    intersect = .false.
                    call find_intersection_between_segments( &
                        fuzzycomb( act, sed )%xypoints( p, 1:2 ), &
                        fuzzycomb( act, sed )%xypoints( p + 1, 1:2 ), &
                        clean_xy( n, 1:2 ), clean_xy( n+1, 1:2 ), intersect, rslt )
                    if( intersect )then
                        interpts = interpts + 1
                        intersect_xy( interpts, 1:2 ) = rslt( 1:2 )
                    endif
                enddo
            enddo

            ! Work with the first fuzzy output distribution
            ! Find points with a new abscissa
            call add_new_absissa(  clean_xy( 1:40, 1 ), cnb, fuzzycomb( act, sed )%xypoints( 1:40, 1 ), &
                fuzzycomb( act, sed )%nb_points, nx, newp )

            ! Insert new points and compute their fuzzy values
            if( newp > 0 ) call add_new_pts(  fuzzycomb( act, sed )%xypoints, fuzzycomb( act, sed )%nb_points, nx, &
                newp, temp_xy, pt )

            ! Add intersection points to the new set of points
            if( interpts > 0 ) call add_intersection_points( intersect_xy, interpts, temp_xy, pt, new_xy, ptt )

            ! Work with the new fuzzy output distribution
            ! Find points with a new abscissa
            call add_new_absissa(  fuzzycomb( act, sed )%xypoints( 1:40, 1 ), fuzzycomb( act, sed )%nb_points, &
                clean_xy( 1:40, 1 ), cnb, nx, newp )

            ! Insert new points and compute their fuzzy values
            if( newp > 0 ) call add_new_pts(  clean_xy, cnb, nx, newp, temp_xy, pt )

            ! Add intersection points to the new set of points
            if( interpts > 0 ) call add_intersection_points( intersect_xy, interpts, temp_xy, pt, xy_new, ptt2 )

            ! Check egality between the number of points in the two fuzzy functions
            if( ptt2 /= ptt )then
                print*,'Something went wrong when comparing number of points between the fuzzy distribution.'
                print*,ptt2,ptt
                stop
            endif

            ! Combine the two fuzzy output distributions
            if( ptt > 0 )then

                ! Combination uses the logical fuzzy "OR" operator which takes the maximum values
                do n = 1, ptt
                    if( new_xy( n, 1 ) /= xy_new( n, 1 ) )then
                        print*,'Something went wrong when checking abscissa values in both fuzzy distribution.'
                        stop
                    endif
                    temp_xy( n, 1 ) = new_xy( n, 1 )
                    temp_xy( n, 2 ) = max( new_xy( n, 2 ), xy_new( n, 2 ) )
                enddo

                ! Clean the points
                call clean_sorted_pointset( temp_xy, ptt, clean_xy, cnb )

                ! Set the combined fuzzy output distribution
                fuzzycomb( act, sed )%nb_points = cnb
                do n = 1, cnb
                    fuzzycomb( act, sed )%xypoints( n, 1:2 ) = clean_xy( n, 1:2 )
                enddo

            endif

        endif


        return

    end subroutine fuzzy_distribution
    ! ============================================================================
    !> Subroutine add_new_absissa
    !! Subroutine add_new_absissa adds points with new absissa in respect to another
    !! distribution.
    !!
    !<
    ! ============================================================================
    subroutine add_new_absissa(  addSet, addPts, compSet, compPts, newX, newPts )

        logical :: add_the_point

        integer :: addPts, compPts, newPts, n, p

        real( tkind ) :: addSet( 40 ), compSet( 40 ), newX( 40 )

        newPts = 0
        do n = 1, addPts
            add_the_point = .true.
            do p = 1, compPts
                if( addSet( n ) == compSet( p ) ) add_the_point = .false.
            enddo
            if( add_the_point )then
                newPts = newPts + 1
                newX( newPts ) = addSet( n )
            endif
        enddo

        return

    end subroutine add_new_absissa
    ! ============================================================================
    !> Subroutine add_new_pts
    !! Subroutine add_new_pts adds points in respect to another distribution.
    !!
    !<
    ! ============================================================================
    subroutine add_new_pts(  setXY, setPts, addX, addPts, new_XY, newPts )

        integer :: n, p, newPts, addPts, setPts

        real( tkind ), dimension( 40 ) :: addX
        real( tkind ), dimension( 40, 2 ) :: setXY, new_XY


        newPts = 0
        do n = 1, addPts
            if( addX( n ) < setXY( 1, 1 ) )then
                newPts = newPts + 1
                new_XY( newPts, 1 ) = addX( n )
                new_XY( newPts, 2 ) = setXY( 1, 2 )
            elseif( addX( n ) < setXY( setPts, 1 ) )then
                do p = 1, setPts - 1
                    if( addX( n ) > setXY( p, 1 ) .and. addX( n ) <  setXY( p+1, 1 ) )then
                        newPts = newPts + 1
                        new_XY( newPts, 1 ) = addX( n )
                        new_XY( newPts, 2 ) = ( setXY( p + 1, 2 ) - setXY( p, 2 ) ) * ( addX( n ) - setXY( p, 1 ) ) &
                            / ( setXY( p + 1, 1 ) - setXY( p, 1 ) ) + setXY( p, 2 )
                    endif
                enddo
            elseif( addX( n ) > setXY( setPts, 1 ) )then
                newPts = newPts + 1
                new_XY( newPts, 1 ) = addX( n )
                new_XY( newPts, 2 ) = setXY( setPts, 2 )
            endif
        enddo

        if( newPts /= addPts + setPts )then
            print*,'Something went wrong went summing number of points in fuzzy distribution.'
            print*,newPts, addPts, setPts
            stop
        endif

        return

    end subroutine add_new_pts
    ! ============================================================================
    !> Subroutine add_intersection_points
    !! Subroutine add_intersection_points adds the intersection points to the active set of
    !! points.
    !!
    !<
    ! ============================================================================
    subroutine add_intersection_points( interXY, interPts, tempXY, tempPts, newXY, newPts )

        integer :: interPts, tempPts, newPts, p, n

        real( tkind ), dimension( 40, 2 ) :: interXY, tempXY, newXY

        newPts = 0
        do p = 1, interPts
            do n = 1, tempPts - 1
                if( interXY( p, 1 ) > tempXY( n+1, 1 ) )then
                    newPts = newPts + 1
                    newXY( newPts, 1:2 ) = tempXY( n, 1:2 )
                elseif( interXY( p, 1 ) > tempXY( n, 1 ) .and. interXY( p, 1 ) < tempXY( n + 1, 1 ) )then
                    newPts = newPts + 1
                    newXY( newPts, 1:2 ) = interXY( p, 1:2 )
                elseif( interXY( p, 1 ) < tempXY( n, 1 ) )then
                    newPts = newPts + 1
                    newXY( newPts, 1:2 ) = tempXY( n, 1:2 )
                endif
            enddo
        enddo
        if( newPts /= interPts + tempPts )then
            print*,'Something went wrong when adding intersecting points to the fuzzy distribution.'
            print*,newPts, interPts, tempPts
            stop
        endif
        if( newPts > 40 )then
            print*,'Something went wrong when adding intersecting points to fuzzy distribution.'
            print*,'The number of points in the distribution now exceeds 40.'
            stop
        endif

        return

    end subroutine add_intersection_points
    ! ============================================================================
    !> Subroutine clean_sorted_pointset
    !! Subroutine clean_sorted_pointset looks at a set of points and clean the unnecessary
    !! points.
    !!
    !<
    ! ============================================================================
    subroutine clean_sorted_pointset( tempset, temp_pts, cleanset, clean_pts )

        integer :: clean_pts, temp_pts, m

        real( tkind ) :: temp
        real( tkind ), dimension( 40, 2 ) :: cleanset, tempset

        cleanset( 1, 1:2 ) = tempset( 1, 1:2 )
        clean_pts = 1
        do m = 2, temp_pts - 1
            temp = tempset( m, 2 )
            if( tempset( m - 1, 2 ) == temp .and. tempset( m + 1, 2 ) == temp  )then
                goto 13
            else
                clean_pts = clean_pts + 1
                cleanset( clean_pts, 1:2 ) = tempset( m, 1:2 )
            endif
13      continue
        enddo
        clean_pts = clean_pts + 1
        cleanset( clean_pts, 1:2 ) = tempset( temp_pts, 1:2 )

        if( clean_pts > 40 )then
            print*,'Something went wrong when summing number of points in fuzzy distribution.'
            print*,'The number of points in the distribution exceeds 40.'
            stop
        endif

        return

    end subroutine clean_sorted_pointset
    ! ============================================================================

end module mod_buildcarbs
