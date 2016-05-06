! ============================================================================
! Name        : FunctionCarbs.f90
! Author      : tristan salles
! Created on: April 05, 2014
! Copyright (C) 2014 CSIRO
! ============================================================================
!> \file FunctionCarbs.f90
!!
!! File FunctionCarbs is a set of functions to compute reefs evolution. The approach is
!! based on fuzzy logic, inspired by Ulf Nordlund code FUZZIM and is similar to Chris Dyt
!! implementation in Sedsim.
!!
!<
! ============================================================================
module mod_functioncarbs

    use file_data
    use mpidata
    use flux_data
    use time_data
    use error_data
    use TIN_data
    use strata_data
    use param_data
    use forces_data
    use sediment_data
    use mod_compaction
    use kdtree2_module
    use kdtree2_precision_module

    implicit none

    integer, dimension( 25, 2 ) :: stepfct

    real( tkind ), dimension( : ), allocatable :: carborg_parameter

    real( tkind ), dimension( :,: ), allocatable :: distance_mats
    real( tkind ), dimension( : ), allocatable :: distance_river, distance_shore
    real( tkind ), dimension( : ), allocatable :: aslp_carbs, nslp_carbs
    real( tkind ), dimension( :,: ), allocatable :: burial_carbs !, silprop_carbs

    public

contains

    ! ============================================================================
    !> Subroutine get_distance_to_nearest_shore_and_rivers
    !! This subroutine finds distance to nearest shore and rivers within the simulation area.
    !!
    !<
    ! ============================================================================
    subroutine get_distance_to_nearest_shore_and_rivers

        logical :: continental, marine

        integer :: n, f, k, pid( 4 ), shore_nb, river_nb, p

        integer, dimension( G_strat%faces_nb ) :: shore_ids

        real( tkind ) :: pt( 2 ), distance
        real( tkind ),dimension( :, : ), allocatable :: data

        type(kdtree2), pointer :: Ctree
        type(kdtree2_result), dimension( 1 ) :: Rslt

        ! Shoreline distance
        if( carb_shore_flag == 1 )then

            if( .not. allocated( distance_shore ) ) allocate( distance_shore( L_strat%nodes_nb ) )
            ! We use the faces to determine the shorelines points
            n = 0
            shore_nb = 0
            shore_ids = -1

            do f = 1, G_strat%faces_nb
                ! Get global face point IDs
                pid( 1:4 ) = gf_points( f, 1:4 )
                ! Check if this is a shoreline faces
                do k = 1, 4
                    if( elev_record( pid( k ) ) >= gsea%actual_sea + 0.01_8 ) continental = .true.
                    if( elev_record( pid( k ) ) < gsea%actual_sea ) marine = .true.
                enddo
                ! Record shore points iteratively
                if( continental .and. marine )then
                    n = n + 1
                    shore_nb =  shore_nb + 1
                    shore_ids( n ) = f
                endif
            enddo

            ! Now build the kdtree
            if( allocated( data ) ) deallocate( data )
            allocate( data( 2, shore_nb ) )
            do k = 1, shore_nb
                data( 1, k ) =  tinf_centroid( shore_ids( k ), 1 )
                data( 2, k ) =  tinf_centroid( shore_ids( k ), 2 )
            enddo
            Ctree => kdtree2_create( data, sort = .true., rearrange = .true. )

            ! Find the closest shoreline points for each local stratigraphic nodes
            do k = 1, L_strat%nodes_nb

                distance = ( strat_X**2.0_8 + strat_Y**2.0_8 ) * strat_dx**2.0_8

                ! Get the distance from local point to shoreline
                if( elev_record( locv_gid( k ) ) < gsea%actual_sea )then
                    pt( 1:2 ) = locv_coord( k, 1:2 )

                    ! Find closest shoreline point
                    call kdtree2_n_nearest( Ctree, pt, nn=1, results=Rslt )
                    distance = Rslt( 1 )%dis
                endif

                distance_shore( k ) = sqrt( distance )

            enddo

            ! Destroy kdtree
            call kdtree2_destroy(Ctree)
            deallocate( data )

        endif

        ! Rivers distance
        if( carb_river_flag == 1 )then

            if( .not. allocated( distance_river ) ) allocate( distance_river( L_strat%nodes_nb ) )

            ! Find the number of river sources active for the current time step
            river_nb = 0
            do n = 1, num_src
                if( tnow >= fws_tstrt( n ) .and. tnow < fws_tend( n ) ) river_nb = river_nb + 1
            enddo

            ! Loop over the river sources
            p = 0
            if( allocated( data ) ) deallocate( data )
            allocate( data( 2, river_nb ) )
            ! Find position of active rivers
            do n = 1, num_src
                if( tnow >= fws_tstrt( n ) .and. tnow < fws_tend( n ) )then
                    p = p + 1
                    data( 1, p ) = fws_xposition( n )
                    data( 2, p ) = fws_yposition( n )
                endif
            enddo

            ! Find the closest river points for each local stratigraphic nodes
            do k = 1, L_strat%nodes_nb

                distance = ( strat_X**2.0_8 + strat_Y**2.0_8 ) * strat_dx**2.0_8

                ! Get the distance from local point to river
                if( elev_record( locv_gid( k ) ) < gsea%actual_sea )then
                    pt( 1:2 ) = locv_coord( k, 1:2 )

                    ! Loop over active rivers
                    do n = 1, river_nb
                        distance = min( distance, ( pt( 1 ) - data( 1, n ) )**2.0_8 + &
                            ( pt( 2 ) - data( 2, n ) )**2.0_8 )
                    enddo
                endif
                distance_river( k ) = sqrt( distance )
            enddo

        endif

        return

    end subroutine get_distance_to_nearest_shore_and_rivers
    ! ============================================================================
    !> Subroutine get_distance_to_nearest_materials
    !! This subroutine  finds distance to nearest materials the simulation area.
    !!
    !<
    ! ============================================================================
    subroutine get_distance_to_nearest_materials( sedID )

        integer :: n, k, layers_nb, sedID
        integer :: local_carb_nb, carb_nb
        integer, dimension( nproc ) :: num_carbs, carbs_displacement

        integer, dimension( L_strat%nodes_nb ) :: temp_carb_id
        integer, dimension( : ), allocatable :: carb_id

        real( tkind ) :: pt( 2 ), distance, distance1
        real( tkind ),dimension( :, : ), allocatable :: data

        type(kdtree2), pointer :: Ctree
        type(kdtree2_result), dimension( 1 ) :: cRslt

        if( .not. allocated( distance_mats ) ) allocate( distance_mats( L_strat%nodes_nb, totgrn ) )

        ! Determine if considered material type is present in the top layer
        local_carb_nb = 0

        do n = 1, L_strat%nodes_nb
            distance_mats( n, sedID ) = 0.0_8
            layers_nb = locv_NbLay( n )
            if( layers_nb == 1 ) goto 20

           if( strat_sedh(  n, layers_nb ,  sedID  ) > 0.0_8 )then
                local_carb_nb = local_carb_nb + 1
                temp_carb_id( local_carb_nb ) = locv_gid( n )
            endif

20      continue

        enddo

        ! Send information globally for specified material
        call mpi_allreduce( local_carb_nb, carb_nb, 1, int_type, sum_type, SPModel_comm_world, ierr )
        if( carb_nb > 0 )then
            call mpi_allgather( local_carb_nb, 1, int_type, num_carbs, 1, int_type, SPModel_comm_world, ierr )
            carbs_displacement( 1 ) = 0
            do k = 2, nproc
                carbs_displacement( k ) = carbs_displacement( k - 1 ) + num_carbs( k - 1 )
            enddo
            allocate( carb_id( carb_nb ) )
            call mpi_allgatherv( temp_carb_id( 1:local_carb_nb), local_carb_nb, int_type, carb_id, num_carbs, &
                carbs_displacement, int_type, SPModel_comm_world, ierr )

            ! Now build kdtrees
            if( allocated( data ) ) deallocate( data )
            allocate( data( 2, carb_nb ) )
            do k = 1, carb_nb
                data( 1, k ) =  gv_coord( carb_id( k ), 1 )
                data( 2, k ) =  gv_coord( carb_id( k ), 2 )
            enddo
            deallocate( carb_id )
            Ctree => kdtree2_create( data, sort = .true., rearrange = .true. )


            ! Find the closest carbonate points for each local stratigraphic nodes
            do k = 1, L_strat%nodes_nb
                distance = ( strat_X**2.0_8 + strat_Y**2.0_8 ) * strat_dx**2.0_8
                distance1 = distance

                ! Get the distance from local point to closest material
                pt( 1:2 ) = locv_coord( k, 1:2 )
                call kdtree2_n_nearest( Ctree, pt, nn=1, results=cRslt )
                distance = cRslt( 1 )%dis
                distance_mats( k, sedID ) = sqrt( distance )
            enddo

            ! Destroy kdtrees and data arrays
            call kdtree2_destroy(Ctree)
            deallocate( data )
        endif

        return

    end subroutine get_distance_to_nearest_materials
    ! ============================================================================
    !> Subroutine get_top_layer_siliciclastic_proportion
    !! This subroutine returns the proportion of siliciclastic sediment present in the top
    !! layer.
    !!
    !<
    ! ============================================================================
    subroutine get_top_layer_siliciclastic_proportion

!        integer :: n, k, layers_nb
!
!        if( carb_sediment_flag == 1 )then
!            if( .not. allocated( silprop_carbs ) ) allocate( silprop_carbs( L_strat%nodes_nb, silgrn ) )
!        endif
!
!        if( carb_sediment_flag == 0 ) return
!
!        ! Loop over stratal nodes
!        do n = 1, L_strat%nodes_nb
!
!            silprop_carbs( n, 1:silgrn ) = 0.0_8
!            ! Top layer ID
!            layers_nb = locv_NbLay( n )
!            if( layers_nb == 1 ) goto 20
!
!            if( strat_thick( n, layers_nb ) > 0.0_8 )then
!                do k = 1, silgrn
!                    silprop_carbs( n, k ) =  strat_sedh(  n, layers_nb ,  k  ) / strat_thick( n, layers_nb )
!                enddo
!            endif
!20      continue
!        enddo

        return

    end subroutine get_top_layer_siliciclastic_proportion
    ! ============================================================================
    !> Subroutine get_carborg_sedimentation
    !! This subroutine returns the thickness of sediment deposited since the previous call.
    !!
    !<
    ! ============================================================================
    subroutine get_carborg_sedimentation( nid, param )

        integer :: nid, layers_nb
        real( tkind ) :: param

        ! Top layer ID
        layers_nb = locv_NbLay( nid )
        param = strat_zelev(  nid, layers_nb ) - carb_sedimentation( nid )
        param = param / time_carborg

        return

    end subroutine get_carborg_sedimentation
    ! ============================================================================
    !> Subroutine get_carborg_exposed_time
    !! This subroutine returns the time since last deposition occurs on the considered node.
    !!
    !<
    ! ============================================================================
    subroutine get_carborg_exposed_time( nid, param )

        integer :: nid, layers_nb, layid
        real( tkind ) :: param

        layers_nb = locv_NbLay( nid )
        layid = strat_layID( nid, layers_nb )
        param = ( L_strat%layersID - layid ) * time_display

        return

    end subroutine get_carborg_exposed_time
    ! ============================================================================
    !> Subroutine get_carborg_buriedlayer_time
    !! This subroutine returns the time difference between considered layer and current layer.
    !!
    !<
    ! ============================================================================
    subroutine get_carborg_buriedlayer_time( nid, layid, param )

        integer :: nid, layid, lay_id
        real( tkind ) :: param

        lay_id = strat_layID( nid, layid )
        param = ( L_strat%layersID - lay_id ) * time_display

        return

    end subroutine get_carborg_buriedlayer_time
    ! ============================================================================
    !> Subroutine get_averaged_slope_from_neighbours
    !! This subroutine calculates slope based on neighbours elevations.
    !! nslp_carbs defines position of the point on a hill (<0) or in a valley (>0).
    !! aslp_carbs defines average gradient ( >0 for most nodes below current node )
    !<
    ! ============================================================================
    subroutine get_averaged_slope_from_neighbours

        integer :: n, k, gid, pid, Xid, Yid, points_nb

        real( tkind ) :: x, y, xmin, xmax, ymin, ymax, xn, yn
        real( tkind ) :: dist, delta_z, sum_up, sum_grad

        if( carb_slope_flag == 1 )then
            if( .not. allocated( aslp_carbs ) ) allocate( aslp_carbs( L_strat%nodes_nb ) )
        endif
        if( carb_valley_flag == 1 )then
            if( .not. allocated( nslp_carbs ) ) allocate( nslp_carbs( L_strat%nodes_nb ) )
        endif

        if( carb_slope_flag == 0  .and. carb_valley_flag == 0 ) return

        xmin = strat_xo
        ymin = strat_yo
        xmax = gv_coord( G_strat%nodes_nb, 1 )
        ymax = gv_coord( G_strat%nodes_nb, 2 )

        ! Loop over stratal nodes
        do n = 1, L_strat%nodes_nb

            ! Initialise values
            points_nb = 0
            sum_up = 0.0_8
            sum_grad = 0.0_8
            if( carb_slope_flag == 1 ) aslp_carbs( n ) = 0.0_8
            if( carb_valley_flag == 1 ) nslp_carbs( n ) = 0.0_8

            ! Get the global node ID and coordinates
            gid = locv_gid( n )
            x = gv_coord( gid, 1 )
            y = gv_coord( gid, 2 )

            ! Find neighbours elevation and compute slope attributes
            do k = 1, 25
                if( k /= 13 )then
                    xn = x + strat_dx * stepfct( k, 1 )
                    yn = y + strat_dx * stepfct( k, 2 )

                    ! Distance to point
                    dist = sqrt( ( x - xn )**2.0_8 + ( y - yn )**2.0_8 )

                    ! Check if file is in the simulation domain
                    if( xn >= xmin .and. xn <= xmax .and. yn >= ymin .and. yn <= ymax )then
                        points_nb = points_nb + 1

                        ! Find point global ID
                        Xid = int( ( xn - xmin ) / strat_dx ) + 1
                        Yid = int( ( yn - ymin ) / strat_dx )
                        pid = Xid + strat_X * Yid

                        ! Get slope attributes
                        delta_z = elev_record( gid ) - elev_record( pid )
                        if( delta_z > 0.0_8 ) sum_up = sum_up - 1.0_8
                        if( delta_z < 0.0_8 ) sum_up = sum_up + 1.0_8
                        sum_grad = sum_grad + delta_z / dist

                    endif

                endif
            enddo

            if( carb_slope_flag == 1 ) aslp_carbs( n ) = sum_grad / points_nb
            if( carb_valley_flag == 1 ) nslp_carbs( n ) = sum_up / points_nb

        enddo

        return

    end subroutine get_averaged_slope_from_neighbours
    ! ============================================================================
    !> Subroutine get_sedimentary_layer_burial_depth
    !! This subroutine is used to calculate the current depth of burial of each layer.
    !!
    !<
    ! ============================================================================
    subroutine get_sedimentary_layer_burial_depth

        integer :: n, k, layers_nb

        real( tkind ) :: top_elev

        if( carb_burial_flag == 1 )then
            if( .not. allocated( burial_carbs ) ) allocate( burial_carbs( L_strat%nodes_nb, G_strat%total_layer ) )
        endif

        if( carb_burial_flag == 0 ) return

        ! Loop over stratal nodes
        do n = 1, L_strat%nodes_nb

            ! Initialise parameters
            layers_nb = locv_NbLay( n )
            burial_carbs( n, 1:layers_nb ) = 0.0_8
            top_elev = strat_zelev( n, layers_nb )

            ! Get burial depth
            do k = 1, layers_nb
                burial_carbs( n, k ) = top_elev - strat_zelev( n, k )
            enddo

        enddo

        return

    end subroutine get_sedimentary_layer_burial_depth
    ! ============================================================================
    !> Subroutine add_carborg_deposit
    !! This subroutine update the deposit layer based on carbonate growth.
    !!
    !<
    ! ============================================================================
    subroutine add_carborg_deposit( nid, carborg_sed, sumdep )

        integer :: nid, gid, layID, lb, ks

        real( tkind ) :: carborg_sed( totgrn ), sumdep, newdep, perc

        gid = locv_gid( nid )
        lb = locv_NbLay( nid )
        layID = strat_layID(  nid, lb )

        ! Keep carbonates under water
        if( carbgrn > 0 .and. strat_zelev(  nid,  lb  ) < gsea%actual_sea )then
            if( orggrn == 0 .and. sumdep + strat_zelev(  nid,  lb  ) > gsea%actual_sea )then
                newdep = gsea%actual_sea - strat_zelev(  nid,  lb  ) - 0.01_8
                perc = newdep / sumdep
                sumdep = 0.0
                do ks = silgrn + 1, silgrn + carbgrn
                    carborg_sed( ks ) = carborg_sed( ks ) * perc
                    sumdep = sumdep + carborg_sed( ks )
                enddo
            endif
        endif

        ! Create a new layer in case a current layer doesn't exist
        if( layID < L_strat%layersID )then
            lb = lb + 1
            strat_sedh( nid, lb, 1:max_grn ) = 0.0_8
            locv_NbLay( nid ) = lb
            strat_layID(  nid,  lb  ) = L_strat%layersID
            strat_thick(  nid,  lb  ) = sumdep
            strat_zelev(  nid,  lb  ) = sumdep + &
                strat_zelev(  nid,  lb - 1  )
            strat_porosity(  nid,  lb, 1:max_grn  ) = 0.0_8
            strat_hardness(  nid,  lb  ) = 1.0_8
            do ks = silgrn+1, totgrn
                strat_sedh(  nid,  lb ,  ks  ) = carborg_sed( ks )
                if( carborg_sed( ks ) > 0.0_8 .and. gporo%compaction ) &
                    strat_porosity(  nid,  lb, ks  ) = porosity( ks, 1 )
            enddo
        ! Otherwise add on an existing layer
        else
            strat_thick(  nid,  lb  ) = sumdep + &
                strat_thick(  nid,  lb )
            strat_zelev(  nid,  lb  ) = sumdep + &
                strat_zelev(  nid,  lb )
            do ks = silgrn+1, totgrn
                strat_sedh(  nid,  lb ,  ks  ) = strat_sedh(  nid,  lb ,  ks  ) + &
                    carborg_sed( ks )
                if( carborg_sed( ks ) > 0.0_8 .and. gporo%compaction .and. &
                    strat_porosity( nid,  lb, ks  ) == 0.0_8 ) &
                    strat_porosity(  nid,  lb, ks  ) = porosity( ks, 1 )
            enddo
        endif

        ! Update top elevation
        proc_elev( gid ) = strat_zelev( nid,  lb  )

        return

    end subroutine add_carborg_deposit
    ! ============================================================================
    !> Subroutine erode_carborg_deposit
    !! This subroutine update the deposit layer based on carbonate erosion.
    !!
    !<
    ! ============================================================================
    subroutine erode_carborg_deposit( nid, carborg_sed )

        integer :: nid, gid, layID, lb, ks

        real( tkind ) :: carborg_sed( totgrn ), erode_sed

        gid = locv_gid( nid )
        lb = locv_NbLay( nid )
        layID = strat_layID(  nid, lb )

        if( lb == 1 )return
        if( strat_thick(  nid ,  lb  ) < tor )then
            strat_thick(  nid ,  lb  ) = 0.0_8
            strat_zelev(  nid ,  lb  ) = strat_zelev(  nid ,  lb - 1  )
            strat_sedh( nid, lb, 1:max_grn ) = 0.0_8
            strat_porosity(  nid ,  lb, 1:max_grn  ) = 0.0_8
            locv_NbLay( nid ) = lb - 1
            lb = lb - 1
            if( lb == 1 ) return
            if( strat_thick(  nid ,  lb  ) < tor )&
                print*,'Something went wrong when eroding carbonate/organic sediment.'
        endif

        do ks = silgrn+1, totgrn

            ! Take into account the layer hardness
            erode_sed = carborg_sed( ks ) / strat_hardness(  nid ,  lb  )

            ! If erosion greater than sediment available in the layer
            ! take it all
            if(  strat_sedh(  nid ,  lb , ks  ) < erode_sed )&
                erode_sed = strat_sedh(  nid ,  lb ,  ks  )

            ! Adjust sediment layer
            strat_sedh(  nid ,  lb ,  ks  ) = -erode_sed + &
                strat_sedh(  nid ,  lb ,  ks  )
            strat_thick(  nid ,  lb  ) = -erode_sed + strat_thick(  nid ,  lb  )
            strat_zelev(  nid ,  lb  ) = -erode_sed + strat_zelev(  nid ,  lb  )
            if( strat_sedh(  nid ,  lb ,  ks  ) > 0.0_8 .and. gporo%compaction .and. &
                strat_porosity( nid ,  lb, ks  ) == 0.0_8 ) &
                strat_porosity(  nid ,  lb, ks  ) = porosity( ks, 1 )

            ! If not sufficient sediment in the layer just delete the sediment layer
            if( strat_thick(  nid ,  lb  ) < tor )then
                strat_thick(  nid ,  lb  ) = 0.0_8
                strat_zelev(  nid ,  lb  ) = strat_zelev(  nid ,  lb - 1  )
                strat_sedh( nid, lb, 1:max_grn ) = 0.0_8
                strat_porosity(  nid ,  lb, 1:max_grn  ) = 0.0_8
                locv_NbLay( nid ) = lb - 1
                lb = lb - 1
                if( lb == 1 ) return
                if( strat_thick(  nid ,  lb  ) < tor )&
                    print*,'Something went wrong when eroding carbonate/organic sediment2.'
            endif

        enddo

        ! Update top elevation
        proc_elev( gid ) = strat_zelev( nid,  lb  )


        return

    end subroutine erode_carborg_deposit
    ! ============================================================================
    !> Subroutine allocate_step_function
    !! This subroutine calculates step function to look at neighbours around 2 cells
    !! distance from current node
    !!
    !<
    ! ============================================================================
    subroutine allocate_step_function

        ! Define step function around 2 cells distance
        stepfct( 1, 1:2 ) = -2

        stepfct( 2, 1 ) = -1
        stepfct( 2, 2 ) = -2

        stepfct( 3, 1 ) = 0
        stepfct( 3, 2 ) = -2

        stepfct( 4, 1 ) = 1
        stepfct( 4, 2 ) = -2

        stepfct( 5, 1 ) = 2
        stepfct( 5, 2 ) = -2

        stepfct( 6, 1 ) = -2
        stepfct( 6, 2 ) = -1

        stepfct( 7, 1 ) = -1
        stepfct( 7, 2 ) = -1

        stepfct( 8, 1 ) = 0
        stepfct( 8, 2 ) = -1

        stepfct( 9, 1 ) = 1
        stepfct( 9, 2 ) = -1

        stepfct( 10, 1 ) = 2
        stepfct( 10, 2 ) = -1

        stepfct( 11, 1 ) = -2
        stepfct( 11, 2 ) = 0

        stepfct( 12, 1 ) = -1
        stepfct( 12, 2 ) = 0

        stepfct( 13, 1 ) = 0
        stepfct( 13, 2 ) = 0

        stepfct( 14, 1 ) = 1
        stepfct( 14, 2 ) = 0

        stepfct( 15, 1 ) = 2
        stepfct( 15, 2 ) = 0

        stepfct( 16, 1 ) = -2
        stepfct( 16, 2 ) = 1

        stepfct( 17, 1 ) = -1
        stepfct( 17, 2 ) = 1

        stepfct( 18, 1 ) = 0
        stepfct( 18, 2 ) = 1

        stepfct( 19, 1 ) = 1
        stepfct( 19, 2 ) = 1

        stepfct( 20, 1 ) = 2
        stepfct( 20, 2 ) = 1

        stepfct( 21, 1 ) = -2
        stepfct( 21, 2 ) = 2

        stepfct( 22, 1 ) = -1
        stepfct( 22, 2 ) = 2

        stepfct( 23, 1 ) = 0
        stepfct( 23, 2 ) = 2

        stepfct( 24, 1 ) = 1
        stepfct( 24, 2 ) = 2

        stepfct( 25, 1:2 ) = 2

        return

    end subroutine allocate_step_function
    ! ============================================================================

end module mod_functioncarbs
