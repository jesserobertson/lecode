! ============================================================================
! Name        : EDTransfer.f90
! Author      : tristan salles
! Created on: April 02, 2014
! Copyright (C) 2014 CSIRO
! ============================================================================
!> \file EDTransfer.f90
!!
!! EDTransfer transfer flow walkers for erosion deposition calculation.
!!
!<
! ============================================================================
module mod_edtransfer

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
    use mod_rungekutta
    use mod_direction
    use mod_slopefw
    use sediment_data

    implicit none

    type flowTransfer
        integer :: procid = -1
        integer :: refid = 0
        integer :: bound = 0
        integer :: type = 0
        real(tkind) :: zpos = 0
        real(tkind) :: height = 0
        real(tkind) :: volume = 0
        real(tkind) :: density = 0
        real(tkind) :: velocity = 0
        integer :: AInb = 0
        integer, dimension(20) :: AIpts
        real(tkind), dimension(max_grn) :: sedcharge
    end type flowTransfer

    type flowWalker
        integer :: refid = 0
        integer :: inbound = 0
        integer :: count = 0
        integer :: type = 0
        real(tkind) :: xpos = 0
        real(tkind) :: ypos = 0
        real(tkind) :: zpos = 0
        real(tkind) :: xslp = 0
        real(tkind) :: yslp = 0
        real(tkind) :: pgrad = 0
        real(tkind) :: height = 0
        real(tkind) :: width = 0
        real(tkind) :: xvel = 0
        real(tkind) :: yvel = 0
        real(tkind) :: volume = 0
        real(tkind) :: density = 0
        real(tkind) :: pzfe = 0
        real(tkind) :: pcoeff = 0
        real(tkind) :: pdist = 0
        real(tkind) :: pvelx = 0
        real(tkind) :: pvely = 0
        integer :: AInb = 0
        integer, dimension(20) :: AIpts
        real(tkind), dimension(max_grn) :: sedcharge
        integer, dimension(numold) :: faceid
    end type flowWalker

    integer :: mpi_flowwalker, mpi_flowtransfer

    public

contains

    ! ============================================================================
    !> Subroutine exchange_flow_walkers
    !! Subroutine exchange_flow_walkers checks each local flow walker for location, and sends
    !! it to another process if necessary.
    !<
    ! ============================================================================
    subroutine exchange_flow_walkers

        integer :: k, m, p, gid, temp
        integer :: faceLocation, nodeLocation
        integer :: outgoingCount, incomingCount
        integer :: outgoingCountn, incomingCountn

        integer, dimension( nproc ) :: sendCounts, sendOffsets, sendIdx
        integer, dimension( nproc ) :: recvCounts, recvOffsets
        integer, dimension( num_fw ) :: sendTargets

        integer, dimension( nproc ) :: sendCountsn, sendOffsetsn, sendIdxn
        integer, dimension( nproc ) :: recvCountsn, recvOffsetsn
        integer, dimension( tot_fw ) :: sendTargetsn

        type( flowWalker ), dimension(:), allocatable :: outgoingFWs
        type( flowWalker ), dimension(:), allocatable :: incomingFWs
        type( flowWalker ) :: tempFW

        type( flowTransfer ), dimension(:), allocatable :: outgoingFWsn
        type( flowTransfer ), dimension(:), allocatable :: incomingFWsn
        type( flowTransfer ) :: tempFWn

        ! Run over all flow walkers, and check their location on the stratigraphy.
        ! Then work out which processor they belong to, keeping track of how many will
        ! be going to that processor in total.
        sendCounts = 0
        sendTargets = -1
        sendCountsn = 0
        sendTargetsn = -1
        do k = 1, num_fw

            if( fa_nbPtAI( k ) > 20 ) print*,'Something went wrong allocation ero/dep points number.'

            if( fa_inbound( k ) /= 1 )then
                ! Extract the FW's face and node.
                call closest_stratal_node( fa_xpos( k ), fa_ypos( k ), faceLocation, nodeLocation )

                ! If the FW belongs somewhere else, take note.
                if( gf_pid( faceLocation ) /= iam )then
                    sendTargets( k ) = gf_pid( faceLocation )
                    sendCounts( gf_pid( faceLocation ) + 1 ) = sendCounts( gf_pid( faceLocation ) + 1 ) + 1
                endif

                ! Loop over nodes processors to find possible transfer ones
                ! No need to transfer the information if the processor already contains
                ! the information and if the fw is going to be transfer to the considered
                ! partition where the node belongs.
                do m = 1, nproc

                    ! Loop over the number of nodes in the area of influence of the flow walker
                    do p = 1, fa_nbPtAI( k )

                        ! Get influenced points global ID
                        gid = fa_ptsAIfw( k , p )

                        ! If the nodes belongs somewhere else, and if the nodes doesn't sit in the
                        ! new transfer partition, record it for transfer.
                        if( gv_SharenID( gid, m ) > 0 .and. sendTargets( k ) /= m - 1 .and. &
                            iam /= m - 1 .and. sendTargetsn( k ) == -1 )then
                            sendTargetsn( k ) = m - 1
                            sendCountsn( m ) = sendCountsn( m ) + 1
                        elseif( gv_SharenID( gid, m ) > 0 .and. iam == m - 1 )then
                            ! Populate FWs for local computation
                            ed_procid( k ) = gf_pid( faceLocation )
                            ed_refid( k ) = fa_refid( k )
                            ed_bound( k ) = fa_inbound( k )
                            ed_type( k ) = fa_type( k )
                            ed_zpos( k ) = fa_zpos( k )
                            ed_h( k ) = fa_h( k )
                            ed_vel( k ) = sqrt( fa_xvel( k )**2 + fa_yvel( k )**2 )
                            ed_vol( k ) = fa_volume( k )
                            ed_dens( k ) = fa_density( k )
                            ed_AInb( k ) = fa_nbPtAI( k )
                            ed_AIpts( k, 1:fa_nbPtAI( k ) ) = fa_ptsAIfw( k, 1:fa_nbPtAI( k ) )
                            ed_sed( k, 1:totgrn ) = fa_sedcharge( k,1:totgrn )
                        endif
                    enddo
                enddo
            endif
        enddo

        ! Count how many particles I'm sending, and allocate my outgoing FW buffer.
        outgoingCount = 0
        outgoingCountn = 0
        do k = 1, nproc
            outgoingCount = outgoingCount + sendCounts( k )
            outgoingCountn = outgoingCountn + sendCountsn( k )
        enddo
        allocate( outgoingFWs( outgoingCount ) )
        allocate( outgoingFWsn( outgoingCountn ) )

        ! Make sure that every processor knows how many FWs it will be receiving from
        ! every other processor.
        recvCounts = 0
        recvCountsn = 0
        call mpi_alltoall( sendCounts, 1, int_type, recvCounts, &
            1, int_type, SPModel_comm_world, ierr )
        call mpi_alltoall( sendCountsn, 1, int_type, recvCountsn, &
            1, int_type, SPModel_comm_world, ierr )

        ! Count how many new FWs I'll be handled, calculate offsets,
        ! and allocate my incoming FW buffer.
        incomingCount = 0
        recvOffsets = 0
        incomingCountn = 0
        recvOffsetsn = 0
        do k = 1, nproc
            ! Set the offset for data flowing from this process.
            recvOffsets( k ) = incomingCount
            recvOffsetsn( k ) = incomingCountn

            ! Incorporate the incoming particles into our count.
            incomingCount = incomingCount + recvCounts( k )
            incomingCountn = incomingCountn + recvCountsn( k )

            ! Sanity-check.
            if( k == iam + 1 .and. ( recvCounts( k ) /= 0 .or. &
                recvCountsn( k ) /= 0 ) )then
                print *, "Something went wrong: trying to receive particles from myself!"
                stop
            endif
        enddo
        allocate( incomingFWs( incomingCount ) )
        allocate( incomingFWsn( incomingCountn ) )

        ! Set up our outgoing offset array.
        sendOffsets = 0
        outgoingCount = 0
        sendOffsetsn = 0
        outgoingCountn = 0
        do k = 1, nproc
            sendOffsets( k ) = outgoingCount
            outgoingCount = outgoingCount + sendCounts( k )
            sendOffsetsn( k ) = outgoingCountn
            outgoingCountn = outgoingCountn + sendCountsn( k )
        end do

        ! Marshall the FWs into the output array.
        sendIdx = sendOffsets
        sendIdxn = sendOffsetsn
        do k = 1, num_fw

            if( sendTargetsn( k ) >= 0 )then

                ! Move our offset along by one.
                sendIdxn(sendTargetsn(k)+1) = sendIdxn(sendTargetsn(k)+1) + 1

                ! Populate the FW datatype.
                tempFWn%procid = iam
                tempFWn%refid = fa_refid( k )
                tempFWn%bound = fa_inbound( k )
                tempFWn%type = fa_type( k )
                tempFWn%zpos = fa_zpos( k )
                tempFWn%height = fa_h( k )
                tempFWn%volume = fa_volume( k )
                tempFWn%density = fa_density( k )
                tempFWn%velocity = sqrt( fa_xvel( k )**2 + fa_yvel( k )**2 )
                tempFWn%AInb = fa_nbPtAI( k )
                tempFWn%AIpts = 0
                tempFWn%AIpts( 1:fa_nbPtAI( k ) ) = fa_ptsAIfw( k, 1:fa_nbPtAI( k ) )
                tempFWn%sedcharge = 0.0
                tempFWn%sedcharge( 1:totgrn ) = fa_sedcharge( k, 1:totgrn )

                ! Copy the FW into the marshalling area.
                outgoingFWsn( sendIdxn( sendTargetsn( k ) + 1 ) ) = tempFWn

            endif

            if( sendTargets( k ) >= 0 )then

                ! Move our offset along by one.
                sendIdx( sendTargets( k ) + 1 ) = sendIdx( sendTargets( k ) + 1 ) + 1

                ! Populate the FW datatype.
                tempFW%refid = fa_refid( k )
                tempFW%inbound = fa_inbound( k )
                tempFW%count = fa_count( k )
                tempFW%type = fa_type( k )
                tempFW%xpos = fa_xpos( k )
                tempFW%ypos = fa_ypos( k )
                tempFW%zpos = fa_zpos( k )
                tempFW%xslp = fa_xslp( k )
                tempFW%yslp = fa_yslp( k )
                tempFW%pgrad = fa_pgrad( k )
                tempFW%height = fa_h( k )
                tempFW%width = fa_width( k )
                tempFW%xvel = fa_xvel( k )
                tempFW%yvel = fa_yvel( k )
                tempFW%volume = fa_volume( k )
                tempFW%density = fa_density( k )
                tempFW%pzfe = fa_pzfe( k )
                tempFW%pcoeff = fa_pcoeff( k )
                tempFW%pdist = fa_pdist( k )
                tempFW%pvelx = fa_pvelx( k )
                tempFW%pvely = fa_pvely( k )
                tempFW%AInb = fa_nbPtAI( k )
                tempFW%AIpts = 0
                tempFW%AIpts( 1:fa_nbPtAI( k ) ) = fa_ptsAIfw( k, 1:fa_nbPtAI( k ) )
                tempFW%sedcharge = 0.0
                tempFW%sedcharge( 1:totgrn ) = fa_sedcharge( k, 1:totgrn )
                tempFW%faceid( 1:numold ) = fa_faceid( k, 1:numold )

                ! Mark it for removal from this processor.
                fa_inbound( k ) = 1

                ! Copy the FW into the marshalling area.
                outgoingFWs( sendIdx( sendTargets( k ) + 1 ) ) = tempFW

            endif

        enddo

        ! Exchange FWs.
        call mpi_alltoallv( outgoingFWs, sendCounts, sendOffsets, mpi_flowwalker, &
            incomingFWs, recvCounts, recvOffsets, mpi_flowwalker, SPModel_comm_world, ierr )
        call mpi_alltoallv( outgoingFWsn, sendCountsn, sendOffsetsn, mpi_flowtransfer, &
            incomingFWsn, recvCountsn, recvOffsetsn, mpi_flowtransfer, SPModel_comm_world, ierr )

        ! Unpack our received particles into our FW array.
        temp = num_fw + 1
        do k = 1, incomingCount
            fa_refid( temp ) = incomingFWs( k )%refid
            fa_inbound( temp ) = incomingFWs( k )%inbound
            fa_count( temp ) = incomingFWs( k )%count
            fa_type( temp ) = incomingFWs( k )%type
            fa_xpos( temp ) = incomingFWs( k )%xpos
            fa_ypos( temp ) = incomingFWs( k )%ypos
            fa_zpos( temp ) = incomingFWs( k )%zpos
            fa_xslp( temp ) = incomingFWs( k )%xslp
            fa_yslp( temp ) = incomingFWs( k )%yslp
            fa_pgrad( temp ) = incomingFWs( k )%pgrad
            fa_h( temp ) = incomingFWs( k )%height
            fa_width( temp ) = incomingFWs( k )%width
            fa_xvel( temp ) = incomingFWs( k )%xvel
            fa_yvel( temp ) = incomingFWs( k )%yvel
            fa_volume( temp ) = incomingFWs( k )%volume
            fa_density( temp ) = incomingFWs( k )%density
            fa_pzfe( temp ) = incomingFWs( k )%pzfe
            fa_pcoeff( temp ) = incomingFWs( k )%pcoeff
            fa_pdist( temp ) = incomingFWs( k )%pdist
            fa_pvelx( temp ) = incomingFWs( k )%pvelx
            fa_pvely( temp ) = incomingFWs( k )%pvely
            fa_nbPtAI( temp ) = incomingFWs( k )%AInb
            fa_ptsAIfw( temp, 1:fa_nbPtAI( temp ) ) = incomingFWs( k )%AIpts( 1:fa_nbPtAI( temp ) )
            fa_sedcharge( temp,1:totgrn ) = incomingFWs( k )%sedcharge( 1:totgrn )
            fa_faceid( temp,1:numold ) = incomingFWs( k )%faceid

            ! Populate FWs for local erosion/deposition computation
            ed_procid( temp ) = iam
            ed_refid( temp ) = fa_refid( temp )
            ed_bound( temp ) = fa_inbound( temp )
            ed_type( temp ) = fa_type( temp )
            ed_zpos( temp ) = fa_zpos( temp )
            ed_h( temp ) = fa_h( temp )
            ed_vel( temp ) = sqrt( fa_xvel( temp )**2 + fa_yvel( temp )**2 )
            ed_vol( temp ) = fa_volume( temp )
            ed_dens( temp ) = fa_density( temp )
            ed_AInb( temp ) = fa_nbPtAI( temp )
            ed_AIpts( temp, 1:fa_nbPtAI( temp ) ) = fa_ptsAIfw( temp, 1:fa_nbPtAI( temp ) )
            ed_sed( temp, 1:totgrn ) = fa_sedcharge( temp,1:totgrn )

            temp = temp + 1
        enddo
        num_fw = num_fw + incomingCount

        ! Unpack our received particles into our FW array.
        temp = num_fw + 1
        do k = 1, incomingCountn
            ed_procid( temp ) = incomingFWsn( k )%procid
            ed_refid( temp ) = incomingFWsn( k )%refid
            ed_bound( temp ) = incomingFWsn( k )%bound
            ed_type( temp ) = incomingFWsn( k )%type
            ed_zpos( temp ) = incomingFWsn( k )%zpos
            ed_h( temp ) = incomingFWsn( k )%height
            ed_vel( temp ) = incomingFWsn( k )%velocity
            ed_vol( temp ) = incomingFWsn( k )%volume
            ed_dens( temp ) = incomingFWsn( k )%density
            ed_AInb( temp ) = incomingFWsn( k )%AInb
            ed_AIpts( temp, 1:ed_AInb( temp ) ) = incomingFWsn( k )%AIpts( 1:ed_AInb( temp ) )
            ed_sed( temp, 1:totgrn ) = incomingFWsn( k )%sedcharge( 1:totgrn )

            temp = temp + 1
        enddo
        num_fwed = num_fw + incomingCountn

        ! Check for duplicates, and mark for deletion if found.
!        temp = 0
!        do m = 1, num_fw
!            if( fa_inbound( m ) /= 1 )then
!                temp = temp + 1
!                do k = m + 1 , num_fw
!                    if( fa_refid( m ) == fa_refid( k ) ) fa_inbound( k ) = 1
!                enddo
!            endif
!        enddo

            call mpi_barrier( SPModel_comm_world,ierr )
        return

    end subroutine exchange_flow_walkers
    ! ============================================================================
    !> Subroutine points_area_of_influence
    !! Subroutine points_area_of_influence is used to compute the nodes in the area of
    !! influence of each FW.
    !<
    ! ============================================================================
    subroutine points_area_of_influence

        integer :: k, fcS, nmax, gid, l, m, nk, p, idS, linePts, step, topPt, botPt

        real( tkind ) :: rectx( 4 ), recty( 4 ), fwelev

        integer, dimension( : ), allocatable :: GdIdx, Gdbot, Gdtop

        do k = 1, num_fw

            ! Initialise points in the region of influence of the flow walker
            fa_ptsAIfw( k ,1:maxPtsFW) = -1

            if( fa_inbound( k ) /= 1 )then

                ! Define top flow walker elevation
                fwelev = fa_h( k ) + fa_zpos( k )

                ! Get the flow walker's closest stratal point
                call closest_stratal_node( fa_xpos( k ), fa_ypos( k ), fcS, idS )

                ! Find maximum number of points in the area of influence based
                ! on flow walkers width
                step = int( 0.5_8 * fa_width( k ) / strat_dx ) + 1
                nmax = ( step + step + 1 )**2

                ! Allocate local arrays
                allocate( GdIdx( nmax ) )
                allocate( Gdbot( nmax ) )
                allocate( Gdtop( nmax ) )
                GdIdx = -1
                Gdbot = -1
                Gdtop = -1

                ! First ID is the closest point
                GdIdx( 1 ) = idS

                ! Get the IDs of the point around closest point ID using
                ! 'step' value
                p = 1
                ! Right side nodes
                m = 1
                do nk = 1, step
                    if( gv_coord( GdIdx( m ), 1 ) < strat_xo + ( strat_X - 1 ) * strat_dx )then
                        p = p + 1
                        GdIdx( p ) = GdIdx( m ) + 1
                        m = p
                    endif
                enddo

                ! Left side nodes
                m = 1
                do nk = 1, step
                    if( gv_coord( GdIdx( m ), 1 ) > strat_xo )then
                        p = p + 1
                        GdIdx( p ) = GdIdx( m ) - 1
                        m = p
                    endif
                enddo
                linePts = p

                ! Top nodes
                p = 1
                m = linePts
                Gdtop( 1:linePts ) = GdIdx( 1:linePts )
                lp00: do nk = 1, step
                    if( Gdtop( p ) > 0 .and. Gdtop( p ) <= G_strat%nodes_nb )then
                        if( gv_coord( Gdtop( p ), 2 ) <= strat_yo + ( strat_Y - 1 ) * strat_dx )then
                            do l = 1, linePts
                                Gdtop( l + m ) = Gdtop( p ) + strat_X
                                p = p + 1
                            enddo
                            m = m + linePts
                         else
                            exit lp00
                         endif
                    endif
                enddo lp00
                topPt = m

                ! Bottom nodes
                p = 1
                m = linePts
                Gdbot( 1:linePts ) = GdIdx( 1:linePts )
                lp01: do nk = 1, step
                    if( Gdbot( p ) > 0 .and. Gdbot( p ) <= G_strat%nodes_nb )then
                        if( gv_coord( Gdbot( p ), 2 ) >= strat_yo )then
                            do l = 1, linePts
                                Gdbot( l + m ) = Gdbot( p ) - strat_X
                                p = p + 1
                            enddo
                            m = m + linePts
                        else
                            exit lp01
                        endif
                    endif
                 enddo lp01
                 botPt = m

                 p = linePts
                 do nk = linePts + 1, botPt
                    if( Gdbot( nk ) > 0 .and. Gdbot( nk ) < G_strat%nodes_nb )then
                        p = p + 1
                        GdIdx( p ) = Gdbot( nk )
                    endif
                 enddo
                 do nk = linePts + 1, topPt
                    if( Gdtop( nk ) > 0 .and. Gdtop( nk ) < G_strat%nodes_nb )then
                        p = p + 1
                        GdIdx( p ) = Gdtop( nk )
                    endif
                 enddo

                if( p == 0 .or. p > nmax )print*,'Something went wrong when looking in the region of influence of the FW.'
                if( p == 0 .or. p > nmax )stop

                ! We have all the points IDs now, we need to pick the ones perpendicular
                ! to the flow walker direction. First find the 4 points making the rectangle
                call build_rectangle_around_fw( k, rectx, recty )

                ! Find the nodes within the rectangular shape
                p = 0
                do nk = 1, nmax
                    gid = GdIdx( nk )
                    if( gid > 0 )then

                        ! Check if the point lies within the rectangle
                        call find_points_rectangular_shape( gv_coord( gid,  1 ), &
                            gv_coord( gid,  2 ), rectx, recty, 4, l, m )

                        ! In case this is the case
                        if( l >= 0 )then

                            ! Check its elevation to assure we are below the flow height
                            if( wdem( gid ) < fwelev )then
                                ! Store the node for erosion/deposition
                                p = p + 1
                                fa_ptsAIfw( k , p ) = gid
                            endif
                        endif
                     endif
                enddo

                ! In case there is no point just take the closest one
                if( p == 0 )then
                    p = 1
                    fa_ptsAIfw( k , p ) = idS
                endif
                fa_nbPtAI( k ) = p

                ! Deallocate local arrays
                deallocate( GdIdx, Gdbot, Gdtop )

            endif

        enddo

        return

    end subroutine points_area_of_influence
    ! ============================================================================
    !> Subroutine build_rectangle_around_fw()
    !! Function build_rectangle_around_fw is used to build the rectangle
    !<
    ! ============================================================================
    subroutine build_rectangle_around_fw( k, x, y )

        integer :: k
        real( tkind ) :: fwx, fwy, d, e, ab, bc, hyp, ex, ey, dx, dy
        real( tkind ) :: sinalpha, cosalpha, x( 4 ), y ( 4 )

        ! Position of the flow walker
        fwx = fa_xpos( k )
        fwy = fa_ypos( k )

        ! Half of rectangle's lenght
        d = fa_width( k ) * 0.5_8
        ! Half og rectangle's width
        e = strat_dx * 0.5_8

        ! In case FW's velocity is null along X axis
        if( fa_xvel( k ) == 0.0_8 )then
            x( 1 ) = fwx - d
            x( 2 ) = fwx + d
            x( 3 ) = x( 2 )
            x( 4 ) = x( 1 )
            y( 1 ) = fwy - e
            y( 2 ) = y ( 1 )
            y( 3 ) = fwy + e
            y( 4 ) = y ( 3 )
        ! In case FW's velocity is null along Y axis
        elseif( fa_yvel( k ) == 0.0_8 )then
            x( 1 ) = fwx - e
            x( 2 ) = fwx + e
            x( 3 ) = x( 2 )
            x( 4 ) = x( 1 )
            y( 1 ) = fwy - d
            y( 2 ) = y ( 1 )
            y( 3 ) = fwy + d
            y( 4 ) = y ( 3 )
        ! Otherwise
        else
            ab = 1.0e8_8 * fa_xvel( k )
            bc = 1.0e8_8 * fa_yvel( k )
            hyp = sqrt( bc**2.0_8 + ab**2.0_8 )
            sinalpha = bc / hyp
            cosalpha = ab / hyp
            dx = fwx + sinalpha * d
            dy = fwy - cosalpha * d
            ex = fwx - sinalpha * d
            ey = fwy + cosalpha * d
            x( 1 ) = dx - e * cosalpha
            y( 1 ) = dy - e * sinalpha
            x( 2 ) = dx + e * cosalpha
            y( 2 ) = dy + e * sinalpha
            x( 3 ) = ex + e * cosalpha
            y( 3 ) = ey + e * sinalpha
            x( 4 ) = ex - e * cosalpha
            y( 4 ) = ey - e * sinalpha
        endif

        return

    end subroutine build_rectangle_around_fw
    ! ============================================================================

end module mod_edtransfer
