! ============================================================================
! Name        : ODE_solver.f90
! Author      : tristan salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file ODE_solver.f90
!!
!! ODE_solver describes the method used to solve the ODE for flow walkers evolution
!! The method is based on an adaptation of the Cash Karp Runge Kutta method with automatic time stepping.
!!
!<
! ============================================================================
module mod_rungekutta

    use mpidata
    use time_data
    use error_data
    use fwalker_data
    use param_data
    use strata_data

    implicit none

    public

    real(tkind), parameter :: accuracy = 1.e5_8
    real(tkind), parameter :: a2 = 0.2_8
    real(tkind), parameter :: a3 = 0.3_8
    real(tkind), parameter :: a4 = 0.6_8
    real(tkind), parameter :: a5 = 1._8
    real(tkind), parameter :: a6 = 0.875_8
    real(tkind), parameter :: c1 = 37._8/378._8
    real(tkind), parameter :: c3 = 250._8/621._8
    real(tkind), parameter :: c4 = 125._8/594._8
    real(tkind), parameter :: c6 = 512._8/1771._8
    real(tkind), parameter :: d1 = 2825._8/27648._8
    real(tkind), parameter :: d3 = 18575._8/48384._8
    real(tkind), parameter :: d4 = 13525._8/55296._8
    real(tkind), parameter :: d5 = 277._8/14336._8
    real(tkind), parameter :: d6 = 0.25_8

    real(tkind), parameter :: b21 = 0.2_8
    real(tkind), parameter :: b31 = 0.075_8
    real(tkind), parameter :: b32 = 0.225_8
    real(tkind), parameter :: b41 = 0.3_8
    real(tkind), parameter :: b42 = -0.9_8
    real(tkind), parameter :: b43 = 1.2_8
    real(tkind), parameter :: b51 = -11._8/54._8
    real(tkind), parameter :: b52 = 2.5_8
    real(tkind), parameter :: b53 = -70._8/27._8
    real(tkind), parameter :: b54 = 35._8/27._8
    real(tkind), parameter :: b61 = 1631._8/55296._8
    real(tkind), parameter :: b62 = 175._8/512._8
    real(tkind), parameter :: b63 = 575._8/13824._8
    real(tkind), parameter :: b64 = 44275._8/110592._8
    real(tkind), parameter :: b65 = 253._8/4096._8

contains

    ! ============================================================================
    !> Function runge_kutta
    !! Function runge_kutta determines the next time step.
    !<
    ! ============================================================================
    function runge_kutta( ) result( dtnew )

        real( tkind ) :: dtnew, err

        err = 3.0_8

        errloop: do
            if( err <= 1.0_8 ) exit errloop
            call Cash_Karp_RKmethod( err )
            if( err == 0 )then
                dtnew = dble( 5.0_8 * dt )
            else
                if( err > 1.0_8) then
                    dtnew = dble( 0.90_8 * dt * err**(-0.25_8) )
                    if( dtnew < 0.2_8 * dt ) dtnew = dble( 0.2_8 * dt )
                    dt = dtnew
                else
                    dtnew = dble( 0.90_8 * dt * err**(-0.20_8) )
                    if( dtnew > 5.0_8 * dt ) dtnew = dble( 5.0_8 * dt )
                endif
            endif
            if( dtnew < 0 )attempt = DTNNEG
        enddo errloop

        call completion

    end function runge_kutta
    ! ============================================================================
    !> Function plume_step
    !! Function plume_step determines the next time step for plume.
    !<
    ! ============================================================================
    function plume_step( ) result( dtplume )

        integer :: k

        real( tkind ) :: dtplume, dtp, frac, fv

        dtp = 1.0e8_8

        do k = 1, num_fw
            if( fa_type( k ) >= tpplum1 )then

                ! Zone of Flow Establishement (ZFE)
                if( fa_pzfe( k ) > fa_pdist( k ) )then
                    fv = sqrt( fa_xvel( k )**2 + fa_yvel( k )**2 )
                    if( fv > 0.0_8 ) dtp = min( dtp, 0.5_8 * strat_dx / fv )
                    fa_type( k ) = tpplum1

                ! Zone of Established Flow (ZEF)
                else
                    fa_type( k ) = tpplum2
                    frac = ( fa_pdist( k ) - fa_pzfe( k ) ) / ( plume_lenght )
                    if( frac > 1.0_8 ) frac = 1.0_8
                    fv = sqrt( ( fa_xvel( k )+ frac * ocean_flow( 1 ) )**2 + &
                        ( fa_yvel( k ) + frac * ocean_flow( 2 ) )**2 )
                    if( fv > 0.0_8 ) dtp = min( dtp, 0.5_8 * strat_dx / fv )
                endif
            endif
        enddo

        call mpi_allreduce( dtp, dtplume, 1, dbl_type, min_type, SPModel_comm_world,ierr)

    end function plume_step
    ! ============================================================================
    !> Subroutine Cash_Karp_RKmethod
    !! Subroutine Cash_Karp_RKmethod utilises an adaptation of the Cash Karp Runge Kutta method
    !! with automatic time stepping.
    !! \param all_err
    !<
    ! ============================================================================
    subroutine Cash_Karp_RKmethod( all_err )

        integer :: k

        real( tkind ) :: all_err, err
        real( tkind ) :: aqz, vxt, vyt
        real( tkind ) :: k1x, k2x, k3x, k4x, k5x, k6x
        real( tkind ) :: k1y, k2y, k3y, k4y, k5y, k6y
        real( tkind ) :: xer, yer, xva, yva, xvb, yvb
        real( tkind ) :: frx6, fry6, frx4, fry4, frx4a, fry4a

        all_err = 0.0_8
        err = 0.0_8
        do k = 1, num_fw
            if( fa_inbound( k ) == 0 .and. fa_type( k ) < tpplum1 ) then
                ! First step
                aqz = fa_Cfric( k ) * sqrt( fa_xvel( k )**2 + fa_yvel( k )**2 )
                k1x = -aqz * fa_xvel( k )
                k1y = -aqz * fa_yvel( k )
                ! Second step
                vxt = fa_xvel( k ) + dt * ( b21 * k1x + a2 * fa_accx( k ) )
                vyt = fa_yvel( k ) + dt * ( b21 * k1y + a2 * fa_accy( k ) )
                aqz = fa_Cfric( k ) * sqrt( vxt**2 + vyt**2 )
                k2x = -aqz * vxt
                k2y = -aqz * vyt
                ! Third step
                vxt = fa_xvel( k ) + dt * ( b31 * k1x + b32 * k2x + a3 * fa_accx( k ) )
                vyt = fa_yvel( k ) + dt * ( b31 * k1y + b32 * k2y + a3 * fa_accy( k ) )
                aqz = fa_Cfric( k ) * sqrt( vxt**2 + vyt**2 )
                k3x = -aqz * vxt
                k3y = -aqz * vyt
                ! Fourth step
                vxt = fa_xvel( k ) + dt * ( b41 * k1x + b42 * k2x + b43 * k3x + a4 * fa_accx( k ) )
                vyt = fa_yvel( k ) + dt * ( b41 * k1y + b42 * k2y + b43 * k3y + a4 * fa_accy( k ) )
                aqz = fa_Cfric( k ) * sqrt( vxt**2 + vyt**2 )
                k4x = -aqz * vxt
                k4y = -aqz * vyt
                ! Fifth step
                vxt = fa_xvel( k ) + dt * ( b51 * k1x + b52 * k2x + b53 * k3x + b54 * k4x + a5 * fa_accx( k ) )
                vyt = fa_yvel( k ) + dt * ( b51 * k1y + b52 * k2y + b53 * k3y + b54 * k4y + a5 * fa_accy( k ) )
                aqz = fa_Cfric( k ) * sqrt( vxt**2 + vyt**2 )
                k5x = -aqz * vxt
                k5y = -aqz * vyt
                ! Sixth step
                vxt = fa_xvel( k ) + dt * ( b61 * k1x + b62 * k2x + b63 * k3x + &
                    b64 * k4x + b65 * k5x + a6 * fa_accx( k ) )
                vyt = fa_yvel( k ) + dt * ( b61 * k1y + b62 * k2y + b63 * k3y + &
                    b64 * k4y + b65 * k5y + a6 * fa_accy( k ) )
                aqz = fa_Cfric( k ) * sqrt( vxt**2 + vyt**2 )
                k6x = -aqz * vxt
                k6y = -aqz * vyt
                ! Calculate sixth and 4th order solution
                frx6 = c1 * k1x + c3 * k3x + c4 * k4x + c6 * k6x
                fry6 = c1 * k1y + c3 * k3y + c4 * k4y + c6 * k6y
                frx4 = d1 * k1x + d3 * k3x + d4 * k4x + d5 * k5x + d6 * k6x
                fry4 = d1 * k1y + d3 * k3y + d4 * k4y + d5 * k5y + d6 * k6y
                ! Calculate a second 4th order solution- standard runge kutta second step
                vxt = fa_xvel( k ) + 0.5_8 * dt * ( k1x + fa_accx( k ) )
                vyt = fa_yvel( k ) + 0.5_8 * dt * ( k1y + fa_accy( k ) )
                aqz = fa_Cfric( k ) * sqrt( vxt**2 + vyt**2 )
                k2x = -aqz * vxt
                k2y = -aqz * vyt
                ! Third step
                vxt = fa_xvel( k ) + 0.5_8 * dt * ( k2x + fa_accx( k ) )
                vyt = fa_yvel( k ) + 0.5_8 * dt * ( k2y + fa_accy( k ) )
                aqz = fa_Cfric( k ) * sqrt( vxt**2 + vyt**2 )
                k3x = -aqz * vxt
                k3y = -aqz * vyt
                ! Fourth step
                vxt = fa_xvel( k ) + dt * ( k3x + fa_accx( k ) )
                vyt = fa_yvel( k ) + dt * ( k3y + fa_accy( k ) )
                aqz = fa_Cfric( k ) * sqrt( vxt**2 + vyt**2 )
                k4x = -aqz * vxt
                k4y = -aqz * vyt
                frx4a = ( k1x + 2.0_8 * ( k2x + k3x ) + k4x ) / 6._8
                fry4a = ( k1y + 2.0_8 * ( k2y + k3y ) + k4y ) / 6._8
                ! Use the 6th order solution as the final guess.
                fa_xvelO( k ) = dble( fa_xvel( k ) + dt * ( frx6 + fa_accx( k ) ) )
                fa_yvelO( k ) = dble( fa_yvel( k ) + dt * ( fry6 + fa_accy( k ) ) )
                xva = abs( fa_xvel( k ) )
                yva = abs( fa_yvel( k ) )
                xvb = abs( fa_xvelO( k ) )
                yvb = abs( fa_yvelO( k ) )
                xva = min( xva, xvb )
                yva = min( yva, yvb )
                xva = max( xva, 0.001_8 )
                yva = max( yva, 0.001_8 )
                ! Check the error size of the solution
                xer = max( abs( frx6 - frx4 ), abs( frx6 - frx4a ) )
                xer = max( xer, abs( frx4 - frx4a ) )
                xer = accuracy * dt * xer / min( xva, 0.7_8 )
                yer = max( abs( fry6 - fry4 ), abs( fry6 - fry4a ) )
                yer = max( yer, abs( fry4 - fry4a ) )
                yer = accuracy * dt * yer / min( yva, 0.7_8 )
                xer = max( xer, yer )
                err = max( err, xer )
                ! Don't step more than a cell in one step (to ensure
                ! even distribution of sediment)
                yer = 2._8 * dt * max( xvb, yvb ) / (strat_dx)
                err = max( err, yer )
            endif

        enddo

        call mpi_allreduce( err, all_err, 1, dbl_type, max_type, SPModel_comm_world,ierr)

        return

    end subroutine Cash_Karp_RKmethod
  ! ============================================================================

end module mod_rungekutta
