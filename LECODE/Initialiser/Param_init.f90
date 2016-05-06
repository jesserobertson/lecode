! ============================================================================
! Name        : Param_init.f90
! Author      : Tristan Salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file Param_init.f90
!!
!! Param_init is called at the beginning of the simulation to
!! initialize some parameters.
!!
!<
! ============================================================================
module params_init

   use mpidata
   use time_data
   use fwalker_data
   use rain_initial
   use strata_data
   use forces_data
   use flux_data
   use mod_sealevel
   use mod_tempsal
   use mod_displacements

   implicit none

   real( tkind ), parameter :: fvc1=18.0_8
   real( tkind ), parameter :: fvc2=1.0_8

contains

   ! ============================================================================
   !> Subroutine spm_parameters_initialise
   !! Subroutine spm_parameters_initialise initialises simulation paramters.
   ! ============================================================================
   subroutine spm_parameters_initialise

      integer :: i, ks

      real( tkind ) :: dens, topfv, fv, flowtime

      dtnext = 5.e2_8
      gdisp%lastdisp = time_start

      ! Define time tolerance accuracy
      i = 0
      if( time_start /= 0.0_8 ) i = int( log10( abs( time_start ) ) )
      if( time_end /= 0.0_8 ) i = max( i , int( log10( abs( time_end ) ) ) )
      i = i - 13
      time_tolerance = max( tor, 10.0_8**i )

      ! Time flags initialisation
      next_displacement = 0.0
      next_display = 0.0
      next_output = 0.0

      ! Source section
      topfv = 0.0
      tot_fw = 0
      num_fw = 0
      num_fw_o = 0
      flowtime = flow_int

      if( num_src > 0 )then
          do i = 1, num_src

             ! Calculate the Euclidean norm of the FW's current velocity.
             fv = sqrt( fws_xvel( i )**2 + fws_yvel( i )**2 )
             topfv = max( topfv, fv )

             ! Determine number of flow walkers released at source
             ! by checking against max. flow walker's depth
             if( fws_srch( i ) <= max_elev_src ) then
                fws_num( i ) = 1
             else
                fws_num( i ) = int( fws_srch( i ) / max_elev_src )
                if( mod( fws_srch( i ), max_elev_src ) > 0.0_8 ) fws_num( i ) = fws_num( i ) + 1
                fws_srch( i ) = fws_srch( i ) / real( fws_num( i ), 8 )
                fws_qfl( i ) = fws_qfl( i ) / real( fws_num( i ), 8 )
                fws_sedconc( i ) = fws_sedconc( i ) / real( fws_num( i ), 8 )
             endif

             ! Transform sed_conc from Mt/year to cubic metres per second
             dens = 0.0_8
             do ks = 1, silgrn
                dens = dens + fws_sedperc( i , ks ) * sediment( ks )%density * 0.01_8
             enddo
             fws_sedconc( i ) = fws_sedconc( i ) * 1.e9_8 / ( dens * secyear )

             ! Sedcharge is the volume of each sediment present in the source cubic metres per second
             do ks = 1, silgrn
                fws_sedcharge( i, ks ) = fws_sedperc( i, ks ) * 0.01_8 * fws_sedconc( i )
             enddo

             if(fws_srch( i ) > max_elev_src ) print*,'WARNING: preset flow height higher than max allowed'

             flowtime = min( flowtime, 1.e6_8 )
          enddo
      endif

      ! Update flow interval
      if( topfv > 0.0_8 ) dtnext = 0.5_8 * strat_dx / topfv
      if( time_display < flowtime ) flowtime = time_display
      flow_int = flowtime
      do  i = 1, num_src
        if( fws_tend( i ) - fws_tstrt( i ) < flow_int ) then
            flowtime = fws_tend( i ) - fws_tstrt( i )
        else
            flowtime = flow_int
        endif

        ! Define source duration in second based of percentage of recurrence and flow interval
        flowtime = flowtime * secyear * fws_perc( i ) * 0.01_8

        ! According to flow duration compute for each sediment class the volume of materials
        ! initially transported ( in cubic metres )
        do ks = 1, silgrn
            fws_sedcharge( i, ks ) = flowtime * fws_sedcharge( i, ks )
        enddo
         fws_volume( i ) = flowtime * fws_qfl( i )
      enddo

      ! Convert rain flow rate in m3/s for the rain interval
      if( geomorphology ) rainstreamQ = rainstreamQ * rain_int * secyear

      ! Sea level section
      gsea%actual_sea = 0.0_8
      if( gsea%sealevel ) call read_sealevel_file

      ! Temperature and salinity section
      if( gtempsal%tpsl ) call read_temperature_salinity_file

      ! Initialise coarse dt to sampling interval time
      coarse_dt = sampling_int

      ! Make sure display interval is a product of sampling interval
      if( int( time_display / coarse_dt ) * coarse_dt /= time_display )&
        time_display = ( int( time_display / coarse_dt ) + 1 ) * coarse_dt

      ! Make sure output interval is a product of sampling interval
      if( int( time_output / coarse_dt ) * coarse_dt /= time_output )&
        time_output = ( int( time_output / coarse_dt ) + 1 ) * coarse_dt

      ! Change units assigned to grain diameters from millimeters to meters.
      do ks = 1, totgrn
         sediment( ks )%diameter = sediment( ks )%diameter * 0.001_8
      enddo

      ! Fall velocities computation based on Zhang
      do ks = 1, totgrn
         if( sediment(ks)%diameter /= 0.0_8 )then
            if( sediment(ks)%diameter > 0.062 * 0.001_8 )then
                sediment(ks)%vfall = settling_velocity_zhang(ks)
            else
                sediment(ks)%vfall = 0.5_8 * 0.001_8
            endif
         else
            sediment(ks)%vfall = 0.0_8
         endif
      end do

      ! Read rainfall file and allocate rain arrays
      if( geomorphology ) call read_rain_file

      ! Order in which the different grain sizes will be deposited,
      ! based on the weight of each grain
      allocate( dep_order( totgrn ) )
      call deposition_order

      ! Allocate displacement fields
      gdisp%disptime = 0.0_8
      if( gdisp%event > 0 .and. .not. udw_plug ) then
         new_disp = .true.
         ! Read and assign vertical displacement rate
         call assign_vertical_displacements_rate( 1 )
         ! Find displacements for considered time
         call update_local_displacements
      endif

      return

   end subroutine spm_parameters_initialise
   ! ============================================================================
   !> Subroutine deposition_order
   !! Subroutine deposition_order computes the order of deposition.
   !<
   ! ============================================================================
   subroutine deposition_order

      integer :: kk, ks, ki, i
      integer :: temp( totgrn )

      real( tkind ) :: heavy
      real( tkind ) :: weight( totgrn )

      ! Compute the weight of each grain considering a spherical shape
      do ks = 1, totgrn
         dep_order( ks ) = 0
         weight( ks ) = sediment( ks )%density
         weight( ks ) = weight( ks ) * sediment( ks )%diameter**3
         weight( ks ) = weight( ks ) * 1.0e9_8
         temp( ks ) = ks
      enddo

      ! Order in which the different grain sizes will be deposited
      do ks = 1, totgrn
         heavy = weight( temp( ks ) )
         kk = temp( ks )
         ki = kk
         do i = ks + 1, totgrn
            if( heavy < weight( temp( i ) ) )then
               kk = temp(i)
               ki = i
               heavy = weight( temp( i ) )
            endif
         enddo
         dep_order( ks ) = kk
         temp( ki ) = temp( ks )
      enddo

      return

   end subroutine deposition_order
   ! ============================================================================
   !> Function settling_velocity_zhang()
   !! Function settling_velocity_zhang is used to compute settling velocity for a considered type of grain
   !! Based on Zhang 1993
   !! \param sed
   !<
   ! ============================================================================
   function settling_velocity_zhang( sed ) result( ws )

      integer, intent( in ) :: sed

      real( tkind ) :: Cvisc, delg, ws

      Cvisc = 1.0e-6_8 * 13.95_8 / sediment( sed )%diameter
      delg=( sediment( sed )%density - fluid_density ) * gravity / fluid_density

      ws = sqrt( Cvisc**2 + 1.09_8 * delg * sediment( sed )%diameter )
      ws = ws - Cvisc

   end function settling_velocity_zhang
   ! ============================================================================

end module params_init
