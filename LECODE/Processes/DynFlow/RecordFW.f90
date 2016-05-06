! ============================================================================
! Name        : RecordFW.f90
! Author      : tristan salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file RecordFW.f90
!!
!! RecordFW records flow walkers position and velocity for output.
!!
!! The file contains the following subroutine:
!! \arg \c record_flow_walkers
!<
! ============================================================================
module mod_recordfw

   use mpidata
   use forces_data
   use fwalker_data
   use time_data

   implicit none

   public

contains

   ! ============================================================================
   !> Subroutine record_flow_walkers
   !! RecordFW records flow walkers position and velocity.
   !<
   ! ============================================================================
   subroutine record_flow_walkers( first )

      integer :: first, k, k1, k2, m, kk, p
      integer,dimension( : ), allocatable:: lasttype
      real( tkind ),dimension( : ), allocatable :: lastxpos, lastypos, lastzpos, lastxvel, &
        lastsedb, lastyvel, lastseds, lasth, lastw

      allocate( lasttype( num_fw_o ), lastxpos( num_fw_o ), lastypos( num_fw_o ), lastzpos( num_fw_o ), &
        lastsedb( num_fw_o ) )
      allocate( lastxvel( num_fw_o ), lastyvel( num_fw_o ), lastseds( num_fw_o ), lasth( num_fw_o ), lastw( num_fw_o ) )

      ! Record flow walkers for output
      k = num_fw_o + num_fw
      k1 = k
      if( first == 0 )then
        allocate( recfw_type( maxfwinfo ) )
        allocate( recfw_xpos( maxfwinfo ) )
        allocate( recfw_ypos( maxfwinfo ) )
        allocate( recfw_zpos( maxfwinfo ) )
        allocate( recfw_xvel( maxfwinfo ) )
        allocate( recfw_yvel( maxfwinfo ) )
        allocate( recfw_srch( maxfwinfo ) )
        allocate( recfw_width( maxfwinfo ) )
        allocate( recfw_volume( maxfwinfo ) )
        allocate( recfw_sedcharge( maxfwinfo ) )
      endif
      if( k < maxfwinfo .and. k > 0 )then
         if( first == 0 )then
            k2 = 1
            do m = 1, k
               if( fa_xvel( m )**2 + fa_yvel( m )**2 > transport%fvmin**2 )then
                  recfw_xpos( k2 ) = real( fa_xpos( m ) )
                  recfw_ypos( k2 ) = real( fa_ypos( m ) )
                  recfw_zpos( k2 ) = real( fa_zpos( m ) )
                  recfw_srch( k2 ) = real( fa_h( m ) )
                  recfw_width( k2 ) = real( fa_width( m ) )
                  recfw_xvel( k2 ) = dble( fa_xvel( m ) )
                  recfw_yvel( k2 ) = dble( fa_yvel( m ) )
                  recfw_type( k2 ) = 0
                  if( fa_type( m ) >= tprain .and. fa_type( m ) < tprain + combine )&
                  recfw_type( k2 ) = 1
                  if( fa_type( m ) >= tpturb )&
                  recfw_type( k2 ) = 2
                  if( fa_type( m ) >= tpplum1 )then
                      recfw_type( k2 ) = 3
                  endif
                  recfw_sedcharge( k2 ) = 0.0_8
                  recfw_volume( k2 ) = fa_volume( m )
                  do p = 1, totgrn
                     recfw_sedcharge( k2 ) = recfw_sedcharge( k2 ) + &
                        fa_sedcharge( m , p ) * sediment( p )%density
                  enddo
                  if( fa_type( m ) >= tprain .and. fa_type( m ) < tpturb )then
                      recfw_sedcharge( k2 ) = recfw_sedcharge( k2 ) / rain_int
                  else
                      recfw_sedcharge( k2 ) = recfw_sedcharge( k2 ) / flow_int
                  endif
                  k2 = k2 + 1
               endif
            enddo
            k1 = k2 - 1
         else
            do m = 1, num_fw_o
               lastxpos( m ) = recfw_xpos( m )
               lastypos( m ) = recfw_ypos( m )
               lastzpos( m ) = recfw_zpos( m )
               lastxvel( m ) = recfw_xvel( m )
               lasth( m ) = recfw_srch( m )
               lastw( m ) = recfw_width( m )
               lastyvel( m ) = recfw_yvel( m )
               lasttype( m ) = recfw_type( m )
               lastseds( m ) = recfw_volume( m )
               lastsedb( m ) = recfw_sedcharge( m )
            enddo
            k2 = 1
            do m = 1, num_fw
               if( fa_xvel( m )**2 + fa_yvel( m )**2 > transport%fvmin**2 )then
                  recfw_xpos( k2 ) = real( fa_xpos( m ) )
                  recfw_ypos( k2 ) = real( fa_ypos( m ) )
                  recfw_zpos( k2 ) = real( fa_zpos( m ) )
                  recfw_srch( k2 ) = real( fa_h( m ) )
                  recfw_width( k2 ) = real( fa_width( m ) )
                  recfw_xvel( k2 ) = dble( fa_xvel( m ) )
                  recfw_yvel( k2 ) = dble( fa_yvel( m ) )
                  recfw_type( k2 ) = 0
                  if( fa_type( m ) >= tprain .and. fa_type( m ) < tprain + combine )&
                  recfw_type( k2 ) = 1
                  if( fa_zpos( m ) < gsea%actual_sea .and. fa_density( m ) >= sea_density )&
                  recfw_type( k2 ) = 2
                  if( fa_type( m ) >= tpturb )&
                  recfw_type( k2 ) = 2
                  if( fa_type( m ) >= tpplum1 )then
                      recfw_type( k2 ) = 3
                  endif
                  recfw_volume( k2 ) = dble( fa_volume( m ) )
                  recfw_sedcharge( k2 ) = 0.0_8
                  do p = 1, totgrn
                     recfw_sedcharge( k2 ) = recfw_sedcharge( k2 ) + &
                        fa_sedcharge( m , p ) * sediment( p )%density
                  enddo
                  if( fa_type( m ) >= tprain .and. fa_type( m ) < tpturb )then
                      recfw_sedcharge( k2 ) = recfw_sedcharge( k2 ) / rain_int
                  else
                      recfw_sedcharge( k2 ) = recfw_sedcharge( k2 ) / flow_int
                  endif
                  k2 = k2 + 1
               endif
            enddo
            k1 = k2 + num_fw_o -1
            do m = k2, k1
               recfw_xpos( m ) = lastxpos( m - k2  + 1 )
               recfw_ypos( m ) = lastypos( m - k2 + 1 )
               recfw_zpos( m ) = lastzpos( m - k2 + 1 )
               recfw_xvel( m ) = lastxvel( m - k2 + 1 )
               recfw_yvel( m ) = lastyvel( m - k2 + 1 )
               recfw_srch( m ) = lasth( m - k2 + 1 )
               recfw_width( m ) = lastw( m - k2 + 1 )
               recfw_type( m ) = lasttype( m - k2 + 1 )
               recfw_volume( m ) = lastseds( m - k2 + 1 )
               recfw_sedcharge( m ) = lastsedb( m - k2 + 1 )
            enddo
         endif
      elseif( k >= maxfwinfo )then
         if( first == 0 )then
            k2 = 1
            do m = 1, maxfwinfo
               if( fa_xvel( m )**2 + fa_yvel( m )**2 > transport%fvmin**2 )then
                  recfw_xpos( k2 ) = real( fa_xpos( m ) )
                  recfw_ypos( k2 ) = real( fa_ypos( m ) )
                  recfw_zpos( k2 ) = real( fa_zpos( m ) )
                  recfw_xvel( k2 ) = dble( fa_xvel( m ) )
                  recfw_yvel( k2 ) = dble( fa_yvel( m ) )
                  recfw_srch( k2 ) = real( fa_h( m ) )
                  recfw_width( k2 ) = real( fa_width( m ) )
                  recfw_volume( k2 ) = dble( fa_volume( m ) )
                  recfw_sedcharge( k2 ) = 0.0_8
                  do p = 1, totgrn
                     recfw_sedcharge( k2 ) = recfw_sedcharge( k2 ) + &
                        fa_sedcharge( m , p ) * sediment( p )%density
                  enddo
                  if( fa_type( m ) >= tprain .and. fa_type( m ) < tpturb )then
                      recfw_sedcharge( k2 ) = recfw_sedcharge( k2 ) / rain_int
                  else
                      recfw_sedcharge( k2 ) = recfw_sedcharge( k2 ) / flow_int
                  endif
                  recfw_type( k2 ) = 0
                  if( fa_type( m ) >= tprain .and. fa_type( m ) < tprain + combine )&
                  recfw_type( k2 ) = 1
                  if( fa_zpos( m ) < gsea%actual_sea .and. fa_density( m ) >= sea_density )&
                  recfw_type( k2 ) = 2
                  if( fa_type( m ) >= tpturb )&
                  recfw_type( k2 ) = 2
                  if( fa_type( m ) >= tpplum1 )then
                      recfw_type( k2 ) = 3
                  endif
                  k2 = k2 + 1
               endif
            enddo
            k1 = k2 - 1
         else
            do m = 1, num_fw_o
               lastxpos( m ) = recfw_xpos( m )
               lastypos( m ) = recfw_ypos( m )
               lastzpos( m ) = recfw_zpos( m )
               lastxvel( m ) = recfw_xvel( m )
               lastyvel( m ) = recfw_yvel( m )
               lasth( m ) = recfw_srch( m )
               lastw( m ) = recfw_width( m )
               lasttype( m ) = recfw_type( m )
               lastseds( m ) = recfw_volume( m )
               lastsedb( m ) = recfw_sedcharge( m )
            enddo
            if( num_fw < maxfwinfo )then
               k2 = 1
               do m = 1, num_fw
                  if( fa_xvel( m )**2 + fa_yvel( m )**2 > transport%fvmin**2 )then
                     recfw_xpos( k2 ) = real( fa_xpos( m ) )
                     recfw_ypos( k2 ) = real( fa_ypos( m ) )
                     recfw_zpos( k2 ) = real( fa_zpos( m ) )
                     recfw_xvel( k2 ) = dble( fa_xvel( m ) )
                     recfw_yvel( k2 ) = dble( fa_yvel( m ) )
                     recfw_srch( k2 ) = real( fa_h( m ) )
                     recfw_width( k2 ) = real( fa_width( m ) )
                     recfw_sedcharge( k2 ) = 0.0_8
                     recfw_volume( k2 ) = real( fa_volume( m ) )
                     do p = 1, totgrn
                         recfw_sedcharge( k2 ) = recfw_sedcharge( k2 ) + &
                            fa_sedcharge( m , p ) * sediment( p )%density
                     enddo
                      if( fa_type( m ) >= tprain .and. fa_type( m ) < tpturb )then
                          recfw_sedcharge( k2 ) = recfw_sedcharge( k2 ) / rain_int
                      else
                          recfw_sedcharge( k2 ) = recfw_sedcharge( k2 ) / flow_int
                      endif
                     recfw_type( k2 ) = 0
                     if( fa_type( m ) >= tprain .and. fa_type( m ) < tprain + combine )&
                     recfw_type( k2 ) = 1
                     if( fa_zpos( m ) < gsea%actual_sea .and. fa_density( m ) >= sea_density )&
                     recfw_type( k2 ) = 2
                     if( fa_type( m ) >= tpturb )&
                     recfw_type( k2 ) = 2
                     if( fa_type( m ) >= tpplum1 )then
                        recfw_type( k2 ) = 3
                     endif
                     k2 = k2 +1
                  endif
               enddo
               k1 = k2 - 1
               kk = 0
               do m = num_fw + 1, maxfwinfo
                  kk = kk + 1
                  recfw_xpos( m ) = lastxpos( kk )
                  recfw_ypos( m ) = lastypos( kk )
                  recfw_zpos( m ) = lastzpos( kk )
                  recfw_xvel( m ) = lastxvel( kk )
                  recfw_yvel( m ) = lastyvel( kk )
                  recfw_srch( m ) = lasth( kk )
                  recfw_width( m ) = lastw( kk )
                  recfw_type( m ) = lasttype( kk )
                  recfw_volume( m ) = lastseds( kk )
                  recfw_sedcharge( m ) = lastsedb( kk )
                  k1 = k1 + 1
               enddo
               k1 = k1 - 1
            else
               k2 = 1
               do m = 1, maxfwinfo
                  if( fa_xvel( m )**2 + fa_yvel( m )**2 < transport%fvmin**2 )then
                     recfw_xpos( k2 ) = real( fa_xpos( m ) )
                     recfw_ypos( k2 ) = real( fa_ypos( m ) )
                     recfw_zpos( k2 ) = real( fa_zpos( m ) )
                     recfw_xvel( k2 ) = dble( fa_xvel( m ) )
                     recfw_srch( k2 ) = real( fa_h( m ) )
                     recfw_width( k2 ) = real( fa_width( m ) )
                     recfw_yvel( k2 ) = dble( fa_yvel( m ) )
                     recfw_volume( k2 ) = dble( fa_volume( m ) )
                     recfw_sedcharge( k2 ) = 0.0_8
                     do p = 1, totgrn
                         recfw_sedcharge( k2 ) = recfw_sedcharge( k2 ) + &
                            fa_sedcharge( m , p ) * sediment( p )%density
                     enddo
                      if( fa_type( m ) >= tprain .and. fa_type( m ) < tpturb )then
                          recfw_sedcharge( k2 ) = recfw_sedcharge( k2 ) / rain_int
                      else
                          recfw_sedcharge( k2 ) = recfw_sedcharge( k2 ) / flow_int
                      endif
                     recfw_type( k2 ) = 0
                     if( fa_type( m ) >= tprain .and. fa_type( m ) < tprain + combine )&
                     recfw_type( k2 ) = 1
                     if( fa_zpos( m ) < gsea%actual_sea .and. fa_density( m ) >= sea_density )&
                     recfw_type( k2 ) = 2
                     if( fa_type( m ) >= tpturb )&
                     recfw_type( k2 ) = 2
                     if( fa_type( m ) >= tpplum1 )then
                        recfw_type( k2 ) = 3
                     endif
                     k2 = k2 + 1
                  endif
               enddo
               k1 = k2 - 1
            endif
         endif
      endif
      num_fw_o = k1

      deallocate( lasttype, lastxpos, lastypos, lastzpos, lastxvel, lastyvel, lastseds, lastsedb, lasth, lastw )

      return

   end subroutine record_flow_walkers
  ! ============================================================================

end module mod_recordfw
