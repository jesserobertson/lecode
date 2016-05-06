! ============================================================================
! Name        : License.f90
! Author      : tristan salles
! Created on: June 07, 2013
! Copyright (C) 2013 CSIRO
!============================================================================
!> \file License.f90
!!
!! File License setup the license and time of validity of the software for a specific date and mac address.
!! This is quite elementary and easy to crack but I have more important things to do then trying to protect
!! the code.
!<
! ============================================================================
module license

    use precision_data

    implicit none

    logical, parameter:: datelimit=.false.
    logical, parameter:: machineid=.false.

    integer, parameter:: numlicen=2, keylen=17

    ! Month and year license validity
    integer, parameter :: endyear = 2012
    integer, parameter :: endmonth = 6

    ! Virtual Box user mac address
    character(len=keylen), parameter:: macad( numlicen )=(/'00:21:9B:67:27:33','00-21-9B-67-27-33'/)

contains

    ! ============================================================================
    !> Subroutine init_license.
    !!
    !<
    ! ============================================================================
    subroutine init_license

        logical:: licmach

        integer :: ios, iu, k, j
        integer, dimension( 8 ) :: values

        character( len=10 ) :: date1,time1,zone1
        character( len=128 ) :: text, cmd, tmpfil

        licmach = .false.

        if( datelimit )then
            call date_and_time(date1,time1,zone1,values)
            if( values(1) > endyear )then
                print*,'LICENCE EXPIRED'
                stop
            else
                if( values(2) > endmonth )then
                    print*,'LICENCE EXPIRED'
                    stop
                endif
            endif
        endif

        if( machineid )then
            iu = 10
            cmd = '/sbin/ifconfig -a | grep eth0 > system._tmp_'
            call term_command( cmd )
            tmpfil = 'system._tmp_'
            open( iu, file=tmpfil, status="old", action="read", iostat = ios )
            do
                text=''
                read( iu,'(a128)', iostat=ios, end=30) text
                do k = 1, numlicen
                    do j = 1, 129 - keylen
                        if( text( j:j + keylen - 1 ) == macad( k )(1 : keylen ) ) licmach = .true.
                    enddo
                enddo
            enddo
30      continue

            close( iu )

            ! Delete the file
            cmd='rm -f system._tmp_'
            call term_command( cmd )

            if( .not. licmach )then
                print*,'Valid licence not found'
                stop
            endif

        endif

        return

    end subroutine init_license
    ! ============================================================================

end module license
