! ============================================================================
! Name        : OceanIO.f90
! Author      : Tristan Salles
! Created on: March 24, 2014
! Copyright (C) 2014 CSIRO
! ============================================================================
!> \file OceanIO.f90
!
! Description : OceanIO.f90 is making the connextion between LECODE and the Ocean
! circulation model.
!
!<
! ============================================================================
module ocean_io

    use file_data
    use mpidata
    use ocean_data
    use FoX_sax
    use FoX_common
    use flux_data

    implicit none


    integer :: hclass, sclass
    integer :: circflag, waveflag

    ! Wave Parameters Region
    logical, save :: in_ocean = .false.
    logical, save :: in_circForcing = .false.
    logical, save :: in_waveForcing = .false.
    logical, save :: in_morphLimit = .false.
    logical, save :: in_morphfac = .false.
    logical, save :: in_osyncfolder = .false.
    logical, save :: in_osynctime = .false.
    logical, save :: in_hindcast = .false.
    logical, save :: in_wwcast = .false.
    logical, save :: in_wtstart = .false.
    logical, save :: in_wtend = .false.
    logical, save :: in_groupnb = .false.
    logical, save :: in_groupperc = .false.
    logical, save :: in_maxdd = .false.
    character(len=128), save :: morphLimit, grid_waveForcing, grid_circForcing, osyncfolder, maxdd
    character(len=128), save :: osynctime, wwcast, wtstart, wtend, groupnb, groupperc, morphfac

contains


    ! ============================================================================
    subroutine startDocument_handler

    end subroutine startDocument_handler
    ! ============================================================================
    subroutine startElement_handler(namespaceURI, localname, name, atts)

        character(len=*), intent(in) :: namespaceURI
        character(len=*), intent(in) :: localname
        character(len=*), intent(in) :: name
        type(dictionary_t), intent(in) :: atts

        ! Ocean Element
        if (name=='OceanPlugin')then
            in_ocean = .true.
            ocean_plug = .true.
        endif
        if(in_ocean) call SoceanElement_handler(name)
        if (name=='forecastClass') in_hindcast = .true.
        if(in_hindcast) call ShindcastElement_handler(name)

    end subroutine startElement_handler
    !============================================================================
    subroutine endElement_handler(namespaceURI, localname, name)

        character(len=*), intent(in) :: namespaceURI
        character(len=*), intent(in) :: localname
        character(len=*), intent(in) :: name

        ! Wave Element
        call EoceanElement_handler(name)
        call EhindcastElement_handler(name)

    end subroutine endElement_handler
    ! ============================================================================
    subroutine characters_handler(chars)

        character(len=*), intent(in) :: chars

        ! Get Ocean Parameters
        if(in_ocean)   call ocean_characters_handler(chars)
        if(in_hindcast)   call hindcast_characters_handler(chars)

    end subroutine characters_handler
    ! ============================================================================
    subroutine SoceanElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='circForcing') in_circForcing = .true.
        if (name=='waveForcing') in_waveForcing = .true.
        if (name=='actionMaxDepth') in_maxdd = .true.
        if (name=='OsyncFolder') in_osyncfolder = .true.
        if (name=='morphLimit') in_morphLimit = .true.
        if (name=='morphAc') in_morphfac = .true.
        if (name=='OsyncTime') in_osynctime = .true.
        if (name=='nbForecast') in_wwcast = .true.
        if (name=='forecastClass') in_hindcast = .true.

    end subroutine SoceanElement_handler
    ! ============================================================================
    subroutine ShindcastElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='Ostart') in_wtstart = .true.
        if (name=='Oend') in_wtend = .true.
        if (name=='subgroupNb') in_groupnb = .true.
        if (name=='subgroupProp') in_groupperc = .true.
        if( in_groupnb )then
            hclass = hclass + 1
            sclass = 0
        endif

    end subroutine ShindcastElement_handler
    ! ============================================================================
    subroutine EoceanElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='circForcing') in_circForcing = .false.
        if (name=='waveForcing') in_waveForcing = .false.
        if (name=='actionMaxDepth') in_maxdd = .false.
        if (name=='OsyncFolder') in_osyncfolder = .false.
        if (name=='OsyncTime') in_osynctime = .false.
        if (name=='morphLimit') in_morphLimit = .false.
        if (name=='morphAc') in_morphfac = .false.
        if (name=='nbForecast') in_wwcast = .false.
        if (name=='forecastClass') in_hindcast = .false.

    end subroutine EoceanElement_handler
    ! ============================================================================
    subroutine EhindcastElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='Ostart') in_wtstart = .false.
        if (name=='Oend') in_wtend = .false.
        if (name=='subgroupNb') in_groupnb = .false.
        if (name=='subgroupProp') in_groupperc = .false.

    end subroutine EhindcastElement_handler
    ! ============================================================================
    subroutine ocean_characters_handler(chars)

        character(len=*), intent(in) :: chars

        if( in_maxdd ) then
            maxdd = chars
            call rts(maxdd,max_WC_action_depth)
        elseif (in_circForcing) then
            grid_circForcing = chars
            call rts(grid_circForcing,circflag)
        elseif (in_waveForcing) then
            grid_waveForcing = chars
            call rts(grid_waveForcing,waveflag)
        elseif ( in_wwcast )then
            wwcast = chars
            call rts(wwcast, forecast_nb )
            allocate( hindcast( forecast_nb ) )
        elseif (in_osyncfolder) then
            outdir4 = chars
        elseif (in_morphLimit) then
            morphLimit = chars
            call rts(morphLimit,ocean_morph_th)
        elseif (in_morphfac) then
            morphfac = chars
            call rts(morphfac,transport%morpho)
        elseif (in_osynctime) then
            osynctime = chars
            call rts(osynctime,ocean_time_step)
        endif

    end subroutine ocean_characters_handler
    ! ============================================================================
    subroutine hindcast_characters_handler(chars)

        character(len=*), intent(in) :: chars
        real :: sbperc

        if (in_wtstart) then
            wtstart = chars
            call rts(wtstart, hindcast( hclass )%tstart )
        elseif (in_wtend) then
            wtend = chars
            call rts(wtend, hindcast( hclass )%tend )
        elseif (in_groupnb) then
            groupnb = chars
            call rts(groupnb, hindcast( hclass )%cnb )
            allocate( hindcast( hclass )%subgroup( hindcast( hclass )%cnb ) )
        elseif( in_groupperc )then
            sclass = sclass + 1
            groupperc = chars
            call rts(groupperc,sbperc)
            hindcast( hclass )%subgroup( sclass )%perc = sbperc
        endif

    end subroutine hindcast_characters_handler
    ! ============================================================================
    subroutine endDocument_handler

    end subroutine endDocument_handler
    ! ============================================================================
    !> Subroutine xml_Ocean()
    !! reads input that defines the Ocean element
    !<
    ! ============================================================================
    subroutine xml_Ocean

        logical :: found

        type(xml_t) :: xf

        integer :: ios, l1

        character(len=128) :: command, stg, fildir

        hclass = 0
        sclass = 0
        circflag = 0
        waveflag = 0
        ocean_plug = .false.
        current_on = .false.
        wave_on = .false.
        max_WC_action_depth = 200.0_8

        ! Open file
        call open_xml_file(xf, finput , ios)
        if (ios/=0) then
            print*,'---------------------'
            print*, 'Error opening input file for parsing Ocean XmL'
            print*,'---------------------'
        endif

        ! Parse file at first and find data
        call parse(xf, &
            startDocument_handler = startDocument_handler, &
            startElement_handler = startElement_handler, &
            endElement_handler = endElement_handler, &
            characters_handler = characters_handler, &
            endDocument_handler = endDocument_handler &
            )

        ! Close file
        call close_xml_t(xf)

        ! Define ocean synchronisation folder and files
        if( ocean_plug )then

            ! Ocean Lecode shared folder
            call noblnk(outdir4)
            fildir = outdir4
            stg = '/topSPM_surface.xyz'
            call noblnk(stg)
            call append_str( fildir,stg )
            call noblnk(fildir)
            inquire(file=fildir,exist=found)
            if( .not. found )then
                command=' '
                command(1:6)='mkdir '
                l1 = len_trim(outdir4)
                command(7:l1+7)=outdir4
                call term_command( command )
            else
                command=' '
                command(1:6) = 'rm -r '
                l1 = len_trim(outdir4)
                command(7:l1+7)=outdir4
                call term_command( command )
                command=' '
                command(1:6)='mkdir '
                l1 = len_trim(outdir4)
                command(7:l1+7)=outdir4
                call term_command( command )
            endif
            poseidon = 'poseidon'
            call noblnk( poseidon )
            call addpath4( poseidon )

            foceanxyz = 'topSPM_surface.xyz'
            call noblnk( foceanxyz )
            call addpath4( foceanxyz )

            if( circflag == 1 ) current_on = .true.
            if( waveflag == 1 ) wave_on = .true.

        endif

        return

    end subroutine xml_Ocean
    ! ============================================================================
    !> Subroutine broadcast_flux()
    !! Broadcast initial parameters which define the flux.
    !<
    ! ============================================================================
    subroutine broadcast_ocean

        integer :: k, sg

        call mpi_bcast( ocean_plug,1,lgc_type,0,SPModel_comm_world,ierr )
        call mpi_bcast( current_on,1,lgc_type,0,SPModel_comm_world,ierr )
        call mpi_bcast( wave_on,1,lgc_type,0,SPModel_comm_world,ierr )
        call mpi_bcast( max_WC_action_depth,1,dbl_type,0,SPModel_comm_world,ierr )

        if( ocean_plug )then

            call mpi_bcast( outdir4,128,mpi_character,0,SPModel_comm_world,ierr )
            call mpi_bcast( poseidon,128,mpi_character,0,SPModel_comm_world,ierr )
            call mpi_bcast( foceanxyz,128,mpi_character,0,SPModel_comm_world,ierr )

            call mpi_bcast( forecast_nb,1,int_type,0,SPModel_comm_world,ierr )
            call mpi_bcast( ocean_time_step,1,dbl_type,0,SPModel_comm_world,ierr )
            call mpi_bcast( ocean_morph_th,1,dbl_type,0,SPModel_comm_world,ierr )
            call mpi_bcast( transport%morpho,1,dbl_type,0,SPModel_comm_world,ierr )

            if( transport%morpho <= 0.0_8 ) transport%morpho = 1.0_8

            if( iam /= 0 ) allocate( hindcast( forecast_nb ) )

            do k = 1, forecast_nb

                call mpi_bcast( hindcast( k )%cnb,1,int_type,0,SPModel_comm_world,ierr )
                call mpi_bcast( hindcast( k )%tstart,1,dbl_type,0,SPModel_comm_world,ierr )
                call mpi_bcast( hindcast( k )%tend,1,dbl_type,0,SPModel_comm_world,ierr )
                if( iam /= 0 ) allocate( hindcast( k )%subgroup( hindcast( k )%cnb ) )

                do sg = 1, hindcast( k )%cnb
                    call mpi_bcast( hindcast( k )%subgroup( sg )%perc,1,dbl_type,0,SPModel_comm_world,ierr )
                enddo

            enddo

        endif

        return

    end subroutine broadcast_ocean
    ! ============================================================================


end module ocean_io
