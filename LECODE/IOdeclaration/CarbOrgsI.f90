! ============================================================================
! Name        : CarbsOrgsI.f90
! Author      : Tristan Salles
! Created on: May 1, 2014
! Copyright (C) 2014 CSIRO
! ============================================================================
!> \file CarbsOrgsI.f90
!
! Description : CarbsOrgsI is used to gather the information within the SPModel XmL input file to
! define the carbonates and organics parameters.
!
!<
! ============================================================================
module carborg_parser

    use file_data
    use FoX_sax
    use mpidata
    use time_data
    use strata_data
    use param_data
    use FoX_common
    use sediment_data
    use forces_data

    implicit none

    public

    integer :: mbshp_Nb, fuzz_Nb, fuzzy_Nb, pts_Nb, coords_Nb, mbs_Nb

    logical, save :: carborgparamsection = .false.
    logical, save :: mbsfctsection = .false.
    logical, save :: fuzzyfctsection = .false.
    logical, save :: fuzzyrulesection = .false.
    logical, save :: in_CarbOrgs = .false.
    logical, save :: in_carborgTime = .false.
    logical, save :: in_tempsalFile = .false.
    logical, save :: in_mbsfct = .false.
    logical, save :: in_fctSubset = .false.
    logical, save :: in_fuzzyfct = .false.
    logical, save :: in_fuzzySubset = .false.
    logical, save :: in_fuzzyrule = .false.
    logical, save :: in_mbsSet = .false.
    logical, save :: in_mbsfctNb = .false.
    logical, save :: in_fuzzyfctNb = .false.
    logical, save :: in_fuzzyruleNb = .false.
    logical, save :: in_fctName = .false.
    logical, save :: in_fctVariable = .false.
    logical, save :: in_fctPts = .false.
    logical, save :: in_mbsfctPoints = .false.
    logical, save :: in_fuzzyName = .false.
    logical, save :: in_fuzzyVariable = .false.
    logical, save :: in_fuzzyPts = .false.
    logical, save :: in_fuzzyCoords = .false.
    logical, save :: in_matName = .false.
    logical, save :: in_sedName = .false.
    logical, save :: in_fuzzyfctName = .false.
    logical, save :: in_varNb = .false.
    logical, save :: in_mbsfctName = .false.

    character(len=128), save :: mbsVar, xyPters, fzzVar, coPters

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

        if (name=='Carbonates-Organics') in_CarbOrgs = .true.
        if(in_CarbOrgs) carborgparamsection = .true.
        if(in_CarbOrgs) call ScarborgElement_handler( name )

        if( name=='membershipFunction' ) in_mbsfct = .true.
        if( in_mbsfct ) mbsfctsection = .true.
        if( in_mbsfct ) call SmbsfctElement_handler( name )
        if( in_fctSubset ) call SmbsfctsubsetElement_handler( name )

        if (name=='fuzzyFunction') in_fuzzyfct = .true.
        if( in_fuzzyfct ) fuzzyfctsection = .true.
        if( in_fuzzyfct ) call SfuzzyfctElement_handler( name )
        if( in_fuzzySubset ) call SfuzzysubsetElement_handler( name )

        if (name=='fuzzyRule') in_fuzzyrule = .true.
        if( in_fuzzyrule ) fuzzyrulesection = .true.
        if( in_fuzzyrule ) call SfuzzyrlElement_handler( name )
        if( in_mbsSet ) call SmbsfuzzyElement_handler( name )

    end subroutine startElement_handler
    !============================================================================
    subroutine endElement_handler(namespaceURI, localname, name)

        character(len=*), intent(in) :: namespaceURI
        character(len=*), intent(in) :: localname
        character(len=*), intent(in) :: name

        call EcarborgElement_handler( name )

        call EmbsfctElement_handler( name )
        call EmbsfctsubsetElement_handler( name )

        call EfuzzyfctElement_handler( name )
        call EfuzzysubsetElement_handler( name )

        call EfuzzyrlElement_handler( name )
        call EmbsfuzzyElement_handler( name )

    end subroutine endElement_handler
    ! ============================================================================
    subroutine ScarborgElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='carborgTime') in_carborgTime = .true.
        if (name=='tempsalFile') in_tempsalFile = .true.
        if (name=='membershipfunctionNb') in_mbsfctNb = .true.
        if (name=='fuzzyfunctionNb') in_fuzzyfctNb = .true.
        if (name=='fuzzyruleNb') in_fuzzyruleNb = .true.

    end subroutine ScarborgElement_handler
    ! ============================================================================
    subroutine SmbsfctElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='functionName') in_fctName = .true.
        if (name=='functionVariable') in_fctVariable = .true.
        if (name=='sedName') in_sedName= .true.
        if (name=='functionPts') in_fctPts = .true.
        if (name=='functionSubset') in_fctSubset = .true.

    end subroutine SmbsfctElement_handler
    ! ============================================================================
    subroutine SmbsfctsubsetElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='xypoints') in_mbsfctPoints = .true.

    end subroutine SmbsfctsubsetElement_handler
    ! ============================================================================
    subroutine SfuzzyfctElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='fuzzyName') in_fuzzyName = .true.
        if (name=='fuzzyVariable') in_fuzzyVariable = .true.
        if (name=='fuzzyPts') in_fuzzyPts = .true.
        if (name=='fuzzySubset') in_fuzzySubset = .true.

    end subroutine SfuzzyfctElement_handler
    ! ============================================================================
    subroutine SfuzzysubsetElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='xycoords') in_fuzzyCoords = .true.

    end subroutine SfuzzysubsetElement_handler
    ! ============================================================================
    subroutine SfuzzyrlElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='matName') in_matName = .true.
        if (name=='fuzzyFctName') in_fuzzyfctName = .true.
        if (name=='variableNb') in_varNb = .true.
        if (name=='membershipSet') in_mbsSet = .true.

    end subroutine SfuzzyrlElement_handler
    ! ============================================================================
    subroutine SmbsfuzzyElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='membershipFctName') in_mbsfctName = .true.

    end subroutine SmbsfuzzyElement_handler
     ! ============================================================================
    subroutine EcarborgElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='carborgTime') in_carborgTime = .false.
        if (name=='tempsalFile') in_tempsalFile = .false.
        if (name=='membershipfunctionNb') in_mbsfctNb = .false.
        if (name=='fuzzyfunctionNb') in_fuzzyfctNb = .false.
        if (name=='fuzzyruleNb') in_fuzzyruleNb = .false.

    end subroutine EcarborgElement_handler
    ! ============================================================================
    subroutine EmbsfctElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='functionName') in_fctName = .false.
        if (name=='functionVariable') in_fctVariable = .false.
        if (name=='sedName') in_sedName= .false.
        if (name=='functionPts') in_fctPts = .false.
        if (name=='functionSubset') in_fctSubset = .false.

    end subroutine EmbsfctElement_handler
    ! ============================================================================
    subroutine EmbsfctsubsetElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='xypoints') in_mbsfctPoints = .false.

    end subroutine EmbsfctsubsetElement_handler
    ! ============================================================================
    subroutine EfuzzyfctElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='fuzzyName') in_fuzzyName = .false.
        if (name=='fuzzyVariable') in_fuzzyVariable = .false.
        if (name=='fuzzyPts') in_fuzzyPts = .false.
        if (name=='fuzzySubset') in_fuzzySubset = .false.

    end subroutine EfuzzyfctElement_handler
    ! ============================================================================
    subroutine EfuzzysubsetElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='xycoords') in_fuzzyCoords = .false.

    end subroutine EfuzzysubsetElement_handler
    ! ============================================================================
    subroutine EfuzzyrlElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='matName') in_matName = .false.
        if (name=='fuzzyFctName') in_fuzzyfctName = .false.
        if (name=='variableNb') in_varNb = .false.
        if (name=='membershipSet') in_mbsSet = .false.

    end subroutine EfuzzyrlElement_handler
    ! ============================================================================
    subroutine EmbsfuzzyElement_handler(name)

        character(len=*), intent(in) :: name

        if (name=='membershipFctName') in_mbsfctName = .false.

    end subroutine EmbsfuzzyElement_handler
    ! ============================================================================
    subroutine characters_handler(chars)

        character(len=*), intent(in) :: chars

        if(in_CarbOrgs) call carborgs_characters_handler(chars)

        if( in_mbsfct ) call membership_characters_handler( chars )
        if( in_fctSubset ) call membershipset_characters_handler( chars )

        if( in_fuzzyfct ) call fuzzyfct_characters_handler( chars )
        if( in_fuzzySubset ) call fuzzysubset_characters_handler( chars )

        if( in_fuzzyrule ) call fuzzyrule_characters_handler( chars )
        if( in_mbsSet ) call mbsfuzzy_characters_handler( chars )

    end subroutine characters_handler
    ! ============================================================================
    subroutine carborgs_characters_handler( chars )

        character(len=*), intent(in) :: chars

        if( in_carborgTime )then
            call rts(chars,time_carborg)
        elseif( in_tempsalFile )then
            ftempsal = chars
        elseif( in_mbsfctNb )then
            call rts( chars, membership_fct_nb )
            allocate( membership( membership_fct_nb ) )
        elseif( in_fuzzyfctNb )then
            call rts( chars, fuzzy_fct_nb )
            allocate( fuzzyfct( fuzzy_fct_nb ) )
        elseif( in_fuzzyruleNb )then
            call rts( chars, fuzzy_rl_nb )
            allocate( fuzzyfy( fuzzy_rl_nb ) )
        endif

    end subroutine carborgs_characters_handler
    ! ============================================================================
    subroutine membership_characters_handler( chars )

        integer :: id, k

        character(len=*), intent(in) :: chars

        if( in_fctName )then
            mbshp_Nb = mbshp_Nb + 1
            pts_Nb = 0
            membership( mbshp_Nb )%name = chars
            membership( mbshp_Nb )%sedclassID = -1
        elseif( in_fctVariable )then
            mbsVar = chars
            if( trim( mbsVar ) == 'depth' )then
                membership( mbshp_Nb )%variableID = 1
            elseif( trim( mbsVar ) == 'ocean-energy' )then
                membership( mbshp_Nb )%variableID = 2
            elseif( trim( mbsVar ) == 'shore-distance' )then
                membership( mbshp_Nb )%variableID = 3
            elseif( trim( mbsVar ) == 'river-distance' )then
                membership( mbshp_Nb )%variableID = 4
            elseif( trim( mbsVar ) == 'material-distance' )then
                membership( mbshp_Nb )%variableID = 5
            elseif( trim( mbsVar ) == 'sedimentation-rate' )then
                membership( mbshp_Nb )%variableID = 6
            elseif( trim( mbsVar ) == 'temperature' )then
                membership( mbshp_Nb )%variableID = 7
            elseif( trim( mbsVar ) == 'salinity' )then
                membership( mbshp_Nb )%variableID = 8
            elseif( trim( mbsVar ) == 'valley' )then
                membership( mbshp_Nb )%variableID = 9
            elseif( trim( mbsVar ) == 'gradient' )then
                membership( mbshp_Nb )%variableID = 10
            elseif( trim( mbsVar ) == 'exposed-time' )then
                membership( mbshp_Nb )%variableID = 11
            elseif( trim( mbsVar ) == 'burial' )then
                membership( mbshp_Nb )%variableID = 12
            elseif( trim( mbsVar ) == 'age' )then
                membership( mbshp_Nb )%variableID = 13
            else
                print*,'Something went wrong in the declaration of the membershipFunction.'
                print*,'the functionVariable name does not match with predefined variables.'
                print*,mbshp_Nb, trim( mbsVar )
                stop
            endif
        elseif( in_sedName )then
            id = 0
            loop_over_sed_class1: do k = 1, totgrn
                if( trim( chars ) == trim( material_name( k ) ) )then
                    id = k
                    exit loop_over_sed_class1
                endif
            enddo loop_over_sed_class1
            if( id == 0 )then
                print*,'Something went wrong in the declaration of the membership function.'
                print*,'the sedName element does not match with defined material name.'
                print*,mbshp_Nb, trim( chars )
                stop
            endif
            membership( mbshp_Nb )%sedclassID = id
        elseif( in_fctPts )then
            call rts( chars, membership( mbshp_Nb )%nb_points )
        endif

    end subroutine membership_characters_handler
    ! ============================================================================
    subroutine membershipset_characters_handler( chars )

        character(len=*), intent(in) :: chars

        if( in_mbsfctPoints )then
            pts_Nb = pts_Nb + 1
            xyPters = chars
            call rts( xyPters, membership( mbshp_Nb )%xypoints( pts_Nb, 1:2 ) )
        endif

    end subroutine membershipset_characters_handler
    ! ============================================================================
    subroutine fuzzyfct_characters_handler( chars )

        character(len=*), intent(in) :: chars

        if( in_fuzzyName )then
            fuzz_Nb = fuzz_Nb + 1
            coords_Nb = 0
            fuzzyfct( fuzz_Nb )%name = chars
        elseif( in_fuzzyVariable )then
            fzzVar = chars
            if( trim( fzzVar ) == 'growth' )then
                fuzzyfct( fuzz_Nb )%variableID = 1
            elseif( trim( fzzVar ) == 'erosion' )then
                fuzzyfct( fuzz_Nb )%variableID = 2
            elseif( trim( fzzVar ) == 'porosity' )then
                fuzzyfct( fuzz_Nb )%variableID = 3
            elseif( trim( fzzVar ) == 'hardness' )then
                fuzzyfct( fuzz_Nb )%variableID = 4
            else
                print*,'Something went wrong in the declaration of the fuzzyFunction.'
                print*,'the fuzzyVariable name does not match with predefined variables.'
                print*,fuzz_Nb, trim( fzzVar )
                stop
            endif
        elseif( in_fuzzyPts )then
            call rts( chars, fuzzyfct( fuzz_Nb )%nb_points )
        endif

    end subroutine fuzzyfct_characters_handler
    ! ============================================================================
    subroutine fuzzysubset_characters_handler( chars )

        character(len=*), intent(in) :: chars

        if( in_fuzzyCoords )then
            coords_Nb = coords_Nb + 1
            coPters = chars
            call rts( coPters, fuzzyfct( fuzz_Nb )%xypoints( coords_Nb, 1:2 ) )
        endif

    end subroutine fuzzysubset_characters_handler
    ! ============================================================================
    subroutine fuzzyrule_characters_handler( chars )

        integer :: k, id

        character(len=*), intent(in) :: chars

        if( in_matName )then
            fuzzy_Nb = fuzzy_Nb + 1
            mbs_Nb = 0
            id = 0
            loop_over_sed_class: do k = 1, totgrn
                if( trim( chars ) == trim( material_name( k ) ) )then
                    id = k
                    exit loop_over_sed_class
                endif
            enddo loop_over_sed_class
            if( id == 0 )then
                print*,'Something went wrong in the declaration of the fuzzyRule.'
                print*,'the matName element does not match with defined material name.'
                print*,fuzzy_Nb, trim( chars )
                stop
            endif
            fuzzyfy( fuzzy_Nb )%sedclassID = id
        elseif( in_fuzzyfctName )then
            id = 0
            loop_over_fuzzy_fct: do k = 1, fuzzy_fct_nb
                if( trim( chars ) == trim( fuzzyfct( k )%name ) )then
                    id = k
                    exit loop_over_fuzzy_fct
                endif
            enddo loop_over_fuzzy_fct
            if( id == 0 )then
                print*,'Something went wrong in the declaration of the fuzzyRule.'
                print*,'the fuzzyFctName element does not match with defined fuzzy name.'
                print*,fuzzy_Nb, trim( chars )
                stop
            endif
            fuzzyfy( fuzzy_Nb )%fuzzyID = id
        elseif( in_varNb )then
            call rts( chars, fuzzyfy( fuzzy_Nb )%membershipNb )
        endif

    end subroutine fuzzyrule_characters_handler
    ! ============================================================================
    subroutine mbsfuzzy_characters_handler( chars )

        integer :: k, id

        character(len=*), intent(in) :: chars

        if( in_mbsfctName )then
            mbs_Nb = mbs_Nb + 1
            id = 0
            loop_over_mbs_fct: do k = 1, membership_fct_nb
                if( trim( chars ) == trim( membership( k )%name ) )then
                    id = k
                    exit loop_over_mbs_fct
                endif
            enddo loop_over_mbs_fct
            if( id == 0 )then
                print*,'Something went wrong in the declaration of the fuzzyRule.'
                print*,'the membershipFctName element does not match with defined membership name.'
                print*,fuzzy_Nb, mbs_Nb, trim( chars )
                stop
            endif
            fuzzyfy( fuzzy_Nb )%membershipID( mbs_Nb ) = id
        endif

    end subroutine mbsfuzzy_characters_handler
    ! ============================================================================
    subroutine endDocument_handler

    end subroutine endDocument_handler
    ! ============================================================================
    !> Subroutine xml_CarbOrgs()
    !! reads input that defines the Carbonate / Organics parameters.
    !<
    ! ============================================================================
    subroutine xml_CarbOrgs

        type(xml_t) :: xf
        integer :: ios

        mbshp_Nb = 0
        fuzz_Nb = 0
        fuzzy_Nb = 0
        pts_Nb = 0
        coords_Nb = 0
        mbs_Nb = 0

        ! Open file
        call open_xml_file(xf, finput , ios)
        if (ios/=0) then
            print*,'---------------------'
            print*, 'Error opening input file for parsing XmL'
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

    end subroutine xml_CarbOrgs
     ! ============================================================================
    !> Subroutine broadcast_carborgs()
    !! Broadcast initial parameters which define the carbonates organics variables.
    !<
    ! ============================================================================
    subroutine broadcast_carborgs

        integer :: k, m, flag

        real( tkind ), dimension( 20, 2 ) :: temp

        call mpi_bcast( time_carborg,1,dbl_type,0,SPModel_comm_world,ierr )
        call mpi_bcast( ftempsal,128,mpi_character,0,SPModel_comm_world,ierr )

        call mpi_bcast( membership_fct_nb,1,int_type,0,SPModel_comm_world,ierr )
        call mpi_bcast( fuzzy_fct_nb,1,int_type,0,SPModel_comm_world,ierr )
        call mpi_bcast( fuzzy_rl_nb,1,int_type,0,SPModel_comm_world,ierr )

        if( iam /= 0 )then
            allocate( membership( membership_fct_nb ) )
            allocate( fuzzyfct( fuzzy_fct_nb ) )
            allocate( fuzzyfy( fuzzy_rl_nb ) )
        endif
		allocate( fuzzycomb( 4, totgrn ) )

        do k = 1, membership_fct_nb
            call mpi_bcast( membership( k )%name,128,mpi_character,0,SPModel_comm_world,ierr )
            call mpi_bcast( membership( k )%variableID,1,int_type,0,SPModel_comm_world,ierr )
            call mpi_bcast( membership( k )%sedclassID,1,int_type,0,SPModel_comm_world,ierr )
            call mpi_bcast( membership( k )%nb_points,1,int_type,0,SPModel_comm_world,ierr )
            do m = 1, membership( k )%nb_points
                call mpi_bcast( membership( k )%xypoints( m, 1:2 ),2,dbl_type,0,SPModel_comm_world,ierr )
            enddo
        enddo

        do k = 1, fuzzy_fct_nb
            call mpi_bcast( fuzzyfct( k )%name,128,mpi_character,0,SPModel_comm_world,ierr )
            call mpi_bcast( fuzzyfct( k )%variableID,1,int_type,0,SPModel_comm_world,ierr )
            call mpi_bcast( fuzzyfct( k )%nb_points,1,int_type,0,SPModel_comm_world,ierr )
            do m = 1, fuzzyfct( k )%nb_points
                call mpi_bcast( fuzzyfct( k )%xypoints( m, 1:2 ),2,dbl_type,0,SPModel_comm_world,ierr )
            enddo
        enddo

        do k = 1, fuzzy_rl_nb
            call mpi_bcast( fuzzyfy( k )%sedclassID,1,int_type,0,SPModel_comm_world,ierr )
            call mpi_bcast( fuzzyfy( k )%fuzzyID,1,int_type,0,SPModel_comm_world,ierr )
            call mpi_bcast( fuzzyfy( k )%membershipNb,1,int_type,0,SPModel_comm_world,ierr )
            do m = 1, fuzzyfy( k )%membershipNb
                call mpi_bcast( fuzzyfy( k )%membershipID( m ),1,int_type,0,SPModel_comm_world,ierr )
            enddo
        enddo

        call mpi_barrier( SPModel_comm_world,ierr )

        ! Check which flags are on.
        carb_sediment_flag = 0
        carb_exposed_flag = 0
        mat_nearest_flag = 0
        temperature_flag = 0
        carb_valley_flag = 0
        carb_burial_flag = 0
        carb_depth_flag = 0
        carb_ocean_flag = 0
        carb_shore_flag = 0
        carb_slope_flag = 0
        carb_river_flag = 0
        carb_age_flag = 0
        salinity_flag = 0

        do k = 1, membership_fct_nb
            flag = membership( k )%variableID
            ! Depth
            if( flag == 1 )then
                carb_depth_flag = 1
            ! Ocean-energy
            elseif( flag == 2 )then
                carb_ocean_flag = 1
            ! Shore-distance
            elseif( flag == 3 )then
                carb_shore_flag = 1
            ! River-distance
            elseif( flag == 4 )then
                carb_river_flag = 1
            ! Material-distance
            elseif( flag == 5 )then
                mat_nearest_flag = 1
            ! Sedimentation-rate
            elseif( flag == 6 )then
                carb_sediment_flag = 1
            ! Temperature
            elseif( flag == 7 )then
                temperature_flag = 1
            ! Salinity
            elseif( flag == 8 )then
                salinity_flag = 1
            ! Valley
            elseif( flag == 9 )then
                carb_valley_flag = 1
            ! Gradient
            elseif( flag == 10 )then
                carb_slope_flag = 1
            ! Exposed-time
            elseif( flag == 11 )then
                carb_exposed_flag = 1
            ! Burial
            elseif( flag == 12 )then
                carb_burial_flag = 1
            ! Age
            elseif( flag == 13 )then
                carb_age_flag = 1
            endif

            ! Force the membership function to form a closed set.

            ! Check lowest value
            if( membership( k )%xypoints( 1, 2 ) /= 0.0_8 )then
                  temp( 1, 1 ) = membership( k )%xypoints( 1, 1 ) - tor
                  temp( 1, 2 ) =  membership( k )%xypoints( 1, 2 )

                  temp( 1, 1 ) = membership( k )%xypoints( 1, 1 ) - 1.e10_8 - 1
                  temp( 1, 2 ) =  0.0_8
                  temp( 2, 1 ) = membership( k )%xypoints( 1, 1 ) - 1.e10_8
                  temp( 2, 2 ) =  membership( k )%xypoints( 1, 2 )
                  do m = 1, membership( k )%nb_points
                    temp( m + 2, 1 ) = membership( k )%xypoints( m, 1 )
                    temp( m + 2, 2 ) = membership( k )%xypoints( m, 2 )
                  enddo
                  membership( k )%nb_points = membership( k )%nb_points + 2
                  if( membership( k )%nb_points > 20 )then
                    print*,'Warning number of points exceeds max allowed for membership: ', &
                        trim( membership( k )%name )
                    stop
                  endif
                  do m = 1, membership( k )%nb_points
                    membership( k )%xypoints( m, 1 ) = temp( m, 1 )
                    membership( k )%xypoints( m, 2 ) = temp( m, 2 )
                  enddo
            endif

            ! Check highest value
            if( membership( k )%xypoints( membership( k )%nb_points, 2 ) /= 0.0_8 )then
                  m = membership( k )%nb_points + 2
                  if( m > 20 )then
                    print*,'Warning number of points exceeds max allowed for membership: ', &
                        trim( membership( k )%name )
                    stop
                  endif
                  membership( k )%nb_points = m
                  membership( k )%xypoints( m-1, 1 ) = membership( k )%xypoints( m-2, 1 ) + 1.e10_8 + 1
                  membership( k )%xypoints( m-1, 2 ) = membership( k )%xypoints( m-2, 2 )
                  membership( k )%xypoints( m, 1 ) = membership( k )%xypoints( m-1, 1 ) + 1
                  membership( k )%xypoints( m, 2 ) = 0.0_8
            endif

        enddo

        ! Perform additional declarations for fuzzy rules
        do k = 1, fuzzy_fct_nb

            ! Force the fuzzy function to form a closed set.

            ! Check lowest value
            if( fuzzyfct( k )%xypoints( 1, 2 ) /= 0.0_8 )then
                  temp( 1, 1 ) = fuzzyfct( k )%xypoints( 1, 1 ) - tor
                  temp( 1, 2 ) =  0.0_8
                  do m = 1, fuzzyfct( k )%nb_points
                    temp( m + 1, 1 ) = fuzzyfct( k )%xypoints( m, 1 )
                    temp( m + 1, 2 ) = fuzzyfct( k )%xypoints( m, 2 )
                  enddo
                  fuzzyfct( k )%nb_points = fuzzyfct( k )%nb_points + 1
                  if( fuzzyfct( k )%nb_points > 20 )then
                    print*,'Warning number of points exceeds max allowed for fuzzy function: ', &
                        trim( fuzzyfct( k )%name )
                    stop
                  endif
                  do m = 1, fuzzyfct( k )%nb_points
                    fuzzyfct( k )%xypoints( m, 1 ) = temp( m, 1 )
                    fuzzyfct( k )%xypoints( m, 2 ) = temp( m, 2 )
                  enddo
            endif

            ! Check highest value
            if( fuzzyfct( k )%xypoints( fuzzyfct( k )%nb_points, 2 ) /= 0.0_8 )then
                  m = fuzzyfct( k )%nb_points + 1
                  if( m > 20 )then
                    print*,'Warning number of points exceeds max allowed for fuzzy function: ', &
                        trim( fuzzyfct( k )%name )
                    stop
                  endif
                  fuzzyfct( k )%nb_points = m
                  fuzzyfct( k )%xypoints( m, 1 ) = fuzzyfct( k )%xypoints( m-1, 1 ) + tor
                  fuzzyfct( k )%xypoints( m, 2 ) = 0.0_8
            endif

        enddo

        return

    end subroutine broadcast_carborgs
   ! ============================================================================

end module carborg_parser
