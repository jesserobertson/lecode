! ============================================================================
! Name        : FlowWalkerO.f90
! Author      : Tristan Salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file FlowWalkerO.f90
!
! Description : FlowWalkerO is used to generate Hdf5 output of the flow walkers during SPModel simulation.
!
!<
! ============================================================================
module FWLK_out

    use hdf5
    use file_data
    use mpidata
    use time_data
    use FoX_wxml
    use fwalker_data

    implicit none

    public

contains
    ! ============================================================================
    !> Subroutine FWLK_hdf5
    !! Subroutine FWLK_hdf5 uses Hdf5 to output SP Model flow walkers evolution.
    !! It records both connectivity and vertices.
    !<
    ! ============================================================================
    subroutine FWLK_hdf5( iter )

        logical :: compression

        ! Parameters Declaration
        integer :: id, i, rank, iter

        character(len=128) :: stg, file

        integer, dimension( : ), allocatable:: ftype
        real(tkind), dimension( : ), allocatable :: sedcharge, nodes, vels, hgt, wdt

        integer(hid_t) :: file_id, plist_id
        integer(hid_t) :: filespace, dset_id
        integer(hsize_t),dimension(2) :: dims

        allocate( ftype( num_fw_o ), nodes( num_fw_o*3 ), vels( num_fw_o*3 ) )
        allocate( sedcharge( num_fw_o ), hgt( num_fw_o ), wdt( num_fw_o ) )

        file = 'FWpoints.'
        call noblnk( file )
        call append_zero( file,iter )
        stg = '.p'
        call append_str( file,stg )
        call append_nb( file,iam )
        stg = '.h5'
        call append_str( file,stg )
        call addpath1( file )
        ! Create nodes and connectivity arrays
        id = 1
        do i = 1, num_fw_o
            nodes( id ) = recfw_xpos( i )
            nodes( id + 1 ) = recfw_ypos( i )
            nodes( id + 2 ) = recfw_zpos( i )
            if( udw ) nodes( id + 1 ) = recfw_zpos( i )
            if( udw ) nodes( id + 2 ) = recfw_ypos( i )
            id = id + 3
        enddo
        ! Create attributes arrays
        ! Velocity fields
        id = 1
        do i = 1, num_fw_o
            vels( id ) = recfw_xvel( i )
            vels( id + 1 ) = recfw_yvel( i )
            vels( id + 2 ) = 0.0
            id = id + 3
        enddo
        ! Flow type
        id = 1
        do i = 1, num_fw_o
            ftype( id ) = recfw_type( i )
            id = id + 1
        enddo
        ! Flow height
        id = 1
        do i = 1, num_fw_o
            hgt( id ) = recfw_srch( i )
            id = id + 1
        enddo
        ! Flow width
        id = 1
        do i = 1, num_fw_o
            wdt( id ) = recfw_width( i )
            id = id + 1
        enddo
        ! Sediment load concentration
        id = 1
        do i = 1, num_fw_o
            sedcharge( id ) = recfw_sedcharge( i )
            id = id + 1
        enddo
        ! Initialize predefined datatypes
        call h5open_f( ierr )
        call h5zfilter_avail_f( h5z_filter_deflate_f, compression, ierr )

        ! Setup file access property list for MPI-IO access.
        call h5pcreate_f( h5p_file_access_f, plist_id, ierr )
        ! Create the file collectively.
        call h5fcreate_f( file, h5f_acc_trunc_f, file_id, ierr, access_prp = plist_id )
        ! The Coordinates - vertices
        dims( 1 ) = 3
        dims( 2 ) = num_fw_o
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        if( .not. compression .or. dims( 2 ) == 0 )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, "/vertices", h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, nodes, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id,"/vertices", h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, nodes, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Velocity field attribute
        ! ========================
        dims( 1 ) = 3
        dims( 2 ) = num_fw_o
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )
        if( .not. compression .or. dims( 2 ) == 0 )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, "/FW_velocity", h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_double, vels, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id,"/FW_velocity", h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, vels, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Type field attribute
        ! ========================
        dims( 1 ) = 1
        dims( 2 ) = num_fw_o
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )

        if( .not. compression .or. dims( 2 ) == 0 )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, "/FW_type", h5t_native_integer, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the datase
            call h5dwrite_f( dset_id, h5t_native_integer, ftype, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id,"/FW_type", h5t_native_integer , filespace, dset_id, ierr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_integer, ftype, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Flow height attribute
        ! ========================
        dims( 1 ) = 1
        dims( 2 ) = num_fw_o
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )

        if( .not. compression .or. dims( 2 ) == 0 )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, "/FW_height", h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the datase
            dims( 1 ) = 1
            dims( 2 ) = num_fw_o
            call h5dwrite_f( dset_id, h5t_native_double, hgt, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id,"/FW_height", h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, hgt, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Flow width attribute
        ! ========================
        dims( 1 ) = 1
        dims( 2 ) = num_fw_o
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )

        if( .not. compression .or. dims( 2 ) == 0 )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, "/FW_width", h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the datase
            dims( 1 ) = 1
            dims( 2 ) = num_fw_o
            call h5dwrite_f( dset_id, h5t_native_double, wdt, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id,"/FW_width", h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, wdt, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! ========================
        ! Sediment load field attribute
        ! ========================
        dims( 1 ) = 1
        dims( 2 ) = num_fw_o
        rank = 2
        call h5screate_simple_f( rank, dims, filespace, ierr )

        if( .not. compression .or. dims( 2 ) == 0 )then
            ! Create the dataset with default properties
            call h5dcreate_f( file_id, "/FW_sedcharge", h5t_native_double, filespace, dset_id, ierr )
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_xfer_f, plist_id, ierr )
            ! Write the datase
            dims( 1 ) = 1
            dims( 2 ) = num_fw_o
            call h5dwrite_f( dset_id, h5t_native_double, sedcharge, dims, ierr, &
                file_space_id = filespace, xfer_prp = plist_id )
        else
            ! Create property list for collective dataset write
            call h5pcreate_f( h5p_dataset_create_f, plist_id, ierr )
            call h5pset_deflate_f( plist_id, 9, ierr )
            call h5pset_chunk_f( plist_id, rank, dims, ierr )
            ! Create the dataset with default properties
            call h5dcreate_f( file_id,"/FW_sedcharge", h5t_native_double , filespace, dset_id, ierr, plist_id )
            ! Write the dataset collectively
            call h5dwrite_f( dset_id, h5t_native_double, sedcharge, dims, ierr )
            call h5pclose_f( plist_id, ierr )
        endif
        ! Close the dataset
        call h5dclose_f( dset_id, ierr )
        call h5sclose_f( filespace, ierr )

        ! Close the file.
        call h5fclose_f( file_id, ierr )
        ! Close interface
        call h5close_f( ierr )

        deallocate( ftype, sedcharge, nodes, vels, hgt, wdt )

        return

    end subroutine FWLK_hdf5
    ! ============================================================================
    !> Subroutine FW_xmf
    !! Subroutine FW_xmf generates the FW XdmF file for the considered time simulation.
    !! It builds the file which calls all the hdf5 created.
    !<
    ! ============================================================================
    subroutine FW_xmf( iter )

        ! Parameters Declaration
        type(xmlf_t) :: xf

        integer :: iter, totnodes, k
        character(len=128) :: str, filename, filename1, filename2, filename3, filename4, file, h5file, stg
        character(len=128) :: filename5, filename6

        integer,dimension( nproc ) :: Fproc

        totnodes = num_fw_o
        Fproc = 0
        call mpi_gather(totnodes,1,int_type,Fproc,1,int_type,0,SPModel_comm_world,ierr)
        if( iam == 0 )then
            file = 'FWpoints.'
            call noblnk( file )
            call append_zero( file,iter )
            str = '.xmf'
            call append_str( file,str )
            call addpath1( file )
            call xml_OpenFile( file, xf )
            ! Header
            call xml_AddDOCTYPE(xf, "Xdmf", "Xdmf.dtd")
            call xml_DeclareNamespace(xf, "http://www.w3.org/2001/XInclude", "xi")
            call xml_NewElement(xf, "Xdmf" )
            call xml_AddAttribute(xf, "Version", "2.0" )
            call xml_NewElement( xf, "Domain" )
            call xml_NewElement( xf, "Grid" )
            call xml_AddAttribute( xf, "Name", "FWRecord" )
            call xml_AddAttribute( xf, "GridType", "Collection" )
            call xml_NewElement( xf, "Time" )
            call xml_AddAttribute( xf, "Type", "Single" )
            call xml_AddAttribute( xf, "Value", next_display-time_display )
            call xml_EndElement( xf, "Time" )
            ! Loop over processors
            do k = 1, nproc
                if( Fproc( k ) > 0 )then
                    h5file = ''
                    h5file = 'FWpoints'
                    call noblnk( h5file )
                    stg = '.'
                    call append_str( h5file,stg )
                    call append_zero( h5file,iter )
                    stg = '.p'
                    call append_str( h5file,stg )
                    call append_nb( h5file,k-1 )
                    stg = '.h5'
                    call append_str( h5file,stg )
                    filename = h5file
                    call noblnk( filename )
                    filename1 = filename
                    str = ':/vertices'
                    call append_str( filename1,str )
                    filename2 = filename
                    str = ':/FW_velocity'
                    call append_str( filename2,str )
                    filename3 = filename
                    str = ':/FW_type'
                    call append_str( filename3,str )
                    filename4 = filename
                    str = ':/FW_height'
                    call append_str( filename4,str )
                    filename5 = filename
                    str = ':/FW_width'
                    call append_str( filename5,str )
                    filename6 = filename
                    str = ':/FW_sedcharge'
                    call append_str( filename6,str )
                    ! Block begin
                    call xml_NewElement( xf, "Grid" )
                    str = 'FWBlock.'
                    call append_zero( str,iter )
                    stg = '.p'
                    call append_str( str,stg )
                    call append_nb( str,k - 1 )
                    call xml_AddAttribute( xf, "Name", trim( str ) )
                    call xml_NewElement( xf, "Topology" )
                    call xml_AddAttribute( xf, "Type", "POLYVERTEX" )
                    call xml_AddAttribute( xf, "NumberOfElements", Fproc( k ) )
                    call xml_EndElement( xf, "Topology" )
                    ! Geometry
                    call xml_NewElement( xf, "Geometry" )
                    call xml_AddAttribute( xf, "Type", "XYZ" )
                    call xml_NewElement( xf, "DataItem" )
                    call xml_AddAttribute( xf, "Format", "HDF" )
                    call xml_AddAttribute( xf, "NumberType", "Float" )
                    call xml_AddAttribute( xf, "Precision", "8" )
                    str = ' '
                    call append_nb2( str,Fproc( k ) )
                    call append_nb2( str, 3 )
                    call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                    call xml_AddCharacters( xf, trim( filename1 ) )
                    call xml_EndElement( xf, "DataItem" )
                    call xml_EndElement( xf, "Geometry" )
                    ! Velocity
                    call xml_NewElement( xf, "Attribute" )
                    call xml_AddAttribute( xf, "Type", "Vector" )
                    call xml_AddAttribute( xf, "Center", "Node" )
                    call xml_AddAttribute( xf, "Name", "XYVelocity" )
                    call xml_NewElement( xf, "DataItem" )
                    call xml_AddAttribute( xf, "Format", "HDF" )
                    call xml_AddAttribute( xf, "NumberType", "Float" )
                    call xml_AddAttribute( xf, "Precision", "8" )
                    call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                    call xml_AddCharacters( xf, trim(filename2 ) )
                    call xml_EndElement( xf, "DataItem" )
                    call xml_EndElement( xf, "Attribute" )
                    ! Flow type
                    call xml_NewElement( xf, "Attribute" )
                    call xml_AddAttribute( xf, "Type", "Scalar" )
                    call xml_AddAttribute( xf, "Center", "Node" )
                    call xml_AddAttribute( xf, "Name", "FlowType" )
                    call xml_NewElement( xf, "DataItem" )
                    call xml_AddAttribute( xf, "Format", "HDF" )
                    call xml_AddAttribute( xf, "NumberType", "Int" )
                    str = ' '
                    call append_nb2( str, Fproc( k ) )
                    call append_nb2( str, 1 )
                    call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                    call xml_AddCharacters( xf, trim( filename3 ) )
                    call xml_EndElement( xf, "DataItem" )
                    call xml_EndElement( xf, "Attribute" )
                    ! Flow height
                    call xml_NewElement( xf, "Attribute" )
                    call xml_AddAttribute( xf, "Type", "Scalar" )
                    call xml_AddAttribute( xf, "Center", "Node" )
                    call xml_AddAttribute( xf, "Name", "Flow Height" )
                    call xml_NewElement( xf, "DataItem" )
                    call xml_AddAttribute( xf, "Format", "HDF" )
                    call xml_AddAttribute( xf, "NumberType", "Float" )
                    str = ' '
                    call append_nb2( str, Fproc( k ) )
                    call append_nb2( str, 1 )
                    call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                    call xml_AddCharacters( xf, trim( filename4 ) )
                    call xml_EndElement( xf, "DataItem" )
                    call xml_EndElement( xf, "Attribute" )
                    ! Flow width
                    call xml_NewElement( xf, "Attribute" )
                    call xml_AddAttribute( xf, "Type", "Scalar" )
                    call xml_AddAttribute( xf, "Center", "Node" )
                    call xml_AddAttribute( xf, "Name", "Stream Width" )
                    call xml_NewElement( xf, "DataItem" )
                    call xml_AddAttribute( xf, "Format", "HDF" )
                    call xml_AddAttribute( xf, "NumberType", "Float" )
                    str = ' '
                    call append_nb2( str, Fproc( k ) )
                    call append_nb2( str, 1 )
                    call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                    call xml_AddCharacters( xf, trim( filename5 ) )
                    call xml_EndElement( xf, "DataItem" )
                    call xml_EndElement( xf, "Attribute" )
                    ! Sediment load
                    call xml_NewElement( xf, "Attribute" )
                    call xml_AddAttribute( xf, "Type", "Scalar" )
                    call xml_AddAttribute( xf, "Center", "Node" )
                    call xml_AddAttribute( xf, "Name", "Sed. Conc. kg/yr" )
                    call xml_NewElement( xf, "DataItem" )
                    call xml_AddAttribute( xf, "Format", "HDF" )
                    call xml_AddAttribute( xf, "NumberType", "Float" )
                    str = ' '
                    call append_nb2( str, Fproc( k ) )
                    call append_nb2( str, 1 )
                    call xml_AddAttribute( xf, "Dimensions", trim( str ) )
                    call xml_AddCharacters( xf, trim( filename6 ) )
                    call xml_EndElement( xf, "DataItem" )
                    call xml_EndElement( xf, "Attribute" )
                    ! Block end
                    call xml_EndElement( xf, "Grid" )
                endif
            enddo
            ! Footer
            call xml_EndElement( xf, "Grid" )
            call xml_EndElement( xf, "Domain" )
            call xml_EndElement( xf, "Xdmf" )
            call xml_Close( xf )
        endif

        return

    end subroutine FW_xmf
    ! ============================================================================
    !> Subroutine FW_series
    !! Subroutine FW_series generates the XmL file for flow walkers.
    !<
    ! ============================================================================
    subroutine FW_series( iter )

        ! Parameters Declaration
        type(xmlf_t) :: xf

        integer :: i, iter, it0
        character(len=128) :: filename, str, fname

        filename = 'FW_series.xdmf'
        call addpath1(filename)
        ! Header
        call xml_OpenFile( filename, xf )
        call xml_AddDOCTYPE(xf, "Xdmf", "Xdmf.dtd")
        call xml_DeclareNamespace(xf, "http://www.w3.org/2001/XInclude", "xi")
        call xml_NewElement(xf, "Xdmf" )
        call xml_AddAttribute(xf, "Version", "2.0" )
        call xml_NewElement( xf, "Domain" )
        call xml_NewElement( xf, "Grid" )
        call xml_AddAttribute( xf, "GridType", "Collection" )
        call xml_AddAttribute( xf, "CollectionType", "Temporal" )
        it0 = 1
        if( restart ) it0 = restart_iter + 1
        ! Loop over time step
        do i = it0, iter + 1
            ! Grid name
            fname = ''
            fname = 'FWpoints.'
            call append_zero( fname,i - 1)
            str = '.xmf'
            call append_str( fname,str )
            call xml_NewElement( xf, "xi:include" )
            call xml_AddAttribute( xf, "href", trim( fname ) )
            call xml_AddAttribute( xf, "xpointer", "xpointer(//Xdmf/Domain/Grid)" )
            call xml_EndElement( xf, "xi:include" )
        enddo
        ! Footer
        call xml_EndElement( xf, "Grid" )
        call xml_EndElement( xf, "Domain" )
        call xml_EndElement( xf, "Xdmf" )
        call xml_Close( xf )

        return

    end subroutine FW_series
  ! ============================================================================

end module FWLK_out
