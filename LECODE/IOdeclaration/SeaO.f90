! ============================================================================
! Name        : SeaO.f90
! Author      : Tristan Salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file SeaO.f90
!
! Description : SeaO is used to generate Hdf5 output of the water surface during SPModel simulation.
!
!<
! ============================================================================
module SL_out

   use file_data
   use mpidata
   use time_data
   use FoX_wxml
   use strata_data
   use forces_data

   implicit none

   public

   integer ::  totnodes, totelems

contains

   ! ============================================================================
   !> Subroutine SL_xmf
   !! Subroutine SL_xmf generates the water surface XdmF file for the considered time simulation.
   !<
   ! ============================================================================
   subroutine SL_xmf( iter )

      ! Parameters Declaration
      type(xmlf_t) :: xf

      integer :: iter, i, j
      real( tkind ) :: minx,miny,maxx,maxy
      real( tkind ) :: pt(4, 3 )

      character(len=128) :: str, file

      ! Create nodes arrays
      minx = strat_xo
      miny = strat_yo
      maxx = strat_xo + strat_dx * ( strat_X - 1 )
      maxy = strat_yo + strat_dx * ( strat_Y - 1 )
      pt( 1, 1 ) = minx
      pt( 2, 1 ) = minx
      pt( 3, 1 ) = maxx
      pt( 4, 1 ) = maxx
      pt( 1, 2 ) = miny
      pt( 1, 3 ) = gsea%actual_sea
      pt( 2, 2 ) = maxy
      pt( 2, 3 ) = gsea%actual_sea
      pt( 3, 2 ) = miny
      pt( 3, 3 ) = gsea%actual_sea
      pt( 4, 2 ) = maxy
      pt( 4, 3 ) = gsea%actual_sea
      if( udw )then
         pt( 1, 3 ) = miny
         pt( 1, 2 ) = gsea%actual_sea
         pt( 2, 3 ) = maxy
         pt( 2, 2 ) = gsea%actual_sea
         pt( 3, 3 ) = miny
         pt( 3, 2 ) = gsea%actual_sea
         pt( 4, 3 ) = maxy
         pt( 4, 2 ) = gsea%actual_sea
      endif
      file = ''
      file = 'Water'
      call noblnk( file )
      str = '.'
      call append_str( file,str )
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
      call xml_AddAttribute( xf, "GridType", "Uniform" )
      call xml_NewElement( xf, "Time" )
      call xml_AddAttribute( xf, "Type", "Single" )
      call xml_AddAttribute( xf, "Value", next_display-time_display )
      call xml_EndElement( xf, "Time" )
      ! Topology
      call xml_NewElement( xf, "Topology" )
      call xml_AddAttribute( xf, "TopologyType", "Quadrilateral" )
      call xml_AddAttribute( xf, "Dimensions", "1" )
      call xml_NewElement( xf, "DataItem" )
      call xml_AddAttribute( xf, "NumberType", "Int" )
      call xml_AddAttribute( xf, "Dimensions", "1 4" )
      call xml_AddAttribute( xf, "Precision", "4" )
      call xml_AddAttribute( xf, "Format", "XML" )
      call xml_AddNewline( xf )
      call xml_AddCharacters(xf, "0 1 3 2")
      call xml_AddNewline( xf )
      call xml_EndElement( xf, "DataItem" )
      call xml_EndElement( xf, "Topology" )
      ! Geometry
      call xml_NewElement( xf, "Geometry" )
      call xml_AddAttribute( xf, "GeometryType", "XYZ" )
      call xml_NewElement( xf, "DataItem" )
      call xml_AddAttribute( xf, "NumberType", "Float" )
      call xml_AddAttribute( xf, "Dimensions", "4 3" )
      call xml_AddAttribute( xf, "Precision", "4" )
      call xml_AddAttribute( xf, "Format", "XML" )
      call xml_AddNewline( xf )
      do i = 1, 4
         str = ' '
         do j = 1, 3
            call append_nbreal( str, pt( i, j ) )
         enddo
         call xml_AddCharacters( xf, trim( str ) )
         call xml_AddNewline( xf )
      enddo
      call xml_EndElement( xf, "DataItem" )
      call xml_EndElement( xf, "Geometry" )
      ! Footer
      call xml_EndElement( xf, "Grid" )
      call xml_EndElement( xf, "Domain" )
      call xml_EndElement( xf, "Xdmf" )
      call xml_Close( xf )

      return

   end subroutine SL_xmf
   ! ============================================================================
   !> Subroutine SL_series
   !! Subroutine SL_series generates the XmL file for sea surface.
   !<
   ! ============================================================================
   subroutine SL_series( iter )

      ! Parameters Declaration
      type(xmlf_t) :: xf

      integer :: i, iter, it0
      character(len=128) :: filename, str, fname

      filename = 'Water_series.xdmf'
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
         fname = 'Water'
         call noblnk( fname )
         str = '.'
         call append_str( fname,str )
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

   end subroutine SL_series
  ! ============================================================================

end module SL_out
