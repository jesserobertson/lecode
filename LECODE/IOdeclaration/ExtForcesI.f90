! ============================================================================
! Name        : ExtForcesI.f90
! Author      : Tristan Salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file ExtForcesI.f90
!
! Description : ExtForcesI is used to gather the information within the SPModel XmL input file to build
! the external forces used in SPModel.
!
!<
! ============================================================================
module Forces_parser

   use file_data
   use mpidata
   use FoX_sax
   use flux_data
   use forces_data
   use time_data
   use FoX_common
   use strata_data

   implicit none

   public

   integer :: dispn
   ! Water section
   logical, save :: in_OutSL = .false.
   character(len=128), save :: outsl
   ! Sea Parameters
   logical, save :: seasection = .false.
   logical, save :: in_Sea = .false.
   logical, save :: in_SeaFile = .false.
   character(len=128), save :: SeaFile
   ! Hemipelagic Parameters
   logical, save :: hemisection = .false.
   logical, save :: in_Hemi = .false.
   logical, save :: in_HemiFile = .false.
   character(len=128), save :: HemiFile
   ! Displacements Parameters
   logical, save :: vdispsection = .false.
   logical, save :: in_DispF = .false.
   logical, save :: in_DispNb = .false.
   logical, save :: in_DispFile = .false.
   logical, save :: in_DispET = .false.
   logical, save :: in_DispST = .false.
   logical, save :: in_disp = .false.
   character(len=128), save :: DispF,DispNb,DispFile,DispET,DispST
   ! Underworld Parameters
   logical, save :: udwsection = .false.
   logical, save :: in_SynchT = .false.
   logical, save :: in_SynchF = .false.
   logical, save :: in_UDW = .false.
   character(len=128), save :: synchF,synchT
   ! Porosity Parameters
   logical, save :: in_Porosity = .false.
   logical, save :: in_EffP = .false.
   logical, save :: in_pFile = .false.
   character(len=128), save :: effp, pfile
   ! Seismic Parameters
   logical, save :: in_Seismic = .false.
   logical, save :: in_Xmin = .false.
   logical, save :: in_Xmax = .false.
   logical, save :: in_Ymin = .false.
   logical, save :: in_Ymax = .false.
   character(len=128), save :: sxmin, sxmax, symin, symax
   ! RMS Parameters
   logical, save :: in_RMS = .false.
   logical, save :: in_RMSXmin = .false.
   logical, save :: in_RMSXmax = .false.
   logical, save :: in_RMSYmin = .false.
   logical, save :: in_RMSYmax = .false.
   character(len=128), save :: rmsxmin, rmsxmax, rmsymin, rmsymax

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

      ! Output Water
      if (name=='OutputWaterLevel') in_OutSL = .true.
      ! Underworld Element
      if (name=='UnderworldPlugin') in_UDW = .true.
      if (in_UDW) udw_plug = .true.
      if (in_UDW) call SudwElement_handler(name)
      ! Sea Element
      if (name=='SeaLevelFluctuations') in_Sea = .true.
      if (name=='SeaLevelFluctuations') seasection = .true.
      if(in_Sea) call SseaElement_handler(name)
      ! Hemipelagic Element
      if (name=='HemipelagicRates') in_Hemi = .true.
      if (name=='HemipelagicRates') hemisection = .true.
      if(in_Hemi) call SHemiElement_handler(name)
      ! Vertical Displacements Parameters
      if (name=='Displacement') in_DispF = .true.
      if (in_DispF) vdispsection = .true.
      if(in_DispF) call SvdispElement_handler(name)
      if (name=='disp')then
         in_disp = .true.
         dispn = dispn + 1
      endif
      if(in_disp) call SdispElement_handler(name)
      ! Porosity Element
      if (name=='Porosity') in_Porosity = .true.
      if(in_Porosity) call SporosityElement_handler(name)
      ! Seismic Element
      if (name=='SeismicLine') in_Seismic = .true.
      if(in_Seismic) seismic_plug = .true.
      if(in_Seismic) call SseismicElement_handler(name)
      ! Seismic Element
      if (name=='RMSModel') in_RMS = .true.
      if(in_RMS) rms_plug = .true.
      if(in_RMS) call SrmsElement_handler(name)

   end subroutine startElement_handler
   ! ============================================================================
   subroutine endElement_handler(namespaceURI, localname, name)

      character(len=*), intent(in) :: namespaceURI
      character(len=*), intent(in) :: localname
      character(len=*), intent(in) :: name

      ! Output Water
      if (name=='OutputWaterLevel') in_OutSL = .false.
      ! Sea Element
      call EseaElement_handler(name)
      ! Hemipelagic Element
      call EHemiElement_handler(name)
      ! Vertical Displacements Parameters
      call EvdispElement_handler(name)
      call EdispElement_handler(name)
      ! Porosity Element
      call EporosityElement_handler(name)
      ! UDW Element
      call EudwElement_handler(name)
      ! Seismic Element
      call EseismicElement_handler(name)
      ! RMS Element
      call ErmsElement_handler(name)

   end subroutine endElement_handler
   ! ============================================================================
   subroutine characters_handler(chars)

      character(len=*), intent(in) :: chars
      integer :: intsl

      ! Water SL Element
      if(in_OutSL)then
         outsl = ''
         outsl = chars
         call rts(outsl, intsl)
         if( intsl == 0 ) gsea%output = .false.
      endif
      ! Mesh Element
      if(in_Sea) call sea_characters_handler(chars)
      ! Hemipelagic Element
      if(in_Hemi) call hemi_characters_handler(chars)
      ! Vertical Displacements Parameters
      if(in_DispF) call vdisp_characters_handler(chars)
      if(in_disp) call disp_characters_handler(chars)
      ! Get Porosity Element
      if(in_Porosity) call porosity_characters_handler(chars)
      ! Get UDW parameter
      if(in_UDW) call udw_characters_handler(chars)
      ! Get Seismic Element
      if(in_Seismic) call seismic_characters_handler(chars)
      ! Get RMS Element
      if(in_RMS) call rms_characters_handler(chars)

   end subroutine characters_handler
   ! ============================================================================
   subroutine SseaElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='SeaLevelFile') in_SeaFile = .true.

   end subroutine SseaElement_handler
   ! ============================================================================
   subroutine SHemiElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='HemipelagicFile') in_HemiFile = .true.

   end subroutine SHemiElement_handler
   ! ============================================================================
   subroutine SvdispElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='nbDispInterval') in_DispNb = .true.

   end subroutine SvdispElement_handler
   ! ============================================================================
   subroutine SseismicElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='Xmin') in_Xmin = .true.
      if (name=='Xmax') in_Xmax = .true.
      if (name=='Ymin') in_Ymin = .true.
      if (name=='Ymax') in_Ymax = .true.

   end subroutine SseismicElement_handler
   ! ============================================================================
   subroutine SrmsElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='Xmin') in_RMSXmin = .true.
      if (name=='Xmax') in_RMSXmax = .true.
      if (name=='Ymin') in_RMSYmin = .true.
      if (name=='Ymax') in_RMSYmax = .true.

   end subroutine SrmsElement_handler
   ! ============================================================================
   subroutine SporosityElement_handler(name)
      character(len=*), intent(in) :: name

      if (name=='EffPressureNb') in_EffP = .true.
      if (name=='PorosityFile') in_pFile = .true.

   end subroutine SporosityElement_handler
   ! ============================================================================
   subroutine SdispElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='dispFile') in_DispFile = .true.
      if (name=='startDispTime') in_DispST = .true.
      if (name=='endDispTime') in_DispET = .true.

   end subroutine SdispElement_handler
   ! ============================================================================
   subroutine SudwElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='UsyncFolder') in_synchF = .true.
      if (name=='UsyncTime') in_synchT = .true.

   end subroutine SudwElement_handler
   ! ============================================================================
   subroutine EseaElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='SeaLevelFluctuations') in_Sea = .false.
      if (name=='SeaLevelFile') in_SeaFile = .false.

   end subroutine EseaElement_handler
   ! ============================================================================
   subroutine EHemiElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='HemipelagicRates') in_Hemi = .false.
      if (name=='HemipelagicFile') in_HemiFile = .false.

   end subroutine EHemiElement_handler
   ! ============================================================================
   subroutine EvdispElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='nbDispInterval') in_DispNb = .false.

   end subroutine EvdispElement_handler
   ! ============================================================================
   subroutine EseismicElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='Xmin') in_Xmin = .false.
      if (name=='Xmax') in_Xmax = .false.
      if (name=='Ymin') in_Ymin = .false.
      if (name=='Ymax') in_Ymax = .false.

   end subroutine EseismicElement_handler
   ! ============================================================================
   subroutine ErmsElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='Xmin') in_RMSXmin = .false.
      if (name=='Xmax') in_RMSXmax = .false.
      if (name=='Ymin') in_RMSYmin = .false.
      if (name=='Ymax') in_RMSYmax = .false.

   end subroutine ErmsElement_handler
   ! ============================================================================
   subroutine EdispElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='startDispTime') in_DispST = .false.
      if (name=='endDispTime') in_DispET = .false.
      if (name=='dispFile') in_DispFile = .false.

   end subroutine EdispElement_handler
   ! ============================================================================
   subroutine EudwElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='UsyncFolder') in_synchF = .false.
      if (name=='UsyncTime') in_synchT = .false.

   end subroutine EudwElement_handler
   ! ============================================================================
   subroutine EporosityElement_handler(name)
      character(len=*), intent(in) :: name

      if (name=='EffPressureNb') in_EffP = .false.
      if (name=='PorosityFile') in_pFile = .false.

   end subroutine EporosityElement_handler
   ! ============================================================================
   subroutine sea_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if( in_SeaFile)then
         fsea = chars
         gsea%sealevel = .true.
      endif

   end subroutine sea_characters_handler
   ! ============================================================================
   subroutine hemi_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if( in_HemiFile)then
         fhemi = chars
         hemi_flag = .true.
      endif

   end subroutine hemi_characters_handler
   ! ============================================================================
   subroutine vdisp_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if(in_DispNb ) then
         DispNb = chars
         call rts(DispNb, gdisp%event )
         allocate( fdisp( gdisp%event ) )
         allocate( gdisp_time( gdisp%event, 2 ) )
      endif

   end subroutine vdisp_characters_handler
   ! ============================================================================
   subroutine seismic_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_Xmin) then
         sxmin = chars
         call rts(sxmin,seismicX_min)
      elseif (in_Xmax) then
         sxmax = chars
         call rts(sxmax,seismicX_max)
      elseif (in_Ymin) then
         symin = chars
         call rts(symin,seismicY_min)
      elseif (in_Ymax) then
         symax = chars
         call rts(symax,seismicY_max)
      endif

   end subroutine seismic_characters_handler
   ! ============================================================================
   subroutine rms_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_RMSXmin) then
         rmsxmin = chars
         call rts(rmsxmin,RMSX_min)
      elseif (in_RMSXmax) then
         rmsxmax = chars
         call rts(rmsxmax,RMSX_max)
      elseif (in_RMSYmin) then
         rmsymin = chars
         call rts(rmsymin,RMSY_min)
      elseif (in_RMSYmax) then
         rmsymax = chars
         call rts(rmsymax,RMSY_max)
      endif

   end subroutine rms_characters_handler
   ! ============================================================================
   subroutine udw_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if(in_synchF ) then
         outdir3 = chars
      endif
      if(in_synchT ) then
         synchT = chars
         call rts(synchT, udw_time )
      endif

   end subroutine udw_characters_handler
   ! ============================================================================
   subroutine disp_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_DispFile ) fdisp( dispn ) = chars
      if(in_DispST ) then
         DispST = chars
         call rts(DispST, gdisp_time( dispn, 1 ) )
      endif
      if(in_DispET ) then
         DispET = chars
         call rts(DispET, gdisp_time( dispn, 2 ) )
      endif

   end subroutine disp_characters_handler
   ! ============================================================================
   subroutine porosity_characters_handler(chars)
      character(len=*), intent(in) :: chars

      if (in_EffP) then
         gporo%compaction = .true.
         effp = chars
         call rts(effp,gporo%ePnb)
      elseif (in_pFile) then
         fporosity = chars
      endif

   end subroutine porosity_characters_handler
   ! ============================================================================
   subroutine endDocument_handler

   end subroutine endDocument_handler
   ! ============================================================================
   !> Subroutine xml_Forces()
   !! reads input that defines the SP Model experiment (XmL)
   !<
   ! ============================================================================
   subroutine xml_Forces

      type(xml_t) :: xf
      integer :: ios

      dispn = 0
      gdisp%event = 0
      gsea%output = .true.
      gsea%sealevel = .false.
      hemi_flag = .false.
      udw_plug = .false.
      seismic_plug = .false.
      rms_plug = .false.

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

   end subroutine xml_Forces
   ! ============================================================================
   !> Subroutine broadcast_forces()
   !! Broadcast initial parameters which define the external forces.
   !<
   ! ============================================================================
   subroutine broadcast_forces

      integer :: k

      character(len=128) :: file

      call mpi_bcast( gporo%compaction,1,lgc_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( gsea%sealevel,1,lgc_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( gsea%output,1,lgc_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( gdisp%event,1,int_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( hemi_flag,1,lgc_type,0,SPModel_comm_world,ierr )
      if( iam /= 0 ) allocate( fdisp( gdisp%event ) )
      if( iam /= 0 ) allocate( gdisp_time( gdisp%event, 2 ) )
      ! Broadcast sea level file name
      if( gsea%sealevel ) call mpi_bcast( fsea,128,mpi_character,0,SPModel_comm_world,ierr )
      ! Broadcast hemipelagic file name
      if( hemi_flag ) call mpi_bcast( fhemi,128,mpi_character,0,SPModel_comm_world,ierr )
      if( gporo%compaction )then
         ! Broadcast porosity file name
         call mpi_bcast( fporosity,128,mpi_character,0,SPModel_comm_world,ierr )
         call mpi_bcast( gporo%ePnb,1,int_type,0,SPModel_comm_world,ierr )
      endif
      if( gdisp%event > 0 )then
         ! Broadcast vertical displacement file name
         do k = 1, gdisp%event
            call mpi_bcast( fdisp( k ),128,mpi_character,0,SPModel_comm_world,ierr )
            call mpi_bcast( gdisp_time( k, 1 ),1,dbl_type,0,SPModel_comm_world,ierr )
            call mpi_bcast( gdisp_time( k, 2 ),1,dbl_type,0,SPModel_comm_world,ierr )
         enddo
      endif
      ! Broadcast UDW plugin
      call mpi_bcast( udw_plug,1,lgc_type,0,SPModel_comm_world,ierr )
      if( udw_plug )then
          call mpi_bcast( outdir3,128,mpi_character,0,SPModel_comm_world,ierr )
          call mpi_bcast( udw_time,1,dbl_type,0,SPModel_comm_world,ierr )
          fudw = 'topsurface.vtk'
          fudisp = 'uw_output.ascii'
          maestro = 'maestro'
          file = ''
          file = fudw
          call addpath3( file )
          fudw = file
          file = ''
          file = maestro
          call addpath3( file )
          maestro = file
          file = ''
          file = fudisp
          call addpath3( file )
          fudisp = file
          ! Assign displacements
          gdisp%event = int( ( time_end - time_start ) / udw_time )
          if( allocated( fdisp ) ) deallocate( fdisp )
          allocate( fdisp( gdisp%event ) )
          if( allocated( gdisp_time ) ) deallocate( gdisp_time )
          allocate( gdisp_time( gdisp%event, 2 ) )
          ! Broadcast vertical displacement parameter
          do k = 1, gdisp%event
            fdisp( k ) = fudisp
            gdisp_time( k, 1 ) = time_start + udw_time * ( k - 1 )
            gdisp_time( k, 2 ) = gdisp_time( k, 1 ) +  udw_time
          enddo
      endif
      ! Broadcast Seismic plugin
      call mpi_bcast( seismic_plug,1,lgc_type,0,SPModel_comm_world,ierr )
      if( seismic_plug )then
          call mpi_bcast( seismicX_min,1,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( seismicX_max,1,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( seismicY_min,1,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( seismicY_max,1,dbl_type,0,SPModel_comm_world,ierr )
      endif
      ! Broadcast RMS plugin
      call mpi_bcast( rms_plug,1,lgc_type,0,SPModel_comm_world,ierr )
      if( rms_plug )then
          call mpi_bcast( RMSX_min,1,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( RMSX_max,1,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( RMSY_min,1,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( RMSY_max,1,dbl_type,0,SPModel_comm_world,ierr )
      endif
      gsea%actual_sea = 0.0_8

      return

   end subroutine broadcast_forces
  ! ============================================================================

end module Forces_parser
