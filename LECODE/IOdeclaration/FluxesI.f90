! ============================================================================
! Name        : FluxesI.f90
! Author      : Tristan Salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file FluxesI.f90
!
! Description : FluxesI is used to gather the information within the SPModel XmL input file to build
! the external influxes variables used in SPModel.
!
!<
! ============================================================================
module Flux_parser

   use file_data
   use mpidata
   use FoX_sax
   use flux_data
   use error_data
   use fwalker_data
   use flux_data
   use strata_data
   use FoX_common

   implicit none

   public

   integer :: srcNb, rainn

   ! Landscape Parameters Region
   logical, save :: rainsection = .false.
   logical, save :: in_RainFields = .false.
   logical, save :: in_RainRegion = .false.
   logical, save :: in_RainElev = .false.
   logical, save :: in_RainAlt = .false.
   logical, save :: in_RainFacc = .false.
   logical, save :: in_RainFile = .false.
   logical, save :: in_RainComb = .false.
   logical, save :: in_RainSQ = .false.
   character(len=128), save :: RainRegion,RainFields,RainElev,RainFile,RainAlt, RainFacc, RainComb, RainSQ
   ! Landscape Parameters Grid
   logical, save :: rainmapsection = .false.
   logical, save :: in_RainGrid = .false.
   logical, save :: in_RainNb = .false.
   logical, save :: in_RainMFile = .false.
   logical, save :: in_RainET = .false.
   logical, save :: in_RainST = .false.
   logical, save :: in_rain = .false.
   character(len=128), save :: RainGrid,RainNb,RainMFile,RainET,RainST
   ! Sources Parameters
   logical, save :: sourcesparamsection = .false.
   logical, save :: sourcesection = .false.
   logical, save :: in_Sources = .false.
   logical, save :: in_SourcesNb = .false.
   logical, save :: in_maxSourcesE = .false.
   logical, save :: in_source = .false.
   logical, save :: in_t1 = .false.
   logical, save :: in_t2 = .false.
   logical, save :: in_x = .false.
   logical, save :: in_y = .false.
   logical, save :: in_xrange = .false.
   logical, save :: in_yrange = .false.
   logical, save :: in_flowHeight = .false.
   logical, save :: in_seaPerc = .false.
   logical, save :: in_Vx = .false.
   logical, save :: in_Vy = .false.
   logical, save :: in_Q = .false.
   logical, save :: in_Qs = .false.
   logical, save :: in_sedConcentration = .false.
   logical, save :: in_flowType = .false.
   character(len=128), save :: srcnum, srcel, srct1, srct2, srcx, srcy, srcxr, srcyr, srcvx, srcvy, srcfh
   character(len=128), save :: srcQ, srcQs, srcsC, srcT, seaPerc

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

      ! Landscape Parameters Region
      if (name=='RainfallRegion') in_RainFields = .true.
      if (name=='RainfallRegion') rainsection = .true.
      if(in_RainFields) call SrainElement_handler(name)
      ! Landscape Parameters Grid
      if (name=='RainfallGrid') in_RainGrid = .true.
      if (in_RainGrid) rainmapsection = .true.
      if(in_RainGrid) call SrainmapElement_handler(name)
      if (name=='rain')then
         in_rain = .true.
         rainn = rainn + 1
      endif
      if(in_rain) call SrainfieldElement_handler(name)
      ! Sources Element
      if (name=='Sources') in_Sources = .true.
      if(in_Sources) sourcesparamsection = .true.
      if(in_Sources) call SsourcesElement_handler(name)
      if (name=='sourceVal') in_source = .true.
      if(in_source) sourcesection = .true.
      if(in_source) call SsourceElement_handler(name)

   end subroutine startElement_handler
   ! ============================================================================
   subroutine endElement_handler(namespaceURI, localname, name)

      character(len=*), intent(in) :: namespaceURI
      character(len=*), intent(in) :: localname
      character(len=*), intent(in) :: name

      ! Landscape Parameters Region
      call ErainElement_handler(name)
      ! Landscape Parameters Grid
      call ErainmapElement_handler(name)
      call ErainfieldElement_handler(name)
      ! Sources Element
      call EsourcesElement_handler(name)
      call EsourceElement_handler(name)

   end subroutine endElement_handler
   ! ============================================================================
   subroutine characters_handler(chars)

      character(len=*), intent(in) :: chars

      ! Landscape Parameters Region
      if(in_RainFields) call rain_characters_handler(chars)
      ! Landscape Parameters Grid
      if(in_RainGrid) call rainmap_characters_handler(chars)
      if(in_rain) call rainfield_characters_handler(chars)
      ! Get Sources Element
      if(in_Sources)  call sources_characters_handler(chars)
      if(in_source)   call source_characters_handler(chars)

   end subroutine characters_handler
   ! ============================================================================
   subroutine SsourcesElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='totalSourcesNb') in_SourcesNb = .true.
      if (name=='maxDepth') in_maxSourcesE = .true.

   end subroutine SsourcesElement_handler
   ! ============================================================================
   subroutine SsourceElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='ts') in_t1 = .true.
      if (name=='te') in_t2 = .true.
      if (name=='x') in_x = .true.
      if (name=='y') in_y = .true.
      if (name=='xRange') in_xrange = .true.
      if (name=='yRange') in_yrange = .true.
      if (name=='streamDepth') in_flowHeight = .true.
      if (name=='recurencePerc') in_seaPerc = .true.
      if (name=='Vx') in_Vx = .true.
      if (name=='Vy') in_Vy = .true.
      if (name=='Q') in_Q = .true.
      if (name=='Qs') in_Qs = .true.
      if(name=='sedConcentration') in_sedConcentration = .true.
      if (name=='flowType') in_flowType = .true.

   end subroutine SsourceElement_handler
   ! ============================================================================
   subroutine SrainElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='rainfallCyclicity') in_RainFile = .true.
      if (name=='rainElevation') in_RainElev = .true.
      if (name=='rainMaxFW') in_RainAlt = .true.
      if (name=='rainFlowAccu') in_RainFacc = .true.
      if (name=='rainCombination') in_RainComb = .true.
      if (name=='rainStreamQ') in_RainSQ = .true.

   end subroutine SrainElement_handler
   ! ============================================================================
   subroutine SrainmapElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='nbRainfallInterval') in_RainNb = .true.
      if (name=='rainElevation') in_RainElev = .true.
      if (name=='rainMaxFW') in_RainAlt = .true.
      if (name=='rainFlowAccu') in_RainFacc = .true.
      if (name=='rainCombination') in_RainComb = .true.
      if (name=='rainStreamQ') in_RainSQ = .true.

   end subroutine SrainmapElement_handler
   ! ============================================================================
   subroutine SrainfieldElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='rainFile') in_RainMFile = .true.
      if (name=='startTime') in_RainST = .true.
      if (name=='endTime') in_RainET = .true.

   end subroutine SrainfieldElement_handler
   ! ============================================================================
   subroutine EsourcesElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='totalSourcesNb') in_SourcesNb = .false.
      if (name=='maxDepth') in_maxSourcesE = .false.

   end subroutine EsourcesElement_handler
   ! ============================================================================
   subroutine EsourceElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='ts') in_t1 = .false.
      if (name=='te') in_t2 = .false.
      if (name=='x')  in_x = .false.
      if (name=='y')  in_y = .false.
      if (name=='xRange')  in_xrange = .false.
      if (name=='yRange')  in_yrange = .false.
      if (name=='streamDepth')  in_flowHeight = .false.
      if (name=='recurencePerc')  in_seaPerc = .false.
      if (name=='Vx') in_Vx = .false.
      if (name=='Vy') in_Vy = .false.
      if (name=='Q') in_Q = .false.
      if (name=='Qs') in_Qs = .false.
      if (name=='sedConcentration') in_sedConcentration = .false.
      if (name=='flowType') in_flowType = .false.

   end subroutine EsourceElement_handler
   ! ============================================================================
   subroutine ErainElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='rainElevation') in_RainElev = .false.
      if (name=='rainfallCyclicity') in_RainFile = .false.
      if (name=='rainMaxFW') in_RainAlt = .false.
      if (name=='rainFlowAccu') in_RainFacc = .false.
      if (name=='rainCombination') in_RainComb = .false.
      if (name=='rainStreamQ') in_RainSQ = .false.

   end subroutine ErainElement_handler
   ! ============================================================================
   subroutine ErainmapElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='rainElevation') in_RainElev = .false.
      if (name=='nbRainfallInterval') in_RainNb = .false.
      if (name=='rainMaxFW') in_RainAlt = .false.
      if (name=='rainFlowAccu') in_RainFacc = .false.
      if (name=='rainCombination') in_RainComb = .false.
      if (name=='rainStreamQ') in_RainSQ = .false.

   end subroutine ErainmapElement_handler
   ! ============================================================================
   subroutine ErainfieldElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='startTime') in_RainST = .false.
      if (name=='endTime') in_RainET = .false.
      if (name=='rainFile') in_RainMFile = .false.

   end subroutine ErainfieldElement_handler
   ! ============================================================================
   subroutine sources_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_SourcesNb) then
         srcnum = chars
         call rts(srcnum, num_src )
         allocate( fws_num(num_src) )
         allocate( fws_type(num_src) )
         allocate( fws_tstrt(num_src) )
         allocate( fws_tend(num_src) )
         allocate( fws_perc(num_src) )
         allocate( fws_xposition(num_src) )
         allocate( fws_yposition(num_src) )
         allocate( fws_xrange(num_src) )
         allocate( fws_yrange(num_src) )
         allocate( fws_srch(num_src) )
         allocate( fws_xvel(num_src) )
         allocate( fws_yvel(num_src) )
         allocate( fws_qfl(num_src) )
         allocate( fws_volume(num_src) )
         allocate( fws_sedconc(num_src) )
         allocate( fws_sedperc(num_src, silgrn ) )
         allocate( fws_sedcharge(num_src, silgrn ) )
         fws_sedperc=0.0_8
         fws_sedcharge=0.0_8
         fws_xrange=0.0_8
         fws_yrange=0.0_8
     elseif( in_maxSourcesE )then
         srcel = chars
         call rts(srcel, max_elev_src )
     endif

   end subroutine sources_characters_handler
   ! ============================================================================
   subroutine source_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_t1) then
         srcNb = srcNb + 1
         if( srcNb > num_src ) attempt = SRCNBS
         srct1 = chars
         call rts(srct1,fws_tstrt(srcNb))
      elseif(in_t2) then
         srct2 = chars
         call rts(srct2,fws_tend(srcNb))
      elseif(in_x)then
         srcx = chars
         call rts(srcx,fws_xposition(srcNb))
         fws_xposition(srcNb) = fws_xposition(srcNb)
      elseif(in_y) then
         srcy = chars
         call rts(srcy,fws_yposition(srcNb))
         fws_yposition(srcNb) = fws_yposition(srcNb)
      elseif(in_xrange)then
         srcxr = chars
         call rts(srcxr,fws_xrange(srcNb))
         fws_xrange(srcNb) = fws_xrange(srcNb) * 0.5_8
      elseif(in_yrange) then
         srcyr = chars
         call rts(srcyr,fws_yrange(srcNb))
         fws_yrange(srcNb) = fws_yrange(srcNb) * 0.5_8
      elseif(in_flowHeight) then
         srcfh = chars
         call rts(srcfh,fws_srch(srcNb))
      elseif(in_seaPerc) then
         seaPerc = chars
         call rts(seaPerc,fws_perc(srcNb))
      elseif(in_Vx)then
         srcvx = chars
         call rts(srcvx,fws_xvel(srcNb))
      elseif(in_Vy) then
         srcvy = chars
         call rts(srcvy,fws_yvel(srcNb))
      elseif(in_Q) then
         srcQ = chars
         call rts(srcQ,fws_qfl(srcNb))
      elseif(in_Qs)then
         srcQs = chars
         call rts(srcQs,fws_sedconc(srcNb))
      elseif(in_sedConcentration) then
         srcsC = chars
         call rts(srcsC,fws_sedperc(srcNb,1:silgrn))
      elseif(in_flowType)then
         srcT = chars
         call rts(srcT,fws_type(srcNb))
      endif

   end subroutine source_characters_handler
   ! ============================================================================
   subroutine rain_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_RainFile) then
         rain_region = .true.
         geomorphology = .true.
         frain = chars
      elseif(in_RainElev) then
         RainElev= chars
         call rts(RainElev, rain_elev )
      elseif(in_RainAlt) then
         RainAlt= chars
         call rts(RainAlt, rain_nb )
      elseif(in_RainFacc) then
         RainFacc= chars
         call rts(RainFacc, rain_facc )
      elseif(in_RainComb) then
         RainComb= chars
         call rts(RainComb, combine )
      elseif(in_RainSQ) then
         RainSQ= chars
         call rts(RainSQ, rainstreamQ )
      endif

   end subroutine rain_characters_handler
   ! ============================================================================
   subroutine rainmap_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_RainNb) then
         geomorphology = .true.
         rain_region = .false.
         RainNb= chars
         call rts(RainNb, rain_event )
         allocate( frainmap( rain_event ) )
         allocate( rain_areaNb( rain_event ) )
         allocate( rain_tend( rain_event ) )
         allocate( rain_tstart( rain_event ) )
         allocate( rain_xmin( rain_event, G_strat%nodes_nb ) )
         allocate( rain_xmax( rain_event, G_strat%nodes_nb ) )
         allocate( rain_ymin( rain_event, G_strat%nodes_nb ) )
         allocate( rain_ymax( rain_event, G_strat%nodes_nb ) )
         allocate( rain_h( rain_event, G_strat%nodes_nb ) )
         allocate( rain_min( rain_event, G_strat%nodes_nb ) )
         allocate( rain_max( rain_event, G_strat%nodes_nb ) )
     elseif(in_RainElev) then
         RainElev= chars
         call rts(RainElev, rain_elev )
     elseif(in_RainAlt) then
         RainAlt= chars
         call rts(RainAlt, rain_nb )
      elseif(in_RainFacc) then
         RainFacc= chars
         call rts(RainFacc, rain_facc )
      elseif(in_RainComb) then
         RainComb= chars
         call rts(RainComb, combine )
      elseif(in_RainSQ) then
         RainSQ= chars
         call rts(RainSQ, rainstreamQ )
      endif

   end subroutine rainmap_characters_handler
   ! ============================================================================
   subroutine rainfield_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_RainMFile) then
         frainmap( rainn ) = chars
      endif
      if(in_RainST ) then
         RainST = chars
         call rts(RainST, rain_tstart( rainn ) )
      endif
      if(in_RainET ) then
         RainET = chars
         call rts(RainET, rain_tend( rainn ) )
      endif

   end subroutine rainfield_characters_handler
   ! ============================================================================
   subroutine endDocument_handler

   end subroutine endDocument_handler
   ! ============================================================================
   !> Subroutine xml_Flux()
   !! reads input that defines the SP Model experiment (XmL)
   !<
   ! ============================================================================
   subroutine xml_Flux

      type(xml_t) :: xf
      integer :: ios

      srcNb = 0
      rainn = 0
      rain_event = 0
      max_elev_src = 1.0_8
      G_strat%nodes_nb = strat_X * strat_Y
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
      if( srcNb /= num_src )   attempt = SRCNBS

   end subroutine xml_Flux
   ! ============================================================================
   !> Subroutine broadcast_flux()
   !! Broadcast initial parameters which define the flux.
   !<
   ! ============================================================================
   subroutine broadcast_flux

      integer :: k

      call mpi_bcast( geomorphology,1,lgc_type,0,SPModel_comm_world,ierr )

      G_strat%nodes_nb = strat_X * strat_Y
      ! Broadcast geomorphology parameters
      if( geomorphology )then
         call mpi_bcast( rain_region,1,lgc_type,0,SPModel_comm_world,ierr )
         if( rain_region )then
            call mpi_bcast( frain,128,mpi_character,0,SPModel_comm_world,ierr )
         else
            call mpi_bcast( rain_event,1,int_type,0,SPModel_comm_world,ierr )
            if( iam /= 0 )then
                 allocate( frainmap( rain_event ) )
                 allocate( rain_areaNb( rain_event ) )
                 allocate( rain_tend( rain_event ) )
                 allocate( rain_tstart( rain_event ) )
                 allocate( rain_xmin( rain_event, G_strat%nodes_nb ) )
                 allocate( rain_xmax( rain_event, G_strat%nodes_nb ) )
                 allocate( rain_ymin( rain_event, G_strat%nodes_nb ) )
                 allocate( rain_ymax( rain_event, G_strat%nodes_nb ) )
                 allocate( rain_h( rain_event, G_strat%nodes_nb ) )
                 allocate( rain_min( rain_event, G_strat%nodes_nb ) )
                 allocate( rain_max( rain_event, G_strat%nodes_nb ) )
            endif
            do k = 1, rain_event
               call mpi_bcast( rain_tstart( k ),1,dbl_type,0,SPModel_comm_world,ierr )
               call mpi_bcast( rain_tend( k ),1,dbl_type,0,SPModel_comm_world,ierr )
               call mpi_bcast( frainmap( k ),128,mpi_character,0,SPModel_comm_world,ierr )
            enddo
         endif
         call mpi_bcast( rain_nb,1,int_type,0,SPModel_comm_world,ierr )
         call mpi_bcast( rain_facc,1,int_type,0,SPModel_comm_world,ierr )
         call mpi_bcast( combine,1,int_type,0,SPModel_comm_world,ierr )
         call mpi_bcast( rain_elev,1,dbl_type,0,SPModel_comm_world,ierr )
         call mpi_bcast( rainstreamQ,1,dbl_type,0,SPModel_comm_world,ierr )
         combine = combine + 1
      else
         combine = 1
      endif
      ! Broadcast source parameters
      call mpi_bcast( num_src,1,int_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( max_elev_src,1,dbl_type,0,SPModel_comm_world,ierr )
      if( iam /= 0 .and. num_src > 0 )then
         allocate( fws_num(num_src) )
         allocate( fws_type(num_src) )
         allocate( fws_tstrt(num_src) )
         allocate( fws_tend(num_src) )
         allocate( fws_perc(num_src) )
         allocate( fws_xposition(num_src) )
         allocate( fws_yposition(num_src) )
         allocate( fws_xrange(num_src) )
         allocate( fws_yrange(num_src) )
         allocate( fws_srch(num_src) )
         allocate( fws_xvel(num_src) )
         allocate( fws_yvel(num_src) )
         allocate( fws_qfl(num_src) )
         allocate( fws_volume(num_src) )
         allocate( fws_sedconc(num_src) )
         allocate( fws_sedperc(num_src, silgrn ) )
         allocate( fws_sedcharge(num_src, silgrn ) )
      endif
      if( num_src > 0 )then
          call mpi_bcast( fws_tstrt,num_src,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( fws_tend,num_src,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( fws_xposition,num_src,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( fws_yposition,num_src,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( fws_srch,num_src,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( fws_perc,num_src,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( fws_xvel,num_src,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( fws_yvel,num_src,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( fws_qfl,num_src,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( fws_sedconc,num_src,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( fws_type,num_src,int_type,0,SPModel_comm_world,ierr )
          do k = 1, silgrn
              call mpi_bcast( fws_sedperc( 1:num_src , k ),num_src,dbl_type,0,SPModel_comm_world,ierr )
          enddo
      endif

      ! Pts in the vicinity of a FW and subject to ero/dep
      if( strat_X > strat_Y )then
          if( 2.0_8 * ( strat_X - 1 ) * strat_dx / nproc <  transport%wdtmax )then
            transport%wdtmax = 2.0_8 * strat_dx * ( strat_X - 2 ) / nproc
            if( iam == 0 ) print*,'Note: the maximum width of streams has been reset to:',&
                2.0_8 * ( strat_X - 2 ) * strat_dx / nproc
          endif
          maxPtsFW = max( 4 * int( transport%wdtmax / strat_dx + 1 ), 1 )
      else
          if( 2.0_8 * ( strat_Y - 1 ) * strat_dx / nproc <  transport%wdtmax )then
            transport%wdtmax = 2.0_8 * strat_dx * ( strat_Y - 2 ) / nproc
            if( iam == 0 ) print*,'Note: the maximum width of streams has been reset to:',&
                2.0_8 * ( strat_Y - 2 ) * strat_dx / nproc
          endif
          maxPtsFW = max( 4 * int( transport%wdtmax / strat_dx + 1 ), 1 )
      endif

      ! Allocate flow walker type
      allocate( fa_refid( maxfw ) )
      allocate( fa_nbPtAI( maxfw ) )
      allocate( fa_inbound( maxfw ) )
      allocate( fa_type( maxfw ) )
      allocate( fa_count( maxfw ) )
      allocate( fa_xpos( maxfw ) )
      allocate( fa_ypos( maxfw ) )
      allocate( fa_zpos( maxfw ) )
      allocate( fa_xslp( maxfw ) )
      allocate( fa_yslp( maxfw ) )
      allocate( fa_width( maxfw ) )
      allocate( fa_h( maxfw ) )
      allocate( fa_xvel( maxfw ) )
      allocate( fa_yvel( maxfw ) )
      allocate( fa_xvelO( maxfw ) )
      allocate( fa_yvelO( maxfw ) )
      allocate( fa_volume( maxfw ) )
      allocate( fa_discharge( maxfw ) )
      allocate( fa_density( maxfw ) )
      allocate( fa_pgrad( maxfw ) )
      allocate( fa_Cfric( maxfw ) )
      allocate( fa_accx( maxfw ) )
      allocate( fa_accy( maxfw ) )
      allocate( fa_pzfe( maxfw ) )
      allocate( fa_pcoeff( maxfw ) )
      allocate( fa_pdist( maxfw ) )
      allocate( fa_pvelx( maxfw ) )
      allocate( fa_pvely( maxfw ) )
      allocate( fa_faceID( maxfw, numold ) )
      allocate( fa_ptsAIfw( maxfw, maxPtsFW ) )
      allocate( fa_sedcharge( maxfw, totgrn ) )

      if( rain_nb > G_strat%nodes_nb )then
        print*,'The number of maximum rain FWs is greater than the number of points in the grid'
        print*,'The value is reset to:',G_strat%nodes_nb
        rain_nb = G_strat%nodes_nb
      endif

      return

   end subroutine broadcast_flux
  ! ============================================================================

end module Flux_parser
