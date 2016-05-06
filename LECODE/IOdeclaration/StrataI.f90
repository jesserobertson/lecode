! ============================================================================
! Name        : StrataI.f90
! Author      : Tristan Salles
! Created on: Aug 14, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
!> \file StrataI.f90
!
! Description : StrataI is used to gather the information within the SPModel XmL input file to build
! the stratigraphic layers used in SPModel.
!
!<
! ============================================================================
module strata_parser

   use file_data
   use FoX_sax
   use mpidata
   use time_data
   use strata_data
   use fwalker_data
   use FoX_common
   use sediment_data

   implicit none

   public

   integer :: depF, sedNb, siNb, caNb, orNb, deform, claNb

   ! Output Directory
   logical, save :: in_OutDir = .false.
   ! Strata Parameters
   logical, save :: stratasection = .false.
   logical, save :: in_Strata = .false.
   logical, save :: in_StrataGrid = .false.
   logical, save :: in_GridXs = .false.
   logical, save :: in_GridYs = .false.
   logical, save :: in_GridXos = .false.
   logical, save :: in_GridYos = .false.
   logical, save :: in_GridSpaces = .false.
   character(len=128), save :: StrataGrid,grid_dxs,grid_dys
   character(len=128), save :: grid_xos,grid_yos,grid_spaces
   ! Sediment Parameters
   logical, save :: sedparamsection = .false.
   logical, save :: siparamsection = .false.
   logical, save :: caparamsection = .false.
   logical, save :: orparamsection = .false.
   logical, save :: in_Sediments = .false.
   logical, save :: in_Minslp = .false.
   logical, save :: in_Nbsiliciclastics = .false.
   logical, save :: in_Nborganics = .false.
   logical, save :: in_Nbcarbonates = .false.
   logical, save :: in_si = .false.
   logical, save :: in_ca = .false.
   logical, save :: in_or = .false.
   logical, save :: in_Matname = .false.
   logical, save :: in_Diameter = .false.
   logical, save :: in_Density = .false.
   logical, save :: in_Suspension = .false.
   logical, save :: in_Mar = .false.
   logical, save :: in_Aer = .false.
   character(len=128), save :: sinum, minslp, ornum, canum
   character(len=128), save :: SedDiameter,SedDensity,SedSuspension, Mar, Aer
   ! Stream Parameters
   logical, save :: streamparamsection = .false.
   logical, save :: clparamsection = .false.
   logical, save :: in_StreamClass = .false.
   logical, save :: in_NbClass = .false.
   logical, save :: in_cl = .false.
   logical, save :: in_BedSlope = .false.
   logical, save :: in_CDepth = .false.
   logical, save :: in_WDratio = .false.
   character(len=128), save :: clnum, BedSlope, CDepth, WDratio
   ! Mass movement Parameters
   logical, save :: masmovparamsection = .false.
   logical, save :: marineparamsection = .false.
   logical, save :: aerialparamsection = .false.
   logical, save :: in_MassMov = .false.
   logical, save :: in_MassFacc = .false.
   logical, save :: in_MassCurv = .false.
   logical, save :: in_AerialClass = .false.
   logical, save :: in_MarineClass = .false.
   logical, save :: in_uSA = .false.
   logical, save :: in_bUS = .false.
   logical, save :: in_aUS = .false.
   logical, save :: in_ssI = .false.
   logical, save :: in_cI = .false.
   character(len=128), save :: uSA, sI, cI, mfacc, mcurv
   ! Deposit Parameters
   logical, save :: depositsection = .false.
   logical, save :: in_Deposit = .false.
   logical, save :: in_DepoFile = .false.
   logical, save :: in_DepNb = .false.
   logical, save :: in_dep = .false.
   character(len=128), save :: Deposit, DepFile, DepNb, depN
   ! Check pointing Parameters
   logical, save :: checksection = .false.
   logical, save :: in_Check = .false.
   logical, save :: in_CheckFreq = .false.
!   logical, save :: in_udw = .false.
!   logical, save :: in_deformableStrata = .false.
!   logical, save :: in_fineSpace = .false.
   character(len=128), save :: checkfr !, udwr, fineSpace, deformableStrata
   ! Restart Parameters
   logical, save :: restartsection = .false.
   logical, save :: in_RestartFields = .false.
   logical, save :: in_RestartFile = .false.
   logical, save :: in_RestartIt = .false.
   logical, save :: in_RestartProc = .false.
   character(len=128), save :: rstit, rstproc

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

      ! Output Element
      if (name=='OutputDirectory') in_OutDir = .true.
      ! Strata Element
      if (name=='Strata_grid') in_Strata = .true.
      if (name=='Strata_grid') stratasection = .true.
      if(in_Strata) call SstrataElement_handler(name)
      ! Sediment Element
      if (name=='Materials') in_Sediments = .true.
      if(in_Sediments) sedparamsection = .true.
      if(in_Sediments) call SsedElement_handler(name)
      if (name=='si') in_si = .true.
      if (name=='carb') in_ca = .true.
      if (name=='org') in_or = .true.
      if(in_si) siparamsection = .true.
      if(in_si) call SsiElement_handler(name)
      if(in_ca) caparamsection = .true.
      if(in_ca) call ScaElement_handler(name)
      if(in_or) orparamsection = .true.
      if(in_or) call SorElement_handler(name)
      ! Stream Element
      if (name=='StreamClass') in_StreamClass = .true.
      if(in_StreamClass) streamparamsection = .true.
      if(in_StreamClass) call SstreamElement_handler(name)
      if (name=='cl') in_cl = .true.
      if(in_cl) clparamsection = .true.
      if(in_cl) call SclElement_handler(name)
      ! Mass movement Element
      if (name=='MassMovement') in_MassMov = .true.
      if(in_MassMov) masmovparamsection = .true.
      if(in_MassMov) call SmassmovElement_handler(name)
      if (name=='aerial') in_AerialClass = .true.
      if(in_AerialClass) aerialparamsection = .true.
      if(in_AerialClass) call SaerialElement_handler(name)
      if( in_bUS ) call SbaerialslopeElement_handler(name)
      if( in_aUS ) call SaaerialslopeElement_handler(name)
      if (name=='marine') in_MarineClass = .true.
      if(in_MarineClass) marineparamsection = .true.
      if(in_MarineClass) call SmarineElement_handler(name)
      if( in_bUS ) call SbmarineslopeElement_handler(name)
      if( in_aUS ) call SamarineslopeElement_handler(name)
      ! Deposit Element
      if (name=='Init_deposit') in_Deposit = .true.
      if(in_Deposit) depositsection = .true.
      if(in_Deposit) call SdepositElement_handler(name)
      if (name=='dep') in_dep = .true.
      if(in_dep) call SdepElement_handler(name)
      ! Restart Parameters
      if (name=='Restart') in_RestartFields = .true.
      if (name=='Restart') restartsection = .true.
      if(in_RestartFields) call SrestartElement_handler(name)
      ! Check pointing Element
      if (name=='CheckPointing') in_Check = .true.
      if(in_Check) checksection = .true.
      if(in_Check) call ScheckElement_handler(name)

   end subroutine startElement_handler
   !============================================================================
   subroutine endElement_handler(namespaceURI, localname, name)

      character(len=*), intent(in) :: namespaceURI
      character(len=*), intent(in) :: localname
      character(len=*), intent(in) :: name

      ! Output Element
      if (name=='OutputDirectory') in_OutDir = .false.
      ! Mesh Element
      call EstrataElement_handler(name)
      ! Sediment Element
      call EsedElement_handler(name)
      call EsiElement_handler(name)
      call EcaElement_handler(name)
      call EorElement_handler(name)
      ! Stream Element
      call EstreamElement_handler(name)
      call EclElement_handler(name)
      ! Mass Movement Element
      call EmassmovElement_handler(name)
      call EaerialElement_handler(name)
      call EbaerialslopeElement_handler( name )
      call EaaerialslopeElement_handler( name )
      call EmarineElement_handler(name)
      call EbmarineslopeElement_handler(name)
      call EamarineslopeElement_handler(name)
      ! Deposit Element
      call EdepositElement_handler(name)
      call EdepElement_handler(name)
      ! Restart Parameters
      call ErestartElement_handler(name)
      ! Check pointing Element
      call EcheckElement_handler(name)

   end subroutine endElement_handler
   ! ============================================================================
   subroutine characters_handler(chars)

      character(len=*), intent(in) :: chars

      ! Get Output Dircetory Name
      if (in_OutDir) then
         outdir=''
         outdir = chars
      endif
      ! Strata Element
      if(in_Strata) call strata_characters_handler(chars)
      ! Get Deposit File Name
      if (in_Deposit)   call deposit_characters_handler(chars)
      if (in_dep)   call dep_characters_handler(chars)
      ! Get Sediment Element
      if(in_Sediments) call sed_characters_handler(chars)
      if( in_si .or. in_ca .or. in_or )then
        if( .not. allocated( material_name ) ) allocate( material_name( totgrn ) )
        if( .not. allocated( sediment ) ) allocate(sediment(totgrn))
      endif
      if(in_si .and. silgrn > 0 ) call si_characters_handler(chars)
      if(in_ca .and. carbgrn > 0 ) call ca_characters_handler(chars)
      if(in_or .and. orggrn > 0 ) call or_characters_handler(chars)
      ! Get Stream Element
      if(in_StreamClass) call stream_characters_handler(chars)
      if(in_cl) call cl_characters_handler(chars)
      ! Get Mass Movement Element
      if(in_MassMov)then
          massmove = .true.
          if( in_MassFacc ) call massfacc_characters_handler(chars)
          if( in_MassCurv ) call masscurv_characters_handler(chars)
          if(in_AerialClass) call aerial_characters_handler(chars)
          if(in_MarineClass) call marine_characters_handler(chars)
      endif
      ! Restart Parameters
      if(in_RestartFields) call restart_characters_handler(chars)
      ! Get Check Pointing Element
      if (in_Check)   call check_characters_handler(chars)

   end subroutine characters_handler
   ! ============================================================================
   subroutine SstrataElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='StrataGrid') in_StrataGrid = .true.
      if (name=='GridX') in_GridXs = .true.
      if (name=='GridY') in_GridYs = .true.
      if (name=='GridXo') in_GridXos = .true.
      if (name=='GridYo') in_GridYos = .true.
      if (name=='GridSpace') in_GridSpaces = .true.

   end subroutine SstrataElement_handler
   ! ============================================================================
   subroutine SdepositElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='LayersNb') in_DepNb = .true.

   end subroutine SdepositElement_handler
   ! ============================================================================
   subroutine ScheckElement_handler(name)
      character(len=*), intent(in) :: name

      if (name=='Frequency') in_CheckFreq = .true.
!      if (name=='underworldCoupling') in_udw = .true.
!      if (name=='deformableStrata') in_deformableStrata = .true.
!      if (name=='fineSpace')        in_fineSpace = .true.

   end subroutine ScheckElement_handler
   ! ============================================================================
   subroutine SrestartElement_handler(name)
      character(len=*), intent(in) :: name

      if (name=='RestartFolder') in_RestartFile = .true.
      if (name=='RestartIt') in_RestartIt = .true.
      if (name=='ProcNb') in_RestartProc = .true.

   end subroutine SrestartElement_handler
   ! ============================================================================
   subroutine SdepElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='FileName') in_DepoFile = .true.

   end subroutine SdepElement_handler
   ! ============================================================================
   subroutine SsedElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='NbSiliciclastics') in_Nbsiliciclastics = .true.
      if (name=='NbOrganics') in_Nborganics = .true.
      if (name=='NbCarbonates') in_Nbcarbonates = .true.
      if (name=='MinSlope') in_Minslp = .true.

   end subroutine SsedElement_handler
   ! ============================================================================
   subroutine SstreamElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='NbClass') in_NbClass = .true.

   end subroutine SstreamElement_handler
   ! ============================================================================
   subroutine SmassmovElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='massFlowAccu') in_MassFacc = .true.
      if (name=='massCurvature') in_MassCurv = .true.
      if (name=='Aerial') in_AerialClass = .true.
      if (name=='Marine') in_MarineClass = .true.

   end subroutine SmassmovElement_handler
   ! ============================================================================
   subroutine SaerialElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='upSlopeArea') in_uSA = .true.
      if (name=='belowUpSlope') in_bUS = .true.
      if (name=='aboveUpSlope') in_aUS = .true.

   end subroutine SaerialElement_handler
   ! ============================================================================
   subroutine SaaerialslopeElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='steepnessIndex') in_ssI = .true.
      if (name=='concavityIndex') in_cI = .true.

   end subroutine SaaerialslopeElement_handler
   ! ============================================================================
   subroutine SbaerialslopeElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='steepnessIndex') in_ssI = .true.
      if (name=='concavityIndex') in_cI = .true.

   end subroutine SbaerialslopeElement_handler
   ! ============================================================================
   subroutine SamarineslopeElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='steepnessIndex') in_ssI = .true.
      if (name=='concavityIndex') in_cI = .true.

   end subroutine SamarineslopeElement_handler
   ! ============================================================================
   subroutine SbmarineslopeElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='steepnessIndex') in_ssI = .true.
      if (name=='concavityIndex') in_cI = .true.

   end subroutine SbmarineslopeElement_handler
   ! ============================================================================
   subroutine SmarineElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='upSlopeArea') in_uSA = .true.
      if (name=='belowUpSlope') in_bUS = .true.
      if (name=='aboveUpSlope') in_aUS = .true.

   end subroutine SmarineElement_handler
   ! ============================================================================
   subroutine SclElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='BedSlope') in_BedSlope = .true.
      if (name=='ChannelDepth') in_CDepth = .true.
      if (name=='WDRatio') in_WDratio = .true.

   end subroutine SclElement_handler
   ! ============================================================================
   subroutine SsiElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='MaterialName') in_Matname = .true.
      if (name=='Diameter') in_Diameter = .true.
      if (name=='Density') in_Density = .true.
      if (name=='Suspension') in_Suspension = .true.
      if (name=='SlopeMarine') in_Mar = .true.
      if (name=='SlopeAerial') in_Aer = .true.

   end subroutine SsiElement_handler
   ! ============================================================================
   subroutine ScaElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='MaterialName') in_Matname = .true.
      if (name=='Diameter') in_Diameter = .true.
      if (name=='Density') in_Density = .true.
      if (name=='SlopeMarine') in_Mar = .true.
      if (name=='SlopeAerial') in_Aer = .true.

   end subroutine ScaElement_handler
   ! ============================================================================
   subroutine SorElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='MaterialName') in_Matname = .true.
      if (name=='Diameter') in_Diameter = .true.
      if (name=='Density') in_Density = .true.
      if (name=='SlopeMarine') in_Mar = .true.
      if (name=='SlopeAerial') in_Aer = .true.

   end subroutine SorElement_handler
   ! ============================================================================
   subroutine EstrataElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='Strata_grid') in_Strata = .false.
      if (name=='StrataGrid') in_StrataGrid = .false.
      if (name=='GridX') in_GridXs = .false.
      if (name=='GridY') in_GridYs = .false.
      if (name=='GridXo') in_GridXos = .false.
      if (name=='GridYo') in_GridYos = .false.
      if (name=='GridSpace') in_GridSpaces = .false.

   end subroutine EstrataElement_handler
   ! ============================================================================
   subroutine EsedElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='NbSiliciclastics') in_Nbsiliciclastics = .false.
      if (name=='NbOrganics') in_Nborganics = .false.
      if (name=='NbCarbonates') in_Nbcarbonates = .false.
      if (name=='MinSlope') in_Minslp = .false.

   end subroutine EsedElement_handler
   ! ============================================================================
   subroutine EdepositElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='LayersNb') in_DepNb = .false.

   end subroutine EdepositElement_handler
   ! ============================================================================
   subroutine EcheckElement_handler(name)
      character(len=*), intent(in) :: name

      if (name=='Frequency') in_CheckFreq = .false.
!      if (name=='underworldCoupling') in_udw = .false.
!      if (name=='deformableStrata') in_deformableStrata = .false.
!      if (name=='fineSpace') in_fineSpace = .false.

   end subroutine EcheckElement_handler
   ! ============================================================================
   subroutine ErestartElement_handler(name)
      character(len=*), intent(in) :: name

      if (name=='RestartFolder') in_RestartFile = .false.
      if (name=='RestartIt') in_RestartIt = .false.
      if (name=='ProcNb') in_RestartProc = .false.

   end subroutine ErestartElement_handler
   ! ============================================================================
   subroutine EdepElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='FileName') in_DepoFile = .false.

   end subroutine EdepElement_handler
   ! ============================================================================
   subroutine EsiElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='MaterialName') in_Matname = .false.
      if (name=='Diameter') in_Diameter = .false.
      if (name=='Density') in_Density = .false.
      if (name=='Suspension') in_Suspension = .false.
      if (name=='SlopeMarine') in_Mar = .false.
      if (name=='SlopeAerial')then
        in_Aer = .false.
        in_si = .false.
      endif

   end subroutine EsiElement_handler
   ! ============================================================================
   subroutine EcaElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='MaterialName') in_Matname = .false.
      if (name=='Diameter') in_Diameter = .false.
      if (name=='Density') in_Density = .false.
      if (name=='SlopeMarine') in_Mar = .false.
      if (name=='SlopeAerial')then
        in_Aer = .false.
        in_ca = .false.
      endif

   end subroutine EcaElement_handler
   ! ============================================================================
   subroutine EorElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='MaterialName') in_Matname = .false.
      if (name=='Diameter') in_Diameter = .false.
      if (name=='Density') in_Density = .false.
      if (name=='SlopeMarine') in_Mar = .false.
      if (name=='SlopeAerial')then
        in_Aer = .false.
        in_or = .false.
      endif

   end subroutine EorElement_handler
   ! ============================================================================
   subroutine EstreamElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='NbClass') in_NbClass = .false.

   end subroutine EstreamElement_handler
   ! ============================================================================
   subroutine EclElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='BedSlope') in_BedSlope = .false.
      if (name=='ChannelDepth') in_CDepth = .false.
      if (name=='WDRatio') in_WDratio = .false.

   end subroutine EclElement_handler
   ! ============================================================================
   subroutine EmassmovElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='massFlowAccu') in_MassFacc = .false.
      if (name=='massCurvature') in_MassCurv = .false.
      if (name=='Aerial') in_AerialClass = .false.
      if (name=='Marine') in_MarineClass = .false.

   end subroutine EmassmovElement_handler
   ! ============================================================================
   subroutine EaerialElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='upSlopeArea') in_uSA = .false.
      if (name=='belowUpSlope') in_bUS = .false.
      if (name=='aboveUpSlope') in_aUS = .false.

   end subroutine EaerialElement_handler
   ! ============================================================================
   subroutine EaaerialslopeElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='steepnessIndex') in_ssI = .false.
      if (name=='concavityIndex') in_cI = .false.

   end subroutine EaaerialslopeElement_handler
   ! ============================================================================
   subroutine EbaerialslopeElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='steepnessIndex') in_ssI = .false.
      if (name=='concavityIndex') in_cI = .false.

   end subroutine EbaerialslopeElement_handler
   ! ============================================================================
   subroutine EamarineslopeElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='steepnessIndex') in_ssI = .false.
      if (name=='concavityIndex') in_cI = .false.

   end subroutine EamarineslopeElement_handler
   ! ============================================================================
   subroutine EbmarineslopeElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='steepnessIndex') in_ssI = .false.
      if (name=='concavityIndex') in_cI = .false.

   end subroutine EbmarineslopeElement_handler
   ! ============================================================================
   subroutine EmarineElement_handler(name)

      character(len=*), intent(in) :: name

      if (name=='upSlopeArea') in_uSA = .false.
      if (name=='belowUpSlope') in_bUS = .false.
      if (name=='aboveUpSlope') in_aUS = .false.

   end subroutine EmarineElement_handler
   ! ============================================================================
   subroutine strata_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if( in_StrataGrid)then
         fstrata = chars
      elseif (in_GridXs) then
         grid_dxs = chars
         call rts(grid_dxs,strat_X)
      elseif (in_GridYs) then
         grid_dys = chars
         call rts(grid_dys,strat_Y)
      elseif (in_GridXos) then
         grid_xos = chars
         call rts(grid_xos,strat_xo)
      elseif (in_GridYos) then
         grid_yos = chars
         call rts(grid_yos,strat_yo)
      elseif (in_GridSpaces) then
         grid_spaces = chars
         call rts(grid_spaces,strat_dx)
      endif

   end subroutine strata_characters_handler
   ! ============================================================================
   subroutine deposit_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if(in_DepNb)then
         depN = chars
         call rts(depN,InitDep)
         allocate( fdep( InitDep ) )
      endif

   end subroutine deposit_characters_handler
   ! ============================================================================
   subroutine dep_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if(in_DepoFile)then
         depF = depF + 1
         if( depF <= InitDep ) fdep( depF ) = chars
      endif

   end subroutine dep_characters_handler
   ! ============================================================================
   subroutine check_characters_handler(chars)
      character(len=*), intent(in) :: chars

      if (in_CheckFreq) then
         checkfr = chars
         checkpointing = .true.
         call rts(checkfr,checkfreq)
!      elseif(in_udw)then
!         udwr = chars
!         call rts(udwr,udw_coupling)
!      elseif(in_deformableStrata)then
!         deformableStrata = chars
!         call rts(deformableStrata,deform)
!         if( deform == 1 ) deformlayer = .true.
!      elseif(in_fineSpace)then
!         fineSpace = chars
!         call rts(fineSpace,fine_dx)
      endif

   end subroutine check_characters_handler
   ! ============================================================================
   subroutine restart_characters_handler(chars)
      character(len=*), intent(in) :: chars
      if (in_RestartFile) then
         restartfolder = chars
         restart = .true.
      elseif(in_RestartIt) then
         rstit = chars
         call rts(rstit, restart_iter )
      elseif(in_RestartProc) then
         rstproc = chars
         call rts(rstproc, restart_proc )
      endif

   end subroutine restart_characters_handler
   !============================================================================
   subroutine sed_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if(in_Nbsiliciclastics)then
         sinum = chars
         call rts(sinum, silgrn )
         totgrn = totgrn + silgrn
      endif
      if(in_Nbcarbonates)then
         canum = chars
         call rts(canum, carbgrn )
         totgrn = totgrn + carbgrn
      endif
      if(in_Nborganics)then
         ornum = chars
         call rts(ornum, orggrn )
         totgrn = totgrn + orggrn
      endif
      if(in_Minslp)then
         minslp = chars
         call rts(minslp, minimum_slp )
      endif

   end subroutine sed_characters_handler
   !============================================================================
   subroutine stream_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if(in_NbClass)then
         clnum = chars
         call rts(clnum, stcl_nb )
         allocate( streamclass( stcl_nb ) )
      endif

   end subroutine stream_characters_handler
   !============================================================================
   subroutine cl_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_BedSlope) then
         claNb = claNb + 1
         BedSlope = chars
         call rts(BedSlope,streamclass(claNb)%slope)
      elseif (in_WDratio) then
         WDratio = chars
         call rts(WDratio,streamclass(claNb)%wd_ratio)
      elseif (in_CDepth) then
         CDepth = chars
         call rts(CDepth,streamclass(claNb)%cdepth)
      endif

   end subroutine cl_characters_handler
   !============================================================================
   subroutine massfacc_characters_handler(chars)

      character(len=*), intent(in) :: chars

      mfacc = chars
      call rts(mfacc,mass_acc)

   end subroutine massfacc_characters_handler
   !============================================================================
   subroutine masscurv_characters_handler(chars)

      character(len=*), intent(in) :: chars

      mcurv = chars
      call rts(mcurv,mass_curv)

   end subroutine masscurv_characters_handler
   !============================================================================
   subroutine aerial_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_uSA) then
         uSA = chars
         call rts(uSA,mass_mov(1)%upslope)
      elseif (in_bUS .and. in_ssI ) then
         sI = chars
         call rts(sI,mass_mov(1)%steepness(1))
      elseif (in_bUS .and. in_cI) then
         cI = chars
         call rts(cI,mass_mov(1)%concavity(1))
      elseif (in_aUS .and. in_ssI ) then
         sI = chars
         call rts(sI,mass_mov(1)%steepness(2))
      elseif (in_aUS .and. in_cI) then
         cI = chars
         call rts(cI,mass_mov(1)%concavity(2))
      endif

   end subroutine aerial_characters_handler
   !============================================================================
   subroutine marine_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_uSA) then
         uSA = chars
         call rts(uSA,mass_mov(2)%upslope)
      elseif (in_bUS .and. in_ssI ) then
         sI = chars
         call rts(sI,mass_mov(2)%steepness(1))
      elseif (in_bUS .and. in_cI) then
         cI = chars
         call rts(cI,mass_mov(2)%concavity(1))
      elseif (in_aUS .and. in_ssI ) then
         sI = chars
         call rts(sI,mass_mov(2)%steepness(2))
      elseif (in_aUS .and. in_cI) then
         cI = chars
         call rts(cI,mass_mov(2)%concavity(2))
      endif

   end subroutine marine_characters_handler
   !============================================================================
   subroutine si_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_Matname) then
         siNb = siNb + 1
         sedNb = siNb
         material_name(sedNb) = chars
      elseif (in_Diameter) then
         SedDiameter = chars
         call rts(SedDiameter,sediment(sedNb)%diameter)
      elseif(in_Density) then
         SedDensity = chars
         call rts(SedDensity,sediment(sedNb)%density)
      elseif(in_Suspension)then
         SedSuspension = chars
         call rts(SedSuspension,sediment(sedNb)%transport)
      elseif(in_Mar)then
         Mar = chars
         call rts(Mar,sediment(sedNb)%slp_marine)
      elseif(in_Aer)then
         Aer = chars
         call rts(Aer,sediment(sedNb)%slp_aerial)
      endif

   end subroutine si_characters_handler
   !============================================================================
   subroutine ca_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_Matname) then
         caNb = caNb + 1
         sedNb = silgrn + caNb
         material_name(sedNb) = chars
      elseif (in_Diameter) then
         SedDiameter = chars
         call rts(SedDiameter,sediment(sedNb)%diameter)
      elseif(in_Density) then
         SedDensity = chars
         call rts(SedDensity,sediment(sedNb)%density)
      elseif(in_Suspension)then
         SedSuspension = chars
         call rts(SedSuspension,sediment(sedNb)%transport)
      elseif(in_Mar)then
         Mar = chars
         call rts(Mar,sediment(sedNb)%slp_marine)
      elseif(in_Aer)then
         Aer = chars
         call rts(Aer,sediment(sedNb)%slp_aerial)
      endif

   end subroutine ca_characters_handler
   !============================================================================
   subroutine or_characters_handler(chars)

      character(len=*), intent(in) :: chars

      if (in_Matname) then
         orNb = orNb + 1
         sedNb = silgrn + carbgrn + orNb
         material_name(sedNb) = chars
      elseif (in_Diameter) then
         SedDiameter = chars
         call rts(SedDiameter,sediment(sedNb)%diameter)
      elseif(in_Density) then
         SedDensity = chars
         call rts(SedDensity,sediment(sedNb)%density)
      elseif(in_Suspension)then
         SedSuspension = chars
         call rts(SedSuspension,sediment(sedNb)%transport)
      elseif(in_Mar)then
         Mar = chars
         call rts(Mar,sediment(sedNb)%slp_marine)
      elseif(in_Aer)then
         Aer = chars
         call rts(Aer,sediment(sedNb)%slp_aerial)
      endif

   end subroutine or_characters_handler
   ! ============================================================================
   subroutine endDocument_handler

   end subroutine endDocument_handler
   ! ============================================================================
   !> Subroutine xml_Strata()
   !! reads input that defines the SP Model experiment (XmL)
   !<
   ! ============================================================================
   subroutine xml_Strata

      type(xml_t) :: xf
      integer :: ios

      depF = 0
      sedNb = 0
      siNb = 0
      caNb = 0
      orNb = 0
      claNb = 0
      InitDep = 0
      silgrn = 0
      orggrn = 0
      carbgrn = 0
      totgrn = 0
      stcl_nb = 0
      massmove = .false.
      checkpointing = .false.

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

   end subroutine xml_Strata
   ! ============================================================================
   !> Subroutine broadcast_strata()
   !! Broadcast initial parameters which define the stratigraphy.
   !<
   ! ============================================================================
   subroutine broadcast_strata

      ! Parameters Declaration
      integer :: i, j

      ! Broadcast restart parameters
      call mpi_bcast( restart,1,lgc_type,0,SPModel_comm_world,ierr )
      if( restart )then
         call mpi_bcast( restartfolder,128,mpi_character,0,SPModel_comm_world,ierr )
         call mpi_bcast( restart_iter,1,int_type,0,SPModel_comm_world,ierr )
         call mpi_bcast( restart_proc,1,int_type,0,SPModel_comm_world,ierr )
      endif
      ! Broadcast checkpointing parameters
      call mpi_bcast( checkpointing,1,lgc_type,0,SPModel_comm_world,ierr )
      if( checkpointing )then
         call mpi_bcast( checkfreq,1,int_type,0,SPModel_comm_world,ierr )
         udw_coupling = 0
         deformlayer = .false.
         fine_dx = 100.0_8
!         call mpi_bcast( udw_coupling,1,int_type,0,SPModel_comm_world,ierr )
!         call mpi_bcast( deformlayer,1,lgc_type,0,SPModel_comm_world,ierr )
!         call mpi_bcast( fine_dx,1,dbl_type,0,SPModel_comm_world,ierr )
      endif
      if( udw_coupling == 1 ) udw = .true.
      ! Broadcast deposit layers
      call mpi_bcast( InitDep,1,int_type,0,SPModel_comm_world,ierr )
      if( InitDep < 1 )then
        if( iam == 0 )print*,'An initial deposit is required to run a simulation'
        stop
      endif
      ! Broadcast file names
      call mpi_bcast( outdir,128,mpi_character,0,SPModel_comm_world,ierr )
      call mpi_bcast( outdir1,128,mpi_character,0,SPModel_comm_world,ierr )
      call mpi_bcast( outdir2,128,mpi_character,0,SPModel_comm_world,ierr )
      call mpi_bcast( fstrata,128,mpi_character,0,SPModel_comm_world,ierr )
      call mpi_bcast( outputs,128,mpi_character,0,SPModel_comm_world,ierr )
      call mpi_bcast( runfiles,128,mpi_character,0,SPModel_comm_world,ierr )
      if( iam /= 0 ) allocate( fdep( InitDep ) )
      do i = 1, InitDep
         call mpi_bcast( fdep( i ),128,mpi_character,0,SPModel_comm_world,ierr )
      enddo
      ! Grid space parameters
      call mpi_bcast( strat_X,1,int_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( strat_Y,1,int_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( strat_xo,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( strat_yo,1,dbl_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( strat_dx,1,dbl_type,0,SPModel_comm_world,ierr )
      ! Broadcast sediment parameters
      call mpi_bcast( silgrn,1,int_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( carbgrn,1,int_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( orggrn,1,int_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( totgrn,1,int_type,0,SPModel_comm_world,ierr )
      call mpi_bcast( minimum_slp,1,dbl_type,0,SPModel_comm_world,ierr )
      if( totgrn < 1 )then
        if( iam == 0 ) print*,'At least one sediment class  is required to run a simulation.'
        stop
      endif
      if( totgrn /= silgrn + carbgrn + orggrn )then
        if( iam == 0 ) print*,'Something went wrong when defining material types.'
        stop
      endif
      if( iam /= 0 ) allocate( sediment(totgrn) )
      if( iam /= 0 ) allocate( material_name(totgrn) )
      do i = 1, totgrn
         sediment(i)%transport = 1
         call mpi_bcast( material_name(i),128,mpi_character,0,SPModel_comm_world,ierr )
         call mpi_bcast( sediment(i)%diameter,1,dbl_type,0,SPModel_comm_world,ierr )
         call mpi_bcast( sediment(i)%density,1,dbl_type,0,SPModel_comm_world,ierr )
!         call mpi_bcast( sediment(i)%transport,1,int_type,0,SPModel_comm_world,ierr )
         call mpi_bcast( sediment(i)%slp_marine,1,dbl_type,0,SPModel_comm_world,ierr )
         call mpi_bcast( sediment(i)%slp_aerial,1,dbl_type,0,SPModel_comm_world,ierr )
      enddo
      ! Broadcast stream classification parameters
      call mpi_bcast( stcl_nb,1,int_type,0,SPModel_comm_world,ierr )
      if( stcl_nb < 2 )then
        if( iam == 0 )print*,'Stream classification should contained at least 2 parameters'
        stop
      endif
      if( iam /= 0 ) allocate( streamclass(stcl_nb) )
      do i = 1, stcl_nb
         call mpi_bcast( streamclass(i)%slope,1,dbl_type,0,SPModel_comm_world,ierr )
         call mpi_bcast( streamclass(i)%cdepth,1,dbl_type,0,SPModel_comm_world,ierr )
         call mpi_bcast( streamclass(i)%wd_ratio,1,dbl_type,0,SPModel_comm_world,ierr )
      enddo
      ! Broadcast mass movement parameters
      call mpi_bcast( massmove,1,lgc_type,0,SPModel_comm_world,ierr )
      if( massmove )then
          call mpi_bcast( mass_curv,1,dbl_type,0,SPModel_comm_world,ierr )
          call mpi_bcast( mass_acc,1,int_type,0,SPModel_comm_world,ierr )
          do i = 1, 2
             call mpi_bcast( mass_mov(i)%upslope,1,dbl_type,0,SPModel_comm_world,ierr )
             do j = 1, 2
                 call mpi_bcast( mass_mov(i)%steepness(j),1,dbl_type,0,SPModel_comm_world,ierr )
                 call mpi_bcast( mass_mov(i)%concavity(j),1,dbl_type,0,SPModel_comm_world,ierr )
            enddo
          enddo
      endif
      call mpi_barrier( SPModel_comm_world,ierr )

      return

   end subroutine broadcast_strata
   ! ============================================================================

end module strata_parser
