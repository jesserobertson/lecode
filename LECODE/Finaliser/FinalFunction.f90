! ============================================================================
! Name        : FinalFunction.f90
! Author      : tristan salles
! Created on: Aug 16, 2012
! Copyright (C) 2012 CSIRO
! ============================================================================
module final

    use file_data
    use mpidata
    use time_data
    use TIN_out
    use flux_data
    use fwalker_data
    use error_data
    use strata_out
    use strata_data
    use forces_data
    use param_data
    use mod_diffuse
    use TIN_function
    use diffuse_fct
    use ocean_data
    use mod_compaction
    use mod_displacements
    use mod_filldepression
    use mod_edtransfer
    use mod_functionsfw

    implicit none

    public

contains

    ! ============================================================================
    !> Subroutine spm_final
    !! Subroutine spm_final finalises the SPModel code execution.
    !<
    ! ============================================================================
    subroutine spm_final

        ! Deallocation
        if( allocated( wdem ) ) deallocate( wdem )
        if( allocated( facc ) ) deallocate( facc )
        if( allocated( halo_proc ) ) deallocate( halo_proc )
        if( allocated( halo_lid ) ) deallocate( halo_lid, halo_rcv2, halo_rcv1 )
        if( allocated( halo_gid ) ) deallocate( halo_gid, halo1_elev, halo2_elev )

        deallocate( locv_NbLay, locv_gid, locv_coord, locv_vrate, locv_vdisp )
        deallocate( gv_halo, gv_ShareID, gv_SharenID, gv_ngbID )
        deallocate( gv_vdisp, gv_coord, gv_orientation )

        if( allocated( strat_layID ) ) deallocate( strat_layID )
        if( allocated( strat_hardness ) ) deallocate( strat_hardness )
        if( allocated( strat_zelev ) ) deallocate( strat_zelev )
        if( allocated( strat_porosity ) ) deallocate( strat_porosity )
        if( allocated( strat_thick ) ) deallocate( strat_thick )
        if( allocated( strat_fac ) ) deallocate( strat_fac )
        if( allocated( strat_sedh ) ) deallocate( strat_sedh )
        if( allocated( hemipelagic ) ) deallocate( hemipelagic )

        deallocate( disp, eromax )
        deallocate( pick_refine, elev_record, proc_elev , uplift )
        deallocate( material_name )
        deallocate( tinf_boundary, tinf_nghb, tinf_pid, fids )
        deallocate( tinf_nids, tinf_centroid, tinf_length )
        deallocate( tinv_pid, ngbID, tinv_ngb, tinv_coord, tinv_boundary )
        deallocate( gf_pid, gf_points, gf_XYcoord )
        deallocate( locf_gid, locf_points, locf_XYcoord )

        if( allocated( fdep ) ) deallocate( fdep )
        if( allocated( rain_areaNb ) )then
            deallocate( rain_areaNb, rain_tend, rain_tstart )
            deallocate( rain_xmin, rain_xmax, rain_ymin, rain_ymax )
            deallocate( rain_h, rain_min, rain_max )
        endif

        if( allocated( frain_node ) )then
            deallocate( frain_node, frain_xposition, frain_yposition, frain_qfl )
        endif

        if( allocated( fa_refid ) )then
          deallocate( fa_refid, fa_nbPtAI, fa_inbound, fa_type, fa_count, fa_width )
          deallocate( fa_xpos, fa_ypos, fa_zpos, fa_xslp, fa_yslp, fa_h, fa_xvel )
          deallocate( fa_yvel, fa_xvelO, fa_yvelO, fa_volume, fa_discharge, fa_density )
          deallocate( fa_pgrad, fa_Cfric, fa_accx, fa_accy, fa_pzfe, fa_pcoeff, fa_pdist )
          deallocate( fa_pvelx, fa_pvely, fa_faceID, fa_ptsAIfw, fa_sedcharge )
        endif

        if( allocated( frainmap ) ) deallocate( frainmap )

        if( num_src > 0 .and. allocated( fws_num ) )then
            deallocate( fws_num, fws_type, fws_tstrt, fws_tend, fws_perc, fws_xposition )
            deallocate( fws_yposition, fws_xrange, fws_yrange, fws_srch, fws_xvel, fws_yvel )
            deallocate( fws_qfl, fws_volume, fws_sedconc, fws_sedperc, fws_sedcharge )
        endif
        if( allocated( sealvl ) ) deallocate( sealvl )
        if( allocated( fdisp ) ) deallocate( fdisp )
        if( allocated( gdisp_time ) ) deallocate( gdisp_time )

        if( allocated( difp ) )then
            deallocate( difp, difo, topmax, cdif, depo, dstart )
            deallocate( top_sedh, top_sedprev )
            deallocate( tdifp, bdifp, bh, th )
            deallocate( difps1, difps2 )
            deallocate( diff_ID, diff_ref, diff_coord, diff_ngb )
        endif

        if( ocean_plug )then
            if( allocated( hindcast ) ) deallocate( hindcast )
            if( allocated( ocean_current ) ) deallocate( ocean_current )
            if( allocated( ocean_wave ) ) deallocate( ocean_wave )
            if( allocated( total_load ) ) deallocate( total_load )
            if( allocated( soulsby_total_load ) ) deallocate( soulsby_total_load )
            if( allocated( transportX ) ) deallocate( transportX )
            if( allocated( transportY ) ) deallocate( transportY )
        endif

        if( allocated( slp ) ) deallocate( slp )
        if( allocated( curv ) ) deallocate( curv )

        ! Deallocate erosion deposition arrays
        if( allocated( ed_refid ) )then
            deallocate( ed_refid, ed_procid, ed_bound, ed_type )
            deallocate( ed_AIpts, ed_AInb, ed_zpos, ed_h, ed_vel )
            deallocate( ed_vol, ed_dens, ed_sed )
        endif

        ! Free MPI types
        call mpi_type_free( mpi_concfw, ierr )
        call mpi_type_free( mpi_flowwalker, ierr )
        call mpi_type_free( mpi_flowtransfer, ierr )

        return

    end subroutine spm_final
    ! ============================================================================

end module final
