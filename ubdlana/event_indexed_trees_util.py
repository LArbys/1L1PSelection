import os,sys

STANDARD_EVENT_INDEXED_TREES_STR = """
    image2d_wire_tree
    chstatus_wire_tree
    image2d_ubspurn_plane0_tree
    image2d_ubspurn_plane1_tree
    image2d_ubspurn_plane2_tree
    image2d_thrumu_tree
    pgraph_test_tree
    pixel2d_test_ctor_tree
    pixel2d_test_img_tree
    pixel2d_test_super_ctor_tree
    pixel2d_test_super_img_tree
    image2d_test_inputimg_tree
    partroi_croimerge_clip_union_tree
    image2d_wcshower_tpc_tree
    image2d_wctrack_tpc_tree
    image2d_segment_tree
    image2d_instance_tree
    image2d_ancestor_tree
    image2d_larflow_tree
    partroi_segment_tree
    NuFilterTree
    MCTree
    EventVertexTree
    pgraph_inter_par_tree
    pixel2d_inter_par_pixel_tree
    pixel2d_inter_img_pixel_tree
    pixel2d_inter_int_pixel_tree
    larlite_id_tree
    daqheadertimeuboone_daq_tree
    gtruth_generator_tree
    mctruth_corsika_tree
    mctruth_generator_tree
    mcflux_generator_tree
    mcshower_mcreco_tree
    mctrack_mcreco_tree
    hit_gaushit_tree
    hit_portedThresholdhit_tree
    hit_dl_tree
    hit_dlraw_tree
    cluster_dl_tree
    cluster_dlraw_tree
    crthit_crthitcorr_tree
    crttrack_crttrack_tree
    ophit_ophitBeam_tree
    ophit_ophitBeamCalib_tree
    ophit_ophitCosmic_tree
    ophit_ophitCosmicCalib_tree
    opflash_opflashBeam_tree
    opflash_opflashCosmic_tree
    opflash_portedFlash_tree
    opflash_simpleFlashBeam_tree
    opflash_simpleFlashBeam::DLWCDeploy_tree
    opflash_simpleFlashCosmic_tree
    opflash_simpleFlashCosmic::DLWCDeploy_tree
    sps_portedSpacePointsThreshold_tree
    track_inter_track_tree
    track_trackReco_tree
    track_trackReco_sceadded_tree
    track_dl_tree
    shower_showerreco_tree
    shower_ssnetshowerreco_tree
    shower_dl_tree
    vertex_inter_vertex_tree
    vertex_trackReco_tree
    vertex_dl_tree
    vertex_dlraw_tree
    pfpart_dl_tree
    pfpart_dlraw_tree
    trigger_daq_tree
    trigger_triggersim_tree
    ass_inter_ass_tree
    ass_opflashBeam_tree
    ass_opflashCosmic_tree
    ass_portedFlash_tree
    ass_portedSpacePointsThreshold_tree
    ass_simpleFlashBeam_tree
    ass_simpleFlashBeam::DLWCDeploy_tree
    ass_simpleFlashCosmic::DLWCDeploy_tree  
    ass_trackReco_tree
    ass_trackReco_sceadded_tree
    ass_dl_tree
    ass_dlraw_tree
    ass_showerreco_tree
    swtrigger_swtrigger_tree
    larflowcluster_ssnetshowerreco_tree
    clustermask_mrcnn_masks_tree
    sparseimg_larflow_tree
    sparseimg_sparseuresnetout_tree
    shower_ssnetshowerrecov2ana_tree
    shower_ssnetshowerrecov2ana_sec_tree
    larflowcluster_ssnetshowerrecov2ana_tree
    track_dqdx_U_tree
    track_dqdx_V_tree
    track_dqdx_Y_tree
    ssnetshowerrecov2ana_anatree
    crtveto_tree
    ubdlana_id_tree
    mceventweight_eventweight4to4aFix_tree
"""

def make_event_indexed_trees_list():
    """ so gross """
    treelist = []
    info = STANDARD_EVENT_INDEXED_TREES_STR.split()
    for i in info:
        i = i.strip()
        if len(i)==0:
            continue
        treelist.append(i)
    return treelist

STANDARD_EVENT_INDEXED_TREES_LIST = make_event_indexed_trees_list()

def is_tree_event_indexed( treename ):
    strtreename = str(treename)

    if strtreename in STANDARD_EVENT_INDEXED_TREES_LIST:
        return True

    # other checks
    if len(strtreename.split("_"))>1 and strtreename.split("_")[0] in ["ass","ophit","opflash","mceventweight"]:
        return True

    return False
    