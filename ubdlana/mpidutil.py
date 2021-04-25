import os, sys, gc
import ROOT
from larcv import larcv
import numpy as np

from lib.config import config_loader
from lib_mpid_torch.rootdata_pid import ROOTData

#larcv.LArbysLoader()

import torch
import torch.nn as nn
from mpid_data import mpid_data
from mpid_net import mpid_net, mpid_func

p_type = {0:"eminus", 1:"gamma", 2:"muon", 3:"piminus",4:"proton"}
train_device = 'cuda' if torch.cuda.is_available() else 'cpu'


def image_np2torch(img_np, cfg, train_device):
    img_np = image_modify(img_np, cfg)
    img_np = torch.from_numpy(img_np.copy())
    img_np = img_np.clone().detach()
    img_torch = img_np.to(train_device).view((-1,1,512,512))
    return img_torch
    

def image_modify (img, cfg):
    img_arr = np.array(img.as_vector())
    img_arr = np.where(img_arr<cfg.adc_lo,         0,img_arr)
    img_arr = np.where(img_arr>cfg.adc_hi,cfg.adc_hi,img_arr)
    img_arr = img_arr.reshape(cfg.batch,img_arr.size).astype(np.float32)

    return img_arr

def make_mpid_anatree(tfile):
    tfile.cd()
    mpid_dir = tfile.mkdir("mpid")
    mpid_dir.cd()
    
    rd = ROOTData()
    tree  = ROOT.TTree("multipid_tree","UB MPID score tree")
    rd.init_tree(tree)
    rd.reset()
    return rd,tree

def load_mpid_model(CFG):
    cfg  = config_loader(CFG)
    assert cfg.batch == 1
    mpid = mpid_net.MPID()
    weight_file = ""
    plane=2
    exec("weight_file = cfg.weight_file_mpid_%i" % plane)
    weight_file = weight_file.replace("__WEIGHT_FOLDER__",os.environ["UBMPIDNET_WEIGHT_DIR"])
    print "MPID MODEL loaded with WEIGHT FILE: ",weight_file
    mpid.load_state_dict(torch.load(weight_file, map_location=train_device))
    mpid.eval()
    return mpid,cfg

def run_mpid_on_larcv_entry( cfg, mpid, iom, rd, outtree, return_result_dict=False ):
    """ run the network on current enetry in larcv iomanager 
    inputs
    ------
    cfg  [config class] created by load_mpid_model
    mpid [torch model]  created by load_mpid_model
    iom  [larcv::IOManager]
    
    modified
    --------
    rd [ROOTData]   created by make_mpid_anatree (use this one to make sure array addresses tied to tree branches)
    outtree [TTree] saves scores per vertex, index via (run,subrun,event,vtxid)
    """
    
    ev_pgr = iom.get_data(larcv.kProductPGraph,"inter_par")
    ev_par = iom.get_data(larcv.kProductPixel2D,"inter_par_pixel")
    ev_pix = iom.get_data(larcv.kProductPixel2D,"inter_img_pixel")
    ev_int = iom.get_data(larcv.kProductPixel2D,"inter_int_pixel")
        
    print '========================>>>>>>>>>>>>>>>>>>>>'
    print 'MPID running on run, subrun, event',ev_pix.run(),ev_pix.subrun(),ev_pix.event()

    rd.run[0]    = int(ev_pix.run())
    rd.subrun[0] = int(ev_pix.subrun())
    rd.event[0]  = int(ev_pix.event())
    rd.entry[0]  = int(iom.current_entry())

    rd.num_vertex[0] = int(ev_pgr.PGraphArray().size())

    print 'num of vertices, ',rd.num_vertex[0]
    print 'pgrapgh size, ',int(ev_pgr.PGraphArray().size())
    nfilled = 0

    result_dict = {}
        
    for ix,pgraph in enumerate(ev_pgr.PGraphArray()):
        print "@pgid=%d" % ix
        #if (ix != 2): continue
        rd.vtxid[0] = int(ix)
   
        rsev = (rd.run[0], rd.subrun[0], rd.event[0], rd.vtxid[0])
        result_dict[rsev] = {}

        pgr = ev_pgr.PGraphArray().at(ix)
        cluster_array = pgr.ClusterIndexArray()
        if cluster_array.size()==0:
            cindex_v = []
        else:
            cindex_v = np.array(pgr.ClusterIndexArray())
        
        pixel2d_par_vv = ev_par.Pixel2DClusterArray()
        pixel2d_pix_vv = ev_pix.Pixel2DClusterArray()
        pixel2d_int_vv = ev_int.Pixel2DClusterArray()
        
        for plane in xrange(3):
            print "@plane=%d" % plane

            if pixel2d_par_vv.size() != 0:
                rd.npar[plane] = 0
                pixel2d_par_v = pixel2d_par_vv.at(plane)
                for cidx in cindex_v:
                    print "@cindex_v=",cidx," of ",pixel2d_par_v.size()
                    try:
                        pixel2d_par = pixel2d_par_v.at(cidx)
                        if pixel2d_par.size()>0:
                            rd.npar[plane] += 1;
                    except:
                        print "error opening particle pixel2d"
				                  

            pixel2d_pix_v = pixel2d_pix_vv.at(plane)
            pixel2d_pix = pixel2d_pix_v.at(ix)

            pixel2d_int_v = pixel2d_int_vv.at(plane)
            pixel2d_int = pixel2d_int_v.at(ix)
            
            # nothing on this plane
            #if pixel2d_pix.empty() == True: continue

            rd.inferred[0] = 1
                
            img_pix = larcv.cluster_to_image2d(pixel2d_pix,cfg.xdim,cfg.ydim)
            img_int = larcv.cluster_to_image2d(pixel2d_int,cfg.xdim,cfg.ydim)

            img_pix_arr=image_np2torch(img_pix, cfg, train_device)
            img_int_arr=image_np2torch(img_int, cfg, train_device)
            '''
            img_pix_arr = image_modify(img_pix, cfg)
            img_int_arr = image_modify(img_int, cfg)
            
            img_pix_arr = torch.from_numpy(img_pix_arr.copy())
            img_int_arr = torch.from_numpy(img_int_arr.copy())
            
            img_pix_arr = img_pix_arr.clone().detach()
            img_int_arr = img_int_arr.clone().detach()
            
            img_pix_arr = img_pix_arr.to(train_device).view((-1,1,512,512))
            img_int_arr = img_int_arr.to(train_device).view((-1,1,512,512))
            '''
            #score

            if (img_pix_arr.sum().cpu() < 100):
                score_pix_v = torch.zeros([5])
            else:
                score_pix_v = nn.Sigmoid()(mpid(img_pix_arr)).cpu().detach().numpy()[0]                    

            if (img_int_arr.sum().cpu() < 100):
                score_int_v = torch.zeros([5])
            else:
                score_int_v = nn.Sigmoid()(mpid(img_int_arr)).cpu().detach().numpy()[0]

            '''
            #Plot the image from pgraph
            fig, (ax1,ax2) = plt.subplots(nrows=1, ncols=2, figsize = (18, 6))
            
            ax1.imshow(image_modify(img_pix,cfg).reshape(512,512), cmap='jet')
            ax2.imshow(image_modify(img_int,cfg).reshape(512,512), cmap='jet')
        
            plt.savefig("out/image/%i_%i_%i_%i_graph_plane_%i"%(ev_pix.run(), ev_pix.subrun(), ev_pix.event(), ix, plane))
            '''
            print 'pix scores are ',score_pix_v
            print 'int scores are ',score_int_v

            rd.eminus_pix_score_torch[plane] = score_pix_v[0]
            rd.gamma_pix_score_torch[plane]  = score_pix_v[1]
            rd.muon_pix_score_torch[plane]   = score_pix_v[2]
            rd.pion_pix_score_torch[plane]   = score_pix_v[3]
            rd.proton_pix_score_torch[plane] = score_pix_v[4]

            rd.eminus_int_score_torch[plane] = score_int_v[0]
            rd.gamma_int_score_torch[plane]  = score_int_v[1]
            rd.muon_int_score_torch[plane]   = score_int_v[2]
            rd.pion_int_score_torch[plane]   = score_int_v[3]
            rd.proton_int_score_torch[plane] = score_int_v[4]

            result_dict[rsev][(plane,"pix")] = score_pix_v
            result_dict[rsev][(plane,"int")] = score_int_v
            
            continue
        nfilled += 1
        outtree.Fill()
        rd.reset_vertex()

    print "number of entries filled: ",nfilled
    print "-------------------------------------------"
    if return_result_dict:
        return nfilled,result_dict

    return nfilled
        

def main(IMAGE_FILE,OUT_DIR,CFG,FILEID=0):
    #
    # initialize
    #
    cfg  = config_loader(CFG)
    assert cfg.batch == 1

    print "*********************************"
    print "INFERENCE PID TORCH DLMERGE WC"
    print cfg
    print "*********************************"

    rd = ROOTData()

    FOUT = os.path.join(OUT_DIR,"multipid_out_%s_WC.root" % (FILEID))
    tfile = ROOT.TFile.Open(FOUT,"RECREATE")
    tfile.cd()
    print "OPEN %s"%FOUT

    tree  = ROOT.TTree("multipid_tree","")
    rd.init_tree(tree)
    rd.reset()

    #
    # initialize PyTorch
    #    

    mpid = mpid_net.MPID()
    #mpid.cuda()

    weight_file = ""
    plane=2
    exec("weight_file = cfg.weight_file_mpid_%i" % plane)
    weight_file = weight_file.replace("__WEIGHT_FOLDER__",os.environ["UBMPIDNET_WEIGHT_DIR"])
    print "WEIGHT FILE: ",weight_file
    mpid.load_state_dict(torch.load(weight_file, map_location=train_device))
    mpid.eval()
    
    #
    # initialize iomanager
    #

    # oiom = larcv.IOManager(larcv.IOManager.kWRITE)
    # oiom.set_out_file("trash.root")
    # oiom.initialize()

    iom  = larcv.IOManager(larcv.IOManager.kREAD)
    iom.add_in_file(IMAGE_FILE)
    #iom.add_in_file(VTX_FILE)
    iom.initialize()

    for entry in xrange(iom.get_n_entries()):
        print "@entry={}".format(entry)

        iom.read_entry(entry)

        ev_pgr = iom.get_data(larcv.kProductPGraph,"inter_par")
        ev_par = iom.get_data(larcv.kProductPixel2D,"inter_par_pixel")
        ev_pix = iom.get_data(larcv.kProductPixel2D,"inter_img_pixel")
        ev_int = iom.get_data(larcv.kProductPixel2D,"inter_int_pixel")
        
        print '========================>>>>>>>>>>>>>>>>>>>>'
        print 'run, subrun, event',ev_pix.run(),ev_pix.subrun(),ev_pix.event()

        rd.run[0]    = int(ev_pix.run())
        rd.subrun[0] = int(ev_pix.subrun())
        rd.event[0]  = int(ev_pix.event())
        rd.entry[0]  = int(iom.current_entry())

        rd.num_vertex[0] = int(ev_pgr.PGraphArray().size())

        print 'num of vertices, ',rd.num_vertex[0]
        print 'pgrapgh size, ',int(ev_pgr.PGraphArray().size())
        
        for ix,pgraph in enumerate(ev_pgr.PGraphArray()):
            print "@pgid=%d" % ix
            #if (ix != 2): continue
            rd.vtxid[0] = int(ix)
   
            pgr = ev_pgr.PGraphArray().at(ix)
            cluster_array = pgr.ClusterIndexArray()
            if cluster_array.size()==0:
                cindex_v = []
            else:
                cindex_v = np.array(pgr.ClusterIndexArray())
            
            pixel2d_par_vv = ev_par.Pixel2DClusterArray()
            pixel2d_pix_vv = ev_pix.Pixel2DClusterArray()
            pixel2d_int_vv = ev_int.Pixel2DClusterArray()

            #parid = pgraph.ClusterIndexArray().front()
            #roi0 = pgraph.ParticleArray().front()
            #x = roi0.X()
            #y = roi0.Y()
            #z = roi0.Z()
            #y_2d_plane_0 = ROOT.Double()

            for plane in xrange(3):
                print "@plane=%d" % plane

                #if not plane==2 : continue
                
                if pixel2d_par_vv.size() != 0:
                    rd.npar[plane] = 0
                    pixel2d_par_v = pixel2d_par_vv.at(plane)
                    for cidx in cindex_v:
                        try:
			    pixel2d_par = pixel2d_par_v.at(cidx)
			    if pixel2d_par.size()>0:
			        rd.npar[plane] += 1;
                        except:
                            print "error opening particle pixel2d container"
				                  
                ### Get 2D vertex Image
                #meta = roi0.BB(plane)
                #x_2d = ROOT.Double()
                #y_2d = ROOT.Double()
                #whole_img = ev_img.at(plane)
                #larcv.Project3D(whole_img.meta(), x, y, z, 0.0, plane, x_2d, y_2d)
                #print 'x2d, ', x_2d, 'y2d, ',y_2d
                #if (plane == 0) : y_2d_plane_0 = y_2d
                #else : y_2d = y_2d_plane_0
                ###

                # nothing
                #if pixel2d_vv.empty()==True: continue

                pixel2d_pix_v = pixel2d_pix_vv.at(plane)
                pixel2d_pix = pixel2d_pix_v.at(ix)

		pixel2d_int_v = pixel2d_int_vv.at(plane)
                pixel2d_int = pixel2d_int_v.at(ix)
            
                # nothing on this plane
                #if pixel2d_pix.empty() == True: continue

                rd.inferred[0] = 1
                
                img_pix = larcv.cluster_to_image2d(pixel2d_pix,cfg.xdim,cfg.ydim)
                img_int = larcv.cluster_to_image2d(pixel2d_int,cfg.xdim,cfg.ydim)

                img_pix_arr=image_np2torch(img_pix, cfg, train_device)
                img_int_arr=image_np2torch(img_int, cfg, train_device)
                '''
                img_pix_arr = image_modify(img_pix, cfg)
                img_int_arr = image_modify(img_int, cfg)

                img_pix_arr = torch.from_numpy(img_pix_arr.copy())
                img_int_arr = torch.from_numpy(img_int_arr.copy())

                img_pix_arr = img_pix_arr.clone().detach()
                img_int_arr = img_int_arr.clone().detach()

                img_pix_arr = img_pix_arr.to(train_device).view((-1,1,512,512))
                img_int_arr = img_int_arr.to(train_device).view((-1,1,512,512))
                '''
                #score

                with torch.no_grad():
                    if (img_pix_arr.sum().cpu() < 100):
                        score_pix_v = torch.zeros([5])
                    else:
                        score_pix_v = nn.Sigmoid()(mpid(img_pix_arr)).cpu().detach().numpy()[0]                    

                    if (img_int_arr.sum().cpu() < 100):
                        score_int_v = torch.zeros([5])
                    else:
                        score_int_v = nn.Sigmoid()(mpid(img_int_arr)).cpu().detach().numpy()[0]

                '''
                #Plot the image from pgraph
                fig, (ax1,ax2) = plt.subplots(nrows=1, ncols=2, figsize = (18, 6))

                ax1.imshow(image_modify(img_pix,cfg).reshape(512,512), cmap='jet')
                ax2.imshow(image_modify(img_int,cfg).reshape(512,512), cmap='jet')
                
                plt.savefig("out/image/%i_%i_%i_%i_graph_plane_%i"%(ev_pix.run(), ev_pix.subrun(), ev_pix.event(), ix, plane))
                '''
                print 'pix scores are ',score_pix_v
                print 'int scores are ',score_int_v

                rd.eminus_pix_score_torch[plane] = score_pix_v[0]
                rd.gamma_pix_score_torch[plane]  = score_pix_v[1]
                rd.muon_pix_score_torch[plane]   = score_pix_v[2]
                rd.pion_pix_score_torch[plane]   = score_pix_v[3]
                rd.proton_pix_score_torch[plane] = score_pix_v[4]

                rd.eminus_int_score_torch[plane] = score_int_v[0]
                rd.gamma_int_score_torch[plane]  = score_int_v[1]
                rd.muon_int_score_torch[plane]   = score_int_v[2]
                rd.pion_int_score_torch[plane]   = score_int_v[3]
                rd.proton_int_score_torch[plane] = score_int_v[4]

                continue
                
            tree.Fill()
            rd.reset_vertex()
        
    tfile.cd()
    tree.Write()
    tfile.Close()
    iom.finalize()


if __name__ == '__main__':
    
    if len(sys.argv) != 4:
        print
        print "\tMERGE_FILE = str(sys.argv[1])"
        print "\tOUT_DIR    = str(sys.argv[2])"
        print "\tCFG        = str(sys.argv[3])"
        print 
        sys.exit(1)
    
    IMAGE_FILE = str(sys.argv[1]) 
    OUT_DIR    = str(sys.argv[2])
    CFG        = str(sys.argv[3])

    #CFG = os.path.join(BASE_PATH,"cfg","simple_config.cfg")

    main(IMAGE_FILE,OUT_DIR,CFG)
    
    print "DONE!"
    sys.exit(0)
