"""
-----------------------------------------------------------------------------------------
post_pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Determine roi h5 files
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number
sys.argv[2]: fit model ('gauss','css')
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
python post_fit/post_pp_roi.py 999999 gauss
-----------------------------------------------------------------------------------------
"""

# Stop warnings
# -------------
import warnings
warnings.filterwarnings("ignore")

# General imports
# ---------------
import os
import sys
import json
import glob
import numpy as np
import matplotlib.pyplot as pl
import ipdb
import platform
import h5py
import scipy.io
opj = os.path.join
deb = ipdb.set_trace

# MRI imports
# -----------
import nibabel as nb
import cortex

# Functions import
# ----------------
from utils import set_pycortex_config_file, mask_gii_2_hdf5 


# Get inputs
# ----------
subject = sys.argv[1]
fit_model = sys.argv[2]

# Define analysis parameters
# --------------------------
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define cluster/server specific parameters
# -----------------------------------------
if 'aeneas' in platform.uname()[1]:
    base_dir = analysis_info['aeneas_base_folder'] 
    main_cmd = '/home/szinte/software/workbench/bin_rh_linux64/wb_command' # put the directory for you wb_command
elif 'local' in platform.uname()[1]:
    base_dir = analysis_info['local_base_folder'] 
    main_cmd = '/Applications/workbench/bin_macosx64/wb_command' # put the directory for you wb_command
 
deriv_dir = opj(base_dir,'pp_data',subject,fit_model,'deriv')
roi_masks_dir = opj(base_dir,'pp_data',subject,fit_model,'roi_masks')
h5_dir = opj(base_dir,'pp_data',subject,fit_model,'h5')
try: os.makedirs(roi_masks_dir)
except OSError: pass

# Determine number of vertex and time_serie
# -----------------------------------------
data = []
data_file  =  sorted(glob.glob(opj(base_dir,'raw_data',subject,'*RETBAR1_7T*.func_bla_psc_av.gii')))
data_file_load = nb.load(data_file[0])
data.append(np.array([data_file_load.darrays[i].data for i in range(len(data_file_load.darrays))]))
data = np.vstack(data) 
ts_num,vox_num = data.shape[0],data.shape[1]

# Change cortex database folder
# -----------------------------
pycortex_folder     =   opj(base_dir,'pp_data','cortex')
set_pycortex_config_file(   project_folder = pycortex_folder)

# Create mask from overlay.svg
# ----------------------------
print('creating roi masks from overlay.svg')
masks = cortex.utils.get_roi_verts( subject = 'fsaverage', 
                                    roi = analysis_info['rois'], 
                                    mask = True)
mat_masks = []
for roi in analysis_info['rois']:
    mat_masks.append(masks[roi])
mat_masks = np.vstack(mat_masks)
mat_masks = mat_masks.astype('float32')

prf_deriv_L_all_fsaverage = nb.load(opj(deriv_dir,'all','prf_deriv_L_all_fsaverage.func.gii'))
mat_masks_L = mat_masks[:,0:163842]
darrays = [nb.gifti.gifti.GiftiDataArray(d) for d in mat_masks_L]
gii_out = nb.gifti.gifti.GiftiImage(header = prf_deriv_L_all_fsaverage.header,
                                    extra = prf_deriv_L_all_fsaverage.extra,
                                    darrays = darrays)
nb.save(gii_out,opj(roi_masks_dir,"masks_L_fsaverage.func.gii"))

prf_deriv_R_all_fsaverage = nb.load(opj(deriv_dir,'all','prf_deriv_R_all_fsaverage.func.gii'))
mat_masks_R = mat_masks[:,163842:327684]
darrays = [nb.gifti.gifti.GiftiDataArray(d) for d in mat_masks_R]
gii_out = nb.gifti.gifti.GiftiImage(header = prf_deriv_R_all_fsaverage.header, 
                                    extra = prf_deriv_R_all_fsaverage.extra, 
                                    darrays = darrays)
nb.save(gii_out,opj(roi_masks_dir,"masks_R_fsaverage.func.gii"))

resample_cmd = """{main_cmd} -metric-resample {metric_in} {current_sphere} {new_sphere} ADAP_BARY_AREA {metric_out} -area-metrics {current_area} {new_area}"""
for hemi in ['L','R']:

    current_sphere = opj(base_dir,'raw_data/surfaces/resample_fsaverage','fsaverage_std_sphere.{hemi}.164k_fsavg_{hemi}.surf.gii'.format(hemi=hemi))
    new_sphere = opj(base_dir,'raw_data/surfaces/resample_fsaverage','fs_LR-deformed_to-fsaverage.{hemi}.sphere.{num_vox_k}k_fs_LR.surf.gii'.format(hemi=hemi,num_vox_k = int(np.round(vox_num/1000))))
    current_area = opj(base_dir,'raw_data/surfaces/resample_fsaverage','fsaverage.{hemi}.midthickness_va_avg.164k_fsavg_{hemi}.shape.gii'.format(hemi=hemi))
    new_area = opj(base_dir,'raw_data/surfaces/resample_fsaverage','fs_LR.{hemi}.midthickness_va_avg.{num_vox_k}k_fs_LR.shape.gii'.format(hemi=hemi,num_vox_k = int(np.round(vox_num/1000))))

    metric_in = opj(roi_masks_dir,"masks_{hemi}_fsaverage.func.gii".format(hemi = hemi))
    metric_out = opj(roi_masks_dir,"masks_{hemi}.func.gii".format(hemi = hemi))

    os.system(resample_cmd.format(  main_cmd = main_cmd,
                                    metric_in = metric_in, 
                                    current_sphere = current_sphere, 
                                    new_sphere = new_sphere, 
                                    metric_out = metric_out, 
                                    current_area = current_area, 
                                    new_area = new_area))

# Save ROIS data in hdf5
# ----------------------
print('creating h5 files')
for roi_num, roi in enumerate(analysis_info['rois']):
    try: os.makedirs(h5_dir)
    except OSError: pass

    h5_file = opj(h5_dir,'{roi}.h5'.format(roi = roi))
    try: os.system('rm '+ h5_file)
    except: pass

    for hemi in ['L','R']:

        mask_file = opj(roi_masks_dir,"masks_{hemi}.func.gii".format(hemi = hemi))
        
        for mask_dir in ['all','pos','neg']:
            
            in_file = opj(deriv_dir,mask_dir,"prf_deriv_{hemi}_{mask_dir}.gii".format(hemi = hemi, mask_dir = mask_dir))
            folder_alias = '{hemi}_{mask_dir}'.format(hemi = hemi,mask_dir = mask_dir)
            
            mask_gii_2_hdf5(in_file = in_file,
                            mask_file = mask_file,
                            hdf5_file = h5_file,
                            folder_alias = folder_alias,
                            roi_num = roi_num)