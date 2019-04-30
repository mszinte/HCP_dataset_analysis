"""
-----------------------------------------------------------------------------------------
pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Region of interests pre-processing
Compute pRF derivatives and plot on pycortex overlay.svg to determine visual ROI
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number
sys.argv[2]: fit model ('gauss','css')
sys.argv[3]: voxels per fit (e.g 2500)
sys.argv[4]: draw roi in pycortex
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
cd /home/szinte/projects/retino_HCP
python post_fit/pp_roi.py 999999 gauss 2500
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
opj = os.path.join
deb = ipdb.set_trace

# MRI imports
# -----------
import nibabel as nb
import cortex

# Functions import
# ----------------
from utils import set_pycortex_config_file, convert_fit_results

# Check system
# ------------
sys.exit('Popeye error with Python 2. Use Python 3 Aborting.') if sys.version_info[0] < 3 else None

# Get inputs
# ----------
subject = sys.argv[1]
fit_model = sys.argv[2]
job_vox = float(sys.argv[3])
draw_roi = float(sys.argv[4])

if fit_model == 'gauss': fit_val = 6
elif fit_model == 'css': fit_val = 7
base_file_name = 'tfMRI_RETBAR1_7T_AP_Atlas_MSMAll_hp2000_clean.dtseries'

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

# determine number of vertex and time_serie
data = []
data_file  =  sorted(glob.glob(opj(base_dir,'raw_data',subject,'*RETBAR1_7T*.func_bla_psc_av.gii')))
data_file_load = nb.load(data_file[0])
data.append(np.array([data_file_load.darrays[i].data for i in range(len(data_file_load.darrays))]))
data = np.vstack(data) 
ts_num,vox_num = data.shape[0],data.shape[1]

# Check if all slices are present
# -------------------------------
start_idx =  np.arange(0,vox_num,job_vox)
end_idx = start_idx+job_vox
end_idx[-1] = vox_num
num_miss_part = 0
fit_est_files_L = []
fit_est_files_R = []
for hemi in ['L','R']:
    for iter_job in np.arange(0,start_idx.shape[0],1):
        fit_est_file = opj(base_dir,'pp_data',subject,fit_model,'fit', '%s_%s.func_bla_psc_est_%s_to_%s.gii' %(base_file_name,hemi,str(int(start_idx[iter_job])),str(int(end_idx[iter_job]))))
        if os.path.isfile(fit_est_file):
            if os.path.getsize(fit_est_file) == 0:
                num_miss_part += 1 
            else:
                exec('fit_est_files_{hemi}.append(fit_est_file)'.format(hemi = hemi))
        else:
            num_miss_part += 1

if num_miss_part != 0:
    sys.exit('%i missing files, analysis stopped'%num_miss_part)


# Combine fit files
# -----------------
print('combining fit files')
for hemi in ['L','R']:
    data_hemi = np.zeros((fit_val,vox_num))
    exec('fit_est_files_hemi = fit_est_files_{hemi}'.format(hemi=hemi))    
    for fit_filename_hemi in fit_est_files_hemi:
        data_fit_hemi = []
        data_fit_file_hemi = nb.load(fit_filename_hemi)
        data_fit_hemi.append(np.array([data_fit_file_hemi.darrays[i].data for i in range(len(data_fit_file_hemi.darrays))]))
        data_fit_hemi = np.vstack(data_fit_hemi)
        data_hemi = data_hemi + data_fit_hemi

    darrays_est_hemi = [nb.gifti.gifti.GiftiDataArray(d) for d in data_hemi]
    exec('gii_out_{hemi} = nb.gifti.gifti.GiftiImage(header = data_fit_file_hemi.header, extra = data_fit_file_hemi.extra,darrays = darrays_est_hemi)'.format(hemi=hemi))
    exec('nb.save(gii_out_{hemi}, opj(base_dir,"pp_data",subject,fit_model,"fit","{bfn}_{hemi}.func_bla_psc_est.gii"))'.format(hemi=hemi,bfn =base_file_name))

# Compute derived measures from prfs
# ----------------------------------
print('extracting pRF derivatives')
for hemi in ['L','R']:
    prf_filename = sorted(glob.glob(opj(base_dir,'pp_data',subject,fit_model,'fit','%s_%s.func_bla_psc_est.gii'%(base_file_name, hemi))))
    convert_fit_results(prf_filename = prf_filename,
                        output_dir = deriv_dir,
                        stim_radius = analysis_info['stim_radius'],
                        hemi = hemi,
                        fit_model = fit_model)


# Resample gii to fsaverage
# -------------------------
print('converting derivative files to fsaverage')
resample_cmd = """{main_cmd} -metric-resample {metric_in} {current_sphere} {new_sphere} ADAP_BARY_AREA {metric_out} -area-metrics {current_area} {new_area}"""
for hemi in ['L','R']:

    current_sphere = opj(base_dir,'raw_data/surfaces/resample_fsaverage','fs_LR-deformed_to-fsaverage.{hemi}.sphere.{num_vox_k}k_fs_LR.surf.gii'.format(hemi=hemi, num_vox_k = int(np.round(vox_num/1000))))
    new_sphere = opj(base_dir,'raw_data/surfaces/resample_fsaverage','fsaverage_std_sphere.{hemi}.164k_fsavg_{hemi}.surf.gii'.format(hemi=hemi))
    current_area = opj(base_dir,'raw_data/surfaces/resample_fsaverage','fs_LR.{hemi}.midthickness_va_avg.{num_vox_k}k_fs_LR.shape.gii'.format(hemi=hemi,num_vox_k = int(np.round(vox_num/1000))))
    new_area = opj(base_dir,'raw_data/surfaces/resample_fsaverage','fsaverage.{hemi}.midthickness_va_avg.164k_fsavg_{hemi}.shape.gii'.format(hemi=hemi))

    for mask_dir in ['all','pos','neg']:
        
        metric_in = opj(deriv_dir,mask_dir,"prf_deriv_{hemi}_{mask_dir}.gii".format(hemi = hemi, mask_dir = mask_dir))
        metric_out = opj(deriv_dir,mask_dir,"prf_deriv_{hemi}_{mask_dir}_fsaverage.func.gii".format(hemi = hemi, mask_dir = mask_dir))

        os.system(resample_cmd.format(  main_cmd = main_cmd,
                                        metric_in = metric_in, 
                                        current_sphere = current_sphere, 
                                        new_sphere = new_sphere, 
                                        metric_out = metric_out, 
                                        current_area = current_area, 
                                        new_area = new_area))

