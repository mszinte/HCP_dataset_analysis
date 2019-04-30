"""
-----------------------------------------------------------------------------------------
add_dmn_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Import to pycortex the Yeo atlas to draw the DMN (default mode network) on flatmap
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
python post_fit/add_dmn_roi.py
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
from utils import set_pycortex_config_file, draw_cortex_vertex

# Check system
# ------------
sys.exit('Drawing Flatmaps only works with Python 2. Aborting.') if sys.version_info[0] > 2 else None

# Define analysis parameters
# --------------------------
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define cluster/server specific parameters
# -----------------------------------------
if 'aeneas' in platform.uname()[1]:
    base_dir = analysis_info['aeneas_base_folder'] 
elif 'local' in platform.uname()[1]:
    base_dir = analysis_info['local_base_folder'] 

# Load Yeo2011 atlas
# ------------------
yeo_mat = []
for hemi in ['lh','rh']:

    yeo_hemi_file = opj(base_dir,'pp_data','freesurfer','fsaverage','yeo_label','{hemi}.7Networks_7.gii'.format(hemi = hemi))
    yeo_hemi_load = nb.load(yeo_hemi_file)
    yeo_mat_all = np.array([yeo_hemi_load.darrays[i].data for i in range(len(yeo_hemi_load.darrays))])
    yeo_mat.append(yeo_mat_all[2])


yeo_mat = np.hstack(yeo_mat)
alpha = np.ones((yeo_mat.shape))

# Change cortex database folder
# -----------------------------
pycortex_folder     =   opj(base_dir,'pp_data','cortex')
set_pycortex_config_file(project_folder     =   pycortex_folder)

# Draw maks on fsaverage overlay.svg
# ----------------------------------
param_mask = {	'data': yeo_mat.T, 'cmap': 'Reds', 'alpha': alpha.T,'vmin': 0, 'vmax': 1, 'cbar': 'discrete',\
				'subject': 'fsaverage','add_roi': True, 'roi_name': 'Yeo_7Networks_7'}

vertex_rgb = draw_cortex_vertex(**param_mask)
