import os
import nibabel as nb
import numpy as np
import matplotlib.pyplot as pl
import scipy as sp
from scipy.signal import savgol_filter
from skimage.transform import rotate
from math import *
import cortex

def set_pycortex_config_file(project_folder):

    
    # Import necessary modules
    import os
    import cortex
    from pathlib import Path

    import ipdb

    # Define the new database and colormaps folder
    pycortex_db_folder = project_folder + '/db/'
    pycortex_cm_folder = project_folder + '/colormaps/'

    # Get pycortex config file location
    pycortex_config_file  = cortex.options.usercfg

    # Create name of new config file that will be written
    new_pycortex_config_file = pycortex_config_file[:-4] + '_new.cfg'

    # Create the new config file
    Path(new_pycortex_config_file).touch()

    # Open the config file in read mode and the newly created one in write mode.
    # Loop over every line in the original file and copy it into the new one.
    # For the lines containing either 'filestore' or 'colormap', it will
    # change the saved folder path to the newly created one above (e.g. pycortex_db_folder)
    with open(pycortex_config_file, 'r') as fileIn:
        with open(new_pycortex_config_file, 'w') as fileOut:

            for line in fileIn:

                if 'filestore' in line:
                    newline = 'filestore=' + pycortex_db_folder
                    fileOut.write(newline)
                    newline = '\n'

                elif 'Colormaps' in line:
                    newline = 'Colormaps=' + pycortex_cm_folder
                    fileOut.write(newline)
                    newline = '\n'

                else:
                    newline = line

                fileOut.write(newline)

    
    # Renames the original config file als '_old' and the newly created one to the original name
    os.rename(pycortex_config_file, pycortex_config_file[:-4] + '_old.cfg')
    os.rename(new_pycortex_config_file, pycortex_config_file)
    return None

def convert_fit_results(prf_filename,
                        output_dir,
                        stim_radius,
                        hemi,
                        fit_model):
    """
    Convert pRF fitting value in different parameters for following analysis
   
    Parameters
    ----------
    prf_filename: absolute paths to prf result files.
    output_dir: absolute path to directory into which to put the resulting files.
    stim_radius: stimulus radius in deg
    hemi: brain hemisphere
    fit_model: fit model ('gauss','css')

    Returns
    -------
    prf_deriv_L_all and  prf_deriv_R_all: derivative of pRF analysis for all pRF voxels
    prf_deriv_L_neg and  prf_deriv_R_neg : derivative of pRF analysis for all negative pRF voxels
    prf_deriv_L_pos and  prf_deriv_R_pos : derivative of pRF analysis for all positive pRF voxels

    stucture output:
    columns: 1->32492
    row00 : sign
    row01 : R2
    row02 : eccentricity in deg
    row03 : polar angle real component in deg
    row04 : polar angle imaginary component in deg
    row05 : size in deg
    row06 : non-linerity or nans
    row07 : amplitude
    row08 : baseline
    row09 : coverage
    row10 : x
    row11 : y
    
    ['prf_sign','prf_rsq','prf_ecc','prf_polar_real','prf_polar_imag','prf_size','prf_non_lin','prf_amp','prf_baseline','prf_cov','prf_x','prf_y']

    """

    # Imports
    # -------
    # General imports
    import os
    import nibabel as nb
    import glob
    import numpy as np
    import ipdb
    deb = ipdb.set_trace

    # Popeye imports
    from popeye.spinach import generate_og_receptive_fields


    # Create folders
    # --------------
    try:
        os.makedirs(os.path.join(output_dir,'all'))
        os.makedirs(os.path.join(output_dir,'pos'))
        os.makedirs(os.path.join(output_dir,'neg'))
    except:
        pass

    # Get data details
    # ----------------
    prf_data = []
    prf_data_load = nb.load(prf_filename[0])
    prf_data.append(np.array([prf_data_load.darrays[i].data for i in range(len(prf_data_load.darrays))]))
    prf_data = np.vstack(prf_data)    
    ext = prf_data_load.extra
    hdr = prf_data_load.header

    # Compute derived measures from prfs
    # ----------------------------------
    # get data index
    if fit_model == 'gauss':
        x_idx, y_idx, sigma_idx, beta_idx, baseline_idx, rsq_idx = 0, 1, 2, 3, 4, 5
    elif fit_model == 'css':
        x_idx, y_idx, sigma_idx, non_lin_idx, beta_idx, baseline_idx, rsq_idx = 0, 1, 2, 3, 4, 5, 6
    
    # pRF sign
    prf_sign_all = np.sign((prf_data[beta_idx,:]))
    pos_mask = prf_sign_all > 0.0
    neg_mask = prf_sign_all < 0.0
    all_mask = pos_mask | neg_mask
    
    # r-square
    prf_rsq_all = prf_data[rsq_idx,:]

    # pRF eccentricity
    prf_ecc_all = np.nan_to_num(np.sqrt(prf_data[x_idx,:]**2 + prf_data[y_idx,:]**2))

    # pRF polar angle
    complex_polar = prf_data[x_idx,:] + 1j * prf_data[y_idx,:]
    normed_polar = complex_polar / np.abs(complex_polar)
    prf_polar_real_all = np.real(normed_polar)
    prf_polar_imag_all = np.imag(normed_polar)
    
    # pRF size
    prf_size_all = prf_data[sigma_idx,:].astype(np.float64)
    prf_size_all[prf_size_all<1e-4] = 1e-4

    # pRF non-linearity
    if fit_model == 'gauss':
        prf_non_lin_all = np.zeros((prf_size_all.shape))*np.nan
    elif fit_model == 'css':
        prf_non_lin_all = prf_data[non_lin_idx,:]

    # pRF amplitude
    if fit_model == 'gauss':
        prf_amp_all = prf_data[beta_idx,:]
    elif fit_model == 'css':
        prf_amp_all = prf_data[beta_idx,:]

    # pRF baseline
    if fit_model == 'gauss':
        prf_baseline_all = prf_data[baseline_idx,:]
    elif fit_model == 'css':
        prf_baseline_all = prf_data[baseline_idx,:]

    # pRF coverage
    deg_x, deg_y = np.meshgrid(np.linspace(-30, 30, 121), np.linspace(-30, 30, 121))         # define prfs in visual space
    
    rfs = generate_og_receptive_fields( prf_data[x_idx,:],
                                        prf_data[y_idx,:],
                                        prf_size_all,
                                        np.ones(np.prod(prf_data[0,:].shape[0])),
                                        deg_x,
                                        deg_y)
    if fit_model == 'css':
        rfs = rfs ** prf_data[non_lin_idx,:]

    total_prf_content = rfs.reshape((-1, prf_data.shape[1])).sum(axis=0)
    stim_vignet = np.sqrt(deg_x ** 2 + deg_y**2) < stim_radius    
    prf_cov_all = rfs[stim_vignet, :].sum(axis=0) / total_prf_content

    # pRF x
    prf_x_all = prf_data[x_idx,:]

    # pRF y
    prf_y_all = prf_data[y_idx,:]

    # Saving
    # ------
    for mask_dir in ['all','pos','neg']:
        print('saving: %s'%('os.path.join(output_dir,"{mask_dir}","prf_deriv_{hemi}_{mask_dir}.gii")'.format(hemi = hemi, mask_dir = mask_dir)))
        for output_type in ['prf_sign','prf_rsq','prf_ecc','prf_polar_real','prf_polar_imag','prf_size','prf_non_lin','prf_amp','prf_baseline','prf_cov','prf_x','prf_y']:
            exec('{output_type}_{mask_dir} = np.copy({output_type}_all)'.format(mask_dir = mask_dir, output_type = output_type))
            exec('{output_type}_{mask_dir}[~{mask_dir}_mask] = np.nan'.format(mask_dir = mask_dir, output_type = output_type))
        
        exec('prf_deriv_{mask_dir} = np.row_stack((prf_sign_{mask_dir},prf_rsq_{mask_dir},prf_ecc_{mask_dir},prf_polar_real_{mask_dir},\
                prf_polar_imag_{mask_dir},prf_size_{mask_dir},prf_non_lin_{mask_dir},prf_amp_{mask_dir},prf_baseline_{mask_dir},prf_cov_{mask_dir},\
                prf_x_{mask_dir},prf_y_{mask_dir}))'.format(mask_dir = mask_dir))
        
        exec('prf_deriv_{mask_dir} = prf_deriv_{mask_dir}.astype(np.float32)'.format(mask_dir = mask_dir))
        exec('darrays = [nb.gifti.gifti.GiftiDataArray(d) for d in prf_deriv_{mask_dir}]'.format(mask_dir = mask_dir))
        exec('gii_out = nb.gifti.gifti.GiftiImage(header = hdr, extra = ext, darrays = darrays)')
        exec('nb.save(gii_out,os.path.join(output_dir,"{mask_dir}","prf_deriv_{hemi}_{mask_dir}.gii"))'.format(hemi = hemi, mask_dir = mask_dir))

    return None

def mask_gii_2_hdf5(in_file, mask_file, hdf5_file, folder_alias, roi_num):
    """masks data in in_file with mask in mask_file,
    to be stored in an hdf5 file

    Takes a list of 3D or 4D fMRI nifti-files and masks the
    data with all masks in the list of nifti-files mask_files.
    These files are assumed to represent the same space, i.e.
    that of the functional acquisitions. 
    These are saved in hdf5_file, in the folder folder_alias.

    Parameters
    ----------
    in_files : list
        list of absolute path to functional nifti-files.
        all nifti files are assumed to have the same ndim
    mask_file : list
        list of absolute path to mask nifti-files.
        mask_files are assumed to be 3D
    hdf5_file : str
        absolute path to hdf5 file.
    folder_alias : str
                name of the to-be-created folder in the hdf5 file.
    roi_num: roi row number

    Returns
    -------
    hdf5_file : str
        absolute path to hdf5 file.
    """

    import nibabel as nb
    import os.path as op
    import numpy as np
    import h5py
    import ipdb
    deb = ipdb.set_trace


    gii_in_data = nb.load(in_file)
    data_mat = np.array([gii_in_data.darrays[i].data for i in range(len(gii_in_data.darrays))])
    data_name = op.split(in_file)[-1].split('.gii')[0]

    gii_in_mask = nb.load(mask_file)
    mask_mat = np.array([gii_in_mask.darrays[i].data for i in range(len(gii_in_mask.darrays))])

    mask_mat = mask_mat[roi_num,:]
    mask_mat = np.round(mask_mat)
    
    roi_data = data_mat[:, mask_mat==1]
    
    try:
        h5file = h5py.File(hdf5_file, "r+")
    except:
        h5file = h5py.File(hdf5_file, "a")
    
    try:
        g_hemi = h5file.create_group(folder_alias)    
    except:
        None

    h5file.create_dataset('{folder_alias}/{data_name}'.format(folder_alias = folder_alias, data_name = data_name),data = roi_data,dtype='float32')

    return None
