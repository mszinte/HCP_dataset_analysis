## HCP_dataset_analyis

code for analysing HCP dataset
Subject number correspond to their native names

## Analysis specifics
---------------------
- HCP subjects were first pre-processed and averaged
- pRF parameters are extracted using fit/submit_fit_jobs.py on Lisa with gaussian model
- get yeo atlas see post_fit/get_yeo_atlas_gii.txt
- DMN regions (ANG/MED_PAR/SUP_MED_FR/LAT_TEMP) are put on fsaverage flatmap overlay.svg
- DMN regions of interest are drawn in inkscape of fsaverage post_fit/add_dmn_roi.py
- pRF derivatives of subject '999999' are analysed and drawn on fsaverage pycortex map using pp_roi.py
- Vision regions of interest (V1/V2/V3/VO/DO/LO/SUP_PAR/TPJ/sPCS/iPCS/mPCS/INS/DLPFC) are drawn manually in inkscape
- '999999' PRF derivatives summary for each ROI are put in h5 files with post_fit/post_pp_roi.py
- pRF derivatives of all others subject are analysed using post_fit/pp_roi.py
- PRF derivatives summary of all others subject for each ROI are put in h5 files with post_fit/post_pp_roi.py
- Figure 1C is made using post_fit/notebooks/MakeFigure1C.ipynb
- Figure 2A is made using post_fit/notebooks/MakeFigure2A.ipynb
- Figure 2B is made using post_fit/notebooks/MakeFigure2B.ipynb
- Figure 3 and S3 are made using post_fit/notebooks/MakeFigure3.ipynb