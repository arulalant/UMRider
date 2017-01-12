This will produce forecast 00utc files.

$ ./submit_all_bsub_jobs.bash  -> we can submit all bsub jobs via this script

# convert all 2 dimensional IMDAA fields into grib2 file 
bsub < ncum_imdaa_regional_post_reanl_2df_00Z.bash
bsub < ncum_imdaa_regional_post_reanl_2df_06Z.bash
bsub < ncum_imdaa_regional_post_reanl_2df_12Z.bash
bsub < ncum_imdaa_regional_post_reanl_2df_18Z.bash

# convert all 3 dimensional model level IMDAA fields into netcdf4 file 
bsub < ncum_imdaa_regional_post_reanl_mdl_00Z.bash
bsub < ncum_imdaa_regional_post_reanl_mdl_06Z.bash
bsub < ncum_imdaa_regional_post_reanl_mdl_12Z.bash
bsub < ncum_imdaa_regional_post_reanl_mdl_18Z.bash

# convert all 3 dimensional pressure level IMDAA fields into grib2 file 
bsub < ncum_imdaa_regional_post_reanl_prl_00Z.bash
bsub < ncum_imdaa_regional_post_reanl_prl_06Z.bash
bsub < ncum_imdaa_regional_post_reanl_prl_12Z.bash
bsub < ncum_imdaa_regional_post_reanl_prl_18Z.bash



Arulalan.T
