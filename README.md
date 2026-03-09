# CAM6_AODs
Fortran code to assimilate MODIS-VIIRS aerosols and nudge through different re-analisys data

these 2 modules are used to nudge assimilated MODIS AODs 
nudging_AODs_qaer_modis_AODVISdn_public.F90
modal_aer_opt_DL_TV_AODVISdn_public.F90

Module modal_aer_opt_DL_TV_AODVISdn_public.F90 contains the modal aerosol optical properties in SW and LW.
D.Leung (2023) dust escheme is also included.
There is subroutine to make AODVISdn availabe and reusable in other modules, in particular nudging_AODs_qaer_modis_AODVISdn_public.F90 to computed the AOD ratio.


Version DF:

nudging_AODs_qaer_modis_AODVISdn_DF.F90 and modal_aer_opt_AODVISdn_DF.F90 contain DF suggestions for the AODBISdn_computed.

Addfld and outfld have been commented in nudging but added in modal_aer.
