module nudging_AODs
!=====================================================================
!
! Purpose: Implement Nudging of the model state of AODs
!          toward specified values from MODIS-VIIRS. 
!
! Author : Toni Viudez (toni.v.mora@nasa.gov)
!
! Description:
!    
!    This module assumes that the user has AOD values from MODIS-VIIRS 
!    which have been preprocessed onto the current model grid and adjusted 
!    for differences in topography. It is also assumed that these resulting 
!    values and are stored in individual files which are indexed with respect 
!    to year, month, day, and second of the day. When the model is between 
!    the given begining and ending times, a relaxation forcing is added to 
!    nudge the model toward the analyses values determined from the forcing 
!    option specified. After the model passes the ending analyses time, the 
!    forcing discontinues.
!
!    The nudging of the model toward the MODIS-VIIRS data is controlled by 
!    the 'nudging_AOD_nl' namelist in 'user_nl_cam'; whose variables control the
!    time interval over which nudging is applied, the strength of the nudging
!    tendencies, and its spatial distribution. 
!
!    FORCING:
!    --------

!    Nudging tendencies are applied as a relaxation force between the current
!    model state values and target state values derived from the avalilable
!    analyses. The form of the target values is selected by the 'Nudge_AOD_Force_Opt'
!    option, the timescale of the forcing is determined from the given 
!    'nudge_timescale_opt_AOD', and the nudging strength Alpha=[0.,1.] for each 
!    variable is specified by the 'Nudge_Qaer_coef' values. 
!
!           F_nudge = Alpha*((Target-Model(t_curr))/TimeScale
!
!
!    WINDOWING:
!    ----------
!    The region of applied nudging can be limited using Horizontal/Vertical 
!    window functions that are constructed using a parameterization of the 
!    Heaviside step function. 
!
!    The Heaviside window function is the product of separate horizonal and vertical 
!    windows that are controled via 12 parameters:
!
!        Nudge_AOD_Hwin_lat0:     Specify the horizontal center of the window in degrees. 
!        Nudge_AOD_Hwin_lon0:     The longitude must be in the range [0,360] and the 
!                             latitude should be [-90,+90].
!        Nudge_AOD_Hwin_latWidth: Specify the lat and lon widths of the window as positive 
!        Nudge_AOD_Hwin_lonWidth: values in degrees.Setting a width to a large value (e.g. 999) 
!                             renders the window a constant in that direction.
!        Nudge_AOD_Hwin_latDelta: Controls the sharpness of the window transition with a 
!        Nudge_AOD_Hwin_lonDelta: length in degrees. Small non-zero values yeild a step 
!                             function while a large value yeilds a smoother transition.
!        Nudge_AOD_Hwin_Invert  : A logical flag used to invert the horizontal window function 
!                             to get its compliment.(e.g. to nudge outside a given window).
!
!        Nudge_AOD_Vwin_Lindex:   In the vertical, the window is specified in terms of model 
!        Nudge_AOD_Vwin_Ldelta:   level indcies. The High and Low transition levels should 
!        Nudge_AOD_Vwin_Hindex:   range from [0,(NLEV+1)]. The transition lengths are also 
!        Nudge_AOD_Vwin_Hdelta:   specified in terms of model indices. For a window function 
!                             constant in the vertical, the Low index should be set to 0,
!                             the High index should be set to (NLEV+1), and the transition 
!                             lengths should be set to 0.001 
!        Nudge_AOD_Vwin_Invert  : A logical flag used to invert the vertical window function 
!                             to get its compliment.
!
!        EXAMPLE: For a channel window function centered at the equator and independent 
!                 of the vertical (30 levels):
!                        Nudge_AOD_Hwin_lat0     = 0.         Nudge_AOD_Vwin_Lindex = 0.
!                        Nudge_AOD_Hwin_latWidth = 30.        Nudge_AOD_Vwin_Ldelta = 0.001
!                        Nudge_AOD_Hwin_latDelta = 5.0        Nudge_AOD_Vwin_Hindex = 31.
!                        Nudge_AOD_Hwin_lon0     = 180.       Nudge_AOD_Vwin_Hdelta = 0.001 
!                        Nudge_AOD_Hwin_lonWidth = 999.       Nudge_AOD_Vwin_Invert = .false.
!                        Nudge_AOD_Hwin_lonDelta = 1.0
!                        Nudge_AOD_Hwin_Invert   = .false.
!
!                 If on the other hand one wanted to apply nudging at the poles and
!                 not at the equator, the settings would be similar but with:
!                        Nudge_AOD_Hwin_Invert = .true.
!
!    A user can preview the window resulting from a given set of namelist values before 
!    running the model. Lookat_NudgeWindow.ncl is a script avalable in the tools directory 
!    which will read in the values for a given namelist and display the resulting window.
!
!    The module is currently configured for only 1 window function. It can readily be 
!    extended for multiple windows if the need arises.
!
!
! Input/Output Values:
!    Forcing contributions are available for history file output by 
!    the names:    {'Nudge_Qaer'}
!    The target values that the model state is nudged toward are available for history 
!    file output via the variables:  {'Target_Qaer'}
!
!    &nudging_AOD_nl
!      Nudge_Model_AOD         - LOGICAL toggle to activate nudging.
!                              TRUE  -> Nudging is on.
!                              FALSE -> Nudging is off.                            [DEFAULT]
!
!	   ----------------------------------------------------------------------------------------				
!      Nudge_Path_Qaer          - CHAR path to the aerosols files.
!                              (e.g. '/glade/u/home/tvmora/CATM/MODIS_VIIRS_AOD/$sYYYY/all/')
!
!      Nudge_File_Template_Qaer - CHAR Aerosols filename with year, month, day, and second
!                                 values replaced by %y, %m, %d, and %s respectively.
!                              (e.g. 'terra_modis_npp_viirs_3hr_lat192xlon288.%y%m%d.%s.nc')

!	   ----------------------------------------------------------------------------------------
!
!      Nudge_Times_Per_Day_AOD - INT Number of analyses files available per day.
!                              1 --> daily analyses.
!                              4 --> 6 hourly analyses.
!                              8 --> 3 hourly.
!
!      model_times_per_day_AOD - INT Number of times to update the model state (used for nudging) 
!                                each day. The value is restricted to be longer than the 
!                                current model timestep and shorter than the analyses 
!                                timestep. As this number is increased, the nudging
!                                force has the form of newtonian cooling.
!                              48 --> 1800 Second timestep.
!                              96 -->  900 Second timestep.
!
!      nudge_beg_year_AOD      - INT nudging begining year.  [1979- ]
!      nudge_beg_month_AOD     - INT nudging begining month. [1-12]
!      nudge_beg_day_AOD       - INT nudging begining day.   [1-31]
!      nudge_end_year_AOD      - INT nudging ending year.    [1979-]
!      nudge_end_month_AOD     - INT nudging ending month.   [1-12]
!      nudge_end_day_AOD       - INT nudging ending day.     [1-31]
!
!      Nudge_AOD_Force_Opt     - INT Index to select the nudging Target for a relaxation 
!                                forcing of the form: 
!                                where (t'==Analysis times ; t==Model Times)
!
!                              0 -> NEXT-OBS: Target=Anal(t'_next)                 [DEFAULT]
!                              1 -> LINEAR:   Target=(F*Anal(t'_curr) +(1-F)*Anal(t'_next))
!                                                 F =(t'_next - t_curr )/Tdlt_Anal
!
!      nudge_timescale_opt_AOD - INT Index to select the timescale for nudging.
!                                where (t'==Analysis times ; t==Model Times) 
!
!                              0 -->  TimeScale = 1/Tdlt_Anal                      [DEFAULT]
!                              1 -->  TimeScale = 1/(t'_next - t_curr )
!
!
!      Nudge_Qaer_prof     - INT index of profile structure to use for Qaer.  [0,1,2]
!
!
!                                The spatial distribution is specified with a profile index.
!                                 Where:  0 == OFF      (No Nudging of this variable)
!                                         1 == CONSTANT (Spatially Uniform Nudging)
!                                         2 == HEAVISIDE WINDOW FUNCTION
!
!      Nudge_Qaer_coef     - REAL fractional nudging coeffcient for Qaer(27). 
!
!                                 The strength of the nudging is specified as a fractional 
!                                 coeffcient between [0,1].
!           
!      Nudge_AOD_Hwin_lat0     - REAL latitudinal center of window in degrees.
!      Nudge_AOD_Hwin_lon0     - REAL longitudinal center of window in degrees.
!      Nudge_AOD_Hwin_latWidth - REAL latitudinal width of window in degrees.
!      Nudge_AOD_Hwin_lonWidth - REAL longitudinal width of window in degrees.
!      Nudge_AOD_Hwin_latDelta - REAL latitudinal transition length of window in degrees.
!      Nudge_AOD_Hwin_lonDelta - REAL longitudinal transition length of window in degrees.
!      Nudge_AOD_Hwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!      Nudge_AOD_Vwin_Lindex   - REAL LO model index of transition
!      Nudge_AOD_Vwin_Hindex   - REAL HI model index of transition
!      Nudge_AOD_Vwin_Ldelta   - REAL LO transition length 
!      Nudge_AOD_Vwin_Hdelta   - REAL HI transition length 
!      Nudge_AOD_Vwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!    /
!
!================
!
! TO DO:
! -----------
!    ** Assimilate CERES-MODIS aerosol optical depth and nudge CESM-CAM6-MAM4 model
!
!		* Create namelist of variables 
!       * Qaer = Aerosols constituents
!
!
! ----------------------> How to modify CESM code
!
! CESMHOME /CERES/sarb/aviudezm/CESM2.2
! CAMOD /CERES/sarb/aviudezm/CESM2.2/CAM_modifications
!
! ./create_newcase --case $casename --compset FCnudged --res f09_f09_mg17 --run-unsupported FOR NUDGING
! cd $CESMHOME/cime/scripts/$casename
!
! After create new_case cp -f $CAMOD/nudging.F90 $CESMHOME/cime/scripts/$casename/SourceMods/src.cam/
! ./case.setup
! modify user_nl_cam: output frequency
! modify user_nl_cam: nudging_nl
! ./xmlchange DEBUG=TRUE ! this option enables the debugger 
! ./case.build --skip-provenance-check
! ./case.submit
!          
!=====================================================================
! 
! 	◦	 subroutine nudging_AOD_readnl(nlfile)
! 
! 	◦	 subroutine nudging_init_AOD
! 
! 	◦	 subroutine nudging_timestep_init_AOD(phys_state)
! 
! 	◦	 subroutine nudging_timestep_tend_AOD(phys_state,phys_tend)
! 
! 	◦	 subroutine nudging_update_analyses_fv(aer_file)
! 
! 	◦	 nudging_set_profile(rlat,rlon,Nudge_prof,Wprof,nlev)
!
! 	◦	 nudging_set_PSprofile(rlat,rlon,Nudge_PSprof)
!
! 	◦	 subroutine calc_DryStaticEnergy(t, q, phis, ps, dse, ncol)
! 
!	-----> subroutines to mimic modal_aer_opt.F90 in nudging and compute AODVIS and AOD_ratio
!
! 	◦	 subroutine calc_MidpointPressure(lchnk,ncol)
! 	◦	 subroutine relative_humidity(lchnk, ncol, pver, q, t, pmid)
! 	◦	 subroutine Create_Model_Qaer_Mode(ncol,pver,lchnk)	
! 	◦	 subroutine water_uptake_mam4(ncol,pver,lchnk)  
! 	◦	 subroutine aerosol_sw_properties(lchnk,ncol,ncoef,phys_state
! 	◦	 subroutine qaer_target_mam4(lchnk,ncol,begchunk,endchunk,phys_state) 
! 
!=====================================================================

  ! Useful modules
  !------------------
  use shr_kind_mod,   only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use time_manager,   only:timemgr_time_ge,timemgr_time_inc,get_curr_date,get_step_size,get_calday
  use phys_grid   ,   only:scatter_field_to_chunk
  
  !->SARB    namelist parameters for constituent nudging need pcnst (number of constituents)
  !          constituent names, long name of constituents
  use constituents,   only:pcnst,cnst_name,cnst_longname,cnst_get_ind 
  
  use cam_abortutils, only:endrun
  use spmd_utils  ,   only:masterproc
  use cam_logfile ,   only:iulog
  use radconstants,      only: nswbands
  
  use physconst,         only: rhoh2o, rga, rair
  use radconstants,      only: nswbands, nlwbands, idx_sw_diag, idx_uv_diag, idx_nir_diag
  use rad_constituents,  only: n_diag, rad_cnst_get_call_list, rad_cnst_get_info,  &
                             rad_cnst_get_aer_props, rad_cnst_get_mode_props
  use physics_types,     only: physics_state

  use physics_buffer, only : pbuf_get_index,physics_buffer_desc, pbuf_get_field
  use cam_history,       only: addfld, add_default, outfld, horiz_only
  use cam_history_support, only: fillvalue
  
  use modal_aero_wateruptake, only: modal_aero_wateruptake_dr
  use modal_aero_calcsize,    only: modal_aero_calcsize_diag
  
!   use nudging, only :: nudging_set_profile
  
#ifdef SPMD
  use mpishorthand
#endif

  ! Set all Global values and routines to private by default 
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public:: Nudge_Model_AOD
  public:: Nudge_AOD_ON
  public:: nudging_AOD_readnl
  public:: nudging_init_AOD
  public:: nudging_timestep_init_AOD
  public:: nudging_timestep_tend_AOD
  
  ! ---- I make public our subroutines
  ! public:: water_uptake_mam4
  ! public:: aerosol_sw_properties
  public:: qaer_target_mam4
        
  private::nudging_update_analyses_fv
  private::nudging_set_PSprofile
  private::nudging_set_profile
  private::calc_DryStaticEnergy

  ! Nudging Parameters
  !--------------------
  logical          :: Nudge_Model_AOD       =.false.
  logical          :: Nudge_AOD_ON      =.false.  
  logical          :: Nudge_AOD_Initialized=.false.
  
  !--> SARB
  character(len=cl):: Nudge_Path_Qaer
  character(len=cs):: Nudge_File_Qaer,Nudge_File_Template_Qaer !,Nudge_File_Template_Qaer
  
  integer          :: Nudge_AOD_Force_Opt
  integer          :: Nudge_Timescale_Opt_AOD
  integer          :: Nudge_TSmode
  integer          :: Nudge_Times_Per_Day_AOD
  integer          :: Model_Times_Per_Day_AOD
  
  !--> SARB    add constituent nudging parameters
  real(r8),dimension(27)         :: Nudge_Qaer_coef
  integer,dimension(27)          :: Nudge_Qaer_prof
  
!   integer          :: Nudge_Beg_Date  ,Nudge_End_Date
  
  integer          :: Nudge_Beg_Year_AOD ,Nudge_Beg_Month_AOD
  integer          :: Nudge_Beg_Day_AOD ,Nudge_Beg_Sec_AOD
  integer          :: Nudge_End_Year_AOD ,Nudge_End_Month_AOD
  integer          :: Nudge_End_Day_AOD  ,Nudge_End_Sec_AOD
  
  integer          :: Nudge_AOD_Curr_Year,Nudge_AOD_Curr_Month
  integer          :: Nudge_AOD_Curr_Day ,Nudge_AOD_Curr_Sec
  integer          :: Nudge_AOD_Next_Year,Nudge_AOD_Next_Month
  integer          :: Nudge_AOD_Next_Day ,Nudge_AOD_Next_Sec
  
  integer          :: Nudge_AOD_Step
  
  integer          :: Model_AOD_Curr_Year,Model_AOD_Curr_Month
  integer          :: Model_AOD_Curr_Day ,Model_AOD_Curr_Sec
  integer          :: Model_AOD_Next_Year,Model_AOD_Next_Month
  integer          :: Model_AOD_Next_Day ,Model_AOD_Next_Sec
  
  integer          :: Model_AOD_Step
    
  integer          :: Nudge_Beg_CalDay ,Nudge_End_CalDay ! Nudging Beginning and End Calendar Day   

  real(r8)         :: Nudge_AOD_Hwin_lat0
  real(r8)         :: Nudge_AOD_Hwin_latWidth
  real(r8)         :: Nudge_AOD_Hwin_latDelta
  real(r8)         :: Nudge_AOD_Hwin_lon0
  real(r8)         :: Nudge_AOD_Hwin_lonWidth
  real(r8)         :: Nudge_AOD_Hwin_lonDelta
  logical          :: Nudge_AOD_Hwin_Invert = .false.
  real(r8)         :: Nudge_AOD_Hwin_lo
  real(r8)         :: Nudge_AOD_Hwin_hi
  real(r8)         :: Nudge_AOD_Vwin_Hindex
  real(r8)         :: Nudge_AOD_Vwin_Hdelta
  real(r8)         :: Nudge_AOD_Vwin_Lindex
  real(r8)         :: Nudge_AOD_Vwin_Ldelta
  logical          :: Nudge_AOD_Vwin_Invert =.false.
  real(r8)         :: Nudge_AOD_Vwin_lo
  real(r8)         :: Nudge_AOD_Vwin_hi
  real(r8)         :: Nudge_AOD_Hwin_latWidthH
  real(r8)         :: Nudge_AOD_Hwin_lonWidthH
  real(r8)         :: Nudge_AOD_Hwin_max
  real(r8)         :: Nudge_AOD_Hwin_min

  ! Nudging State Arrays:  TARGET - MODEL - TAU - STEP
  !---------------------------------------------------------------------
  integer Nudge_nlon,Nudge_nlat,Nudge_ncol,Nudge_nlev
  
  real(r8),allocatable::Target_AOD     (:,:)  !(pcols,begchunk:endchunk)
  
  real(r8),allocatable::Target_Qaer    (:,:,:,:)  !(pcols,pver,pcnst,begchunk:endchunk) pcnst=naer
  
!   real(r8),allocatable:: Model_S     (:,:,:)  !(pcols,pver,begchunk:endchunk)  
    
  real(r8),allocatable:: Model_Qaer     (:,:,:,:)  !(pcols,pver,pcnst,begchunk:endchunk) pcnst=naer
  
  real(r8),allocatable::  q_so4_a1     (:,:,:)
  real(r8),allocatable::  q_pom_a1     (:,:,:)
  real(r8),allocatable::  q_soa1_a1     (:,:,:)
  real(r8),allocatable::  q_soa2_a1     (:,:,:)
  real(r8),allocatable::  q_soa3_a1     (:,:,:)
  real(r8),allocatable::  q_soa4_a1     (:,:,:)
  real(r8),allocatable::  q_soa5_a1     (:,:,:)  
  real(r8),allocatable::  q_bc_a1     (:,:,:)
  real(r8),allocatable::  q_dst_a1     (:,:,:)
  real(r8),allocatable::  q_ncl_a1     (:,:,:)
  real(r8),allocatable::  q_num_a1     (:,:,:)
  real(r8),allocatable::  q_so4_a2     (:,:,:)
  real(r8),allocatable::  q_dst_a2     (:,:,:)
  real(r8),allocatable::  q_soa1_a2     (:,:,:)
  real(r8),allocatable::  q_soa2_a2     (:,:,:)
  real(r8),allocatable::  q_soa3_a2     (:,:,:)
  real(r8),allocatable::  q_soa4_a2     (:,:,:)
  real(r8),allocatable::  q_soa5_a2     (:,:,:)
  real(r8),allocatable::  q_ncl_a2     (:,:,:)
  real(r8),allocatable::  q_num_a2     (:,:,:)
  real(r8),allocatable::  q_dst_a3     (:,:,:)
  real(r8),allocatable::  q_ncl_a3     (:,:,:)
  real(r8),allocatable::  q_so4_a3     (:,:,:)
  real(r8),allocatable::  q_num_a3     (:,:,:) 
  real(r8),allocatable::  q_pom_a4     (:,:,:)
  real(r8),allocatable::  q_bc_a4     (:,:,:)
  real(r8),allocatable::  q_num_a4     (:,:,:)
  
  real(r8),allocatable:: Nudge_Qaer_tau  (:,:,:,:)  !(pcols,pver,pcnst,begchunk:endchunk) pcsnt=naer
  
  real(r8),allocatable:: Nudge_Qaer_step (:,:,:,:)  !(pcols,pver,pcnst,begchunk:endchunk) pcnst=naer

  ! Nudging Observation Arrays: OBSERVATIONS
  !-----------------------------------------------------------------------

  integer               Nudge_NumObs_AOD
  integer,allocatable:: Nudge_ObsInd_AOD(:)
  

  logical ,allocatable::Nudge_File_Present_Qaer(:)  
  real(r8),allocatable::Nobs_Qaer (:,:,:,:,:) !(pcols,pver,pcnst,begchunk:endchunk,Nudge_NumObs_AOD)  pcnst=naer

  !--> SARB    chunked array from AOD nudging files

  real(r8),allocatable:: Nobs_AOD (:,:,:) !(pcols,begchunk:endchunk,Nudge_NumObs_AOD) --> CERES-MODIS AODs
  
!   real(r8), allocatable :: mlp(:,:,:)     ! (ncol,pver,lchnk)   midpoint pressures  
! !   real(r8), allocatable :: mlp_c(:,:,:)     ! (ncol,pver,lchnk)   midpoint pressures computed    
!   real(r8), allocatable :: RH(:,:,:)     ! (ncol,pver,lchnk) relative humidity     
  
!   real(r8), allocatable :: lnRH(:,:,:)     ! (ncol,pver,lchnk) relative humidity    
  
!   real(r8), allocatable :: vmr(:,:,:,:)     ! (ncol,pver,nmde,lchnk) dry volume mixing ratio 
!   real(r8), allocatable :: B_mode(:,:,:,:)     ! (ncol,pver,nmde,lchnk)  
!   real(r8), allocatable :: rdry_mode(:,:,:,:)     ! (ncol,pver,nmde,lchnk)  
!   real(r8), allocatable :: Brd(:,:,:,:)     ! (ncol,pver,nmde,lchnk)    
!   real(r8), allocatable :: rwet_mode(:,:,:,:)     ! (ncol,pver,nmde,lchnk)  

!   real(r8), allocatable :: vmr_wet(:,:,:,:)     ! (ncol,pver,nmde,lchnk) wet volume mixing ratio
!   real(r8), allocatable :: vmr_water(:,:,:,:)     ! (ncol,pver,nmde,lchnk) water volume mixing ratio  

!   real(r8), allocatable :: wet_ref_real_index(:,:,:,:)     ! (ncol,pver,nmde,lchnk) Wet refractive index  
!   real(r8), allocatable :: wet_ref_im_index(:,:,:,:)     ! (ncol,pver,nmde,lchnk) Wet refractive index  
      
!   real(r8),allocatable:: Model_Qaer_Mode(:,:,:,:)  !(pcols,pver,pcnst,begchunk:endchunk) pcnst=naer
  
!   real(r8),allocatable:: cext(:,:,:)  !(pcols,ncoef,begchunk:endchunk) 
!   real(r8),allocatable:: cabs(:,:,:)  !(pcols,ncoef,begchunk:endchunk) 
!   real(r8),allocatable:: casm(:,:,:)  !(pcols,ncoef,begchunk:endchunk) 

!   real(r8),allocatable:: pext(:,:)  !(pcols,begchunk:endchunk) 
!   real(r8),allocatable:: pabs(:,:)  !(pcols,begchunk:endchunk) 
!   real(r8),allocatable:: pasm(:,:)  !(pcols,begchunk:endchunk) 
!   real(r8),allocatable:: palb(:,:)  !(pcols,begchunk:endchunk)   
      
!   real(r8),allocatable:: dopaer(:,:,:,:)  !(pcols,pver,nmde,begchunk:endchunk) --- adding #mode
!   real(r8),allocatable:: tauxar(:,:,:)      ! pcols,pver,begchunk:endchunk
  real(r8),allocatable:: AODVISdn_computed(:,:)      ! pcols,begchunk:endchunk  
  real(r8),allocatable:: AOD_ratio(:,:)   !pcols,begchunk:endchunk ratio between MODIS_AOD/AODVISdn_computed
   
   ! ----------------------------------> aerosol advected constituents
   integer so4_a1,pom_a1,soa1_a1,soa2_a1,soa3_a1,soa4_a1,soa5_a1,bc_a1,dst_a1,ncl_a1,num_a1
   integer so4_a2,dst_a2,soa1_a2,soa2_a2,soa3_a2,soa4_a2,soa5_a2,ncl_a2,num_a2
   integer dst_a3,ncl_a3,so4_a3,num_a3
   integer pom_a4,bc_a4,num_a4
   
   ! --- define index of advected aerosols MAM4 from 
   integer,dimension(27) :: idxaer
   integer :: naer = 27
   
   complex(r8) :: crefwsw(nswbands) ! complex refractive index for water visible

   integer, parameter :: pnmde(4)  = (/ 1, 2, 3, 4 /) 
   
   ! Dimension sizes in coefficient arrays used to parameterize aerosol radiative properties
   ! in terms of refractive index and wet radius - compute aerosols coefficients
   integer, parameter :: ncoef=5, prefr=7, prefi=10
   
   logical, public :: modal_strat_sulfate = .false.   ! If .true. then MAM sulfate surface area density used in stratospheric heterogeneous chemistry
contains

  !================================================================
  
  subroutine nudging_AOD_readnl(nlfile)
   ! 
   ! nudging_AOD_readnl: Initialize default values controlling the Nudging 
   !                 process. Then read namelist values to override 
   !                 them. (This subroutine is called in runtime_opts.F90)
   !===============================================================
   use ppgrid        ,only: pver
   
   !--> SARB
   use constituents,   only:pcnst,cnst_name,cnst_longname    
   
   use namelist_utils,only:find_group_name
   use units         ,only:getunit,freeunit
   !
   ! Arguments
   !-------------
   character(len=*),intent(in)::nlfile
   !
   ! Local Values
   !---------------
   integer ierr,unitn, time_of_day
   integer indc !index for all the aerosol mode-species available (27)
 
   namelist /nudging_AOD_nl/ Nudge_Model_AOD,Nudge_Path_Qaer,                       &
                         Nudge_File_Template_Qaer,Nudge_AOD_Force_Opt,              &
                         Nudge_TimeScale_Opt_AOD,                                   &
                         Nudge_Times_Per_Day_AOD,Model_Times_Per_Day_AOD,           &
                         Nudge_Qaer_coef,Nudge_Qaer_prof,                           &
                         Nudge_Beg_Year_AOD,Nudge_Beg_Month_AOD,Nudge_Beg_Day_AOD,  &
                         Nudge_End_Year_AOD,Nudge_End_Month_AOD,Nudge_End_Day_AOD,  &
                         Nudge_AOD_Hwin_lat0,Nudge_AOD_Hwin_lon0,                   &
                         Nudge_AOD_Hwin_latWidth,Nudge_AOD_Hwin_lonWidth,           &
                         Nudge_AOD_Hwin_latDelta,Nudge_AOD_Hwin_lonDelta,           &
                         Nudge_AOD_Hwin_Invert,                                     &
                         Nudge_AOD_Vwin_Lindex,Nudge_AOD_Vwin_Hindex,               &
                         Nudge_AOD_Vwin_Ldelta,Nudge_AOD_Vwin_Hdelta,               &
                         Nudge_AOD_Vwin_Invert                                                   

   ! Initialize time_of_day to first hour:
   time_of_day = 0   

   ! Nudging is NOT initialized yet, For now
   ! Nudging will always begin/end at midnight.
   !--------------------------------------------
   Nudge_AOD_Initialized=.false.
   Nudge_AOD_ON      =.false.   
   Nudge_Beg_Sec_AOD=0
   Nudge_End_Sec_AOD=0

   ! Set Default Namelist values
   !-----------------------------
   Nudge_Model_AOD         = .false.

   Nudge_Path_Qaer          = './Data/YOTC_ne30np4_001/'
   
   Nudge_File_Template_Qaer = 'Qaer_advected.%y%calday.%h.nc' ! Default Aerosols file in YYYYcalday (Year + calendar day)
   
   Nudge_AOD_Force_Opt     = 0
   nudge_timescale_opt_AOD = 0
   Nudge_TSmode        = 0
   Nudge_Times_Per_Day_AOD = 4
   Model_Times_Per_Day_AOD = 4
   
   !--> SARB initialize constituent nudging parameters

   do indc = 1, 27
      Nudge_Qaer_coef(indc) = 0._r8
      Nudge_Qaer_prof(indc) = 0
   end do
   
   Nudge_Beg_Year_AOD      = 2008
   Nudge_Beg_Month_AOD     = 5
   Nudge_Beg_Day_AOD      = 1
   Nudge_End_Year_AOD      = 2008
   Nudge_End_Month_AOD     = 9
   Nudge_End_Day_AOD       = 1
   Nudge_AOD_Hwin_lat0     = 0._r8
   Nudge_AOD_Hwin_latWidth = 9999._r8
   Nudge_AOD_Hwin_latDelta = 1.0_r8
   Nudge_AOD_Hwin_lon0     = 180._r8
   Nudge_AOD_Hwin_lonWidth = 9999._r8
   Nudge_AOD_Hwin_lonDelta = 1.0_r8
   Nudge_AOD_Hwin_Invert   = .false.
   Nudge_AOD_Hwin_lo       = 0.0_r8
   Nudge_AOD_Hwin_hi       = 1.0_r8
   Nudge_AOD_Vwin_Hindex   = float(pver+1)
   Nudge_AOD_Vwin_Hdelta   = 0.001_r8
   Nudge_AOD_Vwin_Lindex   = 0.0_r8
   Nudge_AOD_Vwin_Ldelta   = 0.001_r8
   Nudge_AOD_Vwin_Invert   = .false.
   Nudge_AOD_Vwin_lo       = 0.0_r8
   Nudge_AOD_Vwin_hi       = 1.0_r8

! --- write in log file default values
     write(iulog,*) '|~~~~~~~~~~~~~~~~ default values ~~~~~~~~~~~~~~~~~ >', Nudge_Beg_Year_AOD,Nudge_Beg_Month_AOD,Nudge_Beg_Day_AOD, & 
     Nudge_End_Year_AOD,Nudge_End_Month_AOD,Nudge_End_Day_AOD     
     
   ! Read in namelist values
   !------------------------
   if(masterproc) then
     unitn = getunit()
     open(unitn,file=trim(nlfile),status='old')
     write(iulog,*) '|~~~~~~~~~~~~ nudging_AODs.F90 ~~~~~~~~~ > unitn ', unitn      
     write(iulog,*) '|~~~~~~~~~~~~ nudging_AODs.F90 ~~~~~~~~~ > file=trim(nlfile) ', trim(nlfile)    
     
     call find_group_name(unitn,'nudging_AOD_nl',status=ierr)
     if(ierr.eq.0) then
       read(unitn,nudging_AOD_nl,iostat=ierr)
       if(ierr.ne.0) then
         call endrun('nudging_AOD_readnl:: ERROR reading namelist --- REVIEW IT  ')
       endif
     endif
     close(unitn)
     call freeunit(unitn)
   endif

   ! Set hi/lo values according to the given '_Invert' parameters
   !--------------------------------------------------------------
   if(Nudge_AOD_Hwin_Invert) then
     Nudge_AOD_Hwin_lo = 1.0_r8
     Nudge_AOD_Hwin_hi = 0.0_r8
   else
     Nudge_AOD_Hwin_lo = 0.0_r8
     Nudge_AOD_Hwin_hi = 1.0_r8
   endif

   if(Nudge_AOD_Vwin_Invert) then
     Nudge_AOD_Vwin_lo = 1.0_r8
     Nudge_AOD_Vwin_hi = 0.0_r8
   else
     Nudge_AOD_Vwin_lo = 0.0_r8
     Nudge_AOD_Vwin_hi = 1.0_r8
   endif

   ! Check for valid namelist values 
   !----------------------------------
   if((Nudge_AOD_Hwin_lat0.lt.-90._r8).or.(Nudge_AOD_Hwin_lat0.gt.+90._r8)) then
     write(iulog,*) 'NUDGING AODs: Window lat0 must be in [-90,+90]'
     write(iulog,*) 'NUDGING AODs:  Nudge_AOD_Hwin_lat0=',Nudge_AOD_Hwin_lat0
     call endrun('nudging_AOD_readnl:: ERROR in namelist')
   endif

   if((Nudge_AOD_Hwin_lon0.lt.0._r8).or.(Nudge_AOD_Hwin_lon0.ge.360._r8)) then
     write(iulog,*) 'NUDGING AODs: Window lon0 must be in [0,+360)'
     write(iulog,*) 'NUDGING AODs:  Nudge_AOD_Hwin_lon0=',Nudge_AOD_Hwin_lon0
     call endrun('nudging_AOD_readnl:: ERROR in namelist')
   endif

   if((Nudge_AOD_Vwin_Lindex.gt.Nudge_AOD_Vwin_Hindex)                         .or. &
      (Nudge_AOD_Vwin_Hindex.gt.float(pver+1)).or.(Nudge_AOD_Vwin_Hindex.lt.0._r8).or. &
      (Nudge_AOD_Vwin_Lindex.gt.float(pver+1)).or.(Nudge_AOD_Vwin_Lindex.lt.0._r8)   ) then
     write(iulog,*) 'NUDGING AODs: Window Lindex must be in [0,pver+1]'
     write(iulog,*) 'NUDGING AODs: Window Hindex must be in [0,pver+1]'
     write(iulog,*) 'NUDGING AODs: Lindex must be LE than Hindex'
     write(iulog,*) 'NUDGING AODs:  Nudge_AOD_Vwin_Lindex=',Nudge_AOD_Vwin_Lindex
     write(iulog,*) 'NUDGING AODs:  Nudge_AOD_Vwin_Hindex=',Nudge_AOD_Vwin_Hindex
     call endrun('nudging_AOD_readnl:: ERROR in namelist')
   endif

   if((Nudge_AOD_Hwin_latDelta.le.0._r8).or.(Nudge_AOD_Hwin_lonDelta.le.0._r8).or. &
      (Nudge_AOD_Vwin_Hdelta  .le.0._r8).or.(Nudge_AOD_Vwin_Ldelta  .le.0._r8)    ) then
     write(iulog,*) 'NUDGING AODs: Window Deltas must be positive'
     write(iulog,*) 'NUDGING AODs:  Nudge_AOD_Hwin_latDelta=',Nudge_AOD_Hwin_latDelta
     write(iulog,*) 'NUDGING AODs:  Nudge_AOD_Hwin_lonDelta=',Nudge_AOD_Hwin_lonDelta
     write(iulog,*) 'NUDGING AODs:  Nudge_AOD_Vwin_Hdelta=',Nudge_AOD_Vwin_Hdelta
     write(iulog,*) 'NUDGING AODs:  Nudge_AOD_Vwin_Ldelta=',Nudge_AOD_Vwin_Ldelta
     call endrun('nudging_AOD_readnl:: ERROR in namelist')

   endif

   if((Nudge_AOD_Hwin_latWidth.le.0._r8).or.(Nudge_AOD_Hwin_lonWidth.le.0._r8)) then
     write(iulog,*) 'NUDGING AODs: Window widths must be positive'
     write(iulog,*) 'NUDGING AODs:  Nudge_AOD_Hwin_latWidth=',Nudge_AOD_Hwin_latWidth
     write(iulog,*) 'NUDGING AODs:  Nudge_AOD_Hwin_lonWidth=',Nudge_AOD_Hwin_lonWidth
     call endrun('nudging_AOD_readnl:: ERROR in namelist')
   endif

   ! Broadcast namelist variables
   !------------------------------
#ifdef SPMD
   
   !--> SARB   
   call mpibcast(Nudge_Path_Qaer         ,len(Nudge_Path_Qaer)         ,mpichar,0,mpicom)
   call mpibcast(Nudge_File_Template_Qaer,len(Nudge_File_Template_Qaer),mpichar,0,mpicom)
   
   call mpibcast(Nudge_Model_AOD        , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_AOD_Initialized , 1, mpilog, 0, mpicom)
   
   call mpibcast(Nudge_AOD_ON       , 1, mpilog, 0, mpicom)   
   
   call mpibcast(Nudge_AOD_Force_Opt    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Timescale_Opt_AOD, 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Times_Per_Day_AOD, 1, mpiint, 0, mpicom)
   call mpibcast(Model_Times_Per_Day_AOD, 1, mpiint, 0, mpicom)

   call mpibcast(Nudge_Qaer_coef    , 27, mpir8, 0, mpicom)   
   call mpibcast(Nudge_Qaer_prof   , 27, mpiint, 0, mpicom)
   
   call mpibcast(Nudge_Beg_Year_AOD     , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Month_AOD    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Day_AOD      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Sec_AOD      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Year_AOD     , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Month_AOD    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Day_AOD      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Sec_AOD      , 1, mpiint, 0, mpicom)
   
   call mpibcast(Nudge_AOD_Hwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Hwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Hwin_lat0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Hwin_latWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Hwin_latDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Hwin_lon0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Hwin_lonWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Hwin_lonDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Hwin_Invert,   1, mpilog, 0, mpicom)
   call mpibcast(Nudge_AOD_Vwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Vwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Vwin_Hindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Vwin_Hdelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Vwin_Lindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Vwin_Ldelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Vwin_Invert,   1, mpilog, 0, mpicom)

#endif

   
   ! End Routine
   !------------
   return
  end subroutine ! nudging_AOD_readnl
  !================================================================


   !================================================================
  subroutine nudging_init_AOD
   ! 
   ! nudging_init_AOD: Allocate space and initialize Nudging values
   !===============================================================
   use ppgrid        ,only: pver,pcols,begchunk,endchunk
   use error_messages,only: alloc_err
   use dycore        ,only: dycore_is
   use dyn_grid      ,only: get_horiz_grid_dim_d
   use phys_grid     ,only: get_rlat_p,get_rlon_p,get_ncols_p
   use cam_history   ,only: addfld, add_default, horiz_only
   use shr_const_mod ,only: SHR_CONST_PI
   use filenames     ,only: interpret_filename_spec
   use cam_history_support, only: add_hist_coord

   ! Local values
   !----------------
   integer  Year,Month,Day,Sec
   integer  YMD1,YMD
   logical  After_Beg,Before_End
   integer  istat,lchnk,ncol,icol,ilev
   integer  hdim1_d,hdim2_d
   integer  dtime
   real(r8) rlat,rlon
   real(r8) Wprof(pver)
   real(r8) lonp,lon0,lonn,latp,lat0,latn
   real(r8) Val1_p,Val2_p,Val3_p,Val4_p
   real(r8) Val1_0,Val2_0,Val3_0,Val4_0
   real(r8) Val1_n,Val2_n,Val3_n,Val4_n
   integer               nn
   integer :: nmde = 4
   integer, parameter :: n_ncoef=5, ncoef=5, prefr=7, prefi=10
   integer :: kaer

   ! Get the time step size
   !------------------------
   dtime = get_step_size()

   ! Allocate Space for Nudging data arrays
   !-----------------------------------------
!    write(iulog,*) '|----------------->>>>> Number of constituents = ', pcnst
   !--> SARB   
   allocate(Target_AOD(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','Target_AOD',pcols*((endchunk-begchunk)+1))
   allocate(Target_Qaer(pcols,pver,naer,begchunk:endchunk),stat=istat) ! -----------------> 04/09/25 replace pcnst=naer
!    allocate(Target_Qaer(pcols,pver,pcnst,begchunk:endchunk),stat=istat)   
   call alloc_err(istat,'nudging_init_AOD','Target_Qaer',pcols*((endchunk-begchunk)+1))
   
!    allocate(Model_S(pcols,pver,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','Model_S',pcols*pver*((endchunk-begchunk)+1))
   
   ! --- Allocate space for Model Qaer
   allocate(Model_Qaer(pcols,pver,naer,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','Model_Qaer',pcols*((endchunk-begchunk)+1))
   
   !------------------------------------------------------------------------------------
   ! --- Allocate space for aerosol mass mixing ratio
   
   allocate(q_so4_a1(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_so4_a1',pcols*((endchunk-begchunk)+1))
   
   allocate(q_pom_a1(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_pom_a1',pcols*((endchunk-begchunk)+1))
   
   allocate(q_soa1_a1(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_soa1_a1',pcols*((endchunk-begchunk)+1))
   
   allocate(q_soa2_a1(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_soa2_a1',pcols*((endchunk-begchunk)+1))
   
   allocate(q_soa3_a1(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_soa3_a1',pcols*((endchunk-begchunk)+1))
   
   allocate(q_soa4_a1(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_soa4_a1',pcols*((endchunk-begchunk)+1))

   allocate(q_soa5_a1(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_soa5_a1',pcols*((endchunk-begchunk)+1))
   
   allocate(q_bc_a1(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_bc_a1',pcols*((endchunk-begchunk)+1))
   
   allocate(q_dst_a1(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_dst_a1',pcols*((endchunk-begchunk)+1))
   
   allocate(q_ncl_a1(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_ncl_a1',pcols*((endchunk-begchunk)+1))
   
   allocate(q_num_a1(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_num_a1',pcols*((endchunk-begchunk)+1))
   
   allocate(q_so4_a2(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_so4_a2',pcols*((endchunk-begchunk)+1))
   
   allocate(q_dst_a2(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_dst_a2',pcols*((endchunk-begchunk)+1))
   
   allocate(q_soa1_a2(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_soa1_a2',pcols*((endchunk-begchunk)+1))

   allocate(q_soa2_a2(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_soa2_a2',pcols*((endchunk-begchunk)+1))
   
   allocate(q_soa3_a2(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_soa3_a2',pcols*((endchunk-begchunk)+1))
   
   allocate(q_soa4_a2(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_soa4_a2',pcols*((endchunk-begchunk)+1))
   
   allocate(q_soa5_a2(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_soa5_a2',pcols*((endchunk-begchunk)+1))
   
   allocate(q_ncl_a2(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_ncl_a2',pcols*((endchunk-begchunk)+1))
   
   allocate(q_num_a2(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_num_a2',pcols*((endchunk-begchunk)+1))
   
   allocate(q_dst_a3(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_dst_a3',pcols*((endchunk-begchunk)+1))
   
   allocate(q_ncl_a3(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_ncl_a3',pcols*((endchunk-begchunk)+1))
   
   allocate(q_so4_a3(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_so4_a3',pcols*((endchunk-begchunk)+1))
   
   allocate(q_num_a3(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_num_a3',pcols*((endchunk-begchunk)+1))
   
   allocate(q_pom_a4(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_pom_a4',pcols*((endchunk-begchunk)+1))
   
   allocate(q_bc_a4(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_bc_a4',pcols*((endchunk-begchunk)+1))
   
   allocate(q_num_a4(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','q_num_a4',pcols*((endchunk-begchunk)+1))
   
	!------------------------------------------------------------------------------------
	   
!    ! --- Allocate space for nudging Mid-point pressure
!    allocate(mlp(pcols,pver,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','mlp',pcols*((endchunk-begchunk)+1))

!    ! --- Allocate space for nudging Mid-point pressure computed
! !    allocate(mlp_c(pcols,pver,begchunk:endchunk),stat=istat)
! !    call alloc_err(istat,'nudging_init_AOD','mlp_c',pcols*((endchunk-begchunk)+1))
   
!    ! --- Allocate space for nudging Relative humidity
!    allocate(RH(pcols,pver,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','RH',pcols*((endchunk-begchunk)+1))

!    ! --- Allocate space for nudging Relative humidity
!    allocate(lnRH(pcols,pver,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','lnRH',pcols*((endchunk-begchunk)+1))
         
!    ! --- Allocate space for nudging bulk hygroscopicity
!    allocate(B_mode(pcols,pver,nmde,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','B_mode',pcols*((endchunk-begchunk)+1))

!    ! --- Allocate space for nudging vmr dry volume mixing ratio
!    allocate(vmr(pcols,pver,nmde,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','vmr',pcols*((endchunk-begchunk)+1))   

!    ! --- Allocate space for nudging rdry_mode
!    allocate(rdry_mode(pcols,pver,nmde,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','rdry_mode',pcols*((endchunk-begchunk)+1))   

!    ! --- Allocate space for nudging Brd
!    allocate(Brd(pcols,pver,nmde,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','Brd',pcols*((endchunk-begchunk)+1))   

!    ! --- Allocate space for nudging rwet_mode
!    allocate(rwet_mode(pcols,pver,nmde,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','rwet_mode',pcols*((endchunk-begchunk)+1))   

!    ! --- Allocate space for nudging vmr_wet volume mixing ratio
!    allocate(vmr_wet(pcols,pver,nmde,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','vmr_wet',pcols*((endchunk-begchunk)+1))   

!    ! --- Allocate space for nudging vmr_water volume mixing ratio
!    allocate(vmr_water(pcols,pver,nmde,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','vmr_water',pcols*((endchunk-begchunk)+1))   

!    ! --- Allocate space for nudging wet_ref_real_index wet refractive index
!    allocate(wet_ref_real_index(pcols,pver,nmde,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','wet_ref_real_index',pcols*((endchunk-begchunk)+1))   

!    ! --- Allocate space for nudging wet_ref_im_index wet refractive index
!    allocate(wet_ref_im_index(pcols,pver,nmde,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','wet_ref_im_index',pcols*((endchunk-begchunk)+1))   

!    allocate(Model_Qaer_Mode(pcols,pver,nmde,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','Model_Qaer_mode',pcols*((endchunk-begchunk)+1))

!    ! ##################################################################################
!    ! ### Aerosols optical properties from Ghan paper ##################################
!    ! ##################################################################################
      
!    ! --- Allocate space for cext coefficient -------------------------------------
!    allocate(cext(pcols,n_ncoef,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','cext',pcols*((endchunk-begchunk)+1))

!    ! --- Allocate space for cabs coefficient
!    allocate(cabs(pcols,n_ncoef,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','cabs',pcols*((endchunk-begchunk)+1))

!    ! --- Allocate space for casm coefficient
!    allocate(casm(pcols,n_ncoef,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','casm',pcols*((endchunk-begchunk)+1))

!    ! --- Allocate space for pext coefficient --- parameterized specific extinction (m2/kg)
!    allocate(pext(pcols,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','pext',pcols*((endchunk-begchunk)+1))

!    ! --- Allocate space for pabs coefficient --- parameterized specific absorption (m2/kg)
!    allocate(pabs(pcols,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','pabs',pcols*((endchunk-begchunk)+1))

!    ! --- Allocate space for pasm coefficient --- parameterized asymmetry factor
!    allocate(pasm(pcols,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','pasm',pcols*((endchunk-begchunk)+1))

!    ! --- Allocate space for palb coefficient --- parameterized albedo
!    allocate(palb(pcols,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','palb',pcols*((endchunk-begchunk)+1))

!    ! --- Allocate space for dopaer ---  aerosol optical depth in layer
!    allocate(dopaer(pcols,pver,nmde,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','dopaer',pcols*((endchunk-begchunk)+1))

!    ! -----FINAL AEROSOL OPTICAL PROPERTIES -------------------------------------	
!    ! --- Allocate space for tauxar coefficient -------------------------------------
!    allocate(tauxar(pcols,pver,begchunk:endchunk),stat=istat)
!    call alloc_err(istat,'nudging_init_AOD','tauxar',pcols*(pver+1)*((endchunk-begchunk)+1))

   ! --- Allocate space for AODVISdn_computed -------------------------------
   allocate(AODVISdn_computed(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','AODVISdn_computed',pcols*((endchunk-begchunk)+1))
   
   ! --- Allocate space for AOD nudging scale factor -------------------------------
   allocate(AOD_ratio(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','AOD_ratio',pcols*((endchunk-begchunk)+1))

   
   ! Allocate Space for spatial dependence of 
   ! Nudging Coefs and Nudging Forcing.
   !-------------------------------------------
   
   !--> SARB   
   allocate(Nudge_Qaer_tau(pcols,pver,naer,begchunk:endchunk),stat=istat) ! replace pcnst by naer = 27
   call alloc_err(istat,'nudging_init_AOD','Nudge_Qaer_tau',pcols*((endchunk-begchunk)+1))

   !----------------------------------> SARB   
   allocate(Nudge_Qaer_step(pcols,pver,naer,begchunk:endchunk),stat=istat) ! replace pcnst by naer = 27
   call alloc_err(istat,'nudging_init_AOD','Nudge_Qaer_step',pcols*((endchunk-begchunk)+1))

   write(iulog,*), '|~~~~~~~~~~~~ nudging_AODs.F90 ~~~~~~~~~ >  Using subroutine nudging_init_AOD '

   ! Register output fields with the cam history module
   !-----------------------------------------------------
   !--> SARB  
   ! --- How to define new filed in CAM6:
   !    call addfld(fname, dimnames, avgflag, units, long_name, gridname,         &
   !          flag_xyfill, sampling_seq, standard_name, fill_value)    

   call addfld('Target_AOD',horiz_only,'A','','Nudging Target AOD', flag_xyfill=.true.,fill_value=-999.0_r8)
   write(iulog,*), '|~~~~~~~~~~~~ nudging_AODs.F90 ~~~~~~~~~ > call addfld(Target_AOD)'
   
!    call addfld('aod_weighted_mean',horiz_only,'A','','AODs weighted mean MODIS', flag_xyfill=.true.,fill_value=-999.0_r8)  
              
   ! ---> add history coordinate for each advected aerosol species           
   call add_hist_coord('nmde', nmde, 'Number of MAM-4 (Modal Aerosol Model mode)', '1')
   call add_hist_coord('naer', naer, 'Number of aerosol mode', '1')                            

   call addfld('Nudge_Qaer_tau',(/ 'lev ', 'naer'/),'A','/s','Nudging Qaer Tau', flag_xyfill=.true.,fill_value=-999.0_r8) 
   
   call addfld('Nudge_Qaer_step',(/ 'lev ', 'naer'/),'A','/s','Nudging Qaer Step', flag_xyfill=.true.,fill_value=-999.0_r8) 
                            
   call addfld('Nudge_Qaer',(/ 'lev ', 'naer'/),'A','kg/kg/s','Nudging Qaer Tendency', flag_xyfill=.true.,fill_value=-999.0_r8)
   
   call addfld('Model_Qaer',(/ 'lev ', 'naer'/),'A','kg/kg'  ,'Nudging Model Qaer', flag_xyfill=.true.,fill_value=-999.0_r8)  
    
   call addfld('Target_Qaer',(/ 'lev ', 'naer'/),'A','kg/kg'  ,'Nudging Target Qaer', flag_xyfill=.true.,fill_value=-999.0_r8) 

   ! ----- add history field for each aerosol mass mixing ratio---------------------------------------------------------
   !
   ! Reference -Liu, X., Ma, P.-L., Wang, H., Tilmes, S., Singh, B., Easter, R. C., Ghan, S. J., and Rasch, P. J.:
   ! Description and evaluation of a new four-mode version of the Modal Aerosol Module (MAM4) within version 5.3 of 
   ! the Community Atmosphere Model, Geosci. Model Dev., 9, 505–522, https://doi.org/10.5194/gmd-9-505-2016, 2016.   
   
   ! -----  mode 1: Aitken   
     call addfld('q_so4_a1',(/ 'lev '/),'A','kg/kg'  ,'so4_a1', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_pom_a1',(/ 'lev '/),'A','kg/kg'  ,'pom_a1', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_soa1_a1',(/ 'lev '/),'A','kg/kg'  ,'soa1_a1', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_soa2_a1',(/ 'lev '/),'A','kg/kg'  ,'soa2_a1', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_soa3_a1',(/ 'lev '/),'A','kg/kg'  ,'soa3_a1', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_soa4_a1',(/ 'lev '/),'A','kg/kg'  ,'soa4_a1', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_soa5_a1',(/ 'lev '/),'A','kg/kg'  ,'soa5_a1', flag_xyfill=.true.,fill_value=-999.0_r8)       
     call addfld('q_bc_a1',(/ 'lev '/),'A','kg/kg'  ,'bc_a1', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_dst_a1',(/ 'lev '/),'A','kg/kg'  ,'dst_a1', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_ncl_a1',(/ 'lev '/),'A','kg/kg'  ,'ncl_a1', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_num_a1',(/ 'lev '/),'A','num'  ,'num_a1', flag_xyfill=.true.,fill_value=-999.0_r8)  
     
   ! -----  mode 2: Accumulation     
     call addfld('q_so4_a2',(/ 'lev '/),'A','kg/kg'  ,'so4_a2', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_dst_a2',(/ 'lev '/),'A','kg/kg'  ,'dst_a2', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_soa1_a2',(/ 'lev '/),'A','kg/kg'  ,'soa1_a2', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_soa2_a2',(/ 'lev '/),'A','kg/kg'  ,'soa2_a2', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_soa3_a2',(/ 'lev '/),'A','kg/kg'  ,'soa3_a2', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_soa4_a2',(/ 'lev '/),'A','kg/kg'  ,'soa4_a2', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_soa5_a2',(/ 'lev '/),'A','kg/kg'  ,'soa5_a2', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_ncl_a2',(/ 'lev '/),'A','kg/kg'  ,'ncl_a2', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_num_a2',(/ 'lev '/),'A','num'  ,'num_a2', flag_xyfill=.true.,fill_value=-999.0_r8)  
     
   ! -----  mode 3: Coarse     
     call addfld('q_dst_a3',(/ 'lev '/),'A','kg/kg'  ,'dst_a3', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_ncl_a3',(/ 'lev '/),'A','kg/kg'  ,'ncl_a3', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_so4_a3',(/ 'lev '/),'A','kg/kg'  ,'so4_a3', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_num_a3',(/ 'lev '/),'A','num'  ,'num_a3', flag_xyfill=.true.,fill_value=-999.0_r8) 
     
   ! -----  mode 4: Primary Carbon      
     call addfld('q_pom_a4',(/ 'lev '/),'A','kg/kg'  ,'pom_a4', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_bc_a4',(/ 'lev '/),'A','kg/kg'  ,'bc_a4', flag_xyfill=.true.,fill_value=-999.0_r8)  
     call addfld('q_num_a4',(/ 'lev '/),'A','num'  ,'num_a4', flag_xyfill=.true.,fill_value=-999.0_r8)  
     
   ! ------------------------------------------------------------------------------------------------------------------      	
   
! !    call addfld('mlp_c',(/ 'lev ' /),'A','Pa'  ,'Mid-point pressure computed ', flag_xyfill=.true.,fill_value=-999.0_r8)    
!    call addfld('mlp',(/ 'lev ' /),'A','Pa'  ,'Mid-point pressure  ', flag_xyfill=.true.,fill_value=-999.0_r8) 
!    call addfld('RH',(/ 'lev '/),'A','%'  ,'Relative Humidity  ', flag_xyfill=.true.,fill_value=-999.0_r8) 
!    call addfld('lnRH',(/ 'lev '/),'A',''  ,'Logarithm Relative Humidity  ', flag_xyfill=.true.,fill_value=-999.0_r8)    
   
!    call addfld('vmr',(/ 'lev ', 'nmde'/),'A',''  ,'Dry Volume Mixing Ratio per mode', flag_xyfill=.true.,fill_value=-999.0_r8)    
!    call addfld('B_mode',(/ 'lev ', 'nmde'/),'A',''  ,'Bulk Hygroscopicity  ', flag_xyfill=.true.,fill_value=-999.0_r8) 
!    call addfld('rdry_mode',(/ 'lev ', 'nmde'/),'A','meters'  ,'Dry surface mode radius  ', flag_xyfill=.true.,fill_value=-999.0_r8) 
!    call addfld('Brd',(/ 'lev ', 'nmde'/),'A',''  ,'Bmode / log(RH)  ', flag_xyfill=.true.,fill_value=-999.0_r8)   
      
!    call addfld('rwet_mode',(/ 'lev ', 'nmde'/),'A','meters'  ,'Wet Modal Radius  ', flag_xyfill=.true.,fill_value=-999.0_r8)  
!    call addfld('vmr_wet',(/ 'lev ', 'nmde'/),'A',''  ,'Wet Volume Mixing Ratio per mode', flag_xyfill=.true.,fill_value=-999.0_r8)       
!    call addfld('vmr_water',(/ 'lev ', 'nmde'/),'A',''  ,'Water Volume Mixing Ratio per mode', flag_xyfill=.true.,fill_value=-999.0_r8)    
   
!    call addfld('wet_ref_real_index',(/ 'lev ', 'nmde'/),'A',''  ,'Real wet refractive index per mode @ 550nm', flag_xyfill=.true.,fill_value=-999.0_r8)      
!    call addfld('wet_ref_im_index',(/ 'lev ', 'nmde'/),'A',''  ,'Imaginary wet refractive index per mode @ 550nm', flag_xyfill=.true.,fill_value=-999.0_r8)    
      
!    call addfld('Model_Qaer_Mode',(/ 'lev ', 'nmde'/),'A',''  ,'Qaer Nudging Model by mode (MAM4) ', flag_xyfill=.true.,fill_value=-999.0_r8)       
    
!    call addfld ('cext',   horiz_only, 'A','  ','Extinction coefficient', flag_xyfill=.true.,fill_value=-999.0_r8)
!    call addfld ('cabs',   horiz_only, 'A','  ','Absorption coefficient', flag_xyfill=.true.,fill_value=-999.0_r8)
!    call addfld ('casm',   horiz_only, 'A','  ','Assymetry coefficient',  flag_xyfill=.true.,fill_value=-999.0_r8)
                
!    call addfld ('pext',   horiz_only, 'A','  ','Parameterized specific extinction (m2/kg)', flag_xyfill=.true.,fill_value=-999.0_r8)
!    call addfld ('pabs',   horiz_only, 'A','  ','Parameterized specific absorption (m2/kg)', flag_xyfill=.true.,fill_value=-999.0_r8)
!    call addfld ('pasm',   horiz_only, 'A','  ','Parameterized asymmetry factor',  flag_xyfill=.true.,fill_value=-999.0_r8)
!    call addfld ('palb',   horiz_only, 'A','  ','Parameterized single scattering albedo',  flag_xyfill=.true.,fill_value=-999.0_r8)
  !  call addfld ('tauxar',(/ 'lev '/),'A',' '  ,'Layer extinction optical depth  ', flag_xyfill=.true.,fill_value=-999.0_r8)   
!    call addfld ('dopaer',   (/ 'lev ', 'nmde'/), 'A','  ','Aerosol optical depth in layer',  flag_xyfill=.true.,fill_value=-999.0_r8) 
   call addfld ('AODVISdn_computed',   horiz_only, 'A','  ','AODVISdn computed',  flag_xyfill=.true.,fill_value=-999.0_r8)                 
   call addfld ('AOD_ratio',   horiz_only, 'A','  ','Nudging scale factor MODIS~VIIRS/Computed',  flag_xyfill=.true.,fill_value=-999.0_r8)                    
!    write(iulog,*), '|~~~~~~~~~~~~ nudging_AODs.F90 ~~~~~~~~~ > call addfld(AOD_ratio)'

   !-----------------------------------------  
   ! Values initialized only by masterproc
   !-----------------------------------------

!    write(iulog,*) '|~~~~~~~~~~~~ nudging_AODs.F90 ~~~~~~~~~ > Values initialized only by masterproc'
   if(masterproc) then

     write(iulog,*) '|~~~~~~~~~~~~ nudging_AODs.F90 ~~~~~~~~~ > masterproc ', masterproc         
     write(iulog,*) '|~~~~~~~~~~~~ nudging_AODs.F90 ~~~~~~~~~ > Model_Times_Per_Day_AOD ', Model_Times_Per_Day_AOD        
     write(iulog,*) '|~~~~~~~~~~~~ nudging_AODs.F90 ~~~~~~~~~ > Nudge_Times_Per_Day_AOD', Nudge_Times_Per_Day_AOD        

     write(iulog,*) '|~~~~~~~~~~~~ nudging_AODs.F90 ~~~~~~~~~ > Before'
     
     ! Set the Stepping intervals for Model and Nudging values
     ! Ensure that the Model_AOD_Step is not smaller then one timestep
     !  and not larger then the Nudge_AOD_Step.
     !--------------------------------------------------------
     Model_AOD_Step=86400/Model_Times_Per_Day_AOD
     Nudge_AOD_Step=86400/Nudge_Times_Per_Day_AOD
     
     write(iulog,*) '|~~~~~~~~~~~~ nudging_AODs.F90 ~~~~~~~~~ > After'     
     
     write(iulog,*) '|~~~~~~~~~~~~ nudging_AODs.F90 ~~~~~~~~~ > Model_AOD_Step ', Model_AOD_Step         
     write(iulog,*) '|~~~~~~~~~~~~ nudging_AODs.F90 ~~~~~~~~~ > Nudge_AOD_Step', Nudge_AOD_Step        
               
     if(Model_AOD_Step.lt.dtime) then
       write(iulog,*) ' '
       write(iulog,*) 'NUDGING AODs: Model_AOD_Step cannot be less than a model timestep'
       write(iulog,*) 'NUDGING AODs:  Setting Model_AOD_Step=dtime , dtime=',dtime
       write(iulog,*) ' '
       Model_AOD_Step=dtime
     endif
     if(Model_AOD_Step.gt.Nudge_AOD_Step) then
       write(iulog,*) ' '
       write(iulog,*) 'NUDGING AODs: Model_AOD_Step cannot be more than Nudge_AOD_Step'
       write(iulog,*) 'NUDGING AODs:  Setting Model_AOD_Step=Nudge_AOD_Step, Nudge_AOD_Step=',Nudge_AOD_Step
       write(iulog,*) ' '
       Model_AOD_Step=Nudge_AOD_Step
     endif

     ! Initialize column and level dimensions
     !--------------------------------------------------------
     write(iulog,*) '|~~~~~~~~~~~~ nudging_AODs.F90 ~~~~~~~~~ > Initialize column and level dimensions'           
     call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
     Nudge_nlon=hdim1_d
     Nudge_nlat=hdim2_d
     Nudge_ncol=hdim1_d*hdim2_d
     Nudge_nlev=pver

     ! Initialize column and level dimensions
     !------------------------------------------------
     write(iulog,*) '|~~~~~~~~~~~~ nudging_AODs.F90 ~~~~~~~~~ Initialize column and level dimensions'     
     call get_curr_date(Year,Month,Day,Sec)
     YMD=(Year*10000) + (Month*100) + Day
     YMD1=(Nudge_Beg_Year_AOD*10000) + (Nudge_Beg_Month_AOD*100) + Nudge_Beg_Day_AOD
     call timemgr_time_ge(YMD1,Nudge_Beg_Sec_AOD,         &
                          YMD ,Sec          ,After_Beg)
     YMD1=(Nudge_End_Year_AOD*10000) + (Nudge_End_Month_AOD*100) + Nudge_End_Day_AOD
     call timemgr_time_ge(YMD ,Sec          ,          &
                          YMD1,Nudge_End_Sec_AOD,Before_End)
  
     if((After_Beg).and.(Before_End)) then
       ! Set Time indicies so that the next call to 
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Model_AOD_Next_Year =Year
       Model_AOD_Next_Month=Month
       Model_AOD_Next_Day  =Day
       Model_AOD_Next_Sec  =(Sec/Model_AOD_Step)*Model_AOD_Step
       Nudge_AOD_Next_Year =Year
       Nudge_AOD_Next_Month=Month
       Nudge_AOD_Next_Day  =Day
       Nudge_AOD_Next_Sec  =(Sec/Nudge_AOD_Step)*Nudge_AOD_Step
     elseif(.not.After_Beg) then
       ! Set Time indicies to Nudging start,
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Model_AOD_Next_Year =Nudge_Beg_Year_AOD
       Model_AOD_Next_Month=Nudge_Beg_Month_AOD
       Model_AOD_Next_Day  =Nudge_Beg_Day_AOD
       Model_AOD_Next_Sec  =Nudge_Beg_Sec_AOD
       Nudge_AOD_Next_Year =Nudge_Beg_Year_AOD
       Nudge_AOD_Next_Month=Nudge_Beg_Month_AOD
       Nudge_AOD_Next_Day  =Nudge_Beg_Day_AOD
       Nudge_AOD_Next_Sec  =Nudge_Beg_Sec_AOD
     elseif(.not.Before_End) then
       ! Nudging will never occur, so switch it off
       !--------------------------------------------
       Nudge_Model_AOD=.false.
       Nudge_AOD_ON   =.false.       
       write(iulog,*) ' '
       write(iulog,*) 'NUDGING AODs: WARNING - Nudging has been requested by it will'
       write(iulog,*) 'NUDGING AODs:           never occur for the given time values'
       write(iulog,*) ' '
     endif

     ! Initialize values for window function  
     !----------------------------------------
     lonp= 180._r8
     lon0=   0._r8
     lonn=-180._r8
     latp=  90._r8-Nudge_AOD_Hwin_lat0
     lat0=   0._r8
     latn= -90._r8-Nudge_AOD_Hwin_lat0
    
     Nudge_AOD_Hwin_lonWidthH=Nudge_AOD_Hwin_lonWidth/2._r8
     Nudge_AOD_Hwin_latWidthH=Nudge_AOD_Hwin_latWidth/2._r8

     Val1_p=(1._r8+tanh((Nudge_AOD_Hwin_lonWidthH+lonp)/Nudge_AOD_Hwin_lonDelta))/2._r8
     Val2_p=(1._r8+tanh((Nudge_AOD_Hwin_lonWidthH-lonp)/Nudge_AOD_Hwin_lonDelta))/2._r8
     Val3_p=(1._r8+tanh((Nudge_AOD_Hwin_latWidthH+latp)/Nudge_AOD_Hwin_latDelta))/2._r8
     Val4_p=(1._r8+tanh((Nudge_AOD_Hwin_latWidthH-latp)/Nudge_AOD_Hwin_latDelta))/2._r8
     Val1_0=(1._r8+tanh((Nudge_AOD_Hwin_lonWidthH+lon0)/Nudge_AOD_Hwin_lonDelta))/2._r8
     Val2_0=(1._r8+tanh((Nudge_AOD_Hwin_lonWidthH-lon0)/Nudge_AOD_Hwin_lonDelta))/2._r8
     Val3_0=(1._r8+tanh((Nudge_AOD_Hwin_latWidthH+lat0)/Nudge_AOD_Hwin_latDelta))/2._r8
     Val4_0=(1._r8+tanh((Nudge_AOD_Hwin_latWidthH-lat0)/Nudge_AOD_Hwin_latDelta))/2._r8

     Val1_n=(1._r8+tanh((Nudge_AOD_Hwin_lonWidthH+lonn)/Nudge_AOD_Hwin_lonDelta))/2._r8
     Val2_n=(1._r8+tanh((Nudge_AOD_Hwin_lonWidthH-lonn)/Nudge_AOD_Hwin_lonDelta))/2._r8
     Val3_n=(1._r8+tanh((Nudge_AOD_Hwin_latWidthH+latn)/Nudge_AOD_Hwin_latDelta))/2._r8
     Val4_n=(1._r8+tanh((Nudge_AOD_Hwin_latWidthH-latn)/Nudge_AOD_Hwin_latDelta))/2._r8

     Nudge_AOD_Hwin_max=     Val1_0*Val2_0*Val3_0*Val4_0
     Nudge_AOD_Hwin_min=min((Val1_p*Val2_p*Val3_n*Val4_n), &
                        (Val1_p*Val2_p*Val3_p*Val4_p), &
                        (Val1_n*Val2_n*Val3_n*Val4_n), &
                        (Val1_n*Val2_n*Val3_p*Val4_p))

     ! Initialize number of nudging observation values to keep track of.
     ! Allocate and initialize observation indices 
     !-----------------------------------------------------------------
     if((Nudge_AOD_Force_Opt.ge.0).and.(Nudge_AOD_Force_Opt.le.1)) then
       Nudge_NumObs_AOD=2
     else
       ! Additional Options may need OBS values at more times.
       !------------------------------------------------------
       Nudge_NumObs_AOD=2
       write(iulog,*) 'NUDGING AODs: Setting Nudge_NumObs_AOD=2'
       write(iulog,*) 'NUDGING AODs: WARNING: Unknown Nudge_AOD_Force_Opt=',Nudge_AOD_Force_Opt
       call endrun('NUDGING AODs: Unknown Forcing Option')
     endif
     allocate(Nudge_ObsInd_AOD(Nudge_NumObs_AOD),stat=istat)
     call alloc_err(istat,'nudging_init_AOD','Nudge_ObsInd_AOD',Nudge_NumObs_AOD)     
     allocate(Nudge_File_Present_Qaer(Nudge_NumObs_AOD),stat=istat)
     call alloc_err(istat,'nudging_init_AOD','Nudge_File_Present_Qaer',Nudge_NumObs_AOD)
 
     do nn=1,Nudge_NumObs_AOD
       Nudge_ObsInd_AOD(nn) = Nudge_NumObs_AOD+1-nn
     end do
     !---> SARB 
     Nudge_File_Present_Qaer(:)=.false.

     ! Initialization is done, 
     !--------------------------
     Nudge_AOD_Initialized=.true.

     ! Check that this is a valid DYCORE model
     !------------------------------------------
     if((.not.dycore_is('UNSTRUCTURED')).and. &
        (.not.dycore_is('EUL')         ).and. &
        (.not.dycore_is('LR')          )      ) then
!      if(.not.dycore_is('LR')) then        
       call endrun('NUDGING IS CURRENTLY ONLY CONFIGURED FOR CAM-SE, FV, or EUL')
     endif

     ! Informational Output
     !---------------------------
     write(iulog,*) ' '
     write(iulog,*) '--------------------------------------------------------------'
     write(iulog,*) '  MODEL NUDGING AODs INITIALIZED WITH THE FOLLOWING SETTINGS: '
     write(iulog,*) '--------------------------------------------------------------'
     write(iulog,*) 'NUDGING AODs: Nudge_Model_AOD=',Nudge_Model_AOD

     write(iulog,*) 'NUDGING AODs: Nudge_Path_Qaer=',Nudge_Path_Qaer
     write(iulog,*) 'NUDGING AODs: Nudge_File_Template_Qaer =',Nudge_File_Template_Qaer
     
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Force_Opt=',Nudge_AOD_Force_Opt    
     write(iulog,*) 'NUDGING AODs: nudge_timescale_opt_AOD=',nudge_timescale_opt_AOD    
     write(iulog,*) 'NUDGING AODs: Nudge_TSmode=',Nudge_TSmode
     write(iulog,*) 'NUDGING AODs: Nudge_Times_Per_Day_AOD=',Nudge_Times_Per_Day_AOD
     write(iulog,*) 'NUDGING AODs: Model_Times_Per_Day_AOD=',Model_Times_Per_Day_AOD
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Step=',Nudge_AOD_Step
     write(iulog,*) 'NUDGING AODs: Model_AOD_Step=',Model_AOD_Step

     write(iulog,*) 'NUDGING AODs: Nudge_Qaer_coef =',Nudge_Qaer_coef   
     
     write(iulog,*) 'NUDGING AODs: Nudge_Qaer_prof =',Nudge_Qaer_prof 
       
     write(iulog,*) 'NUDGING AODs: nudge_beg_year_AOD =',nudge_beg_year_AOD
     write(iulog,*) 'NUDGING AODs: nudge_beg_month_AOD=',nudge_beg_month_AOD
     write(iulog,*) 'NUDGING AODs: nudge_beg_day_AOD =',Nudge_Beg_Day_AOD
     write(iulog,*) 'NUDGING AODs: nudge_end_year_AOD =',nudge_end_year_AOD
     write(iulog,*) 'NUDGING AODs: nudge_end_month_AOD=',nudge_end_month_AOD
     write(iulog,*) 'NUDGING AODs: nudge_end_day_AOD  =',nudge_end_day_AOD
     
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Hwin_lat0     =',Nudge_AOD_Hwin_lat0
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Hwin_latWidth =',Nudge_AOD_Hwin_latWidth
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Hwin_latDelta =',Nudge_AOD_Hwin_latDelta
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Hwin_lon0     =',Nudge_AOD_Hwin_lon0
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Hwin_lonWidth =',Nudge_AOD_Hwin_lonWidth
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Hwin_lonDelta =',Nudge_AOD_Hwin_lonDelta
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Hwin_Invert   =',Nudge_AOD_Hwin_Invert  
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Hwin_lo       =',Nudge_AOD_Hwin_lo
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Hwin_hi       =',Nudge_AOD_Hwin_hi
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Vwin_Hindex   =',Nudge_AOD_Vwin_Hindex
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Vwin_Hdelta   =',Nudge_AOD_Vwin_Hdelta
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Vwin_Lindex   =',Nudge_AOD_Vwin_Lindex
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Vwin_Ldelta   =',Nudge_AOD_Vwin_Ldelta
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Vwin_Invert   =',Nudge_AOD_Vwin_Invert  
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Vwin_lo       =',Nudge_AOD_Vwin_lo
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Vwin_hi       =',Nudge_AOD_Vwin_hi
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Hwin_latWidthH=',Nudge_AOD_Hwin_latWidthH
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Hwin_lonWidthH=',Nudge_AOD_Hwin_lonWidthH
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Hwin_max      =',Nudge_AOD_Hwin_max
     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Hwin_min      =',Nudge_AOD_Hwin_min

     write(iulog,*) 'NUDGING AODs: Nudge_AOD_Initialized  =',Nudge_AOD_Initialized
     write(iulog,*) ' '
     write(iulog,*) 'NUDGING AODs: Nudge_NumObs_AOD=',Nudge_NumObs_AOD
     write(iulog,*) ' '

   endif ! (masterproc) then

   ! Broadcast other variables that have changed
   !---------------------------------------------
#ifdef SPMD
   call mpibcast(Model_AOD_Step          ,            1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Step          ,            1, mpir8 , 0, mpicom)
   call mpibcast(Model_AOD_Next_Year     ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_AOD_Next_Month    ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_AOD_Next_Day      ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_AOD_Next_Sec      ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_AOD_Next_Year     ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_AOD_Next_Month    ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_AOD_Next_Day      ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_AOD_Next_Sec      ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Model_AOD         ,            1, mpilog, 0, mpicom)
!    call mpibcast(Nudge_ON            ,            1, mpilog, 0, mpicom)

   call mpibcast(Nudge_AOD_ON        ,            1, mpilog, 0, mpicom)
   
   call mpibcast(Nudge_AOD_Initialized  ,            1, mpilog, 0, mpicom)
   call mpibcast(Nudge_ncol          ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_nlev          ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_nlon          ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_nlat          ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_AOD_Hwin_max      ,            1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Hwin_min      ,            1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Hwin_lonWidthH,            1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_AOD_Hwin_latWidthH,            1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_NumObs_AOD        ,            1, mpiint, 0, mpicom)
#endif

   ! All non-masterproc processes also need to allocate space
   ! before the broadcast of Nudge_NumObs_AOD dependent data.
   !------------------------------------------------------------
   if(.not.masterproc) then
     allocate(Nudge_ObsInd_AOD(Nudge_NumObs_AOD),stat=istat)
     call alloc_err(istat,'nudging_init_AOD','Nudge_ObsInd_AOD',Nudge_NumObs_AOD)
!      allocate(Nudge_File_Present(Nudge_NumObs_AOD),stat=istat)
     call alloc_err(istat,'nudging_init_AOD','Nudge_File_Present',Nudge_NumObs_AOD)
     !---> SARB
     allocate(Nudge_File_Present_Qaer(Nudge_NumObs_AOD),stat=istat)
     call alloc_err(istat,'nudging_init_AOD','Nudge_File_Present_Qaer',Nudge_NumObs_AOD)
     
   endif
#ifdef SPMD
   call mpibcast(Nudge_ObsInd_AOD        , Nudge_NumObs_AOD, mpiint, 0, mpicom)
   !--->SARB
   call mpibcast(Nudge_File_Present_Qaer  , Nudge_NumObs_AOD, mpilog, 0, mpicom) 
     
#endif

   ! Allocate Space for Nudging observation arrays, initialize with 0's
   !---------------------------------------------------------------------
   !--> SARB   
   allocate(Nobs_Qaer(pcols,pver,naer,begchunk:endchunk,Nudge_NumObs_AOD),stat=istat)  ! -----------------> 04/09/25 replace pcnst=naer
!    allocate(Nobs_Qaer(pcols,pver,pcnst,begchunk:endchunk,Nudge_NumObs_AOD),stat=istat)  
   call alloc_err(istat,'nudging_init_AOD','Nobs_Qaer',pcols*((endchunk-begchunk)+1)*Nudge_NumObs_AOD)
   allocate(Nobs_AOD(pcols,begchunk:endchunk,Nudge_NumObs_AOD),stat=istat)
   call alloc_err(istat,'nudging_init_AOD','Nobs_AOD',pcols*((endchunk-begchunk)+1)*Nudge_NumObs_AOD)

   !--> SARB   
!    Nobs_Qaer(:pcols,:pver,:pcnst,begchunk:endchunk,:Nudge_NumObs_AOD)=0._r8
   Nobs_Qaer(:pcols,:pver,:naer,begchunk:endchunk,:Nudge_NumObs_AOD)=0._r8 ! -----------------> 04/09/25 replace pcnst=naer
   Nobs_AOD(:pcols,begchunk:endchunk,:Nudge_NumObs_AOD)=0._r8
   
!!DIAG
   if(masterproc) then
     write(iulog,*) 'NUDGING AODs: nudging_init_AOD() OBS arrays allocated and initialized'
     write(iulog,*) 'NUDGING AODs: nudging_init_AOD() SIZE#',(9*pcols*pver*((endchunk-begchunk)+1)*Nudge_NumObs_AOD)
     write(iulog,*) 'NUDGING AODs: nudging_init_AOD() MB:',float(8*9*pcols*pver*((endchunk-begchunk)+1)*Nudge_NumObs_AOD)/(1024._r8*1024._r8)
     write(iulog,*) 'NUDGING AODs: nudging_init_AOD() pcols=',pcols,' pver=',pver
     write(iulog,*) 'NUDGING AODs: nudging_init_AOD() begchunk:',begchunk,' endchunk=',endchunk
     write(iulog,*) 'NUDGING AODs: nudging_init_AOD() chunk:',(endchunk-begchunk+1),' Nudge_NumObs_AOD=',Nudge_NumObs_AOD
     write(iulog,*) 'NUDGING AODs: nudging_init_AOD() Nudge_ObsInd_AOD=',Nudge_ObsInd_AOD
!      write(iulog,*) 'NUDGING AODs: nudging_init_AOD() Nudge_File_Present=',Nudge_File_Present
   endif
!!DIAG

   ! Initialize the analysis filename at the NEXT time for startup.
   !---------------------------------------------------------------
   Nudge_File_Qaer=interpret_filename_spec(Nudge_File_Template_Qaer      , &
                                       yr_spec=Nudge_AOD_Next_Year , &
                                      mon_spec=Nudge_AOD_Next_Month, &
                                      day_spec=Nudge_AOD_Next_Day  , &
                                      sec_spec=Nudge_AOD_Next_Sec    )


   if(masterproc) then
    write(iulog,*) 'NUDGING AODs: Reading AODs files:',trim(Nudge_Path_Qaer)//trim(Nudge_File_Qaer)
   endif

   ! Rotate Nudge_ObsInd_AOD() indices for new data, then update 
   ! the Nudge observation arrays with analysis data at the 
   ! NEXT==Nudge_ObsInd_AOD(1) time.
   !----------------------------------------------------------
!    if(dycore_is('UNSTRUCTURED')) then
!      call nudging_update_analyses_se (trim(Nudge_Path)//trim(Nudge_File))
!    elseif(dycore_is('EUL')) then
!      call nudging_update_analyses_eul(trim(Nudge_Path)//trim(Nudge_File))
!    else !if(dycore_is('LR')) then
!    if(dycore_is('LR')) then
     call nudging_update_analyses_fv ((trim(Nudge_Path_Qaer)//trim(Nudge_File_Qaer)))   
!    endif

   ! Initialize Nudging Coefficient profiles in local arrays
   ! Load zeros into nudging arrays
   !------------------------------------------------------
   do lchnk=begchunk,endchunk
     ncol=get_ncols_p(lchnk)
     do icol=1,ncol
       rlat=get_rlat_p(lchnk,icol)*180._r8/SHR_CONST_PI
       rlon=get_rlon_p(lchnk,icol)*180._r8/SHR_CONST_PI

    !--> SARB       

      do kaer = 1,27
!           call nudging_set_profile(rlat,rlon,1,Wprof,pver)       ! --- original from nudgnig.F90 and alwyas used
!           call nudging_set_profile(rlat,rlon,2,Wprof,pver)       
          call nudging_set_profile(rlat,rlon,nudge_Qaer_prof(kaer),Wprof,pver) ! checking version if reads the namelist (test_Nudge_prof in Derecho)     
	      Nudge_Qaer_tau(icol,:,kaer,lchnk)=Wprof(:)  ! ---last version used
	  end do  ! --- naer
	         
     end do ! --- ncols
          
	do kaer = 1,27 
     Nudge_Qaer_tau(:ncol,:pver,kaer,lchnk) =                             &
          Nudge_Qaer_tau(:ncol,:pver,kaer,lchnk) * Nudge_Qaer_coef(kaer)/float(Nudge_AOD_Step) ! Nudge_Qcoef should be Nudge_Qaer_coef
	end do  ! --- naer
     
     !--> SARB
     Nudge_Qaer_step(:pcols,:pver,:naer,lchnk)=0._r8
          
     !--> SARB
     Target_Qaer(:pcols,:pver,:naer,lchnk)=0._r8
     Target_AOD(:pcols,lchnk)=-999.0_r8

! !      mlp_c(:pcols,:pver, lchnk)=0._r8     
!      mlp(:pcols,:pver, lchnk)=0._r8     
!      RH(:pcols,:pver, lchnk)=0._r8           
!      lnRH(:pcols,:pver, lchnk)=0._r8          
!      B_mode(:pcols,:pver,:nmde, lchnk)=0._r8           
! !      rwet_mode(:pcols,:pver,:naer)=0._r8    

!      wet_ref_real_index(:pcols,:pver,:nmde,lchnk)=0._r8    
!      wet_ref_im_index(:pcols,:pver,:nmde,lchnk)=0._r8    
          
!      Model_Qaer_Mode(:pcols,:pver,:nmde, lchnk)=0._r8    
	 AOD_ratio(:pcols,lchnk)=1._r8           

   end do ! --- lchnk
 
   ! End Routine
   !------------
   return
  end subroutine ! nudging_init_AOD
  !================================================================
                
  !================================================================
  subroutine nudging_timestep_init_AOD(phys_state)
   ! 
   ! NUDGING_TIMESTEP_INIT: 
   !                 Check the current time and update Model/Nudging 
   !                 arrays when necessary. Toggle the Nudging flag
   !                 when the time is withing the nudging window.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state
   use constituents ,only: cnst_get_ind
   use dycore       ,only: dycore_is
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use filenames    ,only: interpret_filename_spec
   use ESMF

   !--> Call subroutine Model_Qaer(mode_idx, spec_idx, idx) in rad_constituents.F90
   !        With this call we can get the index (idx)
   !        get constituent index of mam specie mmr (climate list only)
   !
   ! Return constituent index of mam specie mass mixing ratio for aerosol modes in
   ! the climate list.
   ! This is a special routine to allow direct access to information in the 
   ! constituent array inside physics parameterizations that have been passed,
   ! and are operating over the entire constituent array.  The interstitial phase
   ! is assumed since that's what is contained in the constituent array.
   !
   ! Arguments
   ! mode_idx = mode index
   ! spec_idx = index of specie in the mode
   ! idx      = index of specie in the constituent array   
   
!    use rad_constituents  ,only: rad_cnst_get_mam_mmr_idx
   use constituents      ,only: pcnst,cnst_name,cnst_longname 
   
   ! Call subroutines from modal_aer_opt.F90 (rad_constituents:rad_cnst_get_mode_props)
   ! From rad_constituents.F90 use rad_cnst_get_mode_props
   ! Return requested properties for the mode from the specified climate or diagnostic list.

   use rad_constituents,  only: rad_cnst_get_mode_props
   use ppgrid,            only: pcols, pver, pverp     
   
   ! --- modules used by call aerosol_sw_properties
!    use modal_aer_opt ,only: binterp ! bilinear interpolation of table
  !  use modal_aer_opt, only: get_aodvisdn  
   use ppgrid,            only: pver, pverp    
   use physconst,         only: rhoh2o, rga, rair
   
      
   ! Arguments
   !-----------
   type(physics_state),intent(in):: phys_state(begchunk:endchunk)
  !  real(r8), intent(in) :: pmid(:,:)  ! (ncol,pver) mid-layer pressure
   ! Local values
   !----------------
   integer Year,Month,Day,Sec
   integer YMD1,YMD2,YMD
   logical Update_Model,Update_Nudge,Sync_Error
   logical After_Beg   ,Before_End
   integer lchnk,ncol,indw
   integer j   
!   ! --- define index of advected aerosols MAM4 from 
!   integer,dimension(27) :: idxaer
!   integer naer
   
   ! --- Define variables for rad_cnst_get_mode_props
   integer :: list_idx=0  ! index of the climate or a diagnostic list
   integer :: mode_idx  ! mode index

   real(r8), pointer :: refitabsw(:,:)  ! table of imag refractive indices for aerosols
   real(r8), pointer :: refrtabsw(:,:)  ! table of real refractive indices for aerosols
   real(r8), pointer :: abspsw(:,:,:,:) ! specific absorption
   real(r8), pointer :: asmpsw(:,:,:,:) ! asymmetry factor
   real(r8), pointer :: extpsw(:,:,:,:) ! specific extinction
   real(r8) :: sigma_logr_aer           ! geometric standard deviation of number distribution
   
   integer, parameter :: ncoef=5, prefr=7, prefi=10  
    
   real(r8) :: cext!(:,:)
   real(r8) :: cabs!(:,:)
   real(r8) :: casm!(:,:)   

   real(r8) :: pext!(:,:)
   real(r8) :: pabs!(:,:)
   real(r8) :: pasm!(:,:)   

   real(r8) :: tauxar!(:,:)
      
   type(ESMF_Time)         Date1,Date2
   type(ESMF_TimeInterval) DateDiff
   integer                 DeltaT
   real(r8)                Tscale
   real(r8)                Tfrac
   integer                 rc
   integer                 nn
   integer                 kk
   real(r8)                Sbar,Qbar,Wsum
   integer                 dtime
   integer kaer
   ! --- > local values for indexes in the use of  mode_idx and spec,idx
!    integer indc
   integer aer_mod, aer_spec
   integer cnst ! init for writing constituents, names and long names loop
!    integer mode_idx,spec_idx,idx ! parameters in rad_cnst_get_mam_mmr_idx
   
   ! Check if Nudging is initialized
   !---------------------------------
   if(.not.Nudge_AOD_Initialized) then
     call endrun('nudging_timestep_init_AOD:: Nudging NOT Initialized')
   endif

   ! Get time step size
   !--------------------
   dtime = get_step_size()

   ! Get Current time
   !--------------------
   call get_curr_date(Year,Month,Day,Sec)
   YMD=(Year*10000) + (Month*100) + Day

   !-------------------------------------------------------
   ! Determine if the current time is AFTER the begining time
   ! and if it is BEFORE the ending time.
   !-------------------------------------------------------
   YMD1=(nudge_beg_year_AOD*10000) + (nudge_beg_month_AOD*100) + Nudge_Beg_Day_AOD
   call timemgr_time_ge(YMD1,Nudge_Beg_Sec_AOD,         &
                        YMD ,Sec          ,After_Beg)

   YMD1=(nudge_end_year_AOD*10000) + (nudge_end_month_AOD*100) + nudge_end_day_AOD
   call timemgr_time_ge(YMD ,Sec,                    &
                        YMD1,Nudge_End_Sec_AOD,Before_End)

   !--------------------------------------------------------------
   ! When past the NEXT time, Update Model Arrays and time indices
   !--------------------------------------------------------------
   YMD1=(Model_AOD_Next_Year*10000) + (Model_AOD_Next_Month*100) + Model_AOD_Next_Day
   call timemgr_time_ge(YMD1,Model_AOD_Next_Sec,            &
                        YMD ,Sec           ,Update_Model)

   if((Before_End).and.(Update_Model)) then
     ! Increment the Model times by the current interval
     !---------------------------------------------------
     Model_AOD_Curr_Year =Model_AOD_Next_Year
     Model_AOD_Curr_Month=Model_AOD_Next_Month
     Model_AOD_Curr_Day  =Model_AOD_Next_Day
     Model_AOD_Curr_Sec  =Model_AOD_Next_Sec
     YMD1=(Model_AOD_Curr_Year*10000) + (Model_AOD_Curr_Month*100) + Model_AOD_Curr_Day
     call timemgr_time_inc(YMD1,Model_AOD_Curr_Sec,              &
                           YMD2,Model_AOD_Next_Sec,Model_AOD_Step,0,0)

     ! Check for Sync Error where NEXT model time after the update
     ! is before the current time. If so, reset the next model 
     ! time to a Model_AOD_Step after the current time.
     !--------------------------------------------------------------
     call timemgr_time_ge(YMD2,Model_AOD_Next_Sec,            &
                          YMD ,Sec           ,Sync_Error)
     if(Sync_Error) then
       Model_AOD_Curr_Year =Year
       Model_AOD_Curr_Month=Month
       Model_AOD_Curr_Day  =Day
       Model_AOD_Curr_Sec  =Sec
       call timemgr_time_inc(YMD ,Model_AOD_Curr_Sec,              &
                             YMD2,Model_AOD_Next_Sec,Model_AOD_Step,0,0)
       write(iulog,*) 'NUDGING AODs: WARNING - Model_Time Sync ERROR... CORRECTED'
     endif
     Model_AOD_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(Model_AOD_Next_Year*10000)
     Model_AOD_Next_Month=(YMD2/100)
     Model_AOD_Next_Day  = YMD2-(Model_AOD_Next_Month*100)

     ! Load values at Current into the Model arrays
     !-----------------------------------------------
     call cnst_get_ind('Q',indw)
     !write(*,*) 'index water vapor constituent: ',  indw
     
         
     !--> SARB call subroutine to get aerosols index for specific mode-spec
     !
     ! Aerosols in CAM-Chem are represented using the 4-mode version of the modal aerosol model (MAM-4) from Lui et al.,2016.
     ! Background and a description of the basic aerosol model performance can be found in Liu et al.,2012.  
     ! The 4 modes are Accumulation (a1), Aitken (a2), Coarse (a3) and Primary Carbon (a4).  
     ! The species that are included in the base chemical mechanism are:
     !		Black Carbon (bc)
     !		Primary Organic Matter (pom)
     !		Sulfate (so4)
     !		Dust (dst)
     !		Sea Salt (ncl)
     !		Secondary Organic Aerosol (soa)
     !		Number (num)
     ! Secondary organic aerosol is treated using a volatility basis set (VBS) scheme, with 5 volatility bins, derived from Hodzic et al.,2016,
     ! including updates to the chemical reaction rates from GECKO-A.
     !
     ! Each mode is conformed by the following species:  
     !
            
     !call rad_cnst_get_mam_mmr_idx(###,pcnt)
     !call rad_cnst_get_mam_mmr_idx(mode_idx, spec_idx, idx)

     ! --- aerosols mode 1- Accumulation
     call cnst_get_ind('so4_a1',so4_a1)     
     call cnst_get_ind('pom_a1',pom_a1)         
     call cnst_get_ind('soa1_a1',soa1_a1)    
     call cnst_get_ind('soa2_a1',soa2_a1)          
     call cnst_get_ind('soa3_a1',soa3_a1)         
     call cnst_get_ind('soa4_a1',soa4_a1)         
     call cnst_get_ind('soa5_a1',soa5_a1)              
     call cnst_get_ind('bc_a1',bc_a1)         
     call cnst_get_ind('dst_a1',dst_a1)     
     call cnst_get_ind('ncl_a1',ncl_a1)     
     call cnst_get_ind('num_a1',num_a1)     
     
     ! --- aerosols mode 2 - Aitken    
     call cnst_get_ind('so4_a2',so4_a2)         
     call cnst_get_ind('soa1_a2',soa1_a2)     
     call cnst_get_ind('soa2_a2',soa2_a2)         
     call cnst_get_ind('soa3_a2',soa3_a2)         
     call cnst_get_ind('soa4_a2',soa4_a2)         
     call cnst_get_ind('soa5_a2',soa5_a2)     
     call cnst_get_ind('dst_a2',dst_a2)             
     call cnst_get_ind('ncl_a2',ncl_a2)     
     call cnst_get_ind('num_a2',num_a2)     
       
     ! --- aerosols mode 3 - Coarse      
     call cnst_get_ind('so4_a3',so4_a3)      
     call cnst_get_ind('dst_a3',dst_a3)               
     call cnst_get_ind('ncl_a3',ncl_a3)     
     call cnst_get_ind('num_a3',num_a3)      
     
     ! --- aerosols mode 4 - Primary Carbon      
     call cnst_get_ind('pom_a4',pom_a4)         
     call cnst_get_ind('bc_a4',bc_a4)      
     call cnst_get_ind('num_a4',num_a4)     
     
     ! --- array for all the advected aerosols in MAM4

     idxaer =  [so4_a1,pom_a1,soa1_a1,soa2_a1,soa3_a1,soa4_a1,soa5_a1,bc_a1,dst_a1,ncl_a1,num_a1, &
                so4_a2,dst_a2,soa1_a2,soa2_a2,soa3_a2,soa4_a2,soa5_a2,ncl_a2,num_a2, &
                dst_a3,ncl_a3,so4_a3,num_a3, &
                pom_a4,bc_a4,num_a4]
                
!      write(iulog,*) 'index aerosols in MAM4', idxaer
     
     do lchnk=begchunk,endchunk
       ncol=phys_state(lchnk)%ncol

       ! --- creating mid-point pressure values:
      !  mlp(:ncol,:pver,lchnk)=phys_state(lchnk)%pmid(:ncol,:pver)

       ! --- Populating the Model_Qaer array with values for all the aerosols 
       !     in the q state array based on the index of each aerosol
       !     mixing ratio (kg/kg moist or dry air depending on type (use physics_types, only:: physics_state])
       
       do kaer=1,27
          Model_Qaer(:ncol,:pver,kaer,lchnk)=phys_state(lchnk)%q(:ncol,:pver,idxaer(kaer))       
       end do     
       
       
     end do

!      ! Load Dry Static Energy values for Model
!      !-----------------------------------------
!      if(Nudge_TSmode.eq.0) then
!        ! DSE tendencies from Temperature only
!        !---------------------------------------
!        do lchnk=begchunk,endchunk
!          ncol=phys_state(lchnk)%ncol
!          Model_S(:ncol,:pver,lchnk)=cpair*Model_T(:ncol,:pver,lchnk)
!        end do
!      elseif(Nudge_TSmode.eq.1) then
!        ! Calculate DSE tendencies from Temperature, Water Vapor, and Surface Pressure
!        !------------------------------------------------------------------------------
!        do lchnk=begchunk,endchunk
!          ncol=phys_state(lchnk)%ncol
!          call calc_DryStaticEnergy(Model_T(:,:,lchnk)  , Model_Q(:,:,lchnk), &
!                                  phys_state(lchnk)%phis,  Model_PS(:,lchnk), &
!                                                   Model_S(:,:,lchnk), ncol)
!        end do
!      endif 
   endif ! ((Before_End).and.(Update_Model)) then

   !----------------------------------------------------------------
   ! When past the NEXT time, Update Nudging Arrays and time indices
   !----------------------------------------------------------------
   YMD1=(Nudge_AOD_Next_Year*10000) + (Nudge_AOD_Next_Month*100) + Nudge_AOD_Next_Day
   call timemgr_time_ge(YMD1,Nudge_AOD_Next_Sec,            &
                        YMD ,Sec           ,Update_Nudge)

   if((Before_End).and.(Update_Nudge)) then
     ! Increment the Nudge times by the current interval
     !---------------------------------------------------
     Nudge_AOD_Curr_Year =Nudge_AOD_Next_Year
     Nudge_AOD_Curr_Month=Nudge_AOD_Next_Month
     Nudge_AOD_Curr_Day  =Nudge_AOD_Next_Day
     Nudge_AOD_Curr_Sec  =Nudge_AOD_Next_Sec
     YMD1=(Nudge_AOD_Curr_Year*10000) + (Nudge_AOD_Curr_Month*100) + Nudge_AOD_Curr_Day
     call timemgr_time_inc(YMD1,Nudge_AOD_Curr_Sec,              &
                           YMD2,Nudge_AOD_Next_Sec,Nudge_AOD_Step,0,0)
     Nudge_AOD_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(Nudge_AOD_Next_Year*10000)
     Nudge_AOD_Next_Month=(YMD2/100)
     Nudge_AOD_Next_Day  = YMD2-(Nudge_AOD_Next_Month*100)

     ! Set the analysis filename at the NEXT time.
     !---------------------------------------------------------------
     Nudge_File_Qaer=interpret_filename_spec(Nudge_File_Template_Qaer      , &
                                         yr_spec=Nudge_AOD_Next_Year , &
                                        mon_spec=Nudge_AOD_Next_Month, &
                                        day_spec=Nudge_AOD_Next_Day  , &
                                        sec_spec=Nudge_AOD_Next_Sec    )
     if(masterproc) then
      write(iulog,*) '!------------------------------------------------------------------------------'
      write(iulog,*) 'NUDGING AODs: Reading AOD file:',trim(Nudge_Path_Qaer)//trim(Nudge_File_Qaer)
      write(iulog,*) '!------------------------------------------------------------------------------'      
     endif

     ! Rotate Nudge_ObsInd_AOD() indices for new data, then update 
     ! the Nudge observation arrays with analysis data at the 
     ! NEXT==Nudge_ObsInd_AOD(1) time.
     !----------------------------------------------------------
     if(dycore_is('LR')) then
      !  write(iulog,*) '!---------------------- Calling nudging_update_analyses_fv: Finite-Volume Cubed Sphere --------------------------------'
       call nudging_update_analyses_fv (trim(Nudge_Path_Qaer)//trim(Nudge_File_Qaer))
     endif
   endif ! ((Before_End).and.(Update_Nudge)) then

   
   
   ! --------->>>>>>>>> NEW TEST >>>>>>>>>     
   !----------------------------------------------------------------
   ! Toggle Nudging flag when the time interval is between 
   ! beginning and ending times, and all of the analyses files exist.
   !
   !			AODs from MODIS-VIIRS
   !   
   !----------------------------------------------------------------
   if((After_Beg).and.(Before_End)) then
     if(Nudge_AOD_Force_Opt.eq.0) then
       ! Verify that the NEXT analyses are available
       !---------------------------------------------
       Nudge_AOD_ON=Nudge_File_Present_Qaer(Nudge_ObsInd_AOD(1))
       write (iulog,*) '!-------- > Nudge_AOD_Force_Opt.eq.0 :', Nudge_File_Present_Qaer(Nudge_ObsInd_AOD(1))  
     elseif(Nudge_AOD_Force_Opt.eq.1) then
       ! Verify that the CURR and NEXT analyses are available
       !-----------------------------------------------------
       Nudge_AOD_ON=(Nudge_File_Present_Qaer(Nudge_ObsInd_AOD(1)).and. &
                 Nudge_File_Present_Qaer(Nudge_ObsInd_AOD(2))      )
     else
       ! Verify that the ALL analyses are available
       !---------------------------------------------
       Nudge_AOD_ON=.true.
       do nn=1,Nudge_NumObs_AOD
         if(.not.Nudge_File_Present_Qaer(nn)) Nudge_AOD_ON=.false.
       end do
     endif
     if(.not.Nudge_AOD_ON) then
       if(masterproc) then
         write(iulog,*) 'NUDGING AODs: WARNING - AOD file NOT FOUND. Switching '
         write(iulog,*) 'NUDGING AODs:           nudging OFF to coast thru the gap. '
       endif
     endif
   else
     Nudge_AOD_ON=.false.
   endif   

  !  write(iulog,*) '!------------------------------------------------------------------------------'  
   write(iulog,*) '!-------- > Nudge_AOD_ON is on: ', Nudge_AOD_ON 	   
   write(iulog,*) '!-------- > NUDGING: Reading AOD file:',trim(Nudge_File_Qaer)
  !  write(iulog,*) '!------------------------------------------------------------------------------'      

   !-------------------------------------------------------
   ! HERE Implement time dependence of Nudging Coefs HERE
   !-------------------------------------------------------


   !---------------------------------------------------
   ! If Data arrays have changed update stepping arrays
   !---------------------------------------------------
   if((Before_End).and.((Update_Nudge).or.(Update_Model))) then

     ! Now Load the Target values for nudging tendencies
     !---------------------------------------------------
     
     if(Nudge_AOD_Force_Opt.eq.0) then	! --- Nudge_AOD_Force_Opt is indicated in the user_nl_cam
     
       ! -----------------------------------------> Target is OBS data at NEXT time
	   
	   ! ----> SARB
       do lchnk=begchunk,endchunk 
       
        ncol=phys_state(lchnk)%ncol   

! !         call calc_MidpointPressure(lchnk,ncol)
!         call relative_humidity(lchnk,ncol,pver,phys_state(lchnk)%q(:,:,1),phys_state(lchnk)%t(:,:),phys_state(lchnk)%pmid(:,:))   ! ---- works
!         call water_uptake_mam4(ncol,pver,lchnk)    ! -- commented by SK's suggestion 06/18/25     

!     	  ! --- call aerosol_sw_properties	   
! 	      call aerosol_sw_properties(lchnk,ncol,pver,ncoef,phys_state(lchnk)) 

        ! --- qaer_target_mam4 : this subroutine computes Target_Qaer = AOD_ratio * Model_Qaer  
        call qaer_target_mam4(lchnk,ncol)

       end do

                     
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
                             
         Target_AOD(:ncol,lchnk)=Nobs_AOD(:ncol,lchnk,Nudge_ObsInd_AOD(1)) ! ----- Here is when the Target_AOD are the MODIS-VIIRS files                            
         
       end do

      write(iulog,*), '|~~~~~~~~~~~~~~~~~~> if(Nudge_AOD_Force_Opt.eq.0) then Target_AOD=Nobs_AOD'
       
     elseif(Nudge_AOD_Force_Opt.eq.1) then
     
       !--------------------------------------> Target is linear interpolation of OBS data CURRENT <--> NEXT time    

       call ESMF_TimeSet(Date1,YY=Year,MM=Month,DD=Day,S=Sec)
       call ESMF_TimeSet(Date2,YY=Nudge_AOD_Next_Year,MM=Nudge_AOD_Next_Month, &
                               DD=Nudge_AOD_Next_Day , S=Nudge_AOD_Next_Sec    )
       DateDiff =Date2-Date1
       call ESMF_TimeIntervalGet(DateDiff,S=DeltaT,rc=rc)
       Tfrac= float(DeltaT)/float(Nudge_AOD_Step)
       
	   ! ----> SARB - 03/16/23
       do lchnk=begchunk,endchunk 
       
        ncol=phys_state(lchnk)%ncol   

! !         call calc_MidpointPressure(lchnk,ncol)
!         call relative_humidity(lchnk,ncol,pver,phys_state(lchnk)%q(:,:,1),phys_state(lchnk)%t(:,:),phys_state(lchnk)%pmid(:,:))   ! ---- works
!         call water_uptake_mam4(ncol,pver,lchnk)    ! -- commented by SK's suggestion 06/18/25     

!     	  ! --- call aerosol_sw_properties	   
! 	      call aerosol_sw_properties(lchnk,ncol,pver,ncoef,phys_state(lchnk)) 

        ! --- qaer_target_mam4    
        call qaer_target_mam4(lchnk,ncol)

       end do

       
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol

          !--> SARB : this loop was introduced at the beginning of the project as a continuation to do same as in the original nudging.F90
! 		     Contrarily what it was expected, declaring Target_Qaer is filled and updated in qaer_target_mam4 subroutine, therefore it is unnecessary to create 
!           a interpolation loop for Nobs_Qaer because we do not have observation for the aerosols mass concentration in each step. Everything come from the total AOD from MODIS-VIIRS

          ! --- 03/20/2023
!          do naer=1,27
! 
!             Target_Qaer(:ncol,:pver,naer,lchnk)=(1._r8-Tfrac)*Nobs_Qaer(:ncol,:pver,naer,lchnk,Nudge_ObsInd_AOD(1)) &
!                                              +Tfrac *Nobs_Qaer(:ncol,:pver,naer,lchnk,Nudge_ObsInd_AOD(2))
!          end do                                           


         ! --- 04/07/2025 : SK's suggestion based on what was done in repvious loop: 06/12/25 (re-run it)
		 do j=1,27
 
            Target_Qaer(:ncol,:pver,j,lchnk)=(1._r8-Tfrac)*Nobs_Qaer(:ncol,:pver,j,lchnk,Nudge_ObsInd_AOD(1)) &
                                             +Tfrac *Nobs_Qaer(:ncol,:pver,j,lchnk,Nudge_ObsInd_AOD(2))
         end do          
         
         Target_AOD(:ncol,lchnk)=(1._r8-Tfrac)*Nobs_AOD(:ncol,lchnk,Nudge_ObsInd_AOD(1)) &
                                           +Tfrac *Nobs_AOD(:ncol,lchnk,Nudge_ObsInd_AOD(2))
       end do

      write(iulog,*), '|~~~~~~~~~~~~~~~~~~> if(Nudge_AOD_Force_Opt.eq.1) then  Target_AOD is linear interpolation of OBS data CURR<-->NEXT time '
             
     else
       write(iulog,*) 'NUDGING AODs: Unknown Nudge_AOD_Force_Opt=',Nudge_AOD_Force_Opt
       call endrun('nudging_timestep_init_AOD:: ERROR unknown Nudging_Force_Opt')
     endif

     ! Now load Dry Static Energy values for Target
     !---------------------------------------------
!      if(Nudge_TSmode.eq.0) then
!        ! DSE tendencies from Temperature only
!        !---------------------------------------
!        do lchnk=begchunk,endchunk
!          ncol=phys_state(lchnk)%ncol
!          Target_S(:ncol,:pver,lchnk)=cpair*Target_T(:ncol,:pver,lchnk)
!        end do
!      elseif(Nudge_TSmode.eq.1) then
!        ! Caluculate DSE tendencies from Temperature, Water Vapor, and Surface Pressure
!        !------------------------------------------------------------------------------
!        do lchnk=begchunk,endchunk
!          ncol=phys_state(lchnk)%ncol
!          call calc_DryStaticEnergy(Target_T(:,:,lchnk), Target_Q(:,:,lchnk), &
!                                  phys_state(lchnk)%phis, Target_PS(:,lchnk), &
!                                                   Target_S(:,:,lchnk), ncol)
!        end do
!      endif

     ! Set Tscale for the specified Forcing Option 
     !-----------------------------------------------
     if(nudge_timescale_opt_AOD.eq.0) then
       Tscale=1._r8
     elseif(nudge_timescale_opt_AOD.eq.1) then
       call ESMF_TimeSet(Date1,YY=Year,MM=Month,DD=Day,S=Sec)
       call ESMF_TimeSet(Date2,YY=Nudge_AOD_Next_Year,MM=Nudge_AOD_Next_Month, &
                               DD=Nudge_AOD_Next_Day , S=Nudge_AOD_Next_Sec    )
       DateDiff =Date2-Date1
       call ESMF_TimeIntervalGet(DateDiff,S=DeltaT,rc=rc)
       Tscale=float(Nudge_AOD_Step)/float(DeltaT)
     else
       write(iulog,*) 'NUDGING AODs: Unknown nudge_timescale_opt_AOD=',nudge_timescale_opt_AOD
       call endrun('nudging_timestep_init_AOD:: ERROR unknown Nudging_TimeScale_Opt')
     endif

     ! Update the nudging tendencies
     !--------------------------------
     do lchnk=begchunk,endchunk
       ncol=phys_state(lchnk)%ncol

       !--> SARB

       do kaer=1,27
       
          Nudge_Qaer_step(:ncol,:pver,kaer,lchnk)=(  Target_Qaer(:ncol,:pver,kaer,lchnk)      &
                                            -Model_Qaer(:ncol,:pver,kaer,lchnk))     &
                                         *Tscale*Nudge_Qaer_tau(:ncol,:pver,kaer,lchnk)

		! ---- SK-FR suggestion 03-28-25 run and compare to control
!           Nudge_Qaer_step(:ncol,:pver,kaer,lchnk)=(  Model_Qaer(:ncol,:pver,kaer,lchnk)      &
!                                             -Model_Qaer(:ncol,:pver,kaer,lchnk))     &
!                                          *Tscale*Nudge_Qaer_tau(:ncol,:pver,kaer,lchnk)

       end do

                                             
     end do

     !******************
     ! DIAG
     !******************
!    if(masterproc) then
!      write(iulog,*) 'PFC: Target_T(1,:pver,begchunk)=',Target_T(1,:pver,begchunk)  
!      write(iulog,*) 'PFC:  Model_T(1,:pver,begchunk)=',Model_T(1,:pver,begchunk)
!      write(iulog,*) 'PFC: Target_S(1,:pver,begchunk)=',Target_S(1,:pver,begchunk)  
!      write(iulog,*) 'PFC:  Model_S(1,:pver,begchunk)=',Model_S(1,:pver,begchunk)
!      write(iulog,*) 'PFC:      Target_PS(1,begchunk)=',Target_PS(1,begchunk)  
!      write(iulog,*) 'PFC:       Model_PS(1,begchunk)=',Model_PS(1,begchunk)
!      write(iulog,*) 'PFC: Nudge_Sstep(1,:pver,begchunk)=',Nudge_Sstep(1,:pver,begchunk)
!      write(iulog,*) 'PFC: Nudge_Xstep arrays updated:'
!    endif
   endif ! ((Before_End).and.((Update_Nudge).or.(Update_Model))) then

   ! End Routine
   !------------
   return
  end subroutine ! nudging_timestep_init_AOD
  !================================================================
  
  !================================================================
  subroutine nudging_timestep_tend_AOD(phys_state,phys_tend)
   ! 
   ! NUDGING_TIMESTEP_TEND: 
   !                If Nudging is ON, return the Nudging contributions 
   !                to forcing using the current contents of the Nudge 
   !                arrays. Send output to the cam history module as well.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state,physics_ptend,physics_ptend_init
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use cam_history  ,only: outfld,  horiz_only

   ! Arguments
   !-------------
   type(physics_state), intent(inout) :: phys_state
   type(physics_ptend), intent(out):: phys_tend

   ! Local values
   !--------------------
   integer indw,ncol,lchnk, i,j,k
   logical lq(pcnst)

   ! --- advected constituents
   integer so4_a1,pom_a1,soa1_a1,soa2_a1,soa3_a1,soa4_a1,soa5_a1,bc_a1,dst_a1,ncl_a1,num_a1
   integer so4_a2,dst_a2,soa1_a2,soa2_a2,soa3_a2,soa4_a2,soa5_a2,ncl_a2,num_a2
   integer dst_a3,ncl_a3,so4_a3,num_a3
   integer pom_a4,bc_a4,num_a4
   integer mode
   integer, parameter :: nmde = 4
   integer kaer   
   
   ! -  character array constructor that creates an array of character strings called idxaer_name	
   character(len=7), parameter :: idxaer_name(27) = &
	 ['so4_a1 ','pom_a1 ','soa1_a1','soa2_a1','soa3_a1','soa4_a1','soa5_a1', &
	 'bc_a1  ','dst_a1 ','ncl_a1 ','num_a1 ', &
	 'so4_a2 ','dst_a2 ','soa1_a2','soa2_a2','soa3_a2','soa4_a2','soa5_a2', &
	 'ncl_a2 ','num_a2 ', &
	 'dst_a3 ','ncl_a3 ','so4_a3 ','num_a3 ', &
	 'pom_a4 ','bc_a4  ','num_a4 ']
   
   call cnst_get_ind('Q',indw)
   lq(:)   =.false.
   lq(indw)=.true. !----> commented 06/05/24
!    call physics_ptend_init(phys_tend,phys_state%psetcols,'nudging',lu=.true.,lv=.true.,ls=.true.,lq=lq)

     ! --- aerosols mode 1
     call cnst_get_ind('so4_a1',so4_a1)     
     call cnst_get_ind('pom_a1',pom_a1)     
     call cnst_get_ind('soa1_a1',soa1_a1)     
     call cnst_get_ind('soa2_a1',soa2_a1)     
     call cnst_get_ind('soa3_a1',soa3_a1)     
     call cnst_get_ind('soa4_a1',soa4_a1)     
     call cnst_get_ind('soa5_a1',soa5_a1)     
     call cnst_get_ind('bc_a1',bc_a1)     
     call cnst_get_ind('dst_a1',dst_a1)     
     call cnst_get_ind('ncl_a1',ncl_a1)     
     call cnst_get_ind('num_a1',num_a1)     
     
     ! --- aerosols mode 2     
     call cnst_get_ind('so4_a2',so4_a2)     
     call cnst_get_ind('soa1_a2',soa1_a2)     
     call cnst_get_ind('soa2_a2',soa2_a2)     
     call cnst_get_ind('soa3_a2',soa3_a2)     
     call cnst_get_ind('soa4_a2',soa4_a2)     
     call cnst_get_ind('soa5_a2',soa5_a2)     
     call cnst_get_ind('dst_a2',dst_a2)     
     call cnst_get_ind('ncl_a2',ncl_a2)     
     call cnst_get_ind('num_a2',num_a2)     

     ! --- aerosols mode 3       
     call cnst_get_ind('so4_a3',so4_a3)     
     call cnst_get_ind('dst_a3',dst_a3)     
     call cnst_get_ind('ncl_a3',ncl_a3)     
     call cnst_get_ind('num_a3',num_a3)     

     ! --- aerosols mode 4       
     call cnst_get_ind('pom_a4',pom_a4)     
     call cnst_get_ind('bc_a4',bc_a4)     
     call cnst_get_ind('num_a4',num_a4)     
     
     ! --- array for all the advected aerosols indexes in q for the MAM4

     idxaer =  [so4_a1,pom_a1,soa1_a1,soa2_a1,soa3_a1,soa4_a1,soa5_a1,bc_a1,dst_a1,ncl_a1,num_a1, &
                so4_a2,dst_a2,soa1_a2,soa2_a2,soa3_a2,soa4_a2,soa5_a2,ncl_a2,num_a2, &
                dst_a3,ncl_a3,so4_a3,num_a3, &
                pom_a4,bc_a4,num_a4]
                
   
   ! ------ added this similar to Q -- above
   
   if(Nudge_AOD_ON) then
    do  kaer = 1,27            
	   lq(idxaer(kaer))=.true.
    end do
   endif
    
	call physics_ptend_init(phys_tend,phys_state%psetcols,'nudging',lu=.true.,lv=.true.,ls=.true.,lq=lq)   
                	
   !-------------------------------> SARB
   if(Nudge_AOD_ON) then     
     lchnk=phys_state%lchnk
     ncol =phys_state%ncol
   
	 write (iulog,*) '!-------- >  CAM6 is using nudigng_AODs and is forcing AODs: ', trim(Nudge_File_Qaer)  
	 
     do kaer=1,27
     	phys_tend%q(:ncol,:pver,idxaer(kaer))=Nudge_Qaer_step(:ncol,:pver,kaer,lchnk) ! --- this option is to NUDGE target    
!      	phys_tend%q(:ncol,:pver,8+idxaer(kaer))=Nudge_Qaer_step(:ncol,:pver,kaer,lchnk) ! --- this option is to NUDGE target    
     	
!      	phys_tend%q(:ncol,:pver,idxaer(kaer))=Target_Qaer(:ncol,:pver,kaer,lchnk) ! --- this option tendency forced to be target_Qaer
   
!      	phys_state%q(:ncol,:pver,idxaer(kaer))=Target_Qaer(:ncol,:pver,kaer,lchnk) ! --- this option state is FORCED to target
     	
!      	phys_state(lchnk)%q(:ncol,:pver,idxaer(kaer))=Target_Qaer(:ncol,:pver,kaer,lchnk) ! --- this option is to FORCE target: 07/03/24 -Suggested by Seiji & Fred
     end do

     !--> SARB

	 call outfld( 'Nudge_Qaer_tau',Nudge_Qaer_tau(:ncol,:pver,:naer,lchnk) ,pcols,lchnk)    ! -- -added for check	 	 
	 call outfld( 'Nudge_Qaer_step',Nudge_Qaer_step(:ncol,:pver,:naer,lchnk) ,pcols,lchnk)    

	! --- term that generates in our file tendencies for each aerosols mode in MAM4 : 
     do kaer=1,27  
        call outfld( 'Nudge_Qaer',phys_tend%q(1,1,idxaer(kaer)),pcols,lchnk) 
         
         ! --- print name index  
        !  write(iulog,'(A,A,A,I0)') 'Aerosol mode ', trim(idxaer_name(kaer)), &
        !                   ' index value in q-array is: ', idxaer(kaer)
     end do
     
     call outfld('Target_AOD',Target_AOD(:,lchnk),pcols,lchnk) ! ---- AODs from MODIS-VIIRS
     !call outfld('aod_weighted_mean',Target_AOD(:,lchnk),pcols,lchnk) ! AOD weighted mean              
     call outfld('Model_Qaer',Model_Qaer(:,:,:,lchnk),pcols,lchnk) 
     call outfld('Target_Qaer',Target_Qaer(:,:,:,lchnk),pcols,lchnk) 

	! -----> outfield all the Q for the aerosols mass mixing ratios --------------------------------------
     call outfld('q_so4_a1',Model_Qaer(:,:,so4_a1,lchnk),pcols,lchnk)      
     call outfld('q_pom_a1',Model_Qaer(:,:,pom_a1,lchnk),pcols,lchnk)      
     call outfld('q_soa1_a1',Model_Qaer(:,:,soa1_a1,lchnk),pcols,lchnk)      
     call outfld('q_soa2_a1',Model_Qaer(:,:,soa2_a1,lchnk),pcols,lchnk)      
     call outfld('q_soa3_a1',Model_Qaer(:,:,soa3_a1,lchnk),pcols,lchnk)      
     call outfld('q_soa4_a1',Model_Qaer(:,:,soa4_a1,lchnk),pcols,lchnk)      
     call outfld('q_soa5_a1',Model_Qaer(:,:,soa5_a1,lchnk),pcols,lchnk)          
     call outfld('q_bc_a1',Model_Qaer(:,:,bc_a1,lchnk),pcols,lchnk)      
     call outfld('q_dst_a1 ',Model_Qaer(:,:,dst_a1,lchnk),pcols,lchnk)     
     call outfld('q_ncl_a1 ',Model_Qaer(:,:,ncl_a1,lchnk),pcols,lchnk)     
     call outfld('q_num_a1',Model_Qaer(:,:,num_a1,lchnk),pcols,lchnk)      
     call outfld('q_so4_a2',Model_Qaer(:,:,so4_a2,lchnk),pcols,lchnk)      
     call outfld('q_dst_a2',Model_Qaer(:,:,dst_a2,lchnk),pcols,lchnk)      
     call outfld('q_soa1_a2',Model_Qaer(:,:,soa1_a2,lchnk),pcols,lchnk)      
     call outfld('q_soa2_a2',Model_Qaer(:,:,soa2_a2,lchnk),pcols,lchnk)      
     call outfld('q_soa3_a2',Model_Qaer(:,:,soa3_a2,lchnk),pcols,lchnk)      
     call outfld('q_soa4_a2',Model_Qaer(:,:,soa4_a2,lchnk),pcols,lchnk)      
     call outfld('q_soa5_a2',Model_Qaer(:,:,soa5_a2,lchnk),pcols,lchnk)      
     call outfld('q_ncl_a2',Model_Qaer(:,:,ncl_a2,lchnk),pcols,lchnk)      
     call outfld('q_num_a2',Model_Qaer(:,:,num_a2,lchnk),pcols,lchnk)      
     call outfld('q_dst_a3',Model_Qaer(:,:,dst_a3,lchnk),pcols,lchnk)      
     call outfld('q_ncl_a3',Model_Qaer(:,:,ncl_a3,lchnk),pcols,lchnk)      
     call outfld('q_so4_a3',Model_Qaer(:,:,so4_a3,lchnk),pcols,lchnk)      
     call outfld('q_num_a3',Model_Qaer(:,:,num_a3,lchnk),pcols,lchnk)       
     call outfld('q_pom_a4',Model_Qaer(:,:,pom_a4,lchnk),pcols,lchnk)      
     call outfld('q_bc_a4',Model_Qaer(:,:,bc_a4,lchnk),pcols,lchnk)      
     call outfld('q_num_a4',Model_Qaer(:,:,num_a4,lchnk),pcols,lchnk)     
   
	!------------------------------------------------------------------------------------------
	
! !      call outfld('! mlp_c',mlp_c(:,:,lchnk),pcols,lchnk) ! Mid-point pressure computed
!      call outfld('mlp',mlp(:,:,lchnk),pcols,lchnk)   ! Mid-point pressure  
!      call outfld('RH',RH(:,:,lchnk),pcols,lchnk)	          
     
!      call outfld('lnRH',lnRH(:,:,lchnk),pcols,lchnk)	               
     
!      call outfld('vmr',vmr(:,:,:,lchnk),pcols,lchnk)  ! Dry Volume Mixing Ratio per mode            
!      call outfld('B_mode',B_mode(:,:,:,lchnk),pcols,lchnk)
!      call outfld('rdry_mode',rdry_mode(:,:,:,lchnk),pcols,lchnk)  ! Dry Modal Radius   
!      call outfld('Brd',Brd(:,:,:,lchnk),pcols,lchnk)       
!      call outfld('rwet_mode',rwet_mode(:,:,:,lchnk),pcols,lchnk) ! Wet Modal Radius    
     
!      call outfld('vmr_wet',vmr_wet(:,:,:,lchnk),pcols,lchnk)   ! Wet Volume Mixing Ratio per mode    
!      call outfld('vmr_water',vmr_water(:,:,:,lchnk),pcols,lchnk) ! Water Volume Mixing Ratio per mode      

!      call outfld('wet_ref_real_index',wet_ref_real_index(:,:,:,lchnk),pcols,lchnk) ! Real Wet refractive index per mode      
!      call outfld('wet_ref_im_index',wet_ref_im_index(:,:,:,lchnk),pcols,lchnk) ! Imaginary Wet refractive index per mode      
                    
!      call outfld('Model_Qaer_Mode',Model_Qaer_Mode(:,:,:,lchnk),pcols,lchnk) 

      !---- aerosol_sw_properties
!       call outfld('cext',          cext(:,1,lchnk),    pcols, lchnk)      
!       call outfld('cabs',          cabs(:,1,lchnk),    pcols, lchnk)      
!       call outfld('casm',          casm(:,1,lchnk),    pcols, lchnk)      

      ! --- parameterized optical properties
!       call outfld('pext',          pext(:,lchnk),    pcols, lchnk)      
!       call outfld('pabs',          pabs(:,lchnk),    pcols, lchnk)      
!       call outfld('pasm',          pasm(:,lchnk),    pcols, lchnk)      
!       call outfld('palb',          palb(:,lchnk),    pcols, lchnk)           
      
      ! call outfld('dopaer',        dopaer(:,:,:,lchnk),    pcols, lchnk)      
      call outfld('AODVISdn_computed',  AODVISdn_computed(:,lchnk),    pcols, lchnk)            
      call outfld('AOD_ratio',     AOD_ratio(:,lchnk),    pcols, lchnk)       

   endif
   
   ! End Routine
   !------------
   return
  end subroutine ! nudging_timestep_tend_AOD
  !================================================================

  !=====================================================================================!
  !  																					!
  ! --- HERE IT WAS nudging_update_analyses_se AND nudging_update_analyses_eul          !
  !																						!
  !=====================================================================================!
  

  !================================================================
!   subroutine nudging_update_analyses_fv(anal_file,aer_file,phys_state)
!   subroutine nudging_update_analyses_fv(anal_file,aer_file)
  subroutine nudging_update_analyses_fv(aer_file)
   ! 
   ! NUDGING_UPDATE_ANALYSES_FV: 
   !                 Open the given AOD data file, read in 
   !                 aod_weighted_mean values and then distribute
   !                 the values to all of the chunks.
   !
   ! ---> This version is intended to includes MODIS-VIIRS AODs
   !===============================================================
   use ppgrid ,only: pver,begchunk,endchunk
   use netcdf
   use physics_types,only: physics_state 
   
   ! Arguments
   !-------------
   character(len=*),intent(in) :: aer_file

   ! Local values
   !-------------
   integer lev
   integer nlon,nlat,plev,istat
   integer ncid,varid
   integer ilat,ilon,ilev
   real(r8) Xanal(Nudge_nlon,Nudge_nlat,Nudge_nlev)
   real(r8) Qanal(Nudge_nlon,Nudge_nlat)
   real(r8) PSanal(Nudge_nlon,Nudge_nlat)
   real(r8) Lat_anal(Nudge_nlat)
   real(r8) Lon_anal(Nudge_nlon)
   real(r8) Xtrans(Nudge_nlon,Nudge_nlev,Nudge_nlat)
   real(r8) Qtrans(Nudge_nlon,Nudge_nlat)
   integer  nn,Nindex
!    integer lchnk,ncol, l
!    type(physics_state),intent(in):: phys_state(begchunk:endchunk)

   ! Rotate Nudge_ObsInd_AOD() indices, then check the existence of the analyses 
   ! file; broadcast the updated indices and file status to all the other MPI nodes. 
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     Nindex=Nudge_ObsInd_AOD(Nudge_NumObs_AOD)
     do nn=Nudge_NumObs_AOD,2,-1
       Nudge_ObsInd_AOD(nn)=Nudge_ObsInd_AOD(nn-1)
     end do
     Nudge_ObsInd_AOD(1)=Nindex

     inquire(FILE=trim(aer_file),EXIST=Nudge_File_Present_Qaer(Nudge_ObsInd_AOD(1)))
     write(iulog,*)'NUDGING AODs: Nudge_ObsInd_AOD=',Nudge_ObsInd_AOD
     write(iulog,*)'NUDGING AODs: Nudge_File_Present_Qaer=',Nudge_File_Present_Qaer

   endif
#ifdef SPMD

   call mpibcast(Nudge_File_Present_Qaer, Nudge_NumObs_AOD, mpilog, 0, mpicom)
   call mpibcast(Nudge_ObsInd_AOD      , Nudge_NumObs_AOD, mpiint, 0, mpicom)
#endif
   if(.not.Nudge_File_Present_Qaer(Nudge_ObsInd_AOD(1))) return



   ! masterporc does all of the work here
   !-----------------------------------------
   ! ---> SARB 
   ! --- This section of code read and creates a variable for MODIS AODs
     
     if (masterproc) then
    
     ! Open the given file
     !-----------------------
    print*,'Aerosol file ',trim(aer_file)
     istat=nf90_open(trim(aer_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(aer_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif


     ! Read in Dimensions
     !--------------------
     istat=nf90_inq_dimid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlon)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_dimid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlat)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif


     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Lon_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_varid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Lat_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

!     if((Nudge_nlon.ne.nlon).or.(Nudge_nlat.ne.nlat).or.(plev.ne.pver)) then
     if((Nudge_nlon.ne.nlon).or.(Nudge_nlat.ne.nlat)) then
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: nlon=',nlon,' Nudge_nlon=',Nudge_nlon
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: nlat=',nlat,' Nudge_nlat=',Nudge_nlat
!      write(iulog,*) 'ERROR: nudging_update_analyses_fv: plev=',plev,' pver=',pver
      call endrun('nudging_update_analyses_fv: analyses dimension mismatch')
     endif


	! NOTE here it reads the variable in the netCDF file for the AOD
	
!      istat=nf90_inq_varid(ncid,'aod_weighted_mean_modis',varid)
    istat=nf90_inq_varid(ncid,'aod_weighted_mean_modis_old',varid)     
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_get_var(ncid,varid,Qanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     do ilat=1,nlat
     do ilon=1,nlon
!       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)

     Qtrans(ilon,ilat)=Qanal(ilon,ilat)       

!     IF ((Qtrans(ilon,ilat) > 0) .AND. (Qtrans(ilon,ilat) <= 1)) THEN
!         print*,'Reading MODIS file: ',aer_file
!	 print*,'MODIS nudge read: ',Qtrans(ilon,ilat),ilon,ilat
!         print*,'***************'
!     ENDIF

     end do
     end do

  !    do ilat=1,nlat
!       print*,'Aerosol nudging latitude = ',Lat_anal(ilat),ilat
!      enddo
!      do ilon=1,nlon
!       print*,'Aerosol nudging longitude = ',Lon_anal(ilon),ilon
!      enddo

     ! Close the analyses file
     !-----------------------
     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

    endif ! (masterproc) 

    call scatter_field_to_chunk(1,1,1,Nudge_nlon,Qtrans,   &
                                Nobs_AOD(1,begchunk,Nudge_ObsInd_AOD(1)))   
                                
!     print*, '------> Begchunk in Nobs_AOD after using scatter_field_to_chunk : ', begchunk                             
                                
   
   ! End Routine
   !------------
   return
  end subroutine ! nudging_update_analyses_fv
  !================================================================

  !================================================================
  subroutine nudging_set_profile(rlat,rlon,Nudge_prof,Wprof,nlev)
   ! 
   ! NUDGING_SET_PROFILE: for the given lat,lon, and Nudging_prof, set
   !                      the verical profile of window coeffcients.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on nudging strength.
   !===============================================================

   ! Arguments
   !--------------
   integer  nlev,Nudge_prof
   real(r8) rlat,rlon
   real(r8) Wprof(nlev)

   ! Local values
   !----------------
   integer  ilev
   real(r8) Hcoef,latx,lonx,Vmax,Vmin
   real(r8) lon_lo,lon_hi,lat_lo,lat_hi,lev_lo,lev_hi

   !---------------
   ! set coeffcient
   !---------------
   if(Nudge_prof.eq.0) then
     ! No Nudging
     !-------------
     Wprof(:)=0.0_r8
   elseif(Nudge_prof.eq.1) then
     ! Uniform Nudging
     !-----------------
     Wprof(:)=1.0_r8
   elseif(Nudge_prof.eq.2) then
     ! Localized Nudging with specified Heaviside window function
     !------------------------------------------------------------
     if(Nudge_AOD_Hwin_max.le.Nudge_AOD_Hwin_min) then
       ! For a constant Horizontal window function, 
       ! just set Hcoef to the maximum of Hlo/Hhi.
       !--------------------------------------------
       Hcoef=max(Nudge_AOD_Hwin_lo,Nudge_AOD_Hwin_hi)
     else
       ! get lat/lon relative to window center
       !------------------------------------------
       latx=rlat-Nudge_AOD_Hwin_lat0
       lonx=rlon-Nudge_AOD_Hwin_lon0
       if(lonx.gt. 180._r8) lonx=lonx-360._r8
       if(lonx.le.-180._r8) lonx=lonx+360._r8

       ! Calcualte RAW window value
       !-------------------------------
       lon_lo=(Nudge_AOD_Hwin_lonWidthH+lonx)/Nudge_AOD_Hwin_lonDelta
       lon_hi=(Nudge_AOD_Hwin_lonWidthH-lonx)/Nudge_AOD_Hwin_lonDelta
       lat_lo=(Nudge_AOD_Hwin_latWidthH+latx)/Nudge_AOD_Hwin_latDelta
       lat_hi=(Nudge_AOD_Hwin_latWidthH-latx)/Nudge_AOD_Hwin_latDelta
       Hcoef=((1._r8+tanh(lon_lo))/2._r8)*((1._r8+tanh(lon_hi))/2._r8) &
            *((1._r8+tanh(lat_lo))/2._r8)*((1._r8+tanh(lat_hi))/2._r8)

       ! Scale the horizontal window coef for specfied range of values.
       !--------------------------------------------------------
       Hcoef=(Hcoef-Nudge_AOD_Hwin_min)/(Nudge_AOD_Hwin_max-Nudge_AOD_Hwin_min)
       Hcoef=(1._r8-Hcoef)*Nudge_AOD_Hwin_lo + Hcoef*Nudge_AOD_Hwin_hi
     endif

     ! Load the RAW vertical window
     !------------------------------
     do ilev=1,nlev
       lev_lo=(float(ilev)-Nudge_AOD_Vwin_Lindex)/Nudge_AOD_Vwin_Ldelta
       lev_hi=(Nudge_AOD_Vwin_Hindex-float(ilev))/Nudge_AOD_Vwin_Hdelta
       Wprof(ilev)=((1._r8+tanh(lev_lo))/2._r8)*((1._r8+tanh(lev_hi))/2._r8)
     end do 

     ! Scale the Window function to span the values between Vlo and Vhi:
     !-----------------------------------------------------------------
     Vmax=maxval(Wprof)
     Vmin=minval(Wprof)
     if((Vmax.le.Vmin).or.((Nudge_AOD_Vwin_Hindex.ge.(nlev+1)).and. &
                           (Nudge_AOD_Vwin_Lindex.le. 0      )     )) then
       ! For a constant Vertical window function, 
       ! load maximum of Vlo/Vhi into Wprof()
       !--------------------------------------------
       Vmax=max(Nudge_AOD_Vwin_lo,Nudge_AOD_Vwin_hi)
       Wprof(:)=Vmax
     else
       ! Scale the RAW vertical window for specfied range of values.
       !--------------------------------------------------------
       Wprof(:)=(Wprof(:)-Vmin)/(Vmax-Vmin)
       Wprof(:)=Nudge_AOD_Vwin_lo + Wprof(:)*(Nudge_AOD_Vwin_hi-Nudge_AOD_Vwin_lo)
     endif

     ! The desired result is the product of the vertical profile 
     ! and the horizontal window coeffcient.
     !----------------------------------------------------
     Wprof(:)=Hcoef*Wprof(:)
   else
     call endrun('nudging_set_profile:: Unknown Nudge_prof value')
   endif

   ! End Routine
   !------------
   return
  end subroutine ! nudging_set_profile
  !================================================================


  !================================================================
  real(r8) function nudging_set_PSprofile(rlat,rlon,Nudge_PSprof)
   ! 
   ! NUDGING_SET_PSPROFILE: for the given lat and lon set the surface
   !                      pressure profile value for the specified index.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on nudging strength.
   !===============================================================

   ! Arguments
   !--------------
   real(r8) rlat,rlon
   integer  Nudge_PSprof

   ! Local values
   !----------------

   !---------------
   ! set coeffcient
   !---------------
   if(Nudge_PSprof.eq.0) then
     ! No Nudging
     !-------------
     nudging_set_PSprofile=0.0_r8
   elseif(Nudge_PSprof.eq.1) then
     ! Uniform Nudging
     !-----------------
     nudging_set_PSprofile=1.0_r8
   else
     call endrun('nudging_set_PSprofile:: Unknown Nudge_prof value')
   endif

   ! End Routine
   !------------
   return
  end function ! nudging_set_PSprofile
  !================================================================


  !================================================================
  subroutine calc_DryStaticEnergy(t, q, phis, ps, dse, ncol)
   ! 
   ! calc_DryStaticEnergy: Given the temperature, specific humidity, surface pressure,
   !                       and surface geopotential for a chunk containing 'ncol' columns, 
   !                       calculate and return the corresponding dry static energy values.
   !--------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pver, pverp
   use dycore,       only: dycore_is
   use hycoef,       only: hyai, hybi, ps0, hyam, hybm
   use physconst,    only: zvir, gravit, cpair, rair
   !
   ! Input/Output arguments
   !-----------------------
   integer , intent(in) :: ncol      ! Number of columns in chunk
   real(r8), intent(in) :: t(:,:)    ! (pcols,pver) - temperature
   real(r8), intent(in) :: q(:,:)    ! (pcols,pver) - specific humidity
   real(r8), intent(in) :: ps(:)     ! (pcols)      - surface pressure 
   real(r8), intent(in) :: phis(:)   ! (pcols)      - surface geopotential
   real(r8), intent(out):: dse(:,:)  ! (pcols,pver)  - dry static energy
   !--->SARB
    !real(r8), intent(in) :: Qaer(:,:,pcnst)    ! (pcols,pver,pcnst) - aerosol constituent
   !
   ! Local variables
   !------------------
   logical  :: fvdyn                 ! finite volume dynamics
   integer  :: ii,kk                 ! Lon, level, level indices
   real(r8) :: tvfac                 ! Virtual temperature factor
   real(r8) :: hkk(ncol)             ! diagonal element of hydrostatic matrix
   real(r8) :: hkl(ncol)             ! off-diagonal element
   real(r8) :: pint(ncol,pverp)      ! Interface pressures
   real(r8) :: pmid(ncol,pver )      ! Midpoint pressures
   real(r8) :: zi(ncol,pverp)        ! Height above surface at interfaces
   real(r8) :: zm(ncol,pver )        ! Geopotential height at mid level

   ! Set dynamics flag
   !-------------------
   fvdyn = dycore_is ('LR')

   ! Load Pressure values and midpoint pressures 
   !----------------------------------------------
   do kk=1,pverp
     do ii=1,ncol
       pint(ii,kk)=(hyai(kk)*ps0)+(hybi(kk)*ps(ii))
     end do
   end do
   do kk=1,pver
     do ii=1,ncol
       pmid(ii,kk)=(hyam(kk)*ps0)+(hybm(kk)*ps(ii))
     end do
   end do

   ! The surface height is zero by definition.
   !-------------------------------------------
   do ii = 1,ncol
     zi(ii,pverp) = 0.0_r8
   end do

   ! Compute the dry static energy, zi, zm from bottom up
   ! Note, zi(i,k) is the interface above zm(i,k)
   !---------------------------------------------------------
   do kk=pver,1,-1

     ! First set hydrostatic elements consistent with dynamics
     !--------------------------------------------------------
     if(fvdyn) then
       do ii=1,ncol
         hkl(ii)=log(pint(ii,kk+1))-log(pint(ii,kk))
         hkk(ii)=1._r8-(hkl(ii)*pint(ii,kk)/(pint(ii,kk+1)-pint(ii,kk)))
       end do
     else
       do ii=1,ncol
         hkl(ii)=(pint(ii,kk+1)-pint(ii,kk))/pmid(ii,kk)
         hkk(ii)=0.5_r8*hkl(ii)
       end do
     endif

     ! Now compute zm, zi, and dse  (WACCM-X vars rairv/zairv/cpairv not used!)
     !------------------------------------------------------------------------
     do ii=1,ncol
       tvfac=t(ii,kk)*rair*(1._r8+(zvir*q(ii,kk)))/gravit
       zm (ii,kk)=zi(ii,kk+1) + (tvfac*hkk(ii))
       zi (ii,kk)=zi(ii,kk+1) + (tvfac*hkl(ii))
       dse(ii,kk)=(t(ii,kk)*cpair) + phis(ii) + (gravit*zm(ii,kk))
     end do

   end do ! kk=pver,1,-1

   ! End Routine
   !-----------
   return
  end subroutine calc_DryStaticEnergy
  !================================================================

  !==========================================================================!
  !
  !  BELOW THERE ARE ALL THE SUBROUTINES TO COMPUTE AODVIS AND NUDGE AODs    !
  !
  !==========================================================================!  
  
!   !================================================================
!   subroutine calc_MidpointPressure(lchnk,ncol)
! 
!    use ppgrid,       only: pver, pverp
!    use hycoef,       only: ps0, hyam, hybm
!    
! !    Input/Output arguments
! !    -----------------------
!    integer, intent(in) :: ncol  ! number of columns in chunk      
!    integer, intent(in) :: lchnk  ! 
!    
! !    Local variables
! !    ------------------
!    integer  :: ii,kk                 ! Lon, level, level indices
! 
! !    midpoint pressures 
! !    ----------------------------------------------
!    do kk=1,pver
!      do ii=1,ncol
!        mlp_c(ii,kk,lchnk)=(hyam(kk)*ps0)+(hybm(kk)*Model_PS(ii,lchnk))
!      end do
!    end do
! 
!    end subroutine ! calc_MidpointPressure
! !   ================================================================ 

!   Subroutines deleted
!   relative_humidty
!   Create_Model_Qaer_Mode
!   water_uptake_mam4
!   aerosol_sw_properties
              
  !================================================================
  subroutine qaer_target_mam4(lchnk,ncol)  

	! Adjusts MAM4 aerosol mass mixing ratios by scaling with AOD ratio (Target/AODVISdn_computed). 
	! Computes scaling factor from satellite-observed vs. model-calculated AOD, then applies 
	! uniform scaling to all aerosol species and vertical levels to match target AOD values 
	! for data assimilation or model constraint purposes.

   use ppgrid       ,only: pver
   use modal_aer_opt, only: get_aodvisdn  

   integer, intent(in)           :: ncol,lchnk 
   
   integer :: l, kcol
   integer, parameter :: nmde = 4
   real(r8) :: AODVISdn_computed(ncol,lchnk)
     
   call get_aodvisdn(AODVISdn_computed,  ncol, lchnk)
      
  ! Compute the ratio between MODIS AOD and computed from RTM
    
    do kcol =1,ncol
        if 	((AODVISdn_computed(kcol,lchnk) .le. 0.01) .or. (Target_AOD(kcol,lchnk) .lt. 0)) then          	  
          AOD_ratio(kcol,lchnk) = 1
        else      	  	
          AOD_ratio(kcol,lchnk) = Target_AOD(kcol,lchnk)/AODVISdn_computed(kcol,lchnk)
        end if   		
    end do
    
    do kcol = 1,ncol       
        do l = 1, naer
          Target_Qaer(kcol,:pver,l,lchnk)=AOD_ratio(kcol,lchnk)*Model_Qaer(kcol,:pver,l,lchnk)         
        end do
    end do
   
   return
  end subroutine ! qaer_target_mam4
 
  !================================================================  
  
end module nudging_AODs