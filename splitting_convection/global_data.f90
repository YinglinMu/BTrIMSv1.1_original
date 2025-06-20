!%%%%%%% AUTHORSHIP %%%%%%%
! This program was written by Jason P. Evans (UNSW) and modified by Chiara M. Holgate (ANU) and Yinglin Mu(UNSW)
! Please seek permission before publishing all or part of this code, from jason.evans@unsw.edu.au and/or chiara.holgate@anu.edu.au and/or yinglin.mu@unsw.edu.au.


!%%%%%%% PURPOSE %%%%%%%
! This program calculates the back trajectories of water vapor using a method based on Dirmeyer & Brubaker, 1999, "Contrasting evaporative moisture sources during the drought of 1988 and the flood of 1993", Journal of Geophysical Research, 104 D16 pg 19,383-19,397.
!
!%%%%%%% INPUT %%%%%%%
! Input data are taken from ERA5, in this model, each independent file for each month and each variable. 


!%%Model expects:%%
! % ERA5 reanalysis data, hourly.
! % u,v: horizontal wind [m/s], 4d.
! % w: vertical wind [Pa/s], 4d.
! % t: actual temperature [K], 4d.
! % q: specific humidity [kg/kg], 4d.
! % ciwc: cloud ice water content [kg/kg], 4d.
! % clwc: cloud liquid water content [kg/kg], 4d.
! % cswc: cloud snow water content [kg/kg], 4d.
! % crwc: cloud rain water content [kg/kg], 4d.
! % tp: total precipitation,[m], 3d
! % e: evapotranspiration [m], 3d.
! % sp: surface pressure [Pa], 3d.
! % blh: boundary layer height [m], 3d.
! % cp: convective precipitation,[m], 3d



!%%%%%%% OUTPUT %%%%%%%
! The program outputs a daily 2d grid of water vapour contribution of each grid cell to each precipitation cell. 
! It also outputs the total rainfall depth in each cell where it rained, and the coordinates of each rain grid cell.
! 
! The output file is dimensioned (rec,lat,lon) where rec is the number of grid cells where it rained that day, 
! lat & lon are the sizes of the domain in the lat and lon directions. Attributes include location and time and amount 
! of precip in the pixel of interest (y,x,day,month,year,precip) each of which is a 1D array of length rec.
!

!%%%%%%% EXPLANATION OF TIME-RELATED VARIABLES AND OTHER USER SET VARIABLES%%%%%%%
! "set"   = user defined values
! "calcd" = calculated within program


! totdays =  number of days to run simulation forward, based on set start/end dates(calcd)
! totbtadays = number of days to back-track for (set)
! tstep = number of minutes for back-track time step (set)
! daytsteps = 1440/tstep = 48 = number of simulation time steps in a day (calcd)
! totsteps = daytsteps*(totbtadays+1) = total number of simulation time steps over period (calcd)
! datatstep = 180 = input file time step in minutes (set)
! datadaysteps = 1440/datatstep = 8 = number of input file time steps in a day (calcd)
! indatatsteps = datatstep/tstep = number of simulation time steps per input time step  (calcd)
! datatotsteps = (datadaysteps*(totbtadays+1)) + 1 = total number of input time steps over the back-track period (calcd)
!.............we want the event day + totbtadays before it + the time step after ! (ie 0 hour time step which is the last in each file)
! ttdataday = = ((tt-1)/indatatsteps) + 1 = position of parcel time step in input file (calcd)
! ttdata = datatotsteps - datadaysteps - 1 + ttdataday = input file time step from the beginning of the loaded files (calcd)



!%%%%%%% ASSUMPTIONS %%%%%%%
! All forms of water are included (vapour, rain, ice, snow, cloud water) in the parcel's mixing ratio. 
! When releasing air parcels, time of parcel release is determined randomly through a precipitation-weighted time profile.
! Height of parcel release is determined randomly from a humidity-weighted vertical profile, i.e the vertical distribution of total water mass indicates where the rain forms.
! Air parcels move kinetically, i.e. forced by wind field in 3 dimensions.
! The PBL is split from the upper atmosphere, i.e. ET is well-mixed within PBL, 'ET-mixing' method is used, in the upper troposphere, moisture contribution is represented by change of water mass content in the parcel, WaterSip identification method is used.However, the PBL is not well-mixed, i.e. parcels are released at different heights within the PBL.
! However, when deep convection occurs, the ET is well-mixed in the whole air column, i.e. 'ET-mixing' is used even the air parcel is above PBl.

MODULE global_data

IMPLICIT NONE

SAVE

!
!*******user modified variables**********************
!

INTEGER :: sday,smon,syear    !start day for calculations
INTEGER :: edday,edmon,edyear !end day for calculations (Exclusive. Must be at least one day after start day)
INTEGER :: totdays
INTEGER, PARAMETER :: totbtadays = 15   !number of days of data to keep for bta; i.e. how far back in time to calc.
                                       !must be less than days you have input data for
INTEGER, PARAMETER :: tstep = 15   !number of minutes for back trajectory time step (simultion time step)
                      !must divide evenly into number of minutes in day 1440 and number of minutes in MM5 time step (here 60)
INTEGER, PARAMETER :: nparcels = 50   !set the number of parcels to release per day per grid cell if it rains
REAL, PARAMETER :: minpre = 1   !min daily precip to deal with (mm)

INTEGER, PARAMETER :: bdy = 6   !boundary layers to ignore; trajectories will be tracked to this boundary

CHARACTER(LEN=60), PARAMETER :: diri_era5 = "/scratch/w40/ym7079/"  !directory of mask file for the study region 

CHARACTER(LEN=100) :: diro  

CHARACTER(LEN=100), PARAMETER :: dirdata_era5 = "/g/data/rt52/era5/"

INTEGER, PARAMETER :: numthreads = 104   !set the number of parallel openmp threads

LOGICAL, PARAMETER :: peak = .FALSE. !does the daylist indicate storm peaks (TRUE) or whole days (FALSE)

LOGICAL, PARAMETER :: wshed = .TRUE. !only calculate trajectories for watershed

CHARACTER(LEN=50), PARAMETER :: fwshed_era5 = "Scotland_mask.nc" ! mask file for the study region 
                                !set to "" if no watershed
                                !0 outside watershed, >0 inside

REAL, PARAMETER :: min_del_q = 0.0001    !the minimum change in parcel mixing ratio (kg/kg) to be considered "real"
REAL, PARAMETER :: delta_coord = 0.0001  ! 1/10000th degree - for floating point calculations

LOGICAL, PARAMETER :: eachParcel = .FALSE.   !output the data along the trajectory of each parcel
REAL, PARAMETER :: domain(6) = (/ -40.0, 70.0, 10.0, 120.0, 1.0, 1000.0 /) ! domains sorrounding the study area forloading input data, tracking moisture and outputing result. !!! Note that values in the domain array MUST match coord points in the data   

!****************************************************

INTEGER :: daytsteps,totsteps,indatatsteps,datadaysteps,datatotsteps
INTEGER :: dim_i,dim_j,dim_k,fdim_i,fdim_j,ssdim
INTEGER :: dim_i_start, dim_j_start, dim_k_start
INTEGER :: dim_i_end, dim_j_end, dim_k_end
INTEGER :: mon,year,dd,totpts
INTEGER :: day
! Additional variable to define the number of time intervals in the input data, now that the input data is monthly instead of daily
INTEGER :: datansteps


REAL, PARAMETER :: Lv = 2.25E6   !latent heat of vaporization of water (Jkg-1)
REAL, PARAMETER :: g = 9.8     !gravity (m.s-2)
REAL, PARAMETER :: P0 = 100000   !reference surface pressure (Pa)
REAL, PARAMETER :: Rd = 287.053   !ideal gas constant for dry air (J/kgK)
REAL, PARAMETER :: Cp = 1004.67   !heat capacity of air at constant pressure (J/kgK)
REAL, PARAMETER :: Rv = 461.5    !gas constant of water vapor (J/kgK)
REAL, PARAMETER :: Cl = 4400     !heat capacity of liquid water at ~-20C (J/kgK)
REAL, PARAMETER :: pi = 3.14159265
REAL, PARAMETER :: deg_dist = 111.   !average distance of 1 degree lat is assumed to be 111km
REAL, PARAMETER :: water_density = 1000 ! density of water (kg/m3)

END MODULE global_data