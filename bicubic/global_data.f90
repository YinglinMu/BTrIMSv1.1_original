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
INTEGER, PARAMETER :: nparcels = 200   !set the number of parcels to release per day per grid cell if it rains
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