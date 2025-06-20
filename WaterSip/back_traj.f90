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
! nparcels = number of parcels to be tracked per day(set)



!%%%%%%% ASSUMPTIONS %%%%%%%
! All forms of water are included (vapour, rain, ice, snow, cloud water) in the parcel's mixing ratio. 
! When releasing air parcels, time of parcel release is determined randomly through a precipitation-weighted time profile.
! Height of parcel release is determined randomly from a humidity-weighted vertical profile, i.e the vertical distribution of total water mass indicates where the rain forms.
! Air parcels move kinetically, i.e. forced by wind field in 3 dimensions.
! The PBL is NOT split from the upper atmosphere, i.e. in the whole air column, moisture contribution is represented by change of water mass content in the parcel.


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PROGRAM back_traj

	USE netcdf
	USE util
	USE global_data
	USE bt_subs
	USE omp_lib
	USE input_data_handling_era5, ONLY: get_grid_data, get_data, get_watershed


	IMPLICIT NONE

	!
	!netcdf id variables
	!
	INTEGER :: status
	INTEGER :: outncid,wvcid,xlocid,ylocid,dayid,opreid
	!
	!data variables
	!
	INTEGER :: par_lev
	REAL :: ptop,delx,par_lat,par_lon,par_pres,par_q,new_par_q,end_precip
	INTEGER :: datatstep
	REAL,ALLOCATABLE,DIMENSION(:,:) :: lat2d,lon2d
	REAL,ALLOCATABLE,DIMENSION(:,:) :: terrain,WV_cont,WV_cont_day
	REAL,ALLOCATABLE,DIMENSION(:,:,:) :: precip
	REAL,ALLOCATABLE,DIMENSION(:,:,:) :: evap,tpw,pbl_hgt,pstar,psfc
	REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: u,v,w,temp,mix,pw,mixcld,mixtot,pres
	!REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: pot_temp !
	REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: unow,vnow,wnow
	REAL,ALLOCATABLE,DIMENSION(:,:,:) :: pres_then,tempnow
	!REAL,ALLOCATABLE,DIMENSION(:,:,:) :: pot_temp_then !
	REAL,ALLOCATABLE,DIMENSION(:,:) :: psfc_then,evapnow
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: pbl_lev

	INTEGER,ALLOCATABLE,DIMENSION(:) :: par_release
	INTEGER :: xx,yy,tt,nn,mm,npar,orec,x,y,ttdata,nnMM5,ttdataday
	INTEGER :: xx_omp,threadnum,torec
	REAL :: ttfac,nnfac,precip_here,qfac_evap,qfac_pbl,fac

	INTEGER,ALLOCATABLE,DIMENSION(:,:) :: wsmask

	INTEGER :: ssx,ssy

	!I want dd and totpts to persist inside and outside parallel regions and subroutines
	!so they have been added to the global_data module and declared threadprivate

	!for outputting the parcel stats for each parcel
	REAL,ALLOCATABLE,DIMENSION(:,:) :: parcel_stats


	!----------------------------------------------------------------
	! Retrieve simulation start and end dates, and output directory, from the command line input

	INTEGER :: num_args
	character(len=100), dimension(:), allocatable :: args
	num_args = command_argument_count()
	allocate(args(num_args))  ! I've omitted checking the return status of the allocation

	call get_command_argument(1,args(1))
	sday = string_to_int(args(1))
	call get_command_argument(2,args(2))
	smon = string_to_int(args(2))
	call get_command_argument(3,args(3))
	syear = string_to_int(args(3))
	call get_command_argument(4,args(4))
	edday = string_to_int(args(4))
	call get_command_argument(5,args(5))
	edmon = string_to_int(args(5))
	call get_command_argument(6,args(6))
	edyear = string_to_int(args(6))
	call get_command_argument(7,args(7))
	diro = args(7)

	print *,"Results saved to: ",diro

	! Find total number of days to run simulation for, given input start and end dates
	totdays=simlength(sday,smon,syear,edday,edmon,edyear)

	print *,"Total days to run analysis=",totdays
	print *,"Total number of back-track days=",totbtadays
	print *,"Number of parcels=",nparcels
	print *,"Simulation time step (mins)=",tstep

	!----------------------------------------------------------------
	! Get header info from first input file
	!!! Note that values in the extents array MUST match coord points in the data       
        !!! extents = (/ start_lat, end_lat, start_lon, end_lon, start_level, end_level /) levels in hPa
        !Pakistan case
    call get_grid_data(ptop, delx, datatstep, lat2d, lon2d, domain)
    
	!--------------------------------------------------------
    print *,"dim_j, dim_i, dim_k",dim_j,dim_i, dim_k
    !print *, 'lat2d(1,:)',lat2d(1,:)
    !print *, 'lon2d(:,1)',lon2d(:,1)

	!
	! Calculate the number of trajectory time steps in a day and in input file time step
	!
	daytsteps = 1440/tstep                ! number of sub-daily simulation time steps
	indatatsteps = datatstep/tstep          ! divide input file time step by the number of simulation time steps, as they may differ
	totsteps = daytsteps*(totbtadays+1)   ! total number of simulation data time steps to remember

	datadaysteps = 1440/datatstep           ! number of input file time steps in day
	datatotsteps = (datadaysteps*(totbtadays+1)) + 1 ! total number of input file time steps over the back-track period
    
    print *,'simulation time step [mins] (tstep): ',tstep
    print *,'input data timestep [mins] (datatstep): ',datatstep
    print *,'no. of time intervals per daily file (datadaysteps): ',datadaysteps
    print *,'no. of simulation timesteps per input file (daytsteps): ',daytsteps
    print *, 'no. of simulation timesteps per day (daytsteps): ', daytsteps
    print *,'no. of simulation timesteps per input file time interval (indatatsteps): ',indatatsteps
    print *,'total no. of back-track simulation timesteps to remember (totsteps): ',totsteps
    print *,'total no. of back-track input file time intervals (datatotsteps): ',datatotsteps
    print *, 'datansteps', datansteps
	!
	! Read in watershed mask
	!
	if (wshed) then
		call get_watershed(wsmask)
	end if

	! Total number of grid pts inside boundaries
	totpts = (dim_j-2)*(dim_i-2)

	! Set the number of threads to use in the parallel sections
	call OMP_SET_NUM_THREADS(numthreads)

				! Allocate the variable arrays
	ALLOCATE( precip(dim_j,dim_i,datadaysteps), &
		evap(dim_j,dim_i,datatotsteps),   &
		tpw(dim_j,dim_i,datatotsteps),   &
		pbl_hgt(dim_j,dim_i,datatotsteps),   &
		pbl_lev(dim_j,dim_i,datatotsteps),   &
		psfc(dim_j,dim_i,datatotsteps),   &
		u(dim_j,dim_i,dim_k,datatotsteps), &
		v(dim_j,dim_i,dim_k,datatotsteps), &
		w(dim_j,dim_i,dim_k,datatotsteps), &
		temp(dim_j,dim_i,dim_k,datatotsteps), &
		mix(dim_j,dim_i,dim_k,datatotsteps), &
		mixtot(dim_j,dim_i,dim_k,datatotsteps), &
		pw(dim_j,dim_i,dim_k,daytsteps+1), &
		mixcld(dim_j,dim_i,dim_k,datatotsteps), &
		pres(dim_j,dim_i,dim_k,datatotsteps), &
		STAT = status )

	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	! FOR EVERY DAY OF THE SIMULATION PERIOD
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do dd = 1, totdays

		orec = 0

		!Get date to open correct input file
		call day_month_year(dd)
		print *,"day,mon,year",day,mon,year

		! Create output file (will create empty file even if it didn't rain anywhere in the domain on that day)
		call new_out_file(outncid,wvcid,xlocid,ylocid,dayid,opreid,day,lat2d,lon2d)
		print *,'created new out file'

		! Get the variables required for back trajectory calculations for the current day
		call get_data(precip,evap,u,v,w,temp,mix,mixcld,mixtot,pres,pbl_hgt,psfc,tpw)

		print *,'got data...'

		! ERA5 gives total precipitation in m. Multiply by 1000 to get it in mm as model expects.
		precip = precip * 1000

  
		! Calculate the precipitable water accumulated from the ground up on day of interest (lat,lon,height,time). 
		! This is used to determine the parcel initial height.
		! Note that pw has an extra timestep in length, to allow the lin_interp_inMM5tsteps(pw) to interpolate between 2 values at the end.
		call calc_pw(mixtot(:,:,:,datatotsteps-datadaysteps:),pres(:,:,:,datatotsteps-datadaysteps:),psfc(:,:,datatotsteps-datadaysteps:),ptop,pw)

		! Calculate the total precipitable water (lat,lon,time).
		call calc_tpw(mixtot,pres,psfc,ptop,tpw)

		! Calculate the subsection x & y dimensions, based on the max distance a parcel can travel in the sim timestep
		ssdim = (ceiling((sqrt(maxval(u)**2+maxval(v)**2)*tstep*60)/delx) *2) + 1

		!loop over x and y grid points
		!

		!parallelize over the grid points
		!ie each grid point will be sent to a new thread
		!program will wait until all grid points are complete before
		!moving on to next day


		print *, 'Starting parallelisation'
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(pw,tpw,u,v,w,pres,temp,psfc,evap,precip,mix,mixtot,lat2d,lon2d,orec,outncid,wvcid,xlocid,ylocid,dayid,opreid,wsmask,daytsteps,totsteps,indatatsteps,datadaysteps,datatotsteps,dim_i,dim_j,dim_k,sday,smon,syear,mon,year,day,dd,totpts,ssdim)
		!allocate these arrays for each thread
		ALLOCATE( WV_cont(dim_j,dim_i),WV_cont_day(dim_j,dim_i), &
				unow(ssdim,ssdim,dim_k,2),vnow(ssdim,ssdim,dim_k,2), &
				par_release(daytsteps), &
				!pot_temp_then(ssdim,ssdim,dim_k), &
				pres_then(ssdim,ssdim,dim_k),wnow(ssdim,ssdim,dim_k,2), &
				psfc_then(ssdim,ssdim),tempnow(ssdim,ssdim,dim_k), &
				STAT = status)
        if (eachParcel) then
            ALLOCATE(parcel_stats(14,totsteps), STAT = status)
        end if
  
!$OMP DO &
!$OMP SCHEDULE (DYNAMIC)
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		! FOR EVERY POINT IN THE GRID
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		do xx_omp = 0, totpts-1

			xx = 2 + AINT(xx_omp*1./(dim_i-2))
			yy = 2 + (xx_omp - (xx-2)*(dim_i-2))

			threadnum = OMP_GET_THREAD_NUM()

			! Only do something if within watershed, if we care about the watershed
			if (wshed) then
				if (wsmask(xx,yy)==0) CYCLE
			end if

			! Only do something if rain fell at this point on this day
			if (SUM(precip(xx,yy,:))>minpre) then

				!$OMP CRITICAL (output_index)
				orec = orec + 1
				torec = orec
    
				! *Output results per parcel can be specified here.*
                if (eachParcel) then
         	        OPEN(unit=threadnum+10,file=TRIM(diro)//"parcel"//TRIM(int_to_string(dd))//"_"//TRIM(int_to_string(orec)), &
         	  	        form="UNFORMATTED",status="REPLACE") 
         	        !print *,threadnum+10
         	    end if 
                
				!$OMP END CRITICAL (output_index)
    
				WV_cont_day = 0.

				!
				! Determine how many parcels to release today and use precip
				! distribution to determine when to release parcels.
				! The globally set nparcels is just a maximum number of parcels
				! to release in each data timestep. We release at least one parcel
				! per simulation timestep within a data time step when it rained.
				! npar calculates how many parcels to release that day. parcel_release_time
				! spreads that number of parcels out of the 144 timesteps depending on
				! when it rained.
				!
				par_release = 0

				!$OMP CRITICAL (par_rel_time)
				if (COUNT(MASK = precip(xx,yy,:)>0.)<(nparcels/indatatsteps)) then
					npar = COUNT(MASK = precip(xx,yy,:)>0.) * indatatsteps
					call parcel_release_time(precip(xx,yy,:),npar,par_release)
				else
					npar = nparcels
					call parcel_release_time(precip(xx,yy,:),npar,par_release)
				end if
				!npar = 10*COUNT(MASK = precip(xx,yy,:)>0.)

				!print*,'npar',npar

				!$OMP END CRITICAL (par_rel_time)
	

				! * Parcel release height can be set here if you want to remove randomness.*

				!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				! FOR EVERY SIMULATION SUB-DAILY TIMESTEP
				!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				do tt = 1, daytsteps
					if (par_release(tt)==0) then
						CYCLE
					end if

					!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					! FOR EVERY LOT OF PARCELS
					!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					do mm = 1, par_release(tt)

						WV_cont = 0.
						fac = 1.
						x = xx
						y = yy

						!the input data time step before this parcel time step on the rain day
						ttdataday = INT((tt-1)/indatatsteps) + 1

						!the input data time step from the beginning of the loaded files
						ttdata = datatotsteps - datadaysteps - 1 + ttdataday

						!factor for linear interpolation to parcel time step
						ttfac = MOD(tt,indatatsteps)*1./indatatsteps

						!the precip produced here at this parcel time step
						end_precip = precip(xx,yy,ttdataday)/indatatsteps

						!determine model level from which to release parcel
						call parcel_release_height(pw(xx,yy,:,tt),par_lev)

						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						!release the parcel and track the back trajectory
						!until all of initial precipitable water is accounted for
						!or the parcel leaves the domain
						!or the user specified time to calculate the back trajectory runs out
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						!
						! We always release parcels from the centre of the grid cell. I think this differs to D&B: See D&B 2007 Figure 1.
						par_lat = lat2d(xx,yy)
						par_lon = lon2d(xx,yy)

						! Calculate the parcel mixing ratio. This is used in the calculation of new parcel level in new_parcel_level_w.
						par_q = lin_interp(mixtot(xx,yy,par_lev,ttdata:ttdata+1),ttfac)

						! * Parcel potential temperature was calculated here.*

						! Calculate parcel pressure.This is used in the calculation of new parcel level in new_parcel_level_w.
						par_pres = lin_interp(pres(xx,yy,par_lev,ttdata:ttdata+1),ttfac)

						!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						! FOR EACH PARCEL RELEASE TIME, FOR EACH SIMULATION TIME STEP IN THE
						! WHOLE BACK-TRACK PERIOD
						!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						do nn = totsteps-daytsteps+tt, 2, -1
							!
							!advect the parcel back in time one step
							!
                            !current parcel stats
                            if (eachParcel) then
                                !print *,"nn ",nn,threadnum,par_lev
                                parcel_stats(1,totsteps-daytsteps+tt+1-nn) = nn*1.
                                parcel_stats(2,totsteps-daytsteps+tt+1-nn) = xx
                                parcel_stats(3,totsteps-daytsteps+tt+1-nn) = yy
                                parcel_stats(4,totsteps-daytsteps+tt+1-nn) = par_lon
                                parcel_stats(5,totsteps-daytsteps+tt+1-nn) = par_lat
                                parcel_stats(6,totsteps-daytsteps+tt+1-nn) = par_pres
                                parcel_stats(7,totsteps-daytsteps+tt+1-nn) = par_lev
                                parcel_stats(8,totsteps-daytsteps+tt+1-nn) = par_q
                                parcel_stats(9,totsteps-daytsteps+tt+1-nn) = u(x,y,par_lev,ttdata)
                                parcel_stats(10,totsteps-daytsteps+tt+1-nn) = v(x,y,par_lev,ttdata)
                                parcel_stats(11,totsteps-daytsteps+tt+1-nn) = w(x,y,par_lev,ttdata)
                            end if
							!
							!calculate the lower left location for the subsection
							!
							if (x+floor(ssdim/2.)>dim_j) then
								ssx = dim_j - ssdim + 1
							else
								ssx = max(x-floor(ssdim/2.),1)
							end if

							if (y+floor(ssdim/2.)>dim_i) then
								ssy = dim_i - ssdim + 1
							else
								ssy = max(y-floor(ssdim/2.),1)
							end if

							! Find where you are in the simlength 3-hourly timeseries
							nnMM5 = INT(nn/indatatsteps) + 1
							nnfac = MOD(nn,indatatsteps)*1./indatatsteps

							!
							!get u,v and temp for this and the previous parcel time step (just subsection)
							!

							unow(:,:,:,2) = lin_interp3D(u(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,:,nnMM5:nnMM5+1),nnfac)

							vnow(:,:,:,2) = lin_interp3D(v(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,:,nnMM5:nnMM5+1),nnfac)

							wnow(:,:,:,2) = lin_interp3D(w(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,:,nnMM5:nnMM5+1),nnfac)

							! The temperature (now) is used to determine the temperature of the parcel before it's advected. (The initial pressure of the parcel was already calculated before nn. Subsequent parcel pressures, as the parcel is moved backward in each time step, are determined within the back-trajectory routine, or more specifically, during the routine to determine the parcel's new height (i.e. pressure).)
							tempnow(:,:,:) = lin_interp3D(temp(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,:,nnMM5:nnMM5+1),nnfac)


							! Find where you are in the nn timeseries
							nnMM5 = INT((nn-1)/indatatsteps) + 1
							nnfac = MOD(nn-1,indatatsteps)*1./indatatsteps

							unow(:,:,:,1) = lin_interp3D(u(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,:,nnMM5:nnMM5+1),nnfac)

							vnow(:,:,:,1) = lin_interp3D(v(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,:,nnMM5:nnMM5+1),nnfac)

							wnow(:,:,:,1) = lin_interp3D(w(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,:,nnMM5:nnMM5+1),nnfac)

							pres_then(:,:,:) = lin_interp3D(pres(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,:,nnMM5:nnMM5+1),nnfac)

							psfc_then(:,:) = lin_interp2D(psfc(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1,nnMM5:nnMM5+1),nnfac)

							! SPECIFY WHICH VERSION OF THE BACK-TRAJECTORY YOU WANT TO USE
							! Here parcels move with vertical wind speed (w) and have their new pressures calculated using actual temp

							call implicit_back_traj_w(unow,vnow,wnow,tempnow,pres_then,psfc_then,lon2d(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1),lat2d(ssx:ssx+ssdim-1,ssy:ssy+ssdim-1),par_lon,par_lat,par_lev,par_pres,par_q,threadnum)

							! Find the grid cell nearest the new lat,lon of the parcel

							! While in the first time step the parcel x,y may be the same cell as the parcel was released from, as you back-track that parcel in time the x,y will change.
							call near_pt(lon2d,lat2d,par_lon,par_lat,x,y)

							! Find the water mass contribution of the new grid square, at this time	      !
							new_par_q = lin_interp(mixtot(x,y,par_lev,nnMM5:nnMM5+1),nnfac)

							!
							!adjust the q reduction factor if we had a decrease in q
							!so long as it isn't the first time step.
							! i.e. If the amount of water in the atmosphere at the parcel position decreases backward in time, then the parcel q at the current time step must not have come from the cell evap...maybe from some other process like convection.
							!

                            !! TO DO - create user setting whether you want to split the PBL or not
                            
                            !was moisture contributed to the parcel?
							!is the parcel in the pbl?
                            ! Unlike WRF, ERA5 evap and twp units are consistent, so no need to divde by indatatsteps.   

							if (par_q > new_par_q+min_del_q) then  
								WV_cont(x,y) = WV_cont(x,y) + ((par_q - new_par_q)/par_q)*fac
								if (nn < totsteps-daytsteps+tt) then										
									fac = fac*(1-(par_q - new_par_q)/par_q)
								end if
							end if

							par_q = new_par_q

                            !saving parcel stats
                            if (eachParcel) then
                                parcel_stats(12,totsteps-daytsteps+tt+1-nn) = evap(x,y,ttdata)
                                parcel_stats(13,totsteps-daytsteps+tt+1-nn) = tpw(x,y,ttdata)
                                parcel_stats(14,totsteps-daytsteps+tt+1-nn) = WV_cont(x,y)
                            end if
              
							!
							!if we have accounted for all the precip  then go to next parcel
							!
							if (SUM(WV_cont)>=1.) then
    							WV_cont(x,y) = WV_cont(x,y) - (SUM(WV_cont) - 1)
								EXIT
							end if

							!
							!if qfac = 0 then all of the increases in the water vapor
							!further back along the trajectory are lost before reaching
							!the end precipitation point so they don't contribute and
							!there is no need to continue
							!
							!the water not accounted for must have come from convection
							!or some other process that remains unaccounted for
							!
							! if (qfac==0) then
							! 	EXIT
							! end if

							! if (qfac_evap<=0.00001) then
							! 	EXIT
							! end if

							!if we have left the domain then assign the remaining precip to
							!outside and go to next parcel

							if (x<2) then
								WV_cont(1,y) = 1. - SUM(WV_cont)
								EXIT
							else if (x>dim_j-2) then
								WV_cont(dim_j,y) = 1. - SUM(WV_cont)
								EXIT

							end if

							if (y<2) then
								WV_cont(x,1) = 1. - SUM(WV_cont)
								EXIT
							else if (y>dim_i-2) then
								WV_cont(x,dim_i) = 1. - SUM(WV_cont)
								EXIT
							end if

							!
							!if we have reached here and nn=2 then we have neither left the
							!domain nor acounted for the precip in the allocated back-track period (didn't go back far enough in time).
							!
						end do   !nn loop

						! wv_cont(x,y) is a 2d grid of E/TPW values. The grid is added to for every nn parcel back-track. E.g. in one 10min daytstep, we might release 1 parcel. This parcel will calculate the contribution from every cell in the grid. However we could release more, like 5 parcels. The contribution from the grid should be the same no matter how many parcels we release. So we take the average grid contribution per parcel released.
						!$OMP CRITICAL (wvcont_day1)
						WV_cont_day = WV_cont_day + WV_cont/npar
						!$OMP END CRITICAL (wvcont_day1)

						if (par_lev==0) then
							write(*,*) "par_lev==0"
							STOP
						end if

                        !if keeping track of each parcel
                        if (eachParcel) then
                            !print *,"output",threadnum+10
                            WRITE(threadnum+10) parcel_stats
                            CLOSE(threadnum+10)
                        end if

                        !print *, 'parcel_stats(:,:10)', parcel_stats(:,:2)

					end do  !mm loop

				end do  !tt loop

				!
				!write output to netcdf file
				!
				print*,MAXVAL(WV_cont_day)
				!$OMP CRITICAL (output)
				status = nf90_put_var(outncid,wvcid,WV_cont_day,start=(/1,1,torec/),count=(/dim_j,dim_i,1/))
				if(status /= nf90_NoErr) call handle_err(status)
				status = nf90_put_var(outncid,xlocid,xx,start=(/torec/))
				if(status /= nf90_NoErr) call handle_err(status)
				status = nf90_put_var(outncid,ylocid,yy,start=(/torec/))
				if(status /= nf90_NoErr) call handle_err(status)
				status = nf90_put_var(outncid,dayid,dd,start=(/torec/))
				if(status /= nf90_NoErr) call handle_err(status)
				status = nf90_put_var(outncid,opreid,SUM(precip(xx,yy,:)),start=(/torec/))
				if(status /= nf90_NoErr) call handle_err(status)
				!$OMP END CRITICAL (output)

			end if

		end do   !xx_omp loop
		!$OMP END DO NOWAIT
		!$OMP END PARALLEL

		status = nf90_close(outncid)
		if(status /= nf90_NoErr) call handle_err(status)

	end do !! dd loop

END PROGRAM back_traj