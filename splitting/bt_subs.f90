MODULE bt_subs

	IMPLICIT NONE

	CONTAINS

!***********************************************************************

	SUBROUTINE new_out_file(outncid,wvcid,wvc2id,xlocid,ylocid,dayid,opreid,daynum,lat2d,lon2d)
	!----------------------------------------------
	!create the output netcdf file and prepare it to accept data
	!--------------------------------------------------------

		USE global_data
		USE util
		USE netcdf

		IMPLICIT NONE

		INTEGER, INTENT(INOUT) :: outncid,wvcid,wvc2id,xlocid,ylocid,dayid,opreid
		INTEGER, INTENT(IN) :: daynum
		REAL,INTENT(IN),DIMENSION(:,:) :: lat2d,lon2d

		INTEGER :: status,jdimid,idimid,gwvcdimid,latid,lonid

		!
		!create the file
		!
		!differentiate whether we are doing whole days or around storm peaks
		if (peak) then
			print *,'we are doing peaks here!'
			if (mon<10) then
				status = nf90_create(TRIM(diro)//"bt."//TRIM(int_to_string(year))//"0" &
					//TRIM(int_to_string(mon))//"_"//TRIM(real_to_string(daynum))// &
					".nc",nf90_clobber,outncid)
				if (status /= NF90_NOERR) call handle_err(status)
			else
				status = nf90_create(TRIM(diro)//"bt."//TRIM(int_to_string(year)) &
					//TRIM(int_to_string(mon))//"_"//TRIM(real_to_string(daynum))// &
					".nc",nf90_clobber,outncid)
					if (status /= NF90_NOERR) call handle_err(status)
			end if
		else
			if (mon<10) then
				status = nf90_create(TRIM(diro)//"bt."//TRIM(int_to_string(year))//"0" &
				//TRIM(int_to_string(mon))//"_"//TRIM(int_to_string(INT(daynum)))// &
					".nc",nf90_clobber,outncid)
				if (status /= NF90_NOERR) call handle_err(status)
			else
				status = nf90_create(TRIM(diro)//"bt."//TRIM(int_to_string(year)) &
		  			//TRIM(int_to_string(mon))//"_"//TRIM(int_to_string(INT(daynum)))// &
					".nc",nf90_clobber,outncid)
				if (status /= NF90_NOERR) call handle_err(status)
			end if
		end if

		print *,'outfile=',TRIM(diro)//"bt."//TRIM(int_to_string(year))//"0" &
				//TRIM(int_to_string(mon))//"_"//TRIM(int_to_string(INT(daynum)))// &
				".nc"

		!
		!define dimensions
		!
		status = nf90_def_dim(outncid,"j_cross",dim_j,jdimid)
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_def_dim(outncid,"i_cross",dim_i,idimid)
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_def_dim(outncid,"gridcell_wvc",nf90_unlimited,gwvcdimid)
		if (status /= NF90_NOERR) call handle_err(status)

		!
		!define the variable
		!
		status = nf90_def_var(outncid,"wv_cont",nf90_float,(/jdimid,idimid,gwvcdimid/),wvcid)
		if (status /= NF90_NOERR) call handle_err(status)
		!turning off apbl output if necessary
        status = nf90_def_var(outncid,"wv_cont_apbl",nf90_float,(/jdimid,idimid,gwvcdimid/),wvc2id)
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_def_var(outncid,"x_loc",nf90_int,(/gwvcdimid/),xlocid)
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_def_var(outncid,"y_loc",nf90_int,(/gwvcdimid/),ylocid)
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_def_var(outncid,"day",nf90_float,(/gwvcdimid/),dayid)
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_def_var(outncid,"pre",nf90_float,(/gwvcdimid/),opreid)
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_def_var(outncid,"latitcrs",nf90_float,(/jdimid,idimid/),latid)
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_def_var(outncid,"longicrs",nf90_float,(/jdimid,idimid/),lonid)
		if (status /= NF90_NOERR) call handle_err(status)

		!
		!define attributes
		!
		status = nf90_put_att(outncid,wvcid,"long_name","Water Vapor Contribution")
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_put_att(outncid,wvcid,"units","proportion of precipitation")
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_put_att(outncid,wvcid,"num_boundary_layers",bdy)
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_put_att(outncid,wvcid,"parcels_per_grid_point",nparcels)
		if (status /= NF90_NOERR) call handle_err(status)

		!turn off apbl output if necessary
        status = nf90_put_att(outncid,wvc2id,"long_name","Water Vapor Contribution above PBL")
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_put_att(outncid,wvc2id,"units","proportion of precipitation")
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_put_att(outncid,wvc2id,"num_boundary_layers",bdy)
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_put_att(outncid,wvc2id,"parcels_per_grid_point",nparcels)
		if (status /= NF90_NOERR) call handle_err(status)

		status = nf90_put_att(outncid,xlocid,"long_name","x index location of precipitation (from 0)")
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_put_att(outncid,ylocid,"long_name","y index location of precipitation (from 0)")
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_put_att(outncid,dayid,"long_name","days since "// &
			TRIM(int_to_string(sday))//"/"//TRIM(int_to_string(smon))//"/"// &
			TRIM(int_to_string(syear)))
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_put_att(outncid,opreid,"long_name","precipitation")
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_put_att(outncid,opreid,"units","mm")
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_put_att(outncid,latid,"long_name","LATITUDE (SOUTH NEGATIVE)")
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_put_att(outncid,latid,"units","degrees")
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_put_att(outncid,lonid,"long_name","LONGITUDE (WEST NEGATIVE)")
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_put_att(outncid,lonid,"units","degrees")
		if (status /= NF90_NOERR) call handle_err(status)


		!
		!leave define mode
		!
		status = nf90_enddef(outncid)


		status = nf90_put_var(outncid,latid,lat2d,start=(/1,1/),count=(/dim_j,dim_i/))
		if(status /= nf90_NoErr) call handle_err(status)
		status = nf90_put_var(outncid,lonid,lon2d,start=(/1,1/),count=(/dim_j,dim_i/))
		if(status /= nf90_NoErr) call handle_err(status)

	END SUBROUTINE new_out_file

	!***********************************************************************


	!***********************************************************************  

	REAL FUNCTION lin_interp(var,fac)
	!---------------------------------------
	!linearly interpolate between the values of var
	!fac is the proporational distance from the first value
	!------------------------------------------

		IMPLICIT NONE

		REAL, INTENT(IN), DIMENSION(2) :: var
		REAL, INTENT(IN) :: fac

		lin_interp = var(1)*(1-fac) + var(2)*fac

	END FUNCTION lin_interp

	!***********************************************************************

	FUNCTION lin_interp2D(var,fac)
	!---------------------------------------
	!linearly interpolate between the values of var (last dimension must have size 2)
	!fac is the proporational distance from the first value
	!------------------------------------------

		IMPLICIT NONE

		REAL, INTENT(IN), DIMENSION(:,:,:) :: var
		REAL, INTENT(IN) :: fac

		REAL, DIMENSION(SIZE(var,1),SIZE(var,2)):: lin_interp2D

		lin_interp2D = var(:,:,1)*(1-fac) + var(:,:,2)*fac

	END FUNCTION lin_interp2D

	!***********************************************************************

	FUNCTION lin_interp3D(var,fac)
	!---------------------------------------
	!linearly interpolate between the values of var (last dimension must have size 2)
	!fac is the proporational distance from the first value
	!------------------------------------------

		IMPLICIT NONE

		REAL, INTENT(IN), DIMENSION(:,:,:,:) :: var
		REAL, INTENT(IN) :: fac

		REAL, DIMENSION(SIZE(var,1),SIZE(var,2),SIZE(var,3)) :: lin_interp3D

		lin_interp3D = var(:,:,:,1)*(1-fac) + var(:,:,:,2)*fac

	END FUNCTION lin_interp3D

	!***********************************************************************

	SUBROUTINE parcel_release_time(precip,npar,par_release)
	!---------------------------------------------------
	! here we calculate the times of day to release our parcels
	! based on a random precipitation weighted sampling
	!-----------------------------------------------------

		USE global_data

		IMPLICIT NONE

		REAL, DIMENSION(:), INTENT(IN) :: precip
		INTEGER, DIMENSION(:), INTENT(OUT) :: par_release
		INTEGER, INTENT(IN) :: npar

		REAL, DIMENSION(SIZE(precip)*indatatsteps) :: cumm_precip
		REAL, DIMENSION(npar) :: rand_nums
		INTEGER :: tt,rr,ss,rec


		par_release = 0

		call RANDOM_NUMBER(rand_nums)
		cumm_precip = 0.

		rec = 0

		do tt = 1,SIZE(precip)
			do ss = 1,indatatsteps
				rec = rec + 1
				if (rec==1) then
					cumm_precip(1) = precip(1)/indatatsteps
				else
					cumm_precip(rec) = cumm_precip(rec-1) + precip(tt)/indatatsteps
				end if
			end do
		end do

		cumm_precip = cumm_precip/cumm_precip(SIZE(cumm_precip))

		do rr = 1,npar
			do tt = 1,SIZE(cumm_precip)
				if (cumm_precip(tt)>rand_nums(rr)) then
					par_release(tt) = par_release(tt) + 1
					EXIT
				end if
			end do
		end do


	END SUBROUTINE parcel_release_time

	!***********************************************************************

	SUBROUTINE parcel_release_height(pw,par_lev)
	!----------------------------------------------
	! calculate the height to release the parcel from
	! based on precipitable water weighted random sampling
	!-----------------------------------------------

		IMPLICIT NONE

		REAL, INTENT(IN), DIMENSION(:) :: pw ! This is the precipitable water, accumulated from the ground up, in the column at that point in time.

		INTEGER, INTENT(OUT) :: par_lev

		REAL :: rand_num
		INTEGER :: kk

		call RANDOM_NUMBER(rand_num)

		! Take random number as a random proportion of the total pw in the column at that time; pw(1) is the pw at the top of the atm column, which represents the accumulated pw over the column below it (same as TPW).
		rand_num = rand_num*pw(1)

		do kk = SIZE(pw),1,-1 
			if (pw(kk)>rand_num) then
				par_lev = kk
				EXIT
			end if
		end do


		! For testing purposes only: take random number as a purely random model level, not weighted by pw.
		!rand_num = 1 + FLOOR(size(pw)*rand_num)
		!par_lev = rand_num

	END SUBROUTINE parcel_release_height

	!***********************************************************************

	SUBROUTINE lin_interp_inMM5tsteps(var)
	!------------------------------------------
	!linearly interpolate through time inside ERA5 time steps
	!----------------------------------------------

		USE global_data

		IMPLICIT NONE

		REAL, INTENT(INOUT), DIMENSION(:,:,:,:) :: var

		INTEGER :: i

		do i = 1,indatatsteps-1
		  var(:,:,:,i+1::indatatsteps) = (1-(i*1./indatatsteps))*var(:,:,:,1:datadaysteps+1-indatatsteps:indatatsteps) &
		                 & + (i*1./indatatsteps)*var(:,:,:,1+indatatsteps::indatatsteps)
		end do


	END SUBROUTINE lin_interp_inMM5tsteps

	!***********************************************************************

	SUBROUTINE calc_pw(mix,pres,surf_pres,ptop,pw)
	!------------------------------------------
	! calculate the precipitable water from input data fields
	! only for the day of current interest
	! save as accumulated field from the ground up (4D field)
	! This is used to calc parcel initial height.
	!-------------------------------------------------

		USE global_data

		IMPLICIT NONE

		REAL, INTENT(IN), DIMENSION(:,:,:,:) :: mix,pres
		REAL, INTENT(IN), DIMENSION(:,:,:) :: surf_pres
		REAL, INTENT(OUT), DIMENSION(:,:,:,:) :: pw

		REAL, DIMENSION(dim_j,dim_i,dim_k,SIZE(mix,4)) :: dp
		INTEGER :: k

		REAL, INTENT(IN) :: ptop

		!
		! calculate the change in pressure (Pa) represented by each point
		!
		!for highest level
		dp(:,:,1,:) = SUM(pres(:,:,:2,:),3)/2. - ptop
	
        !need to account for posibility that pressure levels go below the ground
        !for the middle levels
		do k = 2,dim_k-1
            where (pres(:,:,k+1,:) <= surf_pres(:,:,:)) 
                dp(:,:,k,:) = (pres(:,:,k+1,:) - pres(:,:,k-1,:)) /2. 
            elsewhere (pres(:,:,k,:) <= surf_pres(:,:,:))
                !for the lowest level above surface
                dp(:,:,k,:) = surf_pres(:,:,:) - (pres(:,:,k,:) + pres(:,:,k-1,:))/2.
            elsewhere
                dp(:,:,k,:) = 0.
            end where

		end do

		!for the lowest level
        where (pres(:,:,dim_k,:) <= surf_pres(:,:,:))
        dp(:,:,dim_k,:) = surf_pres(:,:,:) - SUM(pres(:,:,dim_k-1:,:),3)/2.
        elsewhere
                dp(:,:,dim_k,:) = 0.
        end where

		!mass in mm
		pw(:,:,:,::indatatsteps) = dp*mix/g

		!interpolate inside input data time steps
		call lin_interp_inMM5tsteps(pw)

		!accumulate from the bottom up. The precipitable water is then the total moisture in the column below it.
		do k = dim_k-1,1,-1             ! i.e. from level 28 to 1
			!pw(:,:,k,2:) = pw(:,:,k+1,2:) + pw(:,:,k,2:)! THE PW IS ACCUMULATED FROM THE SECOND TS ON, SO THE FIRST 10MIN IS NOT ACCUMULATED. WHY?? This only matters if tt=1. Here I change it to do all timesteps.
			pw(:,:,k,:) = pw(:,:,k+1,:) + pw(:,:,k,:)
		end do
		! print *,shape(pw)
		! print *,"pw ",pw(1,1,dim_k,:),surf_pres(1,1,1),ptop

	END SUBROUTINE calc_pw

	!***********************************************************************

	SUBROUTINE calc_tpw(mix,pres,surf_pres,ptop,tpw)
	!------------------------------------------
	! calculate the total precipitable water from input data fields
	! at every level and time (3D field)
	!-------------------------------------------------

		USE global_data

		IMPLICIT NONE

		REAL, INTENT(IN), DIMENSION(:,:,:,:) :: mix,pres
		REAL, INTENT(IN), DIMENSION(:,:,:) :: surf_pres
		REAL, INTENT(OUT), DIMENSION(:,:,:) :: tpw

		REAL, DIMENSION(dim_j,dim_i,dim_k,datatotsteps) :: dp
		INTEGER :: j,i,k,t

		REAL, INTENT(IN) :: ptop
		!
		! calculate the change in pressure (Pa) represented by each point
		!
		!for highest level
		dp(:,:,1,:) = SUM(pres(:,:,:2,:),3)/2. - ptop

        !need to account for posibility that pressure levels go below the ground
        !for the middle levels
        do k = 2,dim_k-1
			where (pres(:,:,k+1,:) <= surf_pres(:,:,:))
				dp(:,:,k,:) = (pres(:,:,k+1,:) - pres(:,:,k-1,:)) /2. !dp(:,:,k,:) = SUM(pres(:,:,k-1:k+1:2,:),3)/2.
			elsewhere (pres(:,:,k,:) <= surf_pres(:,:,:))
				!for the lowest level above surface
				dp(:,:,k,:) = surf_pres(:,:,:) - (pres(:,:,k,:) + pres(:,:,k-1,:))/2.
			elsewhere
				dp(:,:,k,:) = 0.
			end where

        end do

        !for the lowest level also accounting for possibility of being below ground
        where (pres(:,:,dim_k,:) <= surf_pres(:,:,:))
            dp(:,:,dim_k,:) = surf_pres(:,:,:) - SUM(pres(:,:,dim_k-1:,:),3)/2.
        elsewhere
            dp(:,:,dim_k,:) = 0.
        end where

		!mass in mm tpw
		do j = 1,dim_j
			do i = 1,dim_i
				do t = 1,datatotsteps
					tpw(j,i,t) = SUM(dp(j,i,:,t)*mix(j,i,:,t)/g)
				end do
			end do
		end do

	END SUBROUTINE calc_tpw

	!***********************************************************************

	SUBROUTINE calc_tpw_pbl(mix,pres,surf_pres,tpw,pbl_lev)
	!------------------------------------------
	! calculate the total precipitable water in the pbl from MM5 fields
	! at every level and time

     ! pw = 1/g.rho * integral(q.dp) over p
     ! units: pw [in meters] s-2/m . m-3/kg * kg/kg*kg/(m.s-2) recalling 1 Pa = 1 kg/(m.s-2)
	!-------------------------------------------------

		USE global_data

		IMPLICIT NONE

		REAL, INTENT(IN), DIMENSION(:,:,:,:) :: mix,pres
		REAL, INTENT(IN), DIMENSION(:,:,:) :: surf_pres
		INTEGER, INTENT(IN), DIMENSION(:,:,:) :: pbl_lev
		REAL, INTENT(OUT), DIMENSION(:,:,:) :: tpw

		REAL, DIMENSION(dim_j,dim_i,dim_k,datatotsteps) :: dp
		INTEGER :: j,i,k,t


		!
		! calculate the change in pressure (Pa) represented by each point
		!
		!for highest level
		dp(:,:,1,:) = SUM(pres(:,:,:2,:),3)/2.

        !need to account for posibility that pressure levels go below the ground
        !for the middle levels
        do k = 2,dim_k-1
                where (pres(:,:,k+1,:) <= surf_pres(:,:,:))
                    dp(:,:,k,:) = (pres(:,:,k+1,:) - pres(:,:,k-1,:)) /2. 
                elsewhere (pres(:,:,k,:) <= surf_pres(:,:,:))
					!for the lowest level above surface
					dp(:,:,k,:) = surf_pres(:,:,:) - (pres(:,:,k,:) + pres(:,:,k-1,:))/2.
                elsewhere
                    dp(:,:,k,:) = 0.
                end where

        end do

        !for the lowest level
        where (pres(:,:,dim_k,:) <= surf_pres(:,:,:))
            dp(:,:,dim_k,:) = surf_pres(:,:,:) - SUM(pres(:,:,dim_k-1:,:),3)/2.
        elsewhere
            dp(:,:,dim_k,:) = 0.
        end where

		!mass in mm
		do j = 1,dim_j
			do i = 1,dim_i
				do t = 1,datatotsteps
					tpw(j,i,t) = SUM(dp(j,i,pbl_lev(j,i,t):,t)*mix(j,i,pbl_lev(j,i,t):,t)/g)
				end do
			end do
		end do

    print *, 'mix(1,1,:,2)',mix(1,1,:,2)
    print *, 'dp(1,1,:,2)',dp(1,1,:,2)


	END SUBROUTINE calc_tpw_pbl

	!***********************************************************************

	SUBROUTINE calc_pot_temp(temp,pres,pot_temp)
	!------------------------------------------------

	! calculate the potential temperature at every level and time
	!-----------------------------------------------------

		USE global_data

		IMPLICIT NONE

		REAL, INTENT(IN), DIMENSION(:,:,:,:) :: temp,pres

		REAL, INTENT(OUT), DIMENSION(:,:,:,:) :: pot_temp

		!
		! calculate the potential temperature
		!
		pot_temp = temp * ((P0/pres)**(Rd/Cp))

	END SUBROUTINE calc_pot_temp

	!***********************************************************************

	SUBROUTINE calc_actual_temp(temp,pres,act_temp)
	!------------------------------------------------
    ! SUBROUTINE UNUSED

	! calculate the actual temperature at every level and time,
	! given the input perturbation potential temperature.
	!-----------------------------------------------------

		USE global_data

		IMPLICIT NONE

		REAL, INTENT(IN), DIMENSION(:,:,:,:) :: temp,pres

		REAL, INTENT(OUT), DIMENSION(:,:,:,:) :: act_temp

		!
		! Calculate the actual temperature [K] from potential temperature.
		! Need to add 300K since wrfout gives *pertubation* potential temperature as temp T.
		!
		act_temp = temp * ((P0/pres)**(-Rd/Cp))

	END SUBROUTINE calc_actual_temp

	!***********************************************************************

	SUBROUTINE calc_eq_pot_temp(mix,mixtot,temp,pres,eq_pot_temp)
	!------------------------------------------------
	! SUBROUTINE UNUSED

	! calculate the equivalent potential temperature at every level and time
	!-----------------------------------------------------

		USE global_data

		IMPLICIT NONE

		REAL, INTENT(IN), DIMENSION(:,:,:,:) :: mix,mixtot,temp,pres

		REAL, INTENT(OUT), DIMENSION(:,:,:,:) :: eq_pot_temp


		REAL, DIMENSION(dim_j,dim_i,dim_k,datatotsteps) :: es,mix_s,pot_temp

		!
		! use Teten's formula to calculate the saturation vapor pressure (in stull)
		!
		es = 611 * exp((17.2694*(temp - 273.16))/(temp - 35.86))

        !
		! calculate the saturated mixing ratio everywhere
		!
		mix_s = (Rd/Rv)*es/(pres - es)

		!
		! calculate the equivalent potential temperature (in atmospheric convection, Emanuel)
		!
		eq_pot_temp = temp * ((P0/pres)**(Rd/(Cp+Cl*mixtot))) * &
				mix_s**(-1.*mix*Rv/(Cp+Cl*mixtot)) * &
				exp(Lv*mix/((Cp+Cl*mixtot)*temp))

		!
		!calculate the equivalent potential temperature (in atmospheric science..., Wallace & Hobbs)
		!
		!this isn't right - need to know mix_s at the lifting condensation level!!
		!
		!call calc_pot_temp(temp,pres,pot_temp)

		!eq_pot_temp = pot_temp*exp(Lv*mix_s/(Cp*temp))

	END SUBROUTINE calc_eq_pot_temp

	!***********************************************************************

	SUBROUTINE calc_pbl_lev(pbl_hgt,pres,surf_pres,pbl_lev)
	!------------------------------------------------
	! SUBROUTINE UNUSED

	! calculate the model level just above the pbl height
	!-----------------------------------------------------

	USE global_data

		IMPLICIT NONE

		REAL, INTENT(IN), DIMENSION(:,:,:,:) :: pres
		REAL, INTENT(IN), DIMENSION(:,:,:) ::pbl_hgt,surf_pres

		INTEGER, INTENT(OUT), DIMENSION(:,:,:) :: pbl_lev

		REAL, DIMENSION(dim_j,dim_i,datatotsteps) :: pbl_pres
		INTEGER :: j,i,t
		INTEGER,DIMENSION(1) :: dummy_min

		!
		! calculate pressure at the pbl height using the hydrostatic equation
		! here I assume that the density averages 1kg m-3 in the pbl
		!
		pbl_pres = -1*pbl_hgt*g + surf_pres        

		!
		! calculate the model level just above the pbl height
		! also add the level above that as it gains moisture by detrainment from the
		! PBL and so this moisture can also be associated with the current
		! location
		!
		do j = 1,dim_j
			do i = 1,dim_i
				do t = 1,datatotsteps
					dummy_min = MINLOC(abs(pbl_pres(j,i,t) - pres(j,i,:,t)))
					if ((pbl_pres(j,i,t) - pres(j,i,dummy_min(1),t)) < 0.) then
						pbl_lev(j,i,t) = dummy_min(1) - 2
					else
						pbl_lev(j,i,t) = dummy_min(1) - 1
					end if
				end do
			end do
		end do

	END SUBROUTINE calc_pbl_lev

	!***********************************************************************

	SUBROUTINE near_pt(lon2d,lat2d,lon,lat,x,y)
	!---------------------------------------------------------
	!calculate the grid point nearest the lat and lon location
	!------------------------------------------------------

		USE global_data

		IMPLICIT NONE

		REAL, INTENT(IN), DIMENSION(:,:) :: lon2d,lat2d
		REAL, INTENT(IN) :: lon,lat
		INTEGER, INTENT(OUT) :: x,y

		REAL, DIMENSION(SIZE(lon2d(:,1)),SIZE(lon2d(1,:))) :: dist
		INTEGER, DIMENSION(2) :: loc

		REAL, DIMENSION(SIZE(lon2d(:,1)),SIZE(lon2d(1,:))) :: lcos ! --svetlana

		!
		!calculate the distance from the parcel location to every grid point
		!must account for changing distance between longitude lines as latitude changes
		!

		lcos=cos(lat2d*pi/180)
		dist = sqrt((lat2d-lat)**2 + (lcos*(lon2d-lon))**2) ! --svetlana

		loc = MINLOC(dist)

		x = loc(1)
		y = loc(2)

	END SUBROUTINE near_pt

	!***********************************************************************

	SUBROUTINE bilin_interp(var2d,lon2d,lat2d,x,y,par_lon,par_lat,var)
	!------------------------------------------------------------------
	! find the bi-linearly interpolated value at par_lon,par_lat
	! lon and lat are not regularly spaced grids
	!------------------------------------------------------------------

		USE global_data

		IMPLICIT NONE

		REAL, INTENT(IN), DIMENSION(:,:) :: var2d,lon2d,lat2d
		REAL, INTENT(IN) :: par_lon,par_lat
		REAL, INTENT(OUT) :: var

		REAL :: fac,t,u
		INTEGER :: xx,yy,xbdy,ybdy,x,y
		LOGICAL :: changex,changey

		changex = .FALSE.
		changey = .FALSE.
		xbdy = 0
		ybdy = 0

		! Find what xx,yy is in the subgrid
		! call near_pt(lon2d,lat2d,par_lon,par_lat,xx,yy)
		! Instead of calling near_pt every time you do bilin_interp, just give bilin_interp the xx and yy
		xx=x
		yy=y

		!
		!check if we are currently exactly on a grid pt
		! ..where lon/lat2d are the subgrids
		! The initial call to bilin_interp in the implicit_back_traj_w will have the parcel on a grid point, as this is the initial location of the parcel when released (par_lat,par_lon) before entering the back_traj routine. Parcel lat,lon may not be on a cell centre after the parcel has been advected.
		If (lon2d(xx,yy)==par_lon.AND.lat2d(xx,yy)==par_lat) then
			var = var2d(xx,yy)
			RETURN
		end if

		!
		!get indices of closest grid value to south and west
		!be careful of boundary
		! If you are in the bottom corner of the subgrid:
		if (xx==1) xbdy = -1
		if (yy==1) ybdy = -1

		! If you're not in the bottom corner of the subgrid:
		! If the nearest cell point is further east than the advected parcel, take x position as xx-1, i.e. move it west.
		if (lon2d(xx,yy)>par_lon.AND.xbdy>-1) then
			changex = .TRUE.
		end if
		! If the nearest cell point is further north than the advected parcel, take y position as yy-1, i.e. move it south.
		if (lat2d(xx,yy)>par_lat.AND.ybdy>-1) then
			changey = .TRUE.
		end if

		! What about case where xx,yy need to increase, not decrease?

		if (changex) then
			xx = xx-1
		end if
		if (changey) then
			yy = yy-1
		end if

		!
		!if we are at the top or right boundary
		!
		if (xx==ssdim) then
			xbdy = 1
		end if
		if (yy==ssdim) then
			ybdy = 1
		end if

		!
		!check to see if point in inside lower and left boundaries
		!
		if (xbdy<0.AND.lon2d(xx,yy)<par_lon) xbdy = 0
		if (ybdy<0.AND.lat2d(xx,yy)<par_lat) ybdy = 0

		!
		!calculate the t and u weights
		!account for possibility of being outside the boundaries
		!even though this should never happen
		!if outside a boundary use the value at the boundary
		!
		if ((xbdy/=0.AND.ybdy/=0)) then
			var = var2d(xx,yy)
		else if (xbdy/=0) then
			fac = (par_lat - lat2d(xx,yy))/(lat2d(xx,yy+1)-lat2d(xx,yy))
			var = lin_interp(var2d(xx,yy:yy+1),fac)
		else if (ybdy/=0) then
			fac = (par_lon - lon2d(xx,yy))/(lon2d(xx+1,yy)-lon2d(xx,yy))
			var = lin_interp(var2d(xx:xx+1,yy),fac)
		else
			t = (par_lon - lon2d(xx,yy))/(lon2d(xx+1,yy)-lon2d(xx,yy))
			u = (par_lat - lat2d(xx,yy))/(lat2d(xx,yy+1)-lat2d(xx,yy))
			var = (1-t)*(1-u)*var2d(xx,yy) + (1-t)*u*var2d(xx,yy+1) + &
			       t*u*var2d(xx+1,yy+1) + t*(1-u)*var2d(xx+1,yy)
		end if

	END SUBROUTINE bilin_interp

	!***********************************************************************

	SUBROUTINE new_parcel_level_w(par_pres,pres,w,temp,mix,lev,psfc)
	!-------------------------------------------------------------------------
	!calculate the new parcel level given w at this location
	!--------------------------------------------------------------------------

		USE global_data

		IMPLICIT NONE

		REAL, INTENT(IN) :: w,temp,mix,psfc
		REAL, INTENT(IN), DIMENSION(:) :: pres
		REAL, INTENT(INOUT) :: par_pres
		INTEGER, INTENT(OUT) :: lev

		INTEGER, DIMENSION(1) :: dummy_lev

                !If w is in ms-1 then use this
                !
		! Here I use the hydrostatic eqn to calculate the change in pressure given w.
		! deltaP = rho*g*deltaz when in hydrostatic equilibrium
		! Note that "(1+0.61*mix)*temp" is the virtual temp. See p80 Wallace & Hobbs.
		!

		!par_pres = par_pres + -1.*(par_pres/(Rd*(1+0.61*mix)*temp))*g*w*tstep*60

		!If w is in Pas-1 (so it is really omega) then use this
		par_pres = par_pres + w*tstep*60

		!if the parcel is below the surface pressure then move it to 5hPa above the surface
		if (par_pres > psfc) par_pres = psfc - 5. 

		! Find the model level where the difference in pressure between the parcel
		! and the atmosphere is the smallest, i.e. which height in pres does the
		! smallest difference occur, where pres dims are (lat,lon,height).
		dummy_lev = MINLOC(ABS(pres - par_pres))

		lev = dummy_lev(1)

		!if the parcel is below the lowest model level then set it to the lowest level
		!if (par_pres > MAXVAL(pres)) par_pres = MAXVAL(pres)
                
		!make sure the level used is above the surface pressure (not underground)
		if (pres(lev) > psfc) lev = lev - 1

		! if (lev==0) then
		!   print *,'par_lev_w - pres_dis',(pres - par_pres),temp,w
		!   print *,'par_lev_w -',par_pres,pres
		! end if

	END SUBROUTINE new_parcel_level_w

	!***********************************************************************

	SUBROUTINE new_parcel_level_pt(par_pot_temp,pot_temp,par_lev,lev)
	!-------------------------------------------------------------------------
	! SUBROUTINE UNUSED

	!calculate what level has the parcel potential temperature at this location
	!--------------------------------------------------------------------------

		USE global_data

		IMPLICIT NONE

		REAL, INTENT(IN), DIMENSION(:) :: pot_temp
		INTEGER, INTENT(IN) :: par_lev
		REAL, INTENT(INOUT) :: par_pot_temp
		INTEGER, INTENT(OUT) :: lev

		REAL, DIMENSION(dim_k) :: pot_dist
		REAL :: pot_min
		INTEGER :: xx,yy,kk
		INTEGER, DIMENSION(1) :: lev_dummy
		INTEGER, DIMENSION(2) :: lev_dis
		LOGICAL :: getdownmin, getupmin


		getdownmin = .TRUE.
		getupmin = .TRUE.


		!
		!adjust parcel potential temperature if required (ie. driven into ground)
		!
		!print *,'1',par_pot_temp,xx,yy
		pot_min = MINVAL(pot_temp)
		par_pot_temp = MAX(par_pot_temp,pot_min)
		!print *,'2',pot_temp(xx,yy,:)
		!
		!calculate nearest potential temperature in vertical column
		!
		pot_dist = pot_temp - par_pot_temp

		!
		!since equivalent pot temperature is not monotonic I need
		!to find the level closest to previous level with right potential temp
		!
		lev_dis = dim_k*2

		do kk = 1, dim_k
			if (par_lev-kk>0 .AND. getdownmin) then
				if (pot_dist(par_lev-kk)<0.AND.pot_dist(par_lev-kk+1)>0.OR. &
						pot_dist(par_lev-kk)>0.AND.pot_dist(par_lev-kk+1)<0.) then
					lev_dummy = MINLOC(ABS(pot_dist(par_lev-kk:par_lev-kk+1)))
					lev_dis(1) = -kk-1+lev_dummy(1)
					getdownmin = .FALSE.
				end if
			end if

			if (par_lev+kk<dim_k+1 .AND. getupmin) then
				if (pot_dist(par_lev+kk)<0.AND.pot_dist(par_lev+kk-1)>0.OR. &
						pot_dist(par_lev+kk)>0.AND.pot_dist(par_lev+kk-1)<0.) then
					lev_dummy = MINLOC(ABS(pot_dist(par_lev+kk-1:par_lev+kk)))
					lev_dis(2) = kk-2+lev_dummy(1)
					getupmin = .FALSE.
				end if
			end if
		end do

		!print *,lev_dis,SUM(lev_dis),par_pot_temp

		if (SUM(lev_dis)==dim_k*4) then
			lev_dummy = MINLOC(ABS(pot_dist)) !this is when 'par_pot_temp' equals 'pot_temp' at some level
			lev = lev_dummy(1)
		else
			lev_dummy = MINLOC(ABS(lev_dis))
			lev = lev_dis(lev_dummy(1)) + par_lev
		end if

		if (lev==0) then
			print *,"new_parcel_pt - pot dist",pot_dist
			print *,"new_parcel_pt - par_lev,lev",par_lev,lev
		end if

	END SUBROUTINE new_parcel_level_pt

	!***********************************************************************

	SUBROUTINE advect(u,v,lon,lat)
	!-----------------------------------------------
	!linear advection in direction given by u and v
	!u and v are m/s
	!lon and lat are in degrees east and north
	!time step is given by tstep from global_data
	!-----------------------------------------------

		USE global_data

		IMPLICIT NONE

		REAL, INTENT(IN) :: u,v
		REAL, INTENT(INOUT) :: lon,lat

		!
		!calculate the new lat
		!
		lat = lat + v*tstep*60/(deg_dist*1000)

		!
		!calculate new lon
		!
		lon = lon + u*tstep*60/(cos(lat*pi/180)*deg_dist*1000)
		if (lon>360) then
			lon = lon - 360
		else if (lon<0) then
			lon = lon + 360
		end if


	END SUBROUTINE advect

	!***********************************************************************

	SUBROUTINE implicit_back_traj(u,v,w,temp,pbl_lev,pot_temp,pres,psfc,lon2d,lat2d, &
					par_lon,par_lat,par_lev, &
					par_pot_temp,par_pres,par_q,thread)
	!-------------------------------------------------------------------------------
	! SUBROUTINE UNUSED

	! Using Merrill's fully implicit isentropic technique
	! calculate the parcels position one time step before
	!
	!u,v,w should only have 2 time steps in them
	!pot_temp & pres should only have the end time step
	!--------------------------------------------------------------------------

		USE global_data

		IMPLICIT NONE

		REAL, INTENT(IN), DIMENSION(:,:,:,:) :: u,v,w
		REAL, INTENT(IN), DIMENSION(:,:,:) :: pot_temp,pres,temp
		REAL, INTENT(IN), DIMENSION(:,:) :: lon2d,lat2d,psfc
		INTEGER, INTENT(IN), DIMENSION(:,:) :: pbl_lev
		REAL, INTENT(INOUT) :: par_lon,par_lat,par_pot_temp,par_pres
		REAL, INTENT(IN) :: par_q
		INTEGER, INTENT(INOUT) :: par_lev
		INTEGER, INTENT(IN) :: thread

		INTEGER :: xx,yy,ll,lev
		REAL :: lon,lat,u_back,v_back,w_back,temp_back,u_for,v_for,w_for,pr
		REAL :: pt1,pt2,vfac,pr1,pr2
		INTEGER, DIMENSION(1) :: dummy_lev

		!print *,'1st',par_pres,par_pot_temp,par_lon,par_lat,par_lev,thread
		!
		!get u and v at parcel location
		!
		lon = par_lon
		lat = par_lat
		call near_pt(lon2d,lat2d,lon,lat,xx,yy)
		call bilin_interp(u(:,:,par_lev,2),lon2d,lat2d,xx,yy,lon,lat,u_back)
		call bilin_interp(v(:,:,par_lev,2),lon2d,lat2d,xx,yy,lon,lat,v_back)

		!
		!find lat and lon after advecting back in time
		!
		call advect(-1.*u_back,-1.*v_back,lon,lat)

		!
		!calculate which vertical level has correct potential temperature
		!at new location or that w moves us to
		!
		call near_pt(lon2d,lat2d,lon,lat,xx,yy)

		!print *,'par_lev,pbl_lev',par_lev,pbl_lev(xx,yy),xx,yy,thread

		call bilin_interp(temp(:,:,par_lev),lon2d,lat2d,xx,yy,par_lon,par_lat,temp_back)

		!if (par_lev >= pbl_lev(xx,yy)) then
		if (.TRUE.) then
			call new_parcel_level_pt(par_pot_temp,pot_temp(xx,yy,:),par_lev,lev)
		else
			pr = par_pres
			call bilin_interp(w(:,:,par_lev,2),lon2d,lat2d,xx,yy,lon,lat,w_back)
                        call new_parcel_level_w(pr,pres(xx,yy,:),w_back,temp_back,par_q,lev,psfc(xx,yy))
		end if

		!print *,'2nd',par_pres,par_pot_temp,par_lon,par_lat,par_lev,thread

		!
		!get u and v at new location
		!
		call bilin_interp(u(:,:,lev,1),lon2d,lat2d,xx,yy,lon,lat,u_for)
		call bilin_interp(v(:,:,lev,1),lon2d,lat2d,xx,yy,lon,lat,v_for)

		!
		!find new location of parcel as mean location given by back and forward trajectory
		!
		call advect(-1.*(u_back+u_for)/2.,-1.*(v_back+v_for)/2.,par_lon,par_lat)

		!
		!find final parcel potential temperature and level
		!
		call near_pt(lon2d,lat2d,par_lon,par_lat,xx,yy)

		!print *,'par_lev,pbl_lev',par_lev,pbl_lev(xx,yy),xx,yy,thread

		!if (par_lev >= pbl_lev(xx,yy)) then
		if (.TRUE.) then
			call new_parcel_level_pt(par_pot_temp,pot_temp(xx,yy,:),par_lev,lev)

			!print *,'lev_pt1',par_pres,par_pot_temp,pot_temp(xx,yy,lev-1:lev+1)

			!
			!need to calculate the new parcel pressure
			!need to be extra careful as pot_temp may not be monotonic
			!
			call bilin_interp(pres(:,:,lev),lon2d,lat2d,xx,yy,par_lon,par_lat,pr1)
			if (lev == dim_k .OR. lev == 1 .OR. par_pot_temp == pot_temp(xx,yy,lev)) then
				par_pres = pr1
			else
				!
				!in case this level is a min or max in pot_temp
				!
				if (par_pot_temp>pot_temp(xx,yy,lev).AND.par_pot_temp<pot_temp(xx,yy,lev-1).AND.par_pot_temp<pot_temp(xx,yy,lev+1)) then

					dummy_lev = MINLOC(ABS(par_pres - pres(xx,yy,lev-1:lev+1:2)))
					if (dummy_lev(1)==1) then
						dummy_lev(1) = lev-1
					else
						dummy_lev(1) = lev+1
					end if
					call bilin_interp(pres(:,:,dummy_lev(1)),lon2d,lat2d,xx,yy,par_lon,par_lat,pr2)
					vfac = (pot_temp(xx,yy,dummy_lev(1))-par_pot_temp)/(pot_temp(xx,yy,dummy_lev(1))-pot_temp(xx,yy,lev))
					par_pres = exp((1-vfac)*log(pr2) + vfac*log(pr1))

				else if (par_pot_temp<pot_temp(xx,yy,lev).AND.par_pot_temp>pot_temp(xx,yy,lev-1).AND.par_pot_temp>pot_temp(xx,yy,lev+1)) then

					dummy_lev = MINLOC(ABS(par_pres - pres(xx,yy,lev-1:lev+1:2)))
					if (dummy_lev(1)==1) then
						dummy_lev(1) = lev-1
					else
						dummy_lev(1) = lev+1
					end if
					call bilin_interp(pres(:,:,dummy_lev(1)),lon2d,lat2d,xx,yy,par_lon,par_lat,pr2)
					vfac = (pot_temp(xx,yy,lev)-par_pot_temp)/(pot_temp(xx,yy,lev)-pot_temp(xx,yy,dummy_lev(1)))
					par_pres = exp((1-vfac)*log(pr1) + vfac*log(pr2))

				!
				!in other cases
				!
				else if (par_pot_temp > pot_temp(xx,yy,lev) .AND. par_pot_temp < pot_temp(xx,yy,lev-1)) then

					call bilin_interp(pres(:,:,lev-1),lon2d,lat2d,xx,yy,par_lon,par_lat,pr2)
					vfac = (pot_temp(xx,yy,lev-1)-par_pot_temp)/(pot_temp(xx,yy,lev-1)-pot_temp(xx,yy,lev))
					par_pres = exp((1-vfac)*log(pr2) + vfac*log(pr1))

				else if (par_pot_temp < pot_temp(xx,yy,lev) .AND. par_pot_temp > pot_temp(xx,yy,lev+1)) then

					call bilin_interp(pres(:,:,lev+1),lon2d,lat2d,xx,yy,par_lon,par_lat,pr2)
					vfac = (pot_temp(xx,yy,lev)-par_pot_temp)/(pot_temp(xx,yy,lev)-pot_temp(xx,yy,lev+1))
					par_pres = exp((1-vfac)*log(pr1) + vfac*log(pr2))

				else if (par_pot_temp > pot_temp(xx,yy,lev) .AND. par_pot_temp < pot_temp(xx,yy,lev+1)) then

					call bilin_interp(pres(:,:,lev+1),lon2d,lat2d,xx,yy,par_lon,par_lat,pr2)
					vfac = (pot_temp(xx,yy,lev+1)-par_pot_temp)/(pot_temp(xx,yy,lev+1)-pot_temp(xx,yy,lev))
					par_pres = exp((1-vfac)*log(pr2) + vfac*log(pr1))

				else if (par_pot_temp < pot_temp(xx,yy,lev) .AND. par_pot_temp > pot_temp(xx,yy,lev-1)) then

					call bilin_interp(pres(:,:,lev-1),lon2d,lat2d,xx,yy,par_lon,par_lat,pr2)
					vfac = (pot_temp(xx,yy,lev)-par_pot_temp)/(pot_temp(xx,yy,lev)-pot_temp(xx,yy,lev-1))
					par_pres = exp((1-vfac)*log(pr1) + vfac*log(pr2))

				end if
			end if
  
		!print *,'lev_pt2',par_pres,par_pot_temp,pot_temp(xx,yy,lev-1:lev+1)

		else
			call bilin_interp(w(:,:,par_lev,1),lon2d,lat2d,xx,yy,lon,lat,w_for)
			call new_parcel_level_w(par_pres,pres(xx,yy,:),(w_back+w_for)/2.,temp_back,par_q,lev,psfc(xx,yy))

			!need to calculate the new parcel potential temperature
			call bilin_interp(pot_temp(:,:,lev),lon2d,lat2d,xx,yy,par_lon,par_lat,pt1)
			if (lev == dim_k .OR. lev == 1 .OR. par_pres == pres(xx,yy,lev)) then
				par_pot_temp = pt1
			else
				if (par_pres > pres(xx,yy,lev)) then
					call bilin_interp(pot_temp(:,:,lev+1),lon2d,lat2d,xx,yy,par_lon,par_lat,pt2)
					vfac = (log(pres(xx,yy,lev+1))-log(par_pres))/(log(pres(xx,yy,lev+1))-log(pres(xx,yy,lev)))
					par_pot_temp = (1-vfac)*pt2 + vfac*pt1
				else
					call bilin_interp(pot_temp(:,:,lev-1),lon2d,lat2d,xx,yy,par_lon,par_lat,pt2)
					vfac = (log(pres(xx,yy,lev))-log(par_pres))/(log(pres(xx,yy,lev))-log(pres(xx,yy,lev-1)))
					par_pot_temp = (1-vfac)*pt1 + vfac*pt2
				end if
			end if
		end if


		if (lev==0) then
		  print *,'L2389, ',par_lev,lev,par_pres,par_pot_temp,temp_back,thread
		  STOP
		end if

		par_lev = lev

	END SUBROUTINE implicit_back_traj

	!***********************************************************************

	SUBROUTINE implicit_back_traj_w(u,v,w,temp,pres,psfc,lon2d,lat2d, &
					par_lon,par_lat,par_lev, &
					par_pres,par_q,thread)
	!-------------------------------------------------------------------------------
	! Using Merrill's fully implicit  technique
	! calculate the parcels position one time step before
	!
	! u,v,w should only have 2 time steps in them
	! pres should only have the end time step

	! Output parcel lat,lon, height and pressure
	!--------------------------------------------------------------------------

		USE global_data

		IMPLICIT NONE

		REAL, INTENT(IN), DIMENSION(:,:,:,:) :: u,v,w
		REAL, INTENT(IN), DIMENSION(:,:,:) :: pres,temp
		REAL, INTENT(IN), DIMENSION(:,:) :: psfc,lon2d,lat2d
		REAL, INTENT(INOUT) :: par_lon,par_lat,par_pres
		REAL, INTENT(IN) :: par_q
		INTEGER, INTENT(INOUT) :: par_lev
		INTEGER, INTENT(IN) :: thread

		INTEGER :: xx,yy,lev !ll
		REAL :: lon,lat,u_back,v_back,w_par,temp_par,u_for,v_for,pr

		! Get u and v at parcel location, by interpolating in time (current and previous parcel timestep), and space (parcel lat/lon and lat/lon of the nearest grid point).
		lon = par_lon
		lat = par_lat

		call near_pt(lon2d,lat2d,lon,lat,xx,yy)

		call bilin_interp(u(:,:,par_lev,2),lon2d,lat2d,xx,yy,lon,lat,u_back) ! where lon2d/lat2d are the subgrids
		call bilin_interp(v(:,:,par_lev,2),lon2d,lat2d,xx,yy,lon,lat,v_back)


		!find lat and lon after advecting back in time. Reverse the wind directions (make negative) as you're going back in time.
		call advect(-1.*u_back,-1.*v_back,lon,lat)

		!calculate which vertical level that w moves us to
		call near_pt(lon2d,lat2d,lon,lat,xx,yy)

        !print *,'par_pres,psfc,pres',par_pres,psfc(xx,yy),pres(xx,yy,:)
		pr = par_pres
		call bilin_interp(temp(:,:,par_lev),lon2d,lat2d,xx,yy,lon,lat,temp_par)
		call bilin_interp(w(:,:,par_lev,2),lon2d,lat2d,xx,yy,lon,lat,w_par)
		! Reverse vertical wind direction as you're going back in time.
		call new_parcel_level_w(pr,pres(xx,yy,:),-1.*w_par,temp_par,par_q,lev,psfc(xx,yy))

		!get u and v at new lon/lat found using first advect call above
		call bilin_interp(u(:,:,lev,1),lon2d,lat2d,xx,yy,lon,lat,u_for)
		call bilin_interp(v(:,:,lev,1),lon2d,lat2d,xx,yy,lon,lat,v_for)

		!
		!find new location of parcel as mean location given by back and forward trajectory
		call advect(-1.*(u_back+u_for)/2.,-1.*(v_back+v_for)/2.,par_lon,par_lat)

		!
		!find final parcel level
		call near_pt(lon2d,lat2d,par_lon,par_lat,xx,yy)

		call bilin_interp(temp(:,:,par_lev),lon2d,lat2d,xx,yy,par_lon,par_lat,temp_par)
		call bilin_interp(w(:,:,par_lev,1),lon2d,lat2d,xx,yy,par_lon,par_lat,w_par)
		call new_parcel_level_w(par_pres,pres(xx,yy,:),-1.*w_par,temp_par,par_q,lev,psfc(xx,yy))

		if (lev==0) then
		  print *,'parcel lev=0',par_lev,lev,par_pres,temp_par,thread
		  STOP
		end if

		par_lev = lev
                !print *,'parvel_lev=',par_lev

	END SUBROUTINE implicit_back_traj_w

	!***********************************************************************

END MODULE bt_subs