MODULE input_data_handling_era5

	IMPLICIT NONE

	INTERFACE get_era5_field
		MODULE PROCEDURE :: get_era5_field_r1, get_era5_field_r2, get_era5_field_r3, get_era5_field_r4
	END INTERFACE

	CONTAINS


	!!! This series of functions assumes you've already worked out how
	!!! big these fields are going to be
	SUBROUTINE get_era5_field_r1(d,m,y,field,out,start,count)

		USE netcdf
		USE global_data, ONLY: dirdata_era5, fdim_j, dim_j_start, dim_j_end
		USE util, ONLY: handle_err

		IMPLICIT NONE

		INTEGER, INTENT(IN) :: d,m,y
		CHARACTER(len=*), INTENT(IN) :: field
		REAL, DIMENSION(:) :: out
		INTEGER, OPTIONAL, INTENT(IN) :: start, count

		!!! Locals
		CHARACTER(len=100) :: fn
		INTEGER :: fid, vid
		INTEGER :: status
		REAL    :: scale_factor, add_offset
		INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: dimids
		CHARACTER(len=NF90_MAX_NAME) :: dim_name
		INTEGER :: lon_idx = -1

		call get_filename(d,m,y,field,fn)

		status = NF90_OPEN(TRIM(fn),NF90_NOWRITE,fid)
		if (status /= NF90_NOERR) call handle_err(status)

		status = nf90_inq_varid(fid, trim(field), vid)
		if(status /= nf90_NoErr) call handle_err(status)

		!!! Are any of our variable's dimensions longitude?
		status = nf90_inquire_variable(fid,vid,dimids=dimids)
		if(status /= nf90_NoErr) call handle_err(status)

		status = nf90_inquire_dimension(fid,dimids(1),name=dim_name)
		if(status /= nf90_NoErr) call handle_err(status)

		if ( TRIM(dim_name) == "longitude" ) then
			lon_idx = 1
		endif

		if(PRESENT(start)) then
			if(PRESENT(count)) then
				!!! If these conditions are true we're going to cross a periodic boundary
				if( lon_idx > 0 .and. start + count > fdim_j ) then
					status = nf90_get_var(fid,vid,out(:fdim_j-dim_j_start+1),start=(/dim_j_start/),count=(/fdim_j - dim_j_start + 1/))
					if(status /= nf90_NoErr) call handle_err(status)
					status = nf90_get_var(fid,vid,out(fdim_j-dim_j_start+2:),start=(/1/),count=(/dim_j_end/))
				else
					status = nf90_get_var(fid,vid,out,start=(/start/),count=(/count/))
				end if
			else
				status = nf90_get_var(fid,vid,out,start=(/start/))
			end if
		else
			if(PRESENT(count)) then
				status = nf90_get_var(fid,vid,out,count=(/count/))
			else
				status = nf90_get_var(fid,vid,out)
			end if
		end if
		if(status /= nf90_NoErr) call handle_err(status)

		!!! Get scale_factor and add_offset (if any)
    status = nf90_get_att(fid, vid, "scale_factor", scale_factor)
    if( status /= nf90_NoErr ) then
      !!! OK if attribute not found - set it to 1.0
      if ( status /= nf90_eNotAtt ) call handle_err(status)
      scale_factor = 1.0
    endif

    status = nf90_get_att(fid, vid, "add_offset", add_offset)
    if( status /= nf90_NoErr ) then
      !!! OK if attribute not found - set it to 0.0
      if ( status /= nf90_eNotAtt ) call handle_err(status)
      add_offset = 0.0
    endif

		out = scale_factor * out + add_offset

	END SUBROUTINE

	SUBROUTINE get_era5_field_r2(d,m,y,field,out,starts,counts)

		USE netcdf
		USE global_data, ONLY: dirdata_era5, fdim_j, dim_j_start, dim_j_end
		USE util, ONLY: handle_err

		IMPLICIT NONE

		INTEGER, INTENT(IN) :: d,m,y
		CHARACTER(len=*), INTENT(IN) :: field
		REAL, DIMENSION(:,:) :: out
		INTEGER, OPTIONAL, DIMENSION(2), INTENT(IN) :: starts, counts

		!!! Locals
		CHARACTER(len=100) :: fn
		INTEGER :: fid, vid
		INTEGER :: status
		REAL    :: scale_factor, add_offset
		INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: dimids
		CHARACTER(len=NF90_MAX_NAME) :: dim_name
		INTEGER :: lon_idx = -1
		INTEGER :: idim

		INTEGER, DIMENSION(2) :: temp_starts
		INTEGER, DIMENSION(2) :: temp_counts
		INTEGER, DIMENSION(2,2) :: read_bounds !!! (/ (/ end1, start2 /), (/ end1, start2 /) /)

		call get_filename(d,m,y,field,fn)

		status = NF90_OPEN(TRIM(fn),NF90_NOWRITE,fid)
		if (status /= NF90_NOERR) call handle_err(status)

		status = nf90_inq_varid(fid, trim(field), vid)
		if(status /= nf90_NoErr) call handle_err(status)

		!!! Are any of our variable's dimensions longitude?
		status = nf90_inquire_variable(fid,vid,dimids=dimids)
		if(status /= nf90_NoErr) call handle_err(status)

		do idim=1,2
			status = nf90_inquire_dimension(fid,dimids(idim),name=dim_name)
			if(status /= nf90_NoErr) call handle_err(status)

			if ( TRIM(dim_name) == "longitude" ) then
				lon_idx = idim
				read_bounds(:,idim) = (/ fdim_j-dim_j_start+1, fdim_j-dim_j_start+2  /)
			else
				read_bounds(:,idim) = (/ size(out,dim=idim),1 /)
			endif
		end do

		if(PRESENT(starts)) then
			if(PRESENT(counts)) then
				if( lon_idx > 0 .and. starts(lon_idx) + counts(lon_idx) > fdim_j ) then
					temp_starts = starts
					temp_counts = counts
					temp_starts(lon_idx) = dim_j_start
					temp_counts(lon_idx) = fdim_j - dim_j_start + 1
					status = nf90_get_var(fid,vid,out(:read_bounds(1,1),:read_bounds(1,2)),start=temp_starts,count=temp_counts)
					if(status /= nf90_NoErr) call handle_err(status)
					temp_starts(lon_idx) = 1
					temp_counts(lon_idx) = dim_j_end
					status = nf90_get_var(fid,vid,out(read_bounds(2,1):,read_bounds(2,2):),start=temp_starts,count=temp_counts)
				else
					status = nf90_get_var(fid,vid,out,start=starts,count=counts)
				end if
			else
				status = nf90_get_var(fid,vid,out,start=starts)
			end if
		else
			if(PRESENT(counts)) then
				status = nf90_get_var(fid,vid,out,count=counts)
			else
				status = nf90_get_var(fid,vid,out)
			end if
		end if
		if(status /= nf90_NoErr) call handle_err(status)

		!!! Get scale_factor and add_offset (if any)
    status = nf90_get_att(fid, vid, "scale_factor", scale_factor)
    if( status /= nf90_NoErr ) then
      !!! OK if attribute not found - set it to 1.0
      if ( status /= nf90_eNotAtt ) call handle_err(status)
      scale_factor = 1.0
    endif

    status = nf90_get_att(fid, vid, "add_offset", add_offset)
    if( status /= nf90_NoErr ) then
      !!! OK if attribute not found - set it to 0.0
      if ( status /= nf90_eNotAtt ) call handle_err(status)
      add_offset = 0.0
    endif

		out = scale_factor * out + add_offset

	END SUBROUTINE

	SUBROUTINE get_era5_field_r3(d,m,y,field,out,starts,counts)

		USE netcdf
		USE global_data, ONLY: dirdata_era5, fdim_j, dim_j_start, dim_j_end
		USE util, ONLY: handle_err

		IMPLICIT NONE

		INTEGER, INTENT(IN) :: d,m,y
		CHARACTER(len=*), INTENT(IN) :: field
		REAL, DIMENSION(:,:,:) :: out
		INTEGER, OPTIONAL, DIMENSION(3), INTENT(IN) :: starts, counts

		!!! Locals
		CHARACTER(len=100) :: fn
		INTEGER :: fid, vid
		INTEGER :: status
		REAL    :: scale_factor, add_offset
		INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: dimids
		CHARACTER(len=NF90_MAX_NAME) :: dim_name
		INTEGER :: lon_idx = -1
		INTEGER :: idim

		INTEGER, DIMENSION(3) :: temp_starts
		INTEGER, DIMENSION(3) :: temp_counts
		INTEGER, DIMENSION(2,3) :: read_bounds !!! (/ (/ end1, start2 /), (/ end1, start2 /), (/ end1, start2 /) /)

		call get_filename(d,m,y,field,fn)

		status = NF90_OPEN(TRIM(fn),NF90_NOWRITE,fid)
		if (status /= NF90_NOERR) call handle_err(status)

		status = nf90_inq_varid(fid, trim(field), vid)
		if(status /= nf90_NoErr) call handle_err(status)

		!!! Are any of our variable's dimensions longitude?
		status = nf90_inquire_variable(fid,vid,dimids=dimids)
		if(status /= nf90_NoErr) call handle_err(status)

		do idim=1,3
			status = nf90_inquire_dimension(fid,dimids(idim),name=dim_name)
			if(status /= nf90_NoErr) call handle_err(status)

			if ( TRIM(dim_name) == "longitude" ) then
				lon_idx = idim
				read_bounds(:,idim) = (/ fdim_j-dim_j_start+1, fdim_j-dim_j_start+2  /)
			else
				read_bounds(:,idim) = (/ size(out,dim=idim),1 /)
			endif
		end do

		if(PRESENT(starts)) then
			if(PRESENT(counts)) then
				if( lon_idx > 0 .and. starts(lon_idx) + counts(lon_idx) > fdim_j ) then
					temp_starts = starts
					temp_counts = counts
					temp_starts(lon_idx) = dim_j_start
					temp_counts(lon_idx) = fdim_j - dim_j_start + 1
					status = nf90_get_var(fid,vid,out(:read_bounds(1,1),:read_bounds(1,2),:read_bounds(1,3)),start=temp_starts,count=temp_counts)
					if(status /= nf90_NoErr) call handle_err(status)
					temp_starts(lon_idx) = 1
					temp_counts(lon_idx) = dim_j_end
					status = nf90_get_var(fid,vid,out(read_bounds(2,1):,read_bounds(2,2):,read_bounds(2,3):),start=temp_starts,count=temp_counts)
				else
					status = nf90_get_var(fid,vid,out,start=starts,count=counts)
				end if
			else
				status = nf90_get_var(fid,vid,out,start=starts)
			end if
		else
			if(PRESENT(counts)) then
				status = nf90_get_var(fid,vid,out,count=counts)
			else
				status = nf90_get_var(fid,vid,out)
			end if
		end if
		if(status /= nf90_NoErr) call handle_err(status)

		!!! Get scale_factor and add_offset (if any)
		status = nf90_get_att(fid, vid, "scale_factor", scale_factor)
		if( status /= nf90_NoErr ) then
		  !!! OK if attribute not found - set it to 1.0
		  if ( status /= nf90_eNotAtt ) call handle_err(status)
		  scale_factor = 1.0
		endif

		status = nf90_get_att(fid, vid, "add_offset", add_offset)
		if( status /= nf90_NoErr ) then
		  !!! OK if attribute not found - set it to 0.0
		  if ( status /= nf90_eNotAtt ) call handle_err(status)
		  add_offset = 0.0
		endif

		out = scale_factor * out + add_offset

	END SUBROUTINE

	SUBROUTINE get_era5_field_r4(d,m,y,field,out,starts,counts)

		USE netcdf
		USE global_data, ONLY: dirdata_era5, fdim_j, dim_j_start, dim_j_end, dim_i, dim_j
		USE util, ONLY: handle_err

		IMPLICIT NONE

		INTEGER, INTENT(IN) :: d,m,y
		CHARACTER(len=*), INTENT(IN) :: field
		REAL, DIMENSION(:,:,:,:) :: out
		INTEGER, OPTIONAL, DIMENSION(4), INTENT(IN) :: starts, counts

		!!! Locals
		CHARACTER(len=100) :: fn
		INTEGER :: fid, vid
		INTEGER :: status
		REAL    :: scale_factor, add_offset
		INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: dimids
		CHARACTER(len=NF90_MAX_NAME) :: dim_name
		INTEGER :: lon_idx = -1
		INTEGER :: idim

		INTEGER, DIMENSION(4) :: temp_starts
		INTEGER, DIMENSION(4) :: temp_counts
		INTEGER, DIMENSION(2,4) :: read_bounds !!! (/ (/ end1, start2 /), (/ end1, start2 /), (/ end1, start2 /), (/ end1, start2 /) /)

		INTEGER :: latid,lonid,vlatid,vlonid, stat

		call get_filename(d,m,y,field,fn)

		status = NF90_OPEN(TRIM(fn),NF90_NOWRITE,fid)
		if (status /= NF90_NOERR) call handle_err(status)

		status = nf90_inq_varid(fid, trim(field), vid)
		if(status /= nf90_NoErr) call handle_err(status)

				!!! Are any of our variable's dimensions longitude?
		status = nf90_inquire_variable(fid,vid,dimids=dimids)
		if(status /= nf90_NoErr) call handle_err(status)

		do idim=1,4
			status = nf90_inquire_dimension(fid,dimids(idim),name=dim_name)
			if(status /= nf90_NoErr) call handle_err(status)

			if ( TRIM(dim_name) == "longitude" ) then
				lon_idx = idim
				read_bounds(:,idim) = (/ fdim_j-dim_j_start+1, fdim_j-dim_j_start+2  /)
			else
				read_bounds(:,idim) = (/ size(out,dim=idim),1 /)
			endif
		end do

		if(PRESENT(starts)) then
			if(PRESENT(counts)) then
				if( lon_idx > 0 .and. starts(lon_idx) + counts(lon_idx) > fdim_j ) then
					temp_starts = starts
					temp_counts = counts
					temp_starts(lon_idx) = dim_j_start
					temp_counts(lon_idx) = fdim_j - dim_j_start + 1
					status = nf90_get_var(fid,vid,out(:read_bounds(1,1),:read_bounds(1,2),:read_bounds(1,3),:read_bounds(1,4)),start=temp_starts,count=temp_counts)
					if(status /= nf90_NoErr) call handle_err(status)
					temp_starts(lon_idx) = 1
					temp_counts(lon_idx) = dim_j_end
					status = nf90_get_var(fid,vid,out(read_bounds(2,1):,read_bounds(2,2):,read_bounds(2,3):,read_bounds(2,4):),start=temp_starts,count=temp_counts)
				else
					status = nf90_get_var(fid,vid,out,start=starts,count=counts)
				end if
			else
				status = nf90_get_var(fid,vid,out,start=starts)
			end if
		else
			if(PRESENT(counts)) then
				status = nf90_get_var(fid,vid,out,count=counts)
			else
				status = nf90_get_var(fid,vid,out)
			end if
		end if
		if(status /= nf90_NoErr) call handle_err(status)

		!!! Get scale_factor and add_offset (if any)
    status = nf90_get_att(fid, vid, "scale_factor", scale_factor)
    if( status /= nf90_NoErr ) then
      !!! OK if attribute not found - set it to 1.0
      if ( status /= nf90_eNotAtt ) call handle_err(status)
      scale_factor = 1.0
    endif

    status = nf90_get_att(fid, vid, "add_offset", add_offset)
    if( status /= nf90_NoErr ) then
      !!! OK if attribute not found - set it to 0.0
      if ( status /= nf90_eNotAtt ) call handle_err(status)
      add_offset = 0.0
    endif

		out = scale_factor * out + add_offset

	END SUBROUTINE

	SUBROUTINE get_filename(d,mn,yr,field,fn)

		!-----------------------------------------------
		! given the month and year get the filename extension string
		!---------------------------------------------------

		USE global_data, ONLY: dirdata_era5
		USE util, ONLY: month_end, to_iso_date

		IMPLICIT NONE

		INTEGER, INTENT(IN) :: d
		INTEGER, INTENT(IN) :: mn, yr
		CHARACTER(LEN=*), INTENT(IN) :: field
		CHARACTER(LEN=100), INTENT(OUT) :: fn
		! len('YYYYMMDD-YYYYMMDD.nc') = 20
		CHARACTER(LEN=20) :: suffix
		CHARACTER(LEN=5)  :: era5_field

		suffix = to_iso_date(yr,mn,1)//"-"//to_iso_date(yr,mn,month_end(yr,mn))//".nc"
		era5_field = field

		SELECT CASE(trim(era5_field))
		CASE("tp","mslhf","sp","blh","tcw")
			write(fn,'(a,i4.4,a)') TRIM(dirdata_era5)//"single-levels/reanalysis/"//trim(era5_field)//"/",yr,"/"//trim(era5_field)//"_era5_oper_sfc_"//suffix   ! CHECKED
		CASE DEFAULT
			write(fn,'(a,i4.4,a)') TRIM(dirdata_era5)//"pressure-levels/reanalysis/"//trim(era5_field)//"/",yr,"/"//trim(era5_field)//"_era5_oper_pl_"//suffix   ! CHECKED
		fn = ADJUSTL(fn)
		END SELECT

	END SUBROUTINE get_filename

	SUBROUTINE get_grid_data(ptop, delx, datatstep, lat2d, lon2d, extents)

		USE global_data, ONLY: syear, smon, sday, dirdata_era5, totpts, bdy, datansteps, dim_i_start, dim_j_start, dim_k_start, dim_i, dim_j, dim_k, fdim_i, fdim_j, dim_i_end, dim_j_end, dim_k_end
		USE util, ONLY: int_to_string, handle_err, all_positive_longitude, to_iso_date, month_end, array_extents
		USE netcdf

		IMPLICIT NONE

		REAL,INTENT(OUT)                             :: ptop, delx
		INTEGER, INTENT(OUT)                         :: datatstep
		REAL,ALLOCATABLE,DIMENSION(:,:), INTENT(OUT) :: lat2d,lon2d
		REAL,DIMENSION(6), INTENT(IN),OPTIONAL       :: extents

		!!! Locals
		INTEGER :: status
		INTEGER :: headncid, tstepid, latcrsid, loncrsid, levelid
		INTEGER :: fdim_k
		INTEGER :: idim
		REAL, ALLOCATABLE, DIMENSION(:) :: lat1d, lon1d, levels

		REAL, DIMENSION(2)              :: input_timeseries
		CHARACTER(len=100) :: fname

		write(fname,'(a,i4.4,a)') TRIM(dirdata_era5)//"pressure-levels/reanalysis/q/",syear,"/q_era5_oper_pl_"//to_iso_date(syear,smon,1)//"-"//to_iso_date(syear,smon,month_end(syear,smon))//".nc"
		print *,'using header from',fname
		status = NF90_OPEN(fname, NF90_NOWRITE, headncid)
		if (status /= NF90_NOERR) call handle_err(status)

		! Read the length of each dimension of input data
		status = nf90_inq_dimid(headncid, "time", tstepid)  ! number of time steps in file
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_inq_dimid(headncid, "latitude", latcrsid)  !latitudes of grid points (degrees)
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_inq_dimid(headncid, "longitude", loncrsid)  !longitudes of grid points (degrees)
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_inq_dimid(headncid, "level", levelid)  ! pressure levels of dataset
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_inquire_dimension(headncid, tstepid, len = datansteps)
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_inquire_dimension(headncid, latcrsid, len = fdim_i)
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_inquire_dimension(headncid, loncrsid, len = fdim_j)
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_inquire_dimension(headncid, levelid, len = fdim_k)
		if (status /= NF90_NOERR) call handle_err(status)

		allocate(lat1d(fdim_i),lon1d(fdim_j),levels(fdim_k))

		status = nf90_inq_varid(headncid, "latitude", latcrsid)
		if(status /= nf90_NoErr) call handle_err(status)
		status = nf90_get_var(headncid, latcrsid, lat1d)
		if(status /= nf90_NoErr) call handle_err(status)

		status = nf90_inq_varid(headncid, "longitude", loncrsid)
		if(status /= nf90_NoErr) call handle_err(status)
		status = nf90_get_var(headncid, loncrsid, lon1d)
		if(status /= nf90_NoErr) call handle_err(status)

		status = nf90_inq_varid(headncid, "level", levelid)
		if(status /= nf90_NoErr) call handle_err(status)
		status = nf90_get_var(headncid, levelid, levels)
		if(status /= nf90_NoErr) call handle_err(status)

		status = nf90_inq_varid(headncid, "time", tstepid)
		if(status /= nf90_NoErr) call handle_err(status)

		! Time dimension - just take the difference between the first two values as the file timestep
		status = nf90_get_var(headncid, tstepid, input_timeseries, start=(/1/), count=(/2/))
		if(status /= nf90_NoErr) call handle_err(status)
  !!! Future TO-DO - This could be automated to read units of time from header file
		datatstep = 60*(input_timeseries(2) - input_timeseries(1)) ! model expects timestep in minutes


		if( present(extents) ) then
			call array_extents(lat1d, extents(1),extents(2),dim_i_start,dim_i_end,reverse=.true.)
			call array_extents(lon1d, extents(3),extents(4),dim_j_start,dim_j_end,periodic=.true.)
			call array_extents(levels,extents(5),extents(6),dim_k_start,dim_k_end)
			dim_i_start = dim_i_start + bdy
			dim_j_start = dim_j_start + bdy
			dim_i_end   = dim_i_end - bdy
			dim_j_end   = dim_j_end - bdy
		else
			!!! extents not present, must want global grid
			dim_i_end  = fdim_i - bdy
			dim_j_end  = fdim_j - bdy
			dim_k_end  = fdim_k
			dim_i_start = 1 + bdy
			dim_j_start = 1 + bdy
			dim_k_start = 1
		endif

		dim_i = dim_i_end - dim_i_start + 1
		!!! dim_j is longitude and is periodic
		if ( dim_j_end > dim_j_start ) then
			dim_j = dim_j_end - dim_j_start + 1
		else
			dim_j = fdim_j - dim_j_start + dim_j_end + 1
		end if
		dim_k = dim_k_end - dim_k_start + 1

		allocate(lat2d(dim_j,dim_i))
		allocate(lon2d(dim_j,dim_i))

		do idim = 1,dim_j
			lat2d(idim,:) = lat1d(dim_i_start:dim_i_end)
		end do

		!!! Handle periodicity in longitude
		if ( dim_j_start < dim_j_end ) then
			do idim = 1,dim_i
				lon2d(:,idim) = lon1d(dim_j_start:dim_j_end)
			end do
		else
			do idim = 1,dim_i
				lon2d(1:fdim_j-dim_j_start+1,idim) = lon1d(dim_j_start:fdim_j)
				lon2d(fdim_j-dim_j_start+2:dim_j,idim) = lon1d(1:dim_j_end)
			end do
		end if

		call all_positive_longitude(lon2d,lon2d)

		!! mBar -> Pa
		ptop = levels(dim_k_start)*100.0

		!!! TESTING - THIS NEEDS TO BE AUTOMATED
		! Calculation of ssdim requires delx of the grid.
		delx=25202.112430956095



print *, 'dim_k_start,dim_k_end,ptop',dim_k_start,dim_k_end,ptop

  

	END SUBROUTINE

	SUBROUTINE get_watershed(wsmask)

		USE global_data, ONLY: diri_era5, fwshed_era5,dim_i,dim_j,dim_i_start, dim_j_start, dim_j_end, fdim_j
		USE util, ONLY: handle_err
		USE netcdf

		IMPLICIT NONE

		INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: wsmask

		!!! Locals
		INTEGER :: wsncid, wsid, status
		CHARACTER(len=100) :: fname
		INTEGER :: stat

		fname=TRIM(diri_era5)//TRIM(fwshed_era5)

		print *,'using wshed from',fname
		ALLOCATE( wsmask(dim_j,dim_i), STAT = status )

		status = NF90_OPEN(fname, NF90_NOWRITE, wsncid)
		if (status /= NF90_NOERR) call handle_err(status)

		status = nf90_inq_varid(wsncid, "wsmask", wsid)  !watershed mask
		if(status /= nf90_NoErr) call handle_err(status)

		if ( dim_j_start < dim_j_end ) then
			status = nf90_get_var(wsncid, wsid, wsmask,start=(/dim_j_start, dim_i_start/), count=(/dim_j, dim_i/))
		else
			status = nf90_get_var(wsncid, wsid, wsmask(1:fdim_j-dim_j_start+1,:),start=(/dim_j_start, dim_i_start/), count=(/fdim_j-dim_j_start+1, dim_i/))
			if(status /= nf90_NoErr) call handle_err(status)
			status = nf90_get_var(wsncid, wsid, wsmask(fdim_j-dim_j_start+2:dim_j,:),start=(/1, dim_i_start/), count=(/dim_j_end, dim_i/))
		end if
  	if(status /= nf90_NoErr) call handle_err(status)

		status = nf90_close(wsncid)
		if(status /= nf90_NoErr) call handle_err(status)

	END SUBROUTINE

	!!! Future TO-DO - rewrite this and get_data to open each era5 file only once
	!!! Future TO-DO - rewrite this to use a rolling window to avoid re-reading over
	!!!                back trajectory days
	SUBROUTINE get_data_mixtot(qc,qt)

		USE global_data, ONLY: peak, datatotsteps, datadaysteps, dim_i, dim_j, dim_k, dim_i_start, dim_j_start, dim_k_start, day, mon, year, water_density, totbtadays, sday
		USE util, ONLY: julian, gregorian, month_end

		IMPLICIT NONE

		REAL, DIMENSION(:,:,:,:) :: qc,qt

		!!! Locals
		REAL, DIMENSION(SIZE(qt,1),SIZE(qt,2),SIZE(qt,3),SIZE(qt,4)) :: clw,rnw,snow,ice
		INTEGER :: jd_today, jd_before
		INTEGER :: new_y, new_m, new_d
		INTEGER :: i
		REAL    :: dayend

		ice = 9999.0

		if (peak) then
			dayend = day + 0.5
		else
			dayend = day
		end if

		if (totbtadays>1) then

			call get_era5_field(day, mon, year, "clwc",  clw(:,:,:,(datatotsteps-datadaysteps):(datatotsteps-1)), starts=(/dim_j_start, dim_i_start, dim_k_start,(day-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,datadaysteps/))
			call get_era5_field(day, mon, year, "crwc",  rnw(:,:,:,(datatotsteps-datadaysteps):(datatotsteps-1)), starts=(/dim_j_start, dim_i_start, dim_k_start,(day-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,datadaysteps/))
			call get_era5_field(day, mon, year, "cswc", snow(:,:,:,(datatotsteps-datadaysteps):(datatotsteps-1)), starts=(/dim_j_start, dim_i_start, dim_k_start,(day-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,datadaysteps/))
			call get_era5_field(day, mon, year, "ciwc",  ice(:,:,:,(datatotsteps-datadaysteps):(datatotsteps-1)), starts=(/dim_j_start, dim_i_start, dim_k_start,(day-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,datadaysteps/))

			jd_today = julian(year,mon,day)

			do i = 1,totbtadays
				jd_before = jd_today-i
				call gregorian(jd_before,new_y,new_m,new_d)

				call get_era5_field(new_d, new_m, new_y, "clwc",  clw(:,:,:,(datatotsteps-datadaysteps*(i+1)):(datatotsteps-(datadaysteps*i)-1)), starts=(/dim_j_start, dim_i_start, dim_k_start,(new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,datadaysteps/))
				call get_era5_field(new_d, new_m, new_y, "crwc",  rnw(:,:,:,(datatotsteps-datadaysteps*(i+1)):(datatotsteps-(datadaysteps*i)-1)), starts=(/dim_j_start, dim_i_start, dim_k_start,(new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,datadaysteps/))
				call get_era5_field(new_d, new_m, new_y, "cswc", snow(:,:,:,(datatotsteps-datadaysteps*(i+1)):(datatotsteps-(datadaysteps*i)-1)), starts=(/dim_j_start, dim_i_start, dim_k_start,(new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,datadaysteps/))
				call get_era5_field(new_d, new_m, new_y, "ciwc",  ice(:,:,:,(datatotsteps-datadaysteps*(i+1)):(datatotsteps-(datadaysteps*i)-1)), starts=(/dim_j_start, dim_i_start, dim_k_start,(new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,datadaysteps/))
			end do

			jd_before = jd_today+1
			call gregorian(jd_before,new_y,new_m,new_d)

			call get_era5_field(new_d, new_m, new_y, "clwc",  clw(:,:,:,datatotsteps:datatotsteps), starts=(/dim_j_start, dim_i_start, dim_k_start,(new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,1/))
			call get_era5_field(new_d, new_m, new_y, "crwc",  rnw(:,:,:,datatotsteps:datatotsteps), starts=(/dim_j_start, dim_i_start, dim_k_start,(new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,1/))
			call get_era5_field(new_d, new_m, new_y, "cswc", snow(:,:,:,datatotsteps:datatotsteps), starts=(/dim_j_start, dim_i_start, dim_k_start,(new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,1/))
			call get_era5_field(new_d, new_m, new_y, "ciwc",  ice(:,:,:,datatotsteps:datatotsteps), starts=(/dim_j_start, dim_i_start, dim_k_start,(new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,1/))

			print *,'Input mixtot file of next day (1st time step) loaded successfully'
		else
			print*, 'If you only want to back-track for one day, must change how input data is retrieved.'
		end if

		qc = clw + rnw + snow + ice
		qt = qc

		print *, 'finished getting data mixtot'
		!print *, 'ice(1,1,1,:)', ice(1,1,1,:)

	END SUBROUTINE

	SUBROUTINE get_data(precip,evap,u,v,w,q,qc,qt,pres,pbl_hgt,psfc,tpw,cpre)

		USE global_data, ONLY: day, mon, year, dim_i, dim_j, dim_k, dim_i_start, dim_j_start, dim_k_start, sday, datadaysteps, totbtadays, datatotsteps, sday, syear, smon, peak, water_density, Lv, dirdata_era5
		USE util, ONLY: julian, gregorian, to_iso_date, month_end, handle_err
		USE netcdf

		IMPLICIT NONE

		REAL, DIMENSION(:,:,:) :: precip,evap,pbl_hgt,psfc,tpw,cpre
		REAL, DIMENSION(:,:,:,:) :: u,v,w,q,qc,qt,pres
		!REAL, DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3),datadaysteps) :: temp


		!!! Locals
		INTEGER :: jd_today, jd_before, new_y, new_m, new_d
		INTEGER :: i, ik, it
		INTEGER :: status, headncid, levelid, fdim_k
		REAL :: dayend
		CHARACTER(len=100) :: fname
		REAL, ALLOCATABLE, DIMENSION(:) :: levels

		call get_data_mixtot(qc,qt)

		!if this is a day around a storm peak we want the half day after as well
		if (peak) then
		  dayend = day + 0.5
		else
		  dayend = day
		end if

		print *, 'up to start of loading get data'

		if (totbtadays>1) then

			call get_era5_field(day, mon, year, "tp", precip, starts=(/dim_j_start, dim_i_start, (day-1)*datadaysteps+1/), counts=(/dim_j,dim_i,datadaysteps/))

			call get_era5_field(day, mon, year,  "e",   evap(:,:,(datatotsteps-datadaysteps):(datatotsteps-1)), starts=(/dim_j_start, dim_i_start, (day-1)*datadaysteps+1/), counts=(/dim_j,dim_i,datadaysteps/))
			call get_era5_field(day, mon, year,"sp",   psfc(:,:,(datatotsteps-datadaysteps):(datatotsteps-1)), starts=(/dim_j_start, dim_i_start, (day-1)*datadaysteps+1/), counts=(/dim_j,dim_i,datadaysteps/))
			call get_era5_field(day, mon, year,"blh",pbl_hgt(:,:,(datatotsteps-datadaysteps):(datatotsteps-1)), starts=(/dim_j_start, dim_i_start, (day-1)*datadaysteps+1/), counts=(/dim_j,dim_i,datadaysteps/))
			call get_era5_field(day, mon, year, "tcw",    tpw(:,:,(datatotsteps-datadaysteps):(datatotsteps-1)), starts=(/dim_j_start, dim_i_start, (day-1)*datadaysteps+1/), counts=(/dim_j,dim_i,datadaysteps/))
			call get_era5_field(day, mon, year,  "cp",   cpre(:,:,(datatotsteps-datadaysteps):(datatotsteps-1)), starts=(/dim_j_start, dim_i_start, (day-1)*datadaysteps+1/), counts=(/dim_j,dim_i,datadaysteps/))

			call get_era5_field(day, mon, year,"q",q(:,:,:,(datatotsteps-datadaysteps):(datatotsteps-1)), starts=(/dim_j_start, dim_i_start,dim_k_start,(day-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,datadaysteps/))
			call get_era5_field(day, mon, year,"u",u(:,:,:,(datatotsteps-datadaysteps):(datatotsteps-1)), starts=(/dim_j_start, dim_i_start,dim_k_start,(day-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,datadaysteps/))
			call get_era5_field(day, mon, year,"v",v(:,:,:,(datatotsteps-datadaysteps):(datatotsteps-1)), starts=(/dim_j_start, dim_i_start,dim_k_start,(day-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,datadaysteps/))
			call get_era5_field(day, mon, year,"w",w(:,:,:,(datatotsteps-datadaysteps):(datatotsteps-1)), starts=(/dim_j_start, dim_i_start,dim_k_start,(day-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,datadaysteps/))


			! Get julian day of current day
			jd_today = julian(year,mon,day)
			print *,'jd_today=',jd_today

			do i = 1,totbtadays
				jd_before = jd_today-i
				call gregorian(jd_before,new_y,new_m,new_d)

				call get_era5_field(new_d, new_m, new_y,   "e",   evap(:,:,(datatotsteps-(datadaysteps*(i+1))):(datatotsteps-(datadaysteps*i)-1)), starts=(/dim_j_start, dim_i_start, (new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,datadaysteps/))
				call get_era5_field(new_d, new_m, new_y, "sp",   psfc(:,:,(datatotsteps-(datadaysteps*(i+1))):(datatotsteps-(datadaysteps*i)-1)), starts=(/dim_j_start, dim_i_start, (new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,datadaysteps/))
				call get_era5_field(new_d, new_m, new_y, "blh",pbl_hgt(:,:,(datatotsteps-(datadaysteps*(i+1))):(datatotsteps-(datadaysteps*i)-1)), starts=(/dim_j_start, dim_i_start, (new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,datadaysteps/))
				call get_era5_field(new_d, new_m, new_y,  "tcw",    tpw(:,:,(datatotsteps-(datadaysteps*(i+1))):(datatotsteps-(datadaysteps*i)-1)), starts=(/dim_j_start, dim_i_start, (new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,datadaysteps/))
				call get_era5_field(new_d, new_m, new_y,  "cp",    cpre(:,:,(datatotsteps-(datadaysteps*(i+1))):(datatotsteps-(datadaysteps*i)-1)), starts=(/dim_j_start, dim_i_start, (new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,datadaysteps/))

				call get_era5_field(new_d, new_m, new_y,"q",q(:,:,:,(datatotsteps-(datadaysteps*(i+1))):(datatotsteps-(datadaysteps*i)-1)), starts=(/dim_j_start, dim_i_start,dim_k_start,(new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,datadaysteps/))
				call get_era5_field(new_d, new_m, new_y,"u",u(:,:,:,(datatotsteps-(datadaysteps*(i+1))):(datatotsteps-(datadaysteps*i)-1)), starts=(/dim_j_start, dim_i_start,dim_k_start,(new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,datadaysteps/))
				call get_era5_field(new_d, new_m, new_y,"v",v(:,:,:,(datatotsteps-(datadaysteps*(i+1))):(datatotsteps-(datadaysteps*i)-1)), starts=(/dim_j_start, dim_i_start,dim_k_start,(new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,datadaysteps/))
				call get_era5_field(new_d, new_m, new_y,"w",w(:,:,:,(datatotsteps-(datadaysteps*(i+1))):(datatotsteps-(datadaysteps*i)-1)), starts=(/dim_j_start, dim_i_start,dim_k_start,(new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,datadaysteps/))


				print *,'Input files of days to be back-tracked loaded successfully'
			end do

			jd_before = jd_today+1
			call gregorian(jd_before,new_y,new_m,new_d)

			call get_era5_field(new_d, new_m, new_y,   "e",   evap(:,:,datatotsteps:datatotsteps), starts=(/dim_j_start, dim_i_start, (new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,1/))
			call get_era5_field(new_d, new_m, new_y, "sp",   psfc(:,:,datatotsteps:datatotsteps), starts=(/dim_j_start, dim_i_start, (new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,1/))
			call get_era5_field(new_d, new_m, new_y, "blh",pbl_hgt(:,:,datatotsteps:datatotsteps), starts=(/dim_j_start, dim_i_start, (new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,1/))
			call get_era5_field(new_d, new_m, new_y, "tcw" ,    tpw(:,:,datatotsteps:datatotsteps), starts=(/dim_j_start, dim_i_start, (new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,1/))
			call get_era5_field(new_d, new_m, new_y, "cp" ,    cpre(:,:,datatotsteps:datatotsteps), starts=(/dim_j_start, dim_i_start, (new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,1/))

			call get_era5_field(new_d, new_m, new_y,"q",q(:,:,:,datatotsteps:datatotsteps), starts=(/dim_j_start, dim_i_start,dim_k_start,(new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,1/))
			call get_era5_field(new_d, new_m, new_y,"u",u(:,:,:,datatotsteps:datatotsteps), starts=(/dim_j_start, dim_i_start,dim_k_start,(new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,1/))
			call get_era5_field(new_d, new_m, new_y,"v",v(:,:,:,datatotsteps:datatotsteps), starts=(/dim_j_start, dim_i_start,dim_k_start,(new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,1/))
			call get_era5_field(new_d, new_m, new_y,"w",w(:,:,:,datatotsteps:datatotsteps), starts=(/dim_j_start, dim_i_start,dim_k_start,(new_d-1)*datadaysteps+1/), counts=(/dim_j,dim_i,dim_k,1/))


			print *,'Input file of next day (1st time step) loaded successfully'
		else
			print*, 'If you only want to back-track for one day, must change how input data is retrieved.'
		end if


		qt = qt + q ! i.e. SUM(QCLD,QRAIN,QSNOW,QICE) + QVAPOUR


		!!! Pressure is special, derive it from a coordinate
		write(fname,'(a,i4.4,a)') TRIM(dirdata_era5)//"pressure-levels/reanalysis/q/",syear,"/q_era5_oper_pl_"//to_iso_date(syear,smon,1)//"-"//to_iso_date(syear,smon,month_end(syear,smon))//".nc"
		print *,'using header from',fname
		status = NF90_OPEN(fname, NF90_NOWRITE, headncid)
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_inq_dimid(headncid, "level", levelid)  ! pressure levels of dataset
		if (status /= NF90_NOERR) call handle_err(status)
		status = nf90_inquire_dimension(headncid, levelid, len = fdim_k)
		if (status /= NF90_NOERR) call handle_err(status)
		allocate(levels(fdim_k))
		status = nf90_inq_varid(headncid, "level", levelid)
		if(status /= nf90_NoErr) call handle_err(status)
		status = nf90_get_var(headncid, levelid, levels)
		if(status /= nf90_NoErr) call handle_err(status)
		status = nf90_close(headncid)
		if(status /= nf90_NoErr) call handle_err(status)

	print *,'dim_k_start',dim_k_start
	print *,'levels',levels
	print *,'levels(1)',levels(1)
	print *,'levels(dim_k)',levels(dim_k)
  

		do it = 1,datatotsteps
			do ik = 1,dim_k
				!! mBar -> Pa
                pres(:,:,ik,it) = levels(dim_k_start - 1 + ik)*100.0
			end do
		end do

  	print *, 'pres(1,1,:,1)', pres(1,1,:,1)



	END SUBROUTINE

END MODULE input_data_handling_era5