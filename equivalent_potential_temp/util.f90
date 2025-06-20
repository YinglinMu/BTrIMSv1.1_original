MODULE util

	IMPLICIT NONE

	CONTAINS

	SUBROUTINE handle_err(status)
	!---------------------------------
	! handle any errors from the fortran 90 netCDF interface
	!---------------------------------

		USE netcdf

		IMPLICIT NONE

		integer, intent (in) :: status

		print *,'status = ',status

		if(status /= nf90_noerr) then
			print *, trim(nf90_strerror(status))
			stop "Stopped"
		end if

	END SUBROUTINE handle_err

	!***********************************************************************

	INTEGER FUNCTION string_to_int(string)
	!------------------------------------------
	! converts a character string to an integer
	!--------------------------------------------

		IMPLICIT NONE

		character (len=*), intent(in) :: string

		! local constant
		integer, parameter :: zero = iachar("0")
		integer :: i,  sign, integ
		character (len=50) :: str

		str = trim(string)
		integ = 0

		select case (str(1:1))
		case ("-")
		  sign = -1
		  str = str(2:)
		case ("+")
		  sign = 1
		  str = str(2:)
		case ("0":"9")
		  sign = 1
		end select

		do i=len(trim(str)),1,-1
			select case (str(i:i))
				case ("0":"9")
					integ = integ + (iachar(string(i:i))-zero)*10**(len(trim(str))-i)
				case default
					print *, "cannot convert a non-integer character to an integer!!"
					return
			end select
		end do

		string_to_int = integ

	end FUNCTION string_to_int

	!***********************************************************************

	CHARACTER(LEN=50) FUNCTION int_to_string(integ)
	!----------------------------------------------
	! converts an integer to a character string
	!----------------------------------------------

		IMPLICIT NONE

		integer, intent(in) :: integ

		! local constant
		integer, parameter :: zero = iachar("0")
		character (len=50) :: str
		integer :: i, inte

		str="                                                  "
		inte = integ

		do i=1,50
			str(50-i:50-i) = achar(mod(abs(inte),10)+zero)
			inte = int(inte/10)
			if (abs(inte)<1) exit
		end do

		if (integ<0) str(50-i-1:50-i) = "-"

		int_to_string = adjustl(str)

	end FUNCTION int_to_string

	!***********************************************************************

	CHARACTER(LEN=50) FUNCTION real_to_string(num)
	!----------------------------------------------
	! converts an real to a character string
	! here I allow a maximum of 10 decimal places
	!----------------------------------------------

		IMPLICIT NONE

		!real, intent(in) :: num
		INTEGER, intent(in) :: num

		! local
		integer :: i, whole, frac
		real :: fracnum


		whole = num !AINT(num)
		fracnum = num - whole

		do i=1,10
			if (MOD(fracnum,1.)==0) exit
			fracnum = fracnum * 10.
		end do

		frac = AINT(fracnum)

		real_to_string = TRIM(int_to_string(whole))//"."//TRIM(int_to_string(frac))

	end FUNCTION real_to_string

	!***********************************************************************

	INTEGER FUNCTION julian(year,mon,day)

	! From http://aa.usno.navy.mil/faq/docs/JD_Formula.php

		IMPLICIT NONE
		INTEGER, INTENT(IN) :: year,mon,day

		julian= day-32075+1461*(year+4800+(mon-14)/12)/4+367*(mon-2-(mon-14)/12*12)/12-3*((year+4900+(mon-14)/12)/100)/4

	END FUNCTION julian

		!***********************************************************************

	SUBROUTINE GREGORIAN(JD,YEAR,MONTH,DAY)

	! From http://aa.usno.navy.mil/faq/docs/JD_Formula.php

		IMPLICIT NONE
		!REAL, INTENT(IN) :: JD
		INTEGER :: JD
		INTEGER, INTENT(OUT) :: YEAR, MONTH, DAY!, HOUR, MINUTE, SECOND
		REAL :: JT
		INTEGER :: I,J,K,L,N

		L = INT(JD)+68569!JD= K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12)/12-3*((I+4900+(J-14)/12)/100)/4
		N = 4*L/146097
		L = L-(146097*N+3)/4
		I = 4000*(L+1)/1461001
		L = L-1461*I/4+31
		J = 80*L/2447
		K = L-2447*J/80
		L = J/11
		J = J+2-12*L
		I = 100*(N-49)+I+L

		YEAR = I
		MONTH = J
		DAY = K

	END SUBROUTINE GREGORIAN

	!***********************************************************************

	INTEGER FUNCTION simlength(startday,startmon,startyear,endday,endmon,endyear)
	!-----------------------------------------------------
	! given the start and end dates of the desired simulation
	! period (top of global data), calculate the number of days
	! in the period

	! Functions: don't need to specify what comes out, e.g. fn_name(in,in,in)
	! Subroutines: do specify in and out, e.g. sbrtn_name(in,in,in,out,out)

	!-----------------------------------------------------
		USE global_data

		IMPLICIT NONE

		INTEGER, intent(in) :: startday,startmon,startyear,endday,endmon,endyear
		INTEGER :: start_jd, end_jd

		start_jd = julian(startyear,startmon,startday)
		end_jd = julian(endyear,endmon,endday)
		simlength = end_jd - start_jd

	END FUNCTION simlength

	!***********************************************************************

	SUBROUTINE day_month_year(simday)
	!----------------------------------------------
	! given the simulation day, calculate the corresponding month and year
	! simulation day one is the first day
	!-----------------------------------------------------

		USE global_data

		IMPLICIT NONE

		INTEGER,INTENT(IN) :: simday	!simulation day

		INTEGER :: jday,simdayyear,simdaymonth,simdayday

		if (simday==1) then
			year=syear
			mon=smon
			day=sday
		else
			! Find julian day of simday, being the start day + an increment of days
			jday=julian(syear,smon,sday)+(simday-1)
			! Convert the julian day to a gregorian day
		call gregorian(jday,simdayyear,simdaymonth,simdayday)
			year=simdayyear
			mon=simdaymonth
			day=simdayday
		end if

	END SUBROUTINE day_month_year

	!***********************************************************************

	SUBROUTINE all_positive_longitude(lon2d,lon2d_corrected)
		!------------------------------------------------
		! ERA5 longitude 2d grid is negative east of dateline.
		! This subroutine converts the negative values to positive.
		!------------------------------------------------

			REAL, DIMENSION(:,:) :: lon2d,lon2d_corrected

			lon2d_corrected=lon2d

			WHERE (lon2d < 0)
				lon2d_corrected = lon2d+360
			END WHERE

		END SUBROUTINE all_positive_longitude

	!***********************************************************************

	FUNCTION to_iso_date(y,m,d)
	!----------------------------------------------------
	! Convert year month day integers to iso date string
	!----------------------------------------------------

		INTEGER, INTENT(IN) :: y,m,d
		CHARACTER(len=8)    :: to_iso_date

		write(to_iso_date,'(i4.4,i2.2,i2.2)') y,m,d

	END FUNCTION to_iso_date

	FUNCTION month_end(y,m)
	!----------------------------------------------------
	! Lookup table for ends of months
	!----------------------------------------------------
		INTEGER, INTENT(IN) :: y,m
		INTEGER :: month_end

		SELECT CASE	(m)
			CASE(1,3,5,7,8,10,12)
				month_end = 31
				RETURN
			CASE(4,6,9,11)
				month_end = 30
				RETURN
		END SELECT

		IF ( MODULO(y,4) == 0 ) then
			IF ( MODULO(y,100) == 0 .and. MODULO(y,400) /= 0 ) then
				month_end = 28
				RETURN
			ELSE
				month_end = 29
				RETURN
			END IF
		ELSE
			month_end = 28
			RETURN
		END IF

	END FUNCTION month_end

	SUBROUTINE array_extents(arr,in_start,in_end,i_start,i_end,reverse,periodic)

		USE global_data, ONLY: delta_coord

		IMPLICIT NONE

		REAL, DIMENSION(:), INTENT(IN) :: arr
		REAL, INTENT(IN)               :: in_start, in_end
		INTEGER, INTENT(OUT)           :: i_start, i_end
		LOGICAL,OPTIONAL,INTENT(IN)    :: reverse
		LOGICAL,OPTIONAL,INTENT(IN)    :: periodic

		!!! Locals
		INTEGER :: i
		REAL :: start, end

		i_start = -1
		i_end   = -1
		if( .not.present(periodic) .or. .not.periodic ) then
			if ( in_start > in_end ) then
				!!! Swap bounds if necessary
				end   = in_start
				start = in_end
			else
				start = in_start
				end   = in_end
			end if
		else
			start = in_start
			end   = in_end
		end if


		if( present(reverse) .and. reverse ) then
			do i=1,SIZE(arr)
				if ( i_start == -1 ) then
					if ( abs(arr(i) - end) < delta_coord ) i_start = i
				else if ( i_end == -1 ) then
					if ( abs(arr(i) - start) < delta_coord ) i_end = i
				else
					exit
				end if
			end do
		else
			do i=1,SIZE(arr)
				if ( i_start == -1 ) then
					if ( abs(arr(i) - start) < delta_coord ) i_start = i
				else if ( i_end == -1 ) then
					if ( abs(arr(i) - end) < delta_coord ) i_end = i
				else
					exit
				end if
			end do
		end if

		if ( i_start == -1 ) i_start = 1
		if ( i_end == -1 ) then
			if ( present(periodic) .and. periodic ) then
				!!! If we haven't found and 'end', check if we need to wrap around
				if ( end < start ) then
					!!! Redo the loop from the start
					do i=1,SIZE(arr)
						if ( i_end == -1 ) then
							if ( abs(arr(i) - end) < delta_coord ) i_end = i
						else
							exit
						end if
					end do
				else
					i_end = size(arr)
				end if
			else
				i_end = size(arr)
			end if
		end if

	END SUBROUTINE array_extents

END MODULE util