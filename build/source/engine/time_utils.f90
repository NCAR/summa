! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module time_utils_module

! data types
USE nrtype

! model constants
USE multiconst,only:secprday,secprhour,secprmin  ! seconds in an (day, hour, minute)

! privacy
implicit none
private
public::extractTime
public::compjulday
public::compcalday
public::elapsedSec
public::fracDay
contains


 ! ******************************************************************************************
 ! public subroutine extractTime: extract year/month/day/hour/minute/second from units string
 ! as well as any time zone information (hour/minute/second)
 ! ******************************************************************************************
 subroutine extractTime(refdate,iyyy,im,id,ih,imin,dsec,ih_tz,imin_tz,dsec_tz,err,message)
 implicit none
 ! dummy variables
 character(*),intent(in)    :: refdate             ! units string (time since...)
 integer(i4b),intent(out)   :: iyyy,im,id,ih,imin  ! time (year/month/day/hour/minute)
 real(rkind),intent(out)       :: dsec                ! seconds
 integer(i4b),intent(out)   :: ih_tz,imin_tz       ! time zone information (hour/minute)
 real(rkind),intent(out)       :: dsec_tz             ! time zone information (seconds)
 integer(i4b),intent(out)   :: err                 ! error code
 character(*),intent(out)   :: message             ! error message
 ! local variables
 integer(i4b)               :: n,nsub              ! length of the string and substring
 integer(i4b)               :: istart              ! position in string
 ! initialize error control
 err=0; message="extractTime/"

 ! There are up to three space delimited fields: iyyy-im-id ih:imin:dsec ih_tz:imin_tz
 ! we'll parse each of these in order.

 ! Missing ih, imin, dsec, ih_tz, imin_tz and dsec_tz fields will be set to zero without causing an error.
 ih=0; imin=0; dsec=0._rkind; ih_tz=0; imin_tz=0; dsec_tz=0._rkind;

 ! get the length of the string
 n = len_trim(refdate)

 ! FIELD 1: move to a position in string past the time units (seconds since , days since , hours since )
 istart = index(refdate,'since')  ! get the index at the beginning of the word "since"
 if (istart>0) then ! if the word "since" exists
  istart = istart + index(refdate(istart:n)," ")
 else
  istart=1
 end if

 ! eat all the whitespace at the start of the string
 do while (refdate(istart:istart)==" ")
  istart = istart+1
 end do

 ! Find the end of FIELD 1
 nsub = index(refdate(istart:n)," ")
 if(nsub==0)then
  nsub=n             ! not found - read till the end of string
 else
  nsub=nsub+istart-1 ! found - read till the character before the whitespace
 end if

 ! parse the iyyy-im-id string
 call extract_dmy(refdate(istart:nsub),"-",iyyy,im,id,err,message); if (err/=0) return
 if(iyyy < 1900)then; err=20; message=trim(message)//'year < 1900'; return; end if
 if(iyyy > 2100)then; err=20; message=trim(message)//'year > 2100'; return; end if
 if(im   <    1)then; err=20; message=trim(message)//'month < 1';   return; end if
 if(im   >   12)then; err=20; message=trim(message)//'month > 12';  return; end if
 if(id   <    1)then; err=20; message=trim(message)//'day < 1';     return; end if
 if(id   >   31)then; err=20; message=trim(message)//'day > 31';    return; end if

 ! FIELD 2: Advance to the ih:imin:dsec string
 istart=nsub+1
 if(istart>n) return  ! no more time info

 ! eat all the whitespace at the start of the string
 do while (refdate(istart:istart)==" ")
  istart = istart+1
 end do

 if(istart>n) return  ! no more time info

 ! Find the end of FIELD 2
 nsub = index(refdate(istart:n)," ")
 if(nsub==0)then
  nsub=n             ! not found - read till the end of string
 else
  nsub=nsub+istart-1 ! found - read till the character before the whitespace
 end if

 ! parse the ih:imin:dsec string. Anything other than ih is optional
 call extract_hms(refdate(istart:nsub),":",ih,imin,dsec,err,message); if (err/=0) return
 if(ih   <  0)    then; err=20; message=trim(message)//'hour < 0';    return; end if
 if(ih   > 24)    then; err=20; message=trim(message)//'hour > 24';   return; end if
 if(imin <  0)    then; err=20; message=trim(message)//'minute < 0';  return; end if
 if(imin > 60)    then; err=20; message=trim(message)//'minute > 60'; return; end if
 if(dsec <  0._rkind)then; err=20; message=trim(message)//'second < 0';  return; end if
 if(dsec > 60._rkind)then; err=20; message=trim(message)//'second > 60'; return; end if

 ! FIELD 3: Advance to the ih_tz:imin_tz string
 istart=nsub+1
 if(istart>n) return  ! no more time info

 ! eat all the whitespace at the start of the string
 do while (refdate(istart:istart)==" ")
  istart = istart+1
 end do

 if(istart>n) return  ! no more time info

 ! Find the end of FIELD 3
 nsub = index(refdate(istart:n)," ")
 if(nsub==0)then
  nsub=n             ! not found - read till the end of string
 else
  nsub=nsub+istart-1 ! found - read till the character before the whitespace
 end if

 ! parse the ih_tz:imin_tz string.
 call extract_hms(refdate(istart:nsub),":",ih_tz,imin_tz,dsec_tz,err,message); if (err/=0) return
 if(ih_tz   <  -12)    then; err=20; message=trim(message)//'time zone hour < -12';  return; end if
 if(ih_tz   >   12)    then; err=20; message=trim(message)//'time zone hour > 12';   return; end if
 if(imin_tz <    0)    then; err=20; message=trim(message)//'time zone minute < 0';  return; end if
 if(imin_tz >   60)    then; err=20; message=trim(message)//'time zone minute > 60'; return; end if
 if(dsec_tz <    0._rkind)then; err=20; message=trim(message)//'time zone second < 0';  return; end if
 if(dsec_tz >   60._rkind)then; err=20; message=trim(message)//'time zone second > 60'; return; end if

 contains


  ! ******************************************************************************************
  ! internal subroutine extract: extract substring
  ! ******************************************************************************************
  subroutine extract(substring,cdelim,iend,itemp,err,message)
  implicit none
  ! input
  character(*),intent(in)     :: substring  ! sub-string to process
  character(len=1),intent(in) :: cdelim     ! string delimiter
  ! output
  integer(i4b),intent(out)    :: iend       ! index at the end of desired string
  integer(i4b),intent(out)    :: itemp      ! output date
  integer(i4b),intent(out)    :: err        ! error code
  character(*),intent(out)    :: message    ! error message
  ! initialize error code and message
  err=0; message="extract/"
  ! identify end-point of string
  iend = index(substring,cdelim)
  ! if sub-string does not exist, assume end is at end of string
  if (iend==0) iend=len_trim(substring)+1
  ! convert string to integer
  read(substring(1:iend-1),*,iostat=err) itemp
  ! read error
  if (err/=0) then
   err=20; message=trim(message)//"unexpected characters [string='"//trim(substring)//"']"; return
  end if
  end subroutine extract


  ! ******************************************************************************************
  ! internal subroutine extract_dmy: extract iyyy-im-id
  ! ******************************************************************************************
  subroutine extract_dmy(substring,cdelim,iyyy,im,id,err,message)
  implicit none
  ! input
  character(*),intent(in)     :: substring       ! sub-string to process
  character(len=1),intent(in) :: cdelim          ! string delimiter
  ! output
  integer(i4b),intent(out)    :: iyyy            ! year
  integer(i4b),intent(out)    :: im              ! month
  integer(i4b),intent(out)    :: id              ! day
  integer(i4b),intent(out)    :: err             ! error code
  character(*),intent(out)    :: message         ! error message
  ! local variables
  integer(i4b)                :: istart,iend,n   ! position in string

  ! initialize error code and message
  err=0; message="extract_dmy/"

  ! start indices
  n = len_trim(substring)
  istart = 1

  ! extract the year
  call extract(substring(istart:n),cdelim,iend,iyyy,err,message); if (err/=0) return

  ! extract the month
  istart=istart+iend
  call extract(substring(istart:n),cdelim,iend,im,err,message); if (err/=0) return

  ! extract the day
  istart=istart+iend
  call extract(substring(istart:n),cdelim,iend,id,err,message); if (err/=0) return

  end subroutine extract_dmy

  ! ******************************************************************************************
  ! internal subroutine extract_hms: extract hh:mm:sec
  ! ******************************************************************************************
  subroutine extract_hms(substring,cdelim,hh,mm,ss,err,message)
  implicit none
  ! input
  character(*),intent(in)     :: substring      ! sub-string to process
  character(len=1),intent(in) :: cdelim         ! string delimiter
  ! output
  integer(i4b),intent(out)    :: hh             ! hour
  integer(i4b),intent(out)    :: mm             ! minute
  real(rkind)    ,intent(out)    :: ss             ! sec
  integer(i4b),intent(out)    :: err            ! error code
  character(*),intent(out)    :: message        ! error message
  ! local variables
  integer(i4b)                :: istart,iend,n  ! position in string

  ! initialize error code and message
  err=0; message="extract_hms/"

  ! initialize indices
  n = len_trim(substring)
  istart = 1

  ! extract the hour
  call extract(substring(istart:n),cdelim,iend,hh,err,message); if (err/=0) return

  ! extract the minute
  istart=istart+iend
  if(istart > n) return
  call extract(substring(istart:n),cdelim,iend,mm,err,message); if (err/=0) return

  ! extract the second
  istart=istart+iend
  if(istart > n) return
  read(substring(istart:n),*) ss

  end subroutine extract_hms

 end subroutine extractTime



 ! ***************************************************************************************
 ! public subroutine compjulday: convert date to julian day (units of days)
 ! ***************************************************************************************
 subroutine compjulday(iyyy,mm,id,ih,imin,dsec,&  ! input
                       juldayss,err,message)      ! output
 implicit none
 ! input variables
 integer(i4b),intent(in)   :: iyyy,mm,id   ! year, month, day
 integer(i4b),intent(in)   :: ih,imin      ! hour, minute
 real(rkind),intent(in)       :: dsec         ! seconds
 ! output
 real(rkind),intent(out)      :: juldayss
  integer(i4b),intent(out) :: err          ! error code
  character(*),intent(out) :: message      ! error message
 ! local variables
 integer(i4b)              :: julday       ! julian day
 integer(i4b),parameter    :: igreg=15+31*(10+12*1582)  !IGREG = 588829
 integer(i4b)              :: ja,jm,jy
 real(rkind)                  :: jfrac        ! fraction of julian day

 ! initialize errors
 err=0; message="juldayss"

 ! compute julian day
 jy=iyyy
 if (jy.eq.0) then; err=10; message=trim(message)//"noYearZero/"; return; end if
 if (jy.lt.0) jy=jy+1
 if (mm.gt.2) then
  jm=mm+1
 else
  jy=jy-1
  jm=mm+13
 end if
 julday=int(365.25*jy)+int(30.6001*jm)+id+1720995
 if (id+31*(mm+12*iyyy).ge.IGREG) then
  ja=int(0.01*jy)
  julday=julday+2-ja+int(0.25*ja)
 end if

 ! compute fraction of the day
 jfrac = fracDay(ih, imin, dsec)

 ! and return the julian day, expressed in fraction of a day
 juldayss = real(julday,kind(rkind)) + jfrac

 end subroutine compjulday

 ! ***************************************************************************************
 ! public subroutine compgregcal: convert julian day (units of days) to calendar date
 ! source: https://en.wikipedia.org/wiki/Julian_day#Julian_or_Gregorian_calendar_from_Julian_day_number
 ! ***************************************************************************************

 subroutine compcalday(julday,                              & !input
                       iyyy,mm,id,ih,imin,dsec,err,message)   !output
 implicit none

 ! input variables
 real(rkind), intent(in)          :: julday       ! julian day

 ! output varibles
 integer(i4b), intent(out)     :: iyyy         ! year
 integer(i4b), intent(out)     :: mm           ! month
 integer(i4b), intent(out)     :: id           ! day
 integer(i4b), intent(out)     :: ih           ! hour
 integer(i4b), intent(out)     :: imin         ! minute
 real(rkind),     intent(out)     :: dsec         ! seconds
 integer(i4b), intent(out)     :: err          ! error code
 character(*), intent(out)     :: message      ! error message

 ! local parameters
 integer(i4b),parameter       :: y = 4716
 integer(i4b),parameter       :: j = 1401
 integer(i4b),parameter       :: m = 2
 integer(i4b),parameter       :: n = 12
 integer(i4b),parameter       :: r = 4
 integer(i4b),parameter       :: p = 1461
 integer(i4b),parameter       :: v = 3
 integer(i4b),parameter       :: u = 5
 integer(i4b),parameter       :: s = 153
 integer(i4b),parameter       :: w = 2
 integer(i4b),parameter       :: b = 274277
 integer(i4b),parameter       :: c = -38
 real(rkind),parameter           :: hr_per_day = 24.0_rkind
 real(rkind),parameter           :: min_per_hour = 60.0_rkind

 ! local variables
 integer(i4b)          :: f,e,g,h                            ! various step variables from wikipedia
 integer(i4b)          :: step_1a,step_1b,step_1c,step_1d    ! temporary variables for calendar calculations
 real(rkind)              :: frac_day  ! fractional day
 real(rkind)              :: remainder ! remainder of modulus operation

 ! initialize errors
 err=0; message="compcalday"
 if(julday<=0)then;err=10;message=trim(message)//"no negative julian days/"; return; end if

 ! step 1
 step_1a = 4*int(julday)+b
 step_1b = step_1a/146097
 step_1c = step_1b*3
 step_1d = step_1c/4

 f = int(julday)+j+step_1d+c

 ! step 2
 e = r * f + v

 ! step 3
 g = mod(e,p)/r

 ! step 4
 h = u * g + w

 ! find day
 id = (mod(h,s))/u + 1

 ! find month
 mm = mod(h/s+m,n)+1

 ! find year
 iyyy = (e/p)-y + (n+m-mm)/n

 ! now find hour,min,second

 frac_day = julday - floor(julday)
 ih = floor((frac_day+1e-9)*hr_per_day)

 remainder = (frac_day+1e-9)*hr_per_day - ih
 imin = floor(remainder*min_per_hour)

 remainder = remainder*min_per_hour - imin
 dsec = nint(remainder*secprmin)

 end subroutine compcalday

 ! ***************************************************************************************
 ! public function elapsedSec: calculate difference of two time marks obtained by date_and_time()
 ! ***************************************************************************************
 function elapsedSec(startTime, endTime)
 integer(i4b),intent(in)        :: startTime(8),endTime(8)            ! state time and end time
 real(rkind)                       :: elapsedSec                         ! elapsed time in seconds
 ! local variables
 integer(i4b)                   :: elapsedDay                         ! elapsed full days
 integer(i4b)                   :: yy                                 ! index of year
 ! number of days of each month
 integer(i4b)                   :: days1(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
 integer(i4b)                   :: days2(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)

 ! calculate the elapsed time smaller than a day
 elapsedSec = (endTime(8)-startTime(8))*.001_rkind + (endTime(7)-startTime(7)) + (endTime(6)-startTime(6))*secprmin + (endTime(5)-startTime(5))*secprhour

 ! check if the run is within the same day otherwise calculate how many days
 if (endTime(1) > startTime(1) .or. endTime(2) > startTime(2) .or. endTime(3) > startTime(3)) then

  elapsedDay = 0
  ! diffenece in year
  do yy = startTime(1), endTime(1) - 1
   elapsedDay = elapsedDay + 365
   if ((mod(yy,4)==0 .and. .not. mod(yy,100)==0) .or. (mod(yy,400)==0)) elapsedDay = elapsedDay + 1
  end do
  if ((mod(startTime(1),4)==0 .and. .not. mod(startTime(1),100)==0) .or. (mod(startTime(1),400)==0)) days1(2) = 29
  if ((mod(endTime(1),4)==0 .and. .not. mod(endTime(1),100)==0) .or. (mod(endTime(1),400)==0)) days2(2) = 29
  ! difference in month
  if (startTime(2) > 1) elapsedDay = elapsedDay - sum(days1(1:(startTime(2)-1)))
  elapsedDay = elapsedDay - startTime(3)
  ! difference in day
  if (endTime(2) > 1) elapsedDay = elapsedDay + sum(days2(1:(endTime(2)-1)))
  elapsedDay = elapsedDay + endTime(3)
  ! convert to seconds
  elapsedSec = elapsedSec + elapsedDay * secprday
 end if
 end function elapsedSec

 ! ***************************************************************************************
 ! public function fracDay: calculate fraction of a day
 ! ***************************************************************************************
 function fracDay(ih, imin, dsec)
 integer(i4b),intent(in)   :: ih,imin      ! hour, minute
 real(rkind),intent(in)       :: dsec         ! seconds
 real(rkind)                  :: fracDay      ! fraction of a day
 ! local variable

 fracDay = (real(ih,kind(rkind))*secprhour + real(imin,kind(rkind))*secprmin + dsec) / secprday
 if(ih < 0) fracDay=-fracDay
 return
 end function fracDay

end module time_utils_module
