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

module read_force_module

! data types
USE nrtype                                    ! variable types, etc.

! derived data types
USE data_types,only:gru_hru_double            ! x%gru(:)%hru(:)%var(:)     (dp)

! constants
USE multiconst,only:secprday                  ! number of seconds in a day

! access missing values
USE globalData,only:realMissing               ! real missing value
USE globalData,only:integerMissing            ! integer missing value

! access the mapping betweeen GRUs and HRUs
USE globalData,only:gru_struc                 ! gru-hru mapping structures

! access the minimum and maximum HRUs in the file
USE globalData,only:ixHRUfile_min,ixHRUfile_max

! global data on the forcing file
USE globalData,only:data_step                 ! length of the data step (s)
USE globalData,only:forcFileInfo              ! forcing file info
USE globalData,only:dJulianStart              ! julian day of start time of simulation
USE globalData,only:refJulday                 ! reference time (fractional julian days)
USE globalData,only:refJulday_data            ! reference time for data files (fractional julian days)
USE globalData,only:fracJulDay                ! fractional julian days since the start of year
USE globalData,only:yearLength                ! number of days in the current year
USE globalData,only:nHRUfile                  ! number of days in the data file

! global metadata
USE globalData,only:time_meta,forc_meta       ! metadata structures
USE var_lookup,only:iLookTIME,iLookFORCE      ! named variables to define structure elements
USE var_lookup,only:iLookDECISIONS            ! named variables for elements of the decision structure

! file paths
USE summaFileManager,only:FORCING_PATH        ! path of the forcing data file

! privacy
implicit none
private
public::read_force

! global parameters
real(rkind),parameter  :: verySmall=1e-3_rkind      ! tiny number
real(rkind),parameter  :: smallOffset=1.e-8_rkind   ! small offset (units=days) to force ih=0 at the start of the day

contains


 ! ************************************************************************************************
 ! public subroutine read_force: read in forcing data
 ! ************************************************************************************************
 subroutine read_force(istep,iFile,iRead,ncid,time_data,forcStruct,err,message)
 ! provide access to subroutines
 USE netcdf                                            ! netcdf capability
 USE time_utils_module,only:compJulday                 ! convert calendar date to julian day
 USE time_utils_module,only:compcalday                 ! convert julian day to calendar date
 USE time_utils_module,only:elapsedSec                 ! calculate the elapsed time
 implicit none
 ! define input variables
 integer(i4b),intent(in)           :: istep            ! time index AFTER the start index
 ! define input-output variables
 integer(i4b),intent(inout)        :: iFile            ! index of current forcing file in forcing file list
 integer(i4b),intent(inout)        :: iRead            ! index of read position in time dimension in current netcdf file
 integer(i4b),intent(inout)        :: ncid             ! netcdf file identifier
 ! define output variables
 integer(i4b),intent(out)          :: time_data(:)     ! vector of time data for a given time step
 type(gru_hru_double)              :: forcStruct       ! x%gru(:)%hru(:)%var(:)     -- model forcing data
 integer(i4b),intent(out)          :: err              ! error code
 character(*),intent(out)          :: message          ! error message
 ! define local variables
 integer(i4b)                      :: nHRUlocal        ! number of HRUs in the local simulation
 integer(i4b)                      :: iGRU,iHRU        ! index of GRU and HRU
 character(len=256),save           :: infile           ! filename
 character(len=256)                :: cmessage         ! error message for downwind routine
 real(rkind)                          :: startJulDay      ! julian day at the start of the year
 real(rkind)                          :: currentJulday    ! Julian day of current time step
 logical(lgt),parameter            :: checkTime=.false.  ! flag to check the time
 ! Start procedure here
 err=0; message="read_force/"

 ! get the number of HRUs in the local simulation
 nHRUlocal = sum(gru_struc(:)%hruCount)

 ! determine the julDay of current model step (istep) we need to read
 if(istep==1)then
  currentJulDay = dJulianStart
 else
  currentJulDay = dJulianStart + (data_step*real(iStep-1,dp))/secprday
 end if

 ! **********************************************************************************************
 ! ***** part 0: if initial step, then open first file and find initial model time step
 ! *****         loop through as many forcing files as necessary to find the initial model step
 ! **********************************************************************************************
 ! check if file is open
 if(ncid==integerMissing)then ! file is closed if ncid==integerMissing

  ! identify the first time step
  call getFirstTimestep(currentJulday,iFile,iRead,ncid,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 end if  ! if the file is not yet open

 ! **********************************************************************************************
 ! ***** part 1: if file open, check to see if we've reached the end of the file, if so close it,
 ! *****         and open new file
 ! *****         Then read the data
 ! **********************************************************************************************
 if(ncid>0)then

  ! check to see if we've passed end of netcdf file
  if(iRead>forcFileInfo(iFile)%nTimeSteps)then

   ! close the NetCDF file
   err = nf90_close(ncid)
   if(err/=nf90_noerr)then; message=trim(message)//'problem closing file ['//trim(infile)//']'; return; endif

   ! increment iFile so we open next forcing file
   iFile = iFile+1
   if(size(forcFileInfo)<iFile)then
    message=trim(message)//'files in list do not include desired data'
    err=20; return
   endif

   ! define new forcing filename
   infile=trim(FORCING_PATH)//trim(forcFileInfo(iFile)%filenmData)

   ! open up the forcing file
   call openForcingFile(iFile,trim(infile),ncId,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

   ! reset iRead since we opened a new file
   iRead=1

  end if  ! if we've passed the end of the NetCDF file

  ! read forcing data
  call readForcingData(currentJulday,ncId,iFile,iRead,nHRUlocal,time_data,forcStruct,err,message)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! check that the file was in fact open
 else
  message=trim(message)//'expect the file to be open'
  err=20; return
 end if  ! end ncid open check

 ! **********************************************************************************************
 ! ***** part 2: compute time
 ! **********************************************************************************************

 ! compute the julian day at the start of the year
 call compjulday(time_data(iLookTIME%iyyy),          & ! input  = year
                 1, 1, 1, 1, 0._rkind,                  & ! input  = month, day, hour, minute, second
                 startJulDay,err,cmessage)                 ! output = julian day (fraction of day) + error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! compute the fractional julian day for the current time step
 call compjulday(time_data(iLookTIME%iyyy),           & ! input  = year
                 time_data(iLookTIME%im),             & ! input  = month
                 time_data(iLookTIME%id),             & ! input  = day
                 time_data(iLookTIME%ih),             & ! input  = hour
                 time_data(iLookTIME%imin),0._rkind,     & ! input  = minute/second
                 currentJulday,err,cmessage)            ! output = julian day (fraction of day) + error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
 ! compute the time since the start of the year (in fractional days)
 fracJulday = currentJulday - startJulDay
 ! set timing of current forcing vector (in seconds since reference day)
 ! NOTE: It is a bit silly to have time information for each HRU and GRU
 do iGRU=1,size(gru_struc)
  do iHRU=1,gru_struc(iGRU)%hruCount
   forcStruct%gru(iGRU)%hru(iHRU)%var(iLookFORCE%time) = (currentJulday-refJulday)*secprday
  end do  ! looping through HRUs
 end do  ! looping through GRUs

 ! compute the number of days in the current year
 yearLength = 365
 if(mod(time_data(iLookTIME%iyyy),4) == 0)then
  yearLength = 366
  if(mod(time_data(iLookTIME%iyyy),100) == 0)then
   yearLength = 365
   if(mod(time_data(iLookTIME%iyyy),400) == 0)then
    yearLength = 366
   end if
  end if
 end if

 ! test
 if(checkTime)then
  write(*,'(i4,1x,4(i2,1x),f9.3,1x,i4)') time_data(iLookTIME%iyyy),           & ! year
                                         time_data(iLookTIME%im),             & ! month
                                         time_data(iLookTIME%id),             & ! day
                                         time_data(iLookTIME%ih),             & ! hour
                                         time_data(iLookTIME%imin),           & ! minute
                                         fracJulday,                          & ! fractional julian day for the current time step
                                         yearLength                             ! number of days in the current year
  !pause ' checking time'
 end if

 end subroutine read_force

 ! *******************************************************************************************************************
 ! *******************************************************************************************************************
 ! *******************************************************************************************************************
 ! *******************************************************************************************************************
 ! *******************************************************************************************************************

 ! *************************************************************************
 ! * private subroutine: find first timestep in any of the forcing files...
 ! *************************************************************************
 subroutine getFirstTimestep(currentJulday,iFile,iRead,ncid,err,message)
 USE netcdf                                            ! netcdf capability
 USE nr_utility_module,only:arth                       ! get a sequence of numbers
 implicit none
 ! define input
 real(rkind),intent(in)               :: currentJulday    ! Julian day of current time step
 ! define input-output variables
 integer(i4b),intent(inout)        :: iFile            ! index of current forcing file in forcing file list
 integer(i4b),intent(inout)        :: iRead            ! index of read position in time dimension in current netcdf file
 integer(i4b),intent(inout)        :: ncid             ! netcdf file identifier
 ! define output variables
 integer(i4b),intent(out)          :: err              ! error code
 character(*),intent(out)          :: message          ! error message
 ! ------------------------------------------------------------------------------------------------------------------
 ! netcdf related
 integer(i4b)                      :: varId            ! variable identifier
 integer(i4b)                      :: dimId            ! dimension identifier
 integer(i4b)                      :: dimLen           ! dimension length
 ! other local variables
 character(len=256),save           :: infile           ! filename
 character(len=256)                :: cmessage         ! error message for downwind routine
 integer(i4b)                      :: nFiles           ! number of forcing files
 real(rkind)                          :: timeVal(1)       ! single time value (restrict time read)
 real(rkind),allocatable              :: fileTime(:)      ! array of time from netcdf file
 real(rkind),allocatable              :: diffTime(:)      ! array of time differences
 ! Start procedure here
 err=0; message="getFirstTimestep/"

 ! get the number of forcing files
 nFiles=size(forcFileInfo)  ! number of forcing files

 ! keep going until we find the file containing the first time step
 do iFile=1,nFiles

  ! define new forcing filename
  infile=trim(FORCING_PATH)//trim(forcFileInfo(iFile)%filenmData)

  ! open netCDF file
  call openForcingFile(iFile,trim(infile),ncId,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

  ! how many time steps in current file?
  err = nf90_inq_dimid(ncid,'time',dimId);             if(err/=nf90_noerr)then; message=trim(message)//'trouble finding time dimension/'//trim(nf90_strerror(err)); return; endif
  err = nf90_inquire_dimension(ncid,dimId,len=dimLen); if(err/=nf90_noerr)then; message=trim(message)//'trouble reading time dimension size/'//trim(nf90_strerror(err)); return; endif

  ! how many HRUs in current file?
  err = nf90_inq_dimid(ncid,'hru',dimId);                if(err/=nf90_noerr)then; message=trim(message)//'trouble finding hru dimension/'//trim(nf90_strerror(err)); return; endif
  err = nf90_inquire_dimension(ncid,dimId,len=nHRUfile); if(err/=nf90_noerr)then; message=trim(message)//'trouble reading hru dimension size/'//trim(nf90_strerror(err)); return; endif

  ! allocate space for time vectors
  if(allocated(fileTime)) deallocate(fileTime)
  if(allocated(diffTime)) deallocate(diffTime)
  allocate(fileTime(dimLen),diffTime(dimLen),stat=err)
  if(err/=0)then; message=trim(message)//'problem allocating time vectors'; return; end if

  ! read time vector from current file
  ! NOTE: This could be faster by checking just the start and the end times

  ! (get variable ID)
  err = nf90_inq_varid(ncid,'time',varId)
  if(err/=nf90_noerr)then; message=trim(message)//'trouble finding time variable/'//trim(nf90_strerror(err)); return; endif

  ! (get single time data value)
  err = nf90_get_var(ncid,varId,timeVal,start=(/1/),count=(/1/))
  if(err/=nf90_noerr)then; message=trim(message)//'trouble reading time vector/'//trim(nf90_strerror(err)); return; endif

  ! get time vector & convert units based on offset and data step
  fileTime = arth(0,1,dimLen) * data_step/secprday + refJulday_data &
             + timeVal(1)/forcFileInfo(iFile)%convTime2Days

  ! find difference of fileTime from currentJulday
  diffTime=abs(fileTime-currentJulday)

  ! start time is in the current file
  if(any(diffTime < verySmall))then

   iRead=minloc(diffTime,1)
   exit

  ! time step is not in current file
  else
   ! close file
   err = nf90_close(ncid)
   if(err/=nf90_noerr)then; message=trim(message)//'trouble closing file '//trim(infile); return; endif

   ! check that it is not the last file
   if(iFile==nFiles)then; err=99; message=trim(message)//'first requested simulation timestep not in any forcing file'; return; end if

  end if  ! first time step is not in any forcing files

 end do ! end of search for model first time step in forcing files

 end subroutine getFirstTimestep

 ! *************************************************************************
 ! * open the NetCDF forcing file and get the time information
 ! *************************************************************************
 subroutine openForcingFile(iFile,infile,ncId,err,message)
 USE netcdf                                              ! netcdf capability
 USE netcdf_util_module,only:nc_file_open                ! open netcdf file
 USE time_utils_module,only:fracDay                      ! compute fractional day
 USE time_utils_module,only:extractTime                  ! extract time info from units string
 USE time_utils_module,only:compJulday                   ! convert calendar date to julian day
 USE globalData,only:tmZoneOffsetFracDay                 ! time zone offset in fractional days
 USE globalData,only:ncTime                              ! time zone information from NetCDF file (timeOffset = longitude/15. - ncTimeOffset)
 USE globalData,only:utcTime                             ! all times in UTC (timeOffset = longitude/15. hours)
 USE globalData,only:localTime                           ! all times local (timeOffset = 0)
 USE summafilemanager,only:NC_TIME_ZONE
 ! dummy variables
 integer(i4b),intent(in)           :: iFile              ! index of current forcing file in forcing file list
 character(*) ,intent(in)          :: infile             ! input file
 integer(i4b) ,intent(out)         :: ncId               ! NetCDF ID
 integer(i4b) ,intent(out)         :: err                ! error code
 character(*) ,intent(out)         :: message            ! error message
 ! local variables
 character(len=256)                :: cmessage           ! error message for downwind routine
 integer(i4b)                      :: iyyy,im,id,ih,imin ! date
 integer(i4b)                      :: ih_tz,imin_tz      ! time zone information
 real(rkind)                          :: dsec,dsec_tz       ! seconds
 integer(i4b)                      :: varId              ! variable identifier
 integer(i4b)                      :: mode               ! netcdf file mode
 integer(i4b)                      :: attLen             ! attribute length
 character(len=256)                :: refTimeString      ! reference time string

 ! initialize error control
 err=0; message='openForcingFile/'

 ! open file
 mode=nf90_NoWrite
 call nc_file_open(trim(infile),mode,ncid,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! get definition of time data
 err = nf90_inq_varid(ncid,'time',varId);                       if(err/=nf90_noerr)then; message=trim(message)//'cannot find time variable/'//trim(nf90_strerror(err)); return; endif
 err = nf90_inquire_attribute(ncid,varId,'units',len = attLen); if(err/=nf90_noerr)then; message=trim(message)//'cannot find time units/'//trim(nf90_strerror(err));    return; endif
 err = nf90_get_att(ncid,varid,'units',refTimeString);          if(err/=nf90_noerr)then; message=trim(message)//'cannot read time units/'//trim(nf90_strerror(err));    return; endif

 ! define the reference time for the model simulation
 call extractTime(refTimeString,                         & ! input  = units string for time data
                  iyyy,im,id,ih,imin,dsec,               & ! output = year, month, day, hour, minute, second
                  ih_tz, imin_tz, dsec_tz,               & ! output = time zone information (hour, minute, second)
                  err,cmessage)                            ! output = error code and error message
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 select case(trim(NC_TIME_ZONE))
  case('ncTime'); tmZoneOffsetFracDay = sign(1, ih_tz) * fracDay(ih_tz,   & ! time zone hour
                                                               imin_tz, & ! time zone minute
                                                               dsec_tz)                        ! time zone second
  case('utcTime');   tmZoneOffsetFracDay = 0._rkind
  case('localTime'); tmZoneOffsetFracDay = 0._rkind
  case default; err=20; message=trim(message)//'unable to identify time zone info option'; return
 end select ! (option time zone option)


 ! convert the reference time to days since the beginning of time
 call compjulday(iyyy,im,id,ih,imin,dsec,                & ! output = year, month, day, hour, minute, second
                 refJulday_data,err,cmessage)              ! output = julian day (fraction of day) + error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! get the time multiplier needed to convert time to units of days
 select case( trim( refTimeString(1:index(refTimeString,' ')) ) )
  case('seconds'); forcFileInfo(iFile)%convTime2Days=86400._rkind
  case('minutes'); forcFileInfo(iFile)%convTime2Days=1440._rkind
  case('hours');   forcFileInfo(iFile)%convTime2Days=24._rkind
  case('days');    forcFileInfo(iFile)%convTime2Days=1._rkind
  case default;    message=trim(message)//'unable to identify time units'; err=20; return
 end select

 end subroutine openForcingFile

 ! *************************************************************************
 ! * read the NetCDF forcing data
 ! *************************************************************************
 subroutine readForcingData(currentJulday,ncId,iFile,iRead,nHRUlocal,time_data,forcStruct,err,message)
 USE netcdf                                            ! netcdf capability
 USE time_utils_module,only:compcalday                 ! convert julian day to calendar date
 USE time_utils_module,only:compJulday                 ! convert calendar date to julian day
 USE get_ixname_module,only:get_ixforce                ! identify index of named variable
 ! dummy variables
 real(rkind),intent(in)               :: currentJulday    ! Julian day of current time step
 integer(i4b) ,intent(in)          :: ncId             ! NetCDF ID
 integer(i4b) ,intent(in)          :: iFile            ! index of forcing file
 integer(i4b) ,intent(in)          :: iRead            ! index in data file
 integer(i4b) ,intent(in)          :: nHRUlocal        ! number of HRUs in the local simulation
 integer(i4b),intent(out)          :: time_data(:)     ! vector of time data for a given time step
 type(gru_hru_double)              :: forcStruct       ! x%gru(:)%hru(:)%var(:)     -- model forcing data
 integer(i4b) ,intent(out)         :: err              ! error code
 character(*) ,intent(out)         :: message          ! error message
 ! local variables
 character(len=256)                :: cmessage         ! error message for downwind routine
 integer(i4b)                      :: varId            ! variable identifier
 character(len = nf90_max_name)    :: varName          ! dimenison name
 real(rkind)                          :: varTime(1)       ! time variable of current forcing data step being read
 ! other local variables
 integer(i4b)                      :: iGRU,iHRU        ! index of GRU and HRU
 integer(i4b)                      :: iHRU_global      ! index of HRU in the NetCDF file
 integer(i4b)                      :: iHRU_local       ! index of HRU in the data subset
 integer(i4b)                      :: iline            ! loop through lines in the file
 integer(i4b)                      :: iNC              ! loop through variables in forcing file
 integer(i4b)                      :: iVar             ! index of forcing variable in forcing data vector
 logical(lgt),parameter            :: checkTime=.false.  ! flag to check the time
 real(rkind)                          :: dsec             ! double precision seconds (not used)
 real(rkind)                          :: dataJulDay       ! julian day of current forcing data step being read
 real(rkind),dimension(nHRUlocal)     :: dataVec          ! vector of data
 real(rkind),dimension(1)             :: dataVal          ! single data value
 real(rkind),parameter                :: dataMin=-1._rkind   ! minimum allowable data value (all forcing variables should be positive)
 logical(lgt),dimension(size(forc_meta)) :: checkForce ! flags to check forcing data variables exist
 logical(lgt),parameter            :: simultaneousRead=.true. ! flag to denote reading all HRUs at once
 ! Start procedure here
 err=0; message="readForcingData/"

 ! initialize time and forcing data structures
 time_data(:) = integerMissing

 ! read time data from iRead location in netcdf file
 err = nf90_inq_varid(ncid,'time',varId);                   if(err/=nf90_noerr)then; message=trim(message)//'trouble finding time variable/'//trim(nf90_strerror(err)); return; endif
 err = nf90_get_var(ncid,varId,varTime,start=(/iRead/));    if(err/=nf90_noerr)then; message=trim(message)//'trouble reading time variable/'//trim(nf90_strerror(err)); return; endif

 ! check that the computed julian day matches the time information in the NetCDF file
 dataJulDay = varTime(1)/forcFileInfo(iFile)%convTime2Days + refJulday_data
 if(abs(currentJulday - dataJulDay) > verySmall)then
  write(message,'(a,f18.8,a,f18.8)') trim(message)//'date for time step: ',dataJulDay,' differs from the expected date: ',currentJulDay
  err=40; return
 end if

 ! convert julian day to time vector
 ! NOTE: use small offset to force ih=0 at the start of the day
 call compcalday(dataJulDay+smallOffset,         & ! input  = julian day
                 time_data(iLookTIME%iyyy),      & ! output = year
                 time_data(iLookTIME%im),        & ! output = month
                 time_data(iLookTIME%id),        & ! output = day
                 time_data(iLookTIME%ih),        & ! output = hour
                 time_data(iLookTIME%imin),dsec, & ! output = minute/second
                 err,cmessage)                     ! output = error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! check to see if any of the time data is missing -- note that it is OK if ih_tz or imin_tz are missing
 if((time_data(iLookTIME%iyyy)==integerMissing) .or. (time_data(iLookTIME%im)==integerMissing) .or. (time_data(iLookTIME%id)==integerMissing) .or. (time_data(iLookTIME%ih)==integerMissing) .or. (time_data(iLookTIME%imin)==integerMissing))then
  do iline=1,size(time_data)
   if(time_data(iline)==integerMissing)then; err=40; message=trim(message)//"variableMissing[var='"//trim(time_meta(iline)%varname)//"']"; return; end if
  end do
 end if

 ! initialize flags for forcing data
 checkForce(:) = .false.
 checkForce(iLookFORCE%time) = .true.  ! time is handled separately

 ! loop through forcing data variables
 do iNC=1,forcFileInfo(iFile)%nVars

  ! check variable is desired
  if(forcFileInfo(iFile)%var_ix(iNC)==integerMissing) cycle

  ! get index in forcing structure
  iVar = forcFileInfo(iFile)%var_ix(iNC)
  checkForce(iVar) = .true.
  
  ! get variable name for error reporting
  err=nf90_inquire_variable(ncid,iNC,name=varName)
  if(err/=nf90_noerr)then; message=trim(message)//'problem reading forcing variable name from netCDF: '//trim(nf90_strerror(err)); return; endif

  ! read forcing data for all HRUs
  if(simultaneousRead)then
   err=nf90_get_var(ncid,forcFileInfo(iFile)%data_id(ivar),dataVec,start=(/ixHRUfile_min,iRead/),count=(/nHRUlocal,1/))
   if(err/=nf90_noerr)then; message=trim(message)//'problem reading forcing data: '//trim(varName)//'/'//trim(nf90_strerror(err)); return; endif
  endif

  ! loop through GRUs and HRUs
  do iGRU=1,size(gru_struc)
   do iHRU=1,gru_struc(iGRU)%hruCount

    ! define global HRU
    iHRU_global = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
    iHRU_local  = (iHRU_global - ixHRUfile_min)+1
    !print*, 'iGRU, iHRU, iHRU_global, iHRU_local = ', iGRU, iHRU, iHRU_global, iHRU_local

    ! read forcing data for a single HRU
    if(.not.simultaneousRead)then
     err=nf90_get_var(ncid,forcFileInfo(iFile)%data_id(ivar),dataVal,start=(/iHRU_global,iRead/))
     if(err/=nf90_noerr)then; message=trim(message)//'problem reading forcing data: '//trim(varName)//'/'//trim(nf90_strerror(err)); return; endif
    endif

    ! check the number of HRUs
    if(iHRU_global > nHRUfile)then
     message=trim(message)//'HRU index exceeds the number of HRUs in the forcing data file'
     err=20; return
    endif

    ! get individual data value
    if(simultaneousRead) dataVal(1) = dataVec(iHRU_local)
    !print*, trim(varname)//': ', dataVal(1)

    ! check individual data value
    if(dataVal(1)<dataMin)then
     write(message,'(a,f13.5)') trim(message)//'forcing data for variable '//trim(varname)//' is less than minimum allowable value ', dataMin
     err=20; return
    endif

    ! put the data into structures
    forcStruct%gru(iGRU)%hru(iHRU)%var(ivar) = dataVal(1)

   end do  ! looping through HRUs within a given GRU
  end do  ! looping through GRUs

 end do  ! loop through forcing variables

 ! check if any forcing data is missing
 if(count(checkForce)<size(forc_meta))then
  do iline=1,size(forc_meta)
   if(.not.checkForce(iline))then
    message=trim(message)//"variableMissing[var='"//trim(forc_meta(iline)%varname)//"']"
    err=20; return
   endif    ! if variable is missing
  end do   ! looping through variables
 end if   ! if any variables are missing

 end subroutine readForcingData


end module read_force_module
