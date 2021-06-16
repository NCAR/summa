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

module derivforce_module

! data types
USE nrtype
USE data_types,only:var_dlength                             ! data structure: x%var(:)%dat (dp)

! model constants
USE multiconst,only:Tfreeze                                 ! freezing point of pure water (K)
USE multiconst,only:secprday                                ! number of seconds in a day
USE multiconst,only:secprhour                               ! number of seconds in an hour
USE multiconst,only:minprhour                               ! number of minutes in an hour

! global time information
USE globalData,only:refJulday                               ! reference time (fractional julian days)
USE globalData,only:data_step                               ! length of the data step (s)
USE globalData,only:tmZoneOffsetFracDay                     ! time zone offset in fractional days

! model decisions
USE globalData,only:model_decisions                         ! model decision structure
USE var_lookup,only:iLookDECISIONS                          ! named variables for elements of the decision structure

! named variables for structure elements
USE var_lookup,only:iLookTIME,iLookATTR                     ! named variables for structure elements
USE var_lookup,only:iLookPARAM,iLookFORCE                   ! named variables for structure elements
USE var_lookup,only:iLookPROG,iLookDIAG,iLookFLUX           ! named variables for structure elements

! look-up values for the choice of the time zone information
USE globalData,only:ncTime,utcTime,localTime                ! time zone info: as in NetCDF file, UTC, or local

! look-up values for the choice of snow albedo options
USE mDecisions_module,only:  &
 constDens,              &    ! Constant new snow density
 anderson,               &    ! Anderson 1976
 hedAndPom,              &    ! Hedstrom and Pomeroy (1998), expoential increase
 pahaut_76                    ! Pahaut 1976, wind speed dependent (derived from Col de Porte, French Alps)

! privacy
implicit none
private
public::derivforce
contains

 ! ************************************************************************************************
 ! public subroutine derivforce: compute derived forcing data
 ! ************************************************************************************************
 subroutine derivforce(time_data,forc_data,attr_data,mpar_data,prog_data,diag_data,flux_data,err,message)
 USE sunGeomtry_module,only:clrsky_rad                       ! compute cosine of the solar zenith angle
 USE conv_funcs_module,only:vapPress                         ! compute vapor pressure of air (Pa)
 USE conv_funcs_module,only:SPHM2RELHM,RELHM2SPHM,WETBULBTMP ! conversion functions
 USE snow_utils_module,only:fracliquid,templiquid            ! functions to compute temperature/liquid water
 USE time_utils_module,only:compcalday                       ! convert julian day to calendar date
 USE summaFileManager,only: NC_TIME_ZONE                     ! time zone option from control file
 ! compute derived forcing data variables
 implicit none
 ! input variables
 integer(i4b),     intent(in)    :: time_data(:)             ! vector of time data for a given time step
 real(rkind),         intent(inout) :: forc_data(:)             ! vector of forcing data for a given time step
 real(rkind),         intent(in)    :: attr_data(:)             ! vector of model attributes
 type(var_dlength),intent(in)    :: mpar_data                ! vector of model parameters
 type(var_dlength),intent(in)    :: prog_data                ! data structure of model prognostic variables for a local HRU
 ! output variables
 type(var_dlength),intent(inout) :: diag_data                ! data structure of model diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data                ! data structure of model fluxes for a local HRU
 integer(i4b),intent(out)        :: err                      ! error code
 character(*),intent(out)        :: message                  ! error message
 ! local time
 integer(i4b)                    :: jyyy,jm,jd               ! year, month, day
 integer(i4b)                    :: jh,jmin                  ! hour, minute
 real(rkind)                        :: dsec                     ! double precision seconds (not used)
 real(rkind)                        :: timeOffset               ! time offset from Grenwich (days)
 real(rkind)                        :: julianTime               ! local julian time
 ! cosine of the solar zenith angle
 real(rkind)                        :: ahour                    ! hour at start of time step
 real(rkind)                        :: dataStep                 ! data step (hours)
 real(rkind)                        :: slope                    ! HRU terrain slope (degrees)
 real(rkind)                        :: azimuth                  ! HRU terrain azimuth (degrees)
 real(rkind)                        :: hri                      ! average radiation index over time step DT
 ! general local variables
 character(len=256)              :: cmessage                 ! error message for downwind routine
 integer(i4b),parameter          :: nBands=2                 ! number of spectral bands
 real(rkind),parameter              :: valueMissing=-9999._rkind   ! missing value
 real(rkind),parameter              :: co2Factor=355.e-6_rkind     ! empirical factor to obtain partial pressure of co2
 real(rkind),parameter              :: o2Factor=0.209_rkind        ! empirical factor to obtain partial pressure of o2
 real(rkind),parameter              :: minMeasHeight=1._rkind      ! minimum measurement height (m)
 real(rkind)                        :: relhum                   ! relative humidity (-)
 real(rkind)                        :: fracrain                 ! fraction of precipitation that falls as rain
 real(rkind)                        :: maxFrozenSnowTemp        ! maximum temperature of snow when the snow is predominantely frozen (K)
 real(rkind),parameter              :: unfrozenLiq=0.01_rkind      ! unfrozen liquid water used to compute maxFrozenSnowTemp (-)
 real(rkind),parameter              :: eps=epsilon(fracrain)    ! a number that is almost negligible
 real(rkind)                        :: Tmin,Tmax                ! minimum and maximum wet bulb temperature in the time step (K)
 real(rkind),parameter              :: pomNewSnowDenMax=150._rkind   ! Upper limit for new snow density limit in Hedstrom and Pomeroy 1998. 150 was used because at was the highest observed density at air temperatures used in this study. See Figure 4 of Hedstrom and Pomeroy (1998).
 real(rkind),parameter              :: andersonWarmDenLimit=2._rkind ! Upper air temperature limit in Anderson (1976) new snow density (C)
 real(rkind),parameter              :: andersonColdDenLimit=15._rkind! Lower air temperature limit in Anderson (1976) new snow density (C)
 real(rkind),parameter              :: andersonDenScal=1.5_rkind     ! Scalar parameter in Anderson (1976) new snow density function (-)
 real(rkind),parameter              :: pahautDenWindScal=0.5_rkind   ! Scalar parameter for wind impacts on density using Pahaut (1976) function (-)
 ! ************************************************************************************************
 ! associate local variables with the information in the data structures
 associate(&
 ! model parameters
 Frad_vis                => mpar_data%var(iLookPARAM%Frad_vis)%dat(1)             , & ! fraction radiation absorbed in visible part of spectrum (-)
 directScale             => mpar_data%var(iLookPARAM%directScale)%dat(1)          , & ! scaling factor for fractional driect radiaion parameterization (-)
 Frad_direct             => mpar_data%var(iLookPARAM%Frad_direct)%dat(1)          , & ! maximum fraction direct radiation (-)
 minwind                 => mpar_data%var(iLookPARAM%minwind)%dat(1)              , & ! minimum windspeed (m s-1)
 fc_param                => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)        , & ! freezing curve parameter for snow (K-1)
 tempCritRain            => mpar_data%var(iLookPARAM%tempCritRain)%dat(1)         , & ! critical temperature where precipitation is rain (K)
 tempRangeTimestep       => mpar_data%var(iLookPARAM%tempRangeTimestep)%dat(1)    , & ! temperature range over the time step (K)
 frozenPrecipMultip      => mpar_data%var(iLookPARAM%frozenPrecipMultip)%dat(1)   , & ! frozen precipitation multiplier (-)
 newSnowDenMin           => mpar_data%var(iLookPARAM%newSnowDenMin)%dat(1)        , & ! minimum new snow density (kg m-3)
 newSnowDenMult          => mpar_data%var(iLookPARAM%newSnowDenMult)%dat(1)       , & ! multiplier for new snow density (kg m-3)
 newSnowDenScal          => mpar_data%var(iLookPARAM%newSnowDenScal)%dat(1)       , & ! scaling factor for new snow density (K)
 constSnowDen            => mpar_data%var(iLookPARAM%constSnowDen)%dat(1)         , & ! Constant new snow density (kg m-3)
 newSnowDenAdd           => mpar_data%var(iLookPARAM%newSnowDenAdd)%dat(1)        , & ! Pahaut 1976, additive factor for new snow density (kg m-3)
 newSnowDenMultTemp      => mpar_data%var(iLookPARAM%newSnowDenMultTemp)%dat(1)   , & ! Pahaut 1976, multiplier for new snow density applied to air temperature (kg m-3 K-1)
 newSnowDenMultWind      => mpar_data%var(iLookPARAM%newSnowDenMultWind)%dat(1)   , & ! Pahaut 1976, multiplier for new snow density applied to wind speed (kg m-7/2 s-1/2)
 newSnowDenMultAnd       => mpar_data%var(iLookPARAM%newSnowDenMultAnd)%dat(1)    , & ! Anderson 1976, multiplier for new snow density for Anderson function (K-1)
 newSnowDenBase          => mpar_data%var(iLookPARAM%newSnowDenBase)%dat(1)       , & ! Anderson 1976, base value that is rasied to the (3/2) power (K)
 heightCanopyTop         => mpar_data%var(iLookPARAM%heightCanopyTop)%dat(1)      , & ! height of the top of the canopy layer (m)
 ! radiation geometry variables
 iyyy                    => time_data(iLookTIME%iyyy)                             , & ! year
 im                      => time_data(iLookTIME%im)                               , & ! month
 id                      => time_data(iLookTIME%id)                               , & ! day
 ih                      => time_data(iLookTIME%ih)                               , & ! hour
 imin                    => time_data(iLookTIME%imin)                             , & ! minute
 latitude                => attr_data(iLookATTR%latitude)                         , & ! latitude (degrees north)
 longitude               => attr_data(iLookATTR%longitude)                        , & ! longitude (degrees east)
 tan_slope               => attr_data(iLookATTR%tan_slope)                        , & ! tan HRU ground surface slope (-)
 aspect                  => attr_data(iLookATTR%aspect)                           , & ! mean azimuth of HRU in degrees E of N (degrees)
 cosZenith               => diag_data%var(iLookDIAG%scalarCosZenith)%dat(1)       , & ! average cosine of the zenith angle over time step DT
 ! measurement height
 mHeight                 => attr_data(iLookATTR%mHeight)                          , & ! latitude (degrees north)
 adjMeasHeight           => diag_data%var(iLookDIAG%scalarAdjMeasHeight)%dat(1)   , & ! adjusted measurement height (m)
 scalarSnowDepth         => prog_data%var(iLookPROG%scalarSnowDepth)%dat(1)       , & ! snow depth on the ground surface (m)
 ! model time
 secondsSinceRefTime     => forc_data(iLookFORCE%time)                            , & ! time = seconds since reference time
 ! model forcing data
 SWRadAtm                => forc_data(iLookFORCE%SWRadAtm)                        , & ! downward shortwave radiation (W m-2)
 airtemp                 => forc_data(iLookFORCE%airtemp)                         , & ! air temperature at 2 meter height (K)
 windspd                 => forc_data(iLookFORCE%windspd)                         , & ! wind speed at 10 meter height (m s-1)
 airpres                 => forc_data(iLookFORCE%airpres)                         , & ! air pressure at 2 meter height (Pa)
 spechum                 => forc_data(iLookFORCE%spechum)                         , & ! specific humidity at 2 meter height (g g-1)
 pptrate                 => forc_data(iLookFORCE%pptrate)                         , & ! precipitation rate (kg m-2 s-1)
 ! derived model forcing data
 scalarO2air             => diag_data%var(iLookDIAG%scalarO2air)%dat(1)           , & ! atmospheric o2 concentration (Pa)
 scalarCO2air            => diag_data%var(iLookDIAG%scalarCO2air)%dat(1)          , & ! atmospheric co2 concentration (Pa)
 ! radiation variables
 scalarFractionDirect    => diag_data%var(iLookDIAG%scalarFractionDirect)%dat(1)  , & ! fraction of direct radiation (0-1)
 spectralIncomingDirect  => flux_data%var(iLookFLUX%spectralIncomingDirect)%dat   , & ! downwelling direct shortwave radiation for each waveband (W m-2)
 spectralIncomingDiffuse => flux_data%var(iLookFLUX%spectralIncomingDiffuse)%dat  , & ! downwelling diffuse shortwave radiation for each waveband (W m-2)
 ! snow accumulation variables
 rainfall                => flux_data%var(iLookFLUX%scalarRainfall)%dat(1)        , & ! computed rainfall rate (kg m-2 s-1)
 snowfall                => flux_data%var(iLookFLUX%scalarSnowfall)%dat(1)        , & ! computed snowfall rate (kg m-2 s-1)
 VPair                   => diag_data%var(iLookDIAG%scalarVPair)%dat(1)           , & ! vapor pressure of the air above the vegetation canopy (Pa)
 twetbulb                => diag_data%var(iLookDIAG%scalarTwetbulb)%dat(1)        , & ! wet bulb temperature (K)
 snowfallTemp            => diag_data%var(iLookDIAG%scalarSnowfallTemp)%dat(1)    , & ! computed temperature of fresh snow (K)
 newSnowDensity          => diag_data%var(iLookDIAG%scalarNewSnowDensity)%dat(1)    & ! computed density of new snow (kg m-3)
 ) ! (associating local variables with the information in the data structures)

 ! initialize error control
 err=0; message="derivforce/"

 ! check spectral dimension
 if(size(spectralIncomingDirect) /= nBands .or. size(spectralIncomingDiffuse) /= nBands)then
  write(message,'(a,i0,a)') trim(message)//'expect ', nBands, 'spectral classes for radiation'
  err=20; return
 end if

 ! adjust the measurement height for the vegetation canopy
 ! NOTE: could return an error or a warning
 ! NOTE: this does not need to be done every time step -- doing here for consistency with the snow adjustment
 if(mHeight < heightCanopyTop)then
  adjMeasHeight = heightCanopyTop+minMeasHeight  ! measurement height at least minMeasHeight above the canopy
 else
  adjMeasHeight = mHeight
 endif

 ! adjust the measurement height for snow depth
 if(adjMeasHeight < scalarSnowDepth+minMeasHeight)then
  adjMeasHeight = scalarSnowDepth+minMeasHeight  ! measurement height at least minMeasHeight above the snow surface
 endif

 ! compute the partial pressure of o2 and co2
 scalarCO2air = co2Factor * airpres  ! atmospheric co2 concentration (Pa)
 scalarO2air  = o2Factor * airpres   ! atmospheric o2 concentration (Pa)

 ! determine timeOffset based on tmZoneInfo option number`
 select case(trim(NC_TIME_ZONE))
  ! Time zone information from NetCDF file
  case('ncTime')
   timeOffset = longitude/360._rkind - tmZoneOffsetFracDay ! time offset in days
  ! All times in UTC
  case('utcTime')
   timeOffset = longitude/360._rkind  ! time offset in days
  ! All times local
  case('localTime')
   timeOffset = 0._rkind  ! time offset in days
  case default; message=trim(message)//'unable to identify option for tmZoneInfo'; err=20; return
 end select ! identifying option tmZoneInfo

 ! constrain timeOffset so that it is in the [-0.5, 0.5] range
 if(timeOffset < -0.5)then
  timeOffset = timeOffset + 1
 else if(timeOffset > 0.5)then
  timeOffset = timeOffset - 1
 endif

 ! compute the local time
 julianTime = secondsSinceRefTime/secprday + refJulday ! julian time (days)

 ! convert julian day to year/month/day/hour/minute
 call compcalday(julianTime+timeOffset,          & ! input  = julian day
                 jyyy,jm,jd,jh,jmin,dsec,        & ! output = year, month, day, hour, minute, second
                 err,cmessage)                     ! output = error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! compute the decimal hour at the start of the time step
 dataStep = data_step/secprhour  ! time step (hours)
 ahour    = real(jh,kind(rkind)) + real(jmin,kind(rkind))/minprhour - data_step/secprhour  ! decimal hour (start of the step)

 ! check slope/aspect intent for radiation calculation
 if(aspect == nr_realMissing)then
  azimuth = 0._rkind              ! if aspect is not an input attribute, slope & azimuth = zero (flat Earth)
  slope   = 0._rkind
 else
  azimuth = aspect                               ! in degrees
  slope   = atan(abs(tan_slope))*180._rkind/PI_D    ! convert from m/m to degrees
 endif

 ! compute the cosine of the solar zenith angle
 call clrsky_rad(jm,jd,ahour,dataStep,   &  ! intent(in): time variables
                 slope,azimuth,latitude, &  ! intent(in): location variables
                 hri,cosZenith)             ! intent(out): cosine of the solar zenith angle
 !write(*,'(a,1x,4(i2,1x),5(f9.3,1x))') 'im,id,ih,imin,ahour,dataStep,azimuth,slope,cosZenith = ', &
 !  im,id,ih,imin,ahour,dataStep,azimuth,slope,cosZenith

 ! ensure solar radiation is non-negative
 if(SWRadAtm < 0._rkind) SWRadAtm = 0._rkind
 ! compute the fraction of direct radiation using the parameterization of Nijssen and Lettenmaier (1999)
 if(cosZenith > 0._rkind)then
  scalarFractionDirect = Frad_direct*cosZenith/(cosZenith + directScale)
 else
  scalarFractionDirect = 0._rkind
 end if
 ! compute direct shortwave radiation, in the visible and near-infra-red part of the spectrum
 spectralIncomingDirect(1) = SWRadAtm*scalarFractionDirect*Frad_vis                         ! (direct vis)
 spectralIncomingDirect(2) = SWRadAtm*scalarFractionDirect*(1._rkind - Frad_vis)               ! (direct nir)
 ! compute diffuse shortwave radiation, in the visible and near-infra-red part of the spectrum
 spectralIncomingDiffuse(1) = SWRadAtm*(1._rkind - scalarFractionDirect)*Frad_vis              ! (diffuse vis)
 spectralIncomingDiffuse(2) = SWRadAtm*(1._rkind - scalarFractionDirect)*(1._rkind - Frad_vis)    ! (diffuse nir)

 !print*,'Frad_direct,scalarFractionDirect,directScale,SWRadAtm,Frad_vis,spectralIncomingDirect: ', &
 !  frad_direct,scalarFractionDirect,directScale,SWRadAtm,Frad_vis,spectralIncomingDirect

 ! ensure wind speed is above a prescribed minimum value
 if(windspd < minwind) windspd=minwind

 ! compute relative humidity (-)
 relhum   = SPHM2RELHM(spechum, airpres, airtemp)
 ! if relative humidity exceeds saturation, then set relative and specific humidity to saturation
 if(relhum > 1._rkind)then
  relhum  = 1._rkind
  spechum = RELHM2SPHM(relhum, airpres, airtemp)
 end if

 ! compute vapor pressure of the air above the vegetation canopy (Pa)
 VPair = vapPress(spechum,airpres)
 !print*, 'VPair = ', VPair

 ! compute wet bulb temperature (K)
 twetbulb = WETBULBTMP(airtemp, relhum, airpres)

 ! compute the maximum temperature of snow when the snow is predominantely frozen (K)
 maxFrozenSnowTemp = templiquid(unfrozenLiq,fc_param)

 ! compute fraction of rain and temperature of fresh snow
 Tmin = twetbulb - tempRangeTimestep/2._rkind
 Tmax = twetbulb + tempRangeTimestep/2._rkind
 if(Tmax < tempCritRain)then
  fracrain     = 0._rkind
  snowfallTemp = twetbulb
 elseif(Tmin > tempCritRain)then
  fracrain     = 1._rkind
  snowfallTemp = maxFrozenSnowTemp
 else
  fracrain     = (Tmax - tempCritRain)/(Tmax - Tmin)
  snowfallTemp = 0.5_rkind*(Tmin + maxFrozenSnowTemp)
 end if

 ! ensure that snowfall temperature creates predominantely solid precipitation
 snowfallTemp      = min(maxFrozenSnowTemp,snowfallTemp) ! snowfall temperature

 ! ensure precipitation rate can be resolved by the data model
 if(pptrate<eps)then
  ! set rainfall and snowfall to zero
  rainfall     = 0._rkind
  snowfall     = 0._rkind
 else
  ! compute rainfall and snowfall
  rainfall = fracrain*pptrate
  snowfall = (1._rkind - fracrain)*pptrate*frozenPrecipMultip
 end if

 ! compute density of new snow
 if(snowfall > tiny(fracrain))then
  ! Determine which method to use
  select case(model_decisions(iLookDECISIONS%snowDenNew)%iDecision)
   ! Hedstrom and Pomeroy 1998
   case(hedAndPom)
    newSnowDensity = min(pomNewSnowDenMax,newSnowDenMin + newSnowDenMult*exp((airtemp-Tfreeze)/newSnowDenScal))  ! new snow density (kg m-3)
   ! Pahaut 1976 (Boone et al. 2002)
   case(pahaut_76)
    newSnowDensity = max(newSnowDenMin,newSnowDenAdd + (newSnowDenMultTemp * (airtemp-Tfreeze))+(newSnowDenMultWind*((windspd)**pahautDenWindScal))); ! new snow density (kg m-3)
   ! Anderson 1976
   case(anderson)
    if(airtemp>(Tfreeze+andersonWarmDenLimit))then
     newSnowDensity = newSnowDenMin + newSnowDenMultAnd*(newSnowDenBase)**(andersonDenScal) ! new snow density (kg m-3)
    elseif(airtemp<=(Tfreeze-andersonColdDenLimit))then
     newSnowDensity = newSnowDenMin ! new snow density (kg m-3)
    else
     newSnowDensity = newSnowDenMin + newSnowDenMultAnd*(airtemp-Tfreeze+newSnowDenBase)**(andersonDenScal) ! new snow density (kg m-3)
    end if
   ! Constant new snow density
   case(constDens)
    newSnowDensity = constSnowDen ! new snow density (kg m-3)
   case default; message=trim(message)//'unable to identify option for new snow density'; err=20; return
  end select ! identifying option for new snow density
 else
  newSnowDensity = valueMissing
  rainfall = rainfall + snowfall ! in most cases snowfall will be zero here
  snowfall = 0._rkind
 end if

 ! end association of local variables with the information in the data structures
 end associate

 end subroutine derivforce


end module derivforce_module
