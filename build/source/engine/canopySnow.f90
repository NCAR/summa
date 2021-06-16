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

module canopySnow_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_i,       &  ! data vector (i4b)
                    var_d,       &  ! data vector (dp)
                    var_dlength, &  ! data vector with variable length dimension (dp)
                    model_options   ! defines the model decisions

! physical constants
USE multiconst,only:Tfreeze         ! freezing point of pure water (K)

! named variables defining elements in the data structures
USE var_lookup,only:iLookFORCE,iLookPARAM,iLookDIAG,iLookPROG,iLookFLUX ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS                                      ! named variables for elements of the decision structure

! model decisions
USE mDecisions_module,only:           &
                      stickySnow,     & ! maximum interception capacity an increasing function of temerature
                      lightSnow,      & ! maximum interception capacity an inverse function of new snow density
                      meltDripUnload, & ! Hedstrom and Pomeroy (1998), Storck et al 2002 (snowUnloadingCoeff & ratioDrip2Unloading)
                      windUnload        ! Roesch et al 2001, formulate unloading based on wind and temperature

! privacy
implicit none
private
public::canopySnow

contains


 ! ************************************************************************************************
 ! public subroutine canopySnow: compute change in snow stored on the vegetation canopy
 ! ************************************************************************************************
 subroutine canopySnow(&
                       ! input: model control
                       dt,                          & ! intent(in): time step (seconds)
                       exposedVAI,                  & ! intent(in): exposed vegetation area index (m2 m-2)
                       computeVegFlux,              & ! intent(in): flag to denote if computing energy flux over vegetation
                       ! input/output: data structures
                       model_decisions,             & ! intent(in):    model decisions
                       forc_data,                   & ! intent(in):    model forcing data
                       mpar_data,                   & ! intent(in):    model parameters
                       diag_data,                   & ! intent(in):    model diagnostic variables for a local HRU
                       prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                       flux_data,                   & ! intent(inout): model flux variables
                       ! output: error control
                       err,message)                   ! intent(out): error control
 ! ------------------------------------------------------------------------------------------------
 implicit none
 ! ------------------------------------------------------------------------------------------------
 ! input: model control
 real(rkind),intent(in)             :: dt                  ! time step (seconds)
 real(rkind),intent(in)             :: exposedVAI          ! exposed vegetation area index -- leaf + stem -- after burial by snow (m2 m-2)
 logical(lgt),intent(in)         :: computeVegFlux      ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input/output: data structures
 type(model_options),intent(in)  :: model_decisions(:)  ! model decisions
 type(var_d),intent(in)          :: forc_data           ! model forcing data
 type(var_dlength),intent(in)    :: mpar_data           ! model parameters
 type(var_dlength),intent(in)    :: diag_data           ! model diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: prog_data           ! model prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data           ! model flux variables
 ! output: error control
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! local variables
 real(rkind),parameter            :: valueMissing=-9999._rkind     ! missing value
 integer(i4b)                  :: iter                       ! iteration index
 integer(i4b),parameter        :: maxiter=50                 ! maximum number of iterations
 real(rkind)                      :: unloading_melt             ! unloading associated with canopy drip (kg m-2 s-1)
 real(rkind)                      :: airtemp_degC               ! value of air temperature in degrees Celcius
 real(rkind)                      :: leafScaleFactor            ! scaling factor for interception based on temperature (-)
 real(rkind)                      :: leafInterceptCapSnow       ! storage capacity for snow per unit leaf area (kg m-2)
 real(rkind)                      :: canopyIceScaleFactor       ! capacity scaling factor for throughfall (kg m-2)
 real(rkind)                      :: throughfallDeriv           ! derivative in throughfall flux w.r.t. canopy storage (s-1)
 real(rkind)                      :: unloadingDeriv             ! derivative in unloading flux w.r.t. canopy storage (s-1)
 real(rkind)                      :: scalarCanopyIceIter        ! trial value for mass of ice on the vegetation canopy (kg m-2) (kg m-2)
 real(rkind)                      :: flux                       ! net flux (kg m-2 s-1)
 real(rkind)                      :: delS                       ! change in storage (kg m-2)
 real(rkind)                      :: resMass                    ! residual in mass equation (kg m-2)
 real(rkind)                      :: tempUnloadingFun           ! temperature unloading functions, Eq. 14 in Roesch et al. 2001
 real(rkind)                      :: windUnloadingFun           ! temperature unloading functions, Eq. 15 in Roesch et al. 2001
 real(rkind),parameter            :: convTolerMass=0.0001_rkind    ! convergence tolerance for mass (kg m-2)
 ! -------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='canopySnow/'
 ! ------------------------------------------------------------------------------------------------
 ! associate variables in the data structure
 associate(&

 ! model decisions
 ixSnowInterception        => model_decisions(iLookDECISIONS%snowIncept)%iDecision,        & ! intent(in): [i4b] choice of option to determine maximum snow interception capacity
 ixSnowUnload              => model_decisions(iLookDECISIONS%snowUnload)%iDecision,        & ! intent(in): [i4b] choice of option to determing how snow unloads from canopy

 ! model forcing data
 scalarAirtemp             => forc_data%var(iLookFORCE%airtemp),                           & ! intent(in): [dp] air temperature (K)

 ! model parameters
 refInterceptCapSnow       => mpar_data%var(iLookPARAM%refInterceptCapSnow)%dat(1),        & ! intent(in): [dp] reference canopy interception capacity for snow per unit leaf area (kg m-2)
 ratioDrip2Unloading       => mpar_data%var(iLookPARAM%ratioDrip2Unloading)%dat(1),        & ! intent(in): [dp] ratio of canopy drip to snow unloading (-)
 snowUnloadingCoeff        => mpar_data%var(iLookPARAM%snowUnloadingCoeff)%dat(1),         & ! intent(in): [dp] time constant for unloading of snow from the forest canopy (s-1)
 minTempUnloading          => mpar_data%var(iLookPARAM%minTempUnloading)%dat(1),           & ! constant describing the minimum temperature for snow unloading in windySnow parameterization (K)
 minWindUnloading          => mpar_data%var(iLookPARAM%minWindUnloading)%dat(1),           & ! constant describing the minimum temperature for snow unloading in windySnow parameterization (K)
 rateTempUnloading         => mpar_data%var(iLookPARAM%rateTempUnloading)%dat(1),          & ! constant describing how quickly snow will unload due to temperature in windySnow parameterization (K s)
 rateWindUnloading         => mpar_data%var(iLookPARAM%rateWindUnloading)%dat(1),          & ! constant describing how quickly snow will unload due to wind in windySnow parameterization (K s)

 ! model diagnostic variables
 scalarNewSnowDensity      => diag_data%var(iLookDIAG%scalarNewSnowDensity)%dat(1),        & ! intent(in): [dp] density of new snow (kg m-3)

 ! model prognostic variables (input/output)
 scalarCanopyIce           => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1),             & ! intent(inout): [dp] mass of ice on the vegetation canopy (kg m-2)

 ! model fluxes (input)
 scalarCanairTemp          => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1),            & ! intent(in): [dp] temperature of the canopy air space (k)
 scalarSnowfall            => flux_data%var(iLookFLUX%scalarSnowfall)%dat(1),              & ! intent(in): [dp] computed snowfall rate (kg m-2 s-1)
 scalarCanopyLiqDrainage   => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1),     & ! intent(in): [dp] liquid drainage from the vegetation canopy (kg m-2 s-1)
 scalarWindspdCanopyTop    => flux_data%var(iLookFLUX%scalarWindspdCanopyTop)%dat(1),      & ! intent(in): [dp] windspeed at the top of the canopy (m s-1)
 ! model variables (output)
 scalarThroughfallSnow     => flux_data%var(iLookFLUX%scalarThroughfallSnow)%dat(1),       & ! intent(out): [dp] snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
 scalarCanopySnowUnloading => flux_data%var(iLookFLUX%scalarCanopySnowUnloading)%dat(1)    & ! intent(out): [dp] unloading of snow from the vegetion canopy (kg m-2 s-1)

 )  ! associate variables in the data structures
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------

 ! compute unloading due to melt drip...
 ! *************************************

 if(computeVegFlux)then
  unloading_melt = min(ratioDrip2Unloading*scalarCanopyLiqDrainage, scalarCanopyIce/dt)  ! kg m-2 s-1
 else
  unloading_melt = 0._rkind
 end if
 scalarCanopyIce = scalarCanopyIce - unloading_melt*dt

 ! *****
 ! compute the ice balance due to snowfall and unloading...
 ! ********************************************************
 ! check for early returns
 if(.not.computeVegFlux .or. (scalarSnowfall<tiny(dt) .and. scalarCanopyIce<tiny(dt)))then
  scalarThroughfallSnow     = scalarSnowfall    ! throughfall of snow through the canopy (kg m-2 s-1)
  scalarCanopySnowUnloading = unloading_melt    ! unloading of snow from the canopy (kg m-2 s-1)
  return
 end if

 ! get a trial value for canopy storage
 scalarCanopyIceIter = scalarCanopyIce
 do iter=1,maxiter
     ! ** compute unloading
     if (ixSnowUnload==meltDripUnload) then
         scalarCanopySnowUnloading = snowUnloadingCoeff*scalarCanopyIceIter
         unloadingDeriv            = snowUnloadingCoeff
     else if (ixSnowUnload==windUnload) then
         tempUnloadingFun = max(scalarCanairTemp - minTempUnloading, 0._rkind) / rateTempUnloading   ! (s-1)
         if (scalarWindspdCanopyTop >= minWindUnloading) then
            windUnloadingFun = abs(scalarWindspdCanopyTop) / rateWindUnloading     ! (s-1)
         else
            windUnloadingFun = 0._rkind ! (s-1)
         end if
         ! implement the "windySnow"  Roesch et al. 2001 parameterization, Eq. 13 in Roesch et al. 2001
         scalarCanopySnowUnloading = scalarCanopyIceIter * (tempUnloadingFun + windUnloadingFun)
         unloadingDeriv            = tempUnloadingFun + windUnloadingFun
     end if
     ! no snowfall
     if(scalarSnowfall<tiny(dt))then ! no snow
         scalarThroughfallSnow = scalarSnowfall  ! throughfall (kg m-2 s-1)
         canopyIceScaleFactor  = valueMissing    ! not used
         throughfallDeriv      = 0._rkind
     else
         ! ** process different options for maximum branch snow interception
         select case(ixSnowInterception)
             case(lightSnow)
                 ! (check new snow density is valid)
                 if(scalarNewSnowDensity < 0._rkind)then; err=20; message=trim(message)//'invalid new snow density'; return; end if
                 ! (compute storage capacity of new snow)
                 leafScaleFactor       = 0.27_rkind + 46._rkind/scalarNewSnowDensity
                 leafInterceptCapSnow  = refInterceptCapSnow*leafScaleFactor  ! per unit leaf area (kg m-2)
             case(stickySnow)
                 airtemp_degC = scalarAirtemp - Tfreeze
                 if (airtemp_degC > -1._rkind) then
                    leafScaleFactor = 4.0_rkind
                 elseif(airtemp_degC > -3._rkind) then
                    leafScaleFactor = 1.5_rkind*airtemp_degC + 5.5_rkind
                 else
                    leafScaleFactor = 1.0_rkind
                 end if
                 leafInterceptCapSnow = refInterceptCapSnow*leafScaleFactor
             case default
                 message=trim(message)//'unable to identify option for maximum branch interception capacity'
                 err=20; return
         end select
         ! compute maximum interception capacity for the canopy
         canopyIceScaleFactor = leafInterceptCapSnow*exposedVAI
         ! (compute throughfall)
         scalarThroughfallSnow = scalarSnowfall*(scalarCanopyIceIter/canopyIceScaleFactor)
         throughfallDeriv      = scalarSnowfall/canopyIceScaleFactor
     end if  ! (if snow is falling)
     ! ** compute iteration increment
     flux = scalarSnowfall - scalarThroughfallSnow - scalarCanopySnowUnloading  ! net flux (kg m-2 s-1)
     delS = (flux*dt - (scalarCanopyIceIter - scalarCanopyIce))/(1._rkind + (throughfallDeriv + unloadingDeriv)*dt)
     ! ** check for convergence
     resMass = scalarCanopyIceIter - (scalarCanopyIce + flux*dt)
     if(abs(resMass) < convTolerMass)exit
     ! ** check for non-convengence
     if(iter==maxiter)then; err=20; message=trim(message)//'failed to converge [mass]'; return; end if
     ! ** update value
     scalarCanopyIceIter = scalarCanopyIceIter + delS
 end do  ! iterating

 ! add the unloading associated with melt drip (kg m-2 s-1)
 scalarCanopySnowUnloading = scalarCanopySnowUnloading + unloading_melt

 ! *****
 ! update mass of ice on the canopy (kg m-2)
 scalarCanopyIce = scalarCanopyIceIter
 ! end association to variables in the data structure
 end associate

 end subroutine canopySnow


end module canopySnow_module
