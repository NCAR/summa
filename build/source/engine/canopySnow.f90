module canopySnow_module
! data types
USE nrtype
! physical constants
USE multiconst,only:Tfreeze         ! freezing point of pure water (K)
! model decisions
USE mDecisions_module,only:       &
                      stickySnow, & ! maximum interception capacity an increasing function of temerature
                      lightSnow     ! maximum interception capacity an inverse function of new snow densit
implicit none
private
public::canopySnow

contains

 ! ************************************************************************************************
 ! new subroutine: compute change in snow stored on the vegetation canopy
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
                       mvar_data,                   & ! intent(inout): model variables for a local HRU
                       ! output: error control
                       err,message)                   ! intent(out): error control
 ! ------------------------------------------------------------------------------------------------
 ! provide access to the derived types to define the data structures
 USE data_struc,only:&
                     var_i,            & ! data vector (i4b)
                     var_d,            & ! data vector (dp)
                     var_dlength,      & ! data vector with variable length dimension (dp)
                     model_options       ! defines the model decisions
 ! provide access to named variables defining elements in the data structures
 USE var_lookup,only:iLookTIME,iLookTYPE,iLookATTR,iLookFORCE,iLookPARAM,iLookMVAR,iLookBVAR,iLookINDEX  ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS                               ! named variables for elements of the decision structure
 implicit none
 ! ------------------------------------------------------------------------------------------------
 ! input: model control
 real(dp),intent(in)             :: dt                  ! time step (seconds)
 real(dp),intent(in)             :: exposedVAI          ! exposed vegetation area index -- leaf + stem -- after burial by snow (m2 m-2)
 logical(lgt),intent(in)         :: computeVegFlux      ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input/output: data structures
 type(model_options),intent(in)  :: model_decisions(:)  ! model decisions
 type(var_d),intent(in)          :: forc_data           ! model forcing data
 type(var_d),intent(in)          :: mpar_data           ! model parameters
 type(var_dlength),intent(inout) :: mvar_data           ! model variables for a local HRU
 ! output: error control
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! -------------------------------------------------------------------------------------------------------------------------------
 ! variables in the data structures
 ! input: model decisions
 integer(i4b)                  :: ixSnowInterception         ! choice of option to determine maximum snow interception capacity
 ! input: model forcing data
 real(dp)                      :: scalarAirtemp              ! air temperature (K)
 ! input: model parameters
 real(dp)                      :: refInterceptCapSnow        ! reference canopy interception capacity for snow per unit leaf area (kg m-2)
 real(dp)                      :: ratioDrip2Unloading        ! ratio of canopy drip to snow unloading (-)
 real(dp)                      :: snowUnloadingCoeff         ! time constant for unloading of snow from the forest canopy (s-1)
 ! input: diagnostic variables
 real(dp)                      :: scalarSnowfall             ! computed snowfall rate (kg m-2 s-1)
 real(dp)                      :: scalarNewSnowDensity       ! density of new snow (kg m-3)
 real(dp)                      :: scalarCanopyLiqDrainage    ! liquid drainage from the vegetation canopy (kg m-2 s-1)
 ! input-output: state variables
 real(dp)                      :: scalarCanopyIce            ! mass of ice on the vegetation canopy (kg m-2)
 ! output: diagnostic variables
 real(dp)                      :: scalarThroughfallSnow      ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp)                      :: scalarCanopySnowUnloading  ! unloading of snow from the vegetion canopy (kg m-2 s-1)
 ! ------------------------------------------------------------------------------------------------
 ! local variables
 real(dp),parameter            :: valueMissing=-9999._dp     ! missing value
 integer(i4b)                  :: iter                       ! iteration index
 integer(i4b),parameter        :: maxiter=50                 ! maximum number of iterations
 integer(i4b)                  :: itry                       ! index of loop used for testing
 real(dp)                      :: unloading_melt             ! unloading associated with canopy drip (kg m-2 s-1)
 real(dp)                      :: airtemp_degC               ! value of air temperature in degrees Celcius
 real(dp)                      :: leafScaleFactor            ! scaling factor for interception based on temperature (-)
 real(dp)                      :: leafInterceptCapSnow       ! storage capacity for snow per unit leaf area (kg m-2)
 real(dp)                      :: canopyIceScaleFactor       ! capacity scaling factor for throughfall (kg m-2)
 real(dp)                      :: throughfallDeriv           ! derivative in throughfall flux w.r.t. canopy storage (s-1)
 real(dp)                      :: unloadingDeriv             ! derivative in unloading flux w.r.t. canopy storage (s-1)
 real(dp)                      :: scalarCanopyIceIter        ! trial value for mass of ice on the vegetation canopy (kg m-2) (kg m-2)
 real(dp)                      :: flux                       ! net flux (kg m-2 s-1)
 real(dp)                      :: delS                       ! change in storage (kg m-2)
 real(dp)                      :: resMass                    ! residual in mass equation (kg m-2) 
 real(dp),parameter            :: convTolerMass=0.0001_dp    ! convergence tolerance for mass (kg m-2)
 ! -------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='canopySnow/'
 ! ------------------------------------------------------------------------------------------------
 ! associate variables in the data structure
 associate(&

 ! model decisions
 ixSnowInterception        => model_decisions(iLookDECISIONS%snowIncept)%iDecision,        & ! intent(in): [i4b] choice of option to determine maximum snow interception capacity 

 ! model forcing data
 scalarAirtemp             => forc_data%var(iLookFORCE%airtemp),                           & ! intent(in): [dp] air temperature (K)

 ! model parameters
 refInterceptCapSnow       => mpar_data%var(iLookPARAM%refInterceptCapSnow),               & ! intent(in): [dp] reference canopy interception capacity for snow per unit leaf area (kg m-2)
 ratioDrip2Unloading       => mpar_data%var(iLookPARAM%ratioDrip2Unloading),               & ! intent(in): [dp] ratio of canopy drip to snow unloading (-)
 snowUnloadingCoeff        => mpar_data%var(iLookPARAM%snowUnloadingCoeff),                & ! intent(in): [dp] time constant for unloading of snow from the forest canopy (s-1)

 ! model variables (input)
 scalarSnowfall            => mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1),              & ! intent(in): [dp] computed snowfall rate (kg m-2 s-1)
 scalarNewSnowDensity      => mvar_data%var(iLookMVAR%scalarNewSnowDensity)%dat(1),        & ! intent(in): [dp] density of new snow (kg m-3)
 scalarCanopyLiqDrainage   => mvar_data%var(iLookMVAR%scalarCanopyLiqDrainage)%dat(1),     & ! intent(in): [dp] liquid drainage from the vegetation canopy (kg m-2 s-1)

 ! model variables (input/output)
 scalarCanopyIce           => mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1),             & ! intent(inout): [dp] mass of ice on the vegetation canopy (kg m-2)

 ! model variables (output)
 scalarThroughfallSnow     => mvar_data%var(iLookMVAR%scalarThroughfallSnow)%dat(1),       & ! intent(out): [dp] snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
 scalarCanopySnowUnloading => mvar_data%var(iLookMVAR%scalarCanopySnowUnloading)%dat(1)    & ! intent(out): [dp] unloading of snow from the vegetion canopy (kg m-2 s-1)


 )  ! associate variables in the data structures
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------

 ! compute unloading due to melt drip...
 ! *************************************

 if(computeVegFlux)then
  unloading_melt = min(ratioDrip2Unloading*scalarCanopyLiqDrainage, scalarCanopyIce/dt)  ! kg m-2 s-1
 else
  unloading_melt = 0._dp
 endif
 scalarCanopyIce = scalarCanopyIce - unloading_melt*dt

 ! *****
 ! compute the ice balance due to snowfall and unloading...
 ! ********************************************************

 ! check for early returns
 if(.not.computeVegFlux .or. (scalarSnowfall<tiny(dt) .and. scalarCanopyIce<tiny(dt)))then
  scalarThroughfallSnow     = scalarSnowfall    ! throughfall of snow through the canopy (kg m-2 s-1)
  scalarCanopySnowUnloading = unloading_melt    ! unloading of snow from the canopy (kg m-2 s-1)
  return
 endif

 ! get a trial value for canopy storage
 scalarCanopyIceIter = scalarCanopyIce

 ! iterate
 do iter=1,maxiter

  ! ** compute unloading
  scalarCanopySnowUnloading = snowUnloadingCoeff*scalarCanopyIceIter
  unloadingDeriv            = snowUnloadingCoeff

  ! ** compute throughfall

  ! no snowfall
  if(scalarSnowfall<tiny(dt))then ! no snow
   ! compute throughfall -- note this is effectively zero (no snow case)
   scalarThroughfallSnow = scalarSnowfall  ! throughfall (kg m-2 s-1)
   canopyIceScaleFactor  = valueMissing    ! not used
   throughfallDeriv      = 0._dp

  ! snowfall: compute interception
  else
 
   ! ** process different options for maximum branch snow interception
   select case(ixSnowInterception)

    ! * option 1: maximum interception capacity an inverse function of new snow density (e.g., Mahat and Tarboton, HydProc 2013)
    case(lightSnow)  
     ! (check new snow density is valid)
     if(scalarNewSnowDensity < 0._dp)then; err=20; message=trim(message)//'invalid new snow density'; return; endif
     ! (compute storage capacity of new snow)
     leafScaleFactor       = 0.27_dp + 46._dp/scalarNewSnowDensity
     leafInterceptCapSnow  = refInterceptCapSnow*leafScaleFactor  ! per unit leaf area (kg m-2)

    ! * option 2: maximum interception capacity an increasing function of air temerature
    case(stickySnow)
     airtemp_degC = scalarAirtemp - Tfreeze
     if    (airtemp_degC > -1._dp)then; leafScaleFactor = 4.0_dp
     elseif(airtemp_degC > -3._dp)then; leafScaleFactor = 1.5_dp*airtemp_degC + 5.5_dp
                                  else; leafScaleFactor = 1.0_dp
     endif
     leafInterceptCapSnow = refInterceptCapSnow*leafScaleFactor
     !write(*,'(a,1x,2(f20.10,1x))') 'airtemp_degC, leafInterceptCapSnow = ', airtemp_degC, leafInterceptCapSnow
     !pause 'in stickysnow'
 
    ! check we found the case
    case default
     message=trim(message)//'unable to identify option for maximum branch interception capacity'
     err=20; return

   end select ! identifying option for maximum branch interception capacity

   ! compute maximum interception capacity for the canopy
   canopyIceScaleFactor = leafInterceptCapSnow*exposedVAI

   ! (compute throughfall)
   scalarThroughfallSnow = scalarSnowfall*(scalarCanopyIceIter/canopyIceScaleFactor)
   throughfallDeriv      = scalarSnowfall/canopyIceScaleFactor

   !write(*,'(a,1x,10(e20.10,1x))') 'scalarSnowfall, scalarNewSnowDensity, refInterceptCapSnow, leafScaleFactor, leafInterceptCapSnow, exposedVAI, canopyIceScaleFactor = ', &
   !                                 scalarSnowfall, scalarNewSnowDensity, refInterceptCapSnow, leafScaleFactor, leafInterceptCapSnow, exposedVAI, canopyIceScaleFactor

  endif  ! (if snow is falling)

  !write(*,'(a,1x,10(e20.10,1x))') 'scalarThroughfallSnow, scalarCanopySnowUnloading, unloading_melt = ', &
  !                                 scalarThroughfallSnow, scalarCanopySnowUnloading, unloading_melt

  ! ** compute iteration increment  
  flux = scalarSnowfall - scalarThroughfallSnow - scalarCanopySnowUnloading  ! net flux (kg m-2 s-1)
  delS = (flux*dt - (scalarCanopyIceIter - scalarCanopyIce))/(1._dp + (throughfallDeriv + unloadingDeriv)*dt)
  !write(*,'(a,1x,10(f20.10,1x))') 'scalarCanopyIceIter, flux, delS, scalarSnowfall, scalarThroughfallSnow, scalarCanopySnowUnloading = ',&
  !                                 scalarCanopyIceIter, flux, delS, scalarSnowfall, scalarThroughfallSnow, scalarCanopySnowUnloading 

  ! ** check for convergence
  resMass = scalarCanopyIceIter - (scalarCanopyIce + flux*dt)
  if(abs(resMass) < convTolerMass)exit

  ! ** check for non-convengence
  if(iter==maxiter)then; err=20; message=trim(message)//'failed to converge [mass]'; return; endif

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
