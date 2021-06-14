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

module coupled_em_module

! numerical recipes data types
USE nrtype

! physical constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    LH_sub,       & ! latent heat of sublimation           (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! data types
USE data_types,only:&
                    var_i,               & ! x%var(:)            (i4b)
                    var_d,               & ! x%var(:)            (dp)
                    var_ilength,         & ! x%var(:)%dat        (i4b)
                    var_dlength            ! x%var(:)%dat        (dp)

! named variables for parent structures
USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
USE var_lookup,only:iLookPROG              ! named variables for structure elements
USE var_lookup,only:iLookDIAG              ! named variables for structure elements
USE var_lookup,only:iLookFLUX              ! named variables for structure elements
USE var_lookup,only:iLookPARAM             ! named variables for structure elements
USE var_lookup,only:iLookINDEX             ! named variables for structure elements
USE globalData,only:iname_snow             ! named variables for snow
USE globalData,only:iname_soil             ! named variables for soil

! named variables for child structures
USE var_lookup,only:childFLUX_MEAN

! metadata
USE globalData,only:indx_meta              ! metadata on the model index variables
USE globalData,only:diag_meta              ! metadata on the model diagnostic variables
USE globalData,only:prog_meta              ! metadata on the model prognostic variables
USE globalData,only:averageFlux_meta       ! metadata on the timestep-average model flux structure

! global data
USE globalData,only:data_step              ! time step of forcing data (s)
USE globalData,only:model_decisions        ! model decision structure
USE globalData,only:globalPrintFlag        ! the global print flag

! look-up values for the numerical method
USE mDecisions_module,only:         &
 iterative,                         &      ! iterative
 nonIterative,                      &      ! non-iterative
 iterSurfEnergyBal                         ! iterate only on the surface energy balance

! look-up values for the maximum interception capacity
USE mDecisions_module,only:         &
                      stickySnow,   &      ! maximum interception capacity an increasing function of temerature
                      lightSnow            ! maximum interception capacity an inverse function of new snow density

! look-up values for the groundwater parameterization
USE mDecisions_module,only:         &
                      qbaseTopmodel,&      ! TOPMODEL-ish baseflow parameterization
                      bigBucket    ,&      ! a big bucket (lumped aquifer model)
                      noExplicit           ! no explicit groundwater parameterization

! look-up values for the spatial representation of groundwater
USE mDecisions_module,only:         &
                      localColumn  ,&      ! separate groundwater representation in each local soil column
                      singleBasin          ! single groundwater store over the entire basin

! privacy
implicit none
private
public::coupled_em
! algorithmic parameters
real(rkind),parameter     :: valueMissing=-9999._rkind  ! missing value, used when diagnostic or state variables are undefined
real(rkind),parameter     :: verySmall=1.e-6_rkind   ! used as an additive constant to check if substantial difference among real numbers
real(rkind),parameter     :: mpe=1.e-6_rkind         ! prevents overflow error if division by zero
real(rkind),parameter     :: dx=1.e-6_rkind          ! finite difference increment
contains


 ! ************************************************************************************************
 ! public subroutine coupled_em: run the coupled energy-mass model for one timestep
 ! ************************************************************************************************
 subroutine coupled_em(&
                       ! model control
                       hruId,             & ! intent(in):    hruId
                       dt_init,           & ! intent(inout): used to initialize the size of the sub-step
                       computeVegFlux,    & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                       ! data structures (input)
                       type_data,         & ! intent(in):    local classification of soil veg etc. for each HRU
                       attr_data,         & ! intent(in):    local attributes for each HRU
                       forc_data,         & ! intent(in):    model forcing data
                       mpar_data,         & ! intent(in):    model parameters
                       bvar_data,         & ! intent(in):    basin-average variables
                       ! data structures (input-output)
                       indx_data,         & ! intent(inout): model indices
                       prog_data,         & ! intent(inout): prognostic variables for a local HRU
                       diag_data,         & ! intent(inout): diagnostic variables for a local HRU
                       flux_data,         & ! intent(inout): model fluxes for a local HRU
                       ! error control
                       err,message)         ! intent(out):   error control
 ! structure allocations
 USE allocspace_module,only:allocLocal      ! allocate local data structures
 USE allocspace_module,only:resizeData      ! clone a data structure
 ! preliminary subroutines
 USE vegPhenlgy_module,only:vegPhenlgy      ! compute vegetation phenology
 USE vegNrgFlux_module,only:wettedFrac      ! compute wetted fraction of the canopy (used in sw radiation fluxes)
 USE snowAlbedo_module,only:snowAlbedo      ! compute snow albedo
 USE vegSWavRad_module,only:vegSWavRad      ! compute canopy sw radiation fluxes
 USE canopySnow_module,only:canopySnow      ! compute interception and unloading of snow from the vegetation canopy
 USE volicePack_module,only:newsnwfall      ! compute change in the top snow layer due to throughfall and unloading
 USE volicePack_module,only:volicePack      ! merge and sub-divide snow layers, if necessary
 USE diagn_evar_module,only:diagn_evar      ! compute diagnostic energy variables -- thermal conductivity and heat capacity
 ! the model solver
 USE indexState_module,only:indexState      ! define indices for all model state variables and layers
 USE opSplittin_module,only:opSplittin      ! solve the system of thermodynamic and hydrology equations for a given substep
 USE time_utils_module,only:elapsedSec      ! calculate the elapsed time
 ! additional subroutines
 USE tempAdjust_module,only:tempAdjust      ! adjust snow temperature associated with new snowfall
 USE snwDensify_module,only:snwDensify      ! snow densification (compaction and cavitation)
 USE var_derive_module,only:calcHeight      ! module to calculate height at layer interfaces and layer mid-point
 ! look-up values for the numerical method
 USE mDecisions_module,only:         &
  iterative,                         &      ! iterative
  nonIterative,                      &      ! non-iterative
  iterSurfEnergyBal                         ! iterate only on the surface energy balance
 ! look-up values for the maximum interception capacity
 USE mDecisions_module,only:          &
                       stickySnow,    &      ! maximum interception capacity an increasing function of temerature
                       lightSnow             ! maximum interception capacity an inverse function of new snow density
 implicit none
 ! model control
 integer(8),intent(in)                :: hruId                  ! hruId
 real(rkind),intent(inout)               :: dt_init                ! used to initialize the size of the sub-step
 logical(lgt),intent(inout)           :: computeVegFlux         ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! data structures (input)
 type(var_i),intent(in)               :: type_data              ! type of vegetation and soil
 type(var_d),intent(in)               :: attr_data              ! spatial attributes
 type(var_d),intent(in)               :: forc_data              ! model forcing data
 type(var_dlength),intent(in)         :: mpar_data              ! model parameters
 type(var_dlength),intent(in)         :: bvar_data              ! basin-average model variables
 ! data structures (input-output)
 type(var_ilength),intent(inout)      :: indx_data              ! state vector geometry
 type(var_dlength),intent(inout)      :: prog_data              ! prognostic variables for a local HRU
 type(var_dlength),intent(inout)      :: diag_data              ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout)      :: flux_data              ! model fluxes for a local HRU
 ! error control
 integer(i4b),intent(out)             :: err                    ! error code
 character(*),intent(out)             :: message                ! error message
 ! =====================================================================================================================================================
 ! =====================================================================================================================================================
 ! local variables
 character(len=256)                   :: cmessage               ! error message
 integer(i4b)                         :: nSnow                  ! number of snow layers
 integer(i4b)                         :: nSoil                  ! number of soil layers
 integer(i4b)                         :: nLayers                ! total number of layers
 integer(i4b)                         :: nState                 ! total number of state variables
 real(rkind)                             :: dtSave                 ! length of last input model sub-step (seconds)
 real(rkind)                             :: dt_sub                 ! length of model sub-step (seconds)
 real(rkind)                             :: dt_wght                ! weight applied to model sub-step (dt_sub/data_step)
 real(rkind)                             :: dt_solv                ! seconds in the data step that have been completed
 real(rkind)                             :: dtMultiplier           ! time step multiplier (-) based on what happenned in "opSplittin"
 real(rkind)                             :: minstep,maxstep        ! minimum and maximum time step length (seconds)
 integer(i4b)                         :: nsub                   ! number of substeps
 logical(lgt)                         :: computeVegFluxOld      ! flag to indicate if we are computing fluxes over vegetation on the previous sub step
 logical(lgt)                         :: includeAquifer         ! flag to denote that an aquifer is included
 logical(lgt)                         :: modifiedLayers         ! flag to denote that snow layers were modified
 logical(lgt)                         :: modifiedVegState       ! flag to denote that vegetation states were modified
 type(var_dlength)                    :: flux_mean              ! timestep-average model fluxes for a local HRU
 integer(i4b)                         :: nLayersRoots           ! number of soil layers that contain roots
 real(rkind)                             :: exposedVAI             ! exposed vegetation area index
 real(rkind)                             :: dCanopyWetFraction_dWat ! derivative in wetted fraction w.r.t. canopy total water (kg-1 m2)
 real(rkind)                             :: dCanopyWetFraction_dT   ! derivative in wetted fraction w.r.t. canopy temperature (K-1)
 real(rkind),parameter                   :: varNotUsed1=-9999._rkind  ! variables used to calculate derivatives (not needed here)
 real(rkind),parameter                   :: varNotUsed2=-9999._rkind  ! variables used to calculate derivatives (not needed here)
 integer(i4b)                         :: iSnow                  ! index of snow layers
 integer(i4b)                         :: iLayer                 ! index of model layers
 real(rkind)                             :: massLiquid             ! mass liquid water (kg m-2)
 real(rkind)                             :: superflousSub          ! superflous sublimation (kg m-2 s-1)
 real(rkind)                             :: superflousNrg          ! superflous energy that cannot be used for sublimation (W m-2 [J m-2 s-1])
 integer(i4b)                         :: ixSolution             ! solution method used by opSplitting
 logical(lgt)                         :: firstSubStep           ! flag to denote if the first time step
 logical(lgt)                         :: stepFailure            ! flag to denote the need to reduce length of the coupled step and try again
 logical(lgt)                         :: tooMuchMelt            ! flag to denote that there was too much melt in a given time step
 logical(lgt)                         :: doLayerMerge           ! flag to denote the need to merge snow layers
 logical(lgt)                         :: pauseFlag              ! flag to pause execution
 logical(lgt),parameter               :: backwardsCompatibility=.true.  ! flag to denote a desire to ensure backwards compatibility with previous branches.
 type(var_ilength)                    :: indx_temp              ! temporary model index variables
 type(var_dlength)                    :: prog_temp              ! temporary model prognostic variables
 type(var_dlength)                    :: diag_temp              ! temporary model diagnostic variables
 ! check SWE
 real(rkind)                             :: oldSWE                 ! SWE at the start of the substep
 real(rkind)                             :: newSWE                 ! SWE at the end of the substep
 real(rkind)                             :: delSWE                 ! change in SWE over the subtep
 real(rkind)                             :: effRainfall            ! effective rainfall (kg m-2 s-1)
 real(rkind)                             :: effSnowfall            ! effective snowfall (kg m-2 s-1)
 real(rkind)                             :: sfcMeltPond            ! surface melt pond (kg m-2)
 real(rkind)                             :: massBalance            ! mass balance error (kg m-2)
 ! balance checks
 integer(i4b)                         :: iVar                   ! loop through model variables
 real(rkind)                             :: totalSoilCompress      ! total soil compression (kg m-2)
 real(rkind)                             :: scalarCanopyWatBalError ! water balance error for the vegetation canopy (kg m-2)
 real(rkind)                             :: scalarSoilWatBalError  ! water balance error (kg m-2)
 real(rkind)                             :: scalarInitCanopyLiq    ! initial liquid water on the vegetation canopy (kg m-2)
 real(rkind)                             :: scalarInitCanopyIce    ! initial ice          on the vegetation canopy (kg m-2)
 real(rkind)                             :: balanceCanopyWater0    ! total water stored in the vegetation canopy at the start of the step (kg m-2)
 real(rkind)                             :: balanceCanopyWater1    ! total water stored in the vegetation canopy at the end of the step (kg m-2)
 real(rkind)                             :: balanceSoilWater0      ! total soil storage at the start of the step (kg m-2)
 real(rkind)                             :: balanceSoilWater1      ! total soil storage at the end of the step (kg m-2)
 real(rkind)                             :: balanceSoilInflux      ! input to the soil zone
 real(rkind)                             :: balanceSoilBaseflow    ! output from the soil zone
 real(rkind)                             :: balanceSoilDrainage    ! output from the soil zone
 real(rkind)                             :: balanceSoilET          ! output from the soil zone
 real(rkind)                             :: balanceAquifer0        ! total aquifer storage at the start of the step (kg m-2)
 real(rkind)                             :: balanceAquifer1        ! total aquifer storage at the end of the step (kg m-2)
 ! test balance checks
 logical(lgt), parameter              :: printBalance=.false.   ! flag to print the balance checks
 real(rkind), allocatable                :: liqSnowInit(:)         ! volumetric liquid water conetnt of snow at the start of the time step
 real(rkind), allocatable                :: liqSoilInit(:)         ! soil moisture at the start of the time step
 ! timing information
 real(rkind)                             :: startTime              ! start time (used to compute wall clock time)
 real(rkind)                             :: endTime                ! end time (used to compute wall clock time)
 ! ----------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="coupled_em/"

 ! This is the start of a data step for a local HRU

 ! get the start time
 call cpu_time(startTime)

 ! check that the decision is supported
 if(model_decisions(iLookDECISIONS%groundwatr)%iDecision==bigBucket .and. &
    model_decisions(iLookDECISIONS%spatial_gw)%iDecision/=localColumn)then
  message=trim(message)//'expect "spatial_gw" decision to equal localColumn when "groundwatr" decision is bigBucket'
  err=20; return
 endif

 ! check if the aquifer is included
 includeAquifer = (model_decisions(iLookDECISIONS%groundwatr)%iDecision==bigBucket)

 ! initialize the numerix tracking variables
 indx_data%var(iLookINDEX%numberFluxCalc       )%dat(1) = 0  ! number of flux calculations                     (-)
 indx_data%var(iLookINDEX%numberStateSplit     )%dat(1) = 0  ! number of state splitting solutions             (-)
 indx_data%var(iLookINDEX%numberDomainSplitNrg )%dat(1) = 0  ! number of domain splitting solutions for energy (-)
 indx_data%var(iLookINDEX%numberDomainSplitMass)%dat(1) = 0  ! number of domain splitting solutions for mass   (-)
 indx_data%var(iLookINDEX%numberScalarSolutions)%dat(1) = 0  ! number of scalar solutions                      (-)

 ! link canopy depth to the information in the data structure
 canopy: associate(canopyDepth => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1) )  ! intent(out): [dp] canopy depth (m)

 ! start by NOT pausing
 pauseFlag=.false.

 ! start by assuming that the step is successful
 stepFailure  = .false.
 doLayerMerge = .false.

 ! initialize flags to modify the veg layers or modify snow layers
 modifiedLayers    = .false.    ! flag to denote that snow layers were modified
 modifiedVegState  = .false.    ! flag to denote that vegetation states were modified

 ! define the first step
 firstSubStep = .true.

 ! count the number of snow and soil layers
 ! NOTE: need to re-compute the number of snow and soil layers at the start of each sub-step because the number of layers may change
 !         (nSnow and nSoil are shared in the data structure)
 nSnow = count(indx_data%var(iLookINDEX%layerType)%dat==iname_snow)
 nSoil = count(indx_data%var(iLookINDEX%layerType)%dat==iname_soil)

 ! compute the total number of snow and soil layers
 nLayers = nSnow + nSoil

 ! create temporary data structures for prognostic variables
 call resizeData(prog_meta(:),prog_data,prog_temp,err=err,message=cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! create temporary data structures for diagnostic variables
 call resizeData(diag_meta(:),diag_data,diag_temp,err=err,message=cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! create temporary data structures for index variables
 call resizeData(indx_meta(:),indx_data,indx_temp,err=err,message=cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for the local fluxes
 call allocLocal(averageFlux_meta(:)%var_info,flux_mean,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! initialize compression and surface melt pond
 sfcMeltPond       = 0._rkind  ! change in storage associated with the surface melt pond (kg m-2)
 totalSoilCompress = 0._rkind  ! change in soil storage associated with compression of the matrix (kg m-2)

 ! initialize mean fluxes
 do iVar=1,size(averageFlux_meta)
  flux_mean%var(iVar)%dat(:) = 0._rkind
 end do

 ! associate local variables with information in the data structures
 associate(&
 ! state variables in the vegetation canopy
 scalarCanopyLiq      => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)                 ,&  ! canopy liquid water (kg m-2)
 scalarCanopyIce      => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)                 ,&  ! canopy ice content (kg m-2)
 ! state variables in the soil domain
 mLayerDepth          => prog_data%var(iLookPROG%mLayerDepth)%dat(nSnow+1:nLayers)       ,&  ! depth of each soil layer (m)
 mLayerVolFracIce     => prog_data%var(iLookPROG%mLayerVolFracIce)%dat(nSnow+1:nLayers)  ,&  ! volumetric ice content in each soil layer (-)
 mLayerVolFracLiq     => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(nSnow+1:nLayers)  ,&  ! volumetric liquid water content in each soil layer (-)
 scalarAquiferStorage => prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)            ,&  ! aquifer storage (m)
 scalarTotalSoilIce   => diag_data%var(iLookDIAG%scalarTotalSoilIce)%dat(1)              ,&  ! total ice in the soil column (kg m-2)
 scalarTotalSoilLiq   => diag_data%var(iLookDIAG%scalarTotalSoilLiq)%dat(1)              ,&  ! total liquid water in the soil column (kg m-2)
 scalarTotalSoilWat   => diag_data%var(iLookDIAG%scalarTotalSoilWat)%dat(1)               &  ! total water in the soil column (kg m-2)
 ) ! (association of local variables with information in the data structures

 ! save the liquid water and ice on the vegetation canopy
 scalarInitCanopyLiq = scalarCanopyLiq    ! initial liquid water on the vegetation canopy (kg m-2)
 scalarInitCanopyIce = scalarCanopyIce    ! initial ice          on the vegetation canopy (kg m-2)

 ! compute total soil moisture and ice at the *START* of the step (kg m-2)
 scalarTotalSoilLiq = sum(iden_water*mLayerVolFracLiq(1:nSoil)*mLayerDepth(1:nSoil))
 scalarTotalSoilIce = sum(iden_water*mLayerVolFracIce(1:nSoil)*mLayerDepth(1:nSoil))  ! NOTE: no expansion and hence use iden_water

 ! compute storage of water in the canopy and the soil
 balanceCanopyWater0 = scalarCanopyLiq + scalarCanopyIce
 balanceSoilWater0   = scalarTotalSoilLiq + scalarTotalSoilIce

 ! get the total aquifer storage at the start of the time step (kg m-2)
 balanceAquifer0 = scalarAquiferStorage*iden_water

 ! save liquid water content
 if(printBalance)then
  allocate(liqSnowInit(nSnow), liqSoilInit(nSoil), stat=err)
  if(err/=0)then
   message=trim(message)//'unable to allocate space for the initial vectors'
   err=20; return
  endif
  if(nSnow>0) liqSnowInit = prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)
  liqSoilInit = mLayerVolFracLiq
 endif

 ! end association of local variables with information in the data structures
 end associate

 ! short-cut to the algorithmic control parameters
 ! NOTE - temporary assignment of minstep to foce something reasonable
 minstep = 10._rkind  ! mpar_data%var(iLookPARAM%minstep)%dat(1)  ! minimum time step (s)
 maxstep = mpar_data%var(iLookPARAM%maxstep)%dat(1)  ! maximum time step (s)
 !print*, 'minstep, maxstep = ', minstep, maxstep

 ! compute the number of layers with roots
 nLayersRoots = count(prog_data%var(iLookPROG%iLayerHeight)%dat(nSnow:nLayers-1) < mpar_data%var(iLookPARAM%rootingDepth)%dat(1)-verySmall)
 if(nLayersRoots == 0)then
  message=trim(message)//'no roots within the soil profile'
  err=20; return
 end if

 ! define the foliage nitrogen factor
 diag_data%var(iLookDIAG%scalarFoliageNitrogenFactor)%dat(1) = 1._rkind  ! foliage nitrogen concentration (1.0 = saturated)

 ! save SWE
 oldSWE = prog_data%var(iLookPROG%scalarSWE)%dat(1)
 !print*, 'nSnow = ', nSnow
 !print*, 'oldSWE = ', oldSWE

 ! *** compute phenology...
 ! ------------------------

 ! compute the temperature of the root zone: used in vegetation phenology
 diag_data%var(iLookDIAG%scalarRootZoneTemp)%dat(1) = sum(prog_data%var(iLookPROG%mLayerTemp)%dat(nSnow+1:nSnow+nLayersRoots)) / real(nLayersRoots, kind(rkind))

 ! remember if we compute the vegetation flux on the previous sub-step
 computeVegFluxOld = computeVegFlux

 ! compute the exposed LAI and SAI and whether veg is buried by snow
 call vegPhenlgy(&
                 ! input/output: data structures
                 model_decisions,             & ! intent(in):    model decisions
                 type_data,                   & ! intent(in):    type of vegetation and soil
                 attr_data,                   & ! intent(in):    spatial attributes
                 mpar_data,                   & ! intent(in):    model parameters
                 prog_data,                   & ! intent(in):    model prognostic variables for a local HRU
                 diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                 ! output
                 computeVegFlux,              & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                 canopyDepth,                 & ! intent(out): canopy depth (m)
                 exposedVAI,                  & ! intent(out): exposed vegetation area index (m2 m-2)
                 err,cmessage)                  ! intent(out): error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! check
 if(computeVegFlux)then
  if(canopyDepth < epsilon(canopyDepth))then
   message=trim(message)//'canopy depth is zero when computeVegFlux flag is .true.'
   err=20; return
  endif
 endif

 ! flag the case where number of vegetation states has changed
 modifiedVegState = (computeVegFlux.neqv.computeVegFluxOld)

 ! *** compute wetted canopy area...
 ! ---------------------------------

 ! compute maximum canopy liquid water (kg m-2)
 diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1) = mpar_data%var(iLookPARAM%refInterceptCapRain)%dat(1)*exposedVAI

 ! compute maximum canopy ice content (kg m-2)
 ! NOTE 1: this is used to compute the snow fraction on the canopy, as used in *BOTH* the radiation AND canopy sublimation routines
 ! NOTE 2: this is a different variable than the max ice used in the throughfall (snow interception) calculations
 ! NOTE 3: use maximum per unit leaf area storage capacity for snow (kg m-2)
 select case(model_decisions(iLookDECISIONS%snowIncept)%iDecision)
  case(lightSnow);  diag_data%var(iLookDIAG%scalarCanopyIceMax)%dat(1) = exposedVAI*mpar_data%var(iLookPARAM%refInterceptCapSnow)%dat(1)
  case(stickySnow); diag_data%var(iLookDIAG%scalarCanopyIceMax)%dat(1) = exposedVAI*mpar_data%var(iLookPARAM%refInterceptCapSnow)%dat(1)*4._rkind
  case default; message=trim(message)//'unable to identify option for maximum branch interception capacity'; err=20; return
 end select ! identifying option for maximum branch interception capacity
 !print*, 'diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1) = ', diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1)
 !print*, 'diag_data%var(iLookDIAG%scalarCanopyIceMax)%dat(1) = ', diag_data%var(iLookDIAG%scalarCanopyIceMax)%dat(1)

 ! compute wetted fraction of the canopy
 ! NOTE: assume that the wetted fraction is constant over the substep for the radiation calculations
 if(computeVegFlux)then

  ! compute wetted fraction of the canopy
  call wettedFrac(&
                  ! input
                  .false.,                                                      & ! flag to denote if derivatives are required
                  .false.,                                                      & ! flag to denote if derivatives are calculated numerically
                  (prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1) < Tfreeze), & ! flag to denote if the canopy is frozen
                  varNotUsed1,                                                  & ! derivative in canopy liquid w.r.t. canopy temperature (kg m-2 K-1)
                  varNotUsed2,                                                  & ! fraction of liquid water on the canopy
                  prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1),              & ! canopy liquid water (kg m-2)
                  prog_data%var(iLookPROG%scalarCanopyIce)%dat(1),              & ! canopy ice (kg m-2)
                  diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1),           & ! maximum canopy liquid water (kg m-2)
                  diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1),           & ! maximum canopy ice content (kg m-2)
                  mpar_data%var(iLookPARAM%canopyWettingFactor)%dat(1),         & ! maximum wetted fraction of the canopy (-)
                  mpar_data%var(iLookPARAM%canopyWettingExp)%dat(1),            & ! exponent in canopy wetting function (-)
                  ! output
                  diag_data%var(iLookDIAG%scalarCanopyWetFraction)%dat(1),      & ! canopy wetted fraction (-)
                  dCanopyWetFraction_dWat,                                      & ! derivative in wetted fraction w.r.t. canopy liquid water content (kg-1 m2)
                  dCanopyWetFraction_dT,                                        & ! derivative in wetted fraction w.r.t. canopy liquid water content (kg-1 m2)
                  err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! vegetation is completely buried by snow (or no veg exists at all)
 else
  diag_data%var(iLookDIAG%scalarCanopyWetFraction)%dat(1) = 0._rkind
  dCanopyWetFraction_dWat                                 = 0._rkind
  dCanopyWetFraction_dT                                   = 0._rkind
 end if

 ! *** compute snow albedo...
 ! --------------------------
 ! NOTE: this should be done before the radiation calculations
 ! NOTE: uses snowfall; should really use canopy throughfall + canopy unloading
 call snowAlbedo(&
                 ! input: model control
                 data_step,                   & ! intent(in): model time step (s)
                 (nSnow > 0),                 & ! intent(in): logical flag to denote if snow is present
                 ! input/output: data structures
                 model_decisions,             & ! intent(in):    model decisions
                 mpar_data,                   & ! intent(in):    model parameters
                 flux_data,                   & ! intent(in):    model flux variables
                 diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                 prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                 ! output: error control
                 err,cmessage)                  ! intent(out): error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if


 ! *** compute canopy sw radiation fluxes...
 ! -----------------------------------------
 call vegSWavRad(&
                 data_step,                    & ! intent(in):    time step (s) -- only used in Noah-MP radiation, to compute albedo
                 nSnow,                        & ! intent(in):    number of snow layers
                 nSoil,                        & ! intent(in):    number of soil layers
                 nLayers,                      & ! intent(in):    total number of layers
                 computeVegFlux,               & ! intent(in):    logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                 type_data,                    & ! intent(in):    type of vegetation and soil
                 prog_data,                    & ! intent(inout): model prognostic variables for a local HRU
                 diag_data,                    & ! intent(inout): model diagnostic variables for a local HRU
                 flux_data,                    & ! intent(inout): model flux variables
                 err,cmessage)                   ! intent(out):   error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if


 ! *** compute canopy throughfall and unloading...
 ! -----------------------------------------------
 ! NOTE 1: this needs to be done before solving the energy and liquid water equations, to account for the heat advected with precipitation (and throughfall/unloading)
 ! NOTE 2: the unloading flux is computed using canopy drip (scalarCanopyLiqDrainage) from the previous time step
 call canopySnow(&
                 ! input: model control
                 data_step,                   & ! intent(in): time step (seconds)
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
                 err,cmessage)                  ! intent(out): error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! adjust canopy temperature to account for new snow
 if(computeVegFlux)then ! logical flag to compute vegetation fluxes (.false. if veg buried by snow)
  call tempAdjust(&
                  ! input: derived parameters
                  canopyDepth,                 & ! intent(in): canopy depth (m)
                  ! input/output: data structures
                  mpar_data,                   & ! intent(in):    model parameters
                  prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                  diag_data,                   & ! intent(out):   model diagnostic variables for a local HRU
                  ! output: error control
                  err,cmessage)                  ! intent(out): error control
                  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
  endif ! if computing fluxes over vegetation

 ! initialize drainage and throughfall
 ! NOTE 1: this needs to be done before solving the energy and liquid water equations, to account for the heat advected with precipitation
 ! NOTE 2: this initialization needs to be done AFTER the call to canopySnow, since canopySnow uses canopy drip drom the previous time step
 if(.not.computeVegFlux)then
  flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1)   = flux_data%var(iLookFLUX%scalarRainfall)%dat(1)
  flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1) = 0._rkind
 else
  flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1)   = 0._rkind
  flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1) = 0._rkind
 end if

 ! ****************************************************************************************************
 ! *** MAIN SOLVER ************************************************************************************
 ! ****************************************************************************************************

 ! initialize the length of the sub-step
 dt_solv = 0._rkind   ! length of time step that has been completed (s)
 dt_init = min(data_step,maxstep)  ! initial substep length (s)
 dt_sub  = dt_init                 ! length of substep
 dtSave  = dt_init                 ! length of substep

 ! initialize the number of sub-steps
 nsub=0

 ! loop through sub-steps
 substeps: do  ! continuous do statement with exit clause (alternative to "while")

  ! print progress
  !print*, '*** new substep'
  !write(*,'(a,3(f11.4,1x))') 'dt_sub, dt_init = ', dt_sub, dt_init

  ! print progress
  if(globalPrintFlag)then
   write(*,'(a,1x,4(f13.5,1x))') ' start of step: dt_init, dt_sub, dt_solv, data_step: ', dt_init, dt_sub, dt_solv, data_step
   print*, 'stepFailure = ', stepFailure
   print*, 'before resizeData: nSnow, nSoil = ', nSnow, nSoil
  endif

  ! increment the number of sub-steps
  nsub = nsub+1

  ! resize the "indx_data" structure
  ! NOTE: this is necessary because the length of index variables depends on a given split
  !        --> the resize here is overwritten later (in indexSplit)
  !        --> admittedly ugly, and retained for now
  if(stepFailure)then
   call resizeData(indx_meta(:),indx_temp,indx_data,err=err,message=cmessage)
   if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
  else
   call resizeData(indx_meta(:),indx_data,indx_temp,err=err,message=cmessage)
   if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
  endif

  ! save/recover copies of index variables
  do iVar=1,size(indx_data%var)
   !print*, 'indx_meta(iVar)%varname = ', trim(indx_meta(iVar)%varname)
   select case(stepFailure)
    case(.false.); indx_temp%var(iVar)%dat(:) = indx_data%var(iVar)%dat(:)
    case(.true.);  indx_data%var(iVar)%dat(:) = indx_temp%var(iVar)%dat(:)
   end select
  end do  ! looping through variables

  ! save/recover copies of prognostic variables
  do iVar=1,size(prog_data%var)
   !print*, 'prog_meta(iVar)%varname = ', trim(prog_meta(iVar)%varname)
   select case(stepFailure)
    case(.false.); prog_temp%var(iVar)%dat(:) = prog_data%var(iVar)%dat(:)
    case(.true.);  prog_data%var(iVar)%dat(:) = prog_temp%var(iVar)%dat(:)
   end select
  end do  ! looping through variables

  ! save/recover copies of diagnostic variables
  do iVar=1,size(diag_data%var)
   select case(stepFailure)
    case(.false.); diag_temp%var(iVar)%dat(:) = diag_data%var(iVar)%dat(:)
    case(.true.);  diag_data%var(iVar)%dat(:) = diag_temp%var(iVar)%dat(:)
   end select
  end do  ! looping through variables

  ! re-assign dimension lengths
  nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_snow)
  nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_soil)
  nLayers = nSnow+nSoil

  ! *** merge/sub-divide snow layers...
  ! -----------------------------------
  call volicePack(&
                  ! input/output: model data structures
                  doLayerMerge,                & ! intent(in):    flag to force merge of snow layers
                  model_decisions,             & ! intent(in):    model decisions
                  mpar_data,                   & ! intent(in):    model parameters
                  indx_data,                   & ! intent(inout): type of each layer
                  prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                  diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                  flux_data,                   & ! intent(inout): model fluxes for a local HRU
                  ! output
                  modifiedLayers,              & ! intent(out): flag to denote that layers were modified
                  err,cmessage)                  ! intent(out): error control
  if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if

  ! save the number of snow and soil layers
  nSnow   = indx_data%var(iLookINDEX%nSnow)%dat(1)
  nSoil   = indx_data%var(iLookINDEX%nSoil)%dat(1)
  nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)


  ! compute the indices for the model state variables
  if(firstSubStep .or. modifiedVegState .or. modifiedLayers)then
   call indexState(computeVegFlux,          & ! intent(in):    flag to denote if computing the vegetation flux
                   includeAquifer,          & ! intent(in):    flag to denote if included the aquifer
                   nSnow,nSoil,nLayers,     & ! intent(in):    number of snow and soil layers, and total number of layers
                   indx_data,               & ! intent(inout): indices defining model states and layers
                   err,cmessage)              ! intent(out):   error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  end if

  ! recreate the temporary data structures
  ! NOTE: resizeData(meta, old, new, ..)
  if(modifiedVegState .or. modifiedLayers)then

   ! create temporary data structures for prognostic variables
   call resizeData(prog_meta(:),prog_data,prog_temp,copy=.true.,err=err,message=cmessage)
   if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

   ! create temporary data structures for diagnostic variables
   call resizeData(diag_meta(:),diag_data,diag_temp,copy=.true.,err=err,message=cmessage)
   if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

   ! create temporary data structures for index variables
   call resizeData(indx_meta(:),indx_data,indx_temp,copy=.true.,err=err,message=cmessage)
   if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

   do iVar=1,size(indx_data%var)
    !print*, 'indx_meta(iVar)%varname = ', trim(indx_meta(iVar)%varname)
    select case(stepFailure)
     case(.false.); indx_temp%var(iVar)%dat(:) = indx_data%var(iVar)%dat(:)
     case(.true.);  indx_data%var(iVar)%dat(:) = indx_temp%var(iVar)%dat(:)
    end select
   end do  ! looping through variables

  endif  ! if modified the states

  ! define the number of state variables
  nState = indx_data%var(iLookINDEX%nState)%dat(1)

  ! *** compute diagnostic variables for each layer...
  ! --------------------------------------------------
  ! NOTE: this needs to be done AFTER volicePack, since layers may have been sub-divided and/or merged
  call diagn_evar(&
                  ! input: control variables
                  computeVegFlux,          & ! intent(in): flag to denote if computing the vegetation flux
                  canopyDepth,             & ! intent(in): canopy depth (m)
                  ! input/output: data structures
                  mpar_data,               & ! intent(in):    model parameters
                  indx_data,               & ! intent(in):    model layer indices
                  prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                  diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                  ! output: error control
                  err,cmessage)              ! intent(out): error control
  if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if


  ! *** compute melt of the "snow without a layer"...
  ! -------------------------------------------------
  ! NOTE: forms a surface melt pond, which drains into the upper-most soil layer through the time step
  ! (check for the special case of "snow without a layer")
  if(nSnow==0)then
   call implctMelt(&
                   ! input/output: integrated snowpack properties
                   prog_data%var(iLookPROG%scalarSWE)%dat(1),               & ! intent(inout): snow water equivalent (kg m-2)
                   prog_data%var(iLookPROG%scalarSnowDepth)%dat(1),         & ! intent(inout): snow depth (m)
                   prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1),       & ! intent(inout): surface melt pond (kg m-2)
                   ! input/output: properties of the upper-most soil layer
                   prog_data%var(iLookPROG%mLayerTemp)%dat(nSnow+1),        & ! intent(inout): surface layer temperature (K)
                   prog_data%var(iLookPROG%mLayerDepth)%dat(nSnow+1),       & ! intent(inout): surface layer depth (m)
                   diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat(nSnow+1),& ! intent(inout): surface layer volumetric heat capacity (J m-3 K-1)
                   ! output: error control
                   err,cmessage                                             ) ! intent(out): error control
   if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
  end if

  ! *** solve model equations...
  ! ----------------------------

  ! save input step
  dtSave = dt_sub
  !write(*,'(a,1x,3(f12.5,1x))') trim(message)//'before opSplittin: dt_init, dt_sub, dt_solv = ', dt_init, dt_sub, dt_solv

  ! get the new solution
  call opSplittin(&
                  ! input: model control
                  nSnow,                                  & ! intent(in):    number of snow layers
                  nSoil,                                  & ! intent(in):    number of soil layers
                  nLayers,                                & ! intent(in):    total number of layers
                  nState,                                 & ! intent(in):    total number of layers
                  dt_sub,                                 & ! intent(in):    length of the model sub-step
                  (nsub==1),                              & ! intent(in):    logical flag to denote the first substep
                  computeVegFlux,                         & ! intent(in):    logical flag to compute fluxes within the vegetation canopy
                  ! input/output: data structures
                  type_data,                              & ! intent(in):    type of vegetation and soil
                  attr_data,                              & ! intent(in):    spatial attributes
                  forc_data,                              & ! intent(in):    model forcing data
                  mpar_data,                              & ! intent(in):    model parameters
                  indx_data,                              & ! intent(inout): index data
                  prog_data,                              & ! intent(inout): model prognostic variables for a local HRU
                  diag_data,                              & ! intent(inout): model diagnostic variables for a local HRU
                  flux_data,                              & ! intent(inout): model fluxes for a local HRU
                  bvar_data,                              & ! intent(in):    model variables for the local basin
                  model_decisions,                        & ! intent(in):    model decisions
                  ! output: model control
                  dtMultiplier,                           & ! intent(out):   substep multiplier (-)
                  tooMuchMelt,                            & ! intent(out):   flag to denote that ice is insufficient to support melt
                  stepFailure,                            & ! intent(out):   flag to denote that the coupled step failed
                  ixSolution,                             & ! intent(out):   solution method used in this iteration
                  err,cmessage)                             ! intent(out):   error code and error message

  ! check for all errors (error recovery within opSplittin)
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
  !print*, 'completed step'
  !print*, 'PAUSE: '; read(*,*)

  ! process the flag for too much melt
  if(tooMuchMelt)then
   stepFailure  = .true.
   doLayerMerge = .true.
  else
   doLayerMerge = .false.
  endif

  ! handle special case of the step failure
  ! NOTE: need to revert back to the previous state vector that we were happy with and reduce the time step
  if(stepFailure)then

   ! halve step
   dt_sub = dtSave/2._rkind

   ! check that the step is not tiny
   if(dt_sub < minstep)then
    print*,ixSolution
    print*, 'dtSave, dt_sub', dtSave, dt_sub
    message=trim(message)//'length of the coupled step is below the minimum step length'
    err=20; return
   endif

   ! try again
   cycle substeps

  endif

  ! update first step
  firstSubStep=.false.

  ! ***  remove ice due to sublimation...
  ! --------------------------------------------------------------
  sublime: associate(&
   scalarCanopySublimation => flux_data%var(iLookFLUX%scalarCanopySublimation)%dat(1), & ! sublimation from the vegetation canopy (kg m-2 s-1)
   scalarSnowSublimation   => flux_data%var(iLookFLUX%scalarSnowSublimation)%dat(1),   & ! sublimation from the snow surface (kg m-2 s-1)
   scalarLatHeatCanopyEvap => flux_data%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1), & ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
   scalarSenHeatCanopy     => flux_data%var(iLookFLUX%scalarSenHeatCanopy)%dat(1),     & ! sensible heat flux from the canopy to the canopy air space (W m-2)
   scalarLatHeatGround     => flux_data%var(iLookFLUX%scalarLatHeatGround)%dat(1),     & ! latent heat flux from ground surface below vegetation (W m-2)
   scalarSenHeatGround     => flux_data%var(iLookFLUX%scalarSenHeatGround)%dat(1),     & ! sensible heat flux from ground surface below vegetation (W m-2)
   scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1),         & ! liquid water stored on the vegetation canopy (kg m-2)
   scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1),         & ! ice          stored on the vegetation canopy (kg m-2)
   mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat,           & ! volumetric fraction of ice in the snow+soil domain (-)
   mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat,           & ! volumetric fraction of liquid water in the snow+soil domain (-)
   mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat                 & ! depth of each snow+soil layer (m)
  ) ! associations to variables in data structures

  ! * compute change in canopy ice content due to sublimation...
  ! ------------------------------------------------------------
  if(computeVegFlux)then

   ! remove mass of ice on the canopy
   scalarCanopyIce = scalarCanopyIce + scalarCanopySublimation*dt_sub

   ! if removed all ice, take the remaining sublimation from water
   if(scalarCanopyIce < 0._rkind)then
    scalarCanopyLiq = scalarCanopyLiq + scalarCanopyIce
    scalarCanopyIce = 0._rkind
   endif

   ! modify fluxes if there is insufficient canopy water to support the converged sublimation rate over the time step dt_sub
   if(scalarCanopyLiq < 0._rkind)then
    ! --> superfluous sublimation flux
    superflousSub = -scalarCanopyLiq/dt_sub  ! kg m-2 s-1
    superflousNrg = superflousSub*LH_sub     ! W m-2 (J m-2 s-1)
    ! --> update fluxes and states
    scalarCanopySublimation = scalarCanopySublimation + superflousSub
    scalarLatHeatCanopyEvap = scalarLatHeatCanopyEvap + superflousNrg
    scalarSenHeatCanopy     = scalarSenHeatCanopy - superflousNrg
    scalarCanopyLiq         = 0._rkind
   endif

  end if  ! (if computing the vegetation flux)

  ! * compute change in ice content of the top snow layer due to sublimation...
  ! ---------------------------------------------------------------------------
  ! NOTE: this is done BEFORE densification
  if(nSnow > 0)then ! snow layers exist

   ! try to remove ice from the top layer
   iSnow=1

   ! save the mass of liquid water (kg m-2)
   massLiquid = mLayerDepth(iSnow)*mLayerVolFracLiq(iSnow)*iden_water

   ! add/remove the depth of snow gained/lost by frost/sublimation (m)
   ! NOTE: assume constant density
   mLayerDepth(iSnow) = mLayerDepth(iSnow) + dt_sub*scalarSnowSublimation/(mLayerVolFracIce(iSnow)*iden_ice)

   ! check that we did not remove the entire layer
   if(mLayerDepth(iSnow) < verySmall)then
    stepFailure  = .true.
    doLayerMerge = .true.
    dt_sub      = max(dtSave/2._rkind, minstep)
    cycle substeps
   else
    stepFailure  = .false.
    doLayerMerge = .false.
   endif

   ! update the volumetric fraction of liquid water
   mLayerVolFracLiq(iSnow) = massLiquid / (mLayerDepth(iSnow)*iden_water)

  ! no snow
  else

   ! no snow: check that sublimation is zero
   if(abs(scalarSnowSublimation) > verySmall)then
    message=trim(message)//'sublimation of snow has been computed when no snow exists'
    err=20; return
   end if

  end if  ! (if snow layers exist)

  end associate sublime

  ! *** account for compaction and cavitation in the snowpack...
  ! ------------------------------------------------------------
  if(nSnow>0)then
   call snwDensify(&
                   ! intent(in): variables
                   dt_sub,                                                 & ! intent(in): time step (s)
                   indx_data%var(iLookINDEX%nSnow)%dat(1),                 & ! intent(in): number of snow layers
                   prog_data%var(iLookPROG%mLayerTemp)%dat(1:nSnow),       & ! intent(in): temperature of each layer (K)
                   diag_data%var(iLookDIAG%mLayerMeltFreeze)%dat(1:nSnow), & ! intent(in): volumetric melt in each layer (kg m-3)
                   ! intent(in): parameters
                   mpar_data%var(iLookPARAM%densScalGrowth)%dat(1),        & ! intent(in): density scaling factor for grain growth (kg-1 m3)
                   mpar_data%var(iLookPARAM%tempScalGrowth)%dat(1),        & ! intent(in): temperature scaling factor for grain growth (K-1)
                   mpar_data%var(iLookPARAM%grainGrowthRate)%dat(1),       & ! intent(in): rate of grain growth (s-1)
                   mpar_data%var(iLookPARAM%densScalOvrbdn)%dat(1),        & ! intent(in): density scaling factor for overburden pressure (kg-1 m3)
                   mpar_data%var(iLookPARAM%tempScalOvrbdn)%dat(1),        & ! intent(in): temperature scaling factor for overburden pressure (K-1)
                   mpar_data%var(iLookPARAM%baseViscosity)%dat(1),         & ! intent(in): viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)
                   ! intent(inout): state variables
                   prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow),      & ! intent(inout): depth of each layer (m)
                   prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow), & ! intent(inout):  volumetric fraction of liquid water after itertations (-)
                   prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow), & ! intent(inout):  volumetric fraction of ice after itertations (-)
                   ! output: error control
                   err,cmessage)                     ! intent(out): error control
   if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if
  end if  ! if snow layers exist

  ! update coordinate variables
  call calcHeight(&
                  ! input/output: data structures
                  indx_data,   & ! intent(in): layer type
                  prog_data,   & ! intent(inout): model variables for a local HRU
                  ! output: error control
                  err,cmessage)
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

  ! recompute snow depth and SWE
  if(nSnow > 0)then
   prog_data%var(iLookPROG%scalarSnowDepth)%dat(1) = sum(  prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow))
   prog_data%var(iLookPROG%scalarSWE)%dat(1)       = sum( (prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                           prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                         * prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow) )
  end if

  ! increment fluxes
  dt_wght = dt_sub/data_step ! define weight applied to each sub-step
  do iVar=1,size(averageFlux_meta)
   flux_mean%var(iVar)%dat(:) = flux_mean%var(iVar)%dat(:) + flux_data%var(averageFlux_meta(iVar)%ixParent)%dat(:)*dt_wght
  end do

  ! increment change in storage associated with the surface melt pond (kg m-2)
  if(nSnow==0) sfcMeltPond = sfcMeltPond + prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1)

  ! increment soil compression (kg m-2)
  totalSoilCompress = totalSoilCompress + diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1) ! total soil compression over whole layer (kg m-2)

  ! ****************************************************************************************************
  ! *** END MAIN SOLVER ********************************************************************************
  ! ****************************************************************************************************

  ! increment sub-step
  dt_solv = dt_solv + dt_sub

  ! save the time step to initialize the subsequent step
  if(dt_solv<data_step .or. nsub==1) dt_init = dt_sub

  ! check
  if(globalPrintFlag)&
  write(*,'(a,1x,3(f18.5,1x))') 'dt_sub, dt_solv, data_step: ', dt_sub, dt_solv, data_step

  ! check that we have completed the sub-step
  if(dt_solv >= data_step-verySmall) then
   exit substeps
  endif

  ! adjust length of the sub-step (make sure that we don't exceed the step)
  dt_sub = min(data_step - dt_solv, dt_sub)
  !print*, 'dt_sub = ', dt_sub

 end do  substeps ! (sub-step loop)
 !print*, 'PAUSE: completed time step'; read(*,*)

 ! *** add snowfall to the snowpack...
 ! -----------------------------------

 ! add new snowfall to the snowpack
 ! NOTE: This needs to be done AFTER the call to canopySnow, since throughfall and unloading are computed in canopySnow
 call newsnwfall(&
                ! input: model control
                data_step,                                                 & ! time step (seconds)
                (nSnow > 0),                                               & ! logical flag if snow layers exist
                mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1),            & ! freeezing curve parameter for snow (K-1)
                ! input: diagnostic scalar variables
                diag_data%var(iLookDIAG%scalarSnowfallTemp)%dat(1),        & ! computed temperature of fresh snow (K)
                diag_data%var(iLookDIAG%scalarNewSnowDensity)%dat(1),      & ! computed density of new snow (kg m-3)
                flux_data%var(iLookFLUX%scalarThroughfallSnow)%dat(1),     & ! throughfall of snow through the canopy (kg m-2 s-1)
                flux_data%var(iLookFLUX%scalarCanopySnowUnloading)%dat(1), & ! unloading of snow from the canopy (kg m-2 s-1)
                ! input/output: state variables
                prog_data%var(iLookPROG%scalarSWE)%dat(1),                 & ! SWE (kg m-2)
                prog_data%var(iLookPROG%scalarSnowDepth)%dat(1),           & ! total snow depth (m)
                prog_data%var(iLookPROG%mLayerTemp)%dat(1),                & ! temperature of the top layer (K)
                prog_data%var(iLookPROG%mLayerDepth)%dat(1),               & ! depth of the top layer (m)
                prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1),          & ! volumetric fraction of ice of the top layer (-)
                prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1),          & ! volumetric fraction of liquid water of the top layer (-)
                ! output: error control
                err,cmessage)                                                ! error control
 if(err/=0)then; err=30; message=trim(message)//trim(cmessage); return; end if

 ! re-compute snow depth and SWE
 if(nSnow > 0)then
  prog_data%var(iLookPROG%scalarSnowDepth)%dat(1) = sum(  prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow))
  prog_data%var(iLookPROG%scalarSWE)%dat(1)       = sum( (prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                          prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                        * prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow) )
 end if
 !print*, 'SWE after snowfall = ',  prog_data%var(iLookPROG%scalarSWE)%dat(1)

 ! re-assign dimension lengths
 nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_snow)
 nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_soil)
 nLayers = nSnow+nSoil

 ! update coordinate variables
 call calcHeight(&
                 ! input/output: data structures
                 indx_data,   & ! intent(in): layer type
                 prog_data,   & ! intent(inout): model variables for a local HRU
                 ! output: error control
                 err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! overwrite flux_data with flux_mean (returns timestep-average fluxes for scalar variables)
 do iVar=1,size(averageFlux_meta)
  flux_data%var(averageFlux_meta(iVar)%ixParent)%dat(:) = flux_mean%var(iVar)%dat(:)
 end do

 ! ***********************************************************************************************************************************
 ! ***********************************************************************************************************************************
 ! ***********************************************************************************************************************************
 ! ***********************************************************************************************************************************

 ! ---
 ! *** balance checks...
 ! ---------------------

 ! save the average compression and melt pond storage in the data structures
 prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1)  = sfcMeltPond

 ! associate local variables with information in the data structures
 associate(&
 ! model forcing
 scalarSnowfall             => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarSnowfall)           )%dat(1)     ,&  ! computed snowfall rate (kg m-2 s-1)
 scalarRainfall             => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarRainfall)           )%dat(1)     ,&  ! computed rainfall rate (kg m-2 s-1)
 ! canopy fluxes
 averageThroughfallSnow     => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarThroughfallSnow)    )%dat(1)     ,&  ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
 averageThroughfallRain     => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarThroughfallRain)    )%dat(1)     ,&  ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 averageCanopySnowUnloading => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarCanopySnowUnloading))%dat(1)     ,&  ! unloading of snow from the vegetion canopy (kg m-2 s-1)
 averageCanopyLiqDrainage   => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarCanopyLiqDrainage)  )%dat(1)     ,&  ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 averageCanopySublimation   => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarCanopySublimation)  )%dat(1)     ,&  ! canopy sublimation/frost (kg m-2 s-1)
 averageCanopyEvaporation   => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarCanopyEvaporation)  )%dat(1)     ,&  ! canopy evaporation/condensation (kg m-2 s-1)
 ! snow fluxes
 averageSnowSublimation     => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarSnowSublimation)    )%dat(1)     ,&  ! sublimation from the snow surface (kg m-2 s-1)
 averageSnowDrainage        => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarSnowDrainage)       )%dat(1)     ,&  ! drainage from the bottom of the snowpack (m s-1)
 ! soil fluxes
 averageSoilInflux          => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarInfiltration)       )%dat(1)     ,&  ! influx of water at the top of the soil profile (m s-1)
 averageSoilDrainage        => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarSoilDrainage)       )%dat(1)     ,&  ! drainage from the bottom of the soil profile (m s-1)
 averageSoilBaseflow        => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarSoilBaseflow)       )%dat(1)     ,&  ! total baseflow from throughout the soil profile (m s-1)
 averageGroundEvaporation   => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarGroundEvaporation)  )%dat(1)     ,&  ! soil evaporation (kg m-2 s-1)
 averageCanopyTranspiration => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarCanopyTranspiration))%dat(1)     ,&  ! canopy transpiration (kg m-2 s-1)
 ! state variables in the vegetation canopy
 scalarCanopyLiq            => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)                               ,&  ! canopy liquid water (kg m-2)
 scalarCanopyIce            => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)                               ,&  ! canopy ice content (kg m-2)
 ! state variables in the soil domain
 mLayerDepth                => prog_data%var(iLookPROG%mLayerDepth)%dat(nSnow+1:nLayers)                     ,&  ! depth of each soil layer (m)
 mLayerVolFracIce           => prog_data%var(iLookPROG%mLayerVolFracIce)%dat(nSnow+1:nLayers)                ,&  ! volumetric ice content in each soil layer (-)
 mLayerVolFracLiq           => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(nSnow+1:nLayers)                ,&  ! volumetric liquid water content in each soil layer (-)
 scalarAquiferStorage       => prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)                          ,&  ! aquifer storage (m)
 ! error tolerance
 absConvTol_liquid          => mpar_data%var(iLookPARAM%absConvTol_liquid)%dat(1)                            ,&  ! absolute convergence tolerance for vol frac liq water (-)
 scalarTotalSoilIce         => diag_data%var(iLookDIAG%scalarTotalSoilIce)%dat(1)                            ,&  ! total ice in the soil column (kg m-2)
 scalarTotalSoilLiq         => diag_data%var(iLookDIAG%scalarTotalSoilLiq)%dat(1)                             &  ! total liquid water in the soil column (kg m-2)
 ) ! (association of local variables with information in the data structures

 ! -----
 ! * balance checks for the canopy...
 ! ----------------------------------

 ! if computing the vegetation flux
 if(computeVegFlux)then

  ! canopy water balance
  balanceCanopyWater1 = scalarCanopyLiq + scalarCanopyIce

  ! balance checks for the canopy
  ! NOTE: need to put the balance checks in the sub-step loop so that we can re-compute if necessary
  scalarCanopyWatBalError = balanceCanopyWater1 - (balanceCanopyWater0 + (scalarSnowfall - averageThroughfallSnow)*data_step + (scalarRainfall - averageThroughfallRain)*data_step &
                             - averageCanopySnowUnloading*data_step - averageCanopyLiqDrainage*data_step + averageCanopySublimation*data_step + averageCanopyEvaporation*data_step)
  if(abs(scalarCanopyWatBalError) > absConvTol_liquid*iden_water*10._rkind)then
   print*, '** canopy water balance error:'
   write(*,'(a,1x,f20.10)') 'data_step                                    = ', data_step
   write(*,'(a,1x,f20.10)') 'balanceCanopyWater0                          = ', balanceCanopyWater0
   write(*,'(a,1x,f20.10)') 'balanceCanopyWater1                          = ', balanceCanopyWater1
   write(*,'(a,1x,f20.10)') 'scalarSnowfall                               = ', scalarSnowfall
   write(*,'(a,1x,f20.10)') 'scalarRainfall                               = ', scalarRainfall
   write(*,'(a,1x,f20.10)') '(scalarSnowfall - averageThroughfallSnow)    = ', (scalarSnowfall - averageThroughfallSnow)!*data_step
   write(*,'(a,1x,f20.10)') '(scalarRainfall - averageThroughfallRain)    = ', (scalarRainfall - averageThroughfallRain)!*data_step
   write(*,'(a,1x,f20.10)') 'averageCanopySnowUnloading                   = ', averageCanopySnowUnloading!*data_step
   write(*,'(a,1x,f20.10)') 'averageCanopyLiqDrainage                     = ', averageCanopyLiqDrainage!*data_step
   write(*,'(a,1x,f20.10)') 'averageCanopySublimation                     = ', averageCanopySublimation!*data_step
   write(*,'(a,1x,f20.10)') 'averageCanopyEvaporation                     = ', averageCanopyEvaporation!*data_step
   write(*,'(a,1x,f20.10)') 'scalarCanopyWatBalError                      = ', scalarCanopyWatBalError
   message=trim(message)//'canopy hydrology does not balance'
   err=20; return
  end if

 endif  ! if computing the vegetation flux

 ! -----
 ! * balance checks for SWE...
 ! ---------------------------

 ! recompute snow depth (m) and SWE (kg m-2)
 if(nSnow > 0)then
  prog_data%var(iLookPROG%scalarSnowDepth)%dat(1) = sum(  prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow))
  prog_data%var(iLookPROG%scalarSWE)%dat(1)       = sum( (prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                          prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                        * prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow) )
 end if

 ! check the individual layers
 if(printBalance .and. nSnow>0)then
  write(*,'(a,1x,10(f12.8,1x))') 'liqSnowInit       = ', liqSnowInit
  write(*,'(a,1x,10(f12.8,1x))') 'volFracLiq        = ', prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)
  write(*,'(a,1x,10(f12.8,1x))') 'iLayerLiqFluxSnow = ', flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat*iden_water*data_step
  write(*,'(a,1x,10(f12.8,1x))') 'mLayerLiqFluxSnow = ', flux_data%var(iLookFLUX%mLayerLiqFluxSnow)%dat*data_step
  write(*,'(a,1x,10(f12.8,1x))') 'change volFracLiq = ', prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow) - liqSnowInit
  deallocate(liqSnowInit, stat=err)
  if(err/=0)then
   message=trim(message)//'unable to deallocate space for the initial volumetric liquid water content of snow'
   err=20; return
  endif
 endif

 ! check SWE
 if(nSnow>0)then
  effSnowfall = averageThroughfallSnow + averageCanopySnowUnloading
  effRainfall = averageThroughfallRain + averageCanopyLiqDrainage
  newSWE      = prog_data%var(iLookPROG%scalarSWE)%dat(1)
  delSWE      = newSWE - (oldSWE - sfcMeltPond)
  massBalance = delSWE - (effSnowfall + effRainfall + averageSnowSublimation - averageSnowDrainage*iden_water)*data_step
  if(abs(massBalance) > 1.d-6)then
   print*,                  'nSnow       = ', nSnow
   print*,                  'nSub        = ', nSub
   write(*,'(a,1x,f20.10)') 'data_step   = ', data_step
   write(*,'(a,1x,f20.10)') 'oldSWE      = ', oldSWE
   write(*,'(a,1x,f20.10)') 'newSWE      = ', newSWE
   write(*,'(a,1x,f20.10)') 'delSWE      = ', delSWE
   write(*,'(a,1x,f20.10)') 'effRainfall = ', effRainfall*data_step
   write(*,'(a,1x,f20.10)') 'effSnowfall = ', effSnowfall*data_step
   write(*,'(a,1x,f20.10)') 'sublimation = ', averageSnowSublimation*data_step
   write(*,'(a,1x,f20.10)') 'snwDrainage = ', averageSnowDrainage*iden_water*data_step
   write(*,'(a,1x,f20.10)') 'sfcMeltPond = ', sfcMeltPond
   write(*,'(a,1x,f20.10)') 'massBalance = ', massBalance
   message=trim(message)//'SWE does not balance'
   err=20; return
  endif  ! if failed mass balance check
 endif  ! if snow layers exist

 ! -----
 ! * balance checks for soil...
 ! ----------------------------

 ! compute the liquid water and ice content at the end of the time step
 scalarTotalSoilLiq = sum(iden_water*mLayerVolFracLiq(1:nSoil)*mLayerDepth(1:nSoil))
 scalarTotalSoilIce = sum(iden_water*mLayerVolFracIce(1:nSoil)*mLayerDepth(1:nSoil))   ! NOTE: no expansion of soil, hence use iden_water

 ! get the total water in the soil (liquid plus ice) at the end of the time step (kg m-2)
 balanceSoilWater1 = scalarTotalSoilLiq + scalarTotalSoilIce

 ! get the total aquifer storage at the start of the time step (kg m-2)
 balanceAquifer1 = scalarAquiferStorage*iden_water

 ! get the input and output to/from the soil zone (kg m-2)
 balanceSoilInflux        = averageSoilInflux*iden_water*data_step
 balanceSoilBaseflow      = averageSoilBaseflow*iden_water*data_step
 balanceSoilDrainage      = averageSoilDrainage*iden_water*data_step
 balanceSoilET            = (averageCanopyTranspiration + averageGroundEvaporation)*data_step

 ! check the individual layers
 if(printBalance)then
  write(*,'(a,1x,10(f12.8,1x))') 'liqSoilInit       = ', liqSoilInit
  write(*,'(a,1x,10(f12.8,1x))') 'volFracLiq        = ', mLayerVolFracLiq
  write(*,'(a,1x,10(f12.8,1x))') 'iLayerLiqFluxSoil = ', flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat*iden_water*data_step
  write(*,'(a,1x,10(f12.8,1x))') 'mLayerLiqFluxSoil = ', flux_data%var(iLookFLUX%mLayerLiqFluxSoil)%dat*data_step
  write(*,'(a,1x,10(f12.8,1x))') 'change volFracLiq = ', mLayerVolFracLiq - liqSoilInit
  deallocate(liqSoilInit, stat=err)
  if(err/=0)then
   message=trim(message)//'unable to deallocate space for the initial soil moisture'
   err=20; return
  endif
 endif

 ! check the soil water balance
 scalarSoilWatBalError  = balanceSoilWater1 - (balanceSoilWater0 + (balanceSoilInflux + balanceSoilET - balanceSoilBaseflow - balanceSoilDrainage - totalSoilCompress) )
 if(abs(scalarSoilWatBalError) > absConvTol_liquid*iden_water*10._rkind)then  ! NOTE: kg m-2, so need coarse tolerance to account for precision issues
  write(*,*)               'solution method           = ', ixSolution
  write(*,'(a,1x,f20.10)') 'data_step                 = ', data_step
  write(*,'(a,1x,f20.10)') 'totalSoilCompress         = ', totalSoilCompress
  write(*,'(a,1x,f20.10)') 'scalarTotalSoilLiq        = ', scalarTotalSoilLiq
  write(*,'(a,1x,f20.10)') 'scalarTotalSoilIce        = ', scalarTotalSoilIce
  write(*,'(a,1x,f20.10)') 'balanceSoilWater0         = ', balanceSoilWater0
  write(*,'(a,1x,f20.10)') 'balanceSoilWater1         = ', balanceSoilWater1
  write(*,'(a,1x,f20.10)') 'balanceSoilInflux         = ', balanceSoilInflux
  write(*,'(a,1x,f20.10)') 'balanceSoilBaseflow       = ', balanceSoilBaseflow
  write(*,'(a,1x,f20.10)') 'balanceSoilDrainage       = ', balanceSoilDrainage
  write(*,'(a,1x,f20.10)') 'balanceSoilET             = ', balanceSoilET
  write(*,'(a,1x,f20.10)') 'scalarSoilWatBalError     = ', scalarSoilWatBalError
  write(*,'(a,1x,f20.10)') 'scalarSoilWatBalError     = ', scalarSoilWatBalError/iden_water
  write(*,'(a,1x,f20.10)') 'absConvTol_liquid         = ', absConvTol_liquid
  ! error control
  message=trim(message)//'soil hydrology does not balance'
  err=20; return
 end if

 ! end association of local variables with information in the data structures
 end associate

 ! end association to canopy depth
 end associate canopy

 ! Save the total soil water (Liquid+Ice)
 diag_data%var(iLookDIAG%scalarTotalSoilWat)%dat(1) = balanceSoilWater1
 ! save the surface temperature (just to make things easier to visualize)
 prog_data%var(iLookPROG%scalarSurfaceTemp)%dat(1) = prog_data%var(iLookPROG%mLayerTemp)%dat(1)

 ! overwrite flux data with the timestep-average value
 if(.not.backwardsCompatibility)then
  do iVar=1,size(flux_mean%var)
   flux_data%var(averageFlux_meta(iVar)%ixParent)%dat = flux_mean%var(iVar)%dat
  end do
 end if

 iLayer = nSnow+1
 !print*, 'nsub, mLayerTemp(iLayer), mLayerVolFracIce(iLayer) = ', nsub, mLayerTemp(iLayer), mLayerVolFracIce(iLayer)
 !print*, 'nsub = ', nsub
 if(nsub>50000)then
  write(message,'(a,i0)') trim(cmessage)//'number of sub-steps > 50000 for HRU ', hruID
  err=20; return
 end if

 ! get the end time
 call cpu_time(endTime)

 ! get the elapsed time
 diag_data%var(iLookDIAG%wallClockTime)%dat(1) = endTime - startTime

 end subroutine coupled_em


 ! *********************************************************************************************************
 ! private subroutine implctMelt: compute melt of the "snow without a layer"
 ! *********************************************************************************************************
 subroutine implctMelt(&
                       ! input/output: integrated snowpack properties
                       scalarSWE,         & ! intent(inout): snow water equivalent (kg m-2)
                       scalarSnowDepth,   & ! intent(inout): snow depth (m)
                       scalarSfcMeltPond, & ! intent(inout): surface melt pond (kg m-2)
                       ! input/output: properties of the upper-most soil layer
                       soilTemp,          & ! intent(inout): surface layer temperature (K)
                       soilDepth,         & ! intent(inout): surface layer depth (m)
                       soilHeatcap,       & ! intent(inout): surface layer volumetric heat capacity (J m-3 K-1)
                       ! output: error control
                       err,message        ) ! intent(out): error control
 implicit none
 ! input/output: integrated snowpack properties
 real(rkind),intent(inout)    :: scalarSWE          ! snow water equivalent (kg m-2)
 real(rkind),intent(inout)    :: scalarSnowDepth    ! snow depth (m)
 real(rkind),intent(inout)    :: scalarSfcMeltPond  ! surface melt pond (kg m-2)
 ! input/output: properties of the upper-most soil layer
 real(rkind),intent(inout)    :: soilTemp           ! surface layer temperature (K)
 real(rkind),intent(inout)    :: soilDepth          ! surface layer depth (m)
 real(rkind),intent(inout)    :: soilHeatcap        ! surface layer volumetric heat capacity (J m-3 K-1)
 ! output: error control
 integer(i4b),intent(out)  :: err                ! error code
 character(*),intent(out)  :: message            ! error message
 ! local variables
 real(rkind)                  :: nrgRequired        ! energy required to melt all the snow (J m-2)
 real(rkind)                  :: nrgAvailable       ! energy available to melt the snow (J m-2)
 real(rkind)                  :: snwDensity         ! snow density (kg m-3)
 ! initialize error control
 err=0; message='implctMelt/'

 if(scalarSWE > 0._rkind)then
  ! only melt if temperature of the top soil layer is greater than Tfreeze
  if(soilTemp > Tfreeze)then
   ! compute the energy required to melt all the snow (J m-2)
   nrgRequired     = scalarSWE*LH_fus
   ! compute the energy available to melt the snow (J m-2)
   nrgAvailable    = soilHeatcap*(soilTemp - Tfreeze)*soilDepth
   ! compute the snow density (not saved)
   snwDensity      = scalarSWE/scalarSnowDepth
   ! compute the amount of melt, and update SWE (kg m-2)
   if(nrgAvailable > nrgRequired)then
    scalarSfcMeltPond  = scalarSWE
    scalarSWE          = 0._rkind
   else
    scalarSfcMeltPond  = nrgAvailable/LH_fus
    scalarSWE          = scalarSWE - scalarSfcMeltPond
   end if
   ! update depth
   scalarSnowDepth = scalarSWE/snwDensity
   ! update temperature of the top soil layer (K)
   soilTemp =  soilTemp - (LH_fus*scalarSfcMeltPond/soilDepth)/soilHeatcap
  else  ! melt is zero if the temperature of the top soil layer is less than Tfreeze
   scalarSfcMeltPond = 0._rkind  ! kg m-2
  end if ! (if the temperature of the top soil layer is greater than Tfreeze)
 else  ! melt is zero if the "snow without a layer" does not exist
  scalarSfcMeltPond = 0._rkind  ! kg m-2
 end if ! (if the "snow without a layer" exists)

 end subroutine implctMelt

end module coupled_em_module
