! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
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

module eval8summa_module

! data types
USE nrtype

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number
USE globalData,only:quadMissing     ! missing quadruple precision number

! access the global print flag
USE globalData,only:globalPrintFlag

! define access to state variables to print
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! domain types
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    LH_vap,       & ! latent heat of vaporization          (J kg-1)
                    LH_sub,       & ! latent heat of sublimation           (J kg-1)
                    Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                    iden_air,     & ! intrinsic density of air             (kg m-3)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    zLookup,      & ! data vector with variable length dimension (dp)
                    model_options   ! defines the model decisions

! indices that define elements of the data structures
USE var_lookup,only:iLookDECISIONS               ! named variables for elements of the decision structure
USE var_lookup,only:iLookPARAM                   ! named variables for structure elements
USE var_lookup,only:iLookPROG                    ! named variables for structure elements
USE var_lookup,only:iLookINDEX                   ! named variables for structure elements
USE var_lookup,only:iLookDIAG                    ! named variables for structure elements
USE var_lookup,only:iLookFLUX                    ! named variables for structure elements
USE var_lookup,only:iLookDERIV                   ! named variables for structure elements

! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:  &
 localColumn,                & ! separate groundwater representation in each local soil column
 singleBasin                   ! single groundwater store over the entire basin

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:  &
 qbaseTopmodel,              & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                  & ! a big bucket (lumped aquifer model)
 noExplicit                    ! no explicit groundwater parameterization

! look-up values for the form of Richards' equation
USE mDecisions_module,only:  &
 moisture,                   & ! moisture-based form of Richards' equation
 mixdform                      ! mixed form of Richards' equation

implicit none
private
public::eval8summa

contains

 ! **********************************************************************************************************
 ! public subroutine eval8summa: compute the residual vector and the Jacobian matrix
 ! **********************************************************************************************************
 subroutine eval8summa(&
                       ! input: model control
                       dt,                      & ! intent(in):    length of the time step (seconds)
                       nSnow,                   & ! intent(in):    number of snow layers
                       nSoil,                   & ! intent(in):    number of soil layers
                       nLayers,                 & ! intent(in):    total number of layers
                       nState,                  & ! intent(in):    total number of state variables
                       firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                       firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                       firstSplitOper,          & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                       computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                       scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                       ! input: state vectors
                       stateVecTrial,           & ! intent(in):    model state vector
                       fScale,                  & ! intent(in):    function scaling vector
                       sMul,                    & ! intent(inout): state vector multiplier (used in the residual calculations)
                       ! input: data structures
                       model_decisions,         & ! intent(in):    model decisions
                       lookup_data,             & ! intent(in):    lookup tables
                       type_data,               & ! intent(in):    type of vegetation and soil
                       attr_data,               & ! intent(in):    spatial attributes
                       mpar_data,               & ! intent(in):    model parameters
                       forc_data,               & ! intent(in):    model forcing data
                       bvar_data,               & ! intent(in):    average model variables for the entire basin
                       prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                       ! input-output: data structures
                       indx_data,               & ! intent(inout): index data
                       diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                       flux_data,               & ! intent(inout): model fluxes for a local HRU
                       deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                       ! input-output: baseflow
                       ixSaturation,            & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                       dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                       ! output: flux and residual vectors
                       feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                       fluxVec,                 & ! intent(out):   flux vector
                       resSink,                 & ! intent(out):   additional (sink) terms on the RHS of the state equation
                       resVec,                  & ! intent(out):   residual vector
                       fEval,                   & ! intent(out):   function evaluation
                       err,message)               ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! provide access to subroutines
 USE getVectorz_module, only:varExtract           ! extract variables from the state vector
 USE updateVars_module, only:updateVars           ! update prognostic variables
 USE t2enthalpy_module, only:t2enthalpy           ! compute enthalpy
 USE t2enthalpy_module, only:t2enthalpy_T         ! compute enthalpy
 USE computHeatCap_module,only:computHeatCap      ! compute diagnostic energy variables -- thermal conductivity and heat capacity
 USE computHeatCap_module, only:computStatMult    ! state vector multiplier (used in the residual calculations)
 USE computThermConduct_module,only:computThermConduct    ! compute thermal conductivity
 USE computFlux_module, only:soilCmpres           ! compute soil compression
 USE computFlux_module, only:computFlux           ! compute fluxes given a state vector
 USE computResid_module,only:computResid          ! compute residuals given a state vector
 implicit none
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 real(rkind),intent(in)             :: dt                     ! length of the time step (seconds)
 integer(i4b),intent(in)         :: nSnow                  ! number of snow layers
 integer(i4b),intent(in)         :: nSoil                  ! number of soil layers
 integer(i4b),intent(in)         :: nLayers                ! total number of layers
 integer(i4b),intent(in)         :: nState                 ! total number of state variables
 logical(lgt),intent(in)         :: firstSubStep           ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(inout)      :: firstFluxCall          ! flag to indicate if we are processing the first flux call
 logical(lgt),intent(in)         :: firstSplitOper         ! flag to indicate if we are processing the first flux call in a splitting operation
 logical(lgt),intent(in)         :: computeVegFlux         ! flag to indicate if computing fluxes over vegetation
 logical(lgt),intent(in)         :: scalarSolution         ! flag to denote if implementing the scalar solution
 ! input: state vectors
 real(rkind),intent(in)             :: stateVecTrial(:)       ! model state vector
 real(rkind),intent(in)             :: fScale(:)              ! function scaling vector
 real(rkind),intent(inout)          :: sMul(:)   ! NOTE: qp   ! state vector multiplier (used in the residual calculations)
 ! input: data structures
 type(model_options),intent(in)  :: model_decisions(:)     ! model decisions
 type(zLookup),      intent(in)  :: lookup_data            ! lookup tables
 type(var_i),        intent(in)  :: type_data              ! type of vegetation and soil
 type(var_d),        intent(in)  :: attr_data              ! spatial attributes
 type(var_dlength),  intent(in)  :: mpar_data              ! model parameters
 type(var_d),        intent(in)  :: forc_data              ! model forcing data
 type(var_dlength),  intent(in)  :: bvar_data              ! model variables for the local basin
 type(var_dlength),  intent(in)  :: prog_data              ! prognostic variables for a local HRU
 ! output: data structures
 type(var_ilength),intent(inout) :: indx_data              ! indices defining model states and layers
 type(var_dlength),intent(inout) :: diag_data              ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data              ! model fluxes for a local HRU
 type(var_dlength),intent(inout) :: deriv_data             ! derivatives in model fluxes w.r.t. relevant state variables
 ! input-output: baseflow
 integer(i4b),intent(inout)      :: ixSaturation           ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
 real(rkind),intent(out)            :: dBaseflow_dMatric(:,:) ! derivative in baseflow w.r.t. matric head (s-1)
 ! output: flux and residual vectors
 logical(lgt),intent(out)        :: feasible               ! flag to denote the feasibility of the solution
 real(rkind),intent(out)            :: fluxVec(:)             ! flux vector
 real(rkind),intent(out)            :: resSink(:)             ! sink terms on the RHS of the flux equation
 real(rkind),intent(out)            :: resVec(:) ! NOTE: qp   ! residual vector
 real(rkind),intent(out)            :: fEval                  ! function evaluation
 ! output: error control
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! state variables
 real(rkind)                        :: scalarCanairTempTrial     ! trial value for temperature of the canopy air space (K)
 real(rkind)                        :: scalarCanopyTempTrial     ! trial value for temperature of the vegetation canopy (K)
 real(rkind)                        :: scalarCanopyWatTrial      ! trial value for liquid water storage in the canopy (kg m-2)
 real(rkind),dimension(nLayers)     :: mLayerTempTrial           ! trial value for temperature of layers in the snow and soil domains (K)
 real(rkind),dimension(nLayers)     :: mLayerVolFracWatTrial     ! trial value for volumetric fraction of total water (-)
 real(rkind),dimension(nSoil)       :: mLayerMatricHeadTrial     ! trial value for total water matric potential (m)
 real(rkind),dimension(nSoil)       :: mLayerMatricHeadLiqTrial  ! trial value for liquid water matric potential (m)
 real(rkind)                        :: scalarAquiferStorageTrial ! trial value of storage of water in the aquifer (m)
 ! enthalpy
 logical(lgt),parameter             :: needEnthalpy=.true.      ! flag to compute enthalpy
 ! diagnostic variables
 real(rkind)                        :: scalarCanopyLiqTrial      ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
 real(rkind)                        :: scalarCanopyIceTrial      ! trial value for mass of ice on the vegetation canopy (kg m-2)
 real(rkind),dimension(nLayers)     :: mLayerVolFracLiqTrial     ! trial value for volumetric fraction of liquid water (-)
 real(rkind),dimension(nLayers)     :: mLayerVolFracIceTrial     ! trial value for volumetric fraction of ice (-)
 ! other local variables
 integer(i4b)                    :: iLayer                    ! index of model layer in the snow+soil domain
 integer(i4b)                    :: jState(1)                 ! index of model state for the scalar solution within the soil domain
 integer(i4b)                    :: ixBeg,ixEnd               ! index of indices for the soil compression routine
 integer(i4b),parameter          :: ixVegVolume=1             ! index of the desired vegetation control volumne (currently only one veg layer)
 real(rkind)                        :: xMin,xMax                 ! minimum and maximum values for water content
 real(rkind)                        :: scalarCanopyHydTrial      ! trial value for mass of water on the vegetation canopy (kg m-2)
 real(rkind),parameter              :: canopyTempMax=500._rkind     ! expected maximum value for the canopy temperature (K)
 real(rkind),dimension(nLayers)     :: mLayerVolFracHydTrial     ! trial value for volumetric fraction of water (-), general vector merged from Wat and Liq
 real(rkind),dimension(nState)      :: rVecScaled                ! scaled residual vector
 character(LEN=256)              :: cmessage                  ! error message of downwind routine
 real(rkind)			                   :: scalarCanopyEnthalpyTrial ! enthalpy of the vegetation canopy (J m-3)
 real(rkind)                        :: mLayerEnthalpyTrial(nLayers) ! enthalpy of snow+soil layer (J m-3)
 real(rkind),dimension(nLayers)     :: mLayerHeatCapTrial        ! heat capacity of each snow+soil layer
 real(qp)                        :: heatCapVegTrial           ! heat capacity of vegetation canopy
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! association to variables in the data structures
 ! --------------------------------------------------------------------------------------------------------------------------------
 associate(&
  scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)      ,& ! intent(in): [dp]     temperature of the vegetation canopy (K)
  mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat               ,& ! intent(in): [dp(:)]  temperature of each snow/soil layer (K)
 ! model decisions
 ixRichards              => model_decisions(iLookDECISIONS%f_Richards)%iDecision   ,&  ! intent(in):  [i4b]   index of the form of Richards' equation
 ! snow parameters
 snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)         ,&  ! intent(in):  [dp]    scaling parameter for the snow freezing curve (K-1)
 ! soil parameters
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat                ,&  ! intent(in):  [dp(:)] soil porosity (-)
 specificStorage         => mpar_data%var(iLookPARAM%specificStorage)%dat(1)       ,&  ! intent(in):  [dp]    specific storage coefficient (m-1)
 ! canopy and layer depth
 canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,&  ! intent(in):  [dp   ] canopy depth (m)
 mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,&  ! intent(in):  [dp(:)] depth of each layer in the snow-soil sub-domain (m)
 ! model state variables
 scalarSfcMeltPond       => prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1)      ,&  ! intent(in):  [dp]    ponded water caused by melt of the "snow without a layer" (kg m-2)
 mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,&  ! intent(in):  [dp(:)] volumetric fraction of liquid water (-)
 mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,&  ! intent(in):  [dp(:)] volumetric fraction of ice (-)
 mLayerMatricHeadLiq     => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat       ,&  ! intent(in):  [dp(:)] liquid water matric potential (m)
 ! model diagnostic variables
 scalarFracLiqVeg        => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)       ,&  ! intent(in):  [dp]    fraction of liquid water on vegetation (-)
 mLayerFracLiqSnow       => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat         ,&  ! intent(in):  [dp(:)] fraction of liquid water in each snow layer (-)
! enthalpy
 scalarCanairEnthalpy    => diag_data%var(iLookDIAG%scalarCanairEnthalpy)%dat(1)   ,&  ! intent(in): [dp]    enthalpy of the canopy air space (J m-3)
 scalarCanopyEnthalpy    => diag_data%var(iLookDIAG%scalarCanopyEnthalpy)%dat(1)   ,&  ! intent(in): [dp]    enthalpy of the vegetation canopy (J m-3)
 mLayerEnthalpy          => diag_data%var(iLookDIAG%mLayerEnthalpy)%dat            ,&  ! intent(in): [dp(:)] enthalpy of the snow+soil layers (J m-3)
 ! soil compression
 scalarSoilCompress      => diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1)     ,&  ! intent(in): [dp]    total change in storage associated with compression of the soil matrix (kg m-2)
 mLayerCompress          => diag_data%var(iLookDIAG%mLayerCompress)%dat            ,&  ! intent(in): [dp(:)] change in storage associated with compression of the soil matrix (-)
 ! derivatives
 dVolTot_dPsi0           => deriv_data%var(iLookDERIV%dVolTot_dPsi0)%dat           ,&  ! intent(in): [dp(:)] derivative in total water content w.r.t. total water matric potential
 dCompress_dPsi          => deriv_data%var(iLookDERIV%dCompress_dPsi)%dat          ,&  ! intent(in): [dp(:)] derivative in compressibility w.r.t. matric head (m-1)
 ! mapping
 ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat         ,&  ! intent(in): [i4b(:)] mapping of full state vector to the state subset
 ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat          ,&  ! intent(in): [i4b(:)] index of control volume for different domains (veg, snow, soil)
 ! indices
 ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,&  ! intent(in): [i4b]    index of canopy air space energy state variable (nrg)
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,&  ! intent(in): [i4b]    index of canopy energy state variable (nrg)
 ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,&  ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
 ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow subdomain
 ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,&  ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
 ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,&  ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
 ixHydCanopy             => indx_data%var(iLookINDEX%ixHydCanopy)%dat              ,&  ! intent(in): [i4b(:)] index of the hydrology states in the canopy domain
 ixHydType               => indx_data%var(iLookINDEX%ixHydType)%dat                ,&  ! intent(in): [i4b(:)] index of the type of hydrology states in snow+soil domain
 layerType               => indx_data%var(iLookINDEX%layerType)%dat                 &  ! intent(in): [i4b(:)] layer type (iname_soil or iname_snow)
 ) ! association to variables in the data structures
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="eval8summa/"

 ! check the feasibility of the solution
 feasible=.true.

 ! check that the canopy air space temperature is reasonable
 if(ixCasNrg/=integerMissing)then
  if(stateVecTrial(ixCasNrg) > canopyTempMax) feasible=.false.
 endif

 ! check that the canopy air space temperature is reasonable
 if(ixVegNrg/=integerMissing)then
  if(stateVecTrial(ixVegNrg) > canopyTempMax) feasible=.false.
 endif

 ! check canopy liquid water is not negative
 if(ixVegHyd/=integerMissing)then
  if(stateVecTrial(ixVegHyd) < 0._rkind) feasible=.false.
 end if

 ! check snow temperature is below freezing
 if(count(ixSnowOnlyNrg/=integerMissing)>0)then
  if(any(stateVecTrial( pack(ixSnowOnlyNrg,ixSnowOnlyNrg/=integerMissing) ) > Tfreeze)) feasible=.false.
 endif

 ! loop through non-missing hydrology state variables in the snow+soil domain
 do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)

  ! check the minimum and maximum water constraints
  if(ixHydType(iLayer)==iname_watLayer .or. ixHydType(iLayer)==iname_liqLayer)then

   ! --> minimum
   if (layerType(iLayer) == iname_soil) then
    xMin = theta_sat(iLayer-nSnow)
   else
    xMin = 0._rkind
   endif

   ! --> maximum
   select case( layerType(iLayer) )
    case(iname_snow); xMax = merge(iden_ice,  1._rkind - mLayerVolFracIce(iLayer), ixHydType(iLayer)==iname_watLayer)
    case(iname_soil); xMax = merge(theta_sat(iLayer-nSnow), theta_sat(iLayer-nSnow) - mLayerVolFracIce(iLayer), ixHydType(iLayer)==iname_watLayer)
   end select

   ! --> check
   if(stateVecTrial( ixSnowSoilHyd(iLayer) ) < xMin .or. stateVecTrial( ixSnowSoilHyd(iLayer) ) > xMax) feasible=.false.
   !if(.not.feasible) write(*,'(a,1x,i4,1x,L1,1x,10(f20.10,1x))') 'iLayer, feasible, stateVecTrial( ixSnowSoilHyd(iLayer) ), xMin, xMax = ', iLayer, feasible, stateVecTrial( ixSnowSoilHyd(iLayer) ), xMin, xMax

  endif  ! if water states

 end do  ! loop through non-missing hydrology state variables in the snow+soil domain

 ! early return for non-feasible solutions
 if(.not.feasible)then
  fluxVec(:) = realMissing
  resVec(:)  = quadMissing
  fEval      = realMissing
  return
 end if

 ! get the start and end indices for the soil compression calculations
 if(scalarSolution)then
  jState = pack(ixControlVolume, ixMapFull2Subset/=integerMissing)
  ixBeg  = jState(1)
  ixEnd  = jState(1)
 else
  ixBeg  = 1
  ixEnd  = nSoil
 endif

 ! extract variables from the model state vector
 call varExtract(&
                 ! input
                 stateVecTrial,            & ! intent(in):    model state vector (mixed units)
                 diag_data,                & ! intent(in):    model diagnostic variables for a local HRU
                 prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                 indx_data,                & ! intent(in):    indices defining model states and layers
                 ! output: variables for the vegetation canopy
                 scalarCanairTempTrial,    & ! intent(out):   trial value of canopy air temperature (K)
                 scalarCanopyTempTrial,    & ! intent(out):   trial value of canopy temperature (K)
                 scalarCanopyWatTrial,     & ! intent(out):   trial value of canopy total water (kg m-2)
                 scalarCanopyLiqTrial,     & ! intent(out):   trial value of canopy liquid water (kg m-2)
                 scalarCanopyIceTrial,     & ! intent(out):   trial value of canopy ice content (kg m-2)
                 ! output: variables for the snow-soil domain
                 mLayerTempTrial,          & ! intent(out):   trial vector of layer temperature (K)
                 mLayerVolFracWatTrial,    & ! intent(out):   trial vector of volumetric total water content (-)
                 mLayerVolFracLiqTrial,    & ! intent(out):   trial vector of volumetric liquid water content (-)
                 mLayerVolFracIceTrial,    & ! intent(out):   trial vector of volumetric ice water content (-)
                 mLayerMatricHeadTrial,    & ! intent(out):   trial vector of total water matric potential (m)
                 mLayerMatricHeadLiqTrial, & ! intent(out):   trial vector of liquid water matric potential (m)
                 ! output: variables for the aquifer
                 scalarAquiferStorageTrial,& ! intent(out):   trial value of storage of water in the aquifer (m)
                 ! output: error control
                 err,cmessage)               ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

 ! update diagnostic variables
 call updateVars(&
                 ! input
                 .false.,                                   & ! intent(in):    logical flag to adjust temperature to account for the energy used in melt+freeze
                 lookup_data,                               & ! intent(in):    lookup tables for a local HRU
                 mpar_data,                                 & ! intent(in):    model parameters for a local HRU
                 indx_data,                                 & ! intent(in):    indices defining model states and layers
                 prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                 diag_data,                                 & ! intent(inout): model diagnostic variables for a local HRU
                 deriv_data,                                & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                 ! output: variables for the vegetation canopy
                 scalarCanopyTempTrial,                     & ! intent(inout): trial value of canopy temperature (K)
                 scalarCanopyWatTrial,                      & ! intent(inout): trial value of canopy total water (kg m-2)
                 scalarCanopyLiqTrial,                      & ! intent(inout): trial value of canopy liquid water (kg m-2)
                 scalarCanopyIceTrial,                      & ! intent(inout): trial value of canopy ice content (kg m-2)
                 ! output: variables for the snow-soil domain
                 mLayerTempTrial,                           & ! intent(inout): trial vector of layer temperature (K)
                 mLayerVolFracWatTrial,                     & ! intent(inout): trial vector of volumetric total water content (-)
                 mLayerVolFracLiqTrial,                     & ! intent(inout): trial vector of volumetric liquid water content (-)
                 mLayerVolFracIceTrial,                     & ! intent(inout): trial vector of volumetric ice water content (-)
                 mLayerMatricHeadTrial,                     & ! intent(inout): trial vector of total water matric potential (m)
                 mLayerMatricHeadLiqTrial,                  & ! intent(inout): trial vector of liquid water matric potential (m)
                 ! output: error control
                 err,cmessage)                                ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)
 
 ! print the states in the canopy domain
 !print*, 'dt = ', dt
 !write(*,'(a,1x,10(f20.10,1x))') 'scalarCanopyTempTrial    = ', scalarCanopyTempTrial
 !write(*,'(a,1x,10(f20.10,1x))') 'scalarCanopyWatTrial     = ', scalarCanopyWatTrial
 !write(*,'(a,1x,10(f20.10,1x))') 'scalarCanopyLiqTrial     = ', scalarCanopyLiqTrial
 !write(*,'(a,1x,10(f20.10,1x))') 'scalarCanopyIceTrial     = ', scalarCanopyIceTrial

 ! print the states in the snow+soil domain
 !write(*,'(a,1x,10(f20.10,1x))') 'mLayerTempTrial          = ', mLayerTempTrial(iJac1:min(nLayers,iJac2))
 !write(*,'(a,1x,10(f20.10,1x))') 'mLayerVolFracWatTrial    = ', mLayerVolFracWatTrial(iJac1:min(nLayers,iJac2))
 !write(*,'(a,1x,10(f20.10,1x))') 'mLayerVolFracLiqTrial    = ', mLayerVolFracLiqTrial(iJac1:min(nLayers,iJac2))
 !write(*,'(a,1x,10(f20.10,1x))') 'mLayerVolFracIceTrial    = ', mLayerVolFracIceTrial(iJac1:min(nLayers,iJac2))
 !write(*,'(a,1x,10(f20.10,1x))') 'mLayerMatricHeadTrial    = ', mLayerMatricHeadTrial(iJac1:min(nSoil,iJac2))
 !write(*,'(a,1x,10(f20.10,1x))') 'mLayerMatricHeadLiqTrial = ', mLayerMatricHeadLiqTrial(iJac1:min(nSoil,iJac2))

 ! print the water content
 if(globalPrintFlag)then
  if(iJac1<nSnow) write(*,'(a,10(f16.10,1x))') 'mLayerVolFracWatTrial = ', mLayerVolFracWatTrial(iJac1:min(iJac2,nSnow))
  if(iJac1<nSnow) write(*,'(a,10(f16.10,1x))') 'mLayerVolFracLiqTrial = ', mLayerVolFracLiqTrial(iJac1:min(iJac2,nSnow))
  if(iJac1<nSnow) write(*,'(a,10(f16.10,1x))') 'mLayerVolFracIceTrial = ', mLayerVolFracIceTrial(iJac1:min(iJac2,nSnow))
 endif


  ! *******************************************************************************************************
 ! *******************************************************************************************************
 ! ******************************************************************************************************* 
 if(needEnthalpy)then
  ! compute H_T
  call t2enthalpy_T(&
                  ! input: data structures
                  diag_data,                   & ! intent(in):  model diagnostic variables for a local HRU
                  mpar_data,                   & ! intent(in):  parameter data structure
                  indx_data,                   & ! intent(in):  model indices
                  lookup_data,                 & ! intent(in):  lookup table data structure
                  ! input: state variables for the vegetation canopy
                  scalarCanairTempTrial,       & ! intent(in):  trial value of canopy air temperature (K)
                  scalarCanopyTempTrial,       & ! intent(in):  trial value of canopy temperature (K)
                  scalarCanopyWatTrial,        & ! intent(in):  trial value of canopy total water (kg m-2)
                  scalarCanopyIceTrial,        & ! intent(in):  trial value for canopy ice content (kg m-2)
                  ! input: variables for the snow-soil domain
                  mLayerTempTrial,             & ! intent(in):  trial vector of layer temperature (K)
                  mLayerVolFracWatTrial,       & ! intent(in):  trial vector of volumetric total water content (-)
                  mLayerMatricHeadTrial,       & ! intent(in):  trial vector of total water matric potential (m)
                  mLayerVolFracIceTrial,       & ! intent(in):  trial vector of volumetric fraction of ice (-)
                  ! output: enthalpy
                  scalarCanairEnthalpy,        & ! intent(out):  enthalpy of the canopy air space (J m-3)
                  scalarCanopyEnthalpyTrial,   & ! intent(out):  enthalpy of the vegetation canopy (J m-3)
                  mLayerEnthalpyTrial,         & ! intent(out):  enthalpy of each snow+soil layer (J m-3)
                  ! output: error control
                  err,cmessage)                  ! intent(out): error control
     if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
     
     
   ! *** compute volumetric heat capacity C_p = dH_T/dT
   call computHeatCap(&
                       ! input: control variables
                       nLayers,                 	 & ! intent(in): number of layers (soil+snow)
                       computeVegFlux,         		 & ! intent(in): flag to denote if computing the vegetation flux
                       canopyDepth,             	 & ! intent(in): canopy depth (m)
                       ! input data structures
                       mpar_data,               	 & ! intent(in): model parameters
                       indx_data,               	 & ! intent(in): model layer indices
                       diag_data,               	 & ! intent(in): model diagnostic variables for a local HRU
                       ! input: state variables
                       scalarCanopyIceTrial,       & ! intent(in): trial value for canopy ice content (kg m-2)
                       scalarCanopyLiqTrial,       & ! intent(in): trial value for the liquid water on the vegetation canopy (kg m-2)
                       scalarCanopyTempTrial,   	 & ! intent(in): trial value of canopy temperature (K)
                       scalarCanopyTemp,  	    	 & ! intent(in): previous value of canopy temperature (K)
                       scalarCanopyEnthalpyTrial,  & ! intent(in): trial enthalpy of the vegetation canopy (J m-3)
                       scalarCanopyEnthalpy,       & ! intent(in): previous enthalpy of the vegetation canopy (J m-3)
                       mLayerVolFracIceTrial,      & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                       mLayerVolFracLiqTrial,      & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
                       mLayerTempTrial,          	 & ! intent(in): trial temperature
                       mLayerTemp,		           	 & ! intent(in): previous temperature
                       mLayerEnthalpyTrial,      	 & ! intent(in): trial enthalpy for snow and soil
                       mLayerEnthalpy,		       	 & ! intent(in): previous enthalpy for snow and soil
                       ! output
                       heatCapVegTrial,            & ! intent(out): volumetric heat capacity of vegetation canopy
                       mLayerHeatCapTrial,         & ! intent(out): heat capacity for snow and soil
                       ! output: error control
                       err,message)                  ! intent(out): error control
                       
   ! compute multiplier of state vector
   call computStatMult(&
                 ! input
                 heatCapVegTrial,                  & ! intent(in): volumetric heat capacity of vegetation canopy
                 mLayerHeatCapTrial,               & ! intent(in): volumetric heat capacity of soil and snow
                 diag_data,                        & ! intent(in): model diagnostic variables for a local HRU
                 indx_data,                        & ! intent(in): indices defining model states and layers
                 ! output
                 sMul,                             & ! intent(out): multiplier for state vector (used in the residual calculations)
                 err,cmessage)                       ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
 
   diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1) = heatCapVegTrial
   diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat(:) = mLayerHeatCapTrial(:)
   
   
   ! update thermal conductivity
   call computThermConduct(&
                       ! input: control variables
                       computeVegFlux,               & ! intent(in): flag to denote if computing the vegetation flux
                       canopyDepth,                  & ! intent(in): canopy depth (m)
                       ! input: state variables
                       scalarCanopyIceTrial,         & ! intent(in): trial value for canopy ice content (kg m-2)
                       scalarCanopyLiqTrial,         & ! intent(in): trial value for the liquid water on the vegetation canopy (kg m-2)
                       mLayerVolFracIceTrial,        & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                       mLayerVolFracLiqTrial,        & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
                        ! input/output: data structures
                       mpar_data,                    & ! intent(in): model parameters
                       indx_data,                    & ! intent(in): model layer indices
                       prog_data,                    & ! intent(in): model prognostic variables for a local HRU
                       diag_data,                    & ! intent(inout): model diagnostic variables for a local HRU
                       err,message)               ! intent(out): error control
   if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if
                       
  endif  ! if computing enthalpy

 ! save the number of flux calls per time step
 indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) = indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) + 1

 ! compute the fluxes for a given state vector
 call computFlux(&
                 ! input-output: model control
                 nSnow,                     & ! intent(in):    number of snow layers
                 nSoil,                     & ! intent(in):    number of soil layers
                 nLayers,                   & ! intent(in):    total number of layers
                 firstSubStep,              & ! intent(in):    flag to indicate if we are processing the first sub-step
                 firstFluxCall,             & ! intent(inout): flag to denote the first flux call
                 firstSplitOper,            & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                 computeVegFlux,            & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                 scalarSolution,            & ! intent(in):    flag to indicate the scalar solution
                 scalarSfcMeltPond/dt,      & ! intent(in):    drainage from the surface melt pond (kg m-2 s-1)
                 ! input: state variables
                 scalarCanairTempTrial,     & ! intent(in):    trial value for the temperature of the canopy air space (K)
                 scalarCanopyTempTrial,     & ! intent(in):    trial value for the temperature of the vegetation canopy (K)
                 mLayerTempTrial,           & ! intent(in):    trial value for the temperature of each snow and soil layer (K)
                 mLayerMatricHeadLiqTrial,  & ! intent(in):    trial value for the liquid water matric potential in each soil layer (m)
                 scalarAquiferStorageTrial, & ! intent(in):    trial value of storage of water in the aquifer (m)
                 ! input: diagnostic variables defining the liquid water and ice content
                 scalarCanopyLiqTrial,      & ! intent(in):    trial value for the liquid water on the vegetation canopy (kg m-2)
                 scalarCanopyIceTrial,      & ! intent(in):    trial value for the ice on the vegetation canopy (kg m-2)
                 mLayerVolFracLiqTrial,     & ! intent(in):    trial value for the volumetric liquid water content in each snow and soil layer (-)
                 mLayerVolFracIceTrial,     & ! intent(in):    trial value for the volumetric ice in each snow and soil layer (-)
                 ! input: data structures
                 model_decisions,           & ! intent(in):    model decisions
                 type_data,                 & ! intent(in):    type of vegetation and soil
                 attr_data,                 & ! intent(in):    spatial attributes
                 mpar_data,                 & ! intent(in):    model parameters
                 forc_data,                 & ! intent(in):    model forcing data
                 bvar_data,                 & ! intent(in):    average model variables for the entire basin
                 prog_data,                 & ! intent(in):    model prognostic variables for a local HRU
                 indx_data,                 & ! intent(in):    index data
                 ! input-output: data structures
                 diag_data,                 & ! intent(inout): model diagnostic variables for a local HRU
                 flux_data,                 & ! intent(inout): model fluxes for a local HRU
                 deriv_data,                & ! intent(out):   derivatives in model fluxes w.r.t. relevant state variables
                 ! input-output: flux vector and baseflow derivatives
                 ixSaturation,              & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                 dBaseflow_dMatric,         & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                 fluxVec,                   & ! intent(out):   flux vector (mixed units)
                 ! output: error control
                 err,cmessage)                ! intent(out):   error code and error message
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

 ! compute soil compressibility (-) and its derivative w.r.t. matric head (m)
 ! NOTE: we already extracted trial matrix head and volumetric liquid water as part of the flux calculations
 call soilCmpres(&
                 ! input:
                 ixRichards,                             & ! intent(in): choice of option for Richards' equation
                 ixBeg,ixEnd,                            & ! intent(in): start and end indices defining desired layers
                 mLayerMatricHeadLiq(1:nSoil),           & ! intent(in): matric head at the start of the time step (m)
                 mLayerMatricHeadLiqTrial(1:nSoil),      & ! intent(in): trial value of matric head (m)
                 mLayerVolFracLiqTrial(nSnow+1:nLayers), & ! intent(in): trial value for the volumetric liquid water content in each soil layer (-)
                 mLayerVolFracIceTrial(nSnow+1:nLayers), & ! intent(in): trial value for the volumetric ice content in each soil layer (-)
                 specificStorage,                        & ! intent(in): specific storage coefficient (m-1)
                 theta_sat,                              & ! intent(in): soil porosity (-)
                 ! output:
                 mLayerCompress,                         & ! intent(inout): compressibility of the soil matrix (-)
                 dCompress_dPsi,                         & ! intent(inout): derivative in compressibility w.r.t. matric head (m-1)
                 err,cmessage)                             ! intent(out): error code and error message
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

 ! compute the total change in storage associated with compression of the soil matrix (kg m-2)
 scalarSoilCompress = sum(mLayerCompress(1:nSoil)*mLayerDepth(nSnow+1:nLayers))*iden_water

 ! vegetation domain: get the correct water states (total water, or liquid water, depending on the state type)
 if(computeVegFlux)then
  scalarCanopyHydTrial = merge(scalarCanopyWatTrial, scalarCanopyLiqTrial, (ixStateType( ixHydCanopy(ixVegVolume) )==iname_watCanopy) )
 else
  scalarCanopyHydTrial = realMissing
 endif

 ! snow+soil domain: get the correct water states (total water, or liquid water, depending on the state type)
 mLayerVolFracHydTrial = merge(mLayerVolFracWatTrial, mLayerVolFracLiqTrial, (ixHydType==iname_watLayer .or. ixHydType==iname_matLayer) )
 
 
 ! *******************************************************************************************************
 ! ******************************************************************************************************* 
 ! ******************************************************************************************************* 
 
 ! compute the residual vector
 call computResid(&
                  ! input: model control
                  dt,                        & ! intent(in):    length of the time step (seconds)
                  nSnow,                     & ! intent(in):    number of snow layers
                  nSoil,                     & ! intent(in):    number of soil layers
                  nLayers,                   & ! intent(in):    total number of layers
                  ! input: flux vectors
                  sMul,                      & ! intent(in):    state vector multiplier (used in the residual calculations)
                  fluxVec,                   & ! intent(in):    flux vector
                  ! input: state variables (already disaggregated into scalars and vectors)
                  scalarCanairTempTrial,     & ! intent(in):    trial value for the temperature of the canopy air space (K)
                  scalarCanopyTempTrial,     & ! intent(in):    trial value for the temperature of the vegetation canopy (K)
                  scalarCanopyHydTrial,      & ! intent(in):    trial value of canopy hydrology state variable (kg m-2)
                  mLayerTempTrial,           & ! intent(in):    trial value for the temperature of each snow and soil layer (K)
                  mLayerVolFracHydTrial,     & ! intent(in):    trial vector of volumetric water content (-)
                  scalarAquiferStorageTrial, & ! intent(in):    trial value of storage of water in the aquifer (m)
                  ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
                  scalarCanopyIceTrial,      & ! intent(in):    trial value for the ice on the vegetation canopy (kg m-2)
                  mLayerVolFracIceTrial,     & ! intent(in):    trial value for the volumetric ice in each snow and soil layer (-)
                  ! input: data structures
                  prog_data,                 & ! intent(in):    model prognostic variables for a local HRU
                  diag_data,                 & ! intent(in):    model diagnostic variables for a local HRU
                  flux_data,                 & ! intent(in):    model fluxes for a local HRU
                  indx_data,                 & ! intent(in):    index data
                  ! output
                  resSink,                   & ! intent(out):   additional (sink) terms on the RHS of the state equation
                  resVec,                    & ! intent(out):   residual vector
                  err,cmessage)                ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

 ! compute the function evaluation
 rVecScaled = fScale(:)*real(resVec(:), rkind)   ! scale the residual vector (NOTE: residual vector is in quadruple precision)
 fEval      = 0.5_rkind*dot_product(rVecScaled,rVecScaled)

 ! end association with the information in the data structures
 end associate

 end subroutine eval8summa
end module eval8summa_module
