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

module opSplittin_module

! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing   ! missing integer
USE globalData,only:realMissing      ! missing double precision number
USE globalData,only:quadMissing      ! missing quadruple precision number

! access matrix information
USE globalData,only: nBands          ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix    ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix    ! named variable for the band diagonal matrix
USE globalData,only: iJac1           ! first layer of the Jacobian to print
USE globalData,only: iJac2           ! last layer of the Jacobian to print

! domain types
USE globalData,only:iname_cas        ! named variables for the canopy air space
USE globalData,only:iname_veg        ! named variables for vegetation
USE globalData,only:iname_snow       ! named variables for snow
USE globalData,only:iname_soil       ! named variables for soil

! state variable type
USE globalData,only:iname_nrgCanair  ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy  ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy  ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy  ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer   ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer   ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer   ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer   ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer   ! named variable defining the liquid matric potential state variable for soil layers
USE globalData,only:iname_watAquifer ! named variable defining the water storage in the aquifer

! global metadata
USE globalData,only:flux_meta                        ! metadata on the model fluxes
USE globalData,only:diag_meta                        ! metadata on the model diagnostic variables
USE globalData,only:prog_meta                        ! metadata on the model prognostic variables
USE globalData,only:deriv_meta                       ! metadata on the model derivatives
USE globalData,only:flux2state_orig                  ! metadata on flux-to-state mapping (original state variables)
USE globalData,only:flux2state_liq                   ! metadata on flux-to-state mapping (liquid water state variables)

! constants
USE multiconst,only:&
                    gravity,       & ! acceleration of gravity              (m s-2)
                    Tfreeze,       & ! temperature at freezing              (K)
                    LH_fus,        & ! latent heat of fusion                (J kg-1)
                    LH_vap,        & ! latent heat of vaporization          (J kg-1)
                    LH_sub,        & ! latent heat of sublimation           (J kg-1)
                    Cp_air,        & ! specific heat of air                 (J kg-1 K-1)
                    iden_air,      & ! intrinsic density of air             (kg m-3)
                    iden_ice,      & ! intrinsic density of ice             (kg m-3)
                    iden_water       ! intrinsic density of liquid water    (kg m-3)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookATTR        ! named variables for structure elements
USE var_lookup,only:iLookTYPE        ! named variables for structure elements
USE var_lookup,only:iLookPROG        ! named variables for structure elements
USE var_lookup,only:iLookDIAG        ! named variables for structure elements
USE var_lookup,only:iLookFLUX        ! named variables for structure elements
USE var_lookup,only:iLookFORCE       ! named variables for structure elements
USE var_lookup,only:iLookPARAM       ! named variables for structure elements
USE var_lookup,only:iLookINDEX       ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS   ! named variables for elements of the decision structure

! look up structure for variable types
USE var_lookup,only:iLookVarType

! provide access to the number of flux variables
USE var_lookup,only:nFlux=>maxvarFlux ! number of model flux variables

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_flagVec,  & ! data vector with variable length dimension (i4b)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    zLookup,      & ! data vector with variable length dimension (dp)
                    model_options   ! defines the model decisions

! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:       &
 localColumn,                     & ! separate groundwater representation in each local soil column
 singleBasin                        ! single groundwater store over the entire basin

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:      &
 qbaseTopmodel,                  & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                      & ! a big bucket (lumped aquifer model)
 noExplicit                        ! no explicit groundwater parameterization

! safety: set private unless specified otherwise
implicit none
private
public::opSplittin

! named variables for the coupling method
integer(i4b),parameter  :: fullyCoupled=1             ! 1st try: fully coupled solution
integer(i4b),parameter  :: stateTypeSplit=2           ! 2nd try: separate solutions for each state type

! named variables for the state variable split
integer(i4b),parameter  :: nrgSplit=1                 ! order in sequence for the energy operation
integer(i4b),parameter  :: massSplit=2                ! order in sequence for the mass operation

! named variables for the domain type split
integer(i4b),parameter  :: vegSplit=1                 ! order in sequence for the vegetation split
integer(i4b),parameter  :: snowSplit=2                ! order in sequence for the snow split
integer(i4b),parameter  :: soilSplit=3                ! order in sequence for the soil split
integer(i4b),parameter  :: aquiferSplit=4             ! order in sequence for the aquifer split

! named variables for the solution method
integer(i4b),parameter  :: vector=1                   ! vector solution method
integer(i4b),parameter  :: scalar=2                   ! scalar solution method
integer(i4b),parameter  :: nSolutions=2               ! number of solution methods

! named variables for the switch between states and domains
integer(i4b),parameter  :: fullDomain=1               ! full domain (veg+snow+soil)
integer(i4b),parameter  :: subDomain=2                ! sub domain (veg, snow, and soil separately)

! maximum number of possible splits
integer(i4b),parameter  :: nStateTypes=2              ! number of state types (energy, water)
integer(i4b),parameter  :: nDomains=4                 ! number of domains (vegetation, snow, soil, and aquifer)

! control parameters
real(rkind),parameter      :: valueMissing=-9999._rkind     ! missing value
real(rkind),parameter      :: verySmall=1.e-12_rkind        ! a very small number (used to check consistency)
real(rkind),parameter      :: veryBig=1.e+20_rkind          ! a very big number
real(rkind),parameter      :: dx = 1.e-8_rkind              ! finite difference increment

contains


 ! **********************************************************************************************************
 ! public subroutine opSplittin: run the coupled energy-mass model for one timestep
 !
 ! The logic of the solver is as follows:
 ! (1) Attempt different solutions in the following order: (a) fully coupled; (b) split by state type and by
 !      domain type for a given energy and mass split (vegetation, snow, and soil); and (c) scalar solution
 !      for a given state type and domain subset.
 ! (2) For a given split, compute a variable number of substeps (in varSubstepFida).
 ! **********************************************************************************************************
 subroutine opSplittin(&
                       ! input: model control
                       nSnow,          & ! intent(in):    number of snow layers
                       nSoil,          & ! intent(in):    number of soil layers
                       nLayers,        & ! intent(in):    total number of layers
                       nState,         & ! intent(in):    total number of state variables
                       dt,             & ! intent(inout): time step (s)
                       firstSubStep,   & ! intent(in):    flag to denote first sub-step
                       computeVegFlux, & ! intent(in):    flag to denote if computing energy flux over vegetation
                       ! input/output: data structures
                       type_data,      & ! intent(in):    type of vegetation and soil
                       attr_data,      & ! intent(in):    spatial attributes
                       forc_data,      & ! intent(in):    model forcing data
                       mpar_data,      & ! intent(in):    model parameters
                       indx_data,      & ! intent(inout): index data
                       prog_data,      & ! intent(inout): model prognostic variables for a local HRU
                       diag_data,      & ! intent(inout): model diagnostic variables for a local HRU
                       flux_data,      & ! intent(inout): model fluxes for a local HRU
                       bvar_data,      & ! intent(in):    model variables for the local basin
                       lookup_data,    & ! intent(in):    lookup tables
                       model_decisions,& ! intent(in):    model decisions
                       ! output: model control
                       dtMultiplier,   & ! intent(out):   substep multiplier (-)
                       tooMuchMelt,    & ! intent(out):   flag to denote that ice is insufficient to support melt
                       stepFailure,    & ! intent(out):   flag to denote step failure
                       ixCoupling,     & ! intent(out):   coupling method used in this iteration
                       err,message)      ! intent(out):   error code and error message
 ! ---------------------------------------------------------------------------------------
 ! structure allocations
 USE allocspace_module,only:allocLocal                ! allocate local data structures
 ! simulation of fluxes and residuals given a trial state vector
 USE soil_utils_module,only:matricHead                ! compute the matric head based on volumetric water content
 USE soil_utils_module,only:liquidHead                ! compute the liquid water matric potential
 ! population/extraction of state vectors
 USE indexState_module,only:indexSplit                ! get state indices
 USE varSubstep_module,only:varSubstep                ! complete substeps for a given split
 ! identify name of variable type (for error message)
 USE get_ixName_module,only:get_varTypeName           ! to access type strings for error messages
 implicit none
 ! ---------------------------------------------------------------------------------------
 ! * dummy variables
 ! ---------------------------------------------------------------------------------------
 ! input: model control
 integer(i4b),intent(in)         :: nSnow                          ! number of snow layers
 integer(i4b),intent(in)         :: nSoil                          ! number of soil layers
 integer(i4b),intent(in)         :: nLayers                        ! total number of layers
 integer(i4b),intent(in)         :: nState                         ! total number of state variables
 real(rkind),intent(inout)          :: dt                             ! time step (seconds)
 logical(lgt),intent(in)         :: firstSubStep                   ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)         :: computeVegFlux                 ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input/output: data structures
 type(var_i),intent(in)          :: type_data                      ! type of vegetation and soil
 type(var_d),intent(in)          :: attr_data                      ! spatial attributes
 type(var_d),intent(in)          :: forc_data                      ! model forcing data
 type(var_dlength),intent(in)    :: mpar_data                      ! model parameters
 type(var_ilength),intent(inout) :: indx_data                      ! indices for a local HRU
 type(var_dlength),intent(inout) :: prog_data                      ! prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data                      ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data                      ! model fluxes for a local HRU
 type(var_dlength),intent(in)    :: bvar_data                      ! model variables for the local basin
 type(zLookup),    intent(in)    :: lookup_data                    ! lookup tables
 type(model_options),intent(in)  :: model_decisions(:)             ! model decisions
 ! output: model control
 real(rkind),intent(out)            :: dtMultiplier                   ! substep multiplier (-)
 logical(lgt),intent(out)        :: tooMuchMelt                    ! flag to denote that ice is insufficient to support melt
 logical(lgt),intent(out)        :: stepFailure                    ! flag to denote step failure
 integer(i4b),intent(out)        :: err                            ! error code
 character(*),intent(out)        :: message                        ! error message
 ! *********************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************
 ! ---------------------------------------------------------------------------------------
 ! * general local variables
 ! ---------------------------------------------------------------------------------------
 character(LEN=256)              :: cmessage                       ! error message of downwind routine
 integer(i4b)                    :: minLayer                       ! the minimum layer used in assigning flags for flux aggregations
 integer(i4b)                    :: iOffset                        ! offset to account for different indices in the soil domain
 integer(i4b)                    :: iMin(1),iMax(1)                ! bounds of a given vector
 integer(i4b)                    :: iLayer,jLayer                  ! index of model layer
 integer(i4b)                    :: iSoil                          ! index of soil layer
 integer(i4b)                    :: iVar                           ! index of variables in data structures
 logical(lgt)                    :: firstSuccess                   ! flag to define the first success
 logical(lgt)                    :: firstFluxCall                  ! flag to define the first flux call
 logical(lgt)                    :: reduceCoupledStep              ! flag to define the need to reduce the length of the coupled step
 type(var_dlength)               :: prog_temp                      ! temporary model prognostic variables
 type(var_dlength)               :: diag_temp                      ! temporary model diagnostic variables
 type(var_dlength)               :: flux_temp                      ! temporary model fluxes
 type(var_dlength)               :: deriv_data                     ! derivatives in model fluxes w.r.t. relevant state variables
 real(rkind),dimension(nLayers)     :: mLayerVolFracIceInit           ! initial vector for volumetric fraction of ice (-)
 ! ------------------------------------------------------------------------------------------------------
 ! * operator splitting
 ! ------------------------------------------------------------------------------------------------------
 ! minimum timestep
 real(rkind),parameter              :: dtmin_coupled=1800._rkind         ! minimum time step for the fully coupled solution (seconds)
 real(rkind),parameter              :: dtmin_split=60._rkind             ! minimum time step for the fully split solution (seconds)
 real(rkind),parameter              :: dtmin_scalar=10._rkind            ! minimum time step for the scalar solution (seconds)
 real(rkind)                        :: dt_min                         ! minimum time step (seconds)
 real(rkind)                        :: dtInit                         ! initial time step (seconds)
 ! explicit error tolerance (depends on state type split, so defined here)
 real(rkind),parameter              :: errorTolLiqFlux=0.01_rkind        ! error tolerance in the explicit solution (liquid flux)
 real(rkind),parameter              :: errorTolNrgFlux=10._rkind         ! error tolerance in the explicit solution (energy flux)
 ! number of substeps taken for a given split
 integer(i4b)                    :: nSubsteps                      ! number of substeps taken for a given split
 ! named variables defining the coupling and solution method
 integer(i4b)                    :: ixCoupling                     ! index of coupling method (1,2)
 integer(i4b)                    :: ixSolution                     ! index of solution method (1,2)
 integer(i4b)                    :: ixStateThenDomain              ! switch between the state and domain (1,2)
 integer(i4b)                    :: tryDomainSplit                 ! (0,1) - flag to try the domain split
 ! actual number of splits
 integer(i4b)                    :: nStateTypeSplit                ! number of splits for the state type
 integer(i4b)                    :: nDomainSplit                   ! number of splits for the domain
 integer(i4b)                    :: nStateSplit                    ! number of splits for the states within a given domain
 ! indices for the state type and the domain split
 integer(i4b)                    :: iStateTypeSplit                ! index of the state type split
 integer(i4b)                    :: iDomainSplit                   ! index of the domain split
 integer(i4b)                    :: iStateSplit                    ! index of the state split
 ! flux masks
 logical(lgt)                    :: neededFlux(nFlux)              ! .true. if flux is needed at all
 logical(lgt)                    :: desiredFlux                    ! .true. if flux is desired for a given split
 type(var_ilength)               :: fluxCount                      ! number of times each flux is updated (should equal nSubsteps)
 type(var_flagVec)               :: fluxMask                       ! mask defining model fluxes
 ! state masks
 integer(i4b),dimension(nState)  :: stateCheck                     ! number of times each state variable is updated (should equal 1)
 logical(lgt),dimension(nState)  :: stateMask                      ! mask defining desired state variables
 integer(i4b)                    :: nSubset                        ! number of selected state variables for a given split
 ! flags
 logical(lgt)                    :: failure                        ! flag to denote failure of substepping
 logical(lgt)                    :: doAdjustTemp                   ! flag to adjust temperature after the mass split
 logical(lgt)                    :: failedMinimumStep              ! flag to denote failure of substepping for a given split
 integer(i4b)                    :: ixSaturation                   ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
 integer(i4b),parameter          :: IDA=1
 integer(i4b),parameter          :: BE=2
 integer(i4b)                    :: solver=BE   				   ! BE or IDA
 integer(i4b)                    :: nCoupling
 ! ---------------------------------------------------------------------------------------
 ! point to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 globalVars: associate(&
 ! model decisions
 ixGroundwater           => model_decisions(iLookDECISIONS%groundwatr)%iDecision   ,& ! intent(in):    [i4b]    groundwater parameterization
 ixSpatialGroundwater    => model_decisions(iLookDECISIONS%spatial_gw)%iDecision   ,& ! intent(in):    [i4b]    spatial representation of groundwater (local-column or single-basin)
 ! domain boundary conditions
 airtemp                 => forc_data%var(iLookFORCE%airtemp)                      ,& ! intent(in):    [dp]     temperature of the upper boundary of the snow and soil domains (K)
 ! vector of energy and hydrology indices for the snow and soil domains
 ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
 ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
 nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in):    [i4b]    number of energy state variables in the snow+soil domain
 nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in):    [i4b]    number of hydrology state variables in the snow+soil domain
 ! indices of model state variables
 ixMapSubset2Full        => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat         ,& ! intent(in):    [i4b(:)] list of indices in the state subset (missing for values not in the subset)
 ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,& ! intent(in):    [i4b(:)] indices defining the type of the state (ixNrgState...)
 ixNrgCanair             => indx_data%var(iLookINDEX%ixNrgCanair)%dat              ,& ! intent(in):    [i4b(:)] indices IN THE FULL VECTOR for energy states in canopy air space domain
 ixNrgCanopy             => indx_data%var(iLookINDEX%ixNrgCanopy)%dat              ,& ! intent(in):    [i4b(:)] indices IN THE FULL VECTOR for energy states in the canopy domain
 ixHydCanopy             => indx_data%var(iLookINDEX%ixHydCanopy)%dat              ,& ! intent(in):    [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the canopy domain
 ixNrgLayer              => indx_data%var(iLookINDEX%ixNrgLayer)%dat               ,& ! intent(in):    [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
 ixHydLayer              => indx_data%var(iLookINDEX%ixHydLayer)%dat               ,& ! intent(in):    [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
 ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy air space energy state variable
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy energy state variable
 ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy hydrology state variable (mass)
 ! numerix tracking
 numberStateSplit        => indx_data%var(iLookINDEX%numberStateSplit     )%dat(1) ,& ! intent(inout): [i4b]    number of state splitting solutions             (-)
 numberDomainSplitNrg    => indx_data%var(iLookINDEX%numberDomainSplitNrg )%dat(1) ,& ! intent(inout): [i4b]    number of domain splitting solutions for energy (-)
 numberDomainSplitMass   => indx_data%var(iLookINDEX%numberDomainSplitMass)%dat(1) ,& ! intent(inout): [i4b]    number of domain splitting solutions for mass   (-)
 numberScalarSolutions   => indx_data%var(iLookINDEX%numberScalarSolutions)%dat(1) ,& ! intent(inout): [i4b]    number of scalar solutions                      (-)
 ! domain configuration
 canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in):    [dp]     canopy depth (m)
 mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,& ! intent(in):    [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
 ! snow parameters
 snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)         ,& ! intent(in):    [dp]     scaling parameter for the snow freezing curve (K-1)
 ! depth-varying soil parameters
 vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat               ,& ! intent(in):    [dp(:)]  van Genutchen "m" parameter (-)
 vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)%dat                    ,& ! intent(in):    [dp(:)]  van Genutchen "n" parameter (-)
 vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)%dat                ,& ! intent(in):    [dp(:)]  van Genutchen "alpha" parameter (m-1)
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat                ,& ! intent(in):    [dp(:)]  soil porosity (-)
 theta_res               => mpar_data%var(iLookPARAM%theta_res)%dat                ,& ! intent(in):    [dp(:)]  soil residual volumetric water content (-)
 ! soil parameters
 specificStorage         => mpar_data%var(iLookPARAM%specificStorage)%dat(1)       ,& ! intent(in):    [dp]     specific storage coefficient (m-1)
 ! model diagnostic variables (fraction of liquid water)
 scalarFracLiqVeg        => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)       ,& ! intent(out):   [dp]     fraction of liquid water on vegetation (-)
 mLayerFracLiqSnow       => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat         ,& ! intent(out):   [dp(:)]  fraction of liquid water in each snow layer (-)
 mLayerMeltFreeze        => diag_data%var(iLookDIAG%mLayerMeltFreeze)%dat          ,& ! intent(out):   [dp(:)]  melt-freeze in each snow and soil layer (kg m-3)
 ! model state variables (vegetation canopy)
 scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(out):   [dp]     temperature of the canopy air space (K)
 scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(out):   [dp]     temperature of the vegetation canopy (K)
 scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(out):   [dp]     mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,& ! intent(out):   [dp]     mass of liquid water on the vegetation canopy (kg m-2)
 scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(out):   [dp]     mass of total water on the vegetation canopy (kg m-2)
 ! model state variables (snow and soil domains)
 mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(out):   [dp(:)]  temperature of each snow/soil layer (K)
 mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,& ! intent(out):   [dp(:)]  volumetric fraction of ice (-)
 mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(out):   [dp(:)]  volumetric fraction of liquid water (-)
 mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(out):   [dp(:)]  volumetric fraction of total water (-)
 mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(out):   [dp(:)]  matric head (m)
 mLayerMatricHeadLiq     => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat        & ! intent(out):   [dp(:)]  matric potential of liquid water (m)
 )
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="opSplittin/"
 
 ! we just solve the fully coupled problem by ida
 select case(solver)
 	case(BE); nCoupling = 2
    case(IDA); nCoupling = 1
    case default; err=20; message=trim(message)//'expect case to be IDA or BE'; return  
 end select     

 ! *****
 ! (0) PRELIMINARIES...
 ! ********************

 ! -----
 ! * initialize...
 ! ---------------

 ! set the global print flag
 globalPrintFlag=.false.

 if(globalPrintFlag)&
 print*, trim(message), dt

 ! initialize the first success call
 firstSuccess=.false.

 ! initialize the flags
 tooMuchMelt=.false.  ! too much melt (merge snow layers)
 stepFailure=.false.  ! step failure

 ! initialize flag for the success of the substepping
 failure=.false.

 ! initialize the flux check
 neededFlux(:) = .false.

 ! initialize the state check
 stateCheck(:) = 0

 ! compute the total water content in the vegetation canopy
 scalarCanopyWat = scalarCanopyLiq + scalarCanopyIce  ! kg m-2

 ! save volumetric ice content at the start of the step
 ! NOTE: used for volumetric loss due to melt-freeze
 mLayerVolFracIceInit(:) = mLayerVolFracIce(:)

 ! compute the total water content in snow and soil
 ! NOTE: no ice expansion allowed for soil
 if(nSnow>0)&
 mLayerVolFracWat(      1:nSnow  ) = mLayerVolFracLiq(      1:nSnow  ) + mLayerVolFracIce(      1:nSnow  )*(iden_ice/iden_water)
 mLayerVolFracWat(nSnow+1:nLayers) = mLayerVolFracLiq(nSnow+1:nLayers) + mLayerVolFracIce(nSnow+1:nLayers)

 ! compute the liquid water matric potential (m)
 ! NOTE: include ice content as part of the solid porosity - major effect of ice is to reduce the pore size; ensure that effSat=1 at saturation
 ! (from Zhao et al., J. Hydrol., 1997: Numerical analysis of simultaneous heat and mass transfer...)
 do iSoil=1,nSoil
  call liquidHead(mLayerMatricHead(iSoil),mLayerVolFracLiq(nSnow+iSoil),mLayerVolFracIce(nSnow+iSoil),  & ! input:  state variables
                  vGn_alpha(iSoil),vGn_n(iSoil),theta_sat(iSoil),theta_res(iSoil),vGn_m(iSoil),         & ! input:  parameters
                  matricHeadLiq=mLayerMatricHeadLiq(iSoil),                                             & ! output: liquid water matric potential (m)
                  err=err,message=cmessage)                                                               ! output: error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do  ! looping through soil layers (computing liquid water matric potential)

 ! allocate space for the flux mask (used to define when fluxes are updated)
 call allocLocal(flux_meta(:),fluxMask,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for the flux count (used to check that fluxes are only updated once)
 call allocLocal(flux_meta(:),fluxCount,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for the temporary prognostic variable structure
 call allocLocal(prog_meta(:),prog_temp,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for the temporary diagnostic variable structure
 call allocLocal(diag_meta(:),diag_temp,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for the temporary flux variable structure
 call allocLocal(flux_meta(:),flux_temp,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for the derivative structure
 call allocLocal(deriv_meta(:),deriv_data,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! intialize the flux conter
 do iVar=1,size(flux_meta)  ! loop through fluxes
  fluxCount%var(iVar)%dat(:) = 0
 end do

 ! initialize the model fluxes
 do iVar=1,size(flux_meta)  ! loop through fluxes
  if(flux2state_orig(iVar)%state1==integerMissing .and. flux2state_orig(iVar)%state2==integerMissing) cycle ! flux does not depend on state (e.g., input)
  if(flux2state_orig(iVar)%state1==iname_watCanopy .and. .not.computeVegFlux) cycle ! use input fluxes in cases where there is no canopy
  flux_data%var(iVar)%dat(:) = 0._rkind
 end do

 ! initialize derivatives
 do iVar=1,size(deriv_meta)
  deriv_data%var(iVar)%dat(:) = 0._rkind
 end do

 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================

 ! loop through different coupling strategies
 coupling: do ixCoupling=1,nCoupling

  ! initialize the time step
  dtInit = min( merge(dt,            dtmin_coupled, ixCoupling==fullyCoupled), dt) ! initial time step
  dt_min = min( merge(dtmin_coupled, dtmin_split,   ixCoupling==fullyCoupled), dt) ! minimum time step

  ! keep track of the number of state splits
  if(ixCoupling/=fullyCoupled) numberStateSplit = numberStateSplit + 1

  ! define the number of operator splits for the state type
  select case(ixCoupling)
   case(fullyCoupled);   nStateTypeSplit=1
   case(stateTypeSplit); nStateTypeSplit=nStateTypes
   case default; err=20; message=trim(message)//'coupling case not found'; return
  end select  ! operator splitting option

  ! define if we wish to try the domain split
  select case(ixCoupling)
   case(fullyCoupled);   tryDomainSplit=0
   case(stateTypeSplit); tryDomainSplit=1
   case default; err=20; message=trim(message)//'coupling case not found'; return
  end select  ! operator splitting option

  ! state splitting loop
  stateTypeSplit: do iStateTypeSplit=1,nStateTypeSplit

   !print*, 'iStateTypeSplit, nStateTypeSplit = ', iStateTypeSplit, nStateTypeSplit

   ! -----
   ! * identify state-specific variables for a given state split...
   ! --------------------------------------------------------------

   ! flag to adjust the temperature
   doAdjustTemp = (ixCoupling/=fullyCoupled .and. iStateTypeSplit==massSplit)

   ! modify the state type names associated with the state vector
   if(ixCoupling/=fullyCoupled .and. iStateTypeSplit==massSplit)then
    if(computeVegFlux)then
     where(ixStateType(ixHydCanopy)==iname_watCanopy) ixStateType(ixHydCanopy)=iname_liqCanopy
    endif
    where(ixStateType(ixHydLayer) ==iname_watLayer)  ixStateType(ixHydLayer) =iname_liqLayer
    where(ixStateType(ixHydLayer) ==iname_matLayer)  ixStateType(ixHydLayer) =iname_lmpLayer
   endif  ! if modifying state variables for the mass split

   ! first try the state type split, then try the domain split within a given state type
   stateThenDomain: do ixStateThenDomain=1,1+tryDomainSplit ! 1=state type split; 2=domain split within a given state type

    !print*, 'start of stateThenDomain loop'

    ! keep track of the number of domain splits
    if(iStateTypeSplit==nrgSplit  .and. ixStateThenDomain==subDomain) numberDomainSplitNrg  = numberDomainSplitNrg  + 1
    if(iStateTypeSplit==massSplit .and. ixStateThenDomain==subDomain) numberDomainSplitMass = numberDomainSplitMass + 1

    ! define the number of domain splits for the state type
    select case(ixStateThenDomain)
     case(fullDomain); nDomainSplit=1
     case(subDomain);  nDomainSplit=nDomains
     case default; err=20; message=trim(message)//'coupling case not found'; return
    end select

    ! check that we haven't split the domain when we are fully coupled
    if(ixCoupling==fullyCoupled .and. nDomainSplit==nDomains)then
     message=trim(message)//'cannot split domains when fully coupled'
     err=20; return
    endif

    ! domain splitting loop
    domainSplit: do iDomainSplit=1,nDomainSplit

     ! trial with the vector then scalar solution
     solution: do ixSolution=1,nSolutions

      ! initialize error control
      err=0; message="opSplittin/"

      ! refine the time step
      if(ixSolution==scalar)then
       dtInit = min(dtmin_split, dt)    ! initial time step
       dt_min = min(dtmin_scalar, dt)   ! minimum time step
      endif

      ! initialize the first flux call
      firstFluxCall=.true.

      ! get the number of split layers
      select case(ixSolution)
       case(vector); nStateSplit=1
       case(scalar); nStateSplit=count(stateMask)
       case default; err=20; message=trim(message)//'unknown solution method'; return
      end select

      !print*, '*****'
      !print*, 'computeVegFlux = ', computeVegFlux
      !print*, '(ixSolution==scalar) = ', (ixSolution==scalar)
      !print*, 'ixCoupling, iStateTypeSplit, ixStateThenDomain, iDomainSplit, nDomainSplit: ', ixCoupling, iStateTypeSplit, ixStateThenDomain, iDomainSplit, nDomainSplit
      !print*, 'ixSoilOnlyHyd = ', indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat

      ! loop through layers (NOTE: nStateSplit=1 for the vector solution, hence no looping)
      stateSplit: do iStateSplit=1,nStateSplit

       ! -----
       ! * define state subsets for a given split...
       ! -------------------------------------------

       ! get the mask for the state subset
       call stateFilter(ixCoupling,ixSolution,ixStateThenDomain,iStateTypeSplit,iDomainSplit,iStateSplit,&
                        indx_data,stateMask,nSubset,err,cmessage)
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

       ! check that state variables exist
       if(nSubset==0) cycle domainSplit

       ! avoid redundant case where vector solution is of length 1
       if(ixSolution==vector .and. count(stateMask)==1) cycle solution

       ! check
       !print*, 'after stateFilter: stateMask   = ', stateMask
       !print*, 'count(stateMask) = ', count(stateMask)

       !if(ixSolution==scalar)then
       ! print*, 'iStateSplit, nStateSplit = ', iStateSplit, nStateSplit
       ! print*, 'start of scalar solution'
       ! !print*, 'PAUSE'; read(*,*)
       !endif

       ! -----
       ! * assemble vectors for a given split...
       ! ---------------------------------------

       ! get indices for a given split
       call indexSplit(stateMask,                   & ! intent(in)    : logical vector (.true. if state is in the subset)
                       nSnow,nSoil,nLayers,nSubset, & ! intent(in)    : number of snow and soil layers, and total number of layers
                       indx_data,                   & ! intent(inout) : index data structure
                       err,cmessage)                  ! intent(out)   : error control
       if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

       ! -----
       ! * define the mask of the fluxes used...
       ! ---------------------------------------

       ! identify the type of state for the states in the subset
       stateSubset: associate(ixStateType_subset  => indx_data%var(iLookINDEX%ixStateType_subset)%dat ,& ! intent(in): [i4b(:)] indices of state types
                              ixMapFull2Subset    => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat   ,& ! intent(in): [i4b(:)] mapping of full state vector to the state subset
                              ixControlVolume     => indx_data%var(iLookINDEX%ixControlVolume)%dat    ,& ! intent(in): [i4b(:)] index of control volume for different domains (veg, snow, soil)
                              ixLayerActive       => indx_data%var(iLookINDEX%ixLayerActive)%dat      ,& ! intent(in): [i4b(:)] list of indices for all active layers (inactive=integerM
                              ixDomainType        => indx_data%var(iLookINDEX%ixDomainType)%dat       )  ! intent(in): [i4b(:)] indices defining the type of the domain (iname_veg, iname_snow, iname_soil)

       ! loop through flux variables
       do iVar=1,size(flux_meta)

        ! * identify flux mask for the fully coupled solution
        if(ixCoupling==fullyCoupled)then
         desiredFlux            = any(ixStateType_subset==flux2state_orig(iVar)%state1) .or. any(ixStateType_subset==flux2state_orig(iVar)%state2)
         fluxMask%var(iVar)%dat = desiredFlux

        ! * identify flux mask for the split solution
        else

         ! identify the flux mask for a given state split
         select case(iStateTypeSplit)
          case(nrgSplit);  desiredFlux = any(ixStateType_subset==flux2state_orig(iVar)%state1) .or. any(ixStateType_subset==flux2state_orig(iVar)%state2)
          case(massSplit); desiredFlux = any(ixStateType_subset==flux2state_liq(iVar)%state1)  .or. any(ixStateType_subset==flux2state_liq(iVar)%state2)
          case default; err=20; message=trim(message)//'unable to identify split based on state type'; return
         end select

         ! no domain splitting
         if(nDomains==1)then
          fluxMask%var(iVar)%dat = desiredFlux

         ! domain splitting
         else

          ! initialize to .false.
          fluxMask%var(iVar)%dat = .false.

          ! only need to proceed if the flux is desired
          if(desiredFlux)then

           ! different domain splitting operations
           select case(iDomainSplit)

            ! canopy fluxes -- (:1) gets the upper boundary(0) if it exists
            case(vegSplit)

             ! vector solution (should only be present for energy)
             if(ixSolution==vector)then
              fluxMask%var(iVar)%dat(:1) = desiredFlux
              if(ixStateThenDomain>1 .and. iStateTypeSplit/=nrgSplit)then
               message=trim(message)//'only expect a vector solution for the vegetation domain for energy'
               err=20; return
              endif

             ! scalar solution
             else
              fluxMask%var(iVar)%dat(:1) = desiredFlux
             endif

            ! fluxes through snow and soil
            case(snowSplit,soilSplit)

             ! loop through layers
             do iLayer=1,nLayers
              if(ixlayerActive(iLayer)/=integerMissing)then

               ! get the offset (ixLayerActive=1,2,3,...nLayers, and soil vectors nSnow+1, nSnow+2, ..., nLayers)
               iOffset = merge(nSnow, 0, flux_meta(iVar)%vartype==iLookVarType%midSoil .or. flux_meta(iVar)%vartype==iLookVarType%ifcSoil)
               jLayer  = iLayer-iOffset

               ! identify the minimum layer
               select case(flux_meta(iVar)%vartype)
                case(iLookVarType%ifcToto, iLookVarType%ifcSnow, iLookVarType%ifcSoil); minLayer=merge(jLayer-1, jLayer, jLayer==1)
                case(iLookVarType%midToto, iLookVarType%midSnow, iLookVarType%midSoil); minLayer=jLayer
                case default; minLayer=integerMissing
               end select

               ! set desired layers
               select case(flux_meta(iVar)%vartype)
                case(iLookVarType%midToto,iLookVarType%ifcToto);                   fluxMask%var(iVar)%dat(minLayer:jLayer) = desiredFlux
                case(iLookVarType%midSnow,iLookVarType%ifcSnow); if(iLayer<=nSnow) fluxMask%var(iVar)%dat(minLayer:jLayer) = desiredFlux
                case(iLookVarType%midSoil,iLookVarType%ifcSoil); if(iLayer> nSnow) fluxMask%var(iVar)%dat(minLayer:jLayer) = desiredFlux
               end select

               ! add hydrology states for scalar variables
               if(iStateTypeSplit==massSplit .and. flux_meta(iVar)%vartype==iLookVarType%scalarv)then
                select case(iDomainSplit)
                 case(snowSplit); if(iLayer==nSnow)   fluxMask%var(iVar)%dat = desiredFlux
                 case(soilSplit); if(iLayer==nSnow+1) fluxMask%var(iVar)%dat = desiredFlux
                end select
               endif  ! if hydrology split and scalar

              endif    ! if the layer is active
             end do   ! looping through layers

            ! check
            case default; err=20; message=trim(message)//'unable to identify split based on domain type'; return
           end select  ! domain split

          endif  ! if flux is desired

         endif  ! domain splitting
        endif  ! not fully coupled

        ! define if the flux is desired
        if(desiredFlux) neededFlux(iVar)=.true.
        !if(desiredFlux) print*, flux_meta(iVar)%varname, fluxMask%var(iVar)%dat

        ! * check
        if( globalPrintFlag .and. count(fluxMask%var(iVar)%dat)>0 )&
        print*, trim(flux_meta(iVar)%varname)

       end do  ! (loop through fluxes)

       end associate stateSubset

       ! *******************************************************************************************************************************
       ! *******************************************************************************************************************************
       ! *******************************************************************************************************************************
       ! ***** trial with a given solution method...

       ! check that we do not attempt the scalar solution for the fully coupled case
       if(ixCoupling==fullyCoupled .and. ixSolution==scalar)then
        message=trim(message)//'only apply the scalar solution to the fully split coupling strategy'
        err=20; return
       endif

       ! reset the flag for the first flux call
       if(.not.firstSuccess) firstFluxCall=.true.

       ! save/recover copies of prognostic variables
       do iVar=1,size(prog_data%var)
        select case(failure)
         case(.false.); prog_temp%var(iVar)%dat(:) = prog_data%var(iVar)%dat(:)
         case(.true.);  prog_data%var(iVar)%dat(:) = prog_temp%var(iVar)%dat(:)
        end select
       end do  ! looping through variables

       ! save/recover copies of diagnostic variables
       do iVar=1,size(diag_data%var)
        select case(failure)
         case(.false.); diag_temp%var(iVar)%dat(:) = diag_data%var(iVar)%dat(:)
         case(.true.);  diag_data%var(iVar)%dat(:) = diag_temp%var(iVar)%dat(:)
        end select
       end do  ! looping through variables

       ! save/recover copies of model fluxes
       do iVar=1,size(flux_data%var)
        select case(failure)
         case(.false.); flux_temp%var(iVar)%dat(:) = flux_data%var(iVar)%dat(:)
         case(.true.);  flux_data%var(iVar)%dat(:) = flux_temp%var(iVar)%dat(:)
        end select
       end do  ! looping through variables

       ! -----
       ! * solve variable subset for one time step...
       ! --------------------------------------------

       !print*, trim(message)//'before varSubstepFida: nSubset = ', nSubset

       ! keep track of the number of scalar solutions
       if(ixSolution==scalar) numberScalarSolutions = numberScalarSolutions + 1
       
       ! solve variable subset for one full time step
       select case(solver)
        case(IDA)
			 print *, 'IDA Solver not implemented yet'
			 stop 1
        case(BE) 
             call varSubstep(&
                       ! input: model control
                       dt,                         & ! intent(inout) : time step (s)
                       dtInit,                     & ! intent(in)    : initial time step (seconds)
                       dt_min,                     & ! intent(in)    : minimum time step (seconds)
                       nSubset,                    & ! intent(in)    : total number of variables in the state subset
                       doAdjustTemp,               & ! intent(in)    : flag to indicate if we adjust the temperature
                       firstSubStep,               & ! intent(in)    : flag to denote first sub-step
                       firstFluxCall,              & ! intent(inout) : flag to indicate if we are processing the first flux call
                       computeVegFlux,             & ! intent(in)    : flag to denote if computing energy flux over vegetation
                       (ixSolution==scalar),       & ! intent(in)    : flag to denote computing the scalar solution
                       iStateSplit,                & ! intent(in)    : index of the layer in the splitting operation
                       fluxMask,                   & ! intent(in)    : mask for the fluxes used in this given state subset
                       fluxCount,                  & ! intent(inout) : number of times fluxes are updated (should equal nsubstep)
                       ! input/output: data structures
                       model_decisions,            & ! intent(in)    : model decisions
                       lookup_data,                & ! intent(in)    : lookup tables
                       type_data,                  & ! intent(in)    : type of vegetation and soil
                       attr_data,                  & ! intent(in)    : spatial attributes
                       forc_data,                  & ! intent(in)    : model forcing data
                       mpar_data,                  & ! intent(in)    : model parameters
                       indx_data,                  & ! intent(inout) : index data
                       prog_data,                  & ! intent(inout) : model prognostic variables for a local HRU
                       diag_data,                  & ! intent(inout) : model diagnostic variables for a local HRU
                       flux_data,                  & ! intent(inout) : model fluxes for a local HRU
                       deriv_data,                 & ! intent(inout) : derivatives in model fluxes w.r.t. relevant state variables
                       bvar_data,                  & ! intent(in)    : model variables for the local basin
                       ! output: control
                       ixSaturation,               & ! intent(inout) : index of the lowest saturated layer (NOTE: only computed on the first iteration)
                       dtMultiplier,               & ! intent(out)   : substep multiplier (-)
                       nSubsteps,                  & ! intent(out)   : number of substeps taken for a given split
                       failedMinimumStep,          & ! intent(out)   : flag for failed substeps
                       reduceCoupledStep,          & ! intent(out)   : flag to reduce the length of the coupled step
                       tooMuchMelt,                & ! intent(out)   : flag to denote that ice is insufficient to support melt
                       err,cmessage)                 ! intent(out)   : error code and error message 
                     ! check
          case default; err=20; message=trim(message)//'expect case to be ida or be'; return  
        end select  
        
         
             
       if(err/=0)then
        message=trim(message)//trim(cmessage)
        if(err>0) return
       endif  ! (check for errors)

       !print*, trim(message)//'after varSubstepFida: scalarSnowDrainage = ', flux_data%var(iLookFLUX%scalarSnowDrainage)%dat
       !print*, trim(message)//'after varSubstepFida: iLayerLiqFluxSnow  = ', flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat
       !print*, trim(message)//'after varSubstepFida: iLayerLiqFluxSoil  = ', flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat

       ! check
       !if(ixSolution==scalar)then
       ! print*, 'PAUSE: check scalar'; read(*,*)
       !endif

       ! reduce coupled step if failed the minimum step for the scalar solution
       if(failedMinimumStep .and. ixSolution==scalar) reduceCoupledStep=.true.

       ! check
       !if(ixCoupling/=fullyCoupled)then
       ! print*, 'dt = ', dt
       ! print*, 'after varSubstepFida: err              = ', err
       ! print*, 'after varSubstepFida: cmessage         = ', trim(cmessage)
       ! print*, 'after varSubstepFida: computeVegFlux   = ', computeVegFlux
       ! print*, 'after varSubstepFida: stateMask        = ', stateMask
       ! print*, 'after varSubstepFida: coupling         = ', (ixCoupling==fullyCoupled)
       ! print*, 'after varSubstepFida: scalar solve     = ', (ixSolution==scalar)
       ! print*, 'iStateTypeSplit, nStateTypeSplit = ', iStateTypeSplit, nStateTypeSplit
       ! print*, 'iDomainSplit,    nDomainSplit    = ', iDomainSplit,    nDomainSplit
       ! print*, 'nSubset           = ', nSubset
       ! print*, 'tooMuchMelt       = ', tooMuchMelt
       ! print*, 'reduceCoupledStep = ', reduceCoupledStep
       ! print*, 'failedMinimumStep = ', failedMinimumStep, merge('coupled','opSplit',ixCoupling==fullyCoupled)
       ! if(ixSolution==scalar)then; print*, 'PAUSE'; read(*,*); endif
       !endif

       !if(ixSolution==scalar)then
       ! !print*, trim(message)//'stop: checking scalar solution'; stop
       ! print*, trim(message)//'pause: checking scalar solution'; read(*,*)
       !endif

       !print*, 'tooMuchMelt, reduceCoupledStep = ', tooMuchMelt, reduceCoupledStep

       ! if too much melt (or some other need to reduce the coupled step) then return
       ! NOTE: need to go all the way back to coupled_em and merge snow layers, as all splitting operations need to occur with the same layer geometry
       if(tooMuchMelt .or. reduceCoupledStep)then
        stepFailure=.true.
        err=0 ! recovering
        return
       endif

       ! define failure
       failure = (failedMinimumStep .or. err<0)
       if(.not.failure) firstSuccess=.true.

       ! if failed, need to reset the flux counter
       if(failure)then
        !print*, 'failure!'
        do iVar=1,size(flux_meta)
         iMin=lbound(flux_data%var(iVar)%dat)
         iMax=ubound(flux_data%var(iVar)%dat)
         do iLayer=iMin(1),iMax(1)
          if(fluxMask%var(iVar)%dat(iLayer)) fluxCount%var(iVar)%dat(iLayer) = fluxCount%var(iVar)%dat(iLayer) - nSubsteps
         end do
         !if(iVar==iLookFLUX%mLayerTranspire) print*, flux_meta(iVar)%varname, fluxCount%var(iVar)%dat
        end do
       endif

       ! try the fully split solution if failed to converge with a minimum time step in the coupled solution
       if(ixCoupling==fullyCoupled .and. failure) cycle coupling

       ! try the scalar solution if failed to converge with a minimum time step in the split solution
       if(ixCoupling/=fullyCoupled)then
        select case(ixStateThenDomain)
         case(fullDomain); if(failure) cycle stateThenDomain
         case(subDomain);  if(failure) cycle solution
         case default; err=20; message=trim(message)//'unknown ixStateThenDomain case'
        end select
       endif

       ! check that state variables updated
       where(stateMask) stateCheck = stateCheck+1
       if(any(stateCheck>1))then
        message=trim(message)//'state variable updated more than once!'
        err=20; return
       endif

       ! success = exit solution
       if(.not.failure)then
        select case(ixStateThenDomain)
         case(fullDomain); if(iStateSplit==nStateSplit) exit stateThenDomain
         case(subDomain);  if(iStateSplit==nStateSplit) exit solution
         case default; err=20; message=trim(message)//'unknown ixStateThenDomain case'
        end select
       else

        ! check that we did not fail for the scalar solution (last resort)
        if(ixSolution==scalar)then
         message=trim(message)//'failed the minimum step for the scalar solution'
         err=20; return

        ! check for an unexpected failure
        else
         message=trim(message)//'unexpected failure'
         err=20; return
        endif

       endif  ! success check

      end do stateSplit ! solution with split layers
      !print*, 'after stateSplit'

     end do solution ! trial with the full layer solution then the split layer solution

     !print*, 'after solution loop'

     ! ***** trial with a given solution method...
     ! *******************************************************************************************************************************
     ! *******************************************************************************************************************************
     ! *******************************************************************************************************************************

    end do domainSplit ! domain type splitting loop

    !print*, 'ixStateThenDomain = ', ixStateThenDomain
    !print*, 'after domain split loop'

   end do stateThenDomain  ! switch between the state and the domain

   !print*, 'after stateThenDomain switch'

   ! -----
   ! * reset state variables for the mass split...
   ! ---------------------------------------------

   ! modify the state type names associated with the state vector
   if(ixCoupling/=fullyCoupled .and. iStateTypeSplit==massSplit)then
    if(computeVegFlux)then
     where(ixStateType(ixHydCanopy)==iname_liqCanopy) ixStateType(ixHydCanopy)=iname_watCanopy
    endif
    where(ixStateType(ixHydLayer) ==iname_liqLayer)  ixStateType(ixHydLayer) =iname_watLayer
    where(ixStateType(ixHydLayer) ==iname_lmpLayer)  ixStateType(ixHydLayer) =iname_matLayer
   endif  ! if modifying state variables for the mass split

  end do stateTypeSplit ! state type splitting loop

  ! check
  !if(ixCoupling/=fullyCoupled)then
  ! print*, 'PAUSE: end of splitting loop'; read(*,*)
  !endif

  ! ==========================================================================================================================================
  ! ==========================================================================================================================================

  ! success = exit the coupling loop
  if(ixCoupling==fullyCoupled .and. .not.failure) exit coupling

 end do coupling ! coupling method

 ! check that all state variables were updated
 if(any(stateCheck==0))then
  message=trim(message)//'some state variables were not updated!'
  err=20; return
 endif

 ! check that the desired fluxes were computed
 do iVar=1,size(flux_meta)
  if(neededFlux(iVar) .and. any(fluxCount%var(iVar)%dat==0))then
   print*, 'fluxCount%var(iVar)%dat = ', fluxCount%var(iVar)%dat
   message=trim(message)//'flux '//trim(flux_meta(iVar)%varname)//' was not computed'
   err=20; return
  endif
 end do

 ! use step halving if unable to complete the fully coupled solution in one substep
 if(ixCoupling/=fullyCoupled .or. nSubsteps>1) dtMultiplier=0.5_rkind

 ! compute the melt in each snow and soil layer
  if(nSnow>0) diag_data%var(iLookDIAG%mLayerMeltFreeze)%dat(      1:nSnow  ) = -(mLayerVolFracIce(      1:nSnow  ) - mLayerVolFracIceInit(      1:nSnow  ))*iden_ice
  diag_data%var(iLookDIAG%mLayerMeltFreeze)%dat(nSnow+1:nLayers) = -(mLayerVolFracIce(nSnow+1:nLayers) - mLayerVolFracIceInit(nSnow+1:nLayers))*iden_water
             
 ! end associate statements
 end associate globalVars

 end subroutine opSplittin


 ! **********************************************************************************************************
 ! private subroutine stateFilter: get a mask for the desired state variables
 ! **********************************************************************************************************
 subroutine stateFilter(ixCoupling,ixSolution,ixStateThenDomain,iStateTypeSplit,iDomainSplit,iStateSplit,&
                        indx_data,stateMask,nSubset,err,message)

 USE indexState_module,only:indxSubset                            ! get state indices
 implicit none
 ! input
 integer(i4b),intent(in)         :: ixCoupling                    ! index of coupling method (1,2)
 integer(i4b),intent(in)         :: ixSolution                    ! index of solution method (1,2)
 integer(i4b),intent(in)         :: ixStateThenDomain             ! switch between full domain and sub domains
 integer(i4b),intent(in)         :: iStateTypeSplit               ! index of the state type split
 integer(i4b),intent(in)         :: iDomainSplit                  ! index of the domain split
 integer(i4b),intent(in)         :: iStateSplit                   ! index of the layer split
 type(var_ilength),intent(inout) :: indx_data                     ! indices for a local HRU
 ! output
 logical(lgt),intent(out)        :: stateMask(:)                  ! mask defining desired state variables
 integer(i4b),intent(out)        :: nSubset                       ! number of selected state variables for a given split
 integer(i4b),intent(out)        :: err                           ! error code
 character(*),intent(out)        :: message                       ! error message
 ! local
 integer(i4b),allocatable        :: ixSubset(:)                   ! list of indices in the state subset
 character(len=256)              :: cmessage                      ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! data structures
 associate(&
 ! indices of model state variables
 ixStateType  => indx_data%var(iLookINDEX%ixStateType)%dat      ,& ! intent(in): [i4b(:)] indices defining the type of the state (ixNrgState...)
 ixNrgCanair  => indx_data%var(iLookINDEX%ixNrgCanair)%dat      ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in canopy air space domain
 ixNrgCanopy  => indx_data%var(iLookINDEX%ixNrgCanopy)%dat      ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the canopy domain
 ixHydCanopy  => indx_data%var(iLookINDEX%ixHydCanopy)%dat      ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the canopy domain
 ixNrgLayer   => indx_data%var(iLookINDEX%ixNrgLayer)%dat       ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
 ixHydLayer   => indx_data%var(iLookINDEX%ixHydLayer)%dat       ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
 ixWatAquifer => indx_data%var(iLookINDEX%ixWatAquifer)%dat     ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for water storage in the aquifer
 ixAllState   => indx_data%var(iLookINDEX%ixAllState)%dat       ,& ! intent(in): [i4b(:)] list of indices for all model state variables (1,2,3,...nState)
 ! number of layers
 nSnow        => indx_data%var(iLookINDEX%nSnow)%dat(1)         ,& ! intent(in): [i4b]    number of snow layers
 nSoil        => indx_data%var(iLookINDEX%nSoil)%dat(1)         ,& ! intent(in): [i4b]    number of soil layers
 nLayers      => indx_data%var(iLookINDEX%nLayers)%dat(1)        & ! intent(in): [i4b]    total number of layers
 ) ! data structures
 ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='stateFilter/'

 ! identify splitting option
 select case(ixCoupling)

  ! -----
  ! - fully coupled...
  ! ------------------

  ! use all state variables
  case(fullyCoupled); stateMask(:) = .true.

  ! -----
  ! - splitting by state type...
  ! ----------------------------

  ! initial split by state type
  case(stateTypeSplit)

   ! switch between full domain and sub domains
   select case(ixStateThenDomain)

   ! split into energy and mass
   case(fullDomain)
    select case(iStateTypeSplit)
     case(nrgSplit);  stateMask = (ixStateType==iname_nrgCanair .or. ixStateType==iname_nrgCanopy .or. ixStateType==iname_nrgLayer)
     case(massSplit); stateMask = (ixStateType==iname_liqCanopy .or. ixStateType==iname_liqLayer  .or. ixStateType==iname_lmpLayer .or. ixStateType==iname_watAquifer)
     case default; err=20; message=trim(message)//'unable to identify split based on state type'; return
    end select

   ! split into vegetation, snow, and soil
   case(subDomain)

    ! define state mask
    stateMask=.false. ! (initialize state mask)
    select case(iStateTypeSplit)

     ! define mask for energy
     case(nrgSplit)
      select case(iDomainSplit)
       case(vegSplit)
        if(ixNrgCanair(1)/=integerMissing) stateMask(ixNrgCanair) = .true.  ! energy of the canopy air space
        if(ixNrgCanopy(1)/=integerMissing) stateMask(ixNrgCanopy) = .true.  ! energy of the vegetation canopy
        stateMask(ixNrgLayer(1)) = .true.  ! energy of the upper-most layer in the snow+soil domain
       case(snowSplit);   if(nSnow>1) stateMask(ixNrgLayer(2:nSnow)) = .true.    ! NOTE: (2:) because the top layer in the snow+soil domain included in vegSplit
       case(soilSplit);   stateMask(ixNrgLayer(max(2,nSnow+1):nLayers)) = .true. ! NOTE: max(2,nSnow+1) gives second layer unless more than 2 snow layers
       case(aquiferSplit) ! do nothing: no energy state variable for the aquifer domain
       case default; err=20; message=trim(message)//'unable to identify model sub-domain'; return
      end select

     ! define mask for water
     case(massSplit)
      select case(iDomainSplit)
       case(vegSplit);     if(ixHydCanopy(1)/=integerMissing) stateMask(ixHydCanopy) = .true.  ! hydrology of the vegetation canopy
       case(snowSplit);    stateMask(ixHydLayer(1:nSnow)) = .true.  ! snow hydrology
       case(soilSplit);    stateMask(ixHydLayer(nSnow+1:nLayers)) = .true.  ! soil hydrology
       case(aquiferSplit); if(ixWatAquifer(1)/=integerMissing) stateMask(ixWatAquifer) = .true.  ! aquifer storage
       case default; err=20; message=trim(message)//'unable to identify model sub-domain'; return
      end select

     ! check
     case default; err=20; message=trim(message)//'unable to identify the state type'; return
    end select  ! (split based on state type)

    ! check
    case default; err=20; message=trim(message)//'unable to identify the switch between full domains and sub domains'; return
   end select ! (switch between full domains and sub domains)

  ! check
  case default; err=20; message=trim(message)//'unable to identify coupling method'; return
 end select  ! (selecting solution method)

 !print*, 'stateMask = ', stateMask

 ! identify scalar solutions
 if(ixSolution==scalar)then

  ! get the subset of indices
  call indxSubset(ixSubset, ixAllState, stateMask, err, cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get the mask
  stateMask(:) = .false.
  stateMask( ixSubset(iStateSplit) ) = .true.

  ! check
  if(count(stateMask)/=1)then
   message=trim(message)//'expect size=1 (scalar)'
   err=20; return
  endif

 endif

 ! get the number of selected state variables
 nSubset = count(stateMask)

 ! end associations
 end associate

 end subroutine stateFilter

end module opSplittin_module
