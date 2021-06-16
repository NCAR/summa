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

module systemSolv_module

! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number
USE globalData,only:quadMissing     ! missing quadruple precision number

! access matrix information
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! domain types
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! global metadata
USE globalData,only:flux_meta                        ! metadata on the model fluxes

! constants
USE multiconst,only:&
                    LH_fus,       & ! latent heat of fusion                (J K-1)
                    Tfreeze,      & ! temperature at freezing              (K)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookFORCE      ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
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
public::systemSolv

! control parameters
real(rkind),parameter  :: valueMissing=-9999._rkind     ! missing value
real(rkind),parameter  :: verySmall=1.e-12_rkind        ! a very small number (used to check consistency)
real(rkind),parameter  :: veryBig=1.e+20_rkind          ! a very big number
real(rkind),parameter  :: dx = 1.e-8_rkind              ! finite difference increment

contains


 ! **********************************************************************************************************
 ! public subroutine systemSolv: run the coupled energy-mass model for one timestep
 ! **********************************************************************************************************
 subroutine systemSolv(&
                       ! input: model control
                       dt,                & ! intent(in):    time step (s)
                       nState,            & ! intent(in):    total number of state variables
                       firstSubStep,      & ! intent(in):    flag to denote first sub-step
                       firstFluxCall,     & ! intent(inout): flag to indicate if we are processing the first flux call
                       firstSplitOper,    & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                       computeVegFlux,    & ! intent(in):    flag to denote if computing energy flux over vegetation
                       scalarSolution,    & ! intent(in):    flag to denote if implementing the scalar solution
                       ! input/output: data structures
                       lookup_data,       & ! intent(in):    lookup tables
                       type_data,         & ! intent(in):    type of vegetation and soil
                       attr_data,         & ! intent(in):    spatial attributes
                       forc_data,         & ! intent(in):    model forcing data
                       mpar_data,         & ! intent(in):    model parameters
                       indx_data,         & ! intent(inout): index data
                       prog_data,         & ! intent(inout): model prognostic variables for a local HRU
                       diag_data,         & ! intent(inout): model diagnostic variables for a local HRU
                       flux_temp,         & ! intent(inout): model fluxes for a local HRU
                       bvar_data,         & ! intent(in):    model variables for the local basin
                       model_decisions,   & ! intent(in):    model decisions
                       stateVecInit,      & ! intent(in):    initial state vector
                       ! output
                       deriv_data,        & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                       ixSaturation,      & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                       untappedMelt,      & ! intent(out):   un-tapped melt energy (J m-3 s-1)
                       stateVecTrial,     & ! intent(out):   updated state vector
                       reduceCoupledStep, & ! intent(out):   flag to reduce the length of the coupled step
                       tooMuchMelt,       & ! intent(out):   flag to denote that there was too much melt
                       niter,             & ! intent(out):   number of iterations taken
                       err,message)         ! intent(out):   error code and error message
 ! ---------------------------------------------------------------------------------------
 ! structure allocations
 USE allocspace_module,only:allocLocal                ! allocate local data structures
 ! simulation of fluxes and residuals given a trial state vector
 USE eval8summa_module,only:eval8summa                ! simulation of fluxes and residuals given a trial state vector
 USE summaSolve_module,only:summaSolve                ! calculate the iteration increment, evaluate the new state, and refine if necessary
 USE getVectorz_module,only:getScaling                ! get the scaling vectors
 USE convE2Temp_module,only:temp2ethpy                ! convert temperature to enthalpy
 USE t2enthalpy_module, only:t2enthalpy_T             ! compute enthalpy
 implicit none
 ! ---------------------------------------------------------------------------------------
 ! * dummy variables
 ! ---------------------------------------------------------------------------------------
 ! input: model control
 real(rkind),intent(in)             :: dt                            ! time step (seconds)
 integer(i4b),intent(in)         :: nState                        ! total number of state variables
 logical(lgt),intent(in)         :: firstSubStep                  ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(inout)      :: firstFluxCall                 ! flag to define the first flux call
 logical(lgt),intent(in)         :: firstSplitOper                ! flag to indicate if we are processing the first flux call in a splitting operation
 logical(lgt),intent(in)         :: computeVegFlux                ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 logical(lgt),intent(in)         :: scalarSolution                ! flag to denote if implementing the scalar solution
 ! input/output: data structures
 type(zLookup),intent(in)        :: lookup_data                   ! lookup tables
 type(var_i),intent(in)          :: type_data                     ! type of vegetation and soil
 type(var_d),intent(in)          :: attr_data                     ! spatial attributes
 type(var_d),intent(in)          :: forc_data                     ! model forcing data
 type(var_dlength),intent(in)    :: mpar_data                     ! model parameters
 type(var_ilength),intent(inout) :: indx_data                     ! indices for a local HRU
 type(var_dlength),intent(inout) :: prog_data                     ! prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data                     ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_temp                     ! model fluxes for a local HRU
 type(var_dlength),intent(in)    :: bvar_data                     ! model variables for the local basin
 type(model_options),intent(in)  :: model_decisions(:)            ! model decisions
 real(rkind),intent(in)             :: stateVecInit(:)               ! initial state vector (mixed units)
 ! output: model control
 type(var_dlength),intent(inout) :: deriv_data                    ! derivatives in model fluxes w.r.t. relevant state variables
 integer(i4b),intent(inout)      :: ixSaturation                  ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
 real(rkind),intent(out)            :: untappedMelt(:)               ! un-tapped melt energy (J m-3 s-1)
 real(rkind),intent(out)            :: stateVecTrial(:)              ! trial state vector (mixed units)
 logical(lgt),intent(out)        :: reduceCoupledStep             ! flag to reduce the length of the coupled step
 logical(lgt),intent(out)        :: tooMuchMelt                   ! flag to denote that there was too much melt
 integer(i4b),intent(out)        :: niter                         ! number of iterations taken
 integer(i4b),intent(out)        :: err                           ! error code
 character(*),intent(out)        :: message                       ! error message
 ! *********************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************
 ! ---------------------------------------------------------------------------------------
 ! * general local variables
 ! ---------------------------------------------------------------------------------------
 character(LEN=256)              :: cmessage                      ! error message of downwind routine
 integer(i4b)                    :: iter                          ! iteration index
 integer(i4b)                    :: iVar                          ! index of variable
 integer(i4b)                    :: iLayer                        ! index of layer in the snow+soil domain
 integer(i4b)                    :: iState                        ! index of model state
 integer(i4b)                    :: nLeadDim                      ! length of the leading dimension of the Jacobian matrix (nBands or nState)
 integer(i4b)                    :: local_ixGroundwater           ! local index for groundwater representation
 real(rkind)                        :: bulkDensity                   ! bulk density of a given layer (kg m-3)
 real(rkind)                        :: volEnthalpy                   ! volumetric enthalpy of a given layer (J m-3)
 real(rkind),parameter              :: tempAccelerate=0.00_rkind        ! factor to force initial canopy temperatures to be close to air temperature
 real(rkind),parameter              :: xMinCanopyWater=0.0001_rkind     ! minimum value to initialize canopy water (kg m-2)
 real(rkind),parameter              :: tinyStep=0.000001_rkind          ! stupidly small time step (s)
 ! ------------------------------------------------------------------------------------------------------
 ! * model solver
 ! ------------------------------------------------------------------------------------------------------
 logical(lgt),parameter          :: forceFullMatrix=.false.       ! flag to force the use of the full Jacobian matrix
 integer(i4b)                    :: maxiter                       ! maximum number of iterations
 integer(i4b)                    :: ixMatrix                      ! form of matrix (band diagonal or full matrix)
 integer(i4b)                    :: localMaxIter                  ! maximum number of iterations (depends on solution type)
 integer(i4b), parameter         :: scalarMaxIter=100             ! maximum number of iterations for the scalar solution
 type(var_dlength)               :: flux_init                     ! model fluxes at the start of the time step
 real(rkind),allocatable            :: dBaseflow_dMatric(:,:)        ! derivative in baseflow w.r.t. matric head (s-1)  ! NOTE: allocatable, since not always needed
 real(rkind)                        :: stateVecNew(nState)           ! new state vector (mixed units)
 real(rkind)                        :: fluxVec0(nState)              ! flux vector (mixed units)
 real(rkind)                        :: fScale(nState)                ! characteristic scale of the function evaluations (mixed units)
 real(rkind)                        :: xScale(nState)                ! characteristic scale of the state vector (mixed units)
 real(rkind)                        :: dMat(nState)                  ! diagonal matrix (excludes flux derivatives)
 real(rkind)                        :: sMul(nState)    ! NOTE: qp    ! multiplier for state vector for the residual calculations
 real(rkind)                        :: rVec(nState)    ! NOTE: qp    ! residual vector
 real(rkind)                        :: rAdd(nState)                  ! additional terms in the residual vector
 real(rkind)                        :: fOld,fNew                     ! function values (-); NOTE: dimensionless because scaled
 real(rkind)                        :: xMin,xMax                     ! state minimum and maximum (mixed units)
 logical(lgt)                    :: converged                     ! convergence flag
 logical(lgt)                    :: feasible                      ! feasibility flag
 real(rkind)                        :: resSinkNew(nState)            ! additional terms in the residual vector
 real(rkind)                        :: fluxVecNew(nState)            ! new flux vector
 real(rkind)                        :: resVecNew(nState)  ! NOTE: qp ! new residual vector
 ! ---------------------------------------------------------------------------------------
 ! point to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 globalVars: associate(&
 ! model decisions
 ixGroundwater           => model_decisions(iLookDECISIONS%groundwatr)%iDecision   ,& ! intent(in):    [i4b]    groundwater parameterization
 ixSpatialGroundwater    => model_decisions(iLookDECISIONS%spatial_gw)%iDecision   ,& ! intent(in):    [i4b]    spatial representation of groundwater (local-column or single-basin)
 scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(in): [dp]     temperature of the canopy air space (K)
 scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(in): [dp]     temperature of the vegetation canopy (K)
 scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(in): [dp]     mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(in): [dp]     mass of total water on the vegetation canopy (kg m-2)
 ! model state variables (snow and soil domains)
 mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(in): [dp(:)]  temperature of each snow/soil layer (K)
 mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,& ! intent(in): [dp(:)]  volumetric fraction of ice (-)
 mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(in): [dp(:)]  volumetric fraction of liquid water (-)
 mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(in): [dp(:)]  volumetric fraction of total water (-)
 mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(inout): [dp(:)]  matric head (m)
 ! enthalpy
 scalarCanairEnthalpy    => diag_data%var(iLookDIAG%scalarCanairEnthalpy)%dat(1)   ,&  ! intent(out): [dp]    enthalpy of the canopy air space (J m-3)
 scalarCanopyEnthalpy    => diag_data%var(iLookDIAG%scalarCanopyEnthalpy)%dat(1)   ,&  ! intent(out): [dp]    enthalpy of the vegetation canopy (J m-3)
 mLayerEnthalpy          => diag_data%var(iLookDIAG%mLayerEnthalpy)%dat            ,&  ! intent(out): [dp(:)] enthalpy of the snow+soil layers (J m-3)
 ! check the need to merge snow layers
 mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,& ! intent(in):    [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
 snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)         ,& ! intent(in):    [dp]     scaling parameter for the snow freezing curve (K-1)
 ! accelerate solution for temperature
 airtemp                 => forc_data%var(iLookFORCE%airtemp)                      ,& ! intent(in):    [dp]     temperature of the upper boundary of the snow and soil domains (K)
 ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy air space energy state variable
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy energy state variable
 ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy hydrology state variable (mass)
 ! vector of energy and hydrology indices for the snow and soil domains
 ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
 ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
 ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the soil domain
 nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in):    [i4b]    number of energy state variables in the snow+soil domain
 nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in):    [i4b]    number of hydrology state variables in the snow+soil domain
 nSoilOnlyHyd            => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)         ,& ! intent(in):    [i4b]    number of hydrology state variables in the soil domain
 ! mapping from full domain to the sub-domain
 ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat         ,& ! intent(in):    [i4b]    mapping of full state vector to the state subset
 ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat          ,& ! intent(in):    [i4b]    index of control volume for different domains (veg, snow, soil)
 ! type of state and domain for a given variable
 ixStateType_subset      => indx_data%var(iLookINDEX%ixStateType_subset)%dat       ,& ! intent(in):    [i4b(:)] [state subset] type of desired model state variables
 ixDomainType_subset     => indx_data%var(iLookINDEX%ixDomainType_subset)%dat      ,& ! intent(in):    [i4b(:)] [state subset] domain for desired model state variables
 ! layer geometry
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):    [i4b]    number of snow layers
 nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):    [i4b]    number of soil layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)                & ! intent(in):    [i4b]    total number of layers
 )
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="systemSolv/"

 ! *****
 ! (0) PRELIMINARIES...
 ! ********************

 ! -----
 ! * initialize...
 ! ---------------

 ! check
 if(dt < tinyStep)then
  message=trim(message)//'dt is tiny'
  err=20; return
 endif

 ! initialize the flags
 tooMuchMelt        = .false.   ! too much melt
 reduceCoupledStep  = .false.   ! need to reduce the length of the coupled step

 ! define maximum number of iterations
 maxiter = nint(mpar_data%var(iLookPARAM%maxiter)%dat(1))

 ! modify the groundwater representation for this single-column implementation
 select case(ixSpatialGroundwater)
  case(singleBasin); local_ixGroundwater = noExplicit    ! force no explicit representation of groundwater at the local scale
  case(localColumn); local_ixGroundwater = ixGroundwater ! go with the specified decision
  case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
 end select ! (modify the groundwater representation for this single-column implementation)

 ! allocate space for the model fluxes at the start of the time step
 call allocLocal(flux_meta(:),flux_init,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for the baseflow derivatives
 ! NOTE: needs allocation because only used when baseflow sinks are active
 if(ixGroundwater==qbaseTopmodel)then
  allocate(dBaseflow_dMatric(nSoil,nSoil),stat=err)  ! baseflow depends on total storage in the soil column, hence on matric head in every soil layer
 else
  allocate(dBaseflow_dMatric(0,0),stat=err)          ! allocate zero-length dimnensions to avoid passing around an unallocated matrix
 end if
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the baseflow derivatives'; return; end if

 ! identify the matrix solution method
 ! (the type of matrix used to solve the linear system A.X=B)
 if(local_ixGroundwater==qbaseTopmodel .or. scalarSolution .or. forceFullMatrix)then
  nLeadDim=nState         ! length of the leading dimension
  ixMatrix=ixFullMatrix   ! named variable to denote the full Jacobian matrix
 else
  nLeadDim=nBands         ! length of the leading dimension
  ixMatrix=ixBandMatrix   ! named variable to denote the band-diagonal matrix
 endif

 ! initialize the model fluxes (some model fluxes are not computed in the iterations)
 do iVar=1,size(flux_temp%var)
  flux_init%var(iVar)%dat(:) = flux_temp%var(iVar)%dat(:)
 end do

 ! **************************************************************************************************************************
 ! **************************************************************************************************************************
 ! **************************************************************************************************************************
 ! *** NUMERICAL SOLUTION FOR A GIVEN SUBSTEP AND SPLIT *********************************************************************
 ! **************************************************************************************************************************
 ! **************************************************************************************************************************
 ! **************************************************************************************************************************

 ! -----
 ! * get scaling vectors...
 ! ------------------------

 ! initialize state vectors
 call getScaling(&
                 ! input
                 diag_data,                        & ! intent(in):    model diagnostic variables for a local HRU
                 indx_data,                        & ! intent(in):    indices defining model states and layers
                 ! output
                 fScale,                           & ! intent(out):   function scaling vector (mixed units)
                 xScale,                           & ! intent(out):   variable scaling vector (mixed units)
                 sMul,                             & ! intent(out):   multiplier for state vector (used in the residual calculations)
                 dMat,                             & ! intent(out):   diagonal of the Jacobian matrix (excludes fluxes)
                 err,cmessage)                       ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

 ! -----
 ! * compute the initial function evaluation...
 ! --------------------------------------------

 ! initialize the trial state vectors
 stateVecTrial = stateVecInit

 ! need to intialize canopy water at a positive value
 if(ixVegHyd/=integerMissing)then
  if(stateVecTrial(ixVegHyd) < xMinCanopyWater) stateVecTrial(ixVegHyd) = stateVecTrial(ixVegHyd) + xMinCanopyWater
 endif

 ! try to accelerate solution for energy
 if(ixCasNrg/=integerMissing) stateVecTrial(ixCasNrg) = stateVecInit(ixCasNrg) + (airtemp - stateVecInit(ixCasNrg))*tempAccelerate
 if(ixVegNrg/=integerMissing) stateVecTrial(ixVegNrg) = stateVecInit(ixVegNrg) + (airtemp - stateVecInit(ixVegNrg))*tempAccelerate

 ! compute H_T at the beginning of the data step
 call t2enthalpy_T(&
                ! input: data structures
                diag_data,                   & ! intent(in):  model diagnostic variables for a local HRU
                mpar_data,                   & ! intent(in):  parameter data structure
                indx_data,                   & ! intent(in):  model indices
                lookup_data,                 & ! intent(in):  lookup table data structure
                ! input: state variables for the vegetation canopy
                scalarCanairTemp,            & ! intent(in):  value of canopy air temperature (K)
                scalarCanopyTemp,            & ! intent(in):  value of canopy temperature (K)
                scalarCanopyWat,             & ! intent(in):  value of canopy total water (kg m-2)
                scalarCanopyIce,             & ! intent(in):  value for canopy ice content (kg m-2)
                ! input: variables for the snow-soil domain
                mLayerTemp,                  & ! intent(in):  vector of layer temperature (K)
                mLayerVolFracWat,            & ! intent(in):  vector of volumetric total water content (-)
                mLayerMatricHead,            & ! intent(in):  vector of total water matric potential (m)
                mLayerVolFracIce,            & ! intent(in):  vector of volumetric fraction of ice (-)
                ! output: enthalpy
                scalarCanairEnthalpy,        & ! intent(out): temperature component of enthalpy of the canopy air space (J m-3)
                scalarCanopyEnthalpy,        & ! intent(out): temperature component of enthalpy of the vegetation canopy (J m-3)
                mLayerEnthalpy,              & ! intent(out): temperature component of enthalpy of each snow+soil layer (J m-3)
                ! output: error control
                err,cmessage)                  ! intent(out): error control
if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! compute the flux and the residual vector for a given state vector
 ! NOTE 1: The derivatives computed in eval8summa are used to calculate the Jacobian matrix for the first iteration
 ! NOTE 2: The Jacobian matrix together with the residual vector is used to calculate the first iteration increment
 call eval8summa(&
                 ! input: model control
                 dt,                      & ! intent(in):    length of the time step (seconds)
                 nSnow,                   & ! intent(in):    number of snow layers
                 nSoil,                   & ! intent(in):    number of soil layers
                 nLayers,                 & ! intent(in):    number of layers
                 nState,                  & ! intent(in):    number of state variables in the current subset
                 firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                 firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                 firstSplitOper,          & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                 computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                 scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                 ! input: state vectors
                 stateVecTrial,           & ! intent(in):    model state vector
                 fScale,                  & ! intent(in):    function scaling vector
                 sMul,                    & ! intent(in):    state vector multiplier (used in the residual calculations)
                 ! input: data structures
                 model_decisions,         & ! intent(in):    model decisions
                 lookup_data,             & ! intent(in):    lookup tables
                 type_data,               & ! intent(in):    type of vegetation and soil
                 attr_data,               & ! intent(in):    spatial attributes
                 mpar_data,               & ! intent(in):    model parameters
                 forc_data,               & ! intent(in):    model forcing data
                 bvar_data,               & ! intent(in):    average model variables for the entire basin
                 prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                 indx_data,               & ! intent(in):    index data
                 ! input-output: data structures
                 diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                 flux_init,               & ! intent(inout): model fluxes for a local HRU (initial flux structure)
                 deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                 ! input-output: baseflow
                 ixSaturation,            & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                 dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                 ! output
                 feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                 fluxVec0,                & ! intent(out):   flux vector
                 rAdd,                    & ! intent(out):   additional (sink) terms on the RHS of the state equation
                 rVec,                    & ! intent(out):   residual vector
                 fOld,                    & ! intent(out):   function evaluation
                 err,cmessage)              ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
 if(.not.feasible)then; message=trim(message)//'state vector not feasible'; err=20; return; endif

 ! copy over the initial flux structure since some model fluxes are not computed in the iterations
 do concurrent ( iVar=1:size(flux_meta) )
  flux_temp%var(iVar)%dat(:) = flux_init%var(iVar)%dat(:)
 end do

 ! check the need to merge snow layers
 if(nSnow>0)then
  ! compute the energy required to melt the top snow layer (J m-2)
  bulkDensity = mLayerVolFracIce(1)*iden_ice + mLayerVolFracLiq(1)*iden_water
  volEnthalpy = temp2ethpy(mLayerTemp(1),bulkDensity,snowfrz_scale)
  ! set flag and error codes for too much melt
  if(-volEnthalpy < flux_init%var(iLookFLUX%mLayerNrgFlux)%dat(1)*dt)then
   tooMuchMelt=.true.
   message=trim(message)//'net flux in the top snow layer can melt all the snow in the top layer'
   err=-20; return ! negative error code to denote a warning
  endif
 endif

 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================

 ! **************************
 ! *** MAIN ITERATION LOOP...
 ! **************************

 ! correct the number of iterations
 localMaxIter = merge(scalarMaxIter, maxIter, scalarSolution)

 ! iterate
 do iter=1,localMaxIter

  ! print iteration count
  !print*, '*** iter, maxiter, dt = ', iter, localMaxiter, dt
  !print*, trim(message)//'before summaSolve'

  ! keep track of the number of iterations
  niter = iter+1  ! +1 because xFluxResid was moved outside the iteration loop (for backwards compatibility)

  ! compute the next trial state vector
  !  1) Computes the Jacobian matrix based on derivatives from the last flux evaluation
  !  2) Computes the iteration increment based on Jacobian and residuals from the last flux evaluation
  !  3) Computes new fluxes and derivatives, new residuals, and (if necessary) refines the state vector
  call summaSolve(&
                  ! input: model control
                  dt,                            & ! intent(in):    length of the time step (seconds)
                  iter,                          & ! intent(in):    iteration index
                  nSnow,                         & ! intent(in):    number of snow layers
                  nSoil,                         & ! intent(in):    number of soil layers
                  nLayers,                       & ! intent(in):    total number of layers
                  nLeadDim,                      & ! intent(in):    length of the leading dimension of the Jacobian matrix (either nBands or nState)
                  nState,                        & ! intent(in):    total number of state variables
                  ixMatrix,                      & ! intent(in):    type of matrix (full or band diagonal)
                  firstSubStep,                  & ! intent(in):    flag to indicate if we are processing the first sub-step
                  firstFluxCall,                 & ! intent(inout): flag to indicate if we are processing the first flux call
                  computeVegFlux,                & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                  scalarSolution,                & ! intent(in):    flag to indicate the scalar solution
                  ! input: state vectors
                  stateVecTrial,                 & ! intent(in):    trial state vector
                  xMin,xMax,                     & ! intent(inout): state maximum and minimum
                  fScale,                        & ! intent(in):    function scaling vector
                  xScale,                        & ! intent(in):    "variable" scaling vector, i.e., for state variables
                  rVec,                          & ! intent(in):    residual vector
                  sMul,                          & ! intent(in):    state vector multiplier (used in the residual calculations)
                  dMat,                          & ! intent(inout): diagonal matrix (excludes flux derivatives)
                  fOld,                          & ! intent(in):    old function evaluation
                  ! input: data structures
                  model_decisions,               & ! intent(in):    model decisions
                  lookup_data,                   & ! intent(in):    lookup tables
                  type_data,                     & ! intent(in):    type of vegetation and soil
                  attr_data,                     & ! intent(in):    spatial attributes
                  mpar_data,                     & ! intent(in):    model parameters
                  forc_data,                     & ! intent(in):    model forcing data
                  bvar_data,                     & ! intent(in):    average model variables for the entire basin
                  prog_data,                     & ! intent(in):    model prognostic variables for a local HRU
                  indx_data,                     & ! intent(in):    index data
                  ! input-output: data structures
                  diag_data,                     & ! intent(inout): model diagnostic variables for a local HRU
                  flux_temp,                     & ! intent(inout): model fluxes for a local HRU (temporary structure)
                  deriv_data,                    & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                  ! input-output: baseflow
                  ixSaturation,                  & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                  dBaseflow_dMatric,             & ! intent(inout): derivative in baseflow w.r.t. matric head (s-1)
                  ! output
                  stateVecNew,                   & ! intent(out):   new state vector
                  fluxVecNew,                    & ! intent(out):   new flux vector
                  resSinkNew,                    & ! intent(out):   additional (sink) terms on the RHS of the state equa
                  resVecNew,                     & ! intent(out):   new residual vector
                  fNew,                          & ! intent(out):   new function evaluation
                  converged,                     & ! intent(out):   convergence flag
                  err,cmessage)                    ! intent(out):   error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  !print*, err,trim(cmessage)

  ! update function evaluation, residual vector, and states
  ! NOTE 1: The derivatives computed in summaSolve are used to calculate the Jacobian matrix at the next iteration
  ! NOTE 2: The Jacobian matrix together with the residual vector is used to calculate the new iteration increment

  ! save functions and residuals
  fOld          = fNew
  rVec          = resVecNew
  stateVecTrial = stateVecNew

  ! print progress
  !write(*,'(a,10(f16.14,1x))') 'rVec                  = ', rVec           ( min(nState,iJac1) : min(nState,iJac2) )
  !write(*,'(a,10(f16.10,1x))') 'fluxVecNew            = ', fluxVecNew     ( min(nState,iJac1) : min(nState,iJac2) )*dt
  !write(*,'(a,10(f16.10,1x))') 'stateVecTrial         = ', stateVecTrial  ( min(nState,iJac1) : min(nState,iJac2) )
  !print*, 'PAUSE: check states and fluxes'; read(*,*)

  ! exit iteration loop if converged
  if(converged) exit

  ! check convergence
  if(iter==localMaxiter)then
   message=trim(message)//'failed to converge'
   err=-20; return
  endif
  !print*, 'PAUSE: iterating'; read(*,*)

 end do  ! iterating
 !print*, 'PAUSE: after iterations'; read(*,*)

 ! -----
 ! * update states...
 ! ------------------

 ! set untapped melt energy to zero
 untappedMelt(:) = 0._rkind

 ! update temperatures (ensure new temperature is consistent with the fluxes)
 if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   iState = ixSnowSoilNrg(iLayer)
   stateVecTrial(iState) = stateVecInit(iState) + (fluxVecNew(iState)*dt + resSinkNew(iState))/real(sMul(iState), rkind)
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! update volumetric water content in the snow (ensure change in state is consistent with the fluxes)
 ! NOTE: for soil water balance is constrained within the iteration loop
 if(nSnowSoilHyd>0)then
  do concurrent (iLayer=1:nSnow,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing water state variables in the snow domain)
   iState = ixSnowSoilHyd(iLayer)
   stateVecTrial(iState) = stateVecInit(iState) + (fluxVecNew(iState)*dt + resSinkNew(iState))
  end do  ! looping through non-missing water state variables in the soil domain
 endif
 
 if ( allocated(dBaseflow_dMatric) ) deallocate(dBaseflow_dMatric)  

 ! end associate statements
 end associate globalVars

 end subroutine systemSolv

end module systemSolv_module
