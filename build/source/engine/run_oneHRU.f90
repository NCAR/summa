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

module run_oneHRU_module

! numerical recipes data types
USE nrtype

! data types
USE data_types,only:&
               var_i,                    & ! x%var(:)                (i4b)
               var_d,                    & ! x%var(:)                (dp)
               var_ilength,              & ! x%var(:)%dat            (i4b)
               var_dlength,              & ! x%var(:)%dat            (dp)
               zLookup                     ! x%z(:)%var(:)%lookup(:) (dp)

! access vegetation data
USE globalData,only:greenVegFrac_monthly   ! fraction of green vegetation in each month (0-1)
USE globalData,only:overwriteRSMIN         ! flag to overwrite RSMIN
USE globalData,only:maxSoilLayers          ! Maximum Number of Soil Layers

! provide access to Noah-MP constants
USE module_sf_noahmplsm,only:isWater       ! parameter for water land cover type

! provide access to the named variables that describe elements of parameter structures
USE var_lookup,only:iLookTYPE              ! look-up values for classification of veg, soils etc.
USE var_lookup,only:iLookATTR              ! look-up values for local attributes
USE var_lookup,only:iLookPARAM             ! look-up values for local column model parameters

! provide access to the named variables that describe elements of variable structures
USE var_lookup,only:iLookPROG              ! look-up values for local column model prognostic (state) variables
USE var_lookup,only:iLookDIAG              ! look-up values for local column model diagnostic variables
USE var_lookup,only:iLookINDEX             ! look-up values for local column index variables

! provide access to model decisions
USE globalData,only:model_decisions        ! model decision structure
USE var_lookup,only:iLookDECISIONS         ! look-up values for model decisions

! provide access to the named variables that describe model decisions
USE mDecisions_module,only:&               ! look-up values for LAI decisions
 monthlyTable,& ! LAI/SAI taken directly from a monthly table for different vegetation classes
 specified      ! LAI/SAI computed from green vegetation fraction and winterSAI and summerLAI parameters

! -----------------------------------------------------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------------------------------------------------
! ----- global variables that are modified ------------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------------------------------------------------

! Noah-MP parameters
USE NOAHMP_VEG_PARAMETERS,only:SAIM,LAIM   ! 2-d tables for stem area index and leaf area index (vegType,month)
USE NOAHMP_VEG_PARAMETERS,only:HVT,HVB     ! height at the top and bottom of vegetation (vegType)
USE noahmp_globals,only:RSMIN              ! minimum stomatal resistance (vegType)

! urban vegetation category (could be local)
USE globalData,only:urbanVegCategory       ! vegetation category for urban areas

implicit none
private
public::run_oneHRU

contains

 ! ************************************************************************************************
 ! public subroutine run_oneGRU: simulation for a single GRU
 ! ************************************************************************************************

 ! simulation for a single HRU
 subroutine run_oneHRU(&
                       ! model control
                       hruId,               & ! intent(in):    hruId
                       dt_init,             & ! intent(inout): used to initialize the length of the sub-step for each HRU
                       computeVegFlux,      & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
                       nSnow,nSoil,nLayers, & ! intent(inout): number of snow and soil layers
                       ! data structures (input)
                       timeVec,             & ! intent(in):    model time data
                       typeData,            & ! intent(in):    local classification of soil veg etc. for each HRU
                       attrData,            & ! intent(in):    local attributes for each HRU
                       lookupData,          & ! intent(in):    local lookup tables for each HRU
                       bvarData,            & ! intent(in):    basin-average variables
                       ! data structures (input-output)
                       mparData,            & ! intent(inout): local model parameters
                       indxData,            & ! intent(inout): model indices
                       forcData,            & ! intent(inout): model forcing data
                       progData,            & ! intent(inout): prognostic variables for a local HRU
                       diagData,            & ! intent(inout): diagnostic variables for a local HRU
                       fluxData,            & ! intent(inout): model fluxes for a local HRU
                       ! error control
                       err,message)           ! intent(out):   error control

 ! ----- define downstream subroutines -----------------------------------------------------------------------------------

 USE module_sf_noahmplsm,only:redprm          ! module to assign more Noah-MP parameters
 USE derivforce_module,only:derivforce        ! module to compute derived forcing data
 USE coupled_em_module,only:coupled_em        ! module to run the coupled energy and mass model
 implicit none

 ! ----- define dummy variables ------------------------------------------------------------------------------------------
 
 ! model control
 integer(8)        , intent(in)    :: hruId               ! hruId
 real(rkind)          , intent(inout) :: dt_init             ! used to initialize the length of the sub-step for each HRU
 logical(lgt)      , intent(inout) :: computeVegFlux      ! flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
 integer(i4b)      , intent(inout) :: nSnow,nSoil,nLayers ! number of snow and soil layers
 ! data structures (input)
 integer(i4b)      , intent(in)    :: timeVec(:)          ! int vector               -- model time data
 type(var_i)       , intent(in)    :: typeData            ! x%var(:)                 -- local classification of soil veg etc. for each HRU
 type(var_d)       , intent(in)    :: attrData            ! x%var(:)                 -- local attributes for each HRU
 type(zLookup)     , intent(in)    :: lookupData          ! x%z(:)%var(:)%lookup(:)  -- local lookup tables for each HRU
 type(var_dlength) , intent(in)    :: bvarData            ! x%var(:)%dat -- basin-average variables
 ! data structures (input-output)
 type(var_dlength) , intent(inout) :: mparData            ! x%var(:)%dat -- local (HRU) model parameters
 type(var_ilength) , intent(inout) :: indxData            ! x%var(:)%dat -- model indices
 type(var_d)       , intent(inout) :: forcData            ! x%var(:)     -- model forcing data
 type(var_dlength) , intent(inout) :: progData            ! x%var(:)%dat -- model prognostic (state) variables
 type(var_dlength) , intent(inout) :: diagData            ! x%var(:)%dat -- model diagnostic variables
 type(var_dlength) , intent(inout) :: fluxData            ! x%var(:)%dat -- model fluxes
 ! error control
 integer(i4b)      , intent(out)   :: err                 ! error code
 character(*)      , intent(out)   :: message             ! error message

 ! ----- define local variables ------------------------------------------------------------------------------------------

 ! local variables
 character(len=256)                :: cmessage            ! error message
 real(rkind)          , allocatable   :: zSoilReverseSign(:) ! height at bottom of each soil layer, negative downwards (m)

 ! initialize error control
 err=0; write(message, '(A20,I0,A2)' ) 'run_oneHRU (hruId = ',hruId,')/'
 
 ! ----- hru initialization ---------------------------------------------------------------------------------------------

 ! water pixel: do nothing
 if (typeData%var(iLookTYPE%vegTypeIndex) == isWater) return

 ! get height at bottom of each soil layer, negative downwards (used in Noah MP)
 allocate(zSoilReverseSign(nSoil),stat=err)
 if(err/=0)then
  message=trim(message)//'problem allocating space for zSoilReverseSign'
  err=20; return
 endif
 zSoilReverseSign(:) = -progData%var(iLookPROG%iLayerHeight)%dat(nSnow+1:nLayers)

 ! populate parameters in Noah-MP modules
 ! Passing a maxSoilLayer in order to pass the check for NROOT, that is done to avoid making any changes to Noah-MP code.
 !  --> NROOT from Noah-MP veg tables (as read here) is not used in SUMMA
 call REDPRM(typeData%var(iLookTYPE%vegTypeIndex),      & ! vegetation type index
             typeData%var(iLookTYPE%soilTypeIndex),     & ! soil type
             typeData%var(iLookTYPE%slopeTypeIndex),    & ! slope type index
             zSoilReverseSign,                          & ! * not used: height at bottom of each layer [NOTE: negative] (m)
             maxSoilLayers,                             & ! number of soil layers
             urbanVegCategory)                            ! vegetation category for urban areas

 ! deallocate height at bottom of each soil layer(used in Noah MP)
 deallocate(zSoilReverseSign,stat=err)
 if(err/=0)then
  message=trim(message)//'problem deallocating space for zSoilReverseSign'
  err=20; return
 endif

 ! overwrite the minimum resistance
 if(overwriteRSMIN) RSMIN = mparData%var(iLookPARAM%minStomatalResistance)%dat(1)

 ! overwrite the vegetation height
 HVT(typeData%var(iLookTYPE%vegTypeIndex)) = mparData%var(iLookPARAM%heightCanopyTop)%dat(1)
 HVB(typeData%var(iLookTYPE%vegTypeIndex)) = mparData%var(iLookPARAM%heightCanopyBottom)%dat(1)

 ! overwrite the tables for LAI and SAI
 if(model_decisions(iLookDECISIONS%LAI_method)%iDecision == specified)then
  SAIM(typeData%var(iLookTYPE%vegTypeIndex),:) = mparData%var(iLookPARAM%winterSAI)%dat(1)
  LAIM(typeData%var(iLookTYPE%vegTypeIndex),:) = mparData%var(iLookPARAM%summerLAI)%dat(1)*greenVegFrac_monthly
 end if

 ! ----- hru forcing ----------------------------------------------------------------------------------------------------

 ! compute derived forcing variables
 call derivforce(timeVec,          & ! vector of time information
                 forcData%var,     & ! vector of model forcing data
                 attrData%var,     & ! vector of model attributes
                 mparData,         & ! data structure of model parameters
                 progData,         & ! data structure of model prognostic variables
                 diagData,         & ! data structure of model diagnostic variables
                 fluxData,         & ! data structure of model fluxes
                 err,cmessage)       ! error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! ----- run the model --------------------------------------------------------------------------------------------------

 ! initialize the number of flux calls
 diagData%var(iLookDIAG%numFluxCalls)%dat(1) = 0._rkind

 ! run the model for a single HRU
 call coupled_em(&
                 ! model control
                 hruId,            & ! intent(in):    hruId
                 dt_init,          & ! intent(inout): initial time step
                 computeVegFlux,   & ! intent(inout): flag to indicate if we are computing fluxes over vegetation
                 ! data structures (input)
                 typeData,         & ! intent(in):    local classification of soil veg etc. for each HRU
                 attrData,         & ! intent(in):    local attributes for each HRU
                 forcData,         & ! intent(in):    model forcing data
                 mparData,         & ! intent(in):    model parameters
                 bvarData,         & ! intent(in):    basin-average model variables
                 lookupData,       & ! intent(in):    lookup tables
                 ! data structures (input-output)
                 indxData,         & ! intent(inout): model indices
                 progData,         & ! intent(inout): model prognostic variables for a local HRU
                 diagData,         & ! intent(inout): model diagnostic variables for a local HRU
                 fluxData,         & ! intent(inout): model fluxes for a local HRU
                 ! error control
                 err,cmessage)       ! intent(out): error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! update the number of layers
 nSnow   = indxData%var(iLookINDEX%nSnow)%dat(1)     ! number of snow layers
 nSoil   = indxData%var(iLookINDEX%nSoil)%dat(1)     ! number of soil layers
 nLayers = indxData%var(iLookINDEX%nLayers)%dat(1)   ! total number of layers

 end subroutine run_oneHRU

end module run_oneHRU_module
