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

module vegNrgFlux_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_i,            & ! data vector (i4b)
                    var_d,            & ! data vector (dp)
                    var_ilength,      & ! data vector with variable length dimension (i4b)
                    var_dlength,      & ! data vector with variable length dimension (dp)
                    model_options       ! defines the model decisions

! indices that define elements of the data structures
USE var_lookup,only:iLookTYPE           ! named variables for structure elements
USE var_lookup,only:iLookPROG           ! named variables for structure elements
USE var_lookup,only:iLookDIAG           ! named variables for structure elements
USE var_lookup,only:iLookFLUX           ! named variables for structure elements
USE var_lookup,only:iLookFORCE          ! named variables for structure elements
USE var_lookup,only:iLookPARAM          ! named variables for structure elements
USE var_lookup,only:iLookINDEX          ! named variables for structure elements
USE var_lookup,only:iLookBVAR           ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS                               ! named variables for elements of the decision structure

! constants
USE multiconst,only:gravity    ! acceleration of gravity              (m s-2)
USE multiconst,only:vkc        ! von Karman's constant                (-)
USE multiconst,only:w_ratio    ! molecular ratio water to dry air     (-)
USE multiconst,only:R_wv       ! gas constant for water vapor         (Pa K-1 m3 kg-1; J kg-1 K-1)
USE multiconst,only:Cp_air     ! specific heat of air                 (J kg-1 K-1)
USE multiconst,only:Cp_ice     ! specific heat of ice                 (J kg-1 K-1)
USE multiconst,only:Cp_soil    ! specific heat of soil                (J kg-1 K-1)
USE multiconst,only:Cp_water   ! specific heat of liquid water        (J kg-1 K-1)
USE multiconst,only:Tfreeze    ! temperature at freezing              (K)
USE multiconst,only:LH_fus     ! latent heat of fusion                (J kg-1)
USE multiconst,only:LH_vap     ! latent heat of vaporization          (J kg-1)
USE multiconst,only:LH_sub     ! latent heat of sublimation           (J kg-1)
USE multiconst,only:sb         ! Stefan Boltzman constant             (W m-2 K-4)
USE multiconst,only:iden_air   ! intrinsic density of air             (kg m-3)
USE multiconst,only:iden_ice   ! intrinsic density of ice             (kg m-3)
USE multiconst,only:iden_water ! intrinsic density of liquid water    (kg m-3)

! look-up values for method used to compute derivative
USE mDecisions_module,only:  &
 numerical,                  & ! numerical solution
 analytical                    ! analytical solution

! look-up values for choice of boundary conditions for thermodynamics
USE mDecisions_module,only:  &
 prescribedTemp,             & ! prescribed temperature
 energyFlux,                 & ! energy flux
 zeroFlux                      ! zero flux

! look-up values for the choice of parameterization for vegetation roughness length and displacement height
USE mDecisions_module,only:  &
 Raupach_BLM1994,            & ! Raupach (BLM 1994) "Simplified expressions..."
 CM_QJRMS1988,               & ! Choudhury and Monteith (QJRMS 1988) "A four layer model for the heat budget..."
 vegTypeTable                  ! constant parameters dependent on the vegetation type

! look-up values for the choice of parameterization for canopy emissivity
USE mDecisions_module,only:  &
 simplExp,                   & ! simple exponential function
 difTrans                      ! parameterized as a function of diffuse transmissivity

! look-up values for the choice of canopy wind profile
USE mDecisions_module,only:  &
 exponential,                & ! exponential wind profile extends to the surface
 logBelowCanopy                ! logarithmic profile below the vegetation canopy

! look-up values for choice of stability function
USE mDecisions_module,only:  &
 standard,                   & ! standard MO similarity, a la Anderson (1976)
 louisInversePower,          & ! Louis (1979) inverse power function
 mahrtExponential              ! Mahrt (1987) exponential

! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:  &
 localColumn,                & ! separate groundwater representation in each local soil column
 singleBasin                   ! single groundwater store over the entire basin

! -------------------------------------------------------------------------------------------------
! privacy
implicit none
private
public::vegNrgFlux
public::wettedFrac
! dimensions
integer(i4b),parameter        :: nBands=2      ! number of spectral bands for shortwave radiation
! named variables
integer(i4b),parameter        :: ist     = 1   ! Surface type:  IST=1 => soil;  IST=2 => lake
integer(i4b),parameter        :: isc     = 4   ! Soil color type
integer(i4b),parameter        :: ice     = 0   ! Surface type:  ICE=0 => soil;  ICE=1 => sea-ice
! spatial indices
integer(i4b),parameter        :: iLoc    = 1   ! i-location
integer(i4b),parameter        :: jLoc    = 1   ! j-location
! algorithmic parameters
real(rkind),parameter     :: missingValue=-9999._rkind   ! missing value, used when diagnostic or state variables are undefined
real(rkind),parameter     :: verySmall=1.e-6_rkind       ! used as an additive constant to check if substantial difference among real numbers
real(rkind),parameter     :: tinyVal=epsilon(1._rkind)   ! used as an additive constant to check if substantial difference among real numbers
real(rkind),parameter     :: mpe=1.e-6_rkind             ! prevents overflow error if division by zero
real(rkind),parameter     :: dx=1.e-11_rkind             ! finite difference increment
! control
logical(lgt)           :: printflag            ! flag to turn on printing
contains

 ! *******************************************************************************************************
 ! public subroutine vegNrgFlux: muster program to compute energy fluxes at vegetation and ground surfaces
 ! *******************************************************************************************************
 subroutine vegNrgFlux(&
                       ! input: model control
                       firstSubStep,                            & ! intent(in): flag to indicate if we are processing the first sub-step
                       firstFluxCall,                           & ! intent(in): flag to indicate if we are processing the first flux call
                       computeVegFlux,                          & ! intent(in): flag to indicate if we need to compute fluxes over vegetation

                       ! input: model state variables
                       upperBoundTemp,                          & ! intent(in): temperature of the upper boundary (K) --> NOTE: use air temperature
                       canairTempTrial,                         & ! intent(in): trial value of the canopy air space temperature (K)
                       canopyTempTrial,                         & ! intent(in): trial value of canopy temperature (K)
                       groundTempTrial,                         & ! intent(in): trial value of ground temperature (K)
                       canopyIceTrial,                          & ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)
                       canopyLiqTrial,                          & ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)

                       ! input: model derivatives
                       dCanLiq_dTcanopy,                        & ! intent(in): derivative in canopy liquid w.r.t. canopy temperature (kg m-2 K-1)

                       ! input/output: data structures
                       type_data,                               & ! intent(in):    type of vegetation and soil
                       forc_data,                               & ! intent(in):    model forcing data
                       mpar_data,                               & ! intent(in):    model parameters
                       indx_data,                               & ! intent(in):    state vector geometry
                       prog_data,                               & ! intent(in):    model prognostic variables for a local HRU
                       diag_data,                               & ! intent(inout): model diagnostic variables for a local HRU
                       flux_data,                               & ! intent(inout): model fluxes for a local HRU
                       bvar_data,                               & ! intent(in):    model variables for the local basin
                       model_decisions,                         & ! intent(in):    model decisions

                       ! output: liquid water fluxes associated with evaporation/transpiration (needed for coupling)
                       returnCanopyTranspiration,               & ! intent(out): canopy transpiration (kg m-2 s-1)
                       returnCanopyEvaporation,                 & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                       returnGroundEvaporation,                 & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)

                       ! output: fluxes
                       canairNetFlux,                           & ! intent(out): net energy flux for the canopy air space (W m-2)
                       canopyNetFlux,                           & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                       groundNetFlux,                           & ! intent(out): net energy flux for the ground surface (W m-2)

                       ! output: energy flux derivatives
                       dCanairNetFlux_dCanairTemp,              & ! intent(out): derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
                       dCanairNetFlux_dCanopyTemp,              & ! intent(out): derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
                       dCanairNetFlux_dGroundTemp,              & ! intent(out): derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
                       dCanopyNetFlux_dCanairTemp,              & ! intent(out): derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
                       dCanopyNetFlux_dCanopyTemp,              & ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
                       dCanopyNetFlux_dGroundTemp,              & ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
                       dGroundNetFlux_dCanairTemp,              & ! intent(out): derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
                       dGroundNetFlux_dCanopyTemp,              & ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
                       dGroundNetFlux_dGroundTemp,              & ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)

                       ! output liquid water flux derivarives (canopy evap)
                       dCanopyEvaporation_dCanLiq,              & ! intent(out): derivative in canopy evaporation w.r.t. canopy liquid water content (s-1)
                       dCanopyEvaporation_dTCanair,             & ! intent(out): derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
                       dCanopyEvaporation_dTCanopy,             & ! intent(out): derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
                       dCanopyEvaporation_dTGround,             & ! intent(out): derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)

                       ! output: liquid water flux derivarives (ground evap)
                       dGroundEvaporation_dCanLiq,              & ! intent(out): derivative in ground evaporation w.r.t. canopy liquid water content (s-1)
                       dGroundEvaporation_dTCanair,             & ! intent(out): derivative in ground evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
                       dGroundEvaporation_dTCanopy,             & ! intent(out): derivative in ground evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
                       dGroundEvaporation_dTGround,             & ! intent(out): derivative in ground evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)

                       ! output: cross derivative terms
                       dCanopyNetFlux_dCanLiq,                  & ! intent(out): derivative in net canopy fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                       dGroundNetFlux_dCanLiq,                  & ! intent(out): derivative in net ground fluxes w.r.t. canopy liquid water content (J kg-1 s-1)

                       ! output: error control
                       err,message)                               ! intent(out): error control

 ! utilities
 USE expIntegral_module,only:expInt                               ! function to calculate the exponential integral
 ! conversion functions
 USE conv_funcs_module,only:satVapPress                           ! function to compute the saturated vapor pressure (Pa)
 USE conv_funcs_module,only:getLatentHeatValue                    ! function to identify latent heat of vaporization/sublimation (J kg-1)
 ! stomatal resistance
 USE stomResist_module,only:stomResist                            ! subroutine to calculate stomatal resistance
 ! compute energy and mass fluxes for vegetation
 implicit none

 ! ---------------------------------------------------------------------------------------
 ! * dummy variables
 ! ---------------------------------------------------------------------------------------
 ! input: model control
 logical(lgt),intent(in)         :: firstSubStep                    ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)         :: firstFluxCall                   ! flag to indicate if we are processing the first flux call
 logical(lgt),intent(in)         :: computeVegFlux                  ! flag to indicate if computing fluxes over vegetation

 ! input: model state variables
 real(rkind),intent(in)             :: upperBoundTemp                  ! temperature of the upper boundary (K) --> NOTE: use air temperature
 real(rkind),intent(in)             :: canairTempTrial                 ! trial value of canopy air space temperature (K)
 real(rkind),intent(in)             :: canopyTempTrial                 ! trial value of canopy temperature (K)
 real(rkind),intent(in)             :: groundTempTrial                 ! trial value of ground temperature (K)
 real(rkind),intent(in)             :: canopyIceTrial                  ! trial value of mass of ice on the vegetation canopy (kg m-2)
 real(rkind),intent(in)             :: canopyLiqTrial                  ! trial value of mass of liquid water on the vegetation canopy (kg m-2)

 ! input: model derivatives
 real(rkind),intent(in)             :: dCanLiq_dTcanopy                ! intent(in): derivative in canopy liquid w.r.t. canopy temperature (kg m-2 K-1)

 ! input/output: data structures
 type(var_i),intent(in)          :: type_data                       ! type of vegetation and soil
 type(var_d),intent(in)          :: forc_data                       ! model forcing data
 type(var_dlength),intent(in)    :: mpar_data                       ! model parameters
 type(var_ilength),intent(in)    :: indx_data                       ! state vector geometry
 type(var_dlength),intent(in)    :: prog_data                       ! prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data                       ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data                       ! model fluxes for a local HRU
 type(var_dlength),intent(in)    :: bvar_data                       ! model variables for the local basin
 type(model_options),intent(in)  :: model_decisions(:)              ! model decisions

 ! output: liquid water fluxes associated with evaporation/transpiration (needed for coupling)
 real(rkind),intent(out)            :: returnCanopyTranspiration       ! canopy transpiration (kg m-2 s-1)
 real(rkind),intent(out)            :: returnCanopyEvaporation         ! canopy evaporation/condensation (kg m-2 s-1)
 real(rkind),intent(out)            :: returnGroundEvaporation         ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)

 ! output: fluxes
 real(rkind),intent(out)            :: canairNetFlux                   ! net energy flux for the canopy air space (W m-2)
 real(rkind),intent(out)            :: canopyNetFlux                   ! net energy flux for the vegetation canopy (W m-2)
 real(rkind),intent(out)            :: groundNetFlux                   ! net energy flux for the ground surface (W m-2)

 ! output: energy flux derivatives
 real(rkind),intent(out)            :: dCanairNetFlux_dCanairTemp      ! derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
 real(rkind),intent(out)            :: dCanairNetFlux_dCanopyTemp      ! derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
 real(rkind),intent(out)            :: dCanairNetFlux_dGroundTemp      ! derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
 real(rkind),intent(out)            :: dCanopyNetFlux_dCanairTemp      ! derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
 real(rkind),intent(out)            :: dCanopyNetFlux_dCanopyTemp      ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
 real(rkind),intent(out)            :: dCanopyNetFlux_dGroundTemp      ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
 real(rkind),intent(out)            :: dGroundNetFlux_dCanairTemp      ! derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
 real(rkind),intent(out)            :: dGroundNetFlux_dCanopyTemp      ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
 real(rkind),intent(out)            :: dGroundNetFlux_dGroundTemp      ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)

 ! output: liquid flux derivatives (canopy evap)
 real(rkind),intent(out)            :: dCanopyEvaporation_dCanLiq      ! derivative in canopy evaporation w.r.t. canopy liquid water content (s-1)
 real(rkind),intent(out)            :: dCanopyEvaporation_dTCanair     ! derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
 real(rkind),intent(out)            :: dCanopyEvaporation_dTCanopy     ! derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
 real(rkind),intent(out)            :: dCanopyEvaporation_dTGround     ! derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)

 ! output: liquid flux derivatives (ground evap)
 real(rkind),intent(out)            :: dGroundEvaporation_dCanLiq      ! derivative in ground evaporation w.r.t. canopy liquid water content (s-1)
 real(rkind),intent(out)            :: dGroundEvaporation_dTCanair     ! derivative in ground evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
 real(rkind),intent(out)            :: dGroundEvaporation_dTCanopy     ! derivative in ground evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
 real(rkind),intent(out)            :: dGroundEvaporation_dTGround     ! derivative in ground evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)

 ! output: cross derivative terms
 real(rkind),intent(out)            :: dCanopyNetFlux_dCanLiq          ! derivative in net canopy fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
 real(rkind),intent(out)            :: dGroundNetFlux_dCanLiq          ! derivative in net ground fluxes w.r.t. canopy liquid water content (J kg-1 s-1)

 ! output: error control
 integer(i4b),intent(out)        :: err                             ! error code
 character(*),intent(out)        :: message                         ! error message

 ! ---------------------------------------------------------------------------------------
 ! * local variables
 ! ---------------------------------------------------------------------------------------
 ! local (general)
 character(LEN=256)             :: cmessage                         ! error message of downwind routine
 real(rkind)                       :: VAI                              ! vegetation area index (m2 m-2)
 real(rkind)                       :: exposedVAI                       ! exposed vegetation area index (m2 m-2)
 real(rkind)                       :: totalCanopyWater                 ! total water on the vegetation canopy (kg m-2)
 real(rkind)                       :: scalarAquiferStorage             ! aquifer storage (m)

 ! local (compute numerical derivatives)
 integer(i4b),parameter         :: unperturbed=1                    ! named variable to identify the case of unperturbed state variables
 integer(i4b),parameter         :: perturbStateGround=2             ! named variable to identify the case where we perturb the ground temperature
 integer(i4b),parameter         :: perturbStateCanopy=3             ! named variable to identify the case where we perturb the canopy temperature
 integer(i4b),parameter         :: perturbStateCanair=4             ! named variable to identify the case where we perturb the canopy air temperature
 integer(i4b),parameter         :: perturbStateCanLiq=5             ! named variable to identify the case where we perturb the canopy liquid water content
 integer(i4b)                   :: itry                             ! index of flux evaluation
 integer(i4b)                   :: nFlux                            ! number of flux evaluations
 real(rkind)                       :: groundTemp                       ! value of ground temperature used in flux calculations (may be perturbed)
 real(rkind)                       :: canopyTemp                       ! value of canopy temperature used in flux calculations (may be perturbed)
 real(rkind)                       :: canairTemp                       ! value of canopy air temperature used in flux calculations (may be perturbed)
 real(rkind)                       :: try0,try1                        ! trial values to evaluate specific derivatives (testing only)

 ! local (saturation vapor pressure of veg)
 real(rkind)                       :: TV_celcius                       ! vegetaion temperature (C)
 real(rkind)                       :: TG_celcius                       ! ground temperature (C)
 real(rkind)                       :: dSVPCanopy_dCanopyTemp           ! derivative in canopy saturated vapor pressure w.r.t. vegetation temperature (Pa/K)
 real(rkind)                       :: dSVPGround_dGroundTemp           ! derivative in ground saturated vapor pressure w.r.t. ground temperature (Pa/K)

 ! local (wetted canopy area)
 real(rkind)                       :: fracLiquidCanopy                 ! fraction of liquid water in the canopy (-)
 real(rkind)                       :: canopyWetFraction                ! trial value of the canopy wetted fraction (-)
 real(rkind)                       :: dCanopyWetFraction_dWat          ! derivative in wetted fraction w.r.t. canopy total water (kg-1 m2)
 real(rkind)                       :: dCanopyWetFraction_dT            ! derivative in wetted fraction w.r.t. canopy temperature (K-1)

 ! local (longwave radiation)
 real(rkind)                       :: expi                             ! exponential integral
 real(rkind)                       :: scaleLAI                         ! scaled LAI (computing diffuse transmissivity)
 real(rkind)                       :: diffuseTrans                     ! diffuse transmissivity (-)
 real(rkind)                       :: groundEmissivity                 ! emissivity of the ground surface (-)
 real(rkind),parameter             :: vegEmissivity=0.98_rkind            ! emissivity of vegetation (0.9665 in JULES) (-)
 real(rkind),parameter             :: soilEmissivity=0.98_rkind           ! emmisivity of the soil (0.9665 in JULES) (-)
 real(rkind),parameter             :: snowEmissivity=0.99_rkind           ! emissivity of snow (-)
 real(rkind)                       :: dLWNetCanopy_dTCanopy            ! derivative in net canopy radiation w.r.t. canopy temperature (W m-2 K-1)
 real(rkind)                       :: dLWNetGround_dTGround            ! derivative in net ground radiation w.r.t. ground temperature (W m-2 K-1)
 real(rkind)                       :: dLWNetCanopy_dTGround            ! derivative in net canopy radiation w.r.t. ground temperature (W m-2 K-1)
 real(rkind)                       :: dLWNetGround_dTCanopy            ! derivative in net ground radiation w.r.t. canopy temperature (W m-2 K-1)

 ! local (aerodynamic resistance)
 real(rkind)                       :: scalarCanopyStabilityCorrection_old    ! stability correction for the canopy (-)
 real(rkind)                       :: scalarGroundStabilityCorrection_old    ! stability correction for the ground surface (-)

 ! local (turbulent heat transfer)
 real(rkind)                       :: z0Ground                         ! roughness length of the ground (ground below the canopy or non-vegetated surface) (m)
 real(rkind)                       :: soilEvapFactor                   ! soil water control on evaporation from non-vegetated surfaces
 real(rkind)                       :: soilRelHumidity_noSnow           ! relative humidity in the soil pores [0-1]
 real(rkind)                       :: scalarLeafConductance            ! leaf conductance (m s-1)
 real(rkind)                       :: scalarCanopyConductance          ! canopy conductance (m s-1)
 real(rkind)                       :: scalarGroundConductanceSH        ! ground conductance for sensible heat (m s-1)
 real(rkind)                       :: scalarGroundConductanceLH        ! ground conductance for latent heat -- includes soil resistance (m s-1)
 real(rkind)                       :: scalarEvapConductance            ! conductance for evaporation (m s-1)
 real(rkind)                       :: scalarTransConductance           ! conductance for transpiration (m s-1)
 real(rkind)                       :: scalarTotalConductanceSH         ! total conductance for sensible heat (m s-1)
 real(rkind)                       :: scalarTotalConductanceLH         ! total conductance for latent heat (m s-1)
 real(rkind)                       :: dGroundResistance_dTGround       ! derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
 real(rkind)                       :: dGroundResistance_dTCanopy       ! derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
 real(rkind)                       :: dGroundResistance_dTCanair       ! derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
 real(rkind)                       :: dCanopyResistance_dTCanopy       ! derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
 real(rkind)                       :: dCanopyResistance_dTCanair       ! derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
 real(rkind)                       :: turbFluxCanair                   ! total turbulent heat fluxes exchanged at the canopy air space (W m-2)
 real(rkind)                       :: turbFluxCanopy                   ! total turbulent heat fluxes from the canopy to the canopy air space (W m-2)
 real(rkind)                       :: turbFluxGround                   ! total turbulent heat fluxes from the ground to the canopy air space (W m-2)

 ! local (turbulent heat transfer -- compute numerical derivatives)
 ! (temporary scalar resistances when states are perturbed)
 real(rkind)                       :: trialLeafResistance              ! mean leaf boundary layer resistance per unit leaf area (s m-1)
 real(rkind)                       :: trialGroundResistance            ! below canopy aerodynamic resistance (s m-1)
 real(rkind)                       :: trialCanopyResistance            ! above canopy aerodynamic resistance (s m-1)
 real(rkind)                       :: notUsed_RiBulkCanopy             ! bulk Richardson number for the canopy (-)
 real(rkind)                       :: notUsed_RiBulkGround             ! bulk Richardson number for the ground surface (-)
 real(rkind)                       :: notUsed_z0Canopy                 ! roughness length of the vegetation canopy (m)
 real(rkind)                       :: notUsed_WindReductionFactor      ! canopy wind reduction factor (-)
 real(rkind)                       :: notUsed_ZeroPlaneDisplacement    ! zero plane displacement (m)
 real(rkind)                       :: notUsed_scalarCanopyStabilityCorrection  ! stability correction for the canopy (-)
 real(rkind)                       :: notUsed_scalarGroundStabilityCorrection  ! stability correction for the ground surface (-)
 real(rkind)                       :: notUsed_EddyDiffusCanopyTop      ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
 real(rkind)                       :: notUsed_FrictionVelocity         ! friction velocity (m s-1)
 real(rkind)                       :: notUsed_WindspdCanopyTop         ! windspeed at the top of the canopy (m s-1)
 real(rkind)                       :: notUsed_WindspdCanopyBottom      ! windspeed at the height of the bottom of the canopy (m s-1)
 real(rkind)                       :: notUsed_dGroundResistance_dTGround  ! derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
 real(rkind)                       :: notUsed_dGroundResistance_dTCanopy  ! derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
 real(rkind)                       :: notUsed_dGroundResistance_dTCanair  ! derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
 real(rkind)                       :: notUsed_dCanopyResistance_dTCanopy  ! derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
 real(rkind)                       :: notUsed_dCanopyResistance_dTCanair  ! derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)

 ! (fluxes after perturbations in model states -- canopy air space)
 real(rkind)                       :: turbFluxCanair_dStateCanair      ! turbulent exchange from the canopy air space to the atmosphere, after canopy air temperature is perturbed (W m-2)
 real(rkind)                       :: turbFluxCanair_dStateCanopy      ! turbulent exchange from the canopy air space to the atmosphere, after canopy temperature is perturbed (W m-2)
 real(rkind)                       :: turbFluxCanair_dStateGround      ! turbulent exchange from the canopy air space to the atmosphere, after ground temperature is perturbed (W m-2)
 real(rkind)                       :: turbFluxCanair_dStateCanliq      ! turbulent exchange from the canopy air space to the atmosphere, after canopy liquid water content is perturbed (W m-2)
 ! (fluxes after perturbations in model states -- vegetation canopy)
 real(rkind)                       :: turbFluxCanopy_dStateCanair      ! total turbulent heat fluxes from the canopy to the canopy air space, after canopy air temperature is perturbed (W m-2)
 real(rkind)                       :: turbFluxCanopy_dStateCanopy      ! total turbulent heat fluxes from the canopy to the canopy air space, after canopy temperature is perturbed (W m-2)
 real(rkind)                       :: turbFluxCanopy_dStateGround      ! total turbulent heat fluxes from the canopy to the canopy air space, after ground temperature is perturbed (W m-2)
 real(rkind)                       :: turbFluxCanopy_dStateCanLiq      ! total turbulent heat fluxes from the canopy to the canopy air space, after canopy liquid water content is perturbed (W m-2)

 ! (fluxes after perturbations in model states -- ground surface)
 real(rkind)                       :: turbFluxGround_dStateCanair      ! total turbulent heat fluxes from the ground to the canopy air space, after canopy air temperature is perturbed (W m-2)
 real(rkind)                       :: turbFluxGround_dStateCanopy      ! total turbulent heat fluxes from the ground to the canopy air space, after canopy temperature is perturbed (W m-2)
 real(rkind)                       :: turbFluxGround_dStateGround      ! total turbulent heat fluxes from the ground to the canopy air space, after ground temperature is perturbed (W m-2)
 real(rkind)                       :: turbFluxGround_dStateCanLiq      ! total turbulent heat fluxes from the ground to the canopy air space, after canopy liquid water content is perturbed (W m-2)

 ! (fluxes after perturbations in model states -- canopy evaporation)
 real(rkind)                       :: latHeatCanEvap_dStateCanair      ! canopy evaporation after canopy air temperature is perturbed (W m-2)
 real(rkind)                       :: latHeatCanEvap_dStateCanopy      ! canopy evaporation after canopy temperature is perturbed (W m-2)
 real(rkind)                       :: latHeatCanEvap_dStateGround      ! canopy evaporation after ground temperature is perturbed (W m-2)
 real(rkind)                       :: latHeatCanEvap_dStateCanLiq      ! canopy evaporation after canopy liquid water content is perturbed (W m-2)

 ! (flux derivatives -- canopy air space)
 real(rkind)                       :: dTurbFluxCanair_dTCanair         ! derivative in net canopy air space fluxes w.r.t. canopy air temperature (W m-2 K-1)
 real(rkind)                       :: dTurbFluxCanair_dTCanopy         ! derivative in net canopy air space fluxes w.r.t. canopy temperature (W m-2 K-1)
 real(rkind)                       :: dTurbFluxCanair_dTGround         ! derivative in net canopy air space fluxes w.r.t. ground temperature (W m-2 K-1)
 real(rkind)                       :: dTurbFluxCanair_dCanLiq          ! derivative in net canopy air space fluxes w.r.t. canopy liquid water content (J kg-1 s-1)

 ! (flux derivatives -- vegetation canopy)
 real(rkind)                       :: dTurbFluxCanopy_dTCanair         ! derivative in net canopy turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
 real(rkind)                       :: dTurbFluxCanopy_dTCanopy         ! derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
 real(rkind)                       :: dTurbFluxCanopy_dTGround         ! derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
 real(rkind)                       :: dTurbFluxCanopy_dCanLiq          ! derivative in net canopy turbulent fluxes w.r.t. canopy liquid water content (J kg-1 s-1)

 ! (flux derivatives -- ground surface)
 real(rkind)                       :: dTurbFluxGround_dTCanair         ! derivative in net ground turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
 real(rkind)                       :: dTurbFluxGround_dTCanopy         ! derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
 real(rkind)                       :: dTurbFluxGround_dTGround         ! derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
 real(rkind)                       :: dTurbFluxGround_dCanLiq          ! derivative in net ground turbulent fluxes w.r.t. canopy liquid water content (J kg-1 s-1)

 ! (liquid water flux derivatives -- canopy evap)
 real(rkind)                       :: dLatHeatCanopyEvap_dCanLiq       ! derivative in latent heat of canopy evaporation w.r.t. canopy liquid water content (W kg-1)
 real(rkind)                       :: dLatHeatCanopyEvap_dTCanair      ! derivative in latent heat of canopy evaporation w.r.t. canopy air temperature (W m-2 K-1)
 real(rkind)                       :: dLatHeatCanopyEvap_dTCanopy      ! derivative in latent heat of canopy evaporation w.r.t. canopy temperature (W m-2 K-1)
 real(rkind)                       :: dLatHeatCanopyEvap_dTGround      ! derivative in latent heat of canopy evaporation w.r.t. ground temperature (W m-2 K-1)

 ! (liquid water flux derivatives -- ground evap)
 real(rkind)                       :: dLatHeatGroundEvap_dCanLiq       ! derivative in latent heat of ground evaporation w.r.t. canopy liquid water content (J kg-1 s-1)
 real(rkind)                       :: dLatHeatGroundEvap_dTCanair      ! derivative in latent heat of ground evaporation w.r.t. canopy air temperature (W m-2 K-1)
 real(rkind)                       :: dLatHeatGroundEvap_dTCanopy      ! derivative in latent heat of ground evaporation w.r.t. canopy temperature (W m-2 K-1)
 real(rkind)                       :: dLatHeatGroundEvap_dTGround      ! derivative in latent heat of ground evaporation w.r.t. ground temperature (W m-2 K-1)

 ! ---------------------------------------------------------------------------------------
 ! point to variables in the data structure
 ! ---------------------------------------------------------------------------------------
 associate(&

 ! input: model decisions
 ix_bcUpprTdyn                   => model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision,           & ! intent(in): [i4b] choice of upper boundary condition for thermodynamics
 ix_fDerivMeth                   => model_decisions(iLookDECISIONS%fDerivMeth)%iDecision,           & ! intent(in): [i4b] choice of method to compute derivatives
 ix_veg_traits                   => model_decisions(iLookDECISIONS%veg_traits)%iDecision,           & ! intent(in): [i4b] choice of parameterization for vegetation roughness length and displacement height
 ix_canopyEmis                   => model_decisions(iLookDECISIONS%canopyEmis)%iDecision,           & ! intent(in): [i4b] choice of parameterization for canopy emissivity
 ix_windPrfile                   => model_decisions(iLookDECISIONS%windPrfile)%iDecision,           & ! intent(in): [i4b] choice of canopy wind profile
 ix_astability                   => model_decisions(iLookDECISIONS%astability)%iDecision,           & ! intent(in): [i4b] choice of stability function
 ix_soilStress                   => model_decisions(iLookDECISIONS%soilStress)%iDecision,           & ! intent(in): [i4b] choice of function for the soil moisture control on stomatal resistance
 ix_groundwatr                   => model_decisions(iLookDECISIONS%groundwatr)%iDecision,           & ! intent(in): [i4b] choice of groundwater parameterization
 ix_stomResist                   => model_decisions(iLookDECISIONS%stomResist)%iDecision,           & ! intent(in): [i4b] choice of function for stomatal resistance
 ix_spatial_gw                   => model_decisions(iLookDECISIONS%spatial_gw)%iDecision,           & ! intent(in): [i4b] choice of groundwater representation (local, basin)

 ! input: layer geometry
 nSnow                           => indx_data%var(iLookINDEX%nSnow)%dat(1),                         & ! intent(in): [i4b] number of snow layers
 nSoil                           => indx_data%var(iLookINDEX%nSoil)%dat(1),                         & ! intent(in): [i4b] number of soil layers
 nLayers                         => indx_data%var(iLookINDEX%nLayers)%dat(1),                       & ! intent(in): [i4b] total number of layers

 ! input: physical attributes
 vegTypeIndex                    => type_data%var(iLookTYPE%vegTypeIndex),                          & ! intent(in): [i4b] vegetation type index
 soilTypeIndex                   => type_data%var(iLookTYPE%soilTypeIndex),                         & ! intent(in): [i4b] soil type index

 ! input: vegetation parameters
 heightCanopyTop                 => mpar_data%var(iLookPARAM%heightCanopyTop)%dat(1),               & ! intent(in): [dp] height at the top of the vegetation canopy (m)
 heightCanopyBottom              => mpar_data%var(iLookPARAM%heightCanopyBottom)%dat(1),            & ! intent(in): [dp] height at the bottom of the vegetation canopy (m)
 canopyWettingFactor             => mpar_data%var(iLookPARAM%canopyWettingFactor)%dat(1),           & ! intent(in): [dp] maximum wetted fraction of the canopy (-)
 canopyWettingExp                => mpar_data%var(iLookPARAM%canopyWettingExp)%dat(1),              & ! intent(in): [dp] exponent in canopy wetting function (-)
 scalarCanopyIceMax              => diag_data%var(iLookDIAG%scalarCanopyIceMax)%dat(1),             & ! intent(in): [dp] maximum interception storage capacity for ice (kg m-2)
 scalarCanopyLiqMax              => diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1),             & ! intent(in): [dp] maximum interception storage capacity for liquid water (kg m-2)

 ! input: vegetation phenology
 scalarLAI                       => diag_data%var(iLookDIAG%scalarLAI)%dat(1),                      & ! intent(in): [dp] one-sided leaf area index (m2 m-2)
 scalarSAI                       => diag_data%var(iLookDIAG%scalarSAI)%dat(1),                      & ! intent(in): [dp] one-sided stem area index (m2 m-2)
 scalarExposedLAI                => diag_data%var(iLookDIAG%scalarExposedLAI)%dat(1),               & ! intent(in): [dp] exposed leaf area index after burial by snow (m2 m-2)
 scalarExposedSAI                => diag_data%var(iLookDIAG%scalarExposedSAI)%dat(1),               & ! intent(in): [dp] exposed stem area index after burial by snow (m2 m-2)
 scalarGrowingSeasonIndex        => diag_data%var(iLookDIAG%scalarGrowingSeasonIndex)%dat(1),       & ! intent(in): [dp] growing season index (0=off, 1=on)
 scalarFoliageNitrogenFactor     => diag_data%var(iLookDIAG%scalarFoliageNitrogenFactor)%dat(1),    & ! intent(in): [dp] foliage nitrogen concentration (1.0 = saturated)

 ! input: aerodynamic resistance parameters
 z0Snow                          => mpar_data%var(iLookPARAM%z0Snow)%dat(1),                        & ! intent(in): [dp] roughness length of snow (m)
 z0Soil                          => mpar_data%var(iLookPARAM%z0Soil)%dat(1),                        & ! intent(in): [dp] roughness length of soil (m)
 z0CanopyParam                   => mpar_data%var(iLookPARAM%z0Canopy)%dat(1),                      & ! intent(in): [dp] roughness length of the canopy (m)
 zpdFraction                     => mpar_data%var(iLookPARAM%zpdFraction)%dat(1),                   & ! intent(in): [dp] zero plane displacement / canopy height (-)
 critRichNumber                  => mpar_data%var(iLookPARAM%critRichNumber)%dat(1),                & ! intent(in): [dp] critical value for the bulk Richardson number where turbulence ceases (-)
 Louis79_bparam                  => mpar_data%var(iLookPARAM%Louis79_bparam)%dat(1),                & ! intent(in): [dp] parameter in Louis (1979) stability function
 Louis79_cStar                   => mpar_data%var(iLookPARAM%Louis79_cStar)%dat(1),                 & ! intent(in): [dp] parameter in Louis (1979) stability function
 Mahrt87_eScale                  => mpar_data%var(iLookPARAM%Mahrt87_eScale)%dat(1),                & ! intent(in): [dp] exponential scaling factor in the Mahrt (1987) stability function
 windReductionParam              => mpar_data%var(iLookPARAM%windReductionParam)%dat(1),            & ! intent(in): [dp] canopy wind reduction parameter (-)
 leafExchangeCoeff               => mpar_data%var(iLookPARAM%leafExchangeCoeff)%dat(1),             & ! intent(in): [dp] turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
 leafDimension                   => mpar_data%var(iLookPARAM%leafDimension)%dat(1),                 & ! intent(in): [dp] characteristic leaf dimension (m)

 ! input: soil stress parameters
 theta_sat                       => mpar_data%var(iLookPARAM%theta_sat)%dat(1),                     & ! intent(in): [dp] soil porosity (-)
 theta_res                       => mpar_data%var(iLookPARAM%theta_res)%dat(1),                     & ! intent(in): [dp] residual volumetric liquid water content (-)
 plantWiltPsi                    => mpar_data%var(iLookPARAM%plantWiltPsi)%dat(1),                  & ! intent(in): [dp] matric head at wilting point (m)
 soilStressParam                 => mpar_data%var(iLookPARAM%soilStressParam)%dat(1),               & ! intent(in): [dp] parameter in the exponential soil stress function (-)
 critSoilWilting                 => mpar_data%var(iLookPARAM%critSoilWilting)%dat(1),               & ! intent(in): [dp] critical vol. liq. water content when plants are wilting (-)
 critSoilTranspire               => mpar_data%var(iLookPARAM%critSoilTranspire)%dat(1),             & ! intent(in): [dp] critical vol. liq. water content when transpiration is limited (-)
 critAquiferTranspire            => mpar_data%var(iLookPARAM%critAquiferTranspire)%dat(1),          & ! intent(in): [dp] critical aquifer storage value when transpiration is limited (m)
 minStomatalResistance           => mpar_data%var(iLookPARAM%minStomatalResistance)%dat(1),         & ! intent(in): [dp] mimimum stomatal resistance (s m-1)

 ! input: forcing at the upper boundary
 mHeight                         => diag_data%var(iLookDIAG%scalarAdjMeasHeight)%dat(1),            & ! intent(in): [dp] measurement height (m)
 airtemp                         => forc_data%var(iLookFORCE%airtemp),                              & ! intent(in): [dp] air temperature at some height above the surface (K)
 windspd                         => forc_data%var(iLookFORCE%windspd),                              & ! intent(in): [dp] wind speed at some height above the surface (m s-1)
 airpres                         => forc_data%var(iLookFORCE%airpres),                              & ! intent(in): [dp] air pressure at some height above the surface (Pa)
 LWRadAtm                        => forc_data%var(iLookFORCE%LWRadAtm),                             & ! intent(in): [dp] downwelling longwave radiation at the upper boundary (W m-2)
 scalarVPair                     => diag_data%var(iLookDIAG%scalarVPair)%dat(1),                    & ! intent(in): [dp] vapor pressure at some height above the surface (Pa)
 scalarO2air                     => diag_data%var(iLookDIAG%scalarO2air)%dat(1),                    & ! intent(in): [dp] atmospheric o2 concentration (Pa)
 scalarCO2air                    => diag_data%var(iLookDIAG%scalarCO2air)%dat(1),                   & ! intent(in): [dp] atmospheric co2 concentration (Pa)
 scalarTwetbulb                  => diag_data%var(iLookDIAG%scalarTwetbulb)%dat(1),                 & ! intent(in): [dp] wetbulb temperature (K)
 scalarRainfall                  => flux_data%var(iLookFLUX%scalarRainfall)%dat(1),                 & ! intent(in): [dp] computed rainfall rate (kg m-2 s-1)
 scalarSnowfall                  => flux_data%var(iLookFLUX%scalarSnowfall)%dat(1),                 & ! intent(in): [dp] computed snowfall rate (kg m-2 s-1)
 scalarThroughfallRain           => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1),          & ! intent(in): [dp] rainfall through the vegetation canopy (kg m-2 s-1)
 scalarThroughfallSnow           => flux_data%var(iLookFLUX%scalarThroughfallSnow)%dat(1),          & ! intent(in): [dp] snowfall through the vegetation canopy (kg m-2 s-1)

 ! input: water storage
 ! NOTE: soil stress only computed at the start of the substep (firstFluxCall=.true.)
 scalarSWE                       => prog_data%var(iLookPROG%scalarSWE)%dat(1),                      & ! intent(in): [dp]    snow water equivalent on the ground (kg m-2)
 scalarSnowDepth                 => prog_data%var(iLookPROG%scalarSnowDepth)%dat(1),                & ! intent(in): [dp]    snow depth on the ground surface (m)
 mLayerVolFracLiq                => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat,                  & ! intent(in): [dp(:)] volumetric fraction of liquid water in each layer (-)
 mLayerMatricHead                => prog_data%var(iLookPROG%mLayerMatricHead)%dat,                  & ! intent(in): [dp(:)] matric head in each soil layer (m)
 localAquiferStorage             => prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1),           & ! intent(in): [dp]    aquifer storage for the local column (m)
 basinAquiferStorage             => bvar_data%var(iLookBVAR%basin__AquiferStorage)%dat(1),          & ! intent(in): [dp]    aquifer storage for the single basin (m)

 ! input: shortwave radiation fluxes
 scalarCanopySunlitLAI           => diag_data%var(iLookDIAG%scalarCanopySunlitLAI)%dat(1),          & ! intent(in): [dp] sunlit leaf area (-)
 scalarCanopyShadedLAI           => diag_data%var(iLookDIAG%scalarCanopyShadedLAI)%dat(1),          & ! intent(in): [dp] shaded leaf area (-)
 scalarCanopySunlitPAR           => flux_data%var(iLookFLUX%scalarCanopySunlitPAR)%dat(1),          & ! intent(in): [dp] average absorbed par for sunlit leaves (w m-2)
 scalarCanopyShadedPAR           => flux_data%var(iLookFLUX%scalarCanopyShadedPAR)%dat(1),          & ! intent(in): [dp] average absorbed par for shaded leaves (w m-2)
 scalarCanopyAbsorbedSolar       => flux_data%var(iLookFLUX%scalarCanopyAbsorbedSolar)%dat(1),      & ! intent(in): [dp] solar radiation absorbed by canopy (W m-2)
 scalarGroundAbsorbedSolar       => flux_data%var(iLookFLUX%scalarGroundAbsorbedSolar)%dat(1),      & ! intent(in): [dp] solar radiation absorbed by ground (W m-2)

 ! output: fraction of wetted canopy area and fraction of snow on the ground
 scalarCanopyWetFraction         => diag_data%var(iLookDIAG%scalarCanopyWetFraction)%dat(1),        & ! intent(out): [dp] fraction of canopy that is wet
 scalarGroundSnowFraction        => diag_data%var(iLookDIAG%scalarGroundSnowFraction)%dat(1),       & ! intent(out): [dp] fraction of ground covered with snow (-)

 ! output: longwave radiation fluxes
 scalarCanopyEmissivity          => diag_data%var(iLookDIAG%scalarCanopyEmissivity)%dat(1),         & ! intent(out): [dp] effective emissivity of the canopy (-)
 scalarLWRadCanopy               => flux_data%var(iLookFLUX%scalarLWRadCanopy)%dat(1),              & ! intent(out): [dp] longwave radiation emitted from the canopy (W m-2)
 scalarLWRadGround               => flux_data%var(iLookFLUX%scalarLWRadGround)%dat(1),              & ! intent(out): [dp] longwave radiation emitted at the ground surface (W m-2)
 scalarLWRadUbound2Canopy        => flux_data%var(iLookFLUX%scalarLWRadUbound2Canopy)%dat(1),       & ! intent(out): [dp] downward atmospheric longwave radiation absorbed by the canopy (W m-2)
 scalarLWRadUbound2Ground        => flux_data%var(iLookFLUX%scalarLWRadUbound2Ground)%dat(1),       & ! intent(out): [dp] downward atmospheric longwave radiation absorbed by the ground (W m-2)
 scalarLWRadUbound2Ubound        => flux_data%var(iLookFLUX%scalarLWRadUbound2Ubound)%dat(1),       & ! intent(out): [dp] atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)
 scalarLWRadCanopy2Ubound        => flux_data%var(iLookFLUX%scalarLWRadCanopy2Ubound)%dat(1),       & ! intent(out): [dp] longwave radiation emitted from canopy lost thru upper boundary (W m-2)
 scalarLWRadCanopy2Ground        => flux_data%var(iLookFLUX%scalarLWRadCanopy2Ground)%dat(1),       & ! intent(out): [dp] longwave radiation emitted from canopy absorbed by the ground (W m-2)
 scalarLWRadCanopy2Canopy        => flux_data%var(iLookFLUX%scalarLWRadCanopy2Canopy)%dat(1),       & ! intent(out): [dp] canopy longwave reflected from ground and absorbed by the canopy (W m-2)
 scalarLWRadGround2Ubound        => flux_data%var(iLookFLUX%scalarLWRadGround2Ubound)%dat(1),       & ! intent(out): [dp] longwave radiation emitted from ground lost thru upper boundary (W m-2)
 scalarLWRadGround2Canopy        => flux_data%var(iLookFLUX%scalarLWRadGround2Canopy)%dat(1),       & ! intent(out): [dp] longwave radiation emitted from ground and absorbed by the canopy (W m-2)
 scalarLWNetCanopy               => flux_data%var(iLookFLUX%scalarLWNetCanopy)%dat(1),              & ! intent(out): [dp] net longwave radiation at the canopy (W m-2)
 scalarLWNetGround               => flux_data%var(iLookFLUX%scalarLWNetGround)%dat(1),              & ! intent(out): [dp] net longwave radiation at the ground surface (W m-2)
 scalarLWNetUbound               => flux_data%var(iLookFLUX%scalarLWNetUbound)%dat(1),              & ! intent(out): [dp] net longwave radiation at the upper boundary (W m-2)

 ! output: aerodynamic resistance
 scalarZ0Canopy                  => diag_data%var(iLookDIAG%scalarZ0Canopy)%dat(1),                 & ! intent(out): [dp] roughness length of the canopy (m)
 scalarWindReductionFactor       => diag_data%var(iLookDIAG%scalarWindReductionFactor)%dat(1),      & ! intent(out): [dp] canopy wind reduction factor (-)
 scalarZeroPlaneDisplacement     => diag_data%var(iLookDIAG%scalarZeroPlaneDisplacement)%dat(1),    & ! intent(out): [dp] zero plane displacement (m)
 scalarRiBulkCanopy              => diag_data%var(iLookDIAG%scalarRiBulkCanopy)%dat(1),             & ! intent(out): [dp] bulk Richardson number for the canopy (-)
 scalarRiBulkGround              => diag_data%var(iLookDIAG%scalarRiBulkGround)%dat(1),             & ! intent(out): [dp] bulk Richardson number for the ground surface (-)
 scalarEddyDiffusCanopyTop       => flux_data%var(iLookFLUX%scalarEddyDiffusCanopyTop)%dat(1),      & ! intent(out): [dp] eddy diffusivity for heat at the top of the canopy (m2 s-1)
 scalarFrictionVelocity          => flux_data%var(iLookFLUX%scalarFrictionVelocity)%dat(1),         & ! intent(out): [dp] friction velocity (m s-1)
 scalarWindspdCanopyTop          => flux_data%var(iLookFLUX%scalarWindspdCanopyTop)%dat(1),         & ! intent(out): [dp] windspeed at the top of the canopy (m s-1)
 scalarWindspdCanopyBottom       => flux_data%var(iLookFLUX%scalarWindspdCanopyBottom)%dat(1),      & ! intent(out): [dp] windspeed at the height of the bottom of the canopy (m s-1)
 scalarLeafResistance            => flux_data%var(iLookFLUX%scalarLeafResistance)%dat(1),           & ! intent(out): [dp] mean leaf boundary layer resistance per unit leaf area (s m-1)
 scalarGroundResistance          => flux_data%var(iLookFLUX%scalarGroundResistance)%dat(1),         & ! intent(out): [dp] below canopy aerodynamic resistance (s m-1)
 scalarCanopyResistance          => flux_data%var(iLookFLUX%scalarCanopyResistance)%dat(1),         & ! intent(out): [dp] above canopy aerodynamic resistance (s m-1)

 ! input/output: soil resistance -- intent(in) and intent(inout) because only called at the first flux call
 mLayerRootDensity               => diag_data%var(iLookDIAG%mLayerRootDensity)%dat,                 & ! intent(in):    [dp] root density in each layer (-)
 scalarAquiferRootFrac           => diag_data%var(iLookDIAG%scalarAquiferRootFrac)%dat(1),          & ! intent(in):    [dp] fraction of roots below the lowest soil layer (-)
 scalarTranspireLim              => diag_data%var(iLookDIAG%scalarTranspireLim)%dat(1),             & ! intent(inout): [dp] weighted average of the transpiration limiting factor (-)
 mLayerTranspireLim              => diag_data%var(iLookDIAG%mLayerTranspireLim)%dat,                & ! intent(inout): [dp] transpiration limiting factor in each layer (-)
 scalarTranspireLimAqfr          => diag_data%var(iLookDIAG%scalarTranspireLimAqfr)%dat(1),         & ! intent(inout): [dp] transpiration limiting factor for the aquifer (-)
 scalarSoilRelHumidity           => diag_data%var(iLookDIAG%scalarSoilRelHumidity)%dat(1),          & ! intent(inout): [dp] relative humidity in the soil pores [0-1]
 scalarSoilResistance            => flux_data%var(iLookFLUX%scalarSoilResistance)%dat(1),           & ! intent(inout): [dp] resistance from the soil (s m-1)

 ! input/output: stomatal resistance -- intent(inout) because only called at the first flux call
 scalarStomResistSunlit          => flux_data%var(iLookFLUX%scalarStomResistSunlit)%dat(1),         & ! intent(inout): [dp] stomatal resistance for sunlit leaves (s m-1)
 scalarStomResistShaded          => flux_data%var(iLookFLUX%scalarStomResistShaded)%dat(1),         & ! intent(inout): [dp] stomatal resistance for shaded leaves (s m-1)
 scalarPhotosynthesisSunlit      => flux_data%var(iLookFLUX%scalarPhotosynthesisSunlit)%dat(1),     & ! intent(inout): [dp] sunlit photosynthesis (umolco2 m-2 s-1)
 scalarPhotosynthesisShaded      => flux_data%var(iLookFLUX%scalarPhotosynthesisShaded)%dat(1),     & ! intent(inout): [dp] shaded photosynthesis (umolco2 m-2 s-1)

 ! output: turbulent heat fluxes
 scalarLatHeatSubVapCanopy       => diag_data%var(iLookDIAG%scalarLatHeatSubVapCanopy)%dat(1),      & ! intent(inout): [dp] latent heat of sublimation/vaporization for the vegetation canopy (J kg-1)
 scalarLatHeatSubVapGround       => diag_data%var(iLookDIAG%scalarLatHeatSubVapGround)%dat(1),      & ! intent(inout): [dp] latent heat of sublimation/vaporization for the ground surface (J kg-1)
 scalarSatVP_canopyTemp          => diag_data%var(iLookDIAG%scalarSatVP_CanopyTemp)%dat(1),         & ! intent(out):   [dp] saturation vapor pressure at the temperature of the vegetation canopy (Pa)
 scalarSatVP_groundTemp          => diag_data%var(iLookDIAG%scalarSatVP_GroundTemp)%dat(1),         & ! intent(out):   [dp] saturation vapor pressure at the temperature of the ground surface (Pa)
 scalarSenHeatTotal              => flux_data%var(iLookFLUX%scalarSenHeatTotal)%dat(1),             & ! intent(out):   [dp] sensible heat from the canopy air space to the atmosphere (W m-2)
 scalarSenHeatCanopy             => flux_data%var(iLookFLUX%scalarSenHeatCanopy)%dat(1),            & ! intent(out):   [dp] sensible heat flux from the canopy to the canopy air space (W m-2)
 scalarSenHeatGround             => flux_data%var(iLookFLUX%scalarSenHeatGround)%dat(1),            & ! intent(out):   [dp] sensible heat flux from ground surface below vegetation (W m-2)
 scalarLatHeatTotal              => flux_data%var(iLookFLUX%scalarLatHeatTotal)%dat(1),             & ! intent(out):   [dp] latent heat from the canopy air space to the atmosphere (W m-2)
 scalarLatHeatCanopyEvap         => flux_data%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1),        & ! intent(out):   [dp] latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
 scalarLatHeatCanopyTrans        => flux_data%var(iLookFLUX%scalarLatHeatCanopyTrans)%dat(1),       & ! intent(out):   [dp] latent heat flux for transpiration from the canopy to the canopy air space (W m-2)
 scalarLatHeatGround             => flux_data%var(iLookFLUX%scalarLatHeatGround)%dat(1),            & ! intent(out):   [dp] latent heat flux from ground surface below vegetation (W m-2)

 ! output: advective heat fluxes
 scalarCanopyAdvectiveHeatFlux   => flux_data%var(iLookFLUX%scalarCanopyAdvectiveHeatFlux)%dat(1),  & ! intent(out): [dp] heat advected to the canopy surface with rain + snow (W m-2)
 scalarGroundAdvectiveHeatFlux   => flux_data%var(iLookFLUX%scalarGroundAdvectiveHeatFlux)%dat(1),  & ! intent(out): [dp] heat advected to the ground surface with throughfall (W m-2)

 ! output: mass fluxes
 scalarCanopySublimation         => flux_data%var(iLookFLUX%scalarCanopySublimation)%dat(1),        & ! intent(out): [dp] canopy sublimation/frost (kg m-2 s-1)
 scalarSnowSublimation           => flux_data%var(iLookFLUX%scalarSnowSublimation)%dat(1),          & ! intent(out): [dp] snow sublimation/frost -- below canopy or non-vegetated (kg m-2 s-1)

 ! input/output: canopy air space variables
 scalarVP_CanopyAir              => diag_data%var(iLookDIAG%scalarVP_CanopyAir)%dat(1),             & ! intent(inout): [dp] vapor pressure of the canopy air space (Pa)
 scalarCanopyStabilityCorrection => diag_data%var(iLookDIAG%scalarCanopyStabilityCorrection)%dat(1),& ! intent(inout): [dp] stability correction for the canopy (-)
 scalarGroundStabilityCorrection => diag_data%var(iLookDIAG%scalarGroundStabilityCorrection)%dat(1),& ! intent(inout): [dp] stability correction for the ground surface (-)

 ! output: liquid water fluxes
 scalarCanopyTranspiration       => flux_data%var(iLookFLUX%scalarCanopyTranspiration)%dat(1),      & ! intent(out): [dp] canopy transpiration (kg m-2 s-1)
 scalarCanopyEvaporation         => flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1),        & ! intent(out): [dp] canopy evaporation/condensation (kg m-2 s-1)
 scalarGroundEvaporation         => flux_data%var(iLookFLUX%scalarGroundEvaporation)%dat(1),        & ! intent(out): [dp] ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)

 ! output: derived fluxes
 scalarTotalET                   => flux_data%var(iLookFLUX%scalarTotalET)%dat(1),                  & ! intent(out): [dp] total ET (kg m-2 s-1)
 scalarNetRadiation              => flux_data%var(iLookFLUX%scalarNetRadiation)%dat(1)              & ! intent(out): [dp] net radiation (W m-2)
 )
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="vegNrgFlux/"

 ! initialize printflag
 printflag = .false.

 ! identify the type of boundary condition for thermodynamics
 select case(ix_bcUpprTdyn)

  ! *****
  ! (1) DIRICHLET OR ZERO FLUX BOUNDARY CONDITION...
  ! ************************************************

  ! NOTE: Vegetation fluxes are not computed in this case

  ! ** prescribed temperature or zero flux at the upper boundary of the snow-soil system
  case(prescribedTemp,zeroFlux)

   ! derived fluxes
   scalarTotalET             = 0._rkind    ! total ET (kg m-2 s-1)
   scalarNetRadiation        = 0._rkind    ! net radiation (W m-2)
   ! liquid water fluxes associated with evaporation/transpiration
   scalarCanopyTranspiration = 0._rkind    ! canopy transpiration (kg m-2 s-1)
   scalarCanopyEvaporation   = 0._rkind    ! canopy evaporation/condensation (kg m-2 s-1)
   scalarGroundEvaporation   = 0._rkind    ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
   ! solid water fluxes associated with sublimation/frost
   scalarCanopySublimation   = 0._rkind    ! sublimation from the vegetation canopy ((kg m-2 s-1)
   scalarSnowSublimation     = 0._rkind    ! sublimation from the snow surface ((kg m-2 s-1)
   ! set canopy fluxes to zero (no canopy)
   canairNetFlux             = 0._rkind    ! net energy flux for the canopy air space (W m-2)
   canopyNetFlux             = 0._rkind    ! net energy flux for the vegetation canopy (W m-2)
   ! set canopy derivatives to zero
   dCanairNetFlux_dCanairTemp = 0._rkind   ! derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
   dCanairNetFlux_dCanopyTemp = 0._rkind   ! derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
   dCanairNetFlux_dGroundTemp = 0._rkind   ! derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
   dCanopyNetFlux_dCanairTemp = 0._rkind   ! derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
   dCanopyNetFlux_dCanopyTemp = 0._rkind   ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
   dCanopyNetFlux_dGroundTemp = 0._rkind   ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
   dGroundNetFlux_dCanairTemp = 0._rkind   ! derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
   dGroundNetFlux_dCanopyTemp = 0._rkind   ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
   ! set liquid flux derivatives to zero (canopy evap)
   dCanopyEvaporation_dCanLiq = 0._rkind    ! derivative in canopy evaporation w.r.t. canopy liquid water content (s-1)
   dCanopyEvaporation_dTCanair= 0._rkind    ! derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   dCanopyEvaporation_dTCanopy= 0._rkind    ! derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
   dCanopyEvaporation_dTGround= 0._rkind    ! derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
   ! set liquid flux derivatives to zero (ground evap)
   dGroundEvaporation_dCanLiq = 0._rkind    ! derivative in ground evaporation w.r.t. canopy liquid water content (s-1)
   dGroundEvaporation_dTCanair= 0._rkind    ! derivative in ground evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   dGroundEvaporation_dTCanopy= 0._rkind    ! derivative in ground evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
   dGroundEvaporation_dTGround= 0._rkind    ! derivative in ground evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)

   ! compute fluxes and derivatives -- separate approach for prescribed temperature and zero flux
   if(ix_bcUpprTdyn == prescribedTemp)then
    ! compute ground net flux (W m-2)
    groundNetFlux = -diag_data%var(iLookDIAG%iLayerThermalC)%dat(0)*(groundTempTrial - upperBoundTemp)/(prog_data%var(iLookPROG%mLayerDepth)%dat(1)*0.5_rkind)
    ! compute derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
    dGroundNetFlux_dGroundTemp = -diag_data%var(iLookDIAG%iLayerThermalC)%dat(0)/(prog_data%var(iLookPROG%mLayerDepth)%dat(1)*0.5_rkind)
   elseif(model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision == zeroFlux)then
    groundNetFlux              = 0._rkind
    dGroundNetFlux_dGroundTemp = 0._rkind
   else
    err=20; message=trim(message)//'unable to identify upper boundary condition for thermodynamics: expect the case to be prescribedTemp or zeroFlux'; return
   end if

  ! *****
  ! (2) NEUMANN BOUNDARY CONDITION...
  ! *********************************

  ! NOTE 1: This is the main routine for calculating vegetation fluxes
  ! NOTE 2: This routine also calculates surface fluxes for the case where vegetation is buried with snow (or bare soil)

  ! *******************************************************************************************************************************************************************
  ! *******************************************************************************************************************************************************************
  ! ***** PRELIMINARIES  **********************************************************************************************************************************************
  ! *******************************************************************************************************************************************************************
  ! *******************************************************************************************************************************************************************

  ! * flux boundary condition
  case(energyFlux)

   ! identify the appropriate groundwater variable
   select case(ix_spatial_gw)
    case(singleBasin); scalarAquiferStorage = basinAquiferStorage
    case(localColumn); scalarAquiferStorage = localAquiferStorage
    case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
   end select ! (modify the groundwater representation for this single-column implementation)

   ! set canopy stability corrections to the previous values
   scalarCanopyStabilityCorrection_old = scalarCanopyStabilityCorrection       ! stability correction for the canopy (-)
   scalarGroundStabilityCorrection_old = scalarGroundStabilityCorrection       ! stability correction for the ground surface (-)

   ! initialize variables to compute stomatal resistance
   if(firstFluxCall .and. firstSubStep)then
    ! vapor pressure in the canopy air space initialized as vapor pressure of air above the vegetation canopy
    ! NOTE: this is needed for the stomatal resistance calculations
    if(scalarVP_CanopyAir < 0._rkind)then
     scalarVP_CanopyAir    = scalarVPair - 1._rkind    ! "small" offset used to assist in checking initial derivative calculations
    end if
   end if

   ! set latent heat of sublimation/vaporization for canopy and ground surface (Pa/K)
   ! NOTE: variables are constant over the substep, to simplify relating energy and mass fluxes
   if(firstFluxCall)then
    scalarLatHeatSubVapCanopy = getLatentHeatValue(canopyTempTrial)
    ! case when there is snow on the ground (EXCLUDE "snow without a layer" -- in this case, evaporate from the soil)
    if(nSnow > 0)then
     if(groundTempTrial > Tfreeze)then; err=20; message=trim(message)//'do not expect ground temperature > 0 when snow is on the ground'; return; end if
     scalarLatHeatSubVapGround = LH_sub  ! sublimation from snow
     scalarGroundSnowFraction  = 1._rkind
    ! case when the ground is snow-free
    else
     scalarLatHeatSubVapGround = LH_vap  ! evaporation of water in the soil pores: this occurs even if frozen because of super-cooled water
     scalarGroundSnowFraction  = 0._rkind
    end if  ! (if there is snow on the ground)
   end if  ! (if the first flux call)
   !write(*,'(a,1x,10(f30.10,1x))') 'groundTempTrial, scalarLatHeatSubVapGround = ', groundTempTrial, scalarLatHeatSubVapGround

   ! compute the roughness length of the ground (ground below the canopy or non-vegetated surface)
   z0Ground = z0soil*(1._rkind - scalarGroundSnowFraction) + z0Snow*scalarGroundSnowFraction     ! roughness length (m)

   ! compute the total vegetation area index (leaf plus stem)
   VAI        = scalarLAI + scalarSAI  ! vegetation area index
   exposedVAI = scalarExposedLAI + scalarExposedSAI  !  exposed vegetation area index

   ! compute emissivity of the canopy (-)
   if(computeVegFlux)then
    select case(ix_canopyEmis)
     ! *** simple exponential function
     case(simplExp)
      scalarCanopyEmissivity = 1._rkind - exp(-exposedVAI)                                     ! effective emissivity of the canopy (-)
     ! *** canopy emissivity parameterized as a function of diffuse transmissivity
     case(difTrans)
      ! compute the exponential integral
      scaleLAI = 0.5_rkind*exposedVAI
      expi     = expInt(scaleLAI)
      ! compute diffuse transmissivity (-)
      diffuseTrans = (1._rkind - scaleLAI)*exp(-scaleLAI) + (scaleLAI**2._rkind)*expi
      ! compute the canopy emissivity
      scalarCanopyEmissivity = (1._rkind - diffuseTrans)*vegEmissivity
     ! *** check we found the correct option
     case default
      err=20; message=trim(message)//'unable to identify option for canopy emissivity'; return
    end select
   end if

   ! ensure canopy longwave fluxes are zero when not computing canopy fluxes
   if(.not.computeVegFlux) scalarCanopyEmissivity=0._rkind

   ! compute emissivity of the ground surface (-)
   groundEmissivity = scalarGroundSnowFraction*snowEmissivity + (1._rkind - scalarGroundSnowFraction)*soilEmissivity  ! emissivity of the ground surface (-)

   ! compute the fraction of canopy that is wet
   ! NOTE: we either sublimate or evaporate over the entire substep
   if(computeVegFlux)then

    ! compute the fraction of liquid water in the canopy (-)
    totalCanopyWater = canopyLiqTrial + canopyIceTrial
    if(totalCanopyWater > tiny(1.0_rkind))then
     fracLiquidCanopy = canopyLiqTrial / (canopyLiqTrial + canopyIceTrial)
    else
     fracLiquidCanopy = 0._rkind
    end if

    ! get wetted fraction and derivatives
    call wettedFrac(&
                    ! input
                    .true.,                                         & ! flag to denote if derivative is desired
                    (ix_fDerivMeth == numerical),                   & ! flag to denote that numerical derivatives are required (otherwise, analytical derivatives are calculated)
                    (scalarLatHeatSubVapCanopy > LH_vap+verySmall), & ! flag to denote if the canopy is frozen
                    dCanLiq_dTcanopy,                               & ! derivative in canopy liquid w.r.t. canopy temperature (kg m-2 K-1)
                    fracLiquidCanopy,                               & ! fraction of liquid water on the canopy (-)
                    canopyLiqTrial,                                 & ! canopy liquid water (kg m-2)
                    canopyIceTrial,                                 & ! canopy ice (kg m-2)
                    scalarCanopyLiqMax,                             & ! maximum canopy liquid water (kg m-2)
                    scalarCanopyIceMax,                             & ! maximum canopy ice content (kg m-2)
                    canopyWettingFactor,                            & ! maximum wetted fraction of the canopy (-)
                    canopyWettingExp,                               & ! exponent in canopy wetting function (-)
                    ! output
                    scalarCanopyWetFraction,                        & ! canopy wetted fraction (-)
                    dCanopyWetFraction_dWat,                        & ! derivative in wetted fraction w.r.t. canopy total water (kg-1 m2)
                    dCanopyWetFraction_dT,                          & ! derivative in wetted fraction w.r.t. canopy temperature (K-1)
                    err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

   else
    scalarCanopyWetFraction = 0._rkind  ! canopy wetted fraction (-)
    dCanopyWetFraction_dWat = 0._rkind  ! derivative in wetted fraction w.r.t. canopy liquid water (kg-1 m2)
    dCanopyWetFraction_dT   = 0._rkind  ! derivative in wetted fraction w.r.t. canopy temperature (K-1)
   end if
   !write(*,'(a,1x,L1,1x,f25.15,1x))') 'computeVegFlux, scalarCanopyWetFraction = ', computeVegFlux, scalarCanopyWetFraction
   !print*, 'dCanopyWetFraction_dWat = ', dCanopyWetFraction_dWat
   !print*, 'dCanopyWetFraction_dT   = ', dCanopyWetFraction_dT
   !print*, 'canopyLiqTrial = ', canopyLiqTrial
   !print*, 'canopyIceTrial = ', canopyIceTrial
   !print*, 'scalarCanopyLiqMax = ', scalarCanopyLiqMax
   !print*, 'scalarCanopyIceMax = ', scalarCanopyIceMax

   ! *******************************************************************************************************************************************************************
   ! *******************************************************************************************************************************************************************
   ! ***** AERODYNAMIC RESISTANCE *****************************************************************************************************************************************
   ! *******************************************************************************************************************************************************************
   ! *******************************************************************************************************************************************************************

   ! NOTE: compute for all iterations

   ! compute aerodynamic resistances
   ! Refs: Choudhury and Monteith (4-layer model for heat budget of homogenous surfaces; QJRMS, 1988)
   !       Niu and Yang (Canopy effects on snow processes; JGR, 2004)
   !       Mahat et al. (Below-canopy turbulence in a snowmelt model, WRR, 2012)
   call aeroResist(&
                   ! input: model control
                   computeVegFlux,                     & ! intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                   (ix_fDerivMeth == analytical),      & ! intent(in): logical flag if would like to compute analytical derivaties
                   ix_veg_traits,                      & ! intent(in): choice of parameterization for vegetation roughness length and displacement height
                   ix_windPrfile,                      & ! intent(in): choice of canopy wind profile
                   ix_astability,                      & ! intent(in): choice of stability function
                   ! input: above-canopy forcing data
                   mHeight,                            & ! intent(in): measurement height (m)
                   airtemp,                            & ! intent(in): air temperature at some height above the surface (K)
                   windspd,                            & ! intent(in): wind speed at some height above the surface (m s-1)
                   ! input: canopy and ground temperature
                   canairTempTrial,                    & ! intent(in): temperature of the canopy air space (K)
                   groundTempTrial,                    & ! intent(in): temperature of the ground surface (K)
                   ! input: diagnostic variables
                   exposedVAI,                         & ! intent(in): exposed vegetation area index -- leaf plus stem (m2 m-2)
                   scalarSnowDepth,                    & ! intent(in): snow depth (m)
                   ! input: parameters
                   z0Ground,                           & ! intent(in): roughness length of the ground (below canopy or non-vegetated surface [snow]) (m)
                   z0CanopyParam,                      & ! intent(in): roughness length of the canopy (m)
                   zpdFraction,                        & ! intent(in): zero plane displacement / canopy height (-)
                   critRichNumber,                     & ! intent(in): critical value for the bulk Richardson number where turbulence ceases (-)
                   Louis79_bparam,                     & ! intent(in): parameter in Louis (1979) stability function
                   Mahrt87_eScale,                     & ! intent(in): exponential scaling factor in the Mahrt (1987) stability function
                   windReductionParam,                 & ! intent(in): canopy wind reduction parameter (-)
                   leafExchangeCoeff,                  & ! intent(in): turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
                   leafDimension,                      & ! intent(in): characteristic leaf dimension (m)
                   heightCanopyTop,                    & ! intent(in): height at the top of the vegetation canopy (m)
                   heightCanopyBottom,                 & ! intent(in): height at the bottom of the vegetation canopy (m)
                   ! output: stability corrections
                   scalarRiBulkCanopy,                 & ! intent(out): bulk Richardson number for the canopy (-)
                   scalarRiBulkGround,                 & ! intent(out): bulk Richardson number for the ground surface (-)
                   scalarCanopyStabilityCorrection,    & ! intent(out): stability correction for the canopy (-)
                   scalarGroundStabilityCorrection,    & ! intent(out): stability correction for the ground surface (-)
                   ! output: scalar resistances
                   scalarZ0Canopy,                     & ! intent(out): roughness length of the canopy (m)
                   scalarWindReductionFactor,          & ! intent(out): canopy wind reduction factor (-)
                   scalarZeroPlaneDisplacement,        & ! intent(out): zero plane displacement (m)
                   scalarEddyDiffusCanopyTop,          & ! intent(out): eddy diffusivity for heat at the top of the canopy (m2 s-1)
                   scalarFrictionVelocity,             & ! intent(out): friction velocity (m s-1)
                   scalarWindspdCanopyTop,             & ! intent(out): windspeed at the top of the canopy (m s-1)
                   scalarWindspdCanopyBottom,          & ! intent(out): windspeed at the height of the bottom of the canopy (m s-1)
                   scalarLeafResistance,               & ! intent(out): mean leaf boundary layer resistance per unit leaf area (s m-1)
                   scalarGroundResistance,             & ! intent(out): below canopy aerodynamic resistance (s m-1)
                   scalarCanopyResistance,             & ! intent(out): above canopy aerodynamic resistance (s m-1)
                   ! output: derivatives in scalar resistances
                   dGroundResistance_dTGround,         & ! intent(out): derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
                   dGroundResistance_dTCanopy,         & ! intent(out): derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
                   dGroundResistance_dTCanair,         & ! intent(out): derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
                   dCanopyResistance_dTCanopy,         & ! intent(out): derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
                   dCanopyResistance_dTCanair,         & ! intent(out): derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
                   ! output: error control
                   err,cmessage                        ) ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
   !print*,         scalarLeafResistance,    & ! mean leaf boundary layer resistance per unit leaf area (s m-1)
   !                scalarGroundResistance,  & ! below canopy aerodynamic resistance (s m-1)
   !                scalarCanopyResistance,  & ! above canopy aerodynamic resistance (s m-1)
   !                '(leaf, ground, canopy)'

   ! *******************************************************************************************************************************************************************
   ! *******************************************************************************************************************************************************************
   ! ***** STOMATAL RESISTANCE *****************************************************************************************************************************************
   ! *******************************************************************************************************************************************************************
   ! *******************************************************************************************************************************************************************

   ! stomatal resistance is constant over the SUBSTEP
   ! NOTE: This is a simplification, as stomatal resistance does depend on canopy temperature
   !       This "short-cut" made because:
   !         (1) computations are expensive;
   !         (2) derivative calculations are rather complex (iterations within the Ball-Berry routine); and
   !         (3) stomatal resistance does not change rapidly
   if(firstFluxCall)then

    ! compute the saturation vapor pressure for vegetation temperature
    TV_celcius = canopyTempTrial - Tfreeze
    call satVapPress(TV_celcius, scalarSatVP_CanopyTemp, dSVPCanopy_dCanopyTemp)

    ! compute soil moisture factor controlling stomatal resistance
    call soilResist(&
                    ! input (model decisions)
                    ix_soilStress,                     & ! intent(in): choice of function for the soil moisture control on stomatal resistance
                    ix_groundwatr,                     & ! intent(in): groundwater parameterization
                    ! input (state variables)
                    mLayerMatricHead(1:nSoil),         & ! intent(in): matric head in each soil layer (m)
                    mLayerVolFracLiq(nSnow+1:nLayers), & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                    scalarAquiferStorage,              & ! intent(in): aquifer storage (m)
                    ! input (diagnostic variables)
                    mLayerRootDensity(1:nSoil),        & ! intent(in): root density in each layer (-)
                    scalarAquiferRootFrac,             & ! intent(in): fraction of roots below the lowest soil layer (-)
                    ! input (parameters)
                    plantWiltPsi,                      & ! intent(in): matric head at wilting point (m)
                    soilStressParam,                   & ! intent(in): parameter in the exponential soil stress function (-)
                    critSoilWilting,                   & ! intent(in): critical vol. liq. water content when plants are wilting (-)
                    critSoilTranspire,                 & ! intent(in): critical vol. liq. water content when transpiration is limited (-)
                    critAquiferTranspire,              & ! intent(in): critical aquifer storage value when transpiration is limited (m)
                    ! output
                    scalarTranspireLim,                & ! intent(out): weighted average of the transpiration limiting factor (-)
                    mLayerTranspireLim(1:nSoil),       & ! intent(out): transpiration limiting factor in each layer (-)
                    scalarTranspireLimAqfr,            & ! intent(out): transpiration limiting factor for the aquifer (-)
                    err,cmessage                       ) ! intent(out): error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    !print*, 'weighted average of the soil moiture factor controlling stomatal resistance (-) = ', scalarTranspireLim

    !write(*,'(a,1x,10(f20.10,1x))') 'canopyTempTrial, scalarSatVP_CanopyTemp, scalarVP_CanopyAir = ', &
    !                                 canopyTempTrial, scalarSatVP_CanopyTemp, scalarVP_CanopyAir

    ! compute stomatal resistance
    call stomResist(&
                    ! input (state and diagnostic variables)
                    canopyTempTrial,                   & ! intent(in): temperature of the vegetation canopy (K)
                    scalarSatVP_CanopyTemp,            & ! intent(in): saturation vapor pressure at the temperature of the veg canopy (Pa)
                    scalarVP_CanopyAir,                & ! intent(in): canopy air vapor pressure (Pa)
                    ! input: data structures
                    type_data,                         & ! intent(in):    type of vegetation and soil
                    forc_data,                         & ! intent(in):    model forcing data
                    mpar_data,                         & ! intent(in):    model parameters
                    model_decisions,                   & ! intent(in):    model decisions
                    ! input-output: data structures
                    diag_data,                         & ! intent(inout): model diagnostic variables for a local HRU
                    flux_data,                         & ! intent(inout): model fluxes for a local HRU
                    ! output: error control
                    err,cmessage                       ) ! intent(out): error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

   end if  ! (if the first flux call in a given sub-step)


   ! *******************************************************************************************************************************************************************
   ! *******************************************************************************************************************************************************************
   ! ***** LONGWAVE RADIATION  *****************************************************************************************************************************************
   ! *******************************************************************************************************************************************************************
   ! *******************************************************************************************************************************************************************

   ! compute canopy longwave radiation balance
   call longwaveBal(&
                    ! input: model control
                    ix_fDerivMeth,                     & ! intent(in): method used to calculate flux derivatives
                    computeVegFlux,                    & ! intent(in): flag to compute fluxes over vegetation
                    ! input: canopy and ground temperature
                    canopyTempTrial,                   & ! intent(in): temperature of the vegetation canopy (K)
                    groundTempTrial,                   & ! intent(in): temperature of the ground surface (K)
                    ! input: canopy and ground emissivity
                    scalarCanopyEmissivity,            & ! intent(in): canopy emissivity (-)
                    groundEmissivity,                  & ! intent(in): ground emissivity (-)
                    ! input: forcing
                    LWRadAtm,                          & ! intent(in): downwelling longwave radiation at the upper boundary (W m-2)
                    ! output: emitted radiation from the canopy and ground
                    scalarLWRadCanopy,                 & ! intent(out): longwave radiation emitted from the canopy (W m-2)
                    scalarLWRadGround,                 & ! intent(out): longwave radiation emitted at the ground surface (W m-2)
                    ! output: individual fluxes
                    scalarLWRadUbound2Canopy,          & ! intent(out): downward atmospheric longwave radiation absorbed by the canopy (W m-2)
                    scalarLWRadUbound2Ground,          & ! intent(out): downward atmospheric longwave radiation absorbed by the ground (W m-2)
                    scalarLWRadUbound2Ubound,          & ! intent(out): atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)
                    scalarLWRadCanopy2Ubound,          & ! intent(out): longwave radiation emitted from canopy lost thru upper boundary (W m-2)
                    scalarLWRadCanopy2Ground,          & ! intent(out): longwave radiation emitted from canopy absorbed by the ground (W m-2)
                    scalarLWRadCanopy2Canopy,          & ! intent(out): canopy longwave reflected from ground and absorbed by the canopy (W m-2)
                    scalarLWRadGround2Ubound,          & ! intent(out): longwave radiation emitted from ground lost thru upper boundary (W m-2)
                    scalarLWRadGround2Canopy,          & ! intent(out): longwave radiation emitted from ground and absorbed by the canopy (W m-2)
                    ! output: net fluxes
                    scalarLWNetCanopy,                 & ! intent(out): net longwave radiation at the canopy (W m-2)
                    scalarLWNetGround,                 & ! intent(out): net longwave radiation at the ground surface (W m-2)
                    scalarLWNetUbound,                 & ! intent(out): net longwave radiation at the upper boundary (W m-2)
                    ! output: flux derivatives
                    dLWNetCanopy_dTCanopy,             & ! intent(out): derivative in net canopy radiation w.r.t. canopy temperature (W m-2 K-1)
                    dLWNetGround_dTGround,             & ! intent(out): derivative in net ground radiation w.r.t. ground temperature (W m-2 K-1)
                    dLWNetCanopy_dTGround,             & ! intent(out): derivative in net canopy radiation w.r.t. ground temperature (W m-2 K-1)
                    dLWNetGround_dTCanopy,             & ! intent(out): derivative in net ground radiation w.r.t. canopy temperature (W m-2 K-1)
                    ! output: error control
                    err,cmessage                       ) ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
   !print*, 'dLWNetCanopy_dTGround = ', dLWNetCanopy_dTGround


   ! *******************************************************************************************************************************************************************
   ! *******************************************************************************************************************************************************************
   ! ***** TURBULENT HEAT FLUXES  **************************************************************************************************************************************
   ! *******************************************************************************************************************************************************************
   ! *******************************************************************************************************************************************************************

   ! check the need to compute numerical derivatives
   if(ix_fDerivMeth == numerical)then
    nFlux=5  ! compute the derivatives using one-sided finite differences
   else
    nFlux=1  ! compute analytical derivatives
   end if

   ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
   do itry=nFlux,1,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)

    ! -------------------------------------------------------------------------------------
    ! state perturbations for numerical deriavtives with one-sided finite differences
    ! note: no perturbations performed using analytical derivatives (nFlux=1)
    ! -------------------------------------------------------------------------------------

    ! identify the type of perturbation
    select case(itry)

     ! un-perturbed case
     case(unperturbed)
      groundTemp        = groundTempTrial
      canopyTemp        = canopyTempTrial
      canairTemp        = canairTempTrial
      canopyWetFraction = scalarCanopyWetFraction

     ! perturb ground temperature
     case(perturbStateGround)
      groundTemp        = groundTempTrial + dx
      canopyTemp        = canopyTempTrial
      canairTemp        = canairTempTrial
      canopyWetFraction = scalarCanopyWetFraction

     ! perturb canopy temperature
     case(perturbStateCanopy)
      groundTemp        = groundTempTrial
      canopyTemp        = canopyTempTrial + dx
      canairTemp        = canairTempTrial
      canopyWetFraction = scalarCanopyWetFraction

     ! perturb canopy air temperature
     case(perturbStateCanair)
      groundTemp        = groundTempTrial
      canopyTemp        = canopyTempTrial
      canairTemp        = canairTempTrial + dx
      canopyWetFraction = scalarCanopyWetFraction

     ! perturb canopy liquid water content
     case(perturbStateCanLiq)
      groundTemp        = groundTempTrial
      canopyTemp        = canopyTempTrial
      canairTemp        = canairTempTrial

      ! perturbations in canopy liquid water content affect canopy wetted fraction
      if(computeVegFlux)then
       call wettedFrac(&
                       ! input
                       .false.,                                        & ! flag to denote if derivative is desired
                       (ix_fDerivMeth == numerical),                   & ! flag to denote that numerical derivatives are required (otherwise, analytical derivatives are calculated)
                       (scalarLatHeatSubVapCanopy > LH_vap+verySmall), & ! flag to denote if the canopy is frozen
                       dCanLiq_dTcanopy,                               & ! derivative in canopy liquid w.r.t. canopy temperature (kg m-2 K-1)
                       fracLiquidCanopy,                               & ! fraction of liquid water on the canopy (-)
                       canopyLiqTrial+dx,                              & ! canopy liquid water (kg m-2)
                       canopyIceTrial,                                 & ! canopy ice (kg m-2)
                       scalarCanopyLiqMax,                             & ! maximum canopy liquid water (kg m-2)
                       scalarCanopyIceMax,                             & ! maximum canopy ice content (kg m-2)
                       canopyWettingFactor,                            & ! maximum wetted fraction of the canopy (-)
                       canopyWettingExp,                               & ! exponent in canopy wetting function (-)
                       ! output
                       canopyWetFraction,                              & ! canopy wetted fraction (-)
                       dCanopyWetFraction_dWat,                        & ! derivative in wetted fraction w.r.t. canopy liquid water (kg-1 m2)
                       dCanopyWetFraction_dT,                          & ! derivative in wetted fraction w.r.t. canopy temperature (K-1)
                       err,cmessage)
       if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

      else
       canopyWetFraction = 0._rkind
      end if
      !print*, 'wetted fraction derivative = ', (canopyWetFraction - scalarCanopyWetFraction)/dx
      !pause

     ! check for an unknown perturbation
     case default; err=10; message=trim(message)//"unknown perturbation"; return

    end select ! (type of perturbation)

    ! compute the saturation vapor pressure for vegetation temperature
    ! NOTE: saturated vapor pressure derivatives don't seem that accurate....
    TV_celcius = canopyTemp - Tfreeze
    call satVapPress(TV_celcius, scalarSatVP_CanopyTemp, dSVPCanopy_dCanopyTemp)

    ! compute the saturation vapor pressure for ground temperature
    ! NOTE: saturated vapor pressure derivatives don't seem that accurate....
    TG_celcius = groundTemp - Tfreeze
    call satVapPress(TG_celcius, scalarSatVP_GroundTemp, dSVPGround_dGroundTemp)

    ! -------------------------------------------------------------------------------------
    ! calculation block (unperturbed fluxes returned [computed last])
    ! -------------------------------------------------------------------------------------

    ! re-compute aerodynamic resistances for perturbed cases
    ! NOTE: unperturbed fluxes computed earlier, and not over-written
    if(itry /= unperturbed)then
     call aeroResist(&
                     ! input: model control
                     computeVegFlux,                          & ! intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                     .false.,                                 & ! intent(in): logical flag if would like to compute analytical derivaties
                     ix_veg_traits,                           & ! intent(in): choice of parameterization for vegetation roughness length and displacement height
                     ix_windPrfile,                           & ! intent(in): choice of canopy wind profile
                     ix_astability,                           & ! intent(in): choice of stability function
                     ! input: above-canopy forcing data
                     mHeight,                                 & ! intent(in): measurement height (m)
                     airtemp,                                 & ! intent(in): air temperature at some height above the surface (K)
                     windspd,                                 & ! intent(in): wind speed at some height above the surface (m s-1)
                     ! input: temperature (canopy, ground, canopy air space)
                     canairTemp,                              & ! intent(in): temperature of the canopy air space (K)
                     groundTemp,                              & ! intent(in): ground temperature (K)
                     ! input: diagnostic variables
                     exposedVAI,                              & ! intent(in): exposed vegetation area index -- leaf plus stem (m2 m-2)
                     scalarSnowDepth,                         & ! intent(in): snow depth (m)
                     ! input: parameters
                     z0Ground,                                & ! intent(in): roughness length of the ground (below canopy or non-vegetated surface [snow]) (m)
                     z0CanopyParam,                           & ! intent(in): roughness length of the canopy (m)
                     zpdFraction,                             & ! intent(in): zero plane displacement / canopy height (-)
                     critRichNumber,                          & ! intent(in): critical value for the bulk Richardson number where turbulence ceases (-)
                     Louis79_bparam,                          & ! intent(in): parameter in Louis (1979) stability function
                     Mahrt87_eScale,                          & ! intent(in): exponential scaling factor in the Mahrt (1987) stability function
                     windReductionParam,                      & ! intent(in): canopy wind reduction parameter (-)
                     leafExchangeCoeff,                       & ! intent(in): turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
                     leafDimension,                           & ! intent(in): characteristic leaf dimension (m)
                     heightCanopyTop,                         & ! intent(in): height at the top of the vegetation canopy (m)
                     heightCanopyBottom,                      & ! intent(in): height at the bottom of the vegetation canopy (m)
                     ! output: stability corrections
                     notUsed_RiBulkCanopy,                    & ! intent(out): bulk Richardson number for the canopy (-)
                     notUsed_RiBulkGround,                    & ! intent(out): bulk Richardson number for the ground surface (-)
                     notUsed_scalarCanopyStabilityCorrection, & ! intent(out): stability correction for the canopy (-)
                     notUsed_scalarGroundStabilityCorrection, & ! intent(out): stability correction for the ground surface (-)
                     ! output: scalar resistances
                     notUsed_z0Canopy,                        & ! intent(out): roughness length of the canopy (m)
                     notUsed_WindReductionFactor,             & ! intent(out): canopy wind reduction factor (-)
                     notUsed_ZeroPlaneDisplacement,           & ! intent(out): zero plane displacement (m)
                     notUsed_EddyDiffusCanopyTop,             & ! intent(out): eddy diffusivity for heat at the top of the canopy (m2 s-1)
                     notUsed_FrictionVelocity,                & ! intent(out): friction velocity (m s-1)
                     notUsed_WindspdCanopyTop,                & ! intent(out): windspeed at the top of the canopy (m s-1)
                     notUsed_WindspdCanopyBottom,             & ! intent(out): windspeed at the height of the bottom of the canopy (m s-1)
                     trialLeafResistance,                     & ! intent(out): mean leaf boundary layer resistance per unit leaf area (s m-1)
                     trialGroundResistance,                   & ! intent(out): below canopy aerodynamic resistance (s m-1)
                     trialCanopyResistance,                   & ! intent(out): above canopy aerodynamic resistance (s m-1)
                     ! output: derivatives in scalar resistances
                     notUsed_dGroundResistance_dTGround,      & ! intent(out): derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
                     notUsed_dGroundResistance_dTCanopy,      & ! intent(out): derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
                     notUsed_dGroundResistance_dTCanair,      & ! intent(out): derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
                     notUsed_dCanopyResistance_dTCanopy,      & ! intent(out): derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
                     notUsed_dCanopyResistance_dTCanair,      & ! intent(out): derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
                     ! output: error control
                     err,cmessage                             ) ! intent(out): error control
     if(err/=0)then; message=trim(message)//trim(cmessage); return; end if


    ! assign scalar resistances for un-perturbed cases
    else
     trialLeafResistance   = scalarLeafResistance
     trialGroundResistance = scalarGroundResistance
     trialCanopyResistance = scalarCanopyResistance

    end if  ! (re-computing resistances for perturbed cases)
    !print*, 'trialLeafResistance = ', trialLeafResistance
    !print*, 'trialGroundResistance = ', trialGroundResistance
    !print*, 'trialCanopyResistance = ', trialCanopyResistance

    ! compute the relative humidity in the top soil layer and the resistance at the ground surface
    ! NOTE: computations are based on start-of-step values, so only compute for the first flux call
    if(firstFluxCall)then
     ! (soil water evaporation factor [0-1])
     soilEvapFactor = mLayerVolFracLiq(nSnow+1)/(theta_sat - theta_res)
     ! (resistance from the soil [s m-1])
     scalarSoilResistance = scalarGroundSnowFraction*1._rkind + (1._rkind - scalarGroundSnowFraction)*EXP(8.25_rkind - 4.225_rkind*soilEvapFactor)  ! Sellers (1992)
     !scalarSoilResistance = scalarGroundSnowFraction*0._rkind + (1._rkind - scalarGroundSnowFraction)*exp(8.25_rkind - 6.0_rkind*soilEvapFactor)    ! Niu adjustment to decrease resitance for wet soil
     ! (relative humidity in the soil pores [0-1])
     if(mLayerMatricHead(1) > -1.e+6_rkind)then  ! avoid problems with numerical precision when soil is very dry
      soilRelHumidity_noSnow = exp( (mLayerMatricHead(1)*gravity) / (groundTemp*R_wv) )
     else
      soilRelHumidity_noSnow = 0._rkind
     end if ! (if matric head is very low)
     scalarSoilRelHumidity  = scalarGroundSnowFraction*1._rkind + (1._rkind - scalarGroundSnowFraction)*soilRelHumidity_noSnow
     !print*, 'mLayerMatricHead(1), scalarSoilRelHumidity = ', mLayerMatricHead(1), scalarSoilRelHumidity
    end if  ! (if the first flux call)

    ! compute turbulent heat fluxes
    call turbFluxes(&
                    ! input: model control
                    computeVegFlux,                       & ! intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                    ix_fDerivMeth,                        & ! intent(in): method used to calculate flux derivatives
                    ! input: above-canopy forcing data
                    airtemp,                              & ! intent(in): air temperature at some height above the surface (K)
                    airpres,                              & ! intent(in): air pressure of the air above the vegetation canopy (Pa)
                    scalarVPair,                          & ! intent(in): vapor pressure of the air above the vegetation canopy (Pa)
                    ! input: latent heat of sublimation/vaporization
                    scalarLatHeatSubVapCanopy,            & ! intent(in): latent heat of sublimation/vaporization for the vegetation canopy (J kg-1)
                    scalarLatHeatSubVapGround,            & ! intent(in): latent heat of sublimation/vaporization for the ground surface (J kg-1)
                    ! input: canopy/ground temperature and saturated vapor pressure
                    canairTemp,                           & ! intent(in): temperature of the canopy air space (K)
                    canopyTemp,                           & ! intent(in): canopy temperature (K)
                    groundTemp,                           & ! intent(in): ground temperature (K)
                    scalarSatVP_CanopyTemp,               & ! intent(in): saturation vapor pressure at the temperature of the veg canopy (Pa)
                    scalarSatVP_GroundTemp,               & ! intent(in): saturation vapor pressure at the temperature of the ground (Pa)
                    dSVPCanopy_dCanopyTemp,               & ! intent(in): derivative in canopy saturation vapor pressure w.r.t. canopy temperature (Pa K-1)
                    dSVPGround_dGroundTemp,               & ! intent(in): derivative in ground saturation vapor pressure w.r.t. ground temperature (Pa K-1)
                    ! input: diagnostic variables
                    exposedVAI,                           & ! intent(in): exposed vegetation area index -- leaf plus stem (m2 m-2)
                    canopyWetFraction,                    & ! intent(in): trial value for the fraction of canopy that is wet [0-1]
                    dCanopyWetFraction_dWat,              & ! intent(in): derivative in the canopy wetted fraction w.r.t. total water content (kg-1 m-2)
                    dCanopyWetFraction_dT,                & ! intent(in): derivative in wetted fraction w.r.t. canopy temperature (K-1)
                    scalarCanopySunlitLAI,                & ! intent(in): sunlit leaf area (-)
                    scalarCanopyShadedLAI,                & ! intent(in): shaded leaf area (-)
                    scalarSoilRelHumidity,                & ! intent(in): relative humidity in the soil pores [0-1]
                    scalarSoilResistance,                 & ! intent(in): resistance from the soil (s m-1)
                    trialLeafResistance,                  & ! intent(in): mean leaf boundary layer resistance per unit leaf area (s m-1)
                    trialGroundResistance,                & ! intent(in): below canopy aerodynamic resistance (s m-1)
                    trialCanopyResistance,                & ! intent(in): above canopy aerodynamic resistance (s m-1)
                    scalarStomResistSunlit,               & ! intent(in): stomatal resistance for sunlit leaves (s m-1)
                    scalarStomResistShaded,               & ! intent(in): stomatal resistance for shaded leaves (s m-1)
                    ! input: derivatives in scalar resistances
                    dGroundResistance_dTGround,           & ! intent(in): derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
                    dGroundResistance_dTCanopy,           & ! intent(in): derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
                    dGroundResistance_dTCanair,           & ! intent(in): derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
                    dCanopyResistance_dTCanopy,           & ! intent(in): derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
                    dCanopyResistance_dTCanair,           & ! intent(in): derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
                    ! output: conductances (used to check derivative calculations)
                    scalarLeafConductance,                & ! intent(out): leaf conductance (m s-1)
                    scalarCanopyConductance,              & ! intent(out): canopy conductance (m s-1)
                    scalarGroundConductanceSH,            & ! intent(out): ground conductance for sensible heat (m s-1)
                    scalarGroundConductanceLH,            & ! intent(out): ground conductance for latent heat -- includes soil resistance (m s-1)
                    scalarEvapConductance,                & ! intent(out): conductance for evaporation (m s-1)
                    scalarTransConductance,               & ! intent(out): conductance for transpiration (m s-1)
                    scalarTotalConductanceSH,             & ! intent(out): total conductance for sensible heat (m s-1)
                    scalarTotalConductanceLH,             & ! intent(out): total conductance for latent heat (m s-1)
                    ! output: canopy air space variables
                    scalarVP_CanopyAir,                   & ! intent(out): vapor pressure of the canopy air space (Pa)
                    ! output: fluxes from the vegetation canopy
                    scalarSenHeatCanopy,                  & ! intent(out): sensible heat flux from the canopy to the canopy air space (W m-2)
                    scalarLatHeatCanopyEvap,              & ! intent(out): latent heat flux associated with evaporation from the canopy to the canopy air space (W m-2)
                    scalarLatHeatCanopyTrans,             & ! intent(out): latent heat flux associated with transpiration from the canopy to the canopy air space (W m-2)
                    ! output: fluxes from non-vegetated surfaces (ground surface below vegetation, bare ground, or snow covered vegetation)
                    scalarSenHeatGround,                  & ! intent(out): sensible heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
                    scalarLatHeatGround,                  & ! intent(out): latent heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
                    ! output: total heat fluxes to the atmosphere
                    scalarSenHeatTotal,                   & ! intent(out): total sensible heat flux to the atmosphere (W m-2)
                    scalarLatHeatTotal,                   & ! intent(out): total latent heat flux to the atmosphere (W m-2)
                    ! output: net fluxes
                    turbFluxCanair,                       & ! intent(out): net turbulent heat fluxes at the canopy air space (W m-2)
                    turbFluxCanopy,                       & ! intent(out): net turbulent heat fluxes at the canopy (W m-2)
                    turbFluxGround,                       & ! intent(out): net turbulent heat fluxes at the ground surface (W m-2)
                    ! output: energy flux derivatives
                    dTurbFluxCanair_dTCanair,             & ! intent(out): derivative in net canopy air space fluxes w.r.t. canopy air temperature (W m-2 K-1)
                    dTurbFluxCanair_dTCanopy,             & ! intent(out): derivative in net canopy air space fluxes w.r.t. canopy temperature (W m-2 K-1)
                    dTurbFluxCanair_dTGround,             & ! intent(out): derivative in net canopy air space fluxes w.r.t. ground temperature (W m-2 K-1)
                    dTurbFluxCanopy_dTCanair,             & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
                    dTurbFluxCanopy_dTCanopy,             & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
                    dTurbFluxCanopy_dTGround,             & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
                    dTurbFluxGround_dTCanair,             & ! intent(out): derivative in net ground turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
                    dTurbFluxGround_dTCanopy,             & ! intent(out): derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
                    dTurbFluxGround_dTGround,             & ! intent(out): derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
                    ! output: liquid flux derivatives (canopy evap)
                    dLatHeatCanopyEvap_dCanLiq,           & ! intent(out): derivative in latent heat of canopy evaporation w.r.t. canopy liquid water content (W kg-1)
                    dLatHeatCanopyEvap_dTCanair,          & ! intent(out): derivative in latent heat of canopy evaporation w.r.t. canopy air temperature (W m-2 K-1)
                    dLatHeatCanopyEvap_dTCanopy,          & ! intent(out): derivative in latent heat of canopy evaporation w.r.t. canopy temperature (W m-2 K-1)
                    dLatHeatCanopyEvap_dTGround,          & ! intent(out): derivative in latent heat of canopy evaporation w.r.t. ground temperature (W m-2 K-1)
                    ! output: liquid flux derivatives (ground evap)
                    dLatHeatGroundEvap_dCanLiq,           & ! intent(out): derivative in latent heat of ground evaporation w.r.t. canopy liquid water content (J kg-1 s-1)
                    dLatHeatGroundEvap_dTCanair,          & ! intent(out): derivative in latent heat of ground evaporation w.r.t. canopy air temperature
                    dLatHeatGroundEvap_dTCanopy,          & ! intent(out): derivative in latent heat of ground evaporation w.r.t. canopy temperature
                    dLatHeatGroundEvap_dTGround,          & ! intent(out): derivative in latent heat of ground evaporation w.r.t. ground temperature
                    ! output: cross derivatives
                    dTurbFluxCanair_dCanLiq,              & ! intent(out): derivative in net canopy air space fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                    dTurbFluxCanopy_dCanLiq,              & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                    dTurbFluxGround_dCanLiq,              & ! intent(out): derivative in net ground turbulent fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                    ! output: error control
                    err,cmessage                          ) ! intent(out): error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

    !write(*,'(a,f25.15)') 'scalarSenHeatTotal = ', scalarSenHeatTotal
    !write(*,'(a,f25.15)') 'scalarSenHeatCanopy = ', scalarSenHeatCanopy
    !write(*,'(a,f25.15)') 'scalarLatHeatCanopyEvap = ', scalarLatHeatCanopyEvap
    !write(*,'(a,f25.15)') 'scalarLatHeatCanopyTrans = ', scalarLatHeatCanopyTrans

    !print*, 'scalarSenHeatGround = ', scalarSenHeatGround
    !print*, 'scalarLatHeatGround = ', scalarLatHeatGround

    !notUsed_scalarCanopyStabilityCorrection  ! stability correction for the canopy (-)
    !notUsed_scalarGroundStabilityCorrection  ! stability correction for the ground surface (-)
    !notUsed_EddyDiffusCanopyTop              ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
    !notUsed_FrictionVelocity                 ! friction velocity (m s-1)
    !notUsed_WindspdCanopyTop                 ! windspeed at the top of the canopy (m s-1)
    !notUsed_WindspdCanopyBottom              ! windspeed at the height of the bottom of the canopy (m s-1)
    !trialLeafResistance                      ! mean leaf boundary layer resistance per unit leaf area (s m-1)
    !trialGroundResistance                    ! below canopy aerodynamic resistance (s m-1)
    !trialCanopyResistance                    ! above canopy aerodynamic resistance (s m-1)

    ! save perturbed fluxes
    if(ix_fDerivMeth == numerical)then
     select case(itry) ! (select type of perturbation)
      case(unperturbed)
       try0 = turbFluxGround
       exit
      case(perturbStateCanair)
       turbFluxCanair_dStateCanair = turbFluxCanair          ! turbulent exchange from the canopy air space to the atmosphere (W m-2)
       turbFluxCanopy_dStateCanair = turbFluxCanopy          ! total turbulent heat fluxes from the canopy to the canopy air space (W m-2)
       turbFluxGround_dStateCanair = turbFluxGround          ! total turbulent heat fluxes from the ground to the canopy air space (W m-2)
       latHeatCanEvap_dStateCanair = scalarLatHeatCanopyEvap ! perturbed value for the latent heat associated with canopy evaporation (W m-2)
      case(perturbStateCanopy)
       turbFluxCanair_dStateCanopy = turbFluxCanair          ! turbulent exchange from the canopy air space to the atmosphere (W m-2)
       turbFluxCanopy_dStateCanopy = turbFluxCanopy          ! total turbulent heat fluxes from the canopy to the canopy air space (W m-2)
       turbFluxGround_dStateCanopy = turbFluxGround          ! total turbulent heat fluxes from the ground to the canopy air space (W m-2)
       latHeatCanEvap_dStateCanopy = scalarLatHeatCanopyEvap ! perturbed value for the latent heat associated with canopy evaporation (W m-2)
      case(perturbStateGround)
       try1 = turbFluxGround
       turbFluxCanair_dStateGround = turbFluxCanair          ! turbulent exchange from the canopy air space to the atmosphere (W m-2)
       turbFluxCanopy_dStateGround = turbFluxCanopy          ! total turbulent heat fluxes from the canopy to the canopy air space (W m-2)
       turbFluxGround_dStateGround = turbFluxGround          ! total turbulent heat fluxes from the ground to the canopy air space (W m-2)
       latHeatCanEvap_dStateGround = scalarLatHeatCanopyEvap ! perturbed value for the latent heat associated with canopy evaporation (W m-2)
      case(perturbStateCanLiq)
       turbFluxCanair_dStateCanliq = turbFluxCanair          ! turbulent exchange from the canopy air space to the atmosphere (W m-2)
       turbFluxCanopy_dStateCanLiq = turbFluxCanopy          ! total turbulent heat fluxes from the canopy to the canopy air space (W m-2)
       turbFluxGround_dStateCanLiq = turbFluxGround          ! total turbulent heat fluxes from the ground to the canopy air space (W m-2)
       latHeatCanEvap_dStateCanliq = scalarLatHeatCanopyEvap ! perturbed value for the latent heat associated with canopy evaporation (W m-2)
      case default; err=10; message=trim(message)//"unknown perturbation"; return
     end select ! (type of perturbation)
    end if ! (if numerical)

   end do  ! (looping through different flux perturbations)

   ! test derivative
   !if(ix_fDerivMeth == numerical)  print*, 'try0, try1 = ', try0, try1
   !if(ix_fDerivMeth == numerical)  print*, 'derivative = ', (ix_fDerivMeth == numerical), (try1 - try0)/dx
   !if(ix_fDerivMeth == analytical) print*, 'derivative = ', (ix_fDerivMeth == numerical), dTurbFluxGround_dTGround
   !pause

   ! compute numerical derivatives
   if(ix_fDerivMeth == numerical)then
    ! derivatives w.r.t. canopy air temperature
    dTurbFluxCanair_dTCanair    = (turbFluxCanair_dStateCanair - turbFluxCanair) / dx          ! derivative in net canopy air space fluxes w.r.t. canopy air temperature (W m-2 K-1)
    dTurbFluxCanopy_dTCanair    = (turbFluxCanopy_dStateCanair - turbFluxCanopy) / dx          ! derivative in net canopy turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
    dTurbFluxGround_dTCanair    = (turbFluxGround_dStateCanair - turbFluxGround) / dx          ! derivative in net ground turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
    dLatHeatCanopyEvap_dTCanair = (latHeatCanEvap_dStateCanair - scalarLatHeatCanopyEvap) / dx ! derivative in latent heat of canopy evaporation w.r.t. canopy air temperature (W m-2 K-1)
    ! derivatives w.r.t. canopy temperature
    dTurbFluxCanair_dTCanopy    = (turbFluxCanair_dStateCanopy - turbFluxCanair) / dx          ! derivative in net canopy air space fluxes w.r.t. canopy temperature (W m-2 K-1)
    dTurbFluxCanopy_dTCanopy    = (turbFluxCanopy_dStateCanopy - turbFluxCanopy) / dx          ! derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
    dTurbFluxGround_dTCanopy    = (turbFluxGround_dStateCanopy - turbFluxGround) / dx          ! derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
    dLatHeatCanopyEvap_dTCanopy = (latHeatCanEvap_dStateCanopy - scalarLatHeatCanopyEvap) / dx ! derivative in latent heat of canopy evaporation w.r.t. canopy temperature (W m-2 K-1)
    ! derivatives w.r.t. ground temperature
    dTurbFluxCanair_dTGround    = (turbFluxCanair_dStateGround - turbFluxCanair) / dx          ! derivative in net canopy air space fluxes w.r.t. ground temperature (W m-2 K-1)
    dTurbFluxCanopy_dTGround    = (turbFluxCanopy_dStateGround - turbFluxCanopy) / dx          ! derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
    dTurbFluxGround_dTGround    = (turbFluxGround_dStateGround - turbFluxGround) / dx          ! derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
    dLatHeatCanopyEvap_dTGround = (latHeatCanEvap_dStateGround - scalarLatHeatCanopyEvap) / dx ! derivative in latent heat of canopy evaporation w.r.t. ground temperature (W m-2 K-1)
    ! derivatives w.r.t. canopy liquid water content
    dTurbFluxCanair_dCanLiq    = (turbFluxCanair_dStateCanliq  - turbFluxCanair) / dx          ! derivative in net canopy air space fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
    dTurbFluxCanopy_dCanLiq    = (turbFluxCanopy_dStateCanLiq  - turbFluxCanopy) / dx          ! derivative in net canopy turbulent fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
    dTurbFluxGround_dCanLiq    = (turbFluxGround_dStateCanLiq  - turbFluxGround) / dx          ! derivative in net ground turbulent fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
    dLatHeatCanopyEvap_dCanLiq = (latHeatCanEvap_dStateCanliq  - scalarLatHeatCanopyEvap) / dx ! derivative in latent heat of canopy evaporation w.r.t. canopy liquid water content (J kg-1 s-1)
   end if
   !if(heightCanopyBottom < scalarSnowDepth+z0Ground) pause 'bottom of the canopy is covered'

   ! test
   !print*, (ix_fDerivMeth == numerical)
   !print*, 'dTurbFluxCanair_dTCanair = ', dTurbFluxCanair_dTCanair
   !print*, 'dTurbFluxCanair_dTCanopy = ', dTurbFluxCanair_dTCanopy
   !print*, 'dTurbFluxCanair_dTGround = ', dTurbFluxCanair_dTGround
   !print*, 'dTurbFluxCanopy_dTCanair = ', dTurbFluxCanopy_dTCanair
   !print*, 'dTurbFluxCanopy_dTCanopy = ', dTurbFluxCanopy_dTCanopy
   !print*, 'dTurbFluxCanopy_dTGround = ', dTurbFluxCanopy_dTGround
   !print*, 'dTurbFluxGround_dTCanair = ', dTurbFluxGround_dTCanair
   !print*, 'dTurbFluxGround_dTCanopy = ', dTurbFluxGround_dTCanopy
   !print*, 'dTurbFluxGround_dTGround = ', dTurbFluxGround_dTGround
   !print*, 'dLatHeatCanopyEvap_dCanLiq = ', dLatHeatCanopyEvap_dCanLiq
   !print*, 'dLatHeatCanopyEvap_dTCanair = ', dLatHeatCanopyEvap_dTCanair
   !print*, 'dLatHeatCanopyEvap_dTCanopy = ', dLatHeatCanopyEvap_dTCanopy
   !print*, 'dLatHeatCanopyEvap_dTGround = ', dLatHeatCanopyEvap_dTGround
   !print*, 'dTurbFluxCanair_dCanLiq = ', dTurbFluxCanair_dCanLiq
   !print*, 'dTurbFluxCanopy_dCanLiq = ', dTurbFluxCanopy_dCanLiq
   !print*, 'dTurbFluxGround_dCanLiq = ', dTurbFluxGround_dCanLiq
   !print*, '*****'
   !pause

   !print*, 'scalarRainfall, scalarThroughfallRain, scalarSnowfall, scalarThroughfallSnow, canopyTempTrial, scalarTwetbulb = ', &
   !         scalarRainfall, scalarThroughfallRain, scalarSnowfall, scalarThroughfallSnow, canopyTempTrial, scalarTwetbulb

   ! compute the heat advected with precipitation (W m-2)
   ! NOTE: fluxes are in kg m-2 s-1, so no need to use density of water/ice here
   scalarCanopyAdvectiveHeatFlux = -Cp_water*(scalarRainfall - scalarThroughfallRain)*(canopyTempTrial - scalarTwetbulb) + &
                                   (-Cp_ice)*(scalarSnowfall - scalarThroughfallSnow)*(canopyTempTrial - scalarTwetbulb)
   scalarGroundAdvectiveHeatFlux = -Cp_water*scalarThroughfallRain*(groundTempTrial - scalarTwetbulb)         + &
                                   (-Cp_ice)*scalarThroughfallSnow*(groundTempTrial - scalarTwetbulb)         !+ &
   !                                -Cp_water*scalarCanopyLiqDrainage  *(groundTempTrial - canopyTempTrial) + &
   !                                -Cp_ice  *scalarCanopySnowUnloading*(groundTempTrial - canopyTempTrial)
   !print*, 'scalarRainfall, scalarThroughfallRain, scalarSnowfall, scalarThroughfallSnow = ', scalarRainfall, scalarThroughfallRain, scalarSnowfall, scalarThroughfallSnow
   !print*, 'scalarCanopyAdvectiveHeatFlux, scalarGroundAdvectiveHeatFlux = ', scalarCanopyAdvectiveHeatFlux, scalarGroundAdvectiveHeatFlux

   ! compute the mass flux associated with transpiration and evaporation/sublimation (J m-2 s-1 --> kg m-2 s-1)
   ! NOTE: remove water from the snow on the ground in preference to removing water from the water in soil pores
   !print*, 'scalarLatHeatCanopyTrans = ', scalarLatHeatCanopyTrans
   !print*, 'scalarLatHeatGround      = ', scalarLatHeatGround
   ! (canopy transpiration/sublimation)
   if(scalarLatHeatSubVapCanopy > LH_vap+verySmall)then ! sublimation
    scalarCanopyEvaporation = 0._rkind
    scalarCanopySublimation = scalarLatHeatCanopyEvap/LH_sub
    if(scalarLatHeatCanopyTrans > 0._rkind)then ! flux directed towards the veg
     scalarCanopySublimation   = scalarCanopySublimation + scalarLatHeatCanopyTrans/LH_sub ! frost
     scalarCanopyTranspiration = 0._rkind
    else
     scalarCanopyTranspiration = scalarLatHeatCanopyTrans/LH_vap  ! transpiration is always vapor
    end if
   ! (canopy transpiration/evaporation)
   else                                                 ! evaporation
    scalarCanopyEvaporation = scalarLatHeatCanopyEvap/LH_vap
    scalarCanopySublimation = 0._rkind
    if(scalarLatHeatCanopyTrans > 0._rkind)then ! flux directed towards the veg
     scalarCanopyEvaporation   = scalarCanopyEvaporation + scalarLatHeatCanopyTrans/LH_vap
     scalarCanopyTranspiration = 0._rkind
    else
     scalarCanopyTranspiration = scalarLatHeatCanopyTrans/LH_vap
    end if
   end if
   ! (ground evaporation/sublimation)
   if(scalarLatHeatSubVapGround > LH_vap+verySmall)then ! sublimation
    ! NOTE: this should only occur when we have formed snow layers, so check
    if(nSnow == 0)then; err=20; message=trim(message)//'only expect snow sublimation when we have formed some snow layers'; return; end if
    scalarGroundEvaporation = 0._rkind  ! ground evaporation is zero once the snowpack has formed
    scalarSnowSublimation   = scalarLatHeatGround/LH_sub
   else
    ! NOTE: this should only occur when we have no snow layers, so check
    if(nSnow > 0)then; err=20; message=trim(message)//'only expect ground evaporation when there are no snow layers'; return; end if
    scalarGroundEvaporation = scalarLatHeatGround/LH_vap
    scalarSnowSublimation   = 0._rkind  ! no sublimation from snow if no snow layers have formed
   end if
   !print*, 'scalarSnowSublimation, scalarLatHeatGround = ', scalarSnowSublimation, scalarLatHeatGround

   !print*, 'canopyWetFraction, scalarCanopyEvaporation = ', canopyWetFraction, scalarCanopyEvaporation

   ! *******************************************************************************************************************************************************************
   ! *******************************************************************************************************************************************************************
   ! ***** AND STITCH EVERYTHING TOGETHER  *****************************************************************************************************************************
   ! *******************************************************************************************************************************************************************
   ! *******************************************************************************************************************************************************************

   ! compute derived fluxes
   scalarTotalET      = scalarGroundEvaporation + scalarCanopyEvaporation + scalarCanopyTranspiration
   scalarNetRadiation = scalarCanopyAbsorbedSolar + scalarLWNetCanopy + scalarGroundAbsorbedSolar + scalarLWNetGround

   ! compute net fluxes at the canopy and ground surface
   canairNetFlux = turbFluxCanair
   canopyNetFlux = scalarCanopyAbsorbedSolar + scalarLWNetCanopy + turbFluxCanopy + scalarCanopyAdvectiveHeatFlux
   groundNetFlux = scalarGroundAbsorbedSolar + scalarLWNetGround + turbFluxGround + scalarGroundAdvectiveHeatFlux
   !write(*,'(a,1x,10(e17.10,1x))') 'canopyNetFlux, groundNetFlux,  scalarLWNetCanopy, turbFluxCanopy, turbFluxGround, scalarLWNetGround, scalarCanopyAdvectiveHeatFlux = ', &
   !                                 canopyNetFlux, groundNetFlux,  scalarLWNetCanopy, turbFluxCanopy, turbFluxGround, scalarLWNetGround, scalarCanopyAdvectiveHeatFlux
   !write(*,'(a,1x,10(e20.14,1x))') 'groundNetFlux, scalarGroundAbsorbedSolar,  scalarLWNetGround, turbFluxGround, scalarGroundAdvectiveHeatFlux = ', &
   !                                 groundNetFlux, scalarGroundAbsorbedSolar,  scalarLWNetGround, turbFluxGround, scalarGroundAdvectiveHeatFlux

   ! compute the energy derivatives
   dCanairNetFlux_dCanairTemp = dTurbFluxCanair_dTCanair
   dCanairNetFlux_dCanopyTemp = dTurbFluxCanair_dTCanopy
   dCanairNetFlux_dGroundTemp = dTurbFluxCanair_dTGround
   dCanopyNetFlux_dCanairTemp = dTurbFluxCanopy_dTCanair
   dCanopyNetFlux_dCanopyTemp = dLWNetCanopy_dTCanopy + dTurbFluxCanopy_dTCanopy - Cp_water*(scalarRainfall - scalarThroughfallRain) - Cp_ice*(scalarSnowfall - scalarThroughfallSnow)
   dCanopyNetFlux_dGroundTemp = dLWNetCanopy_dTGround + dTurbFluxCanopy_dTGround
   dGroundNetFlux_dCanairTemp = dTurbFluxGround_dTCanair
   dGroundNetFlux_dCanopyTemp = dLWNetGround_dTCanopy + dTurbFluxGround_dTCanopy
   dGroundNetFlux_dGroundTemp = dLWNetGround_dTGround + dTurbFluxGround_dTGround - Cp_water*scalarThroughfallRain - Cp_ice*scalarThroughfallSnow

   ! check if evaporation or sublimation
   if(scalarLatHeatSubVapCanopy < LH_vap+verySmall)then ! evaporation

    ! compute the liquid water derivarives
    dCanopyEvaporation_dCanLiq  = dLatHeatCanopyEvap_dCanLiq/LH_vap    ! (s-1)
    dCanopyEvaporation_dTCanair = dLatHeatCanopyEvap_dTCanair/LH_vap   ! (kg m-2 s-1 K-1)
    dCanopyEvaporation_dTCanopy = dLatHeatCanopyEvap_dTCanopy/LH_vap   ! (kg m-2 s-1 K-1)
    dCanopyEvaporation_dTGround = dLatHeatCanopyEvap_dTGround/LH_vap   ! (kg m-2 s-1 K-1)

   ! sublimation
   else
    dCanopyEvaporation_dCanLiq  = 0._rkind  ! (s-1)
    dCanopyEvaporation_dTCanair = 0._rkind  ! (kg m-2 s-1 K-1)
    dCanopyEvaporation_dTCanopy = 0._rkind  ! (kg m-2 s-1 K-1)
    dCanopyEvaporation_dTGround = 0._rkind  ! (kg m-2 s-1 K-1)
   end if

   ! compute the liquid water derivarives (ground evap)
   dGroundEvaporation_dCanLiq  = dLatHeatGroundEvap_dCanLiq/LH_vap    ! (s-1)
   dGroundEvaporation_dTCanair = dLatHeatGroundEvap_dTCanair/LH_vap   ! (kg m-2 s-1 K-1)
   dGroundEvaporation_dTCanopy = dLatHeatGroundEvap_dTCanopy/LH_vap   ! (kg m-2 s-1 K-1)
   dGroundEvaporation_dTGround = dLatHeatGroundEvap_dTGround/LH_vap   ! (kg m-2 s-1 K-1)

   ! compute the cross derivative terms (only related to turbulent fluxes; specifically canopy evaporation and transpiration)
   dCanopyNetFlux_dCanLiq = dTurbFluxCanopy_dCanLiq  ! derivative in net canopy fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
   dGroundNetFlux_dCanLiq = dTurbFluxGround_dCanLiq  ! derivative in net ground fluxes w.r.t. canopy liquid water content (J kg-1 s-1)

   !print*, (ix_fDerivMeth == numerical)
   !print*, 'dGroundNetFlux_dCanairTemp = ', dGroundNetFlux_dCanairTemp
   !print*, 'dCanopyNetFlux_dCanopyTemp = ', dCanopyNetFlux_dCanopyTemp
   !print*, 'dGroundNetFlux_dCanopyTemp = ', dGroundNetFlux_dCanopyTemp
   !print*, 'dCanopyNetFlux_dGroundTemp = ', dCanopyNetFlux_dGroundTemp
   !print*, 'dGroundNetFlux_dGroundTemp = ', dGroundNetFlux_dGroundTemp
   !print*, 'dLWNetCanopy_dTGround = ', dLWNetCanopy_dTGround
   !print*, 'dTurbFluxCanopy_dTGround = ', dTurbFluxCanopy_dTGround
   !pause

  ! * check
  case default; err=10; message=trim(message)//'unable to identify upper boundary condition for thermodynamics'; return

 ! end case statement
 end select ! upper boundary condition for thermodynamics

 ! return liquid fluxes (needed for coupling)
 returnCanopyTranspiration = scalarCanopyTranspiration    ! canopy transpiration (kg m-2 s-1)
 returnCanopyEvaporation   = scalarCanopyEvaporation      ! canopy evaporation/condensation (kg m-2 s-1)
 returnGroundEvaporation   = scalarGroundEvaporation      ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)

 ! end associations
 end associate


 end subroutine vegNrgFlux


 ! *******************************************************************************************************
 ! public subroutine wettedFrac: compute wetted fraction of the canopy
 ! *******************************************************************************************************
 subroutine wettedFrac(&
                       ! input
                       deriv,                  & ! flag to denote if derivative is desired
                       derNum,                 & ! flag to denote that numerical derivatives are required (otherwise, analytical derivatives are calculated)
                       frozen,                 & ! flag to denote if the canopy is frozen
                       dLiq_dT,                & ! derivative in canopy liquid w.r.t. canopy temperature (kg m-2 K-1)
                       fracLiq,                & ! fraction of liquid water on the canopy (-)
                       canopyLiq,              & ! canopy liquid water (kg m-2)
                       canopyIce,              & ! canopy ice (kg m-2)
                       canopyLiqMax,           & ! maximum canopy liquid water (kg m-2)
                       canopyIceMax,           & ! maximum canopy ice content (kg m-2)
                       canopyWettingFactor,    & ! maximum wetted fraction of the canopy (-)
                       canopyWettingExp,       & ! exponent in canopy wetting function (-)
                       ! output
                       canopyWetFraction,      & ! canopy wetted fraction (-)
                       dCanopyWetFraction_dWat,& ! derivative in wetted fraction w.r.t. canopy total water (kg-1 m2)
                       dCanopyWetFraction_dT,  & ! derivative in wetted fraction w.r.t. canopy temperature (K-1)
                       err,message)              ! error control
 implicit none
 ! input
 logical(lgt),intent(in)       :: deriv                   ! flag to denote if derivative is desired
 logical(lgt),intent(in)       :: derNum                  ! flag to denote that numerical derivatives are required (otherwise, analytical derivatives are calculated)
 logical(lgt),intent(in)       :: frozen                  ! flag to denote if the canopy is frozen
 real(rkind),intent(in)           :: dLiq_dT                 ! derivative in canopy liquid w.r.t. canopy temperature (kg m-2 K-1)
 real(rkind),intent(in)           :: fracLiq                 ! fraction of liquid water on the canopy (-)
 real(rkind),intent(in)           :: canopyLiq               ! canopy liquid water (kg m-2)
 real(rkind),intent(in)           :: canopyIce               ! canopy ice (kg m-2)
 real(rkind),intent(in)           :: canopyLiqMax            ! maximum canopy liquid water (kg m-2)
 real(rkind),intent(in)           :: canopyIceMax            ! maximum canopy ice content (kg m-2)
 real(rkind),intent(in)           :: canopyWettingFactor     ! maximum wetted fraction of the canopy (-)
 real(rkind),intent(in)           :: canopyWettingExp        ! exponent in canopy wetting function (-)
 ! output
 real(rkind),intent(out)          :: canopyWetFraction       ! canopy wetted fraction (-)
 real(rkind),intent(out)          :: dCanopyWetFraction_dWat ! derivative in wetted fraction w.r.t. canopy total water (kg-1 m2)
 real(rkind),intent(out)          :: dCanopyWetFraction_dT   ! derivative in wetted fraction w.r.t. canopy temperature (K-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                     ! error code
 character(*),intent(out)      :: message                 ! error message
 ! local variables
 logical(lgt),parameter        :: smoothing=.true.        ! flag to denote that smoothing is required
 real(rkind)                      :: canopyWetFractionPert   ! canopy wetted fraction after state perturbations (-)
 real(rkind)                      :: canopyWetFractionDeriv  ! derivative in wetted fraction w.r.t. canopy liquid water (kg-1 m2)
 ! -----------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='wettedFrac/'

 ! compute case where the canopy is frozen
 if(frozen)then
  ! compute fraction of liquid water on the canopy
  call wetFraction((deriv .and. .not.derNum),smoothing,canopyIce,canopyIceMax,canopyWettingFactor,canopyWettingExp,canopyWetFraction,canopyWetFractionDeriv)
  ! compute numerical derivative, if derivative is desired
  if(deriv.and.derNum)then
   call wetFraction((deriv .and. .not.derNum),smoothing,canopyIce+dx,canopyIceMax,canopyWettingFactor,canopyWettingExp,canopyWetFractionPert,canopyWetFractionDeriv)
   canopyWetFractionDeriv = (canopyWetFractionPert - canopyWetFraction)/dx
  end if
  ! scale derivative by the fraction of water
  ! NOTE: dIce/dWat = (1._rkind - fracLiq), hence dWet/dWat = dIce/dWat . dWet/dLiq
  dCanopyWetFraction_dWat = canopyWetFractionDeriv*(1._rkind - fracLiq)
  dCanopyWetFraction_dT   = -canopyWetFractionDeriv*dLiq_dT  ! NOTE: dIce/dT = -dLiq/dT
  return
 end if

 ! compute fraction of liquid water on the canopy
 ! NOTE: if(.not.deriv) canopyWetFractionDeriv = 0._rkind
 call wetFraction((deriv .and. .not.derNum),smoothing,canopyLiq,canopyLiqMax,canopyWettingFactor,canopyWettingExp,canopyWetFraction,canopyWetFractionDeriv)

 ! compute numerical derivative
 if(deriv.and.derNum)then
  call wetFraction((deriv .and. .not.derNum),smoothing,canopyLiq+dx,canopyLiqMax,canopyWettingFactor,canopyWettingExp,canopyWetFractionPert,canopyWetFractionDeriv)
  canopyWetFractionDeriv = (canopyWetFractionPert - canopyWetFraction)/dx
 end if

 ! scale derivative by the fraction of water
 ! NOTE: dLiq/dWat = fracLiq, hence dWet/dWat = dLiq/dWat . dWet/dLiq
 dCanopyWetFraction_dWat = canopyWetFractionDeriv*fracLiq
 dCanopyWetFraction_dT   = canopyWetFractionDeriv*dLiq_dT

 ! test
 !write(*,'(a,1x,2(L1,1x),10(f20.10,1x))') 'deriv, derNum, canopyWetFraction, canopyWetFractionDeriv = ', deriv, derNum, canopyWetFraction, canopyWetFractionDeriv
 !if(deriv) pause 'testing canopy wet fraction'

 end subroutine wettedFrac


 ! *******************************************************************************************************
 ! private subroutine wetFraction: compute fraction of canopy covered with liquid water
 ! *******************************************************************************************************
 subroutine wetFraction(derDesire,smoothing,canopyLiq,canopyMax,canopyWettingFactor,canopyWettingExp,canopyWetFraction,canopyWetFractionDeriv)
 implicit none
 ! dummy variables
 logical(lgt),intent(in) :: derDesire              ! flag to denote if analytical derivatives are desired
 logical(lgt),intent(in) :: smoothing              ! flag to denote if smoothing is required
 real(rkind),intent(in)     :: canopyLiq              ! liquid water content (kg m-2)
 real(rkind),intent(in)     :: canopyMax              ! liquid water content (kg m-2)
 real(rkind),intent(in)     :: canopyWettingFactor    ! maximum wetted fraction of the canopy (-)
 real(rkind),intent(in)     :: canopyWettingExp       ! exponent in canopy wetting function (-)

 real(rkind),intent(out)    :: canopyWetFraction      ! canopy wetted fraction (-)
 real(rkind),intent(out)    :: canopyWetFractionDeriv ! derivative in wetted fraction w.r.t. canopy liquid water (kg-1 m2)
 ! local variables
 real(rkind)                :: relativeCanopyWater    ! water stored on vegetation canopy, expressed as a fraction of maximum storage (-)
 real(rkind)                :: rawCanopyWetFraction   ! initial value of the canopy wet fraction (before smoothing)
 real(rkind)                :: rawWetFractionDeriv    ! derivative in canopy wet fraction w.r.t. storage (kg-1 m2)
 real(rkind)                :: smoothFunc             ! smoothing function used to improve numerical stability at times with limited water storage (-)
 real(rkind)                :: smoothFuncDeriv        ! derivative in the smoothing function w.r.t.canopy storage (kg-1 m2)
 real(rkind)                :: verySmall=epsilon(1._rkind) ! a very small number
 ! --------------------------------------------------------------------------------------------------------------

 ! compute relative canopy water
 relativeCanopyWater = canopyLiq/canopyMax
 !write(*,'(a,1x,e20.10,1x,2(f20.10,1x))') 'relativeCanopyWater, canopyLiq, canopyMax = ', relativeCanopyWater, canopyLiq, canopyMax

 ! compute an initial value of the canopy wet fraction
 ! - canopy below value where canopy is 100% wet
 if(relativeCanopyWater < 1._rkind)then
  rawCanopyWetFraction = canopyWettingFactor*(relativeCanopyWater**canopyWettingExp)
  if(derDesire .and. relativeCanopyWater>verySmall)then
   rawWetFractionDeriv = (canopyWettingFactor*canopyWettingExp/canopyMax)*relativeCanopyWater**(canopyWettingExp - 1._rkind)
  else
   rawWetFractionDeriv = 0._rkind
  end if

 ! - canopy is at capacity (canopyWettingFactor)
 else
  rawCanopyWetFraction = canopyWettingFactor
  rawWetFractionDeriv  = 0._rkind
 end if

 ! smooth canopy wetted fraction
 if(smoothing)then
  call logisticSmoother(derDesire,canopyLiq,smoothFunc,smoothFuncDeriv)
  canopyWetFraction = rawCanopyWetFraction*smoothFunc  ! logistic smoother
 else
  canopyWetFraction      = rawCanopyWetFraction
  canopyWetFractionDeriv = rawWetFractionDeriv
 end if

 ! compute derivative (product rule)
 if(derDesire .and. smoothing)then  ! NOTE: raw derivative is used if not smoothing
  canopyWetFractionDeriv = rawWetFractionDeriv*smoothFunc + rawCanopyWetFraction*smoothFuncDeriv
 else
  canopyWetFractionDeriv = 0._rkind
 end if

 end subroutine wetFraction


 ! *******************************************************************************************************
 ! private subroutine logisticSmoother: compute the smoothing function
 ! *******************************************************************************************************
 subroutine logisticSmoother(derDesire,canopyLiq,smoothFunc,smoothFuncDeriv)
 implicit none
 ! dummy variables
 logical(lgt),intent(in) :: derDesire              ! flag to denote if analytical derivatives are desired
 real(rkind),intent(in)     :: canopyLiq              ! liquid water content (kg m-2)
 real(rkind),intent(out)    :: smoothFunc             ! smoothing function (-)
 real(rkind),intent(out)    :: smoothFuncDeriv        ! derivative in smoothing function (kg-1 m-2)
 ! local variables
 real(rkind)                :: xArg                   ! argument used in the smoothing function (-)
 real(rkind)                :: expX                   ! exp(-xArg) -- used multiple times
 real(rkind),parameter      :: smoothThresh=0.01_rkind   ! mid-point of the smoothing function (kg m-2)
 real(rkind),parameter      :: smoothScale=0.001_rkind   ! scaling factor for the smoothing function (kg m-2)
 real(rkind),parameter      :: xLimit=50._rkind          ! don't compute exponents for > xLimit
 ! --------------------------------------------------------------------------------------------------------------
 ! compute argument in the smoothing function
 xArg = (canopyLiq - smoothThresh)/smoothScale

 ! only compute smoothing function for small exponents
 if(xArg > -xLimit .and. xArg < xLimit)then  ! avoid huge exponents
  expX            = exp(-xarg)                                   ! (also used in the derivative)
  smoothFunc      = 1._rkind / (1._rkind + expX)                       ! (logistic smoother)
  if(derDesire)then
   smoothFuncDeriv = expX / (smoothScale * (1._rkind + expX)**2._rkind) ! (derivative in the smoothing function)
  else
   smoothFuncDeriv = 0._rkind
  end if

 ! outside limits: special case of smooth exponents
 else
  if(xArg < 0._rkind)then; smoothFunc = 0._rkind   ! xArg < -xLimit
  else;                 smoothFunc = 1._rkind   ! xArg >  xLimit
  end if
  smoothFuncDeriv = 0._rkind
 end if  ! check for huge exponents

 end subroutine logisticSmoother
 ! --------------------------------------------------------------------------------------------------------------


 ! *******************************************************************************************************
 ! private subroutine longwaveBal: compute longwave radiation balance at the canopy and ground surface
 ! *******************************************************************************************************
 subroutine longwaveBal(&
                        ! input: model control
                        ixDerivMethod,                  & ! intent(in): choice of method used to compute derivative (analytical or numerical)
                        computeVegFlux,                 & ! intent(in): flag to compute fluxes over vegetation
                        ! input: canopy and ground temperature
                        canopyTemp,                     & ! intent(in): canopy temperature (K)
                        groundTemp,                     & ! intent(in): ground temperature (K)
                        ! input: canopy and ground emissivity
                        emc,                            & ! intent(in): canopy emissivity (-)
                        emg,                            & ! intent(in): ground emissivity (-)
                        ! input: forcing
                        LWRadUbound,                    & ! intent(in): downwelling longwave radiation at the upper boundary (W m-2)
                        ! output: sources
                        LWRadCanopy,                    & ! intent(out): longwave radiation emitted from the canopy (W m-2)
                        LWRadGround,                    & ! intent(out): longwave radiation emitted at the ground surface (W m-2)
                        ! output: individual fluxes
                        LWRadUbound2Canopy,             & ! intent(out): downward atmospheric longwave radiation absorbed by the canopy (W m-2)
                        LWRadUbound2Ground,             & ! intent(out): downward atmospheric longwave radiation absorbed by the ground (W m-2)
                        LWRadUbound2Ubound,             & ! intent(out): atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)
                        LWRadCanopy2Ubound,             & ! intent(out): longwave radiation emitted from canopy lost thru upper boundary (W m-2)
                        LWRadCanopy2Ground,             & ! intent(out): longwave radiation emitted from canopy absorbed by the ground (W m-2)
                        LWRadCanopy2Canopy,             & ! intent(out): canopy longwave reflected from ground and absorbed by the canopy (W m-2)
                        LWRadGround2Ubound,             & ! intent(out): longwave radiation emitted from ground lost thru upper boundary (W m-2)
                        LWRadGround2Canopy,             & ! intent(out): longwave radiation emitted from ground and absorbed by the canopy (W m-2)
                        ! output: net fluxes
                        LWNetCanopy,                    & ! intent(out): net longwave radiation at the canopy (W m-2)
                        LWNetGround,                    & ! intent(out): net longwave radiation at the ground surface (W m-2)
                        LWNetUbound,                    & ! intent(out): net longwave radiation at the upper boundary (W m-2)
                        ! output: flux derivatives
                        dLWNetCanopy_dTCanopy,          & ! intent(out): derivative in net canopy radiation w.r.t. canopy temperature (W m-2 K-1)
                        dLWNetGround_dTGround,          & ! intent(out): derivative in net ground radiation w.r.t. ground temperature (W m-2 K-1)
                        dLWNetCanopy_dTGround,          & ! intent(out): derivative in net canopy radiation w.r.t. ground temperature (W m-2 K-1)
                        dLWNetGround_dTCanopy,          & ! intent(out): derivative in net ground radiation w.r.t. canopy temperature (W m-2 K-1)
                        ! output: error control
                        err,message                     ) ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input: model control
 integer(i4b),intent(in)       :: ixDerivMethod            ! choice of method used to compute derivative (analytical or numerical)
 logical(lgt),intent(in)       :: computeVegFlux           ! flag to indicate if computing fluxes over vegetation
 ! input: canopy and ground temperature
 real(rkind),intent(in)           :: canopyTemp               ! canopy temperature (K)
 real(rkind),intent(in)           :: groundTemp               ! ground temperature (K)
 ! input: canopy and ground emissivity
 real(rkind),intent(in)           :: emc                      ! canopy emissivity (-)
 real(rkind),intent(in)           :: emg                      ! ground emissivity (-)
 ! input: forcing
 real(rkind),intent(in)           :: LWRadUbound              ! downwelling longwave radiation at the upper boundary (W m-2)
 ! output: sources
 real(rkind),intent(out)          :: LWRadCanopy              ! longwave radiation emitted from the canopy (W m-2)
 real(rkind),intent(out)          :: LWRadGround              ! longwave radiation emitted at the ground surface (W m-2)
 ! output: individual fluxes
 real(rkind),intent(out)          :: LWRadUbound2Canopy       ! downward atmospheric longwave radiation absorbed by the canopy (W m-2)
 real(rkind),intent(out)          :: LWRadUbound2Ground       ! downward atmospheric longwave radiation absorbed by the ground (W m-2)
 real(rkind),intent(out)          :: LWRadUbound2Ubound       ! atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)
 real(rkind),intent(out)          :: LWRadCanopy2Ubound       ! longwave radiation emitted from canopy lost thru upper boundary (W m-2)
 real(rkind),intent(out)          :: LWRadCanopy2Ground       ! longwave radiation emitted from canopy absorbed by the ground (W m-2)
 real(rkind),intent(out)          :: LWRadCanopy2Canopy       ! canopy longwave reflected from ground and absorbed by the canopy (W m-2)
 real(rkind),intent(out)          :: LWRadGround2Ubound       ! longwave radiation emitted from ground lost thru upper boundary (W m-2)
 real(rkind),intent(out)          :: LWRadGround2Canopy       ! longwave radiation emitted from ground and absorbed by the canopy (W m-2)
 ! output: net fluxes
 real(rkind),intent(out)          :: LWNetCanopy              ! net longwave radiation at the canopy (W m-2)
 real(rkind),intent(out)          :: LWNetGround              ! net longwave radiation at the ground surface (W m-2)
 real(rkind),intent(out)          :: LWNetUbound              ! net longwave radiation at the upper boundary (W m-2)
 ! output: flux derivatives
 real(rkind),intent(out)          :: dLWNetCanopy_dTCanopy    ! derivative in net canopy radiation w.r.t. canopy temperature (W m-2 K-1)
 real(rkind),intent(out)          :: dLWNetGround_dTGround    ! derivative in net ground radiation w.r.t. ground temperature (W m-2 K-1)
 real(rkind),intent(out)          :: dLWNetCanopy_dTGround    ! derivative in net canopy radiation w.r.t. ground temperature (W m-2 K-1)
 real(rkind),intent(out)          :: dLWNetGround_dTCanopy    ! derivative in net ground radiation w.r.t. canopy temperature (W m-2 K-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b),parameter        :: unperturbed=1            ! named variable to identify the case of unperturbed state variables
 integer(i4b),parameter        :: perturbStateCanopy=2     ! named variable to identify the case where we perturb the canopy temperature
 integer(i4b),parameter        :: perturbStateGround=3     ! named variable to identify the case where we perturb the ground temperature
 integer(i4b)                  :: itry                     ! index of flux evaluation
 integer(i4b)                  :: nFlux                    ! number of flux evaluations
 real(rkind)                      :: TCan                     ! value of canopy temperature used in flux calculations (may be perturbed)
 real(rkind)                      :: TGnd                     ! value of ground temperature used in flux calculations (may be perturbed)
 real(rkind)                      :: fluxBalance              ! check energy closure (W m-2)
 real(rkind),parameter            :: fluxTolerance=1.e-10_rkind  ! tolerance for energy closure (W m-2)
 real(rkind)                      :: dLWRadCanopy_dTCanopy    ! derivative in emitted radiation at the canopy w.r.t. canopy temperature
 real(rkind)                      :: dLWRadGround_dTGround    ! derivative in emitted radiation at the ground w.r.t. ground temperature
 real(rkind)                      :: LWNetCanopy_dStateCanopy ! net lw canopy flux after perturbation in canopy temperature
 real(rkind)                      :: LWNetGround_dStateCanopy ! net lw ground flux after perturbation in canopy temperature
 real(rkind)                      :: LWNetCanopy_dStateGround ! net lw canopy flux after perturbation in ground temperature
 real(rkind)                      :: LWNetGround_dStateGround ! net lw ground flux after perturbation in ground temperature
 ! -----------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='longwaveBal/'

 ! check the need to compute numerical derivatives
 if(ixDerivMethod==numerical)then
  nFlux=3  ! compute the derivatives using one-sided finite differences
 else
  nFlux=1  ! compute analytical derivatives
 end if

 ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
 do itry=nFlux,1,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)

  !print*, 'perturbation: ', (itry==unperturbed), (itry==perturbStateCanopy), (itry==perturbStateGround)

  ! -------------------------------------------------------------------------------------
  ! state perturbations for numerical deriavtives with one-sided finite differences
  ! note: no perturbations performed using analytical derivatives (nFlux=1)
  ! -------------------------------------------------------------------------------------

  ! identify the type of perturbation
  select case(itry)

   ! un-perturbed case
   case(unperturbed)
    TCan = canopyTemp
    TGnd = groundTemp

   ! perturb canopy temperature
   case(perturbStateCanopy)
    TCan = canopyTemp + dx
    TGnd = groundTemp

   ! perturb ground temperature
   case(perturbStateGround)
    TCan = canopyTemp
    TGnd = groundTemp + dx

   ! check for an unknown perturbation
   case default; err=10; message=trim(message)//"unknown perturbation"; return

  end select ! (type of perturbation)

  ! -------------------------------------------------------------------------------------
  ! calculation block (unperturbed fluxes returned [computed last])
  ! -------------------------------------------------------------------------------------
  ! NOTE: emc should be set to zero when not computing canopy fluxes

  ! compute longwave fluxes from canopy and the ground
  if(computeVegFlux)then
   LWRadCanopy = emc*sb*TCan**4._rkind                                           ! longwave radiation emitted from the canopy (W m-2)
  else
   LWRadCanopy = 0._rkind
  end if
  LWRadGround = emg*sb*TGnd**4._rkind                                           ! longwave radiation emitted at the ground surface (W m-2)

  ! compute fluxes originating from the atmosphere
  LWRadUbound2Canopy = (emc + (1._rkind - emc)*(1._rkind - emg)*emc)*LWRadUbound   ! downward atmospheric longwave radiation absorbed by the canopy (W m-2)
  LWRadUbound2Ground = (1._rkind - emc)*emg*LWRadUbound                         ! downward atmospheric longwave radiation absorbed by the ground (W m-2)
  LWRadUbound2Ubound = (1._rkind - emc)*(1._rkind - emg)*(1._rkind - emc)*LWRadUbound ! atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)

  ! compute fluxes originating from the canopy
  LWRadCanopy2Ubound = (1._rkind + (1._rkind - emc)*(1._rkind - emg))*LWRadCanopy     ! longwave radiation emitted from canopy lost thru upper boundary (W m-2)
  LWRadCanopy2Ground = emg*LWRadCanopy                                       ! longwave radiation emitted from canopy absorbed by the ground (W m-2)
  LWRadCanopy2Canopy = emc*(1._rkind - emg)*LWRadCanopy                         ! canopy longwave reflected from ground and absorbed by the canopy (W m-2)

  ! compute fluxes originating from the ground surface
  LWRadGround2Ubound = (1._rkind - emc)*LWRadGround                             ! longwave radiation emitted from ground lost thru upper boundary (W m-2)
  LWRadGround2Canopy = emc*LWRadGround                                       ! longwave radiation emitted from ground and absorbed by the canopy (W m-2)

  ! compute net longwave radiation (W m-2)
  LWNetCanopy = LWRadUbound2Canopy + LWRadGround2Canopy + LWRadCanopy2Canopy - 2._rkind*LWRadCanopy  ! canopy
  LWNetGround = LWRadUbound2Ground + LWRadCanopy2Ground - LWRadGround                             ! ground surface
  LWNetUbound = LWRadUbound - LWRadUbound2Ubound - LWRadCanopy2Ubound - LWRadGround2Ubound                             ! upper boundary

  !print*, 'LWRadCanopy = ', LWRadCanopy
  !print*, 'LWRadGround = ', LWRadGround

  !print*, 'LWNetCanopy = ', LWNetCanopy
  !print*, 'LWNetGround = ', LWNetGround
  !print*, 'LWNetUbound = ', LWNetUbound

  ! check the flux balance
  fluxBalance = LWNetUbound - (LWNetCanopy + LWNetGround)
  if(abs(fluxBalance) > fluxTolerance)then
   print*, 'fluxBalance = ', fluxBalance
   print*, 'emg, emc = ', emg, emc
   print*, 'TCan, TGnd = ', TCan, TGnd
   print*, 'LWRadUbound = ', LWRadUbound
   print*, 'LWRadCanopy = ', LWRadCanopy
   print*, 'LWRadGround = ', LWRadGround
   print*, 'LWRadUbound2Canopy = ', LWRadUbound2Canopy
   print*, 'LWRadUbound2Ground = ', LWRadUbound2Ground
   print*, 'LWRadUbound2Ubound = ', LWRadUbound2Ubound
   print*, 'LWRadCanopy2Ubound = ', LWRadCanopy2Ubound
   print*, 'LWRadCanopy2Ground = ', LWRadCanopy2Ground
   print*, 'LWRadCanopy2Canopy = ', LWRadCanopy2Canopy
   print*, 'LWRadGround2Ubound = ', LWRadGround2Ubound
   print*, 'LWRadGround2Canopy = ', LWRadGround2Canopy
   print*, 'LWNetCanopy = ', LWNetCanopy
   print*, 'LWNetGround = ', LWNetGround
   print*, 'LWNetUbound = ', LWNetUbound
   message=trim(message)//'flux imbalance'
   err=20; return
  end if

  ! --------------------------------------------------------------------------------------
  ! save perturbed fluxes to calculate numerical derivatives (one-sided finite difference)
  ! --------------------------------------------------------------------------------------
  if(ixDerivMethod==numerical)then
   select case(itry) ! (select type of perturbation)
    case(unperturbed); exit
    case(perturbStateCanopy)
     LWNetCanopy_dStateCanopy = LWNetCanopy
     LWNetGround_dStateCanopy = LWNetGround
    case(perturbStateGround)
     LWNetCanopy_dStateGround = LWNetCanopy
     LWNetGround_dStateGround = LWNetGround
    case default; err=10; message=trim(message)//"unknown perturbation"; return
   end select ! (type of perturbation)
  end if ! (if numerical)

 end do  ! looping through different perturbations

 ! -------------------------------------------------------------------------------------
 ! compute derivatives
 ! -------------------------------------------------------------------------------------
 select case(ixDerivMethod)

  ! ***** analytical derivatives
  case(analytical)
   ! compute initial derivatives
   dLWRadCanopy_dTCanopy = 4._rkind*emc*sb*TCan**3._rkind
   dLWRadGround_dTGround = 4._rkind*emg*sb*TGnd**3._rkind
   ! compute analytical derivatives
   dLWNetCanopy_dTCanopy = (emc*(1._rkind - emg) - 2._rkind)*dLWRadCanopy_dTCanopy ! derivative in net canopy radiation w.r.t. canopy temperature (W m-2 K-1)
   dLWNetGround_dTGround = -dLWRadGround_dTGround     ! derivative in net ground radiation w.r.t. ground temperature (W m-2 K-1)
   dLWNetCanopy_dTGround = emc*dLWRadGround_dTGround  ! derivative in net canopy radiation w.r.t. ground temperature (W m-2 K-1)
   dLWNetGround_dTCanopy = emg*dLWRadCanopy_dTCanopy  ! derivative in net ground radiation w.r.t. canopy temperature (W m-2 K-1)

  ! ***** numerical derivatives
  case(numerical)
   ! compute numerical derivatives (one-sided finite differences)
   dLWNetCanopy_dTCanopy = (LWNetCanopy_dStateCanopy - LWNetCanopy)/dx  ! derivative in net canopy radiation w.r.t. canopy temperature (W m-2 K-1)
   dLWNetGround_dTGround = (LWNetGround_dStateGround - LWNetGround)/dx  ! derivative in net ground radiation w.r.t. ground temperature (W m-2 K-1)
   dLWNetCanopy_dTGround = (LWNetCanopy_dStateGround - LWNetCanopy)/dx  ! derivative in net canopy radiation w.r.t. ground temperature (W m-2 K-1)
   dLWNetGround_dTCanopy = (LWNetGround_dStateCanopy - LWNetGround)/dx  ! derivative in net ground radiation w.r.t. canopy temperature (W m-2 K-1)

  ! ***** error check
  case default; err=10; message=trim(message)//"unknown method to calculate derivatives"; return

 end select ! (type of method to calculate derivatives)

 end subroutine longwaveBal


 ! *******************************************************************************************************
 ! private subroutine aeroResist: compute aerodynamic resistances
 ! *******************************************************************************************************
 subroutine aeroResist(&
                       ! input: model control
                       computeVegFlux,                & ! intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                       derivDesired,                  & ! intent(in): flag to indicate if analytical derivatives are desired
                       ixVegTraits,                   & ! intent(in): choice of parameterization for vegetation roughness length and displacement height
                       ixWindProfile,                 & ! intent(in): choice of canopy wind profile
                       ixStability,                   & ! intent(in): choice of stability function
                       ! input: above-canopy forcing data
                       mHeight,                       & ! intent(in): measurement height (m)
                       airtemp,                       & ! intent(in): air temperature at some height above the surface (K)
                       windspd,                       & ! intent(in): wind speed at some height above the surface (m s-1)
                       ! input: temperature (canopy, ground, canopy air space)
                       canairTemp,                    & ! intent(in): temperature of the canopy air space (K)
                       groundTemp,                    & ! intent(in): ground temperature (K)
                       ! input: diagnostic variables
                       exposedVAI,                    & ! intent(in): exposed vegetation area index -- leaf plus stem (m2 m-2)
                       snowDepth,                     & ! intent(in): snow depth (m)
                       ! input: parameters
                       z0Ground,                      & ! intent(in): roughness length of the ground (below canopy or non-vegetated surface [snow]) (m)
                       z0CanopyParam,                 & ! intent(in): roughness length of the canopy (m)
                       zpdFraction,                   & ! intent(in): zero plane displacement / canopy height (-)
                       critRichNumber,                & ! intent(in): critical value for the bulk Richardson number where turbulence ceases (-)
                       Louis79_bparam,                & ! intent(in): parameter in Louis (1979) stability function
                       Mahrt87_eScale,                & ! intent(in): exponential scaling factor in the Mahrt (1987) stability function
                       windReductionParam,            & ! intent(in): canopy wind reduction parameter (-)
                       leafExchangeCoeff,             & ! intent(in): turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
                       leafDimension,                 & ! intent(in): characteristic leaf dimension (m)
                       heightCanopyTop,               & ! intent(in): height at the top of the vegetation canopy (m)
                       heightCanopyBottom,            & ! intent(in): height at the bottom of the vegetation canopy (m)
                       ! output: stability corrections
                       RiBulkCanopy,                  & ! intent(out): bulk Richardson number for the canopy (-)
                       RiBulkGround,                  & ! intent(out): bulk Richardson number for the ground surface (-)
                       canopyStabilityCorrection,     & ! intent(out): stability correction for the canopy (-)
                       groundStabilityCorrection,     & ! intent(out): stability correction for the ground surface (-)
                       ! output: scalar resistances
                       z0Canopy,                      & ! intent(out): roughness length of the canopy (m)
                       windReductionFactor,           & ! intent(out): canopy wind reduction factor (-)
                       zeroPlaneDisplacement,         & ! intent(out): zero plane displacement (m)
                       eddyDiffusCanopyTop,           & ! intent(out): eddy diffusivity for heat at the top of the canopy (m2 s-1)
                       frictionVelocity,              & ! intent(out): friction velocity (m s-1)
                       windspdCanopyTop,              & ! intent(out): windspeed at the top of the canopy (m s-1)
                       windspdCanopyBottom,           & ! intent(out): windspeed at the height of the bottom of the canopy (m s-1)
                       leafResistance,                & ! intent(out): mean leaf boundary layer resistance per unit leaf area (s m-1)
                       groundResistance,              & ! intent(out): below canopy aerodynamic resistance (s m-1)
                       canopyResistance,              & ! intent(out): above canopy aerodynamic resistance (s m-1)
                       ! output: derivatives in scalar resistances
                       dGroundResistance_dTGround,    & ! intent(out): derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
                       dGroundResistance_dTCanopy,    & ! intent(out): derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
                       dGroundResistance_dTCanair,    & ! intent(out): derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
                       dCanopyResistance_dTCanopy,    & ! intent(out): derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
                       dCanopyResistance_dTCanair,    & ! intent(out): derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
                       ! output: error control
                       err,message                    ) ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! compute aerodynamic resistances
 ! Refs: Choudhury and Monteith (4-layer model for heat budget of homogenous surfaces; QJRMS, 1988)
 !       Niu and Yang (Canopy effects on snow processes; JGR, 2004)
 !       Mahat et al. (Below-canopy turbulence in a snowmelt model, WRR, 2012)
 implicit none
 ! input: model control
 logical(lgt),intent(in)       :: computeVegFlux                ! logical flag to compute vegetation fluxes (.false. if veg buried by snow)
 logical(lgt),intent(in)       :: derivDesired                  ! logical flag to indicate if analytical derivatives are desired
 integer(i4b),intent(in)       :: ixVegTraits                   ! choice of parameterization for vegetation roughness length and displacement height
 integer(i4b),intent(in)       :: ixWindProfile                 ! choice of canopy wind profile
 integer(i4b),intent(in)       :: ixStability                   ! choice of stability function
 ! input: above-canopy forcing data
 real(rkind),intent(in)           :: mHeight                       ! measurement height (m)
 real(rkind),intent(in)           :: airtemp                       ! air temperature at some height above the surface (K)
 real(rkind),intent(in)           :: windspd                       ! wind speed at some height above the surface (m s-1)
 ! input: temperature (canopy, ground, canopy air space)
 real(rkind),intent(in)           :: canairTemp                    ! temperature of the canopy air space (K)
 real(rkind),intent(in)           :: groundTemp                    ! ground temperature (K)
 ! input: diagnostic variables
 real(rkind),intent(in)           :: exposedVAI                    ! exposed vegetation area index -- leaf plus stem (m2 m-2)
 real(rkind),intent(in)           :: snowDepth                     ! snow depth (m)
 ! input: parameters
 real(rkind),intent(in)           :: z0Ground                      ! roughness length of the ground (below canopy or non-vegetated surface [snow]) (m)
 real(rkind),intent(in)           :: z0CanopyParam                 ! roughness length of the canopy (m)
 real(rkind),intent(in)           :: zpdFraction                   ! zero plane displacement / canopy height (-)
 real(rkind),intent(in)           :: critRichNumber                ! critical value for the bulk Richardson number where turbulence ceases (-)
 real(rkind),intent(in)           :: Louis79_bparam                ! parameter in Louis (1979) stability function
 real(rkind),intent(in)           :: Mahrt87_eScale                ! exponential scaling factor in the Mahrt (1987) stability function
 real(rkind),intent(in)           :: windReductionParam            ! canopy wind reduction parameter (-)
 real(rkind),intent(in)           :: leafExchangeCoeff             ! turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
 real(rkind),intent(in)           :: leafDimension                 ! characteristic leaf dimension (m)
 real(rkind),intent(in)           :: heightCanopyTop               ! height at the top of the vegetation canopy (m)
 real(rkind),intent(in)           :: heightCanopyBottom            ! height at the bottom of the vegetation canopy (m)
 ! output: stability corrections
 real(rkind),intent(out)          :: RiBulkCanopy                  ! bulk Richardson number for the canopy (-)
 real(rkind),intent(out)          :: RiBulkGround                  ! bulk Richardson number for the ground surface (-)
 real(rkind),intent(out)          :: canopyStabilityCorrection     ! stability correction for the canopy (-)
 real(rkind),intent(out)          :: groundStabilityCorrection     ! stability correction for the ground surface (-)
 ! output: scalar resistances
 real(rkind),intent(out)          :: z0Canopy                      ! roughness length of the vegetation canopy (m)
 real(rkind),intent(out)          :: windReductionFactor           ! canopy wind reduction factor (-)
 real(rkind),intent(out)          :: zeroPlaneDisplacement         ! zero plane displacement (m)
 real(rkind),intent(out)          :: eddyDiffusCanopyTop           ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
 real(rkind),intent(out)          :: frictionVelocity              ! friction velocity (m s-1)
 real(rkind),intent(out)          :: windspdCanopyTop              ! windspeed at the top of the canopy (m s-1)
 real(rkind),intent(out)          :: windspdCanopyBottom           ! windspeed at the height of the bottom of the canopy (m s-1)
 real(rkind),intent(out)          :: leafResistance                ! mean leaf boundary layer resistance per unit leaf area (s m-1)
 real(rkind),intent(out)          :: groundResistance              ! below canopy aerodynamic resistance (s m-1)
 real(rkind),intent(out)          :: canopyResistance              ! above canopy aerodynamic resistance (s m-1)
 ! output: derivatives in scalar resistances
 real(rkind),intent(out)          :: dGroundResistance_dTGround    ! derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
 real(rkind),intent(out)          :: dGroundResistance_dTCanopy    ! derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
 real(rkind),intent(out)          :: dGroundResistance_dTCanair    ! derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
 real(rkind),intent(out)          :: dCanopyResistance_dTCanopy    ! derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
 real(rkind),intent(out)          :: dCanopyResistance_dTCanair    ! derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                           ! error code
 character(*),intent(out)      :: message                       ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! local variables: general
 character(LEN=256)            :: cmessage                      ! error message of downwind routine
 ! local variables: vegetation roughness and dispalcement height
 real(rkind),parameter            :: oneThird=1._rkind/3._rkind          ! 1/3
 real(rkind),parameter            :: twoThirds=2._rkind/3._rkind         ! 2/3
 real(rkind),parameter            :: C_r = 0.3                     ! roughness element drag coefficient (-) from Raupach (BLM, 1994)
 real(rkind),parameter            :: C_s = 0.003_rkind                ! substrate surface drag coefficient (-) from Raupach (BLM, 1994)
 real(rkind),parameter            :: approxDragCoef_max = 0.3_rkind   ! maximum value of the approximate drag coefficient (-) from Raupach (BLM, 1994)
 real(rkind),parameter            :: psi_h = 0.193_rkind              ! roughness sub-layer influence function (-) from Raupach (BLM, 1994)
 real(rkind),parameter            :: c_d1 = 7.5_rkind                 ! scaling parameter used to define displacement height (-) from Raupach (BLM, 1994)
 real(rkind),parameter            :: cd_CM = 0.2_rkind                ! mean drag coefficient for individual leaves (-) from Choudhury and Monteith (QJRMS, 1988)
 real(rkind)                      :: funcLAI                       ! temporary variable to calculate zero plane displacement for the canopy
 real(rkind)                      :: fracCanopyHeight              ! zero plane displacement expressed as a fraction of canopy height
 real(rkind)                      :: approxDragCoef                ! approximate drag coefficient used in the computation of canopy roughness length (-)
 ! local variables: resistance
 real(rkind)                      :: canopyExNeut                  ! surface-atmosphere exchange coefficient under neutral conditions (-)
 real(rkind)                      :: groundExNeut                  ! surface-atmosphere exchange coefficient under neutral conditions (-)
 real(rkind)                      :: sfc2AtmExchangeCoeff_canopy   ! surface-atmosphere exchange coefficient after stability corrections (-)
 real(rkind)                      :: groundResistanceNeutral       ! ground resistance under neutral conditions (s m-1)
 real(rkind)                      :: windConvFactor_fv             ! factor to convert friction velocity to wind speed at top of canopy (-)
 real(rkind)                      :: windConvFactor                ! factor to convert wind speed at top of canopy to wind speed at a given height in the canopy (-)
 real(rkind)                      :: referenceHeight               ! z0Canopy+zeroPlaneDisplacement (m)
 real(rkind)                      :: windspdRefHeight              ! windspeed at the reference height (m/s)
 real(rkind)                      :: heightAboveGround             ! height above the snow surface (m)
 real(rkind)                      :: heightCanopyTopAboveSnow      ! height at the top of the vegetation canopy relative to snowpack (m)
 real(rkind)                      :: heightCanopyBottomAboveSnow   ! height at the bottom of the vegetation canopy relative to snowpack (m)
 real(rkind),parameter            :: xTolerance=0.1_rkind             ! tolerance to handle the transition from exponential to log-below canopy
 ! local variables: derivatives
 real(rkind)                      :: dFV_dT                        ! derivative in friction velocity w.r.t. canopy air temperature
 real(rkind)                      :: dED_dT                        ! derivative in eddy diffusivity at the top of the canopy w.r.t. canopy air temperature
 real(rkind)                      :: dGR_dT                        ! derivative in neutral ground resistance w.r.t. canopy air temperature
 real(rkind)                      :: tmp1,tmp2                     ! temporary variables used in calculation of ground resistance
 real(rkind)                      :: dCanopyStabilityCorrection_dRich     ! derivative in stability correction w.r.t. Richardson number for the canopy (-)
 real(rkind)                      :: dGroundStabilityCorrection_dRich     ! derivative in stability correction w.r.t. Richardson number for the ground surface (-)
 real(rkind)                      :: dCanopyStabilityCorrection_dAirTemp  ! (not used) derivative in stability correction w.r.t. air temperature (K-1)
 real(rkind)                      :: dGroundStabilityCorrection_dAirTemp  ! (not used) derivative in stability correction w.r.t. air temperature (K-1)
 real(rkind)                      :: dCanopyStabilityCorrection_dCasTemp  ! derivative in canopy stability correction w.r.t. canopy air space temperature (K-1)
 real(rkind)                      :: dGroundStabilityCorrection_dCasTemp  ! derivative in ground stability correction w.r.t. canopy air space temperature (K-1)
 real(rkind)                      :: dGroundStabilityCorrection_dSfcTemp  ! derivative in ground stability correction w.r.t. surface temperature (K-1)
 real(rkind)                      :: singleLeafConductance         ! leaf boundary layer conductance (m s-1)
 real(rkind)                      :: canopyLeafConductance         ! leaf boundary layer conductance -- scaled up to the canopy (m s-1)
 real(rkind)                      :: leaf2CanopyScaleFactor        ! factor to scale from the leaf to the canopy [m s-(1/2)]
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='aeroResist/'

 ! check that measurement height is above the top of the canopy
 if(mHeight < heightCanopyTop)then
  err=20; message=trim(message)//'measurement height is below the top of the canopy'; return
 end if

 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! * compute vegetation poperties (could be done at the same time as phenology.. does not have to be in the flux routine!)
 if(computeVegFlux) then ! (if vegetation is exposed)

  ! ***** identify zero plane displacement, roughness length, and surface temperature for the canopy (m)
  ! First, calculate new coordinate system above snow - use these to scale wind profiles and resistances
  ! NOTE: the new coordinate system makes zeroPlaneDisplacement and z0Canopy consistent
  heightCanopyTopAboveSnow = heightCanopyTop - snowDepth
  heightCanopyBottomAboveSnow = max(heightCanopyBottom - snowDepth, 0.0_rkind)
  select case(ixVegTraits)

   ! Raupach (BLM 1994) "Simplified expressions..."
   case(Raupach_BLM1994)
    ! (compute zero-plane displacement)
    funcLAI          = sqrt(c_d1*exposedVAI)
    fracCanopyHeight = -(1._rkind - exp(-funcLAI))/funcLAI + 1._rkind
    zeroPlaneDisplacement = fracCanopyHeight*(heightCanopyTopAboveSnow-heightCanopyBottomAboveSnow)+heightCanopyBottomAboveSnow
    ! (coupute roughness length of the veg canopy)
    approxDragCoef   = min( sqrt(C_s + C_r*exposedVAI/2._rkind), approxDragCoef_max)
    z0Canopy         = (1._rkind - fracCanopyHeight) * exp(-vkc*approxDragCoef - psi_h) * (heightCanopyTopAboveSnow-heightCanopyBottomAboveSnow)

   ! Choudhury and Monteith (QJRMS 1988) "A four layer model for the heat budget..."
   case(CM_QJRMS1988)
    funcLAI =  cd_CM*exposedVAI
    zeroPlaneDisplacement = 1.1_rkind*heightCanopyTopAboveSnow*log(1._rkind + funcLAI**0.25_rkind)
    if(funcLAI < 0.2_rkind)then
     z0Canopy = z0Ground + 0.3_rkind*heightCanopyTopAboveSnow*funcLAI**0.5_rkind
    else
     z0Canopy = 0.3_rkind*heightCanopyTopAboveSnow*(1._rkind - zeroPlaneDisplacement/heightCanopyTopAboveSnow)
    end if

   ! constant parameters dependent on the vegetation type
   case(vegTypeTable)
    zeroPlaneDisplacement = zpdFraction*heightCanopyTopAboveSnow  ! zero-plane displacement (m)
    z0Canopy = z0CanopyParam                                      ! roughness length of the veg canopy (m)

   ! check
   case default
    err=10; message=trim(message)//"unknown parameterization for vegetation roughness length and displacement height"; return

  end select  ! vegetation traits (z0, zpd)

  ! check zero plane displacement
  if(zeroPlaneDisplacement < heightCanopyBottomAboveSnow)then
   write(*,'(a,1x,10(f12.5,1x))') 'heightCanopyTop, snowDepth, heightCanopyTopAboveSnow, heightCanopyBottomAboveSnow, exposedVAI = ', &
                                   heightCanopyTop, snowDepth, heightCanopyTopAboveSnow, heightCanopyBottomAboveSnow, exposedVAI
   message=trim(message)//'zero plane displacement is below the canopy bottom'
   err=20; return
  endif

  ! check measurement height
  if(mHeight < zeroPlaneDisplacement)then; err=20; message=trim(message)//'measurement height is below the displacement height'; return; end if
  if(mHeight < z0Canopy)then; err=20; message=trim(message)//'measurement height is below the roughness length'; return; end if

  ! -----------------------------------------------------------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------------------------------------------------------
  ! * compute resistance for the case where the canopy is exposed
  ! compute the stability correction for resistance from canopy air space to air above the canopy (-)
  call aStability(&
                  ! input
                  derivDesired,                                     & ! input: logical flag to compute analytical derivatives
                  ixStability,                                      & ! input: choice of stability function
                  ! input: forcing data, diagnostic and state variables
                  mHeight,                                          & ! input: measurement height (m)
                  airTemp,                                          & ! input: air temperature above the canopy (K)
                  canairTemp,                                       & ! input: temperature of the canopy air space (K)
                  windspd,                                          & ! input: wind speed above the canopy (m s-1)
                  ! input: stability parameters
                  critRichNumber,                                   & ! input: critical value for the bulk Richardson number where turbulence ceases (-)
                  Louis79_bparam,                                   & ! input: parameter in Louis (1979) stability function
                  Mahrt87_eScale,                                   & ! input: exponential scaling factor in the Mahrt (1987) stability function
                  ! output
                  RiBulkCanopy,                                     & ! output: bulk Richardson number (-)
                  canopyStabilityCorrection,                        & ! output: stability correction for turbulent heat fluxes (-)
                  dCanopyStabilityCorrection_dRich,                 & ! output: derivative in stability correction w.r.t. Richardson number for the canopy (-)
                  dCanopyStabilityCorrection_dAirTemp,              & ! output: (not used) derivative in stability correction w.r.t. air temperature (K-1)
                  dCanopyStabilityCorrection_dCasTemp,              & ! output: derivative in stability correction w.r.t. canopy air space temperature (K-1)
                  err, cmessage                                     ) ! output: error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

  ! compute turbulent exchange coefficient (-)
  canopyExNeut = (vkc**2._rkind) / ( log((mHeight - zeroPlaneDisplacement)/z0Canopy))**2._rkind     ! coefficient under conditions of neutral stability
  sfc2AtmExchangeCoeff_canopy = canopyExNeut*canopyStabilityCorrection                        ! after stability corrections

  ! compute the friction velocity (m s-1)
  frictionVelocity = windspd * sqrt(sfc2AtmExchangeCoeff_canopy)

  ! compute the above-canopy resistance (s m-1)
  canopyResistance = 1._rkind/(sfc2AtmExchangeCoeff_canopy*windspd)
  if(canopyResistance < 0._rkind)then; err=20; message=trim(message)//'canopy resistance < 0'; return; end if

  ! compute windspeed at the top of the canopy above snow depth (m s-1)
  ! NOTE: stability corrections cancel out
  windConvFactor_fv = log((heightCanopyTopAboveSnow - zeroPlaneDisplacement)/z0Canopy) / log((mHeight - snowDepth - zeroPlaneDisplacement)/z0Canopy)
  windspdCanopyTop  = windspd*windConvFactor_fv

  ! compute the windspeed reduction
  ! Refs: Norman et al. (Ag. Forest Met., 1995) -- citing Goudriaan (1977 manuscript "crop micrometeorology: a simulation study", Wageningen).
  windReductionFactor = windReductionParam * exposedVAI**twoThirds * (heightCanopyTopAboveSnow - heightCanopyBottomAboveSnow)**oneThird / leafDimension**oneThird

  ! compute windspeed at the height z0Canopy+zeroPlaneDisplacement (m s-1)
  referenceHeight   = z0Canopy+zeroPlaneDisplacement
  windConvFactor    = exp(-windReductionFactor*(1._rkind - (referenceHeight/heightCanopyTopAboveSnow)))
  windspdRefHeight  = windspdCanopyTop*windConvFactor

  ! compute windspeed at the bottom of the canopy relative to the snow depth (m s-1)
  windConvFactor       = exp(-windReductionFactor*(1._rkind - (heightCanopyBottomAboveSnow/heightCanopyTopAboveSnow)))
  windspdCanopyBottom  = windspdCanopyTop*windConvFactor

  ! compute the leaf boundary layer resistance (s m-1)
  singleLeafConductance  = leafExchangeCoeff*sqrt(windspdCanopyTop/leafDimension)
  leaf2CanopyScaleFactor = (2._rkind/windReductionFactor) * (1._rkind - exp(-windReductionFactor/2._rkind)) ! factor to scale from the leaf to the canopy
  canopyLeafConductance  = singleLeafConductance*leaf2CanopyScaleFactor
  leafResistance         = 1._rkind/(canopyLeafConductance)
  if(leafResistance < 0._rkind)then; err=20; message=trim(message)//'leaf resistance < 0'; return; end if

  ! compute eddy diffusivity for heat at the top of the canopy (m2 s-1)
  !  Note: use of friction velocity here includes stability adjustments
  !  Note: max used to avoid dividing by zero
  eddyDiffusCanopyTop = max(vkc*FrictionVelocity*(heightCanopyTopAboveSnow - zeroPlaneDisplacement), mpe)

  ! compute the resistance between the surface and canopy air UNDER NEUTRAL CONDITIONS (s m-1)

  ! case 1: assume exponential profile extends from the snow depth plus surface roughness length to the displacement height plus vegetation roughness
  if(ixWindProfile==exponential .or. heightCanopyBottomAboveSnow<z0Ground+xTolerance)then

   ! compute the neutral ground resistance
   tmp1 = exp(-windReductionFactor* z0Ground/heightCanopyTopAboveSnow)
   tmp2 = exp(-windReductionFactor*(z0Canopy+zeroPlaneDisplacement)/heightCanopyTopAboveSnow)
   groundResistanceNeutral = ( heightCanopyTopAboveSnow*exp(windReductionFactor) / (windReductionFactor*eddyDiffusCanopyTop) ) * (tmp1 - tmp2)   ! s m-1

  ! case 2: logarithmic profile from snow depth plus roughness height to bottom of the canopy
  ! NOTE: heightCanopyBottomAboveSnow>z0Ground+xTolerance
  else

   ! compute the neutral ground resistance
   ! (first, component between heightCanopyBottomAboveSnow and z0Canopy+zeroPlaneDisplacement)
   tmp1  = exp(-windReductionFactor* heightCanopyBottomAboveSnow/heightCanopyTopAboveSnow)
   tmp2  = exp(-windReductionFactor*(z0Canopy+zeroPlaneDisplacement)/heightCanopyTopAboveSnow)
   groundResistanceNeutral = ( heightCanopyTopAboveSnow*exp(windReductionFactor) / (windReductionFactor*eddyDiffusCanopyTop) ) * (tmp1 - tmp2)
   ! (add log-below-canopy component)
   groundResistanceNeutral = groundResistanceNeutral + (1._rkind/(max(0.1_rkind,windspdCanopyBottom)*vkc**2._rkind))*(log(heightCanopyBottomAboveSnow/z0Ground))**2._rkind

  endif  ! switch between exponential profile and log-below-canopy

  ! compute the stability correction for resistance from the ground to the canopy air space (-)
  ! NOTE: here we are interested in the windspeed at height z0Canopy+zeroPlaneDisplacement
  call aStability(&
                  ! input
                  derivDesired,                                     & ! input: logical flag to compute analytical derivatives
                  ixStability,                                      & ! input: choice of stability function
                  ! input: forcing data, diagnostic and state variables
                  referenceHeight,                                  & ! input: height of the canopy air space temperature/wind (m)
                  canairTemp,                                       & ! input: temperature of the canopy air space (K)
                  groundTemp,                                       & ! input: temperature of the ground surface (K)
                  max(0.1_rkind,windspdRefHeight),                     & ! input: wind speed at height z0Canopy+zeroPlaneDisplacement (m s-1)
                  ! input: stability parameters
                  critRichNumber,                                   & ! input: critical value for the bulk Richardson number where turbulence ceases (-)
                  Louis79_bparam,                                   & ! input: parameter in Louis (1979) stability function
                  Mahrt87_eScale,                                   & ! input: exponential scaling factor in the Mahrt (1987) stability function
                  ! output
                  RiBulkGround,                                     & ! output: bulk Richardson number (-)
                  groundStabilityCorrection,                        & ! output: stability correction for turbulent heat fluxes (-)
                  dGroundStabilityCorrection_dRich,                 & ! output: derivative in stability correction w.r.t. Richardson number for the canopy (-)
                  dGroundStabilityCorrection_dCasTemp,              & ! output: derivative in stability correction w.r.t. canopy air space temperature (K-1)
                  dGroundStabilityCorrection_dSfcTemp,              & ! output: derivative in stability correction w.r.t. surface temperature (K-1)
                  err, cmessage                                     ) ! output: error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

  ! compute the ground resistance
  groundResistance = groundResistanceNeutral / groundStabilityCorrection
  if(groundResistance < 0._rkind)then; err=20; message=trim(message)//'ground resistance < 0 [vegetation is present]'; return; end if

 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! * compute resistance for the case without a canopy (bare ground, or canopy completely buried with snow)
 else

  ! no canopy, so set huge resistances (not used)
  canopyResistance = 1.e12_rkind   ! not used: huge resistance, so conductance is essentially zero
  leafResistance   = 1.e12_rkind   ! not used: huge resistance, so conductance is essentially zero

  ! check that measurement height above the ground surface is above the roughness length
  if(mHeight < snowDepth+z0Ground)then; err=20; message=trim(message)//'measurement height < snow depth + roughness length'; return; end if

  ! compute the resistance between the surface and canopy air UNDER NEUTRAL CONDITIONS (s m-1)
  groundExNeut = (vkc**2._rkind) / ( log((mHeight - snowDepth)/z0Ground)**2._rkind) ! turbulent transfer coefficient under conditions of neutral stability (-)
  groundResistanceNeutral = 1._rkind / (groundExNeut*windspd)

  ! define height above the snow surface
  heightAboveGround  = mHeight - snowDepth

  ! check that measurement height above the ground surface is above the roughness length
  if(heightAboveGround < z0Ground)then
   print*, 'z0Ground = ', z0Ground
   print*, 'mHeight  = ', mHeight
   print*, 'snowDepth = ', snowDepth
   print*, 'heightAboveGround = ', heightAboveGround
   message=trim(message)//'height above ground < roughness length [likely due to snow accumulation]'
   err=20; return
  end if

  ! compute ground stability correction
  call aStability(&
                   ! input
                  derivDesired,                                     & ! input: logical flag to compute analytical derivatives
                  ixStability,                                      & ! input: choice of stability function
                  ! input: forcing data, diagnostic and state variables
                  heightAboveGround,                                & ! input: measurement height above the ground surface (m)
                  airtemp,                                          & ! input: temperature above the ground surface (K)
                  groundTemp,                                       & ! input: trial value of surface temperature -- "surface" is either canopy or ground (K)
                  windspd,                                          & ! input: wind speed above the ground surface (m s-1)
                  ! input: stability parameters
                  critRichNumber,                                   & ! input: critical value for the bulk Richardson number where turbulence ceases (-)
                  Louis79_bparam,                                   & ! input: parameter in Louis (1979) stability function
                  Mahrt87_eScale,                                   & ! input: exponential scaling factor in the Mahrt (1987) stability function
                  ! output
                  RiBulkGround,                                     & ! output: bulk Richardson number (-)
                  groundStabilityCorrection,                        & ! output: stability correction for turbulent heat fluxes (-)
                  dGroundStabilityCorrection_dRich,                 & ! output: derivative in stability correction w.r.t. Richardson number for the ground surface (-)
                  dGroundStabilityCorrection_dAirTemp,              & ! output: (not used) derivative in stability correction w.r.t. air temperature (K-1)
                  dGroundStabilityCorrection_dSfcTemp,              & ! output: derivative in stability correction w.r.t. surface temperature (K-1)
                  err, cmessage                                     ) ! output: error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

  ! compute the ground resistance (after stability corrections)
  groundResistance = groundResistanceNeutral/groundStabilityCorrection
  if(groundResistance < 0._rkind)then; err=20; message=trim(message)//'ground resistance < 0 [no vegetation]'; return; end if

  ! set all canopy variables to missing (no canopy!)
  z0Canopy                   = missingValue   ! roughness length of the vegetation canopy (m)
  RiBulkCanopy               = missingValue   ! bulk Richardson number for the canopy (-)
  windReductionFactor        = missingValue   ! canopy wind reduction factor (-)
  zeroPlaneDisplacement      = missingValue   ! zero plane displacement (m)
  canopyStabilityCorrection  = missingValue   ! stability correction for the canopy (-)
  eddyDiffusCanopyTop        = missingValue   ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
  frictionVelocity           = missingValue   ! friction velocity (m s-1)
  windspdCanopyTop           = missingValue   ! windspeed at the top of the canopy (m s-1)
  windspdCanopyBottom        = missingValue   ! windspeed at the height of the bottom of the canopy (m s-1)

 end if  ! (if no canopy)

 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! * compute derivatives
 if(derivDesired)then  ! if analytical derivatives are desired

  ! derivatives for the vegetation canopy
  if(computeVegFlux) then ! (if vegetation is exposed)

   ! ***** compute derivatives w.r.t. canopy temperature
   ! NOTE: derivatives are zero because using canopy air space temperature
   dCanopyResistance_dTCanopy = 0._rkind ! derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
   dGroundResistance_dTCanopy = 0._rkind ! derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)

   ! ***** compute derivatives w.r.t. ground temperature (s m-1 K-1)
   dGroundResistance_dTGround = -(groundResistanceNeutral*dGroundStabilityCorrection_dSfcTemp)/(groundStabilityCorrection**2._rkind)

   ! ***** compute derivatives w.r.t. temperature of the canopy air space (s m-1 K-1)
   ! derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
   dCanopyResistance_dTCanair = -dCanopyStabilityCorrection_dCasTemp/(windspd*canopyExNeut*canopyStabilityCorrection**2._rkind)
   ! derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
   ! (compute derivative in NEUTRAL ground resistance w.r.t. canopy air temperature (s m-1 K-1))
   dFV_dT = windspd*canopyExNeut*dCanopyStabilityCorrection_dCasTemp/(sqrt(sfc2AtmExchangeCoeff_canopy)*2._rkind)                         ! d(frictionVelocity)/d(canopy air temperature)
   dED_dT = dFV_dT*vkc*(heightCanopyTopAboveSnow - zeroPlaneDisplacement)                                                              ! d(eddyDiffusCanopyTop)d(canopy air temperature)
   dGR_dT = -dED_dT*(tmp1 - tmp2)*heightCanopyTopAboveSnow*exp(windReductionFactor) / (windReductionFactor*eddyDiffusCanopyTop**2._rkind) ! d(groundResistanceNeutral)/d(canopy air temperature)
   ! (stitch everything together -- product rule)
   dGroundResistance_dTCanair = dGR_dT/groundStabilityCorrection - groundResistanceNeutral*dGroundStabilityCorrection_dCasTemp/(groundStabilityCorrection**2._rkind)

  ! ***** compute resistances for non-vegetated surfaces (e.g., snow)
  else

   ! set canopy derivatives to zero (non-vegetated, remember)
   dCanopyResistance_dTCanopy = 0._rkind
   dGroundResistance_dTCanopy = 0._rkind

   ! compute derivatives for ground resistance
   dGroundResistance_dTGround = -dGroundStabilityCorrection_dSfcTemp/(windspd*groundExNeut*groundStabilityCorrection**2._rkind)

  end if  ! (switch between vegetated and non-vegetated surfaces)

 ! * analytical derivatives not desired
 else
  dGroundResistance_dTGround = missingValue
  dGroundResistance_dTCanopy = missingValue
  dCanopyResistance_dTCanopy = missingValue
 end if

 ! test
 !print*, 'dGroundResistance_dTGround = ', dGroundResistance_dTGround
 !print*, 'dGroundResistance_dTCanopy = ', dGroundResistance_dTCanopy
 !print*, 'dCanopyResistance_dTCanopy = ', dCanopyResistance_dTCanopy
 !pause 'in aeroResist'

 end subroutine aeroResist


 ! *******************************************************************************************************
 ! private subroutine soilResist: compute soil moisture factor controlling stomatal resistance
 ! *******************************************************************************************************
 subroutine soilResist(&
                       ! input (model decisions)
                       ixSoilResist,             & ! intent(in): choice of function for the soil moisture control on stomatal resistance
                       ixGroundwater,            & ! intent(in): choice of groundwater representation
                       ! input (state variables)
                       mLayerMatricHead,         & ! intent(in): matric head in each layer (m)
                       mLayerVolFracLiq,         & ! intent(in): volumetric fraction of liquid water in each layer
                       scalarAquiferStorage,     & ! intent(in): aquifer storage (m)
                       ! input (diagnostic variables)
                       mLayerRootDensity,        & ! intent(in): root density in each layer (-)
                       scalarAquiferRootFrac,    & ! intent(in): fraction of roots below the lowest unsaturated layer (-)
                       ! input (parameters)
                       plantWiltPsi,             & ! intent(in): matric head at wilting point (m)
                       soilStressParam,          & ! intent(in): parameter in the exponential soil stress function (-)
                       critSoilWilting,          & ! intent(in): critical vol. liq. water content when plants are wilting (-)
                       critSoilTranspire,        & ! intent(in): critical vol. liq. water content when transpiration is limited (-)
                       critAquiferTranspire,     & ! intent(in): critical aquifer storage value when transpiration is limited (m)
                       ! output
                       wAvgTranspireLimitFac,    & ! intent(out): weighted average of the transpiration limiting factor (-)
                       mLayerTranspireLimitFac,  & ! intent(out): transpiration limiting factor in each layer (-)
                       aquiferTranspireLimitFac, & ! intent(out): transpiration limiting factor for the aquifer (-)
                       err,message)                ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 USE mDecisions_module, only: NoahType,CLM_Type,SiB_Type  ! options for the choice of function for the soil moisture control on stomatal resistance
 USE mDecisions_module, only: bigBucket                   ! named variable that defines the "bigBucket" groundwater parameterization
 implicit none
 ! input (model decisions)
 integer(i4b),intent(in)       :: ixSoilResist             ! choice of function for the soil moisture control on stomatal resistance
 integer(i4b),intent(in)       :: ixGroundwater            ! choice of groundwater representation
 ! input (variables)
 real(rkind),intent(in)           :: mLayerMatricHead(:)      ! matric head in each layer (m)
 real(rkind),intent(in)           :: mLayerVolFracLiq(:)      ! volumetric fraction of liquid water in each layer (-)
 real(rkind),intent(in)           :: scalarAquiferStorage     ! aquifer storage (m)
 ! input (diagnostic variables)
 real(rkind),intent(in)           :: mLayerRootDensity(:)     ! root density in each layer (-)
 real(rkind),intent(in)           :: scalarAquiferRootFrac    ! fraction of roots below the lowest unsaturated layer (-)
 ! input (parameters)
 real(rkind),intent(in)           :: plantWiltPsi             ! matric head at wilting point (m)
 real(rkind),intent(in)           :: soilStressParam          ! parameter in the exponential soil stress function (-)
 real(rkind),intent(in)           :: critSoilWilting          ! critical vol. liq. water content when plants are wilting (-)
 real(rkind),intent(in)           :: critSoilTranspire        ! critical vol. liq. water content when transpiration is limited (-)
 real(rkind),intent(in)           :: critAquiferTranspire     ! critical aquifer storage value when transpiration is limited (m)
 ! output
 real(rkind),intent(out)          :: wAvgTranspireLimitFac    ! intent(out): weighted average of the transpiration limiting factor (-)
 real(rkind),intent(out)          :: mLayerTranspireLimitFac(:)  ! intent(out): transpiration limiting factor in each layer (-)
 real(rkind),intent(out)          :: aquiferTranspireLimitFac ! intent(out): transpiration limiting factor for the aquifer (-)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local variables
 real(rkind)                      :: gx                       ! stress function for the soil layers
 real(rkind),parameter            :: verySmall=epsilon(gx)    ! a very small number
 integer(i4b)                  :: iLayer                   ! index of soil layer
 ! initialize error control
 err=0; message='soilResist/'

 ! ** compute the factor limiting transpiration for each soil layer (-)
 wAvgTranspireLimitFac = 0._rkind  ! (initialize the weighted average)
 do iLayer=1,size(mLayerMatricHead)
  ! compute the soil stress function
  select case(ixSoilResist)
   case(NoahType)  ! thresholded linear function of volumetric liquid water content
    gx = (mLayerVolFracLiq(iLayer) - critSoilWilting) / (critSoilTranspire - critSoilWilting)
   case(CLM_Type)  ! thresholded linear function of matric head
    if(mLayerMatricHead(iLayer) > plantWiltPsi)then
     gx = 1._rkind - mLayerMatricHead(iLayer)/plantWiltPsi
    else
     gx = 0._rkind
    end if
   case(SiB_Type)  ! exponential of the log of matric head
    if(mLayerMatricHead(iLayer) < 0._rkind)then  ! (unsaturated)
     gx = 1._rkind - exp( -soilStressParam * ( log(plantWiltPsi/mLayerMatricHead(iLayer)) ) )
    else ! (saturated)
     gx = 1._rkind
    end if
   case default    ! check identified the option
    err=20; message=trim(message)//'cannot identify option for soil resistance'; return
  end select
  ! save the factor for the given layer (ensure between zero and one)
  mLayerTranspireLimitFac(iLayer) = min( max(verySmall,gx), 1._rkind)
  ! compute the weighted average (weighted by root density)
  wAvgTranspireLimitFac = wAvgTranspireLimitFac + mLayerTranspireLimitFac(iLayer)*mLayerRootDensity(iLayer)
 end do ! (looping through soil layers)

 ! ** compute the factor limiting evaporation in the aquifer
 if(scalarAquiferRootFrac > verySmall)then
  ! check that aquifer root fraction is allowed
  if(ixGroundwater /= bigBucket)then
   message=trim(message)//'aquifer evaporation only allowed for the big groundwater bucket -- increase the soil depth to account for roots'
   err=20; return
  end if
  ! compute the factor limiting evaporation for the aquifer
  aquiferTranspireLimitFac = min(scalarAquiferStorage/critAquiferTranspire, 1._rkind)
 else  ! (if there are roots in the aquifer)
  aquiferTranspireLimitFac = 0._rkind
 end if

 ! compute the weighted average (weighted by root density)
 wAvgTranspireLimitFac = wAvgTranspireLimitFac + aquiferTranspireLimitFac*scalarAquiferRootFrac

 end subroutine soilResist


 ! ********************************************************************************
 ! private subroutine turbFluxes: compute turbulent heat fluxes
 ! ********************************************************************************
 subroutine turbFluxes(&
                       ! input: model control
                       computeVegFlux,                & ! intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                       ixDerivMethod,                 & ! intent(in): choice of method used to compute derivative (analytical or numerical)
                       ! input: above-canopy forcing data
                       airtemp,                       & ! intent(in): air temperature at some height above the surface (K)
                       airpres,                       & ! intent(in): air pressure of the air above the vegetation canopy (Pa)
                       VPair,                         & ! intent(in): vapor pressure of the air above the vegetation canopy (Pa)
                       ! input: latent heat of sublimation/vaporization
                       latHeatSubVapCanopy,           & ! intent(in): latent heat of sublimation/vaporization for the vegetation canopy (J kg-1)
                       latHeatSubVapGround,           & ! intent(in): latent heat of sublimation/vaporization for the ground surface (J kg-1)
                       ! input: canopy and ground temperature
                       canairTemp,                    & ! intent(in): temperature of the canopy air space (K)
                       canopyTemp,                    & ! intent(in): canopy temperature (K)
                       groundTemp,                    & ! intent(in): ground temperature (K)
                       satVP_CanopyTemp,              & ! intent(in): saturation vapor pressure at the temperature of the veg canopy (Pa)
                       satVP_GroundTemp,              & ! intent(in): saturation vapor pressure at the temperature of the ground (Pa)
                       dSVPCanopy_dCanopyTemp,        & ! intent(in): derivative in canopy saturation vapor pressure w.r.t. canopy temperature (Pa K-1)
                       dSVPGround_dGroundTemp,        & ! intent(in): derivative in ground saturation vapor pressure w.r.t. ground temperature (Pa K-1)
                       ! input: diagnostic variables
                       exposedVAI,                    & ! intent(in): exposed vegetation area index -- leaf plus stem (m2 m-2)
                       canopyWetFraction,             & ! intent(in): fraction of canopy that is wet [0-1]
                       dCanopyWetFraction_dWat,       & ! intent(in): derivative in the canopy wetted fraction w.r.t. total water content (kg-1 m-2)
                       dCanopyWetFraction_dT,         & ! intent(in): derivative in wetted fraction w.r.t. canopy temperature (K-1)
                       canopySunlitLAI,               & ! intent(in): sunlit leaf area (-)
                       canopyShadedLAI,               & ! intent(in): shaded leaf area (-)
                       soilRelHumidity,               & ! intent(in): relative humidity in the soil pores [0-1]
                       soilResistance,                & ! intent(in): resistance from the soil (s m-1)
                       leafResistance,                & ! intent(in): mean leaf boundary layer resistance per unit leaf area (s m-1)
                       groundResistance,              & ! intent(in): below canopy aerodynamic resistance (s m-1)
                       canopyResistance,              & ! intent(in): above canopy aerodynamic resistance (s m-1)
                       stomResistSunlit,              & ! intent(in): stomatal resistance for sunlit leaves (s m-1)
                       stomResistShaded,              & ! intent(in): stomatal resistance for shaded leaves (s m-1)
                       ! input: derivatives in scalar resistances
                       dGroundResistance_dTGround,    & ! intent(in): derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
                       dGroundResistance_dTCanopy,    & ! intent(in): derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
                       dGroundResistance_dTCanair,    & ! intent(in): derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
                       dCanopyResistance_dTCanopy,    & ! intent(in): derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
                       dCanopyResistance_dTCanair,    & ! intent(in): derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
                       ! output: conductances (used to check derivative calculations)
                       leafConductance,               & ! intent(out): leaf conductance (m s-1)
                       canopyConductance,             & ! intent(out): canopy conductance (m s-1)
                       groundConductanceSH,           & ! intent(out): ground conductance for sensible heat (m s-1)
                       groundConductanceLH,           & ! intent(out): ground conductance for latent heat -- includes soil resistance (m s-1)
                       evapConductance,               & ! intent(out): conductance for evaporation (m s-1)
                       transConductance,              & ! intent(out): conductance for transpiration (m s-1)
                       totalConductanceSH,            & ! intent(out): total conductance for sensible heat (m s-1)
                       totalConductanceLH,            & ! intent(out): total conductance for latent heat (m s-1)
                       ! output: canopy air space variables
                       VP_CanopyAir,                  & ! intent(out): vapor pressure of the canopy air space (Pa)
                       ! output: fluxes from the vegetation canopy
                       senHeatCanopy,                 & ! intent(out): sensible heat flux from the canopy to the canopy air space (W m-2)
                       latHeatCanopyEvap,             & ! intent(out): latent heat flux associated with evaporation from the canopy to the canopy air space (W m-2)
                       latHeatCanopyTrans,            & ! intent(out): latent heat flux associated with transpiration from the canopy to the canopy air space (W m-2)
                       ! output: fluxes from non-vegetated surfaces (ground surface below vegetation, bare ground, or snow covered vegetation)
                       senHeatGround,                 & ! intent(out): sensible heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
                       latHeatGround,                 & ! intent(out): latent heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
                       ! output: total heat fluxes to the atmosphere
                       senHeatTotal,                  & ! intent(out): total sensible heat flux to the atmosphere (W m-2)
                       latHeatTotal,                  & ! intent(out): total latent heat flux to the atmosphere (W m-2)
                       ! output: net fluxes
                       turbFluxCanair,                & ! intent(out): net turbulent heat fluxes at the canopy air space (W m-2)
                       turbFluxCanopy,                & ! intent(out): net turbulent heat fluxes at the canopy (W m-2)
                       turbFluxGround,                & ! intent(out): net turbulent heat fluxes at the ground surface (W m-2)
                       ! output: flux derivatives
                       dTurbFluxCanair_dTCanair,      & ! intent(out): derivative in net canopy air space fluxes w.r.t. canopy air temperature (W m-2 K-1)
                       dTurbFluxCanair_dTCanopy,      & ! intent(out): derivative in net canopy air space fluxes w.r.t. canopy temperature (W m-2 K-1)
                       dTurbFluxCanair_dTGround,      & ! intent(out): derivative in net canopy air space fluxes w.r.t. ground temperature (W m-2 K-1)
                       dTurbFluxCanopy_dTCanair,      & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
                       dTurbFluxCanopy_dTCanopy,      & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
                       dTurbFluxCanopy_dTGround,      & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
                       dTurbFluxGround_dTCanair,      & ! intent(out): derivative in net ground turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
                       dTurbFluxGround_dTCanopy,      & ! intent(out): derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
                       dTurbFluxGround_dTGround,      & ! intent(out): derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
                       ! output: liquid flux derivatives (canopy evap)
                       dLatHeatCanopyEvap_dCanLiq,    & ! intent(out): derivative in latent heat of canopy evaporation w.r.t. canopy liquid water content (J kg-1 s-1)
                       dLatHeatCanopyEvap_dTCanair,   & ! intent(out): derivative in latent heat of canopy evaporation w.r.t. canopy air temperature (W m-2 K-1)
                       dLatHeatCanopyEvap_dTCanopy,   & ! intent(out): derivative in latent heat of canopy evaporation w.r.t. canopy temperature (W m-2 K-1)
                       dLatHeatCanopyEvap_dTGround,   & ! intent(out): derivative in latent heat of canopy evaporation w.r.t. ground temperature (W m-2 K-1)
                       ! output: liquid flux derivatives (ground evap)
                       dLatHeatGroundEvap_dCanLiq,    & ! intent(out): derivative in latent heat of ground evaporation w.r.t. canopy liquid water content (J kg-1 s-1)
                       dLatHeatGroundEvap_dTCanair,   & ! intent(out): derivative in latent heat of ground evaporation w.r.t. canopy air temperature
                       dLatHeatGroundEvap_dTCanopy,   & ! intent(out): derivative in latent heat of ground evaporation w.r.t. canopy temperature
                       dLatHeatGroundEvap_dTGround,   & ! intent(out): derivative in latent heat of ground evaporation w.r.t. ground temperature
                       ! output: cross derivatives
                       dTurbFluxCanair_dCanLiq,       & ! intent(out): derivative in net canopy air space fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                       dTurbFluxCanopy_dCanLiq,       & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                       dTurbFluxGround_dCanLiq,       & ! intent(out): derivative in net ground turbulent fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                       ! output: error control
                       err,message                    ) ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input: model control
 logical(lgt),intent(in)       :: computeVegFlux        ! logical flag to compute vegetation fluxes (.false. if veg buried by snow)
 integer(i4b),intent(in)       :: ixDerivMethod         ! choice of method used to compute derivative (analytical or numerical)
 ! input: above-canopy forcing data
 real(rkind),intent(in)           :: airtemp               ! air temperature at some height above the surface (K)
 real(rkind),intent(in)           :: airpres               ! air pressure of the air above the vegetation canopy (Pa)
 real(rkind),intent(in)           :: VPair                 ! vapor pressure of the air above the vegetation canopy (Pa)
 ! input: latent heat of sublimation/vaporization
 real(rkind),intent(in)           :: latHeatSubVapCanopy   ! latent heat of sublimation/vaporization for the vegetation canopy (J kg-1)
 real(rkind),intent(in)           :: latHeatSubVapGround   ! latent heat of sublimation/vaporization for the ground surface (J kg-1)
 ! input: canopy and ground temperature
 real(rkind),intent(in)           :: canairTemp            ! temperature of the canopy air space (K)
 real(rkind),intent(in)           :: canopyTemp            ! canopy temperature (K)
 real(rkind),intent(in)           :: groundTemp            ! ground temperature (K)
 real(rkind),intent(in)           :: satVP_CanopyTemp      ! saturation vapor pressure at the temperature of the veg canopy (Pa)
 real(rkind),intent(in)           :: satVP_GroundTemp      ! saturation vapor pressure at the temperature of the ground (Pa)
 real(rkind),intent(in)           :: dSVPCanopy_dCanopyTemp  ! derivative in canopy saturation vapor pressure w.r.t. canopy temperature (Pa K-1)
 real(rkind),intent(in)           :: dSVPGround_dGroundTemp  ! derivative in ground saturation vapor pressure w.r.t. ground temperature (Pa K-1)
 ! input: diagnostic variables
 real(rkind),intent(in)           :: exposedVAI            ! exposed vegetation area index -- leaf plus stem (m2 m-2)
 real(rkind),intent(in)           :: canopyWetFraction     ! fraction of canopy that is wet [0-1]
 real(rkind),intent(in)           :: dCanopyWetFraction_dWat ! derivative in the canopy wetted fraction w.r.t. liquid water content (kg-1 m-2)
 real(rkind),intent(in)           :: dCanopyWetFraction_dT   ! derivative in the canopy wetted fraction w.r.t. canopy temperature (K-1)
 real(rkind),intent(in)           :: canopySunlitLAI       ! sunlit leaf area (-)
 real(rkind),intent(in)           :: canopyShadedLAI       ! shaded leaf area (-)
 real(rkind),intent(in)           :: soilRelHumidity       ! relative humidity in the soil pores [0-1]
 real(rkind),intent(in)           :: soilResistance        ! resistance from the soil (s m-1)
 real(rkind),intent(in)           :: leafResistance        ! mean leaf boundary layer resistance per unit leaf area (s m-1)
 real(rkind),intent(in)           :: groundResistance      ! below canopy aerodynamic resistance (s m-1)
 real(rkind),intent(in)           :: canopyResistance      ! above canopy aerodynamic resistance (s m-1)
 real(rkind),intent(in)           :: stomResistSunlit      ! stomatal resistance for sunlit leaves (s m-1)
 real(rkind),intent(in)           :: stomResistShaded      ! stomatal resistance for shaded leaves (s m-1)
 ! input: derivatives in scalar resistances
 real(rkind),intent(in)            :: dGroundResistance_dTGround       ! derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
 real(rkind),intent(in)            :: dGroundResistance_dTCanopy       ! derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
 real(rkind),intent(in)            :: dGroundResistance_dTCanair       ! derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
 real(rkind),intent(in)            :: dCanopyResistance_dTCanopy       ! derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
 real(rkind),intent(in)            :: dCanopyResistance_dTCanair       ! derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
 ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! output: conductances -- used to test derivatives
 real(rkind),intent(out)          :: leafConductance              ! leaf conductance (m s-1)
 real(rkind),intent(out)          :: canopyConductance            ! canopy conductance (m s-1)
 real(rkind),intent(out)          :: groundConductanceSH          ! ground conductance for sensible heat (m s-1)
 real(rkind),intent(out)          :: groundConductanceLH          ! ground conductance for latent heat -- includes soil resistance (m s-1)
 real(rkind),intent(out)          :: evapConductance              ! conductance for evaporation (m s-1)
 real(rkind),intent(out)          :: transConductance             ! conductance for transpiration (m s-1)
 real(rkind),intent(out)          :: totalConductanceSH           ! total conductance for sensible heat (m s-1)
 real(rkind),intent(out)          :: totalConductanceLH           ! total conductance for latent heat (m s-1)
 ! output: canopy air space variables
 real(rkind),intent(out)          :: VP_CanopyAir                 ! vapor pressure of the canopy air space (Pa)
 ! output: fluxes from the vegetation canopy
 real(rkind),intent(out)          :: senHeatCanopy                ! sensible heat flux from the canopy to the canopy air space (W m-2)
 real(rkind),intent(out)          :: latHeatCanopyEvap            ! latent heat flux associated with evaporation from the canopy to the canopy air space (W m-2)
 real(rkind),intent(out)          :: latHeatCanopyTrans           ! latent heat flux associated with transpiration from the canopy to the canopy air space (W m-2)
 ! output: fluxes from non-vegetated surfaces (ground surface below vegetation, bare ground, or snow covered vegetation)
 real(rkind),intent(out)          :: senHeatGround                ! sensible heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
 real(rkind),intent(out)          :: latHeatGround                ! latent heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
 ! output: total heat fluxes to the atmosphere
 real(rkind),intent(out)          :: senHeatTotal                 ! total sensible heat flux to the atmosphere (W m-2)
 real(rkind),intent(out)          :: latHeatTotal                 ! total latent heat flux to the atmosphere (W m-2)
 ! output: net fluxes
 real(rkind),intent(out)          :: turbFluxCanair               ! net turbulent heat fluxes at the canopy air space (W m-2)
 real(rkind),intent(out)          :: turbFluxCanopy               ! net turbulent heat fluxes at the canopy (W m-2)
 real(rkind),intent(out)          :: turbFluxGround               ! net turbulent heat fluxes at the ground surface (W m-2)
 ! output: energy flux derivatives
 real(rkind),intent(out)          :: dTurbFluxCanair_dTCanair     ! derivative in net canopy air space fluxes w.r.t. canopy air temperature (W m-2 K-1)
 real(rkind),intent(out)          :: dTurbFluxCanair_dTCanopy     ! derivative in net canopy air space fluxes w.r.t. canopy temperature (W m-2 K-1)
 real(rkind),intent(out)          :: dTurbFluxCanair_dTGround     ! derivative in net canopy air space fluxes w.r.t. ground temperature (W m-2 K-1)
 real(rkind),intent(out)          :: dTurbFluxCanopy_dTCanair     ! derivative in net canopy turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
 real(rkind),intent(out)          :: dTurbFluxCanopy_dTCanopy     ! derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
 real(rkind),intent(out)          :: dTurbFluxCanopy_dTGround     ! derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
 real(rkind),intent(out)          :: dTurbFluxGround_dTCanair     ! derivative in net ground turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
 real(rkind),intent(out)          :: dTurbFluxGround_dTCanopy     ! derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
 real(rkind),intent(out)          :: dTurbFluxGround_dTGround     ! derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
 ! output: liquid flux derivatives (canopy evap)
 real(rkind),intent(out)          :: dLatHeatCanopyEvap_dCanLiq   ! derivative in latent heat of canopy evaporation w.r.t. canopy liquid water content (W kg-1)
 real(rkind),intent(out)          :: dLatHeatCanopyEvap_dTCanair  ! derivative in latent heat of canopy evaporation w.r.t. canopy air temperature (W m-2 K-1)
 real(rkind),intent(out)          :: dLatHeatCanopyEvap_dTCanopy  ! derivative in latent heat of canopy evaporation w.r.t. canopy temperature (W m-2 K-1)
 real(rkind),intent(out)          :: dLatHeatCanopyEvap_dTGround  ! derivative in latent heat of canopy evaporation w.r.t. ground temperature (W m-2 K-1)
 ! output: liquid flux derivatives (ground evap)
 real(rkind),intent(out)          :: dLatHeatGroundEvap_dCanLiq   ! derivative in latent heat of ground evaporation w.r.t. canopy liquid water content (J kg-1 s-1)
 real(rkind),intent(out)          :: dLatHeatGroundEvap_dTCanair  ! derivative in latent heat of ground evaporation w.r.t. canopy air temperature (W m-2 K-1)
 real(rkind),intent(out)          :: dLatHeatGroundEvap_dTCanopy  ! derivative in latent heat of ground evaporation w.r.t. canopy temperature (W m-2 K-1)
 real(rkind),intent(out)          :: dLatHeatGroundEvap_dTGround  ! derivative in latent heat of ground evaporation w.r.t. ground temperature (W m-2 K-1)
 ! output: cross derivatives
 real(rkind),intent(out)          :: dTurbFluxCanair_dCanLiq      ! derivative in net canopy air space fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
 real(rkind),intent(out)          :: dTurbFluxCanopy_dCanLiq      ! derivative in net canopy turbulent fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
 real(rkind),intent(out)          :: dTurbFluxGround_dCanLiq      ! derivative in net ground turbulent fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                          ! error code
 character(*),intent(out)      :: message                      ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! local variables -- general
 real(rkind)                      :: fpart1,fpart2         ! different parts of a function
 real(rkind)                      :: dPart0,dpart1,dpart2         ! derivatives for different parts of a function
 ! local variables -- "constants"
 real(rkind)                      :: volHeatCapacityAir           ! volumetric heat capacity of air (J m-3)
 real(rkind)                      :: latentHeatConstant           ! latent heat constant (kg m-3 K-1)
 ! local variables -- derivatives for energy conductances
 real(rkind)                      :: dEvapCond_dCanopyTemp        ! derivative in evap conductance w.r.t. canopy temperature
 real(rkind)                      :: dTransCond_dCanopyTemp       ! derivative in trans conductance w.r.t. canopy temperature
 real(rkind)                      :: dCanopyCond_dCanairTemp      ! derivative in canopy conductance w.r.t. canopy air temperature
 real(rkind)                      :: dCanopyCond_dCanopyTemp      ! derivative in canopy conductance w.r.t. canopy temperature
 real(rkind)                      :: dGroundCondSH_dCanairTemp    ! derivative in ground conductance of sensible heat w.r.t. canopy air temperature
 real(rkind)                      :: dGroundCondSH_dCanopyTemp    ! derivative in ground conductance of sensible heat w.r.t. canopy temperature
 real(rkind)                      :: dGroundCondSH_dGroundTemp    ! derivative in ground conductance of sensible heat w.r.t. ground temperature
 ! local variables -- derivatives for mass conductances
 real(rkind)                      :: dGroundCondLH_dCanairTemp    ! derivative in ground conductance w.r.t. canopy air temperature
 real(rkind)                      :: dGroundCondLH_dCanopyTemp    ! derivative in ground conductance w.r.t. canopy temperature
 real(rkind)                      :: dGroundCondLH_dGroundTemp    ! derivative in ground conductance w.r.t. ground temperature
 ! local variables -- derivatives for the canopy air space variables
 real(rkind)                      :: fPart_VP                     ! part of the function for vapor pressure of the canopy air space
 real(rkind)                      :: leafConductanceTr            ! leaf conductance for transpiration (m s-1)
 real(rkind)                      :: dVPCanopyAir_dTCanair        ! derivative in the vapor pressure of the canopy air space w.r.t. temperature of the canopy air space
 real(rkind)                      :: dVPCanopyAir_dTCanopy        ! derivative in the vapor pressure of the canopy air space w.r.t. temperature of the canopy
 real(rkind)                      :: dVPCanopyAir_dTGround        ! derivative in the vapor pressure of the canopy air space w.r.t. temperature of the ground
 real(rkind)                      :: dVPCanopyAir_dWetFrac        ! derivative of vapor pressure in the canopy air space w.r.t. wetted fraction of the canopy
 real(rkind)                      :: dVPCanopyAir_dCanLiq         ! derivative of vapor pressure in the canopy air space w.r.t. canopy liquid water content
 ! local variables -- sensible heat flux derivatives
 real(rkind)                      :: dSenHeatTotal_dTCanair       ! derivative in the total sensible heat flux w.r.t. canopy air temperature
 real(rkind)                      :: dSenHeatTotal_dTCanopy       ! derivative in the total sensible heat flux w.r.t. canopy air temperature
 real(rkind)                      :: dSenHeatTotal_dTGround       ! derivative in the total sensible heat flux w.r.t. ground temperature
 real(rkind)                      :: dSenHeatCanopy_dTCanair      ! derivative in the canopy sensible heat flux w.r.t. canopy air temperature
 real(rkind)                      :: dSenHeatCanopy_dTCanopy      ! derivative in the canopy sensible heat flux w.r.t. canopy temperature
 real(rkind)                      :: dSenHeatCanopy_dTGround      ! derivative in the canopy sensible heat flux w.r.t. ground temperature
 real(rkind)                      :: dSenHeatGround_dTCanair      ! derivative in the ground sensible heat flux w.r.t. canopy air temperature
 real(rkind)                      :: dSenHeatGround_dTCanopy      ! derivative in the ground sensible heat flux w.r.t. canopy temperature
 real(rkind)                      :: dSenHeatGround_dTGround      ! derivative in the ground sensible heat flux w.r.t. ground temperature
 ! local variables -- latent heat flux derivatives
 real(rkind)                      :: dLatHeatCanopyTrans_dTCanair ! derivative in the canopy transpiration flux w.r.t. canopy air temperature
 real(rkind)                      :: dLatHeatCanopyTrans_dTCanopy ! derivative in the canopy transpiration flux w.r.t. canopy temperature
 real(rkind)                      :: dLatHeatCanopyTrans_dTGround ! derivative in the canopy transpiration flux w.r.t. ground temperature
 ! local variables -- wetted fraction derivatives
 real(rkind)                      :: dLatHeatCanopyEvap_dWetFrac  ! derivative in the latent heat of canopy evaporation w.r.t. canopy wet fraction (W m-2)
 real(rkind)                      :: dLatHeatCanopyTrans_dWetFrac ! derivative in the latent heat of canopy transpiration w.r.t. canopy wet fraction (W m-2)
 real(rkind)                      :: dLatHeatCanopyTrans_dCanLiq  ! derivative in the latent heat of canopy transpiration w.r.t. canopy liquid water (J kg-1 s-1)
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='turbFluxes/'

 ! compute constants
 volHeatCapacityAir = iden_air*cp_air           ! volumetric heat capacity of air (J m-3)
 latentHeatConstant = iden_air*w_ratio/airpres  ! latent heat constant for (kg m-3 Pa-1)

 ! *****
 ! * compute conductances, and derivatives...
 ! ******************************************

 ! compute conductances for sensible heat (m s-1)
 if(computeVegFlux)then
  leafConductance    = exposedVAI/leafResistance
  leafConductanceTr  = canopySunlitLAI/(leafResistance+stomResistSunlit) + canopyShadedLAI/(leafResistance+stomResistShaded)
  canopyConductance  = 1._rkind/canopyResistance
 else
  leafConductance    = 0._rkind
  canopyConductance  = 0._rkind
 end if
 groundConductanceSH = 1._rkind/groundResistance

 ! compute total conductance for sensible heat
 totalConductanceSH  = leafConductance + groundConductanceSH + canopyConductance

 ! compute conductances for latent heat (m s-1)
 if(computeVegFlux)then
  evapConductance    = canopyWetFraction*leafConductance
  transConductance   = (1._rkind - canopyWetFraction) * leafConductanceTr
  !write(*,'(a,10(f14.8,1x))') 'canopySunlitLAI, canopyShadedLAI, stomResistSunlit, stomResistShaded, leafResistance, canopyWetFraction = ', &
  !                             canopySunlitLAI, canopyShadedLAI, stomResistSunlit, stomResistShaded, leafResistance, canopyWetFraction
 else
  evapConductance    = 0._rkind
  transConductance   = 0._rkind
 end if
 groundConductanceLH = 1._rkind/(groundResistance + soilResistance)  ! NOTE: soilResistance accounts for fractional snow, and =0 when snow cover is 100%
 totalConductanceLH  = evapConductance + transConductance + groundConductanceLH + canopyConductance

 ! check sensible heat conductance
 if(totalConductanceSH < -tinyVal .or. groundConductanceSH < -tinyVal .or. canopyConductance < -tinyVal)then
  message=trim(message)//'negative conductance for sensible heat'
  err=20; return
 endif

 ! check latent heat conductance
 if(totalConductanceLH < tinyVal .or. groundConductanceLH < -tinyVal)then
  message=trim(message)//'negative conductance for latent heat'
  err=20; return
 endif

 ! * compute derivatives
 ! NOTE: it may be more efficient to compute these derivatives when computing resistances
 if(ixDerivMethod == analytical)then

  ! compute derivatives in individual conductances for sensible heat w.r.t. canopy temperature (m s-1 K-1)
  if(computeVegFlux)then
   dEvapCond_dCanopyTemp     = dCanopyWetFraction_dT*leafConductance                       ! derivative in evap conductance w.r.t. canopy temperature
   dTransCond_dCanopyTemp    = -dCanopyWetFraction_dT*leafConductanceTr                    ! derivative in trans conductance w.r.t. canopy temperature
   dCanopyCond_dCanairTemp   = -dCanopyResistance_dTCanair/canopyResistance**2._rkind         ! derivative in canopy conductance w.r.t. canopy air emperature
   dCanopyCond_dCanopyTemp   = -dCanopyResistance_dTCanopy/canopyResistance**2._rkind         ! derivative in canopy conductance w.r.t. canopy temperature
   dGroundCondSH_dCanairTemp = -dGroundResistance_dTCanair/groundResistance**2._rkind         ! derivative in ground conductance w.r.t. canopy air temperature
   dGroundCondSH_dCanopyTemp = -dGroundResistance_dTCanopy/groundResistance**2._rkind         ! derivative in ground conductance w.r.t. canopy temperature
   dGroundCondSH_dGroundTemp = -dGroundResistance_dTGround/groundResistance**2._rkind         ! derivative in ground conductance w.r.t. ground temperature
  else
   dEvapCond_dCanopyTemp     = 0._rkind  ! derivative in evap conductance w.r.t. canopy temperature
   dTransCond_dCanopyTemp    = 0._rkind  ! derivative in trans conductance w.r.t. canopy temperature
   dCanopyCond_dCanairTemp   = 0._rkind  ! derivative in canopy conductance w.r.t. canopy air emperature
   dCanopyCond_dCanopyTemp   = 0._rkind  ! derivative in canopy conductance w.r.t. canopy temperature
   dGroundCondSH_dCanairTemp = 0._rkind  ! derivative in ground conductance w.r.t. canopy air temperature
   dGroundCondSH_dCanopyTemp = 0._rkind  ! derivative in ground conductance w.r.t. canopy temperature
   dGroundCondSH_dGroundTemp = -dGroundResistance_dTGround/groundResistance**2._rkind         ! derivative in ground conductance w.r.t. ground temperature
  end if

  ! compute derivatives in individual conductances for latent heat w.r.t. canopy temperature (m s-1 K-1)
  if(computeVegFlux)then
   dGroundCondLH_dCanairTemp = -dGroundResistance_dTCanair/(groundResistance+soilResistance)**2._rkind ! derivative in ground conductance w.r.t. canopy air temperature
   dGroundCondLH_dCanopyTemp = -dGroundResistance_dTCanopy/(groundResistance+soilResistance)**2._rkind ! derivative in ground conductance w.r.t. canopy temperature
   dGroundCondLH_dGroundTemp = -dGroundResistance_dTGround/(groundResistance+soilResistance)**2._rkind ! derivative in ground conductance w.r.t. ground temperature
  else
   dGroundCondLH_dCanairTemp = 0._rkind  ! derivative in ground conductance w.r.t. canopy air temperature
   dGroundCondLH_dCanopyTemp = 0._rkind  ! derivative in ground conductance w.r.t. canopy temperature
   dGroundCondLH_dGroundTemp = -dGroundResistance_dTGround/(groundResistance+soilResistance)**2._rkind ! derivative in ground conductance w.r.t. ground temperature
  end if

 end if ! (if computing analytical derivatives)

 ! *****
 ! * compute sensible and latent heat fluxes, and derivatives...
 ! *************************************************************

 ! * compute sensible and latent heat fluxes from the canopy to the canopy air space (W m-2)
 if(computeVegFlux)then

  ! compute the vapor pressure in the canopy air space (Pa)
  fPart_VP     = canopyConductance*VPair + (evapConductance + transConductance)*satVP_CanopyTemp + groundConductanceLH*satVP_GroundTemp*soilRelHumidity
  VP_CanopyAir = fPart_VP/totalConductanceLH
  !write(*,'(a,10(f20.10,1x))') 'canopyConductance, evapConductance, transConductance, groundConductanceLH, soilRelHumidity = ', &
  !                              canopyConductance, evapConductance, transConductance, groundConductanceLH, soilRelHumidity

  ! compute sensible heat flux from the canopy air space to the atmosphere
  ! NOTE: canairTemp is a state variable
  senHeatTotal = -volHeatCapacityAir*canopyConductance*(canairTemp - airtemp)
  !print*, 'canairTemp, airtemp, senHeatTotal = ', canairTemp, airtemp, senHeatTotal

  ! compute fluxes
  senHeatCanopy      = -volHeatCapacityAir*leafConductance*(canopyTemp - canairTemp)        ! (positive downwards)
  latHeatCanopyEvap  = -latHeatSubVapCanopy*latentHeatConstant*evapConductance*(satVP_CanopyTemp - VP_CanopyAir)    ! (positive downwards)
  latHeatCanopyTrans =              -LH_vap*latentHeatConstant*transConductance*(satVP_CanopyTemp - VP_CanopyAir)   ! (positive downwards)
  !write(*,'(a,10(f25.15,1x))') 'latHeatCanopyEvap, VP_CanopyAir = ', latHeatCanopyEvap, VP_CanopyAir
  !write(*,'(a,10(f25.15,1x))') 'latHeatCanopyTrans, VP_CanopyAir = ', latHeatCanopyTrans, VP_CanopyAir
  !write(*,'(a,10(f25.15,1x))') 'transConductance = ', transConductance

  ! check that energy for canopy evaporation does not exhaust the available water
  ! NOTE: do this here, rather than enforcing solution constraints, because energy and mass solutions may be uncoupled
  !if(latHeatSubVapCanopy > LH_vap+verySmall)then ! (sublimation)
  ! maxFlux = -canopyIce*LH_sub/dt       ! W m-2
  !else ! (evaporation)
  ! maxFlux = -canopyLiquid*LH_vap/dt    ! W m-2
  !end if
  ! NOTE: fluxes are positive downwards
  !if(latHeatCanopyEvap < maxFlux) latHeatCanopyEvap = maxFlux
  !write(*,'(a,10(f20.10,1x))') 'maxFlux, latHeatCanopyEvap = ', maxFlux, latHeatCanopyEvap

 ! * no vegetation, so fluxes are zero
 else
  senHeatCanopy      = 0._rkind
  latHeatCanopyEvap  = 0._rkind
  latHeatCanopyTrans = 0._rkind
 end if

 ! compute sensible and latent heat fluxes from the ground to the canopy air space (W m-2)
 if(computeVegFlux)then
  senHeatGround      = -volHeatCapacityAir*groundConductanceSH*(groundTemp - canairTemp)                                          ! (positive downwards)
  latHeatGround      = -latHeatSubVapGround*latentHeatConstant*groundConductanceLH*(satVP_GroundTemp*soilRelHumidity - VP_CanopyAir)  ! (positive downwards)
 else
  senHeatGround      = -volHeatCapacityAir*groundConductanceSH*(groundTemp - airtemp)                                                 ! (positive downwards)
  latHeatGround      = -latHeatSubVapGround*latentHeatConstant*groundConductanceLH*(satVP_GroundTemp*soilRelHumidity - VPair)         ! (positive downwards)
  senHeatTotal       = senHeatGround
 end if
 !write(*,'(a,10(f25.15,1x))') 'latHeatGround = ', latHeatGround

 ! compute latent heat flux from the canopy air space to the atmosphere
 ! NOTE: VP_CanopyAir is a diagnostic variable
 latHeatTotal = latHeatCanopyEvap + latHeatCanopyTrans + latHeatGround

 ! * compute derivatives
 if(ixDerivMethod == analytical)then

  ! differentiate CANOPY fluxes
  if(computeVegFlux)then

   ! compute derivatives of vapor pressure in the canopy air space w.r.t. all state variables
   ! (derivative of vapor pressure in the canopy air space w.r.t. temperature of the canopy air space)
   dPart1 = dCanopyCond_dCanairTemp*VPair + dGroundCondLH_dCanairTemp*satVP_GroundTemp*soilRelHumidity
   dPart2 = -(dCanopyCond_dCanairTemp + dGroundCondLH_dCanairTemp)/(totalConductanceLH**2._rkind)
   dVPCanopyAir_dTCanair = dPart1/totalConductanceLH + fPart_VP*dPart2
   ! (derivative of vapor pressure in the canopy air space w.r.t. temperature of the canopy)
   dPart0 = (evapConductance + transConductance)*dSVPCanopy_dCanopyTemp + (dEvapCond_dCanopyTemp + dTransCond_dCanopyTemp)*satVP_CanopyTemp
   dPart1 = dCanopyCond_dCanopyTemp*VPair + dPart0 + dGroundCondLH_dCanopyTemp*satVP_GroundTemp*soilRelHumidity
   dPart2 = -(dCanopyCond_dCanopyTemp + dEvapCond_dCanopyTemp + dTransCond_dCanopyTemp + dGroundCondLH_dCanopyTemp)/(totalConductanceLH**2._rkind)
   dVPCanopyAir_dTCanopy = dPart1/totalConductanceLH + fPart_VP*dPart2
   ! (derivative of vapor pressure in the canopy air space w.r.t. temperature of the ground)
   dPart1 = dGroundCondLH_dGroundTemp*satVP_GroundTemp*soilRelHumidity + groundConductanceLH*dSVPGround_dGroundTemp*soilRelHumidity
   dPart2 = -dGroundCondLH_dGroundTemp/(totalConductanceLH**2._rkind)
   dVPCanopyAir_dTGround = dPart1/totalConductanceLH + fPart_VP*dPart2
   ! (derivative of vapor pressure in the canopy air space w.r.t. wetted fraction of the canopy)
   dPart1 = (leafConductance - leafConductanceTr)*satVP_CanopyTemp
   dPart2 = -(leafConductance - leafConductanceTr)/(totalConductanceLH**2._rkind)
   dVPCanopyAir_dWetFrac = dPart1/totalConductanceLH + fPart_VP*dPart2
   dVPCanopyAir_dCanLiq  = dVPCanopyAir_dWetFrac*dCanopyWetFraction_dWat
   !write(*,'(a,5(f20.8,1x))') 'dVPCanopyAir_dTCanair, dVPCanopyAir_dTCanopy, dVPCanopyAir_dTGround, dVPCanopyAir_dWetFrac, dVPCanopyAir_dCanLiq = ', &
   !                            dVPCanopyAir_dTCanair, dVPCanopyAir_dTCanopy, dVPCanopyAir_dTGround, dVPCanopyAir_dWetFrac, dVPCanopyAir_dCanLiq

   ! sensible heat from the canopy to the atmosphere
   dSenHeatTotal_dTCanair       = -volHeatCapacityAir*canopyConductance - volHeatCapacityAir*dCanopyCond_dCanairTemp*(canairTemp - airtemp)
   dSenHeatTotal_dTCanopy       = -volHeatCapacityAir*dCanopyCond_dCanopyTemp*(canairTemp - airtemp)
   dSenHeatTotal_dTGround       = 0._rkind
   !write(*,'(a,3(f20.8,1x))') 'dSenHeatTotal_dTCanair, dSenHeatTotal_dTCanopy, dSenHeatTotal_dTGround                   = ', &
   !                            dSenHeatTotal_dTCanair, dSenHeatTotal_dTCanopy, dSenHeatTotal_dTGround

   ! sensible heat from the canopy to the canopy air space
   dSenHeatCanopy_dTCanair      =  volHeatCapacityAir*leafConductance
   dSenHeatCanopy_dTCanopy      = -volHeatCapacityAir*leafConductance
   dSenHeatCanopy_dTGround      = 0._rkind
   !write(*,'(a,3(f20.8,1x))') 'dSenHeatCanopy_dTCanair, dSenHeatCanopy_dTCanopy, dSenHeatCanopy_dTGround                = ', &
   !                            dSenHeatCanopy_dTCanair, dSenHeatCanopy_dTCanopy, dSenHeatCanopy_dTGround

   ! sensible heat from the ground to the canopy air space
   dSenHeatGround_dTCanair      = -volHeatCapacityAir*dGroundCondSH_dCanairTemp*(groundTemp - canairTemp) + volHeatCapacityAir*groundConductanceSH
   dSenHeatGround_dTCanopy      = -volHeatCapacityAir*dGroundCondSH_dCanopyTemp*(groundTemp - canairTemp)
   dSenHeatGround_dTGround      = -volHeatCapacityAir*dGroundCondSH_dGroundTemp*(groundTemp - canairTemp) - volHeatCapacityAir*groundConductanceSH
   !write(*,'(a,3(f20.8,1x))') 'dSenHeatGround_dTCanair, dSenHeatGround_dTCanopy, dSenHeatGround_dTGround                = ', &
   !                            dSenHeatGround_dTCanair, dSenHeatGround_dTCanopy, dSenHeatGround_dTGround

   ! latent heat associated with canopy evaporation
   ! (initial calculations)
   fPart1 = -latHeatSubVapCanopy*latentHeatConstant*evapConductance
   dPart1 = -latHeatSubVapCanopy*latentHeatConstant*dEvapCond_dCanopyTemp
   fPart2 = satVP_CanopyTemp - VP_CanopyAir
   dPart2 = dSVPCanopy_dCanopyTemp - dVPCanopyAir_dTCanopy
   ! (derivatives)
   dLatHeatCanopyEvap_dTCanair  = fPart1*(-dVPCanopyAir_dTCanair)
   dLatHeatCanopyEvap_dTCanopy  = fPart1*dpart2 + fPart2*dPart1
   dLatHeatCanopyEvap_dTGround  = fPart1*(-dVPCanopyAir_dTGround)
   !write(*,'(a,3(f20.8,1x))') 'dLatHeatCanopyEvap_dTCanair, dLatHeatCanopyEvap_dTCanopy, dLatHeatCanopyEvap_dTGround    = ', &
   !                            dLatHeatCanopyEvap_dTCanair, dLatHeatCanopyEvap_dTCanopy, dLatHeatCanopyEvap_dTGround

   ! latent heat associated with canopy transpiration
   ! (initial calculations)
   fPart1 = -LH_vap*latentHeatConstant*transConductance
   dPart1 = -LH_vap*latentHeatConstant*dTransCond_dCanopyTemp
   ! (derivatives)
   dLatHeatCanopyTrans_dTCanair = fPart1*(-dVPCanopyAir_dTCanair)
   dLatHeatCanopyTrans_dTCanopy = fPart1*dPart2 + fPart2*dPart1
   dLatHeatCanopyTrans_dTGround = fPart1*(-dVPCanopyAir_dTGround)
   !write(*,'(a,3(f20.8,1x))') 'dLatHeatCanopyTrans_dTCanair, dLatHeatCanopyTrans_dTCanopy, dLatHeatCanopyTrans_dTGround = ', &
   !                            dLatHeatCanopyTrans_dTCanair, dLatHeatCanopyTrans_dTCanopy, dLatHeatCanopyTrans_dTGround

   ! latent heat flux from the ground
   fPart1 = -latHeatSubVapGround*latentHeatConstant*groundConductanceLH       ! function of the first part
   fPart2 = (satVP_GroundTemp*soilRelHumidity - VP_CanopyAir)                 ! function of the second part
   dLatHeatGroundEvap_dTCanair = -latHeatSubVapGround*latentHeatConstant*dGroundCondLH_dCanairTemp*fPart2 - dVPCanopyAir_dTCanair*fPart1
   dLatHeatGroundEvap_dTCanopy = -latHeatSubVapGround*latentHeatConstant*dGroundCondLH_dCanopyTemp*fPart2 - dVPCanopyAir_dTCanopy*fPart1
   dLatHeatGroundEvap_dTGround = -latHeatSubVapGround*latentHeatConstant*dGroundCondLH_dGroundTemp*fPart2 + (dSVPGround_dGroundTemp*soilRelHumidity - dVPCanopyAir_dTGround)*fPart1
   !write(*,'(a,3(f20.8,1x))') 'dLatHeatGroundEvap_dTCanair, dLatHeatGroundEvap_dTCanopy, dLatHeatGroundEvap_dTGround                = ', &
   !                            dLatHeatGroundEvap_dTCanair, dLatHeatGroundEvap_dTCanopy, dLatHeatGroundEvap_dTGround

   ! latent heat associated with canopy evaporation w.r.t. wetted fraction of the canopy
   dPart1 = -latHeatSubVapCanopy*latentHeatConstant*leafConductance
   fPart1 = dPart1*canopyWetFraction
   dLatHeatCanopyEvap_dWetFrac  = dPart1*(satVP_CanopyTemp - VP_CanopyAir) + fPart1*(-dVPCanopyAir_dWetFrac)

   ! latent heat associated with canopy transpiration w.r.t. wetted fraction of the canopy
   dPart1 = LH_vap*latentHeatConstant*leafConductanceTr  ! NOTE: positive, since (1 - wetFrac)
   fPart1 = -dPart1*(1._rkind - canopyWetFraction)
   dLatHeatCanopyTrans_dWetFrac = dPart1*(satVP_CanopyTemp - VP_CanopyAir) + fPart1*(-dVPCanopyAir_dWetFrac)
   !print*, 'dLatHeatCanopyTrans_dWetFrac = ', dLatHeatCanopyTrans_dWetFrac

   ! latent heat associated with canopy transpiration w.r.t. canopy liquid water
   dLatHeatCanopyTrans_dCanLiq = dLatHeatCanopyTrans_dWetFrac*dCanopyWetFraction_dWat ! (J s-1 kg-1)
   !print*, 'dLatHeatCanopyTrans_dCanLiq = ', dLatHeatCanopyTrans_dCanLiq

  else  ! canopy is undefined

   ! set derivatives for canopy fluxes to zero (no canopy, so fluxes are undefined)
   dSenHeatTotal_dTCanair       = 0._rkind
   dSenHeatTotal_dTCanopy       = 0._rkind
   dSenHeatTotal_dTGround       = 0._rkind
   dSenHeatCanopy_dTCanair      = 0._rkind
   dSenHeatCanopy_dTCanopy      = 0._rkind
   dSenHeatCanopy_dTGround      = 0._rkind
   dLatHeatCanopyEvap_dTCanair  = 0._rkind
   dLatHeatCanopyEvap_dTCanopy  = 0._rkind
   dLatHeatCanopyEvap_dTGround  = 0._rkind
   dLatHeatCanopyTrans_dTCanair = 0._rkind
   dLatHeatCanopyTrans_dTCanopy = 0._rkind
   dLatHeatCanopyTrans_dTGround = 0._rkind

   ! set derivatives for wetted area and canopy transpiration to zero (no canopy, so fluxes are undefined)
   dLatHeatCanopyEvap_dWetFrac  = 0._rkind
   dLatHeatCanopyEvap_dCanLiq   = 0._rkind
   dLatHeatCanopyTrans_dCanLiq  = 0._rkind
   dVPCanopyAir_dCanLiq         = 0._rkind

   ! set derivatives for ground fluxes w.r.t canopy temperature to zero (no canopy, so fluxes are undefined)
   dSenHeatGround_dTCanair     = 0._rkind
   dSenHeatGround_dTCanopy     = 0._rkind
   dLatHeatGroundEvap_dTCanair = 0._rkind
   dLatHeatGroundEvap_dTCanopy = 0._rkind

   ! compute derivatives for the ground fluxes w.r.t. ground temperature
   dSenHeatGround_dTGround     = (-volHeatCapacityAir*dGroundCondSH_dGroundTemp)*(groundTemp - airtemp) + &                                               ! d(ground sensible heat flux)/d(ground temp)
                                 (-volHeatCapacityAir*groundConductanceSH)
   dLatHeatGroundEvap_dTGround = (-latHeatSubVapGround*latentHeatConstant*dGroundCondLH_dGroundTemp)*(satVP_GroundTemp*soilRelHumidity - VPair) + &       ! d(ground latent heat flux)/d(ground temp)
                                 (-latHeatSubVapGround*latentHeatConstant*groundConductanceLH)*dSVPGround_dGroundTemp*soilRelHumidity

   !print*, 'dGroundCondLH_dGroundTemp = ', dGroundCondLH_dGroundTemp

  end if   ! (if canopy is defined)

 end if  ! (if computing analytical derivatives)


 ! *****
 ! * compute net turbulent fluxes, and derivatives...
 ! **************************************************

 ! compute net fluxes
 turbFluxCanair = senHeatTotal - senHeatCanopy - senHeatGround            ! net turbulent flux at the canopy air space (W m-2)
 turbFluxCanopy = senHeatCanopy + latHeatCanopyEvap + latHeatCanopyTrans  ! net turbulent flux at the canopy (W m-2)
 turbFluxGround = senHeatGround + latHeatGround                           ! net turbulent flux at the ground surface (W m-2)
 !write(*,'(a,1x,3(f20.10,1x))') 'senHeatCanopy, latHeatCanopyEvap, latHeatCanopyTrans = ', senHeatCanopy, latHeatCanopyEvap, latHeatCanopyTrans

  ! * compute derivatives
 if(ixDerivMethod == analytical)then
  ! (energy derivatives)
  dTurbFluxCanair_dTCanair = dSenHeatTotal_dTCanair - dSenHeatCanopy_dTCanair - dSenHeatGround_dTCanair            ! derivative in net canopy air space fluxes w.r.t. canopy air temperature (W m-2 K-1)
  dTurbFluxCanair_dTCanopy = dSenHeatTotal_dTCanopy - dSenHeatCanopy_dTCanopy - dSenHeatGround_dTCanopy            ! derivative in net canopy air space fluxes w.r.t. canopy temperature (W m-2 K-1)
  dTurbFluxCanair_dTGround = dSenHeatTotal_dTGround - dSenHeatCanopy_dTGround - dSenHeatGround_dTGround            ! derivative in net canopy air space fluxes w.r.t. ground temperature (W m-2 K-1)
  dTurbFluxCanopy_dTCanair = dSenHeatCanopy_dTCanair + dLatHeatCanopyEvap_dTCanair + dLatHeatCanopyTrans_dTCanair  ! derivative in net canopy turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
  dTurbFluxCanopy_dTCanopy = dSenHeatCanopy_dTCanopy + dLatHeatCanopyEvap_dTCanopy + dLatHeatCanopyTrans_dTCanopy  ! derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
  dTurbFluxCanopy_dTGround = dSenHeatCanopy_dTGround + dLatHeatCanopyEvap_dTGround + dLatHeatCanopyTrans_dTGround  ! derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
  dTurbFluxGround_dTCanair = dSenHeatGround_dTCanair + dLatHeatGroundEvap_dTCanair                                 ! derivative in net ground turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
  dTurbFluxGround_dTCanopy = dSenHeatGround_dTCanopy + dLatHeatGroundEvap_dTCanopy                                 ! derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
  dTurbFluxGround_dTGround = dSenHeatGround_dTGround + dLatHeatGroundEvap_dTGround                                 ! derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
  ! (liquid water derivatives)
  dLatHeatCanopyEvap_dCanLiq = dLatHeatCanopyEvap_dWetFrac*dCanopyWetFraction_dWat                                 ! derivative in latent heat of canopy evaporation w.r.t. canopy liquid water (W kg-1)
  dLatHeatGroundEvap_dCanLiq = latHeatSubVapGround*latentHeatConstant*groundConductanceLH*dVPCanopyAir_dCanLiq     ! derivative in latent heat of ground evaporation w.r.t. canopy liquid water (J kg-1 s-1)
  ! (cross deriavtives)
  dTurbFluxCanair_dCanLiq  = 0._rkind                                                                                 ! derivative in net canopy air space fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
  dTurbFluxCanopy_dCanLiq  = dLatHeatCanopyEvap_dCanLiq + dLatHeatCanopyTrans_dCanLiq                              ! derivative in net canopy turbulent fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
  dTurbFluxGround_dCanLiq  = dLatHeatGroundEvap_dCanLiq                                                            ! derivative in net ground turbulent fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
 else ! (just make sure we return something)
  ! (energy derivatives)
  dTurbFluxCanair_dTCanair = 0._rkind
  dTurbFluxCanair_dTCanopy = 0._rkind
  dTurbFluxCanair_dTGround = 0._rkind
  dTurbFluxCanopy_dTCanair = 0._rkind
  dTurbFluxCanopy_dTCanopy = 0._rkind
  dTurbFluxCanopy_dTGround = 0._rkind
  dTurbFluxGround_dTCanair = 0._rkind
  dTurbFluxGround_dTCanopy = 0._rkind
  dTurbFluxGround_dTGround = 0._rkind
  ! (liquid water derivatives)
  dLatHeatCanopyEvap_dCanLiq   = 0._rkind
  dLatHeatGroundEvap_dCanLiq   = 0._rkind
  ! (cross deriavtives)
  dTurbFluxCanair_dCanLiq  = 0._rkind
  dTurbFluxCanopy_dCanLiq  = 0._rkind
  dTurbFluxGround_dCanLiq  = 0._rkind
 end if

 end subroutine turbFluxes


 ! *******************************************************************************************************
 ! private subroutine aStability: compute stability corrections for turbulent heat fluxes (-)
 ! *******************************************************************************************************
 subroutine aStability(&
                       ! input: control
                       computeDerivative,              & ! input: logical flag to compute analytical derivatives
                       ixStability,                    & ! input: choice of stability function
                       ! input: forcing data, diagnostic and state variables
                       mHeight,                        & ! input: measurement height (m)
                       airTemp,                        & ! input: air temperature (K)
                       sfcTemp,                        & ! input: surface temperature (K)
                       windspd,                        & ! input: wind speed (m s-1)
                       ! input: stability parameters
                       critRichNumber,                 & ! input: critical value for the bulk Richardson number where turbulence ceases (-)
                       Louis79_bparam,                 & ! input: parameter in Louis (1979) stability function
                       Mahrt87_eScale,                 & ! input: exponential scaling factor in the Mahrt (1987) stability function
                       ! output
                       RiBulk,                         & ! output: bulk Richardson number (-)
                       stabilityCorrection,            & ! output: stability correction for turbulent heat fluxes (-)
                       dStabilityCorrection_dRich,     & ! output: derivative in stability correction w.r.t. Richardson number (-)
                       dStabilityCorrection_dAirTemp,  & ! output: derivative in stability correction w.r.t. temperature (K-1)
                       dStabilityCorrection_dSfcTemp,  & ! output: derivative in stability correction w.r.t. temperature (K-1)
                       err, message                    ) ! output: error control
 implicit none
 ! input: control
 logical(lgt),intent(in)       :: computeDerivative      ! flag to compute the derivative
 integer(i4b),intent(in)       :: ixStability            ! choice of stability function
 ! input: forcing data, diagnostic and state variables
 real(rkind),intent(in)           :: mHeight                ! measurement height (m)
 real(rkind),intent(in)           :: airtemp                ! air temperature (K)
 real(rkind),intent(in)           :: sfcTemp                ! surface temperature (K)
 real(rkind),intent(in)           :: windspd                ! wind speed (m s-1)
 ! input: stability parameters
 real(rkind),intent(in)           :: critRichNumber         ! critical value for the bulk Richardson number where turbulence ceases (-)
 real(rkind),intent(in)           :: Louis79_bparam         ! parameter in Louis (1979) stability function
 real(rkind),intent(in)           :: Mahrt87_eScale         ! exponential scaling factor in the Mahrt (1987) stability function
 ! output
 real(rkind),intent(out)          :: RiBulk                 ! bulk Richardson number (-)
 real(rkind),intent(out)          :: stabilityCorrection    ! stability correction for turbulent heat fluxes (-)
 real(rkind),intent(out)          :: dStabilityCorrection_dRich    ! derivative in stability correction w.r.t. Richardson number (-)
 real(rkind),intent(out)          :: dStabilityCorrection_dAirTemp ! derivative in stability correction w.r.t. air temperature (K-1)
 real(rkind),intent(out)          :: dStabilityCorrection_dSfcTemp ! derivative in stability correction w.r.t. surface temperature (K-1)
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
 ! local
 real(rkind), parameter           :: verySmall=1.e-10_rkind    ! a very small number (avoid stability of zero)
 real(rkind)                      :: dRiBulk_dAirTemp       ! derivative in the bulk Richardson number w.r.t. air temperature (K-1)
 real(rkind)                      :: dRiBulk_dSfcTemp       ! derivative in the bulk Richardson number w.r.t. surface temperature (K-1)
 real(rkind)                      :: bPrime                 ! scaled "b" parameter for stability calculations in Louis (1979)
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='aStability/'

 ! compute the bulk Richardson number (-)
 call bulkRichardson(&
                     ! input
                     airTemp,                        & ! input: air temperature (K)
                     sfcTemp,                        & ! input: surface temperature (K)
                     windspd,                        & ! input: wind speed (m s-1)
                     mHeight,                        & ! input: measurement height (m)
                     computeDerivative,              & ! input: flag to compute the derivative
                     ! output
                     RiBulk,                         & ! output: bulk Richardson number (-)
                     dRiBulk_dAirTemp,               & ! output: derivative in the bulk Richardson number w.r.t. air temperature (K-1)
                     dRiBulk_dSfcTemp,               & ! output: derivative in the bulk Richardson number w.r.t. surface temperature (K-1)
                     err,message)                      ! output: error control

 ! set derivative to one if not computing it
 if(.not.computeDerivative)then
  dStabilityCorrection_dRich    = 1._rkind
  dStabilityCorrection_dAirTemp = 1._rkind
  dStabilityCorrection_dSfcTemp = 1._rkind
 end if

 ! ***** process unstable cases
 if(RiBulk<0._rkind)then
  ! compute surface-atmosphere exchange coefficient (-)
  stabilityCorrection = (1._rkind - 16._rkind*RiBulk)**0.5_rkind
  ! compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
  if(computeDerivative)then
   dStabilityCorrection_dRich    = (-16._rkind) * 0.5_rkind*(1._rkind - 16._rkind*RiBulk)**(-0.5_rkind)
   dStabilityCorrection_dAirTemp = dRiBulk_dAirTemp * dStabilityCorrection_dRich
   dStabilityCorrection_dSfcTemp = dRiBulk_dSfcTemp * dStabilityCorrection_dRich
  end if
  return
 end if

 ! ***** process stable cases
 select case(ixStability)

  ! ("standard" stability correction, a la Anderson 1976)
  case(standard)
   ! compute surface-atmosphere exchange coefficient (-)
   if(RiBulk <  critRichNumber) stabilityCorrection = (1._rkind - 5._rkind*RiBulk)**2._rkind
   if(RiBulk >= critRichNumber) stabilityCorrection = verySmall
   ! compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
   if(computeDerivative)then
    if(RiBulk <  critRichNumber) dStabilityCorrection_dRich = (-5._rkind) * 2._rkind*(1._rkind - 5._rkind*RiBulk)
    if(RiBulk >= critRichNumber) dStabilityCorrection_dRich = verySmall
   end if

  ! (Louis 1979)
  case(louisInversePower)
   ! scale the "b" parameter for stable conditions
   bprime = Louis79_bparam/2._rkind
   ! compute surface-atmosphere exchange coefficient (-)
   stabilityCorrection = 1._rkind / ( (1._rkind + bprime*RiBulk)**2._rkind )
   if(stabilityCorrection < epsilon(stabilityCorrection)) stabilityCorrection = epsilon(stabilityCorrection)
   ! compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
   if(computeDerivative)then
    dStabilityCorrection_dRich = bprime * (-2._rkind)*(1._rkind + bprime*RiBulk)**(-3._rkind)
   end if

  ! (Mahrt 1987)
  case(mahrtExponential)
   ! compute surface-atmosphere exchange coefficient (-)
   stabilityCorrection = exp(-Mahrt87_eScale * RiBulk)
   if(stabilityCorrection < epsilon(stabilityCorrection)) stabilityCorrection = epsilon(stabilityCorrection)
   ! compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
   if(computeDerivative)then
    dStabilityCorrection_dRich = (-Mahrt87_eScale) * exp(-Mahrt87_eScale * RiBulk)
   end if

  ! (return error if the stability correction method is not found)
  case default
   err=10; message=trim(message)//"optionNotFound[stability correction]"; return

 end select

 ! get the stability correction with respect to air temperature and surface temperature
 ! NOTE: air temperature is used for canopy air temperature, which is a model state variable
 if(computeDerivative)then
  dStabilityCorrection_dAirTemp = dRiBulk_dAirTemp * dStabilityCorrection_dRich
  dStabilityCorrection_dSfcTemp = dRiBulk_dSfcTemp * dStabilityCorrection_dRich
 end if

 end subroutine aStability


 ! *******************************************************************************************************
 ! private subroutine bulkRichardson: compute bulk Richardson number
 ! *******************************************************************************************************
 subroutine bulkRichardson(&
                           ! input
                           airTemp,                    & ! input: air temperature (K)
                           sfcTemp,                    & ! input: surface temperature (K)
                           windspd,                    & ! input: wind speed (m s-1)
                           mHeight,                    & ! input: measurement height (m)
                           computeDerivative,          & ! input: flag to compute the derivative
                           ! output
                           RiBulk,                     & ! output: bulk Richardson number (-)
                           dRiBulk_dAirTemp,           & ! output: derivative in the bulk Richardson number w.r.t. air temperature (K-1)
                           dRiBulk_dSfcTemp,           & ! output: derivative in the bulk Richardson number w.r.t. surface temperature (K-1)
                           err,message)                  ! output: error control
 implicit none
 ! input
 real(rkind),intent(in)           :: airtemp                ! air temperature (K)
 real(rkind),intent(in)           :: sfcTemp                ! surface temperature (K)
 real(rkind),intent(in)           :: windspd                ! wind speed (m s-1)
 real(rkind),intent(in)           :: mHeight                ! measurement height (m)
 logical(lgt),intent(in)       :: computeDerivative      ! flag to compute the derivative
 ! output
 real(rkind),intent(inout)        :: RiBulk                 ! bulk Richardson number (-)
 real(rkind),intent(out)          :: dRiBulk_dAirTemp       ! derivative in the bulk Richardson number w.r.t. air temperature (K-1)
 real(rkind),intent(out)          :: dRiBulk_dSfcTemp       ! derivative in the bulk Richardson number w.r.t. surface temperature (K-1)
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
 ! local variables
 real(rkind)                      :: T_grad        ! gradient in temperature between the atmosphere and surface (K)
 real(rkind)                      :: T_mean        ! mean of the atmosphere and surface temperature (K)
 real(rkind)                      :: RiMult        ! dimensionless scaling factor (-)
 ! initialize error control
 err=0; message='bulkRichardson/'
 ! compute local variables
 T_grad = airtemp - sfcTemp
 T_mean = 0.5_rkind*(airtemp + sfcTemp)
 RiMult = (gravity*mHeight)/(windspd*windspd)
 ! compute the Richardson number
 RiBulk = (T_grad/T_mean) * RiMult
 ! compute the derivative in the Richardson number
 if(computeDerivative)then
  dRiBulk_dAirTemp =  RiMult/T_mean - RiMult*T_grad/(0.5_rkind*((airtemp + sfcTemp)**2._rkind))
  dRiBulk_dSfcTemp = -RiMult/T_mean - RiMult*T_grad/(0.5_rkind*((airtemp + sfcTemp)**2._rkind))
 else
  dRiBulk_dAirTemp = 1._rkind
  dRiBulk_dSfcTemp = 1._rkind
 end if
 end subroutine bulkRichardson


end module vegNrgFlux_module
