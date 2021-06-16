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

module computHeatCap_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_d,            & ! data vector (rkind)
                    var_ilength,      & ! data vector with variable length dimension (i4b)
                    var_dlength         ! data vector with variable length dimension (rkind)

! named variables defining elements in the data structures
USE var_lookup,only:iLookPARAM,iLookDIAG,iLookINDEX  ! named variables for structure elements

! physical constants
USE multiconst,only:&
                    iden_air,    & ! intrinsic density of air      (kg m-3)
                    iden_ice,    & ! intrinsic density of ice      (kg m-3)
                    iden_water,  & ! intrinsic density of water    (kg m-3)
                    ! specific heat
                    Cp_air,      & ! specific heat of air          (J kg-1 K-1)
                    Cp_ice,      & ! specific heat of ice          (J kg-1 K-1)
                    Cp_soil,     & ! specific heat of soil         (J kg-1 K-1)
                    Cp_water       ! specific heat of liquid water (J kg-1 K-1)
! named variables to describe the state variable type
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

! missing values
USE globalData,only:integerMissing ! missing integer

! named variables that define the layer type
USE globalData,only:iname_snow     ! snow
USE globalData,only:iname_soil     ! soil


! privacy
implicit none
private
public::computHeatCap
public::computStatMult

contains


 ! **********************************************************************************************************
 ! public subroutine computHeatCap: compute diagnostic energy variables (heat capacity)
 ! **********************************************************************************************************
 subroutine computHeatCap(&
                        ! input: control variables
                       nLayers,                 	 & ! intent(in): number of layers (soil+snow)
                       computeVegFlux,         		 & ! intent(in): flag to denote if computing the vegetation flux
                       canopyDepth,             	 & ! intent(in): canopy depth (m)
                       ! input data structures
                       mpar_data,               	 & ! intent(in):    model parameters
                       indx_data,               	 & ! intent(in):    model layer indices
                       diag_data,               	 & ! intent(in):    model diagnostic variables for a local HRU
                       ! input: state variables
                       scalarCanopyIce,         	 & ! intent(in)
                       scalarCanopyLiquid,       	 & ! intent(in)
                       scalarCanopyTempTrial,   	 & ! intent(in):  trial value of canopy temperature (K)
                       scalarCanopyTempPrev,    	 & ! intent(in):  previous value of canopy temperature (K)
                       scalarCanopyEnthalpyTrial,    & ! intent(in):  trial enthalpy of the vegetation canopy (J m-3)
                       scalarCanopyEnthalpyPrev,     & ! intent(in):  previous enthalpy of the vegetation canopy (J m-3)
                       mLayerVolFracIce,        	 & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                       mLayerVolFracLiq,        	 & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
                       mLayerTempTrial,          	 & ! intent(in): trial temperature
                       mLayerTempPrev,           	 & ! intent(in): previous temperature
                       mLayerEnthalpyTrial,      	 & ! intent(in): trial enthalpy for snow and soil
                       mLayerEnthalpyPrev,       	 & ! intent(in): previous enthalpy for snow and soil
                       ! output
                       heatCapVeg,               	 &
                       mLayerHeatCap,            	 & ! intent(out): heat capacity for snow and soil
                       ! output: error control
                       err,message)               ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! input: control variables
 logical(lgt),intent(in)         :: computeVegFlux         ! logical flag to denote if computing the vegetation flux
 real(rkind),intent(in)             :: canopyDepth            ! depth of the vegetation canopy (m)
 ! input/output: data structures
 type(var_dlength),intent(in)    :: mpar_data              ! model parameters
 type(var_ilength),intent(in)    :: indx_data              ! model layer indices
 ! input:
 integer(i4b),intent(in)         :: nLayers
 type(var_dlength),intent(in)    :: diag_data              ! diagnostic variables for a local HRU
 real(rkind),intent(in)             :: scalarCanopyIce        ! trial value of canopy ice content (kg m-2)
 real(rkind),intent(in)             :: scalarCanopyLiquid
 real(rkind),intent(in)		     :: scalarCanopyTempTrial  ! trial value of canopy temperature
 real(rkind),intent(in)             :: scalarCanopyEnthalpyTrial ! trial enthalpy of the vegetation canopy (J m-3)
 real(rkind),intent(in)             :: scalarCanopyEnthalpyPrev  ! intent(in):  previous enthalpy of the vegetation canopy (J m-3)
 real(rkind),intent(in)		     :: scalarCanopyTempPrev   ! Previous value of canopy temperature
 real(rkind),intent(in)             :: mLayerVolFracLiq(:)        ! trial vector of volumetric liquid water content (-)
 real(rkind),intent(in)             :: mLayerVolFracIce(:)        ! trial vector of volumetric ice water content (-)
 real(rkind),intent(in)             :: mLayerTempTrial(:)
 real(rkind),intent(in)             :: mLayerTempPrev(:)
 real(rkind),intent(in)             :: mLayerEnthalpyTrial(:)
 real(rkind),intent(in)             :: mLayerEnthalpyPrev(:)
 ! output:
 real(qp),intent(out)            :: heatCapVeg
 real(qp),intent(out)            :: mLayerHeatCap(:)
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: iLayer                 ! index of model layer
 real(rkind)                        :: delT
 real(rkind)                        :: delEnt
 integer(i4b)                    :: iSoil                  ! index of soil layer
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! associate variables in data structure
 associate(&
 ! input: coordinate variables
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1),                    & ! intent(in): number of snow layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1),                  & ! intent(in): total number of layers
 layerType               => indx_data%var(iLookINDEX%layerType)%dat,                   & ! intent(in): layer type (iname_soil or iname_snow)
 ! input: heat capacity and thermal conductivity
 specificHeatVeg         => mpar_data%var(iLookPARAM%specificHeatVeg)%dat(1),          & ! intent(in): specific heat of vegetation (J kg-1 K-1)
 maxMassVegetation       => mpar_data%var(iLookPARAM%maxMassVegetation)%dat(1),        & ! intent(in): maximum mass of vegetation (kg m-2)
 ! input: depth varying soil parameters
 iden_soil               => mpar_data%var(iLookPARAM%soil_dens_intr)%dat,              & ! intent(in): intrinsic density of soil (kg m-3)
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat                    & ! intent(in): soil porosity (-)
 )  ! end associate statemen
 ! initialize error control
 err=0; message="computHeatCap/"
 
 ! initialize the soil layer
 iSoil=integerMissing

 ! compute the bulk volumetric heat capacity of vegetation (J m-3 K-1)
 if(computeVegFlux)then
   delT = scalarCanopyTempTrial - scalarCanopyTempPrev
   if(abs(delT) <= 1e-14_rkind)then   
	  heatCapVeg =  specificHeatVeg*maxMassVegetation/canopyDepth + & ! vegetation component
                    Cp_water*scalarCanopyLiquid/canopyDepth       + & ! liquid water component
                    Cp_ice*scalarCanopyIce/canopyDepth                ! ice component
   else 
   	  delEnt = scalarCanopyEnthalpyTrial - scalarCanopyEnthalpyPrev
   	  heatCapVeg = delEnt / delT
   end if
 end if

 ! loop through layers
 do iLayer=1,nLayers
   delT = mLayerTempTrial(iLayer) - mLayerTempPrev(iLayer)
   if(abs(delT) <= 1e-14_rkind)then
      ! get the soil layer
      if(iLayer>nSnow) iSoil = iLayer-nSnow
         select case(layerType(iLayer))
            ! * soil
            case(iname_soil)
                  mLayerHeatCap(iLayer) =  iden_soil(iSoil)  * Cp_soil  * ( 1._rkind - theta_sat(iSoil) ) + & ! soil component
                                           iden_ice          * Cp_Ice   * mLayerVolFracIce(iLayer)     + & ! ice component
                                           iden_water        * Cp_water * mLayerVolFracLiq(iLayer)     + & ! liquid water component
                                           iden_air          * Cp_air   * ( theta_sat(iSoil) - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer)) )! air component
            case(iname_snow)
                  mLayerHeatCap(iLayer) = iden_ice          * Cp_ice   * mLayerVolFracIce(iLayer)     + & ! ice component
                                          iden_water        * Cp_water * mLayerVolFracLiq(iLayer)     + & ! liquid water component
                                          iden_air   * Cp_air   * ( 1._rkind - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer)) )   ! air component
            case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute olumetric heat capacity'; return
         end select
   else
      delEnt = mLayerEnthalpyTrial(iLayer) - mLayerEnthalpyPrev(iLayer)
      mLayerHeatCap(iLayer) = delEnt / delT
   endif
 end do  ! looping through layers
 
 end associate

 end subroutine computHeatCap
 
 ! **********************************************************************************************************
 ! public subroutine computStatMult: get scale factors
 ! **********************************************************************************************************
 subroutine computStatMult(&
                       heatCapVeg,              &
                       mLayerHeatCap,             &
                       ! input: data structures
                       diag_data,               & ! intent(in):    model diagnostic variables for a local HRU
                       indx_data,               & ! intent(in):    indices defining model states and layers
                       ! output
                       sMul,                    & ! intent(out):   multiplier for state vector (used in the residual calculations)
                       err,message)               ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 USE nr_utility_module,only:arth                   ! get a sequence of numbers arth(start, incr, count)
 USE f2008funcs_module,only:findIndex              ! finds the index of the first value within a vector
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input: data structures
 real(qp),intent(out)            :: heatCapVeg
 real(qp),intent(out)            :: mLayerHeatCap(:)
 type(var_dlength),intent(in)    :: diag_data              ! diagnostic variables for a local HRU
 type(var_ilength),intent(in)    :: indx_data              ! indices defining model states and layers
 ! output: state vectors
 real(qp),intent(out)            :: sMul(:)    ! NOTE: qp  ! multiplier for state vector (used in the residual calculations)
 ! output: error control
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! state subsets
 integer(i4b)                    :: iLayer                 ! index of layer within the snow+soil domain
 integer(i4b)                    :: ixStateSubset          ! index within the state subset
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 associate(&
 ! model diagnostic variables
 canopyDepth         => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in):  [dp]     canopy depth (m)
 volHeatCapVeg       => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),& ! intent(in) : [dp]     bulk volumetric heat capacity of vegetation (J m-3 K-1)
 ! indices defining specific model states
 ixCasNrg            => indx_data%var(iLookINDEX%ixCasNrg)%dat                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy air space energy state variable
 ixVegNrg            => indx_data%var(iLookINDEX%ixVegNrg)%dat                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy energy state variable
 ixVegHyd            => indx_data%var(iLookINDEX%ixVegHyd)%dat                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy hydrology state variable (mass)
 ! vector of energy and hydrology indices for the snow and soil domains
 ixSnowSoilNrg       => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in) : [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
 ixSnowSoilHyd       => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in) : [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
 nSnowSoilNrg        => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in) : [i4b]    number of energy state variables in the snow+soil domain
 nSnowSoilHyd        => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in) : [i4b]    number of hydrology state variables in the snow+soil domain
 ! type of model state variabless
 ixStateType_subset  => indx_data%var(iLookINDEX%ixStateType_subset)%dat       ,& ! intent(in) : [i4b(:)] [state subset] type of desired model state variables
 ! number of layers
 nSnow               => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in) : [i4b]    number of snow layers
 nSoil               => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in) : [i4b]    number of soil layers
 nLayers             => indx_data%var(iLookINDEX%nLayers)%dat(1)                & ! intent(in) : [i4b]    total number of layers
 )  ! end association with variables in the data structures
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='computStatMult/'

 ! -----
 ! * define components of derivative matrices that are constant over a time step (substep)...
 ! ------------------------------------------------------------------------------------------

 ! define the multiplier for the state vector for residual calculations (vegetation canopy)
 ! NOTE: Use the "where" statement to generalize to multiple canopy layers (currently one canopy layer)

 where(ixStateType_subset==iname_nrgCanair) sMul = Cp_air*iden_air   ! volumetric heat capacity of air (J m-3 K-1)
 where(ixStateType_subset==iname_nrgCanopy) sMul = heatCapVeg     ! volumetric heat capacity of the vegetation (J m-3 K-1)
 where(ixStateType_subset==iname_watCanopy) sMul = 1._rkind             ! nothing else on the left hand side
 where(ixStateType_subset==iname_liqCanopy) sMul = 1._rkind             ! nothing else on the left hand side


 ! define the energy multiplier and diagonal elements for the state vector for residual calculations (snow-soil domain)
 if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   ixStateSubset        = ixSnowSoilNrg(iLayer)      ! index within the state vector
   sMul(ixStateSubset)  = mLayerHeatCap(iLayer)        ! transfer volumetric heat capacity to the state multiplier
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! define the hydrology multiplier and diagonal elements for the state vector for residual calculations (snow-soil domain)
 if(nSnowSoilHyd>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   ixStateSubset        = ixSnowSoilHyd(iLayer)      ! index within the state vector
   sMul(ixStateSubset)  = 1._rkind                      ! state multiplier = 1 (nothing else on the left-hand-side)
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! define the scaling factor and diagonal elements for the aquifer
 where(ixStateType_subset==iname_watAquifer)  sMul = 1._rkind

 ! ------------------------------------------------------------------------------------------
 ! ------------------------------------------------------------------------------------------

 end associate
 ! end association to variables in the data structure where vector length does not change
 end subroutine computStatMult
 
end module computHeatCap_module
