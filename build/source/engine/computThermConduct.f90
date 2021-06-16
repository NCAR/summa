
module computThermConduct_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_d,            & ! data vector (rkind)
                    var_ilength,      & ! data vector with variable length dimension (i4b)
                    var_dlength         ! data vector with variable length dimension (rkind)

! named variables defining elements in the data structures
USE var_lookup,only:iLookPARAM,iLookPROG,iLookDIAG,iLookINDEX  ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS               ! named variables for elements of the decision structure

! physical constants
USE multiconst,only:&
                    iden_air,    & ! intrinsic density of air      (kg m-3)
                    iden_ice,    & ! intrinsic density of ice      (kg m-3)
                    iden_water,  & ! intrinsic density of water    (kg m-3)
                    ! specific heat
                    Cp_air,      & ! specific heat of air          (J kg-1 K-1)
                    Cp_ice,      & ! specific heat of ice          (J kg-1 K-1)
                    Cp_soil,     & ! specific heat of soil         (J kg-1 K-1)
                    Cp_water,    & ! specific heat of liquid water (J kg-1 K-1)
                    ! thermal conductivity
                    lambda_air,  & ! thermal conductivity of air   (J s-1 m-1)
                    lambda_ice,  & ! thermal conductivity of ice   (J s-1 m-1)
                    lambda_water   ! thermal conductivity of water (J s-1 m-1)

! missing values
USE globalData,only:integerMissing ! missing integer
USE globalData,only:realMissing    ! missing real number

! named variables that define the layer type
USE globalData,only:iname_snow     ! snow
USE globalData,only:iname_soil     ! soil

! provide access to named variables for thermal conductivity of soil
USE globalData,only:model_decisions  ! model decision structure

! decisions for thermal conductivity of soil
USE mDecisions_module,only:Smirnova2000    ! option for temporally constant thermal conductivity

! decisions for thermal conductivity of soil
USE mDecisions_module,only: funcSoilWet, & ! function of soil wetness
                            mixConstit,  & ! mixture of constituents
                            hanssonVZJ     ! test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004

! privacy
implicit none
private
public::computThermConduct

! algorithmic parameters
real(rkind),parameter     :: valueMissing=-9999._rkind  ! missing value, used when diagnostic or state variables are undefined
real(rkind),parameter     :: verySmall=1.e-6_rkind   ! used as an additive constant to check if substantial difference among real numbers
real(rkind),parameter     :: mpe=1.e-6_rkind         ! prevents overflow error if division by zero
real(rkind),parameter     :: dx=1.e-6_rkind          ! finite difference increment
contains


 ! **********************************************************************************************************
 ! public subroutine computThermConduct: compute diagnostic energy variables (thermal conductivity and heat capacity)
 ! **********************************************************************************************************
 subroutine computThermConduct(&
                       ! input: control variables
                       computeVegFlux,          & ! intent(in): flag to denote if computing the vegetation flux
                       canopyDepth,             & ! intent(in): canopy depth (m)
                       ! input: state variables
                       scalarCanopyIce,           & ! intent(in): canopy ice content (kg m-2)
                       scalarCanopyLiquid,        & ! intent(in): canopy liquid water content (kg m-2)
                       mLayerVolFracIce,          & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                       mLayerVolFracLiq,          & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
                       ! input/output: data structures
                       mpar_data,               & ! intent(in):    model parameters
                       indx_data,               & ! intent(in):    model layer indices
                       prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                       diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                       ! output: error control
                       err,message)               ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! provide access to external subroutines
 USE snow_utils_module,only:tcond_snow            ! compute thermal conductivity of snow
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 logical(lgt),intent(in)         :: computeVegFlux         ! logical flag to denote if computing the vegetation flux
 real(rkind),intent(in)             :: canopyDepth            ! depth of the vegetation canopy (m)
 real(rkind),intent(in)            :: scalarCanopyIce        ! trial value of canopy ice content (kg m-2)
 real(rkind),intent(in)             :: scalarCanopyLiquid
 real(rkind),intent(in)             :: mLayerVolFracLiq(:)        ! trial vector of volumetric liquid water content (-)
 real(rkind),intent(in)             :: mLayerVolFracIce(:)        ! trial vector of volumetric ice water content (-)
 ! input/output: data structures
 type(var_dlength),intent(in)    :: mpar_data              ! model parameters
 type(var_ilength),intent(in)    :: indx_data              ! model layer indices
 type(var_dlength),intent(in)    :: prog_data              ! model prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data              ! model diagnostic variables for a local HRU
 ! output: error control
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)                :: cmessage               ! error message of downwind routine
 integer(i4b)                      :: iLayer                 ! index of model layer
 integer(i4b)                      :: iSoil                  ! index of soil layer
 real(rkind)                          :: TCn                    ! thermal conductivity below the layer interface (W m-1 K-1)
 real(rkind)                          :: TCp                    ! thermal conductivity above the layer interface (W m-1 K-1)
 real(rkind)                          :: zdn                    ! height difference between interface and lower value (m)
 real(rkind)                          :: zdp                    ! height difference between interface and upper value (m)
 real(rkind)                          :: bulkden_soil           ! bulk density of soil (kg m-3)
 real(rkind)                          :: lambda_drysoil         ! thermal conductivity of dry soil (W m-1)
 real(rkind)                          :: lambda_wetsoil         ! thermal conductivity of wet soil (W m-1)
 real(rkind)                          :: lambda_wet             ! thermal conductivity of the wet material
 real(rkind)                          :: relativeSat            ! relative saturation (-)
 real(rkind)                          :: kerstenNum             ! the Kersten number (-), defining weight applied to conductivity of the wet medium
 real(rkind)                          :: den                    ! denominator in the thermal conductivity calculations
 ! local variables to reproduce the thermal conductivity of Hansson et al. VZJ 2005
 real(rkind),parameter                :: c1=0.55_rkind             ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
 real(rkind),parameter                :: c2=0.8_rkind              ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
 real(rkind),parameter                :: c3=3.07_rkind             ! optimized parameter from Hansson et al. VZJ 2005 (-)
 real(rkind),parameter                :: c4=0.13_rkind             ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
 real(rkind),parameter                :: c5=4._rkind               ! optimized parameter from Hansson et al. VZJ 2005 (-)
 real(rkind),parameter                :: f1=13.05_rkind            ! optimized parameter from Hansson et al. VZJ 2005 (-)
 real(rkind),parameter                :: f2=1.06_rkind             ! optimized parameter from Hansson et al. VZJ 2005 (-)
 real(rkind)                          :: fArg,xArg              ! temporary variables (see Hansson et al. VZJ 2005 for details)
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! associate variables in data structure
 associate(&
 ! input: model decisions
 ixThCondSnow            => model_decisions(iLookDECISIONS%thCondSnow)%iDecision,      & ! intent(in): choice of method for thermal conductivity of snow
 ixThCondSoil            => model_decisions(iLookDECISIONS%thCondSoil)%iDecision,      & ! intent(in): choice of method for thermal conductivity of soil
 ! input: coordinate variables
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1),                    & ! intent(in): number of snow layers
 nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1),                    & ! intent(in): number of soil layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1),                  & ! intent(in): total number of layers
 layerType               => indx_data%var(iLookINDEX%layerType)%dat,                   & ! intent(in): layer type (iname_soil or iname_snow)
 mLayerHeight            => prog_data%var(iLookPROG%mLayerHeight)%dat,                 & ! intent(in): height at the mid-point of each layer (m)
 iLayerHeight            => prog_data%var(iLookPROG%iLayerHeight)%dat,                 & ! intent(in): height at the interface of each layer (m)
 ! input: heat capacity and thermal conductivity
 specificHeatVeg         => mpar_data%var(iLookPARAM%specificHeatVeg)%dat(1),          & ! intent(in): specific heat of vegetation (J kg-1 K-1)
 maxMassVegetation       => mpar_data%var(iLookPARAM%maxMassVegetation)%dat(1),        & ! intent(in): maximum mass of vegetation (kg m-2)
 fixedThermalCond_snow   => mpar_data%var(iLookPARAM%fixedThermalCond_snow)%dat(1),    & ! intent(in): temporally constant thermal conductivity of snow (W m-1 K-1)
 ! input: depth varying soil parameters
 iden_soil               => mpar_data%var(iLookPARAM%soil_dens_intr)%dat,              & ! intent(in): intrinsic density of soil (kg m-3)
 thCond_soil             => mpar_data%var(iLookPARAM%thCond_soil)%dat,                 & ! intent(in): thermal conductivity of soil (W m-1 K-1)
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat,                   & ! intent(in): soil porosity (-)
 frac_sand               => mpar_data%var(iLookPARAM%frac_sand)%dat,                   & ! intent(in): fraction of sand (-)
 frac_silt               => mpar_data%var(iLookPARAM%frac_silt)%dat,                   & ! intent(in): fraction of silt (-)
 frac_clay               => mpar_data%var(iLookPARAM%frac_clay)%dat,                   & ! intent(in): fraction of clay (-)
 ! output: diagnostic variables
 mLayerThermalC          => diag_data%var(iLookDIAG%mLayerThermalC)%dat,               & ! intent(out): thermal conductivity at the mid-point of each layer (W m-1 K-1)
 iLayerThermalC          => diag_data%var(iLookDIAG%iLayerThermalC)%dat,               & ! intent(out): thermal conductivity at the interface of each layer (W m-1 K-1)
 mLayerVolFracAir        => diag_data%var(iLookDIAG%mLayerVolFracAir)%dat              & ! intent(out): volumetric fraction of air in each layer (-)
 )  ! end associate statement
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="computThermConduct/"

 ! initialize the soil layer
 iSoil=integerMissing

 ! loop through layers
 do iLayer=1,nLayers

  ! get the soil layer
  if(iLayer>nSnow) iSoil = iLayer-nSnow

  ! compute the thermal conductivity of dry and wet soils (W m-1)
  ! NOTE: this is actually constant over the simulation, and included here for clarity
  if(ixThCondSoil == funcSoilWet .and. layerType(iLayer)==iname_soil)then
   bulkden_soil   = iden_soil(iSoil)*( 1._rkind - theta_sat(iSoil) )
   lambda_drysoil = (0.135_rkind*bulkden_soil + 64.7_rkind) / (iden_soil(iSoil) - 0.947_rkind*bulkden_soil)
   lambda_wetsoil = (8.80_rkind*frac_sand(iSoil) + 2.92_rkind*frac_clay(iSoil)) / (frac_sand(iSoil) + frac_clay(iSoil))
  end if

  ! *****
  ! * compute the volumetric fraction of air in each layer...
  ! *********************************************************
  select case(layerType(iLayer))
   case(iname_soil); mLayerVolFracAir(iLayer) = theta_sat(iSoil) - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer))
   case(iname_snow); mLayerVolFracAir(iLayer) = 1._rkind - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer))
   case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute volumetric fraction of air'; return
  end select
  
  ! *****
  ! * compute the thermal conductivity of snow and soil at the mid-point of each layer...
  ! *************************************************************************************
  select case(layerType(iLayer))

   ! ***** soil
   case(iname_soil)

    ! select option for thermal conductivity of soil
    select case(ixThCondSoil)

     ! ** function of soil wetness
     case(funcSoilWet)

      ! compute the thermal conductivity of the wet material (W m-1)
      lambda_wet  = lambda_wetsoil**( 1._rkind - theta_sat(iSoil) ) * lambda_water**theta_sat(iSoil) * lambda_ice**(theta_sat(iSoil) - mLayerVolFracLiq(iLayer))
      relativeSat = (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer))/theta_sat(iSoil)  ! relative saturation
      ! compute the Kersten number (-)
      if(relativeSat > 0.1_rkind)then ! log10(0.1) = -1
       kerstenNum = log10(relativeSat) + 1._rkind
      else
       kerstenNum = 0._rkind  ! dry thermal conductivity
      endif
      ! ...and, compute the thermal conductivity
      mLayerThermalC(iLayer) = kerstenNum*lambda_wet + (1._rkind - kerstenNum)*lambda_drysoil

     ! ** mixture of constituents
     case(mixConstit)
      mLayerThermalC(iLayer) = thCond_soil(iSoil) * ( 1._rkind - theta_sat(iSoil) ) + & ! soil component
                               lambda_ice         * mLayerVolFracIce(iLayer)     + & ! ice component
                               lambda_water       * mLayerVolFracLiq(iLayer)     + & ! liquid water component
                               lambda_air         * mLayerVolFracAir(iLayer)         ! air component

     ! ** test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004
     case(hanssonVZJ)
      fArg  = 1._rkind + f1*mLayerVolFracIce(iLayer)**f2
      xArg  = mLayerVolFracLiq(iLayer) + fArg*mLayerVolFracIce(iLayer)
      mLayerThermalC(iLayer) = c1 + c2*xArg + (c1 - c4)*exp(-(c3*xArg)**c5)

     ! ** check
     case default; err=20; message=trim(message)//'unable to identify option for thermal conductivity of soil'; return

    end select  ! option for the thermal conductivity of soil

   ! ***** snow
   case(iname_snow)
    ! temporally constant thermal conductivity
    if(ixThCondSnow==Smirnova2000)then
     mLayerThermalC(iLayer) = fixedThermalCond_snow
    ! thermal conductivity as a function of snow density
    else
     call tcond_snow(mLayerVolFracIce(iLayer)*iden_ice,  & ! input: snow density (kg m-3)
                     mLayerThermalC(iLayer),             & ! output: thermal conductivity (W m-1 K-1)
                     err,cmessage)                         ! output: error control
     if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    endif

   ! * error check
   case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute thermal conductivity'; return

  end select
  !print*, 'iLayer, mLayerThermalC(iLayer) = ', iLayer, mLayerThermalC(iLayer)

 end do  ! looping through layers
 !pause

 ! *****
 ! * compute the thermal conductivity of snow at the interface of each layer...
 ! ****************************************************************************
 do iLayer=1,nLayers-1  ! (loop through layers)
  ! get temporary variables
  TCn = mLayerThermalC(iLayer)    ! thermal conductivity below the layer interface (W m-1 K-1)
  TCp = mLayerThermalC(iLayer+1)  ! thermal conductivity above the layer interface (W m-1 K-1)
  zdn = iLayerHeight(iLayer)   - mLayerHeight(iLayer) ! height difference between interface and lower value (m)
  zdp = mLayerHeight(iLayer+1) - iLayerHeight(iLayer) ! height difference between interface and upper value (m)
  den = TCn*zdp + TCp*zdn  ! denominator
  ! compute thermal conductivity
  if(TCn+TCp > epsilon(TCn))then
   iLayerThermalC(iLayer) = (TCn*TCp*(zdn + zdp)) / den
  else
   iLayerThermalC(iLayer) = (TCn*zdn +  TCp*zdp) / (zdn + zdp)
  endif
  !write(*,'(a,1x,i4,1x,10(f9.3,1x))') 'iLayer, TCn, TCp, zdn, zdp, iLayerThermalC(iLayer) = ', iLayer, TCn, TCp, zdn, zdp, iLayerThermalC(iLayer)
 end do  ! looping through layers

 ! special case of hansson
 if(ixThCondSoil==hanssonVZJ)then
  iLayerThermalC(0) = 28._rkind*(0.5_rkind*(iLayerHeight(1) - iLayerHeight(0)))
 else
  iLayerThermalC(0) = mLayerThermalC(1)
 end if

 ! assume the thermal conductivity at the domain boundaries is equal to the thermal conductivity of the layer
 iLayerThermalC(nLayers) = mLayerThermalC(nLayers)

 ! end association to variables in the data structure
 end associate

 end subroutine computThermConduct


end module computThermConduct_module
