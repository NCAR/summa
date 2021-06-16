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

module t2enthalpy_module

! constants
USE multiconst, only: gravity, &                          ! gravitational acceleration (m s-1)
                      Tfreeze, &                          ! freezing point of water (K)
                      Cp_soil,Cp_water,Cp_ice,Cp_air,&    ! specific heat of soil, water and ice (J kg-1 K-1)
                      iden_water,iden_ice,iden_air,&      ! intrinsic density of water and ice (kg m-3)
                      LH_fus                              ! latent heat of fusion (J kg-1)

! data types
USE nrtype
USE data_types,only:var_iLength                    ! var(:)%dat(:)
USE data_types,only:var_dLength                    ! var(:)%dat(:)
USE data_types,only:zLookup                        ! z(:)%var(:)%lookup(:)

! indices within parameter structure
USE var_lookup,only:iLookPARAM                     ! named variables to define structure element
USE var_lookup,only:iLookINDEX                     ! named variables to define structure element
USE var_lookup,only:iLookLOOKUP                    ! named variables to define structure element
USE var_lookup,only:iLookDIAG       ! named variables for structure elements

! data dimensions
USE var_lookup,only:maxvarLookup                   ! maximum number of variables in the lookup tables

! domain types
USE globalData,only:iname_cas                      ! named variables for canopy air space
USE globalData,only:iname_veg                      ! named variables for vegetation canopy
USE globalData,only:iname_snow                     ! named variables for snow
USE globalData,only:iname_soil                     ! named variables for soil
USE globalData,only:iname_aquifer                  ! named variables for the aquifer

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair                ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy                ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_nrgLayer                 ! named variable defining the energy state variable for snow+soil layers

! missing values
USE globalData,only:integerMissing                 ! missing integer
USE globalData,only:realMissing                    ! missing real number

! privacy
implicit none
private
public::T2E_lookup
public::t2enthalpy
public::t2enthalpy_T

! define the look-up table used to compute temperature based on enthalpy
contains


 ! ************************************************************************************************************************
 ! public subroutine T2E_lookup: define a look-up table to compute enthalpy based on temperature
 ! ************************************************************************************************************************
 subroutine T2E_lookup(nSoil,                       &  ! intent(in):    number of soil layers
                       mpar_data,                   &  ! intent(in):    parameter data structure
                       lookup_data,                 &  ! intent(inout): lookup table data structure
                       err,message)
 USE nr_utility_module,only:arth                       ! use to build vectors with regular increments
 USE spline_int_module,only:spline,splint              ! use for cubic spline interpolation
 USE soil_utils_module,only:volFracLiq                 ! use to compute the volumetric fraction of liquid water
 implicit none
 ! declare dummy variables
 integer(i4b),intent(in)       :: nSoil
 type(var_dlength),intent(in)  :: mpar_data            ! model parameters
 type(zLookup),intent(inout)   :: lookup_data          ! lookup tables
 integer(i4b),intent(out)      :: err                  ! error code
 character(*),intent(out)      :: message              ! error message
 ! declare local variables
 character(len=128)            :: cmessage             ! error message in downwind routine
 logical(lgt),parameter        :: doTest=.false.       ! flag to test
 integer(i4b),parameter        :: nLook=100            ! number of elements in the lookup table
 integer(i4b),parameter        :: nIntegr8=10000       ! number of points used in the numerical integration
 real(rkind),parameter            :: T_lower=260.0_rkind     ! lowest temperature value where all liquid water is assumed frozen (K)
 real(rkind),dimension(nLook)     :: xTemp                ! temporary vector
 real(rkind)                      :: xIncr                ! temporary increment
 real(rkind)                      :: T_incr               ! temperature increment
 real(rkind),parameter            :: T_test=272.9742_rkind   ! test value for temperature (K)
 real(rkind)                      :: E_test               ! test value for enthalpy (J m-3)
 integer(i4b)                  :: iVar                 ! loop through variables
 integer(i4b)                  :: iSoil                ! loop through soil layers
 integer(i4b)                  :: iLook                ! loop through lookup table
 integer(i4b)                  :: jIntegr8             ! index for numerical integration
 logical(lgt)                  :: check                ! flag to check allocation
 real(rkind)                      :: vGn_m                ! van Genuchten "m" parameter (-)
 real(rkind)                      :: vFracLiq             ! volumetric fraction of liquid water (-)
 real(rkind)                      :: vFracIce             ! volumetric fraction of ice (-)
 real(rkind)                      :: matricHead           ! matric head (m)
 ! initialize error control
 err=0; message="T2E_lookup/"

 ! get the values of temperature for the lookup table
 xIncr = 1._rkind/real(nLook-1, kind(rkind))
 xTemp = T_lower + (Tfreeze - T_lower)*arth(0._rkind,xIncr,nLook)**0.25_rkind ! use **0.25 to give more values near freezing

 ! -----
 ! * allocate space for the lookup table...
 ! ----------------------------------------

 ! initialize checks
 check=.false.

 ! allocate space for soil layers
 if(allocated(lookup_data%z))then; check=.true.; else; allocate(lookup_data%z(nSoil), stat=err); endif
 if(check) then; err=20; message=trim(message)//'lookup table z dimension was unexpectedly allocated already'; return; end if
 if(err/=0)then; err=20; message=trim(message)//'problem allocating lookup table z dimension dimension'; return; end if

 ! allocate space for the variables in the lookup table
 do iSoil=1,nSoil
  if(allocated(lookup_data%z(iSoil)%var))then; check=.true.; else; allocate(lookup_data%z(iSoil)%var(maxvarLookup), stat=err); endif
  if(check) then; err=20; message=trim(message)//'lookup table var dimension was unexpectedly allocated already'; return; end if
  if(err/=0)then; err=20; message=trim(message)//'problem allocating lookup table var dimension dimension'; return; end if

  ! allocate space for the values in the lookup table
  do iVar=1,maxvarLookup
   if(allocated(lookup_data%z(iSoil)%var(iVar)%lookup))then; check=.true.; else; allocate(lookup_data%z(iSoil)%var(iVar)%lookup(nLook), stat=err); endif
   if(check) then; err=20; message=trim(message)//'lookup table value dimension was unexpectedly allocated already'; return; end if
   if(err/=0)then; err=20; message=trim(message)//'problem allocating lookup table vaule dimension dimension'; return; end if

  end do ! (looping through variables)
 end do ! (looping through soil layers)

 ! loop through soil layers
 do iSoil=1,nSoil

  ! -----
  ! * make association to variables in the data structures...
  ! ---------------------------------------------------------

  associate(&

  ! associate model parameters
  snowfrz_scale  => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)           , & ! scaling parameter for freezing     (K-1)
  soil_dens_intr => mpar_data%var(iLookPARAM%soil_dens_intr)%dat(iSoil)      , & ! intrinsic soil density             (kg m-3)
  theta_sat      => mpar_data%var(iLookPARAM%theta_sat)%dat(iSoil)           , & ! soil porosity                      (-)
  theta_res      => mpar_data%var(iLookPARAM%theta_res)%dat(iSoil)           , & ! volumetric residual water content  (-)
  vGn_alpha      => mpar_data%var(iLookPARAM%vGn_alpha)%dat(iSoil)           , & ! van Genuchten "alpha" parameter    (m-1)
  vGn_n          => mpar_data%var(iLookPARAM%vGn_n)%dat(iSoil)               , & ! van Genuchten "n" parameter        (-)

  ! associate values in the lookup table
  Tk            => lookup_data%z(iSoil)%var(iLookLOOKUP%temperature)%lookup  , & ! temperature (K)
  Ey            => lookup_data%z(iSoil)%var(iLookLOOKUP%enthalpy)%lookup     , & ! enthalpy (J m-3)
  E2            => lookup_data%z(iSoil)%var(iLookLOOKUP%deriv2)%lookup         & ! second derivative of the interpolating function

  ) ! end associate statement

  ! compute vGn_m
  vGn_m = 1._rkind - 1._rkind/vGn_n

  ! -----
  ! * populate the lookup table...
  ! ------------------------------

  ! initialize temperature and enthalpy
  Tk(nLook) = Tfreeze
  Ey(nLook) = 0._rkind

  ! loop through lookup table
  do iLook=(nLook-1),1,-1

   ! update temperature and enthalpy
   Tk(iLook) = Tk(iLook+1)
   Ey(iLook) = Ey(iLook+1)

   ! get the temperature increment for the numerical integration
   T_incr = (xTemp(iLook)-xTemp(iLook+1))/real(nIntegr8, kind(rkind))

   ! numerical integration between different values of the lookup table
   do jIntegr8=1,nIntegr8

    ! update temperature
    Tk(iLook)  = Tk(iLook) + T_incr

    ! compute the volumetric liquid water and ice content
    ! NOTE: assume saturation
    matricHead = (LH_fus/gravity)*(Tk(iLook) - Tfreeze)/Tfreeze
    vFracLiq   = volFracLiq(matricHead,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
    vFracIce   = theta_sat - vFracLiq

    ! compute enthalpy
    ! NOTE: assume intrrinsic density of ice is the intrinsic density of water
    ! NOTE: kg m-3 J kg-1 K-1 K
    Ey(iLook)  = Ey(iLook) + iden_water*Cp_water*vFracLiq*T_incr + iden_water*Cp_ice*vFracIce*T_incr

   end do  ! numerical integration

  end do  ! loop through lookup table

  ! use cubic spline interpolation to obtain enthalpy values at the desired values of temperature
  call spline(Tk,Ey,1.e30_rkind,1.e30_rkind,E2,err,cmessage)  ! get the second derivatives
  if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

  ! check
  !print*, 'Tk = ', Tk
  !print*, 'Ey = ', Ey
  !print*, 'E2 = ', E2

  ! test
  if(doTest)then

   ! calculate enthalpy
   call splint(Tk,Ey,E2,T_test,E_test,err,cmessage)
   if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

   ! write values
   print*, 'doTest    = ', doTest
   print*, 'T_test    = ', T_test    ! temperature (K)
   print*, 'E_test    = ', E_test    ! enthalpy (J m-3)
   print*, 'theta_sat = ', theta_sat ! soil porosity                      (-)
   print*, 'theta_res = ', theta_res ! volumetric residual water content  (-)
   print*, 'vGn_alpha = ', vGn_alpha ! van Genuchten "alpha" parameter    (m-1)
   print*, 'vGn_n     = ', vGn_n     ! van Genuchten "n" parameter        (-)
   print*, trim(message)//'PAUSE: Set doTest=.false. to complete simulations'
   read(*,*)

  endif  ! if testing

  ! end asssociation to variables in the data structures
  end associate

 end do  ! (looping through soil layers)
 end subroutine T2E_lookup


 ! ************************************************************************************************************************
 ! public subroutine t2enthalpy: compute enthalpy from temperature and total water content
 ! ************************************************************************************************************************
 subroutine t2enthalpy(&
                       ! input: data structures
                       diag_data,                         & ! intent(in):  model diagnostic variables for a local HRU
                       mpar_data,                         & ! intent(in):  parameter data structure
                       indx_data,                         & ! intent(in):  model indices
                       lookup_data,                       & ! intent(in):  lookup table data structure
                       ! input: state variables for the vegetation canopy
                       scalarCanairTempTrial,             & ! intent(in):  trial value of canopy air temperature (K)
                       scalarCanopyTempTrial,             & ! intent(in):  trial value of canopy temperature (K)
                       scalarCanopyWatTrial,              & ! intent(in):  trial value of canopy total water (kg m-2)
                       scalarCanopyIceTrial,              & ! intent(in):  trial value of canopy ice content (kg m-2)
                       ! input: variables for the snow-soil domain
                       mLayerTempTrial,                   & ! intent(in):  trial vector of layer temperature (K)
                       mLayerVolFracWatTrial,             & ! intent(in):  trial vector of volumetric total water content (-)
                       mLayerMatricHeadTrial,             & ! intent(in):  trial vector of total water matric potential (m)
                       mLayerVolFracIceTrial,             & ! intent(in)
                       ! output: enthalpy
                       scalarCanairEnthalpy,              & ! intent(out):  enthalpy of the canopy air space (J m-3)
                       scalarCanopyEnthalpy,              & ! intent(out):  enthalpy of the vegetation canopy (J m-3)
                       mLayerEnthalpy,                    & ! intent(out):  enthalpy of each snow+soil layer (J m-3)
                       ! output: error control
                       err,message)                         ! intent(out): error control
 ! -------------------------------------------------------------------------------------------------------------------------
 ! downwind routines
 USE soil_utils_module,only:crit_soilT     ! compute critical temperature below which ice exists
 USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water
 USE spline_int_module,only:splint         ! use for cubic spline interpolation
 implicit none
 ! delare dummy variables
 ! -------------------------------------------------------------------------------------------------------------------------
 ! input: data structures
 type(var_dlength),intent(in)  :: diag_data                 ! diagnostic variables for a local HRU
 type(var_dlength),intent(in)  :: mpar_data                 ! model parameters
 type(var_ilength),intent(in)  :: indx_data                 ! model indices
 type(zLookup),intent(in)      :: lookup_data               ! lookup tables
 ! input: state variables for the vegetation canopy
 real(rkind),intent(in)           :: scalarCanairTempTrial     ! trial value of canopy air temperature (K)
 real(rkind),intent(in)           :: scalarCanopyTempTrial     ! trial value of canopy temperature (K)
 real(rkind),intent(in)           :: scalarCanopyWatTrial      ! trial value of canopy total water (kg m-2)
 real(rkind),intent(in)           :: scalarCanopyIceTrial      ! trial value of canopy ice content (kg m-2)
 ! input: variables for the snow-soil domain
 real(rkind),intent(in)           :: mLayerTempTrial(:)        ! trial vector of layer temperature (K)
 real(rkind),intent(in)           :: mLayerVolFracWatTrial(:)  ! trial vector of volumetric total water content (-)
 real(rkind),intent(in)           :: mLayerMatricHeadTrial(:)  ! trial vector of total water matric potential (m)
 real(rkind),intent(in)           :: mLayerVolFracIceTrial(:)  ! trial vector of volumetric fraction of Ice (-)
 ! output: enthalpy
 real(rkind),intent(out)          :: scalarCanairEnthalpy      ! enthalpy of the canopy air space (J m-3)
 real(rkind),intent(out)          :: scalarCanopyEnthalpy      ! enthalpy of the vegetation canopy (J m-3)
 real(rkind),intent(out)          :: mLayerEnthalpy(:)         ! enthalpy of each snow+soil layer (J m-3)
 ! output: error control
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! -------------------------------------------------------------------------------------------------------------------------
 ! declare local variables
 character(len=128)            :: cmessage                  ! error message in downwind routine
 integer(i4b)                  :: iState                    ! index of model state variable
 integer(i4b)                  :: iLayer                    ! index of model layer
 integer(i4b)                  :: ixFullVector              ! index within full state vector
 integer(i4b)                  :: ixDomainType              ! name of a given model domain
 integer(i4b)                  :: ixControlIndex            ! index within a given model domain
 real(rkind)                      :: vGn_m                     ! van Genuchten "m" parameter (-)
 real(rkind)                      :: Tcrit                     ! temperature where all water is unfrozen (K)
 real(rkind)                      :: psiLiq                    ! matric head of liquid water (m)
 real(rkind)                      :: vFracWat                  ! volumetric fraction of total water, liquid+ice (-)
 real(rkind)                      :: vFracLiq                  ! volumetric fraction of liquid water (-)
 real(rkind)                      :: vFracIce                  ! volumetric fraction of ice (-)
 real(rkind)                      :: enthTemp                  ! enthalpy at the temperature of the control volume (J m-3)
 real(rkind)                      :: enthTcrit                 ! enthalpy at the critical temperature where all water is unfrozen (J m-3)
 real(rkind)                      :: enthPhase                 ! enthalpy associated with phase change (J m-3)
 real(rkind)                      :: enthWater                 ! enthalpy of total water (J m-3)
 real(rkind)                      :: enthSoil                  ! enthalpy of soil particles (J m-3)
 real(rkind)                      :: enthMix                   ! enthalpy of the mixed region, liquid+ice (J m-3)
 real(rkind)                      :: enthLiq                   ! enthalpy of the liquid region (J m-3)
 real(rkind)                      :: enthAir                   ! enthalpy of air (J m-3)
 real(rkind)                      :: diffT
 real(rkind)                      :: integral
 real(rkind)                      :: enthIce
 real(rkind)                      :: enthVeg
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 generalVars: associate(&
 ! number of model layers, and layer type
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):  [i4b]    total number of snow layers
 nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):  [i4b]    total number of soil layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in):  [i4b]    total number of snow and soil layers
 ! mapping between the full state vector and the state subset
 ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat         ,& ! intent(in):  [i4b(:)] list of indices in the state subset for each state in the full state vector
 ixMapSubset2Full        => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat         ,& ! intent(in):  [i4b(:)] [state subset] list of indices of the full state vector in the state subset
 ! type of domain, type of state variable, and index of control volume within domain
 ixDomainType_subset     => indx_data%var(iLookINDEX%ixDomainType_subset)%dat      ,& ! intent(in):  [i4b(:)] [state subset] id of domain for desired model state variables
 ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat          ,& ! intent(in):  [i4b(:)] index of the control volume for different domains (veg, snow, soil)
 ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,& ! intent(in):  [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
 ! snow parameters
 snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)          & ! intent(in):  [dp] scaling parameter for the snow freezing curve (K-1)
 ) ! end associate statement
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! initialize error control
 err=0; message="t2enthalpy/"

 ! loop through model state variables
 do iState=1,size(ixMapSubset2Full)

  ! -----
  ! - compute indices...
  ! --------------------

  ! get domain type, and index of the control volume within the domain
  ixFullVector   = ixMapSubset2Full(iState)       ! index within full state vector
  ixDomainType   = ixDomainType_subset(iState)    ! named variables defining the domain (iname_cas, iname_veg, etc.)
  ixControlIndex = ixControlVolume(ixFullVector)  ! index within a given domain

  ! check an energy state
  if(ixStateType(ixFullVector)==iname_nrgCanair .or. ixStateType(ixFullVector)==iname_nrgCanopy .or. ixStateType(ixFullVector)==iname_nrgLayer)then

   ! get the layer index
   select case(ixDomainType)
    case(iname_cas);     iLayer = integerMissing
    case(iname_veg);     iLayer = integerMissing
    case(iname_snow);    iLayer = ixControlIndex
    case(iname_soil);    iLayer = ixControlIndex + nSnow
    case(iname_aquifer); cycle ! aquifer: do nothing (no thermodynamics in the aquifer)
    case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil, iname_aquifer'; return
   end select

   ! identify domain
   select case(ixDomainType)

    ! -----
    ! - canopy air space...
    ! ---------------------
    case(iname_cas)
       scalarCanairEnthalpy = Cp_air * iden_air * (scalarCanairTempTrial - Tfreeze)
    ! -----
    ! - vegetation canopy...
    ! -----------------------
    case(iname_veg)
      ! association to necessary variables for vegetation
      vegVars: associate(&     
      canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1),         & ! intent(in): [dp]      canopy depth (m)
      specificHeatVeg         => mpar_data%var(iLookPARAM%specificHeatVeg)%dat(1),          & ! intent(in): specific heat of vegetation (J kg-1 K-1)
      maxMassVegetation       => mpar_data%var(iLookPARAM%maxMassVegetation)%dat(1),        & ! intent(in): maximum mass of vegetation (kg m-2)
      snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)   & ! intent(in):  [dp] scaling parameter for the snow freezing curve (K-1)
      ) 

      diffT = scalarCanopyTempTrial - Tfreeze
      enthVeg = specificHeatVeg * maxMassVegetation * diffT / canopyDepth
      if(diffT>=0._rkind)then
         enthLiq = Cp_water * scalarCanopyWatTrial * diffT / canopyDepth
         enthIce = 0._rkind
         enthPhase = 0._rkind
      else
         integral = (1._rkind/snowfrz_scale) * atan(snowfrz_scale * diffT)
         enthLiq = Cp_water * scalarCanopyWatTrial * integral / canopyDepth
         enthIce = Cp_ice * scalarCanopyWatTrial * ( diffT - integral ) / canopyDepth
         enthPhase = LH_fus * scalarCanopyIceTrial / canopyDepth
      end if
      
       scalarCanopyEnthalpy =  enthVeg + enthLiq + enthIce - enthPhase           

      ! end association 
      end associate vegVars

    ! -----
    ! - snow...
    ! ---------
    case(iname_snow)
      
      ! association to necessary variables for snow
      snowVars: associate(&
      snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)   & ! intent(in):  [dp] scaling parameter for the snow freezing curve (K-1)
      ) 

      diffT = mLayerTempTrial(iLayer) - Tfreeze
      if(diffT>=0._rkind)then
         enthLiq = iden_water * Cp_water * mLayerVolFracWatTrial(iLayer) * diffT
         enthIce = 0._rkind
         enthAir = iden_air * Cp_air * ( 1._rkind - mLayerVolFracWatTrial(iLayer) ) * diffT
         enthPhase = 0._rkind
      else
         integral = (1._rkind/snowfrz_scale) * atan(snowfrz_scale * diffT)
         enthLiq = iden_water * Cp_water * mLayerVolFracWatTrial(iLayer) * integral
         enthIce = iden_water * Cp_ice * mLayerVolFracWatTrial(iLayer) * ( diffT - integral )
         enthAir = iden_air * Cp_air * ( diffT - mLayerVolFracWatTrial(iLayer) * ( (iden_water/iden_ice)*(diffT-integral) + integral ) ) 
         enthPhase = iden_ice * LH_fus * mLayerVolFracIceTrial(iLayer)
      end if
      
       mLayerEnthalpy(iLayer) =  enthLiq + enthIce + enthAir - enthPhase           

      ! end association 
      end associate snowVars

    ! -----
    ! - soil...
    ! ---------
    case(iname_soil)

     ! make association to variables in the data structures...
     soilVars: associate(&

     ! associate model parameters
     soil_dens_intr => mpar_data%var(iLookPARAM%soil_dens_intr)%dat(ixControlIndex)      , & ! intrinsic soil density             (kg m-3)
     theta_sat      => mpar_data%var(iLookPARAM%theta_sat)%dat(ixControlIndex)           , & ! soil porosity                      (-)
     theta_res      => mpar_data%var(iLookPARAM%theta_res)%dat(ixControlIndex)           , & ! volumetric residual water content  (-)
     vGn_alpha      => mpar_data%var(iLookPARAM%vGn_alpha)%dat(ixControlIndex)           , & ! van Genuchten "alpha" parameter    (m-1)
     vGn_n          => mpar_data%var(iLookPARAM%vGn_n)%dat(ixControlIndex)               , & ! van Genuchten "n" parameter        (-)

     ! associate values in the lookup table
     Tk            => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%temperature)%lookup  , & ! temperature (K)
     Ey            => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%enthalpy)%lookup     , & ! enthalpy (J m-3)
     E2            => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%deriv2)%lookup         & ! second derivative of the interpolating function

     ) ! end associate statement

     ! diagnostic variables
     vGn_m    = 1._rkind - 1._rkind/vGn_n
     Tcrit    = crit_soilT( mLayerMatricHeadTrial(ixControlIndex) )
     vFracWat = volFracLiq(mLayerMatricHeadTrial(ixControlIndex),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
     

     ! *** compute enthalpy of water for unfrozen conditions
     if(mlayerTempTrial(iLayer) > Tcrit)then
      enthWater = iden_water*Cp_water*vFracWat*(mlayerTempTrial(iLayer) - Tfreeze) ! valid for temperatures below freezing also
      enthPhase = 0._rkind

     ! *** compute enthalpy of water for frozen conditions
     else
      ! calculate enthalpy at the temperature (cubic spline interpolation)
      call splint(Tk,Ey,E2,mlayerTempTrial(iLayer),enthTemp,err,cmessage)
      if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

      ! calculate enthalpy at the critical temperature (cubic spline interpolation)
      call splint(Tk,Ey,E2,Tcrit,enthTcrit,err,cmessage)
      if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

      ! calculate the enthalpy of water
      enthMix   = enthTemp - enthTcrit ! enthalpy of the liquid+ice mix
      enthLiq   = iden_water*Cp_water*vFracWat*(Tcrit - Tfreeze)
      enthWater = enthMix + enthLiq

      ! *** compute the enthalpy associated with phase change
      psiLiq    = (mLayerTempTrial(iLayer) - Tfreeze)*LH_fus/(gravity*Tfreeze)
      vFracLiq  = volFracLiq(psiLiq,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
      vFracIce  = vFracWat - vFracLiq
      enthPhase = iden_water*LH_fus*vFracIce

     endif ! (if frozen conditions)

     ! *** compute the enthalpy of soil
     enthSoil  = soil_dens_intr*Cp_soil*(1._rkind - theta_sat)*(mlayerTempTrial(iLayer) - Tfreeze)

     ! *** compute the enthalpy of air
     enthAir   = iden_air*Cp_air*(1._rkind - theta_sat - vFracWat)*(mlayerTempTrial(iLayer) - Tfreeze)

     ! *** compute the total enthalpy (J m-3)
     mLayerEnthalpy(iLayer) = enthSoil + enthWater + enthAir - enthPhase

!     print*, iLayer, enthSoil, enthWater, enthAir, enthPhase

     end associate soilVars

    ! -----
    ! - checks...
    ! -----------
    case(iname_aquifer); cycle ! aquifer: do nothing
    case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil, iname_aquifer'; return
   end select

  end if  ! if an energy layer
 end do  ! looping through state variables

 end associate generalVars

 end subroutine t2enthalpy
 

 ! ************************************************************************************************************************
 ! public subroutine t2enthalpy_T: compute enthalpy from temperature and total water content
 ! ************************************************************************************************************************
 subroutine t2enthalpy_T(&
                       ! input: data structures
                       diag_data,                         & ! intent(in):  model diagnostic variables for a local HRU
                       mpar_data,                         & ! intent(in):  parameter data structure
                       indx_data,                         & ! intent(in):  model indices
                       lookup_data,                       & ! intent(in):  lookup table data structure
                       ! input: state variables for the vegetation canopy
                       scalarCanairTempTrial,             & ! intent(in):  trial value of canopy air temperature (K)
                       scalarCanopyTempTrial,             & ! intent(in):  trial value of canopy temperature (K)
                       scalarCanopyWatTrial,              & ! intent(in):  trial value of canopy total water (kg m-2)
                       scalarCanopyIceTrial,              & ! intent(in):  trial value of canopy ice content (kg m-2)
                       ! input: variables for the snow-soil domain
                       mLayerTempTrial,                   & ! intent(in):  trial vector of layer temperature (K)
                       mLayerVolFracWatTrial,             & ! intent(in):  trial vector of volumetric total water content (-)
                       mLayerMatricHeadTrial,             & ! intent(in):  trial vector of total water matric potential (m)
                       mLayerVolFracIceTrial,             & ! intent(in)
                       ! output: enthalpy
                       scalarCanairEnthalpy,              & ! intent(out):  enthalpy of the canopy air space (J m-3)
                       scalarCanopyEnthalpy,              & ! intent(out):  enthalpy of the vegetation canopy (J m-3)
                       mLayerEnthalpy,                    & ! intent(out):  enthalpy of each snow+soil layer (J m-3)
                       ! output: error control
                       err,message)                         ! intent(out): error control
 ! -------------------------------------------------------------------------------------------------------------------------
 ! downwind routines
 USE soil_utils_module,only:crit_soilT     ! compute critical temperature below which ice exists
 USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water
 USE spline_int_module,only:splint         ! use for cubic spline interpolation
 implicit none
 ! delare dummy variables
 ! -------------------------------------------------------------------------------------------------------------------------
 ! input: data structures
 type(var_dlength),intent(in)  :: diag_data                 ! diagnostic variables for a local HRU
 type(var_dlength),intent(in)  :: mpar_data                 ! model parameters
 type(var_ilength),intent(in)  :: indx_data                 ! model indices
 type(zLookup),intent(in)      :: lookup_data               ! lookup tables
 ! input: state variables for the vegetation canopy
 real(rkind),intent(in)           :: scalarCanairTempTrial     ! trial value of canopy air temperature (K)
 real(rkind),intent(in)           :: scalarCanopyTempTrial     ! trial value of canopy temperature (K)
 real(rkind),intent(in)           :: scalarCanopyWatTrial      ! trial value of canopy total water (kg m-2)
 real(rkind),intent(in)           :: scalarCanopyIceTrial      ! trial value of canopy ice content (kg m-2)
 ! input: variables for the snow-soil domain
 real(rkind),intent(in)           :: mLayerTempTrial(:)        ! trial vector of layer temperature (K)
 real(rkind),intent(in)           :: mLayerVolFracWatTrial(:)  ! trial vector of volumetric total water content (-)
 real(rkind),intent(in)           :: mLayerMatricHeadTrial(:)  ! trial vector of total water matric potential (m)
 real(rkind),intent(in)           :: mLayerVolFracIceTrial(:)  ! trial vector of volumetric fraction of Ice (-)
 ! output: enthalpy
 real(rkind),intent(out)          :: scalarCanairEnthalpy      ! enthalpy of the canopy air space (J m-3)
 real(rkind),intent(out)          :: scalarCanopyEnthalpy      ! enthalpy of the vegetation canopy (J m-3)
 real(rkind),intent(out)          :: mLayerEnthalpy(:)         ! enthalpy of each snow+soil layer (J m-3)
 ! output: error control
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! -------------------------------------------------------------------------------------------------------------------------
 ! declare local variables
 character(len=128)            :: cmessage                  ! error message in downwind routine
 integer(i4b)                  :: iState                    ! index of model state variable
 integer(i4b)                  :: iLayer                    ! index of model layer
 integer(i4b)                  :: ixFullVector              ! index within full state vector
 integer(i4b)                  :: ixDomainType              ! name of a given model domain
 integer(i4b)                  :: ixControlIndex            ! index within a given model domain
 real(rkind)                      :: vGn_m                     ! van Genuchten "m" parameter (-)
 real(rkind)                      :: Tcrit                     ! temperature where all water is unfrozen (K)
 real(rkind)                      :: psiLiq                    ! matric head of liquid water (m)
 real(rkind)                      :: vFracWat                  ! volumetric fraction of total water, liquid+ice (-)
 real(rkind)                      :: vFracLiq                  ! volumetric fraction of liquid water (-)
 real(rkind)                      :: vFracIce                  ! volumetric fraction of ice (-)
 real(rkind)                      :: enthTemp                  ! enthalpy at the temperature of the control volume (J m-3)
 real(rkind)                      :: enthTcrit                 ! enthalpy at the critical temperature where all water is unfrozen (J m-3)
 real(rkind)                      :: enthPhase                 ! enthalpy associated with phase change (J m-3)
 real(rkind)                      :: enthWater                 ! enthalpy of total water (J m-3)
 real(rkind)                      :: enthSoil                  ! enthalpy of soil particles (J m-3)
 real(rkind)                      :: enthMix                   ! enthalpy of the mixed region, liquid+ice (J m-3)
 real(rkind)                      :: enthLiq                   ! enthalpy of the liquid region (J m-3)
 real(rkind)                      :: enthAir                   ! enthalpy of air (J m-3)
 real(rkind)                      :: diffT
 real(rkind)                      :: integral
 real(rkind)                      :: enthIce
 real(rkind)                      :: enthVeg
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 generalVars: associate(&
 ! number of model layers, and layer type
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):  [i4b]    total number of snow layers
 nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):  [i4b]    total number of soil layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in):  [i4b]    total number of snow and soil layers
 ! mapping between the full state vector and the state subset
 ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat         ,& ! intent(in):  [i4b(:)] list of indices in the state subset for each state in the full state vector
 ixMapSubset2Full        => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat         ,& ! intent(in):  [i4b(:)] [state subset] list of indices of the full state vector in the state subset
 ! type of domain, type of state variable, and index of control volume within domain
 ixDomainType_subset     => indx_data%var(iLookINDEX%ixDomainType_subset)%dat      ,& ! intent(in):  [i4b(:)] [state subset] id of domain for desired model state variables
 ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat          ,& ! intent(in):  [i4b(:)] index of the control volume for different domains (veg, snow, soil)
 ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,& ! intent(in):  [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
 ! snow parameters
 snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)          & ! intent(in):  [dp] scaling parameter for the snow freezing curve (K-1)
 ) ! end associate statement
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! initialize error control
 err=0; message="t2enthalpy_T/"

 ! loop through model state variables
 do iState=1,size(ixMapSubset2Full)

  ! -----
  ! - compute indices...
  ! --------------------

  ! get domain type, and index of the control volume within the domain
  ixFullVector   = ixMapSubset2Full(iState)       ! index within full state vector
  ixDomainType   = ixDomainType_subset(iState)    ! named variables defining the domain (iname_cas, iname_veg, etc.)
  ixControlIndex = ixControlVolume(ixFullVector)  ! index within a given domain

  ! check an energy state
  if(ixStateType(ixFullVector)==iname_nrgCanair .or. ixStateType(ixFullVector)==iname_nrgCanopy .or. ixStateType(ixFullVector)==iname_nrgLayer)then

   ! get the layer index
   select case(ixDomainType)
    case(iname_cas);     iLayer = integerMissing
    case(iname_veg);     iLayer = integerMissing
    case(iname_snow);    iLayer = ixControlIndex
    case(iname_soil);    iLayer = ixControlIndex + nSnow
    case(iname_aquifer); cycle ! aquifer: do nothing (no thermodynamics in the aquifer)
    case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil, iname_aquifer'; return
   end select

   ! identify domain
   select case(ixDomainType)

    ! -----
    ! - canopy air space...
    ! ---------------------
    case(iname_cas)
       scalarCanairEnthalpy = Cp_air * iden_air * (scalarCanairTempTrial - Tfreeze)
    ! -----
    ! - vegetation canopy...
    ! -----------------------
    case(iname_veg)
      ! association to necessary variables for vegetation
      vegVars: associate(&     
      canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1),         & ! intent(in): [dp]      canopy depth (m)
      specificHeatVeg         => mpar_data%var(iLookPARAM%specificHeatVeg)%dat(1),          & ! intent(in): specific heat of vegetation (J kg-1 K-1)
      maxMassVegetation       => mpar_data%var(iLookPARAM%maxMassVegetation)%dat(1),        & ! intent(in): maximum mass of vegetation (kg m-2)
      snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)   & ! intent(in):  [dp] scaling parameter for the snow freezing curve (K-1)
      ) 

      diffT = scalarCanopyTempTrial - Tfreeze
      enthVeg = specificHeatVeg * maxMassVegetation * diffT / canopyDepth
      if(diffT>=0._rkind)then
         enthLiq = Cp_water * scalarCanopyWatTrial * diffT / canopyDepth
         enthIce = 0._rkind
         enthPhase = 0._rkind
      else
         integral = (1._rkind/snowfrz_scale) * atan(snowfrz_scale * diffT)
         enthLiq = Cp_water * scalarCanopyWatTrial * integral / canopyDepth
         enthIce = Cp_ice * scalarCanopyWatTrial * ( diffT - integral ) / canopyDepth
         enthPhase = LH_fus * scalarCanopyIceTrial / canopyDepth
      end if
      
       scalarCanopyEnthalpy =  enthVeg + enthLiq + enthIce           

      ! end association 
      end associate vegVars

    ! -----
    ! - snow...
    ! ---------
    case(iname_snow)
      
      ! association to necessary variables for snow
      snowVars: associate(&
      snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)   & ! intent(in):  [dp] scaling parameter for the snow freezing curve (K-1)
      ) 

      diffT = mLayerTempTrial(iLayer) - Tfreeze
      if(diffT>=0._rkind)then
         enthLiq = iden_water * Cp_water * mLayerVolFracWatTrial(iLayer) * diffT
         enthIce = 0._rkind
         enthAir = iden_air * Cp_air * ( 1._rkind - mLayerVolFracWatTrial(iLayer) ) * diffT
         enthPhase = 0._rkind
      else
         integral = (1._rkind/snowfrz_scale) * atan(snowfrz_scale * diffT)
         enthLiq = iden_water * Cp_water * mLayerVolFracWatTrial(iLayer) * integral
         enthIce = iden_water * Cp_ice * mLayerVolFracWatTrial(iLayer) * ( diffT - integral )
         enthAir = iden_air * Cp_air * ( diffT - mLayerVolFracWatTrial(iLayer) * ( (iden_water/iden_ice)*(diffT-integral) + integral ) ) 
         enthPhase = iden_ice * LH_fus * mLayerVolFracIceTrial(iLayer)
      end if
      
       mLayerEnthalpy(iLayer) =  enthLiq + enthIce + enthAir           

      ! end association 
      end associate snowVars

    ! -----
    ! - soil...
    ! ---------
    case(iname_soil)

     ! make association to variables in the data structures...
     soilVars: associate(&

     ! associate model parameters
     soil_dens_intr => mpar_data%var(iLookPARAM%soil_dens_intr)%dat(ixControlIndex)      , & ! intrinsic soil density             (kg m-3)
     theta_sat      => mpar_data%var(iLookPARAM%theta_sat)%dat(ixControlIndex)           , & ! soil porosity                      (-)
     theta_res      => mpar_data%var(iLookPARAM%theta_res)%dat(ixControlIndex)           , & ! volumetric residual water content  (-)
     vGn_alpha      => mpar_data%var(iLookPARAM%vGn_alpha)%dat(ixControlIndex)           , & ! van Genuchten "alpha" parameter    (m-1)
     vGn_n          => mpar_data%var(iLookPARAM%vGn_n)%dat(ixControlIndex)               , & ! van Genuchten "n" parameter        (-)

     ! associate values in the lookup table
     Tk            => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%temperature)%lookup  , & ! temperature (K)
     Ey            => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%enthalpy)%lookup     , & ! enthalpy (J m-3)
     E2            => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%deriv2)%lookup         & ! second derivative of the interpolating function

     ) ! end associate statement

     ! diagnostic variables
     vGn_m    = 1._rkind - 1._rkind/vGn_n
     Tcrit    = crit_soilT( mLayerMatricHeadTrial(ixControlIndex) )
     vFracWat = volFracLiq(mLayerMatricHeadTrial(ixControlIndex),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
     

     ! *** compute enthalpy of water for unfrozen conditions
     if(mlayerTempTrial(iLayer) > Tcrit)then
      enthWater = iden_water*Cp_water*vFracWat*(mlayerTempTrial(iLayer) - Tfreeze) ! valid for temperatures below freezing also
      enthPhase = 0._rkind

     ! *** compute enthalpy of water for frozen conditions
     else
      ! calculate enthalpy at the temperature (cubic spline interpolation)
      call splint(Tk,Ey,E2,mlayerTempTrial(iLayer),enthTemp,err,cmessage)
      if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

      ! calculate enthalpy at the critical temperature (cubic spline interpolation)
      call splint(Tk,Ey,E2,Tcrit,enthTcrit,err,cmessage)
      if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

      ! calculate the enthalpy of water
      enthMix   = enthTemp - enthTcrit ! enthalpy of the liquid+ice mix
      enthLiq   = iden_water*Cp_water*vFracWat*(Tcrit - Tfreeze)
      enthWater = enthMix + enthLiq

      ! *** compute the enthalpy associated with phase change
      psiLiq    = (mLayerTempTrial(iLayer) - Tfreeze)*LH_fus/(gravity*Tfreeze)
      vFracLiq  = volFracLiq(psiLiq,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
      vFracIce  = vFracWat - vFracLiq
      enthPhase = iden_water*LH_fus*vFracIce

     endif ! (if frozen conditions)

     ! *** compute the enthalpy of soil
     enthSoil  = soil_dens_intr*Cp_soil*(1._rkind - theta_sat)*(mlayerTempTrial(iLayer) - Tfreeze)

     ! *** compute the enthalpy of air
     enthAir   = iden_air*Cp_air*(1._rkind - theta_sat - vFracWat)*(mlayerTempTrial(iLayer) - Tfreeze)

     ! *** compute the total enthalpy (J m-3)
     mLayerEnthalpy(iLayer) = enthSoil + enthWater + enthAir

!     print*, iLayer, enthSoil, enthWater, enthAir, enthPhase

     end associate soilVars

    ! -----
    ! - checks...
    ! -----------
    case(iname_aquifer); cycle ! aquifer: do nothing
    case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil, iname_aquifer'; return
   end select

  end if  ! if an energy layer
 end do  ! looping through state variables

 end associate generalVars

 end subroutine t2enthalpy_T

end module t2enthalpy_module
