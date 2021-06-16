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

module snwDensify_module

! data types
USE nrtype

! model constants
USE multiconst,only:&
                    Tfreeze,    &     ! freezing point of pure water (K)
                    iden_ice,   &     ! intrinsic density of ice (kg m-3)
                    iden_air,   &     ! intrinsic density of air (kg m-3)
                    iden_water        ! intrinsic density of liquid water (kg m-3)
! privacy
implicit none
private
public::snwDensify
contains

 ! ************************************************************************************************
 ! public subroutine snwDensify: compute change in snow density over the time step
 ! ************************************************************************************************
 subroutine snwDensify(&

                       ! intent(in): variables
                       dt,                             & ! intent(in): time step (s)
                       nSnow,                          & ! intent(in): number of snow layers
                       mLayerTemp,                     & ! intent(in): temperature of each layer (K)
                       mLayerMeltFreeze,               & ! intent(in): volumnetric melt in each layer (kg m-3)

                       ! intent(in): parameters
                       densScalGrowth,                 & ! intent(in): density scaling factor for grain growth (kg-1 m3)
                       tempScalGrowth,                 & ! intent(in): temperature scaling factor for grain growth (K-1)
                       grainGrowthRate,                & ! intent(in): rate of grain growth (s-1)
                       densScalOvrbdn,                 & ! intent(in): density scaling factor for overburden pressure (kg-1 m3)
                       tempScalOvrbdn,                 & ! intent(in): temperature scaling factor for overburden pressure (K-1)
                       baseViscosity,                      & ! intent(in): viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)

                       ! intent(inout): state variables
                       mLayerDepth,                    & ! intent(inout): depth of each layer (m)
                       mLayerVolFracLiqNew,            & ! intent(inout):  volumetric fraction of liquid water after itertations (-)
                       mLayerVolFracIceNew,            & ! intent(inout):  volumetric fraction of ice after itertations (-)

                       ! output: error control
                       err,message)                      ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! compute change in snow density over the time step
 implicit none
 ! intent(in): variables
 real(rkind),intent(in)                 :: dt                       ! time step (seconds)
 integer(i4b),intent(in)             :: nSnow                    ! number of snow layers
 real(rkind),intent(in)                 :: mLayerTemp(:)            ! temperature of each snow layer after iterations (K)
 real(rkind),intent(in)                 :: mLayerMeltFreeze(:)      ! volumetric melt in each layer (kg m-3)
 ! intent(in): parameters
 real(rkind),intent(in)                 :: densScalGrowth           ! density scaling factor for grain growth (kg-1 m3)
 real(rkind),intent(in)                 :: tempScalGrowth           ! temperature scaling factor for grain growth (K-1)
 real(rkind),intent(in)                 :: grainGrowthRate          ! rate of grain growth (s-1)
 real(rkind),intent(in)                 :: densScalOvrbdn           ! density scaling factor for overburden pressure (kg-1 m3)
 real(rkind),intent(in)                 :: tempScalOvrbdn           ! temperature scaling factor for overburden pressure (K-1)
 real(rkind),intent(in)                 :: baseViscosity            ! viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)
 ! intent(inout): state variables
 real(rkind),intent(inout)              :: mLayerDepth(:)           ! depth of each layer (m)
 real(rkind),intent(inout)              :: mLayerVolFracLiqNew(:)   ! volumetric fraction of liquid water in each snow layer after iterations (-)
 real(rkind),intent(inout)              :: mLayerVolFracIceNew(:)   ! volumetric fraction of ice in each snow layer after iterations (-)
 ! intent(out): error control
 integer(i4b),intent(out)            :: err                      ! error code
 character(*),intent(out)            :: message                  ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! define local variables
 integer(i4b)                        :: iSnow                    ! index of snow layers
 real(rkind)                            :: chi1,chi2,chi3,chi4,chi5 ! multipliers in the densification algorithm (-)
 real(rkind)                            :: halfWeight               ! half of the weight of the current snow layer (kg m-2)
 real(rkind)                            :: weightSnow               ! total weight of snow above the current snow layer (kg m-2)
 real(rkind)                            :: CR_grainGrowth           ! compaction rate for grain growth (s-1)
 real(rkind)                            :: CR_ovrvdnPress           ! compaction rate associated with over-burden pressure (s-1)
 real(rkind)                            :: CR_metamorph             ! compaction rate for metamorphism (s-1)
 real(rkind)                            :: massIceOld               ! mass of ice in the snow layer (kg m-2)
 real(rkind)                            :: massLiqOld               ! mass of liquid water in the snow layer (kg m-2)
 real(rkind)                            :: scalarDepthNew           ! updated layer depth (m)
 real(rkind)                            :: scalarDepthMin           ! minimum layer depth (m)
 real(rkind)                            :: volFracIceLoss           ! volumetric fraction of ice lost due to melt and sublimation (-)
 real(rkind), dimension(nSnow)          :: mLayerVolFracAirNew      ! volumetric fraction of air in each layer after compaction (-)
 real(rkind),parameter                  :: snwden_min=100._rkind       ! minimum snow density for reducing metamorphism rate (kg m-3)
 real(rkind),parameter                  :: snwDensityMax=550._rkind    ! maximum snow density for collapse under melt (kg m-3)
 real(rkind),parameter                  :: wetSnowThresh=0.01_rkind    ! threshold to discriminate between "wet" and "dry" snow
 real(rkind),parameter                  :: minLayerDensity=40._rkind   ! minimum snow density allowed for any layer (kg m-3)
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="snwDensify/"

 ! NOTE: still need to process the case of "snow without a layer"
 if(nSnow==0)return

 ! initialize the weight of snow above each layer (kg m-2)
 weightSnow = 0._rkind

 ! loop through snow layers
 do iSnow=1,nSnow

  ! print starting density
  !write(*,'(a,1x,i4,1x,f9.3)') 'b4 compact: iSnow, density = ', iSnow, mLayerVolFracIceNew(iSnow)*iden_ice

  ! save mass of liquid water and ice (mass does not change)
  massIceOld = iden_ice*mLayerVolFracIceNew(iSnow)*mLayerDepth(iSnow)   ! (kg m-2)
  massLiqOld = iden_water*mLayerVolFracLiqNew(iSnow)*mLayerDepth(iSnow) ! (kg m-2)

  ! *** compute the compaction associated with grain growth (s-1)
  ! compute the base rate of grain growth (-)
  if(mLayerVolFracIceNew(iSnow)*iden_ice <snwden_min) chi1=1._rkind
  if(mLayerVolFracIceNew(iSnow)*iden_ice>=snwden_min) chi1=exp(-densScalGrowth*(mLayerVolFracIceNew(iSnow)*iden_ice - snwden_min))
  ! compute the reduction of grain growth under colder snow temperatures (-)
  chi2 = exp(-tempScalGrowth*(Tfreeze - mLayerTemp(iSnow)))
  ! compute the acceleration of grain growth in the presence of liquid water (-)
  if(mLayerVolFracLiqNew(iSnow) > wetSnowThresh)then; chi3=2._rkind  ! snow is "wet"
  else; chi3=1._rkind; end if                                         ! snow is "dry"
  ! compute the compaction associated with grain growth (s-1)
  CR_grainGrowth = grainGrowthRate*chi1*chi2*chi3

  ! **** compute the compaction associated with over-burden pressure (s-1)
  ! compute the weight imposed on the current layer (kg m-2)
  halfWeight = (massIceOld + massLiqOld)/2._rkind  ! there is some over-burden pressure from the layer itself
  weightSnow = weightSnow + halfweight          ! add half of the weight from the current layer
  ! compute the increase in compaction under colder snow temperatures (-)
  chi4 = exp(-tempScalOvrbdn*(Tfreeze - mLayerTemp(iSnow)))
  ! compute the increase in compaction under low density snow (-)
  chi5 = exp(-densScalOvrbdn*mLayerVolFracIceNew(iSnow)*iden_ice)
  ! compute the compaction associated with over-burden pressure (s-1)
  CR_ovrvdnPress = (weightSnow/baseViscosity)*chi4*chi5
  ! update the snow weight with the halfWeight not yet used
  weightSnow = weightSnow + halfweight          ! add half of the weight from the current layer

  ! *** compute the compaction rate associated with snow melt (s-1)
  ! NOTE: loss of ice due to snowmelt is implicit, so can be updated directly
  if(iden_ice*mLayerVolFracIceNew(iSnow) < snwDensityMax)then ! only collapse layers if below a critical density
   ! (compute volumetric losses of ice due to melt and sublimation)
   volFracIceLoss = max(0._rkind,mLayerMeltFreeze(iSnow)/iden_ice)  ! volumetric fraction of ice lost due to melt (-)
   ! (adjust snow depth to account for cavitation)
   scalarDepthNew = mLayerDepth(iSnow) * mLayerVolFracIceNew(iSnow)/(mLayerVolFracIceNew(iSnow) + volFracIceLoss)
   !print*, 'volFracIceLoss = ', volFracIceLoss
  else
   scalarDepthNew = mLayerDepth(iSnow)
  end if
  ! compute the total compaction rate associated with metamorphism
  CR_metamorph = CR_grainGrowth + CR_ovrvdnPress
  ! update depth due to metamorphism (implicit solution)
  ! Ensure that the new depth is in line with the maximum amount of compaction that
  ! can occur given the masses of ice and liquid in the layer
  scalarDepthNew = scalarDepthNew/(1._rkind + CR_metamorph*dt)
  scalarDepthMin = (massIceOld / iden_ice) + (massLiqOld / iden_water)
  mLayerDepth(iSnow) = max(scalarDepthMin, scalarDepthNew)

  ! check that depth is reasonable
  if(mLayerDepth(iSnow) < 0._rkind)then
   write(*,'(a,1x,i4,1x,10(f12.5,1x))') 'iSnow, dt, density,massIceOld, massLiqOld = ', iSnow, dt, mLayerVolFracIceNew(iSnow)*iden_ice, massIceOld, massLiqOld
   write(*,'(a,1x,i4,1x,10(f12.5,1x))') 'iSnow, mLayerDepth(iSnow), scalarDepthNew, mLayerVolFracIceNew(iSnow), mLayerMeltFreeze(iSnow), CR_grainGrowth*dt, CR_ovrvdnPress*dt = ', &
                                         iSnow, mLayerDepth(iSnow), scalarDepthNew, mLayerVolFracIceNew(iSnow), mLayerMeltFreeze(iSnow), CR_grainGrowth*dt, CR_ovrvdnPress*dt
  endif

  ! update volumetric ice and liquid water content
  mLayerVolFracIceNew(iSnow) = massIceOld/(mLayerDepth(iSnow)*iden_ice)
  mLayerVolFracLiqNew(iSnow) = massLiqOld/(mLayerDepth(iSnow)*iden_water)
  mLayerVolFracAirNew(iSnow) = 1.0_rkind - mLayerVolFracIceNew(iSnow) - mLayerVolFracLiqNew(iSnow)
  !write(*,'(a,1x,i4,1x,f9.3)') 'after compact: iSnow, density = ', iSnow, mLayerVolFracIceNew(iSnow)*iden_ice
  !if(mLayerMeltFreeze(iSnow) > 20._rkind) pause 'meaningful melt'

 end do  ! looping through snow layers

 ! check depth
 if(any(mLayerDepth(1:nSnow) < 0._rkind))then
  do iSnow=1,nSnow
   write(*,'(a,1x,i4,1x,4(f12.5,1x))') 'iSnow, mLayerDepth(iSnow)', iSnow, mLayerDepth(iSnow)
  end do
  message=trim(message)//'unreasonable value for snow depth'
  err=20; return
 end if

 ! check for low/high snow density
 if(any(mLayerVolFracIceNew(1:nSnow)*iden_ice + mLayerVolFracLiqNew(1:nSnow)*iden_water + mLayerVolFracAirNew(1:nSnow)*iden_air < minLayerDensity) .or. &
    any(mLayerVolFracIceNew(1:nSnow) + mLayerVolFracLiqNew(1:nSnow) + mLayerVolFracAirNew(1:nSnow) > 1._rkind))then
  do iSnow=1,nSnow
   write(*,*) 'iSnow, volFracIce, density = ', iSnow, mLayerVolFracIceNew(iSnow),  mLayerVolFracIceNew(iSnow)*iden_ice
  end do
  message=trim(message)//'unreasonable value for snow density'
  err=20; return
 end if

 end subroutine snwDensify


end module snwDensify_module
