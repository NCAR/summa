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

module snowLiqFlx_module

! access modules
USE nrtype                                    ! numerical recipes data types
USE multiconst,only:iden_ice,iden_water       ! intrinsic density of ice and water (kg m-3)

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! named variables
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements

! data types
USE data_types,only:var_d           ! x%var(:)       (dp)
USE data_types,only:var_dlength     ! x%var(:)%dat   (dp)
USE data_types,only:var_ilength     ! x%var(:)%dat   (i4b)

! privacy
implicit none
private
public::snowLiqFlx
contains


 ! ************************************************************************************************
 ! public subroutine snowLiqFlx: compute liquid water flux through the snowpack
 ! ************************************************************************************************
 subroutine snowLiqFlx(&
                       ! input: model control
                       nSnow,                   & ! intent(in):    number of snow layers
                       firstFluxCall,           & ! intent(in):    the first flux call
                       scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                       ! input: forcing for the snow domain
                       scalarThroughfallRain,   & ! intent(in):    rain that reaches the snow surface without ever touching vegetation (kg m-2 s-1)
                       scalarCanopyLiqDrainage, & ! intent(in):    liquid drainage from the vegetation canopy (kg m-2 s-1)
                       ! input: model state vector
                       mLayerVolFracLiqTrial,   & ! intent(in):    trial value of volumetric fraction of liquid water at the current iteration (-)
                       ! input-output: data structures
                       indx_data,               & ! intent(in):    model indices
                       mpar_data,               & ! intent(in):    model parameters
                       prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                       diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                       ! output: fluxes and derivatives
                       iLayerLiqFluxSnow,       & ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
                       iLayerLiqFluxSnowDeriv,  & ! intent(inout): derivative in vertical liquid water flux at layer interfaces (m s-1)
                       ! output: error control
                       err,message)               ! intent(out):   error control
 implicit none
 ! input: model control
 integer(i4b),intent(in)         :: nSnow                      ! number of snow layers
 logical(lgt),intent(in)         :: firstFluxCall              ! the first flux call
 logical(lgt),intent(in)         :: scalarSolution             ! flag to denote if implementing the scalar solution
 ! input: forcing for the snow domain
 real(rkind),intent(in)             :: scalarThroughfallRain      ! computed throughfall rate (kg m-2 s-1)
 real(rkind),intent(in)             :: scalarCanopyLiqDrainage    ! computed drainage of liquid water (kg m-2 s-1)
 ! input: model state vector
 real(rkind),intent(in)             :: mLayerVolFracLiqTrial(:)   ! trial value of volumetric fraction of liquid water at the current iteration (-)
 ! input-output: data structures
 type(var_ilength),intent(in)    :: indx_data                  ! model indices
 type(var_dlength),intent(in)    :: mpar_data                  ! model parameters
 type(var_dlength),intent(in)    :: prog_data                  ! prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data                  ! diagnostic variables for a local HRU
 ! output: fluxes and derivatives
 real(rkind),intent(inout)          :: iLayerLiqFluxSnow(0:)      ! vertical liquid water flux at layer interfaces (m s-1)
 real(rkind),intent(inout)          :: iLayerLiqFluxSnowDeriv(0:) ! derivative in vertical liquid water flux at layer interfaces (m s-1)
 ! output: error control
 integer(i4b),intent(out)        :: err                        ! error code
 character(*),intent(out)        :: message                    ! error message
 ! ------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: i                          ! search index for scalar solution
 integer(i4b)                    :: iLayer                     ! layer index
 integer(i4b)                    :: ixTop                      ! top layer in subroutine call
 integer(i4b)                    :: ixBot                      ! bottom layer in subroutine call
 real(rkind)                        :: multResid                  ! multiplier for the residual water content (-)
 real(rkind),parameter              :: residThrs=550._rkind          ! ice density threshold to reduce residual liquid water content (kg m-3)
 real(rkind),parameter              :: residScal=10._rkind           ! scaling factor for residual liquid water content reduction factor (kg m-3)
 real(rkind),parameter              :: maxVolIceContent=0.7_rkind    ! maximum volumetric ice content to store water (-)
 real(rkind)                        :: availCap                   ! available storage capacity [0,1] (-)
 real(rkind)                        :: relSaturn                  ! relative saturation [0,1] (-)
 ! ------------------------------------------------------------------------------------------------------------------------------------------
 ! make association of local variables with information in the data structures
 associate(&
  ! input: layer indices
  ixLayerState     => indx_data%var(iLookINDEX%ixLayerState)%dat,             & ! intent(in): list of indices for all model layers
  ixSnowOnlyHyd    => indx_data%var(iLookINDEX%ixSnowOnlyHyd)%dat,            & ! intent(in): index in the state subset for hydrology state variables in the snow domain
  ! input: snow properties and parameters
  mLayerVolFracIce => prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow), & ! intent(in): volumetric ice content at the start of the time step (-)
  Fcapil           => mpar_data%var(iLookPARAM%Fcapil)%dat(1),                & ! intent(in): capillary retention as a fraction of the total pore volume (-)
  k_snow           => mpar_data%var(iLookPARAM%k_snow)%dat(1),                & ! intent(in): hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
  mw_exp           => mpar_data%var(iLookPARAM%mw_exp)%dat(1),                & ! intent(in): exponent for meltwater flow (-)
  ! input/output: diagnostic variables -- only computed for the first iteration
  mLayerPoreSpace  => diag_data%var(iLookDIAG%mLayerPoreSpace)%dat,           & ! intent(inout): pore space in each snow layer (-)
  mLayerThetaResid => diag_data%var(iLookDIAG%mLayerThetaResid)%dat           & ! intent(inout): esidual volumetric liquid water content in each snow layer (-)
 ) ! association of local variables with information in the data structures
 ! ------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='snowLiqFlx/'

 ! check that the input vectors match nSnow
 if(size(mLayerVolFracLiqTrial)/=nSnow .or. size(mLayerVolFracIce)/=nSnow .or. &
    size(iLayerLiqFluxSnow)/=nSnow+1 .or. size(iLayerLiqFluxSnowDeriv)/=nSnow+1) then
  err=20; message=trim(message)//'size mismatch of input/output vectors'; return
 end if

 ! check the meltwater exponent is >=1
 if(mw_exp<1._rkind)then; err=20; message=trim(message)//'meltwater exponent < 1'; return; end if

 ! get the indices for the snow+soil layers
 ixTop = integerMissing
 if(scalarSolution)then
  ! WARNING: Previously this was implemented as:
  !    ixLayerDesired = pack(ixLayerState, ixSnowOnlyHyd/=integerMissing)
  !    ixTop = ixLayerDesired(1)
  !    ixBot = ixLayerDesired(1)
  ! This implementation can result in a segfault when using JRDN layering.
  ! The segfault occurs when trying to access `mw_exp` in:
  !    iLayerLiqFluxSnow(iLayer)      = k_snow*relSaturn**mw_exp
  ! Debugging found that the `pack` statement caused `mw_exp` to no longer be accessible.
  ! We have not been able to determine the underlying reason for this segfault.
  do i=1,size(ixSnowOnlyHyd)
    if(ixSnowOnlyHyd(i) /= integerMissing)then
      ixTop=ixLayerState(i)
      ixBot=ixTop
      exit  ! break out of loop once found
    endif
  end do
  if(ixTop == integerMissing)then
    err=20; message=trim(message)//'Unable to identify snow layer for scalar solution!'; return
  end if
 else
  ixTop = 1
  ixBot = nSnow
 endif

 ! define the liquid flux at the upper boundary (m s-1)
 iLayerLiqFluxSnow(0)      = (scalarThroughfallRain + scalarCanopyLiqDrainage)/iden_water
 iLayerLiqFluxSnowDeriv(0) = 0._rkind

 ! compute properties fixed over the time step
 if(firstFluxCall)then
  ! loop through snow layers
  do iLayer=1,nSnow
   ! compute the reduction in liquid water holding capacity at high snow density (-)
   multResid = 1._rkind / ( 1._rkind + exp( (mLayerVolFracIce(iLayer)*iden_ice - residThrs) / residScal) )
   ! compute the pore space (-)
   mLayerPoreSpace(iLayer)  = 1._rkind - mLayerVolFracIce(iLayer)
   ! compute the residual volumetric liquid water content (-)
   mLayerThetaResid(iLayer) = Fcapil*mLayerPoreSpace(iLayer) * multResid
  end do  ! (looping through snow layers)
 end if  ! (if the first flux call)

 ! compute fluxes
 do iLayer=ixTop,ixBot  ! (loop through snow layers)
  ! check that flow occurs
  if(mLayerVolFracLiqTrial(iLayer) > mLayerThetaResid(iLayer))then
   ! compute the relative saturation (-)
   availCap  = mLayerPoreSpace(iLayer) - mLayerThetaResid(iLayer)                 ! available capacity
   relSaturn = (mLayerVolFracLiqTrial(iLayer) - mLayerThetaResid(iLayer)) / availCap    ! relative saturation
   iLayerLiqFluxSnow(iLayer)      = k_snow*relSaturn**mw_exp
   iLayerLiqFluxSnowDeriv(iLayer) = ( (k_snow*mw_exp)/availCap ) * relSaturn**(mw_exp - 1._rkind)
   if(mLayerVolFracIce(iLayer) > maxVolIceContent)then ! NOTE: use start-of-step ice content, to avoid convergence problems
     ! ** allow liquid water to pass through under very high ice density
     iLayerLiqFluxSnow(iLayer) = iLayerLiqFluxSnow(iLayer) + iLayerLiqFluxSnow(iLayer-1) !NOTE: derivative may need to be updated in future.
   end if
  else  ! flow does not occur
   iLayerLiqFluxSnow(iLayer)      = 0._rkind
   iLayerLiqFluxSnowDeriv(iLayer) = 0._rkind
  endif  ! storage above residual content
 end do  ! loop through snow layers

 ! end association of local variables with information in the data structures
 end associate

 end subroutine snowLiqFlx


end module snowLiqFlx_module
