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

module soil_utils_fida_module

! data types
USE nrtype

USE multiconst,only: gravity, & ! acceleration of gravity       (m s-2)
                     Tfreeze, & ! temperature at freezing    (K)
                     LH_fus,  & ! latent heat of fusion      (J kg-1, or m2 s-2)
                     R_wv       ! gas constant for water vapor  (J kg-1 K-1; [J = Pa m3])
USE soil_utils_module,only:matricHead
USE soil_utils_module,only:dPsi_dTheta
USE soil_utils_module,only:volFracLiq
USE soil_utils_module,only:dTheta_dPsi

! privacy
implicit none
private

! routines to make public

public::liquidHeadFida
public::d2Theta_dPsi2
public::d2Theta_dTk2

! constant parameters
real(dp),parameter     :: verySmall=epsilon(1.0_dp) ! a very small number (used to avoid divide by zero)
contains


 ! ******************************************************************************************************************************
 ! public subroutine: compute the liquid water matric potential (and the derivatives w.r.t. total matric potential and temperature)
 ! ******************************************************************************************************************************
 subroutine liquidHeadFida(&
                       ! input
                       matricHeadTotal                          ,& ! intent(in)    : total water matric potential (m)
                       matricHeadTotalPrime                     ,& ! intent(in)
                       volFracLiq                               ,& ! intent(in)    : volumetric fraction of liquid water (-)
                       volFracIce                               ,& ! intent(in)    : volumetric fraction of ice (-)
                       vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m,& ! intent(in)    : soil parameters
                       dVolTot_dPsi0                            ,& ! intent(in)    : derivative in the soil water characteristic (m-1)
                       dTheta_dT                                ,& ! intent(in)    : derivative in volumetric total water w.r.t. temperature (K-1)
                       tempPrime                                ,& ! intent(in)
                       volFracLiqPrime                          ,& ! intent(in)
                       volFracIcePrime                          ,& ! intent(in)
                       ! output
                       matricHeadLiq                            ,& ! intent(out)   : liquid water matric potential (m)
                       matricHeadLiqPrime                       ,& ! intent(out)
                       dPsiLiq_dPsi0                            ,& ! intent(out)   : derivative in the liquid water matric potential w.r.t. the total water matric potential (-)
                       dPsiLiq_dTemp                            ,& ! intent(out)   : derivative in the liquid water matric potential w.r.t. temperature (m K-1)
                       err,message)                                ! intent(out)   : error control
 ! computes the liquid water matric potential (and the derivatives w.r.t. total matric potential and temperature)
 implicit none
 ! input
 real(dp),intent(in)            :: matricHeadTotal                           ! total water matric potential (m)
  real(dp),intent(in)           :: matricHeadTotalPrime
 real(dp),intent(in)            :: volFracLiq                                ! volumetric fraction of liquid water (-)
 real(dp),intent(in)            :: volFracIce                                ! volumetric fraction of ice (-)
 real(dp),intent(in)            :: vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m ! soil parameters
 real(dp),intent(in)  ,optional :: dVolTot_dPsi0                             ! derivative in the soil water characteristic (m-1)
 real(dp),intent(in)  ,optional :: dTheta_dT                                 ! derivative in volumetric total water w.r.t. temperature (K-1)
 real(dp),intent(in)            :: TempPrime
 real(dp),intent(in)            :: volFracLiqPrime
 real(dp),intent(in)            :: volFracIcePrime
 ! output
 real(dp),intent(out)           :: matricHeadLiq                             ! liquid water matric potential (m)
  real(dp),intent(out)          :: matricHeadLiqPrime 
 real(dp),intent(out) ,optional :: dPsiLiq_dPsi0                             ! derivative in the liquid water matric potential w.r.t. the total water matric potential (-)
 real(dp),intent(out) ,optional :: dPsiLiq_dTemp                             ! derivative in the liquid water matric potential w.r.t. temperature (m K-1)
 ! output: error control
 integer(i4b),intent(out)       :: err                                       ! error code
 character(*),intent(out)       :: message                                   ! error message
 ! local
 real(dp)                       :: xNum,xDen                                 ! temporary variables (numeratir, denominator)
 real(dp)                       :: effSat                                    ! effective saturation (-)
 real(dp)                       :: dPsiLiq_dEffSat                           ! derivative in liquid water matric potential w.r.t. effective saturation (m)
 real(dp)                       :: dEffSat_dTemp                             ! derivative in effective saturation w.r.t. temperature (K-1)
 real(dp)                       :: dEffSat_dFracLiq
 real(dp)                       :: effSatPrime
 ! ------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='liquidHeadFida/'

 ! ** partially frozen soil
 if(volFracIce > verySmall .and. matricHeadTotal < 0._dp)then  ! check that ice exists and that the soil is unsaturated
 

  ! -----
  ! - compute liquid water matric potential...
  ! ------------------------------------------

  ! - compute effective saturation
  ! NOTE: include ice content as part of the solid porosity - major effect of ice is to reduce the pore size; ensure that effSat=1 at saturation
  ! (from Zhao et al., J. Hydrol., 1997: Numerical analysis of simultaneous heat and mass transfer...)
  xNum   = volFracLiq - theta_res
  xDen   = theta_sat - volFracIce - theta_res
  effSat = xNum/xDen          ! effective saturation

  ! - matric head associated with liquid water
  matricHeadLiq = matricHead(effSat,vGn_alpha,0._dp,1._dp,vGn_n,vGn_m)  ! argument is effective saturation, so theta_res=0 and theta_sat=1
  if (effSat < 1._dp .and. effSat > 0._dp)then
    effSatPrime = (volFracLiqPrime * xDen + volFracIcePrime * xNum) / xDen**2._dp
    matricHeadLiqPrime = -( 1._dp/(vGn_alpha*vGn_n*vGn_m) ) * effSat**(-1._dp-1._dp/vGn_m) * ( effSat**(-1._dp/vGn_m) - 1._dp )**(-1._dp+1._dp/vGn_n) * effSatPrime
  else
   matricHeadLiqPrime = 0._dp
  endif
  

  ! compute derivative in liquid water matric potential w.r.t. effective saturation (m)
  if(present(dPsiLiq_dPsi0).or.present(dPsiLiq_dTemp))then
   dPsiLiq_dEffSat = dPsi_dTheta(effSat,vGn_alpha,0._dp,1._dp,vGn_n,vGn_m)
  endif

  ! -----
  ! - compute derivative in the liquid water matric potential w.r.t. the total water matric potential...
  ! ----------------------------------------------------------------------------------------------------

  ! check if the derivative is desired
  if(present(dPsiLiq_dTemp))then
   ! (check required input derivative is present)
   if(.not.present(dVolTot_dPsi0))then
    message=trim(message)//'dVolTot_dPsi0 argument is missing'
    err=20; return
   endif

   ! (compute derivative in the liquid water matric potential w.r.t. the total water matric potential)
   dPsiLiq_dPsi0 = dVolTot_dPsi0*dPsiLiq_dEffSat*xNum/(xDen**2._dp)
 !  matricHeadLiqPrime = dPsiLiq_dPsi0 * matricHeadTotalPrime

  endif  ! if dPsiLiq_dTemp is desired

  ! -----
  ! - compute the derivative in the liquid water matric potential w.r.t. temperature...
  ! -----------------------------------------------------------------------------------

  ! check if the derivative is desired
  if(present(dPsiLiq_dTemp))then

   ! (check required input derivative is present)
   if(.not.present(dTheta_dT))then
    message=trim(message)//'dTheta_dT argument is missing'
    err=20; return
   endif
   ! (compute the derivative in the liquid water matric potential w.r.t. temperature)
   dEffSat_dTemp = -dTheta_dT*xNum/(xDen**2._dp) + dTheta_dT/xDen
   dPsiLiq_dTemp = dPsiLiq_dEffSat*dEffSat_dTemp
  ! matricHeadLiqPrime = dPsiLiq_dTemp * tempPrime

   

  endif  ! if dPsiLiq_dTemp is desired

 ! ** unfrozen soil
 else   ! (no ice)
  matricHeadLiq = matricHeadTotal
  matricHeadLiqPrime = matricHeadTotalPrime
  if(present(dPsiLiq_dTemp)) dPsiLiq_dPsi0 = 1._dp  ! derivative=1 because values are identical
  if(present(dPsiLiq_dTemp)) dPsiLiq_dTemp = 0._dp  ! derivative=0 because no impact of temperature for unfrozen conditions
 end if  ! (if ice exists)

 end subroutine liquidHeadFida
 
 
  ! ******************************************************************************************************************************
 ! public function dTheta_dPsi: compute the derivative of the soil water characteristic (m-1)
 ! ******************************************************************************************************************************
 function d2Theta_dPsi2(psi,alpha,theta_res,theta_sat,n,m)
 implicit none
 real(dp),intent(in) :: psi         ! soil water suction (m)
 real(dp),intent(in) :: alpha       ! scaling parameter (m-1)
 real(dp),intent(in) :: theta_res   ! residual volumetric water content (-)
 real(dp),intent(in) :: theta_sat   ! porosity (-)
 real(dp),intent(in) :: n           ! vGn "n" parameter (-)
 real(dp),intent(in) :: m           ! vGn "m" parameter (-)
 real(dp)            :: d2Theta_dPsi2 ! derivative of the soil water characteristic (m-1)
 real(dp)            :: mult_fcn
 real(dp)            :: mult_fcnp
 if(psi<0._dp)then
  mult_fcn = (-m*n*alpha*(alpha*psi)**(n-1._dp)) * ( 1._dp + (psi*alpha)**n )**(-1._dp)
  mult_fcnp = -m*n*alpha*(n-1._dp)*alpha*(alpha*psi)**(n-2._dp)*( 1._dp + (psi*alpha)**n )**(-1._dp) - &
              ( n*alpha*(alpha*psi)**(n-1._dp)*(1._dp + (psi*alpha)**n)**(-2._dp) ) * ( -m*n*alpha*(alpha*psi)**(n-1._dp) )
  d2Theta_dPsi2 = mult_fcn * dTheta_dPsi(psi,alpha,theta_res,theta_sat,n,m) + &
                  mult_fcnp * ( volFracLiq(psi,alpha,theta_res,theta_sat,n,m) - theta_res )
 else
  d2Theta_dPsi2 = 0._dp
 end if
 end function d2Theta_dPsi2
 

 ! ******************************************************************************************************************************
 ! public function dTheta_dTk: differentiate the freezing curve w.r.t. temperature
 ! ******************************************************************************************************************************
 function d2Theta_dTk2(Tk,theta_res,theta_sat,alpha,n,m)
 implicit none
 real(dp),intent(in) :: Tk            ! temperature (K)
 real(dp),intent(in) :: theta_res     ! residual liquid water content (-)
 real(dp),intent(in) :: theta_sat     ! porosity (-)
 real(dp),intent(in) :: alpha         ! vGn scaling parameter (m-1)
 real(dp),intent(in) :: n             ! vGn "n" parameter (-)
 real(dp),intent(in) :: m             ! vGn "m" parameter (-)
 real(dp)            :: d2Theta_dTk2    ! derivative of the freezing curve w.r.t. temperature (K-1)
 ! local variables
 real(dp)            :: kappa         ! constant (m K-1)
 real(dp)            :: xtemp         ! alpha*kappa*(Tk-Tfreeze) -- dimensionless variable (used more than once)
 ! compute kappa (m K-1)
 kappa =  LH_fus/(gravity*Tfreeze)    ! NOTE: J = kg m2 s-2
 ! define a tempory variable that is used more than once (-)
 xtemp = alpha*kappa*(Tk-Tfreeze)
 ! differentiate the freezing curve w.r.t. temperature -- making use of the chain rule
 d2Theta_dTk2 = (-alpha*kappa*m*n*alpha*kappa)* (theta_sat - theta_res) * (  (n-1)*xtemp**(n - 2._dp) * (1._dp + xtemp**n)**(-m - 1._dp) &
                                                                        + n*(-m-1)*xtemp**(2*n - 2._dp) * (1._dp + xtemp**n)**(-m - 2._dp) )
 end function d2Theta_dTk2


end module soil_utils_fida_module