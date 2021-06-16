module expIntegral_module
USE nrtype
implicit none
private
public::expint
contains

 ! Numerical recipes routines removed; use code from UEB-Veg

 ! ****************** EXPONENTIAL INTEGRAL FUNCTION *****************************************
 ! From UEB-Veg
 ! Computes the exponential integral function for the given value
 FUNCTION EXPINT (LAI)
 real(rkind) LAI
 real(rkind) EXPINT
 real(rkind) a0,a1,a2,a3,a4,a5,b1,b2,b3,b4
 real(rkind),parameter :: verySmall=tiny(1.0_rkind)     ! a very small number
 IF (LAI < verySmall)THEN
  EXPINT=1._rkind

 ELSEIF (LAI.LE.1.0) THEN
  a0=-.57721566_rkind
  a1=.99999193_rkind
  a2=-.24991055_rkind
  a3=.05519968_rkind
  a4=-.00976004_rkind
  a5=.00107857_rkind

  EXPINT = a0+a1*LAI+a2*LAI**2+a3*LAI**3+a4*LAI**4+a5*LAI**5 - log(LAI)

 ELSE
  a1=8.5733287401_rkind
  a2=18.0590169730_rkind
  a3=8.6347637343_rkind
  a4=.2677737343_rkind
  b1=9.5733223454_rkind
  b2=25.6329561486_rkind
  b3=21.0996530827_rkind
  b4=3.9584969228_rkind

  EXPINT=(LAI**4+a1*LAI**3+a2*LAI**2+a3*LAI+a4)/ &
      ((LAI**4+b1*LAI**3+b2*LAI**2+b3*LAI+b4)*LAI*exp(LAI))

 END IF
 RETURN
 END FUNCTION EXPINT

END MODULE expIntegral_module
