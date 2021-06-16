module spline_int_module
USE nrtype
implicit none
private
public::spline
public::splint
contains

 ! *************************************************************
 ! new subroutine: spline
 ! *************************************************************
 SUBROUTINE spline(x,y,yp1,ypn,y2,err,message)
 ! computes 2nd derivatives of the interpolating function at tabulated points
 IMPLICIT NONE
 ! dummy variables
 real(rkind), DIMENSION(:), INTENT(IN) :: x,y
 real(rkind), INTENT(IN) :: yp1,ypn
 real(rkind), DIMENSION(:), INTENT(OUT) :: y2
 integer(i4b),intent(out)      :: err
 character(*),intent(out)      :: message
 ! local variables
 character(len=128) :: cmessage
 INTEGER(I4B) :: n
 real(rkind), DIMENSION(size(x)) :: a,b,c,r
 ! initialize error control
 err=0; message="f-spline/"
 ! check that the size of the vectors match
 if(size(x)/=size(y) .or. size(y)/=size(y2)) then
  err=20; message="f-spline/sizeMismatch"; return
 else
  n=size(x)
 end if
 ! start procedure
 c(1:n-1)=x(2:n)-x(1:n-1)
 r(1:n-1)=6.0_rkind*((y(2:n)-y(1:n-1))/c(1:n-1))
 r(2:n-1)=r(2:n-1)-r(1:n-2)
 a(2:n-1)=c(1:n-2)
 b(2:n-1)=2.0_rkind*(c(2:n-1)+a(2:n-1))
 b(1)=1.0
 b(n)=1.0
 if (yp1 > 0.99e30_rkind) then
  r(1)=0.0
  c(1)=0.0
 else
  r(1)=(3.0_rkind/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  c(1)=0.5
 end if
 if (ypn > 0.99e30_rkind) then
  r(n)=0.0
  a(n)=0.0
 else
  r(n)=(-3.0_rkind/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
  a(n)=0.5
 end if
 call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n),err,cmessage)
 if (err/=0) message=trim(message)//trim(cmessage); return
 END SUBROUTINE spline

 ! *************************************************************
 ! new subroutine: splint
 ! *************************************************************
 SUBROUTINE splint(xa,ya,y2a,x,y,err,message)
 IMPLICIT NONE
 ! declare dummy variables
 real(rkind), DIMENSION(:), INTENT(IN)  :: xa,ya,y2a
 real(rkind), INTENT(IN)  :: x
 real(rkind), INTENT(OUT) :: y
 integer(i4b),intent(out)      :: err
 character(*),intent(out)      :: message
 ! declare local variables
 INTEGER(I4B) :: khi,klo,n
 real(rkind) :: a,b,h
 ! check size of input vectors
 if (size(xa)==size(ya) .and. size(ya)==size(y2a)) then
  n=size(xa)
 else
  err=20; message="f-splint/sizeMismatch"; return
 end if
 ! start procedure
 klo=max(min(locate(xa,x),n-1),1)
 khi=klo+1
 h=xa(khi)-xa(klo)
 if (h == 0.0_rkind) then; err=20; message="f-splint/badXinput"; return; end if
 a=(xa(khi)-x)/h
 b=(x-xa(klo))/h
 y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_rkind
 END SUBROUTINE splint

 ! *************************************************************
 ! new subroutine: locate
 ! *************************************************************
 FUNCTION locate(xx,x)
 IMPLICIT NONE
 real(rkind), DIMENSION(:), INTENT(IN) :: xx
 real(rkind), INTENT(IN) :: x
 INTEGER(I4B) :: locate
 INTEGER(I4B) :: n,jl,jm,ju
 LOGICAL :: ascnd
 n=size(xx)
 ascnd = (xx(n) >= xx(1))
 jl=0
 ju=n+1
 do
  if (ju-jl <= 1) exit
  jm=(ju+jl)/2
  if (ascnd .eqv. (x >= xx(jm))) then
   jl=jm
  else
   ju=jm
  end if
 end do
 if (x == xx(1)) then
  locate=1
 else if (x == xx(n)) then
  locate=n-1
 else
  locate=jl
 end if
 END FUNCTION locate

 ! *************************************************************
 ! new subroutine: tridag
 ! *************************************************************
 SUBROUTINE tridag(a,b,c,r,u,err,message)
 IMPLICIT NONE
 ! dummy variables
 real(rkind), DIMENSION(:), INTENT(IN) :: a,b,c,r
 real(rkind), DIMENSION(:), INTENT(OUT) :: u
 integer(i4b),intent(out)      :: err
 character(*),intent(out)      :: message
 ! local variables
 real(rkind), DIMENSION(size(b)) :: gam
 INTEGER(I4B) :: n,j
 real(rkind) :: bet
 ! initialize error control
 err=0; message="f-spline/OK"
 ! check that the size of the vectors match
 if ( all((/size(b),size(c)+1,size(r),size(u)/) == size(a)+1) ) then
  n=size(a)+1
 else
  err=20; message="f-tridag/sizeMismatch"; return
 end if
 ! start procedure
 bet=b(1)
 if (bet == 0.0_rkind) then; err=20; message="f-tridag/errorAtCodeStage-1"; return; end if
 u(1)=r(1)/bet
 do j=2,n
  gam(j)=c(j-1)/bet
  bet=b(j)-a(j-1)*gam(j)
  if (bet == 0.0_rkind) then; err=20; message="f-tridag/errorAtCodeStage-2"; return; end if
  u(j)=(r(j)-a(j-1)*u(j-1))/bet
 end do
 do j=n-1,1,-1
  u(j)=u(j)-gam(j+1)*u(j+1)
 end do
 END SUBROUTINE tridag

end module spline_int_module
