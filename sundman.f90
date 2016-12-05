module SUNDMAN
! Description:
!    Contains the event functions for the regularized formulations when
!    integrated with the XRA15 solver, and auxiliary routines.
! 
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
! 
! ==============================================================================

use KINDS, only: dk,qk
implicit none

contains


function TIMETOUT(eqs,neq,tstar,s,x)
! Description:
!    Event function for output at a prescribed time using XRA15 with regularized
!    formulations. It computes the difference between the current time
!    (computed from the state vector of any regularized formulation) and the
!    output time tstar.
! 
! ==============================================================================

use SETTINGS,  only: flag_time_EDr,flag_time_GDr
use TRANSFORM, only: DEDROMO_TE2TIME,DGDROMO_TE2TIME

! INPUTS
integer,intent(in)   ::  eqs,neq
real(dk),intent(in)  ::  s,tstar,x(1:neq)

! OUTPUT
real(dk)  ::  TIMETOUT

! LOCALS
real(dk)  ::  t

! ==============================================================================

select case (eqs)
  case(-2)
    TIMETOUT = x(5) - tstar
  case(2)
    TIMETOUT = x(10) - tstar
  case(3)
    ! Compute time from EDromo state and subtract to tstar
    t = DEDROMO_TE2TIME(x,s,flag_time_EDr)
    TIMETOUT = t - tstar
  case(4)
    ! Compute time from GDromo state and subtract to tstar
    t = DGDROMO_TE2TIME(x,s,flag_time_GDr)
    TIMETOUT = t - tstar
end select

contains

  function CSUM(add)
  ! Kahan's summation algorithm (AKA compensated summation).
  ! Ref: https://en.wikipedia.org/wiki/Kahan_summation_algorithm
  
  real(dk),intent(in) :: add(:)
  real(dk) :: CSUM
  real(dk) :: c,y,t,summ
  integer :: i,length
  
  length = size(add,1)
  c = 0._dk
  summ = 0._dk
  
  do i=1,length
    y = add(i) - c
    t = summ + y
    c = (t - summ) - y
    summ = t
    
  end do
  CSUM = summ
  
  end function CSUM

end function TIMETOUT


function DSUND(eqs,neq,s,x,xdot)
! Description:
!   Computes the value of dt/ds for the regularized formulations. It is a driver
!   subroutine calling DSUND_KS, DSUND_EDROMO, DSUND_GDROMO.
! 
! ==============================================================================

! INPUTS
integer,intent(in)   ::  eqs,neq
real(dk),intent(in)  ::  s,x(1:neq),xdot(1:neq)
! OUTPUTS
real(dk)  ::  DSUND

select case (eqs)
  case(-2)
    DSUND = xdot(5)
  case(2)
    DSUND = DSUND_KS(neq,s,x)
  case(3)
    DSUND = DSUND_EDROMO(neq,s,x)
  case(4)
    DSUND = DSUND_GDROMO(neq,s,x)
end select

end function DSUND


function DDSUND(eqs,neq,s,x,xdot)
! Description:
!    Computes the value of d^2t/ds^2 for the regularized formulations. It is a
!    driver subroutine calling DDSUND_KS, DDSUND_EDROMO, DDSUND_GDROMO.
! 
! ==============================================================================

! INPUTS
integer,intent(in)   ::  eqs,neq
real(dk),intent(in)  ::  s,x(1:neq),xdot(1:neq)
! OUTPUTS
real(dk)  ::  DDSUND

select case (eqs)
  case(-2)
    DDSUND = 2._dk*dot_product(x(1:4),xdot(1:4))
  case(2)
    DDSUND = DDSUND_KS(neq,s,x)
  case(3)
    DDSUND = DDSUND_EDROMO(neq,s,x,xdot)
  case(4)
    DDSUND = DDSUND_GDROMO(neq,s,x,xdot)
end select

end function DDSUND


function DSUND_KS(neq,s,u)
! Gives the value of dt/ds for the KS formulation.

! INPUTS
integer,intent(in)   ::  neq
real(dk),intent(in)  ::  s,u(1:neq)

! OUTPUT
real(dk)  ::  DSUND_KS

! ==============================================================================

DSUND_KS = dot_product(u(1:4),u(1:4))

end function DSUND_KS


function DSUND_EDROMO(neq,s,l)
! Description:
!    Gives the value of dt/ds for the EDromo formulation.
!
! ==============================================================================

! INPUTS
integer,intent(in)  :: neq
real(dk),intent(in) :: s,l(1:neq)

! OUTPUT
real(dk) :: DSUND_EDROMO

! LOCALS
real(dk) :: rho

! ==============================================================================

rho = 1._dk - l(1)*cos(s) - l(2)*sin(s)
DSUND_EDROMO = (l(3)**1.5_dk)*rho

end function DSUND_EDROMO


function DSUND_GDROMO(neq,s,l)
! Description:
!    Gives the value of dt/ds for the GDromo formulation.
! 
! ==============================================================================

! INPUTS
integer,intent(in)  :: neq
real(dk),intent(in) :: s,l(1:neq)

! OUTPUT
real(dk) :: DSUND_GDROMO

! LOCALS
real(dk) :: rho

! ==============================================================================

rho = l(1)*cosh(s) + l(2)*sinh(s) - 1._dk
DSUND_GDROMO = (l(3)**1.5_dk)*rho

end function DSUND_GDROMO


function DDSUND_KS(neq,s,u)
! Description:
!    Gives the value of d^2t/ds^2 for the KS formulation.
!
! ==============================================================================

! INPUTS
integer,intent(in)   ::  neq
real(dk),intent(in)  ::  s,u(1:neq)

! OUTPUT
real(dk)  ::  DDSUND_KS

! ==============================================================================

DDSUND_KS = 2._dk*dot_product(u(1:4),u(5:8))

end function DDSUND_KS


function DDSUND_EDROMO(neq,s,l,ldot)
! Description:
!    Gives the value of d^2t/ds^2 for the EDromo formulation.
!
! ==============================================================================

! INPUTS
integer,intent(in)   ::  neq
real(dk),intent(in)  ::  s,l(1:neq),ldot(1:neq)

! OUTPUT
real(dk)  ::  DDSUND_EDROMO

! LOCALS
real(dk)  ::  cs,ss
real(dk)  ::  rho,zeta
real(dk)  ::  term_en,term_ecc

! ==============================================================================

cs = cos(s); ss = sin(s)

rho = 1._dk - l(1)*cs - l(2)*ss
zeta = l(1)*ss - l(2)*cs
term_en = 1.5_dk*sqrt(l(3))*rho*ldot(3)
term_ecc = (l(3)**1.5_dk)*(zeta - ldot(1)*cs - ldot(2)*ss)
DDSUND_EDROMO = term_en + term_ecc

end function DDSUND_EDROMO


function DDSUND_GDROMO(neq,s,l,ldot)
! Description:
!    Gives the value of d^2t/ds^2 for the GDromo formulation.
!
! ==============================================================================

! INPUTS
integer,intent(in)   ::  neq
real(dk),intent(in)  ::  s,l(1:neq),ldot(1:neq)

! OUTPUT
real(dk)  ::  DDSUND_GDROMO

! LOCALS
real(dk)  ::  chs,shs
real(dk)  ::  rho,zeta
real(dk)  ::  term_en,term_ecc

! ==============================================================================

chs = cosh(s); shs = sinh(s)

rho = l(1)*chs + l(2)*shs - 1._dk
zeta = l(1)*shs + l(2)*chs
term_en = 1.5_dk*sqrt(l(3))*rho*ldot(3)
term_ecc = (l(3)**1.5_dk)*(zeta + ldot(1)*chs + ldot(2)*shs)
DDSUND_GDROMO = term_en + term_ecc

end function DDSUND_GDROMO


end module SUNDMAN
