module THIRDBODIES

use KINDS, only: ik,dk,qk
implicit none

contains

function QACC3B_ND(r3,r2,mu1,mu2)
! Description:
!    Computes non-dimensional perturbation acceleration due to a third body 'mu2' when the main
!    body is 'mu1'.
! 
! ==============================================================================================
!                                  VARIABLES AND DECLARATIONS
! ==============================================================================================

! VARIABLES
implicit none
! Arguments
real(qk),intent(in)  ::  r3(1:3)    ! Position vector of the particle, ND.
real(qk),intent(in)  ::  r2(1:3)    ! main->3body position vector, ND.
real(qk),intent(in)  ::  mu1,mu2    ! Gravitational parameters of main and third body [km,km^3/s].
! Function definition
real(qk)  ::  QACC3B_ND(1:3)
! Locals
real(qk)  ::  r23(1:3)              ! 3body->particle position vector, ND.
real(qk)  ::  r2N,r23N

! ==============================================================================================
!                                            EXECUTION
! ==============================================================================================

! Compute position wrt third body
r23(:) = r3(:) - r2(:)
! Compute ND acceleration
r2N = sqrt(dot_product(r2,r2))
r23N = sqrt(dot_product(r23,r23))
QACC3B_ND = - (mu2/mu1)*(r23/r23N**3 + r2/r2N**3)

end function QACC3B_ND

function DACC3B_ND(r3,r2,mu1,mu2)
! Description:
!    Computes non-dimensional perturbation acceleration due to a third body 'mu2' when the main
!    body is 'mu1'.
! 
! ==============================================================================================
!                                  VARIABLES AND DECLARATIONS
! ==============================================================================================

! VARIABLES
implicit none
! Arguments
real(dk),intent(in)  ::  r3(1:3)    ! Position vector of the particle, ND.
real(dk),intent(in)  ::  r2(1:3)    ! main->3body position vector, ND.
real(qk),intent(in)  ::  mu1,mu2    ! Gravitational parameters of main and third body [km,km^3/s].
! Function definition
real(dk)  ::  DACC3B_ND(1:3)
! Locals
real(dk)  ::  r23(1:3)              ! 3body->particle position vector, ND.
real(dk)  ::  r2N,r23N

! ==============================================================================================
!                                            EXECUTION
! ==============================================================================================

! Compute position wrt third body
r23(:) = r3(:) - r2(:)
! Compute ND acceleration
r2N = sqrt(dot_product(r2,r2))
r23N = sqrt(dot_product(r23,r23))
DACC3B_ND = - (mu2/mu1)*(r23/r23N**3 + r2/r2N**3)

end function DACC3B_ND

function QPOS_VEL_CIRC(t,w,sma)
! Computes the position and velocity of a planet in a circular orbit with mean
! motion "w", semi-major axis "sma" at time "t". Quadruple precision.
real(qk),intent(in)  ::  t,w,sma
real(qk)  ::  QPOS_VEL_CIRC(1:6)
real(qk)  ::  lon

lon = w*t
QPOS_VEL_CIRC(1:3) = sma*[cos(lon),sin(lon),0._qk]
QPOS_VEL_CIRC(4:6) = sma*w*[-sin(lon),cos(lon),0._qk]

end function QPOS_VEL_CIRC

function DPOS_VEL_CIRC(t,w,sma)
! Computes the position and velocity of a planet in a circular orbit with mean
! motion "w", semi-major axis "sma" at time "t". Double precision.
real(dk),intent(in)  ::  t,w,sma
real(dk)  ::  DPOS_VEL_CIRC(1:6)
real(dk)  ::  lon

lon = w*t
DPOS_VEL_CIRC(1:3) = sma*[cos(lon),sin(lon),0._dk]
DPOS_VEL_CIRC(4:6) = sma*w*[-sin(lon),cos(lon),0._dk]

end function DPOS_VEL_CIRC

end module THIRDBODIES
