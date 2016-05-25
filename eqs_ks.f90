module EQS_KS

use KINDS, only: ik,qk,dk
implicit none

contains

subroutine DKS_RHS(neq,s,u,udot)
! Description:
!    Computes the value of the right-hand side of the equations of motion of the Kustaanheimo-
!    Stiefel formulation.
!
! ==============================================================================
!                                  VARIABLES AND DECLARATIONS
! ==============================================================================

! MODULES
use AUXILIARIES,    only: inSoI
use TRANSFORM,      only: DKS2CART,DKSMAT
use THIRDBODIES,    only: DACC3B_ND
use CONSTANTS,      only: DU,TU,wEarth,smaEarth,muSun,muEarth

! VARIABLES
implicit none
! Arguments
integer(ik),intent(in)   ::  neq                  ! Number of equations.
real(dk),intent(in)      ::  s                    ! Value of fictitious time.
real(dk),intent(in)      ::  u(1:neq)             ! KS state vector.
real(dk),intent(out)     ::  udot(1:neq)          ! RHS of EoM's, ND.
! Local variables
! -- State variables
real(dk)    ::  x(1:4),xdot(1:4)    ! Radius and velocity in R^4, ND
real(dk)    ::  r                   ! Radius magnitude, ND
real(dk)    ::  t                   ! Physical time, ND
real(dk)    ::  L(1:4,1:4)          ! KS-matrix
! -- Perturbations
real(dk)    ::  rVpot               ! Perturbing potential term
real(dk)    ::  drVdu(1:4)          ! Derivatives of the perturbing potential term
real(dk)    ::  P(1:4)              ! Non-potential perturbing acceleration, ND
real(dk)    ::  r2(1:3)             ! Position vector of the third body [-]
real(dk)    ::  lon_Earth           ! Current Earth longitude [rad]
real(dk)    ::  Q(1:4)              ! Total acceleration, ND

! Debug
integer(ik)  ::  i


! ==============================================================================


! STATE VECTOR DICTIONARY
! u(1:4)        u1,...,u4; KS-position, in R^4
! u(5:8)        u1',...,u4'; KS-velocity, in R^4
! u(9)          h; (-total energy) = (-Keplerian energy) + (-potential)
! u(10)         t; non-dimensional physical time

! ==============================================================================
! 01. COMPUTE CARTESIAN COORDINATES, ND
! ==============================================================================

call DKS2CART(u,x,xdot)
r = sqrt(dot_product(x,x))
t = u(10)
!!t0 = t0Prop*TU

! ==============================================================================
! 02. COMPUTE PERTURBING POTENTIAL, ND
! ==============================================================================

! Initialize
rVpot = 0._dk; drVdu = 0._dk

! ==============================================================================
! 03. COMPUTE PERTURBING ACCELERATION, ND
! ==============================================================================

! Initialize
P = 0._dk
lon_Earth = wEarth*(t/TU)

if (inSoI) then
  ! Third body: Sun
  r2 = -(smaEarth/DU)*[cos(lon_Earth),sin(lon_Earth),0._dk]
  P(1:3) = DACC3B_ND(x(1:3),r2,muEarth,muSun)
  
else
  ! Third body: Earth
  r2 = (smaEarth/DU)*[cos(lon_Earth),sin(lon_Earth),0._dk]
  P(1:3) = DACC3B_ND(x(1:3),r2,muSun,muEarth)
  
end if

! Compute KS-matrix
L = DKSMAT(u)

! Total acceleration
Q = -drVdu/4._dk + .5_dk*r*matmul(transpose(L),P)

! ==============================================================================
! 04. COMPUTE RHS
! ==============================================================================

! Velocities
udot(1:4) = u(5:8)

! Accelerations
udot(5:8) = -.5_dk*u(9)*u(1:4) + Q

! Total energy
udot(9)   = -2._dk*dot_product(u(5:8),matmul(transpose(L),P))

! Time
udot(10)  = r

end subroutine DKS_RHS
    
end module EQS_KS
