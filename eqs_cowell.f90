module EQS_COWELL
! Description:
!    Contains the subroutines for computing the RHS of the Cowell equations in
!    double- and quadruple-precision.
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

subroutine QCOWELL_RHS(neq,t,y,ydot)
! Description:
!    Computes the value of the right-hand side of the equations of motion of the 
!    Cowell formulation (quad precision).
!
! ==============================================================================
!                                  VARIABLES AND DECLARATIONS
! ==============================================================================

! MODULES
use AUXILIARIES, only: inSOI
use CONSTANTS,   only: smaEarth,wEarth,muEarth,muSun,DU,TU
use THIRDBODIES, only: QACC3B_ND

! VARIABLES

implicit none

! Arguments
integer,intent(in)     ::  neq             ! Number of equations.
real(qk),intent(in)        ::  t               ! Time, ND.
real(qk),intent(in)        ::  y(1:neq)        ! Cartesian state vector, ND.
real(qk),intent(out)       ::  ydot(1:neq)     ! RHS of EoM's, ND.

! Local variables
real(qk)          ::  rMag                      ! Magnitude of position vector. [-]
real(qk)          ::  ap(1:3)                   ! Perturbation acceleration. [-]
real(qk)          ::  r2(1:3)                   ! Position vector of the third body [-]
real(qk)          ::  lon_Earth                 ! Current Earth longitude [rad]

! ==============================================================================
!                                            EXECUTION
! ==============================================================================

rMag = sqrt(dot_product(y(1:3),y(1:3)))

! ==============================================================================
! 01. COMPUTE PERTURBATIONS
! ==============================================================================

ap = 0._qk
lon_Earth = wEarth*(t/TU)
if (inSOI) then
    ! Third body: Sun
    r2 = -(smaEarth/DU)*[cos(lon_Earth),sin(lon_Earth),0._qk]
    ap = QACC3B_ND(y(1:3),r2,muEarth,muSun)
else
    ! Third body: Earth
    r2 = (smaEarth/DU)*[cos(lon_Earth),sin(lon_Earth),0._qk]
    ap = QACC3B_ND(y(1:3),r2,muSun,muEarth)
end if

! ==============================================================================
! 02. EVALUATE RIGHT-HAND SIDE
! ==============================================================================

ydot(1:3) = y(4:6)
ydot(4:6) = -y(1:3)/rMag**3 + ap

end subroutine QCOWELL_RHS

subroutine DCOWELL_RHS(neq,t,y,ydot)
! Description:
!    Computes the value of the right-hand side of the equations of motion of
!    the Cowell formulation (double precision).
!
! ==============================================================================

! MODULES
use AUXILIARIES, only: inSoI
use CONSTANTS,   only: smaEarth,wEarth,muEarth,muSun,DU,TU
use THIRDBODIES, only: DACC3B_ND

! VARIABLES
implicit none

! Arguments
integer,intent(in)     ::  neq                 ! Number of equations.
real(dk),intent(in)        ::  t               ! Time, ND.
real(dk),intent(in)        ::  y(1:neq)        ! Cartesian state vector, ND.
real(dk),intent(out)       ::  ydot(1:neq)     ! RHS of EoM's, ND.

! Local variables
real(dk)          ::  rMag                     ! Magnitude of position vector. [-]
real(dk)          ::  ap(1:3)                  ! Perturbation acceleration. [-]
real(dk)          ::  r2(1:3)                  ! Position vector of the third body [-]
real(dk)          ::  lon_Earth                ! Current Earth longitude [rad]

! ==============================================================================

rMag = sqrt(dot_product(y(1:3),y(1:3)))

! ==============================================================================
! 01. COMPUTE PERTURBATIONS
! ==============================================================================

ap = 0._dk
lon_Earth = wEarth*(t/TU)
if (inSOI) then
    ! Third body: Sun
    r2 = -(smaEarth/DU)*[cos(lon_Earth),sin(lon_Earth),0._dk]
    ap = DACC3B_ND(y(1:3),r2,muEarth,muSun)
else
    ! Third body: Earth
    r2 = (smaEarth/DU)*[cos(lon_Earth),sin(lon_Earth),0._dk]
    ap = DACC3B_ND(y(1:3),r2,muSun,muEarth)
end if

! ==============================================================================
! 02. EVALUATE RIGHT-HAND SIDE
! ==============================================================================
    
ydot(1:3) = y(4:6)
ydot(4:6) = -y(1:3)/rMag**3 + ap

end subroutine DCOWELL_RHS



subroutine DCOWELL_RHS_T2(neq,t,y,ydot,yddot)
! Description:
!    Computes the value of the right-hand side of the equations of motion of the 
!    Cowell formulation. This version provides the right-hand-side of the
!    2nd-order equations of motion. Double precision.
! 
! ==============================================================================
!                                  VARIABLES AND DECLARATIONS
! ==============================================================================
!
! MODULES
use AUXILIARIES, only: inSoI
use CONSTANTS,   only: smaEarth,wEarth,muEarth,muSun,DU,TU
use THIRDBODIES, only: DACC3B_ND

! VARIABLES
implicit none

! Arguments
integer,intent(in)     ::  neq             ! Number of equations.
real(dk),intent(in)    ::  t               ! Time, ND.
real(dk),intent(in)    ::  y(1:neq)        ! Cartesian position, ND.
real(dk),intent(in)    ::  ydot(1:neq)     ! Cartesian velocity, ND.
real(dk),intent(out)   ::  yddot(1:neq)    ! Cartesian acceleration, ND.

! Local variables
real(dk)          ::  rMag                      ! Magnitude of position vector. [-]
real(dk)          ::  ap(1:3)                   ! Perturbation acceleration. [-]
real(dk)          ::  r2(1:3)                   ! Position vector of the third body [-]
real(dk)          ::  lon_Earth                 ! Current Earth longitude [rad]

! ==============================================================================
!                                            EXECUTION
! ==============================================================================

rMag = sqrt(dot_product(y,y))

! ==============================================================================
! 01. COMPUTE PERTURBATIONS
! ==============================================================================

ap = 0._dk
lon_Earth = wEarth*(t/TU)
if (inSOI) then
    ! Third body: Sun
    r2 = -(smaEarth/DU)*[cos(lon_Earth),sin(lon_Earth),0._dk]
    ap = DACC3B_ND(y,r2,muEarth,muSun)
else
    ! Third body: Earth
    r2 = (smaEarth/DU)*[cos(lon_Earth),sin(lon_Earth),0._dk]
    ap = DACC3B_ND(y,r2,muSun,muEarth)
end if

! ==============================================================================
! 02. EVALUATE RIGHT-HAND SIDE
! ==============================================================================
    
yddot = -y/rMag**3 + ap

end subroutine DCOWELL_RHS_T2



end module EQS_COWELL
