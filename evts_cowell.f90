module EVTS_COWELL

use KINDS, only: ik,dk,qk
implicit none

contains

subroutine QCOWELL_EVT(neq,t,y,ng,roots)
! Description:
!    Finds roots to stop the integration for the Cowell formulation.
! 
! ==============================================================================

! MODULES
use AUXILIARIES, only: RSwitch,inSoI
use CONSTANTS,   only: DU,TU,wEarth,smaEarth

! VARIABLES
implicit none
! Arguments IN
integer(ik),intent(in)   ::  neq
integer(ik),intent(in)   ::  ng
real(qk),intent(in)      ::  t
real(qk),intent(in)      ::  y(1:neq)
! Arguments OUT
real(qk),intent(out)     ::  roots(1:ng)

! Locals
real(qk)  ::  r_geo(1:3)        ! Geocentric position vector [-]
real(qk)  ::  r_Earth(1:3)      ! Heliocentric Earth position vector [-]
real(qk)  ::  lon_Earth         ! Heliocentric Earth longitude [rad]
real(qk)  ::  rMag              ! Geocentric distance [-]

! ==============================================================================

! ==============================================================================
! 01. Reaching of switch distance
! ==============================================================================
if (inSoI) then
    r_geo = y(1:3)
else
    lon_Earth = (wEarth/TU)*t
    r_Earth = (smaEarth/DU)*[cos(lon_Earth),sin(lon_Earth),0._qk]
    r_geo = y(1:3) - r_Earth
end if
rMag = sqrt(dot_product(r_geo,r_geo))

roots(1) = rMag - RSwitch/DU

end subroutine QCOWELL_EVT


subroutine DCOWELL_EVT(neq,t,y,ng,roots)
! Find roots to stop the integration at output points and at switch time for
! the Cowell formulation.
! 
! ==============================================================================

! MODULES
use AUXILIARIES, only: JD_next,JD_stop
use CONSTANTS,   only: TU,secsPerDay

! VARIABLES
implicit none
! Arguments IN
integer(ik),intent(in)   ::  neq
integer(ik),intent(in)   ::  ng
real(dk),intent(in)      ::  t
real(dk),intent(in)      ::  y(1:neq)
! Arguments OUT
real(dk),intent(out)     ::  roots(1:ng)

! ==============================================================================

! ==============================================================================
! 01. Next timestep
! ==============================================================================

roots(1) = t - JD_next*secsPerDay*TU

! ==============================================================================
! 02. Stop integration
! ==============================================================================

roots(2) = t - JD_stop*secsPerDay*TU

end subroutine DCOWELL_EVT

end module EVTS_COWELL
