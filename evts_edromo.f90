module EVTS_EDROMO

use KINDS, only: ik,dk,qk

contains

subroutine DEDROMO_EVT(neq,phi,z,ng,roots)
! Description:
!    Finds roots to stop the integration for the Dromo formulation.
! 
! ==============================================================================

! MODULES
use AUXILIARIES, only: JD_next,JD_stop
use CONSTANTS,   only: TU,secsPerDay
use SETTINGS,    only: flag_time_EDr

! VARIABLES
implicit none
! Arguments IN
integer(ik),intent(in)   ::  neq
integer(ik),intent(in)   ::  ng
real(dk),intent(in)      ::  phi
real(dk),intent(in)      ::  z(1:neq)
! Arguments OUT
real(dk),intent(out)     ::  roots(1:ng)

! Locals
real(dk)  ::  zeta
real(dk)  ::  t         ! Current time [-]

! Debug
integer,save  ::  called = 1

! ==============================================================================

! Get time
zeta  = z(1)*sin(phi) - z(2)*cos(phi)
if ( flag_time_EDr == 0_ik ) then
    ! Physical time	
    t = z(8)
elseif  ( flag_time_EDr == 1_ik ) then
    ! Constant time element
    t = z(8) - z(3)**1.5_dk*(zeta - phi)
elseif  ( flag_time_EDr == 2_ik ) then
    ! Linear time element
    t = z(8) - z(3)**1.5_dk*zeta
end if

! ==============================================================================
! 01. Next timestep
! ==============================================================================

roots(1) = t - JD_next*secsPerDay*TU

! ==============================================================================
! 02. Stop integration
! ==============================================================================

roots(2) = t - JD_stop*secsPerDay*TU

end subroutine DEDROMO_EVT

end module EVTS_EDROMO
