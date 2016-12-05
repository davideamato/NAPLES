module EVTS_GDROMO
! Description:
!    Contains event functions for the GDromo formulation in double precision.
!
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
! 
! ==============================================================================

use KINDS, only: ik,dk,qk

contains

subroutine DGDROMO_EVT(neq,phi,z,ng,roots)
! Description:
!    Finds roots to stop the integration for the GDromo formulation.
! 
! ==============================================================================

! MODULES
use AUXILIARIES, only: JD_next,JD_stop
use CONSTANTS,   only: TU,secsPerDay
use SETTINGS,    only: flag_time_GDr

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

! ==============================================================================

! Get time
zeta  = z(1)*sinh(phi) + z(2)*cosh(phi)
if ( flag_time_GDr == 0_ik ) then
	! Physical time	
    t = z(8)
elseif  ( flag_time_GDr == 1_ik ) then
    ! Constant time element
    t = z(8) + z(3)**1.5_dk*(zeta - phi)
elseif  ( flag_time_GDr == 2_ik ) then
    ! Linear time element
    t = z(8) + z(3)**1.5_dk*zeta
end if

! ==============================================================================
! 01. Next timestep
! ==============================================================================

roots(1) = t - JD_next*secsPerDay*TU

! ==============================================================================
! 02. Stop integration
! ==============================================================================

roots(2) = t - JD_stop*secsPerDay*TU

end subroutine DGDROMO_EVT
    
end module EVTS_GDROMO
