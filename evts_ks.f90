module EVTS_KS
! Description:
!    Contains event functions for the K-S formulation in double precision.
!
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
! 
! ==============================================================================

use KINDS, only: ik,dk,qk

contains

subroutine DKS_EVT(neq,s,u,ng,roots)
! Description:
!    Finds roots to stop the integration for the KS formulation.
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
real(dk),intent(in)      ::  s
real(dk),intent(in)      ::  u(1:neq)
! Arguments OUT
real(dk),intent(out)     ::  roots(1:ng)

! Locals
real(dk)  ::  t         ! Current time [-]

! ==============================================================================

t = u(10)

! ==============================================================================
! 01. Next timestep
! ==============================================================================

roots(1) = t - JD_next*secsPerDay*TU

! ==============================================================================
! 02. Stop integration
! ==============================================================================

roots(2) = t - JD_stop*secsPerDay*TU

end subroutine DKS_EVT
    
end module EVTS_KS
