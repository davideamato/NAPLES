module AUXILIARIES
! Description:
!    Contains auxiliary quantities used in the NAPLES program.
! 
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
! 
! ==============================================================================

use KINDS, only: dk,qk
implicit none

! TIMES
! Initial and final times
real(qk)  ::  JD_i,JD_f,JD_stop,JD_CA   ! Initial, final and close approach times [days]
real(dk)  ::  JD_next
! Time steps vector
real(qk),allocatable  ::  time_steps(:)
integer               ::  ind_time

! INTEGRATION DIAGNOSTICS
! LSODAR - BDF method flags
logical      ::  meth_switch

! ONLINE TRAJECTORY MATCHING
logical   ::  inSoI
real(qk)  ::  RSwitch

end module AUXILIARIES
