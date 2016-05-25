module AUXILIARIES

use KINDS, only: dk,qk
implicit none

! TIMES
! Initial and final times
real(qk)  ::  JD_i,JD_f,JD_stop,JD_CA   ! Initial, final and close approach times [days]
real(dk)  ::  JD_next
! Time steps vector
real(qk),allocatable  ::  time_steps(:)
integer               ::  ind_time
! Indices of the time steps vector corresponding to the switch IN and switch OUT
! of the region of influence. They are saved here since they are used by both
! the main and the DINTEGRATION subroutine when eqs = 1.
!integer  ::  ind_in_bwd,ind_out_fwd
! Time step of the Cowell integration. It is used in DINTEGRATION
! when eqs = 1.
!real(qk) :: dt_Cow

! INTEGRATION DIAGNOSTICS
! LSODAR - BDF method flags
! meth_switch:   =.true. if the BDF method was used at least once in the dp
!                propagation.
logical      ::  meth_switch

! ONLINE TRAJECTORY MATCHING
logical   ::  inSoI
real(qk)  ::  RSwitch

end module AUXILIARIES

! MF:            method flag, signalling in which phase of the propagation the
!                switch to BDF was performed. MF = 1,2,3 for phases H1, CE, H2
!                respectively. MF is set in DINTEGRATE (INTEGRATION module) or
!                EXCEPTION_CHECK in the main program.
