module TEST_PROP

use KINDS, only: ik,dk,qk

implicit none

contains

subroutine DP_TRAJECTORY(R_i,V_i,JD_i,JD_f,eqs,integ,tol,Yt,idiag,rdiag)

! MODULES
use AUXILIARIES,  only: inSoI,ind_time
use CONSTANTS,    only: secsPerDay,muSun,muEarth,DU,TU,wEarth,smaEarth
use CONSTANTS,    only: SET_UNITS
use STATE_INIT,   only: DINIT_KS,DFTIME_2D,DINIT_EDROMO,DHYPAN,DINIT_GDROMO
use SETTINGS,     only: flag_time_EDr,flag_time_GDr
use EQS_COWELL,   only: DCOWELL_RHS,DCOWELL_RHS_T2
use EQS_KS,       only: DKS_RHS
use EQS_EDROMO,   only: DEDROMO_RHS
use EQS_GDROMO,   only: DGDROMO_RHS
use EVTS_COWELL,  only: DCOWELL_EVT
use EVTS_KS,      only: DKS_EVT
use EVTS_EDROMO,  only: DEDROMO_EVT
use EVTS_GDROMO,  only: DGDROMO_EVT
use INTEGRATION,  only: DINTEGRATE
use STEPSIZES,    only: DSTEPSIZE
use THIRDBODIES,  only: QPOS_VEL_CIRC

! VARIABLES
implicit none
! Arguments IN
real(qk),intent(in)  ::  R_i(1:3),V_i(1:3)    ! Initial position and velocity [km,km/s]
real(qk),intent(in)  ::  JD_i,JD_f            ! Initial and final times [days]
integer,intent(in)   ::  eqs,integ
real(dk),intent(in)  ::  tol                  ! Integration tolerance
! Arguments OUT
real(dk),allocatable,intent(out) ::  Yt(:,:)  ! Trajectory array [s,km,km/s]
integer,intent(out)  ::  idiag(:)
real(dk),intent(out) ::  rdiag(:)
! States and trajectories
real(dk),allocatable ::  X_i(:),Xdot_i(:)
real(qk),allocatable ::  X_i_q(:),Xdot_i_q(:)
real(dk)             ::  s_i,ds_i
real(dk)             ::  t_D,r_helio(1:3),v_helio(1:3)
! Integrator settings - LSODAR
integer  ::  neq
integer,parameter  ::  itol = 1, itask = 1, iopt = 1
integer,parameter  ::  lrw = 300, liw = 300
integer,parameter  ::  jt = 2
integer            ::  nevts,jroot
real(dk)           ::  rwork(1:lrw)
integer            ::  iwork(1:liw)
! Diagnostics
integer  ::  istate,intstep
! Misc.
real(qk)  ::  mu

! INTERFACES AND PROCEDURE POINTERS
external SLSODAR


! ==============================================================================


! Variable glossary
! eqs:          flag for the equations of motion. eqs = 1 for Cowell,
!               = 2 for K-S, = 3 for EDromo.
!
! integ:        flag for the integrator. integ = 1 for LSODAR, = 2 for Radau.
!
! idiag:        integration diagnostics (integer).
!   (1):        total number of evaluations of the RHS.
!   (2):        istate at the end of propagation. Check the meaning of istate
!               in integration.f90, opksmain.f90, xeverhart.f90.
!   (3):        total number of internal integration steps.
! 
! rdiag:        integration diagnostics (real).
!   (1):        for LSODAR, value of time (JD) which the first method
!               switch to BDF was detected.
!   (2):        Average order of integration. The order is averaged
!               over the output steps instead of the internal integration steps.
!               Thus, the average could be unreliable for a small number of
!               output steps. For XRA15 rdiag(2) = 15.
!   (3):        for LSODAR, in case in which an exception was detected, rdiag(3)
!               is the last value of the physical time at which the integration
!               was successful.


! ==============================================================================
! 01. INITIALIZATIONS
! ==============================================================================

! Set gravitational parameter for non-dimensionalization
if (inSoI) then
  mu = muEarth
else
  mu = muSun
end if

call SET_UNITS(R_i,mu)

! Initialize diagnostics and deallocate stuff
idiag = 0; rdiag = 0._dk
if (allocated(X_i)) deallocate(X_i); if (allocated(X_i_q)) deallocate(X_i_q)
if (allocated(Xdot_i)) deallocate(Xdot_i); if (allocated(Xdot_i_q)) deallocate(Xdot_i_q)

! Select the initialization routine according to the type of variables
! 1: Cowell (n.d. Cartesian coordinates)
! 2: Kustaanheimo-Stiefel
! 3: EDromo
! 4: GDromo

ds_i = DSTEPSIZE((R_i/DU),(V_i/(DU*TU)),1._qk,DU,eqs,inSoI,tol)

! (A): STATE VECTOR AND INDEPENDENT VARIABLE INITIALIZATION
! (B): INTEGRATION LOOP
select case (eqs)

case(-1) ! ONLY RADAU: 3 2nd-order equations
  ! (A)
  neq = 3
  allocate(X_i(1:neq)); allocate(Xdot_i(1:neq))
  X_i = real((R_i/DU),dk)
  Xdot_i = real((V_i/(DU*TU)),dk)
  s_i = JD_i*secsPerDay*TU
  
  ! (B)
  call DINTEGRATE(FAKE_RHS_T1,DCOWELL_RHS_T2,DCOWELL_EVT,X_i,Xdot_i,s_i,ds_i,&
  &neq,real(JD_f,dk),eqs,integ,tol,Yt,idiag,rdiag)

case (1)
  ! LSODAR/Radau: 6 1st-order equations
  ! (A)
  ! Initialize state vector and independent var.
  neq = 6
  allocate(X_i(1:neq)); allocate(Xdot_i(1:neq))
  X_i = real([R_i/DU,V_i/(DU*TU)],dk)
  s_i = JD_i*secsPerDay*TU
  
  ! (B)
  call DINTEGRATE(DCOWELL_RHS,FAKE_RHS_T2,DCOWELL_EVT,X_i,Xdot_i,s_i,ds_i,neq,&
  &real(JD_f,dk),eqs,integ,tol,Yt,idiag,rdiag)

case (2)
  ! K-S is always initialized as a set of 1st-order equations since it is
  ! necessary to integrate the equations for time and ang. momentum, which are
  ! first order.
  ! Initialization is performed in QUAD as to prevent rare round-off problems,
  ! same for EDromo.
  ! (A)
  neq = 10
  allocate(X_i(1:neq)); allocate(Xdot_i(1:neq)); allocate(X_i_q(1:neq))
  call DINIT_KS(R_i,V_i,JD_i*secsPerDay,mu,DU,TU,X_i_q,0._qk)
  X_i = X_i_q
  s_i = DFTIME_2D(R_i,V_i,mu)
  
  ! (B)
  !!! CHECK THIS
  call DINTEGRATE(DKS_RHS,FAKE_RHS_T2,DKS_EVT,X_i,Xdot_i,s_i,ds_i,neq,&
      &real(JD_f,dk),eqs,integ,tol,Yt,idiag,rdiag)
  !!! CHECK THIS
  
case (3)
  ! (A)
  neq = 8
  allocate(X_i(1:neq)); allocate(Xdot_i(1:neq)); allocate(X_i_q(1:neq))
  s_i = DFTIME_2D(R_i,V_i,mu)
  call DINIT_EDROMO(R_i,V_i,JD_i*secsPerDay,DU,TU,X_i_q,real(s_i,qk),0._qk,flag_time_EDr)
  X_i = X_i_q
  ! SANITY CHECK FOR NaNs
  if (any(X_i /= X_i)) then
    idiag(2) = -12
    rdiag(3) = JD_i
    return
  end if
  
  ! (B)
  call DINTEGRATE(DEDROMO_RHS,FAKE_RHS_T2,DEDROMO_EVT,X_i,Xdot_i,s_i,ds_i,neq,&
      &real(JD_f,dk),eqs,integ,tol,Yt,idiag,rdiag)
  
case (4)
  ! (A)
  neq = 8
  allocate(X_i(1:neq)); allocate(Xdot_i(1:neq)); allocate(X_i_q(1:neq))
  s_i = DHYPAN((R_i/DU),(V_i/(DU*TU)))
  call DINIT_GDROMO(R_i,V_i,JD_i*secsPerDay,DU,TU,X_i_q,real(s_i,qk),0._qk,flag_time_GDr)
  X_i = X_i_q
  ! SANITY CHECK FOR NaNs
  if (any(X_i /= X_i)) then
    idiag(2) = -12
    rdiag(3) = JD_i
    return
  end if
  
  ! (B)
  call DINTEGRATE(DGDROMO_RHS,FAKE_RHS_T2,DGDROMO_EVT,X_i,Xdot_i,s_i,ds_i,neq,&
      &real(JD_f,dk),eqs,integ,tol,Yt,idiag,rdiag)
  
end select

contains

subroutine FAKE_RHS_T1(neq,t,y,ydot)
implicit none
integer(ik),intent(in)    ::  neq
real(dk),intent(in)       ::  t
real(dk),intent(in)       ::  y(1:neq)
real(dk),intent(out)      ::  ydot(1:neq)
end subroutine FAKE_RHS_T1

subroutine FAKE_RHS_T2(neq,t,y,ydot,yddot)
implicit none
integer(ik),intent(in)    ::  neq         
real(dk),intent(in)       ::  t           
real(dk),intent(in)       ::  y(1:neq)    
real(dk),intent(in)       ::  ydot(1:neq) 
real(dk),intent(out)      ::  yddot(1:neq)
end subroutine FAKE_RHS_T2

subroutine FAKE_EVT(neq,t,y,ng,roots)
implicit none
integer(ik),intent(in)  ::  neq,ng     
real(dk),intent(in)     ::  t,y(1:neq) 
real(dk),intent(out)    ::  roots(1:ng)
end subroutine FAKE_EVT

end subroutine DP_TRAJECTORY

end module TEST_PROP
