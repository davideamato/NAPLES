module INTEGRATION

use KINDS, only: ik,qk,dk
implicit none

real(dk)  ::  eval_ratio,er_avg
integer   ::  called = 0

contains

subroutine DINTEGRATE(DRHS_T1,DRHS_T2,DEVT,X_i,Xdot_i,s_i,ds_i,neq,JD_f,eqs,&
&integ,tol,Yt,idiag,rdiag)
! Integrates the equations of motion starting from initial values of the state
! vector and independent variable "X_i, s_i" until time "JD_stop".
! It chooses between LSODAR and Radau as integrators and uses any of the
! following formulations:
! eqs = 1: Cowell
! eqs = 2: Kustaanheimo-Stiefel
! eqs = 3: EDromo
! eqs = 4: GDromo
! 
! Double precision version.

! MODULES
use AUXILIARIES, only: meth_switch,time_steps,ind_time,JD_next,JD_stop,inSoI
use CONSTANTS,   only: DU,TU,secsPerDay
use SETTINGS,    only: mxsteps,dt_H,dt_CE
use XEVERHART

! VARIABLES
implicit none
! Arguments IN
integer,intent(in)    ::  neq,eqs,integ
real(dk),intent(in)   ::  X_i(1:neq),Xdot_i(1:neq),s_i,ds_i,tol,JD_f
! Arguments OUT       
integer,intent(out)   ::  idiag(:)
real(dk),intent(out),allocatable  ::  Yt(:,:)
real(dk),intent(out)  ::  rdiag(:)

! INTERFACES AND PROCEDURE POINTERS
external SLSODAR
abstract interface
  subroutine FTYPE1(neq,t,y,ydot)
    import  ::  ik,dk
    implicit none
    integer(ik),intent(in)    ::  neq
    real(dk),intent(in)       ::  t
    real(dk),intent(in)       ::  y(1:neq)
    real(dk),intent(out)      ::  ydot(1:neq)
  end subroutine FTYPE1
  
  subroutine FTYPE2(neq,t,y,ydot,yddot)
    import  ::  ik,dk
    implicit none
    integer(ik),intent(in)    ::  neq
    real(dk),intent(in)       ::  t
    real(dk),intent(in)       ::  y(1:neq)
    real(dk),intent(in)       ::  ydot(1:neq)
    real(dk),intent(out)      ::  yddot(1:neq)
  end subroutine FTYPE2
  
  subroutine EVENTS(neq,t,y,ng,roots)
    import  ::  ik,dk
    implicit none
    integer(ik),intent(in)  ::  neq,ng
    real(dk),intent(in)     ::  t,y(1:neq)
    real(dk),intent(out)    ::  roots(1:ng)
  end subroutine EVENTS
  
end interface
procedure(FTYPE1)  ::  DRHS_T1
procedure(FTYPE2)  ::  DRHS_T2
procedure(EVENTS)  ::  DEVT

! Integration variables
! Integrator settings - LSODAR
integer,parameter  ::  itol = 1, itask = 1, iopt = 1
integer,parameter  ::  lrw = 300, liw = 300
integer,parameter  ::  jt = 2
integer            ::  nevts
integer,allocatable::  jroot(:)
real(dk)           ::  rwork(1:lrw)
integer            ::  iwork(1:liw)
logical            ::  next_t
integer            ::  steps
! Diagnostics
integer  ::  istate,intstep
real(dk) ::  sum_ord
logical  ::  xflag

! RA15 variables
real(dk)  :: s_end,h0,tstar
real(dk),allocatable  ::  tip(:),Xt(:,:)
integer   :: nclass,mip,temp,itaskr,idiagr(1:10)
integer   :: ii

! States and trajectories
real(dk) ::  X(1:neq),Xdot(1:neq),s,ds
real(dk) ::  Y(1:7)
real(dk),allocatable  ::  auxY(:,:)

! Debug
integer  ::  debug_int,it
integer,save :: jj = 1

! ==============================================================================

! Variable glossary
! DRHS_T1:      pointer to Type I equations of motion (y' = F(y,t)).
! 
! DRHS_T2:      pointer to Type II equations of motion (y'' = F(y,ydot,t)).
! 
! eqs:          flag for the equations of motion. eqs = 1 for Cowell,
!               = 2 for K-S, = 3 for EDromo.
!
! integ:        flag for the integrator. integ = 1 for LSODAR, = 2 for Radau.
!
! state_end:    final position, velocity and time [s,km,km/s]
! 
! idiag:        integration diagnostics (integer).
!   (1):        total number of evaluations of the RHS.
!   (2):        istate at the end of propagation.
!   (3):        total number of internal integration steps.
! 
! rdiag:        integration diagnostics (real).
!   (1):        for LSODAR, value of time (JD) which the first method
!               switch to BDF was detected.
!   (2):        average order of integration. The order is averaged over the 
!               output steps instead of the internal integration steps.
!               Thus, the average could be unreliable for a small number of
!               output steps.
!               When using XRA15, rdiag(2) = 15.
!   (3):        for LSODAR, in case in which an exception was detected, rdiag(3) 
!               is the last value of the physical time (JD) at which the 
!               integration was successful.

debug_int = 0

! ==============================================================================
! 01. INITIALIZE INTEGRATOR VARIABLES
! ==============================================================================

! Initialize state vector, i. variable, stepping
X = X_i
Xdot = Xdot_i
s = s_i
ds = huge(0._dk)
! ind_time was set in the main program
JD_next = time_steps(ind_time)/secsPerDay

if (allocated(auxY)) deallocate(auxY)
allocate(auxY(1:2*mxsteps,1:7)); auxY = 0._dk
auxY(1,1:7) = DSTATE_AND_TIME(neq,eqs,X,Xdot,s,DU,TU,0._dk)

! Select the integrator depending on the "integ" flag.
! 1: LSODAR
! 2: RADAU
select case (integ)

case (1)
  ! Initialize LSODAR settings and diagnostics
  istate = 1; iwork = 0; rwork = 0._dk; sum_ord = 0._dk
  
  ! Set maximum number of integration steps per output step and initial
  ! step size
  iwork(6) = 200000
  rwork(5) = ds_i
  
  ! Set counter and method switch flag
  intstep = 0
  meth_switch = .false.
  next_t = .false.
  
  ! The exit condition is always found by root-finding, also for Cowell. This
  ! greatly simplifies implementation since in general we use 4 (slightly)
  ! different time steps due to the backwards propagation.
  nevts = 2
  if (allocated(jroot)) deallocate(jroot)
  allocate(jroot(1:nevts)); jroot = 0
  
  ! ============================================================================
  ! 02. LSODAR INTEGRATION LOOP
  ! ============================================================================
  
  do
    ! Advance the solution
    call SLSODAR(DRHS_T1,neq,X,s,s+ds,itol,tol,tol,itask,istate,iopt,rwork,&
    &lrw,iwork,liw,FAKE,jt,DEVT,nevts,jroot)
!    call DLSODAR(DRHS_T1,neq,X,s,s+ds,itol,tol,tol,itask,istate,iopt,rwork,&
!    &lrw,iwork,liw,FAKE,jt,DEVT,nevts,jroot) ! (for quad)
    if (istate < 0) then
      call CLEANUP()
      idiag(2) = istate
      return
    end if
    
    intstep = intstep + 1
    
    call SANITY_CHECKS(neq,X,intstep,mxsteps,istate,idiag,xflag)
    if (xflag) then 
      call CLEANUP()
      return
    end if
    
    auxY(intstep+1,:) = DSTATE_AND_TIME(neq,eqs,X,Xdot,s,DU,TU,0._dk)
    rdiag(3) = auxY(intstep+1,1)
    
    ! DIAGNOSTICS
    if (iwork(19) == 2 .and. .not.(meth_switch)) then
      meth_switch = .true.
      rdiag(1) = auxY(intstep+1,1)
    end if
    sum_ord = sum_ord + iwork(14)
    
    ! ADVANCE TIMESTEP
    next_t = (istate == 3 .and. (ind_time < size(time_steps,1) - 1))
    if (next_t) then
      ind_time = ind_time + 1
      JD_next = time_steps(ind_time)/secsPerDay
    end if
    
    ! EXIT CONDITIONS
    if (jroot(2) == 1) exit
    
  end do
  
  ! SAVE TO Yt AND DEALLOCATE AUXILIARY TRAJECTORY ARRAY
  if (allocated(Yt)) deallocate(Yt)
  if (allocated(auxY)) then
    allocate(Yt(1:intstep+1,1:7),source=auxY(1:intstep+1,1:7))
    deallocate(auxY)
  end if
  
  ! Save diagnostics
  idiag = [iwork(12),istate,iwork(11)]
  rdiag(2) = sum_ord/intstep
  
case (2)
  ! ============================================================================
  ! 03. XRA15 INTEGRATION (Loop is inside the subroutine).
  ! ============================================================================
  
  if (eqs == -1) then
    ! Cowell equations, Type II
    s_end = JD_f*secsPerDay*TU
    
    mip = 0
    nclass = 2
    itaskr = 1
    tstar = 0._dk
    
    istate = 0
    idiag = 0; idiagr = 0
  
    ! The integration loop is performed inside XRA15.
    call XRA15(nclass,neq,eqs,itaskr,DRHS_T2,tol,mip,X,Xdot,ds_i,s,s_end,tstar&
    &,tip,Xt,istate,idiagr)
    if (istate < 0) then
      call CLEANUP()
      idiag(1) = idiagr(1); idiag(2) = istate; idiag(3) = idiagr(4)
      rdiag(1) = 0._dk; rdiag(2) = 15._dk
      return
    end if
    
    ! Dimensionalization
    if (allocated(Yt)) deallocate(Yt)
    allocate(Yt(1:size(Xt,1),1:7))
    do ii=1,size(Yt,1)
      Yt(ii,:) = DSTATE_AND_TIME(neq,eqs,Xt(ii,2:4),Xt(ii,5:7),Xt(ii,1),DU,TU,0._dk)
    end do
    
  else
    ! KS/EDromo/GDromo equations, Type I
    s_end = huge(0._dk)
    
    mip = 0
    nclass = 1
    itaskr = 3
    tstar = JD_f*secsPerDay*TU
    
    idiag = 0; idiagr = 0
    
    ! The integration loop is performed inside XRA15.
    call XRA15(nclass,neq,eqs,itaskr,DRHS_T1,tol,mip,X,Xdot,ds_i,s,s_end,tstar&
    &,tip,Xt,istate,idiagr)
    if (istate < 0) then
      call CLEANUP()
      idiag(1) = idiagr(1); idiag(2) = istate; idiag(3) = idiagr(4)
      rdiag(1) = 0._dk; rdiag(2) = 15._dk
      return
    end if
    
    ! Conversion to dimensional Cartesian coordinates
    if (allocated(Yt)) deallocate(Yt)
    allocate(Yt(1:size(Xt,1),1:7))
    do ii=1,size(Yt,1)
    ! Xdot is not used here (the equations are Type I), nevertheless we pass it
    ! to comply with the interface of DSTATE_AND_TIME.
      Yt(ii,:) = DSTATE_AND_TIME(neq,eqs,Xt(ii,2:),Xdot,Xt(ii,1),DU,TU,0._dk)
    end do
  
  end if
 
  ! PROCESSING
  ! Move the total number of steps to idiag(3). Overwrite first element of rdiag
  ! to avoid confusion with LSODAR.
  idiag(1) = idiagr(1); idiag(2) = istate; idiag(3) = idiagr(4)
  rdiag(1) = 0._dk; rdiag(2) = 15._dk
  
  ! Compute the maximum ratio between function evaluations due to root-finding
  ! and to stepping.
  called = called + 1
  eval_ratio = real(idiagr(3),dk)/real(idiagr(1),dk)
  er_avg = er_avg + eval_ratio
!  er_avg = er_avg/called
  
end select


contains

    subroutine CLEANUP()
      if (allocated(jroot)) deallocate(jroot)
      if (allocated(auxY))  deallocate(auxY)
      if (allocated(Yt)) deallocate(Yt)
    end subroutine CLEANUP
    
    subroutine FAKE()
    end subroutine FAKE

end subroutine DINTEGRATE


! ==============================================================================
! AUXILIARY PROCEDURES
! ==============================================================================

function DSTATE_AND_TIME(neq,eqs,X,Xdot,s,DU,TU,pot)
! Convert (X,s) to (R,V,t) for any formulation.

use SETTINGS,  only: flag_time_EDr,flag_time_GDr
use TRANSFORM, only: DKS2CART,DEDROMO2CART,DEDROMO_TE2TIME,DGDROMO2CART,&
&DGDROMO_TE2TIME
implicit none
! Inputs
integer,intent(in)   :: neq,eqs
real(dk),intent(in)  :: X(1:neq),Xdot(1:neq),s,pot
real(qk),intent(in)  :: DU,TU
! Function definition
real(dk)  ::  DSTATE_AND_TIME(1:7)

! Locals
real(dk)  ::  r(1:3),v(1:3),t


! ==============================================================================

select case (eqs)
case(-1)
  DSTATE_AND_TIME = [s/TU,X*DU,Xdot*DU*TU]
  
case(1)
  DSTATE_AND_TIME = [s/TU,X(1:3)*DU,X(4:6)*DU*TU]
  
case(2)
  call DKS2CART(X,r,v)
  DSTATE_AND_TIME = [X(10)/TU,r*DU,v*DU*TU]

case(3)
  call DEDROMO2CART(X,s,0._dk,r,v)
  t = DEDROMO_TE2TIME(X,s,flag_time_EDr)
  DSTATE_AND_TIME = [t/TU,r*DU,v*DU*TU]

case(4)
  call DGDROMO2CART(X,s,0._dk,r,v)
  t = DGDROMO_TE2TIME(X,s,flag_time_GDr)
  DSTATE_AND_TIME = [t/TU,r*DU,v*DU*TU]
  
end select

end function DSTATE_AND_TIME


subroutine SANITY_CHECKS(neq,X,i_int,mxsteps,istate,idiag,xflag)

implicit none
! Inputs
integer,intent(in)  ::  neq,istate,i_int,mxsteps
real(dk),intent(in) ::  X(1:neq)
! Outputs
integer,intent(out) ::  idiag(:)
logical  ::  xflag

! ==============================================================================

xflag = .false.
if (any(X /= X)) then
  idiag(2) = -12; xflag = .true.
else if (istate == 2) then
  idiag(2) = -13; xflag = .true.
else if (i_int == mxsteps) then
  idiag(2) = -11; xflag = .true.
end if
end subroutine SANITY_CHECKS



!function EXIT_COND(eqs,integ,intstep,steps,jroot,nevts)

!implicit none
!! Inputs
!integer,intent(in)  ::  eqs,integ,intstep,steps,nevts
!integer,intent(in)  ::  jroot(1:nevts)
!! Function definition
!logical  ::  EXIT_COND

!! ==============================================================================

!select case (eqs)
!case (1)
!  EXIT_COND = (intstep == steps)

!case (2)
!  EXIT_COND = (jroot(2) == 1 )  ! TO DO: or jroot(3) == 1 (NOPE)
!  
!case (3)
!  EXIT_COND = (jroot(2) == 1)   ! TO DO: or jroot(3) == 1 (NOPE)

!end select

!end function EXIT_COND
!select case (eqs)
!case(1)
!!  ds = ((JD_f*secsPerDay*TU - s_i)/steps)
!  if (inSoI) then
!    ds = dt_CE*secsPerDay*TU
!  else
!    ds = dt_H*secsPerDay*TU
!  end if
!  steps = (JD_f*secsPerDay*TU - s_i)/ds

!case(2,3,4)
!end select
!  select case (eqs)
!  case(1)
!    nevts = 0
!  case(2,3,4)
!  end select
!    if (eqs /= 1) then

!    else
!      ! Adjust size of last time step for Cowell
!      if (intstep == steps - 1) ds = JD_f*secsPerDay*TU - s
!    end if
!    if (EXIT_COND(eqs,integ,intstep,steps,jroot,nevts)) exit
!write(*,*) 'DINTEGRATE: X_i = ',X_i
!write(*,*) 'DINTEGRATE: s_i = ',s_i
!  write(*,*) 'DINTEGRATE: jroot = ',jroot
!  write(*,*) 'DINTEGRATE: X = ',X
!  write(*,*) 'DINTEGRATE: s = ',s
!  write(*,*) 'DINTEGRATE: auxY = '
!  write(*,'(7(es22.15,1x))') ( auxY(ii,1:7), ii =1,intstep+1 )
!write(*,*) 'DINTEGRATE: JD_i = ',auxY(1,1)/86400._qk
!write(*,*) 'DINTEGRATE: JD_next = ',JD_next
!    write(*,*) 'istate = ',istate
!    if (istate > 0) then
!      select case (eqs)
!        case(-1)
!    write(*,*) 'istate = ',istate
!    write(*,*) size(xt,2)
!    write(30,'(9(es22.15,1x))') (xt(ii,:), ii=1,size(xt,1))
!    if (istate > 0) then
!      select case (eqs)
!        case(-1)
!          it = 1
!        case(2)
!          it = 11
!        case(3,4)
!          it = 9
!      end select
!      write(*,*) 'Initial JD = ',Xt(1,it)/TU/secsPerDay
!      write(*,*) 'Final JD = ',Xt(size(xt,1),it)/TU/secsPerDay
!      write(*,*) 'tstar (JD) = ',tstar/TU/secsPerDay
!    end if
!          it = 1
!        case(2)
!          it = 11
!        case(3,4)
!          it = 9
!      end select
!      write(*,*) 'Initial JD = ',Xt(1,it)/TU/secsPerDay
!      write(*,*) 'Final JD = ',Xt(size(xt,1),it)/TU/secsPerDay
!      write(*,*) 'tstar (JD) = ',tstar/TU/secsPerDay
!    end if
end module INTEGRATION
