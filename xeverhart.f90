module XEVERHART
! Description:
!    This is a major re-implementation in modern Fortran of Everhart's integrator
!    RA15, with added features such as low-order dense output and event location. 
!    It is based on the implementation in modern Fortran of RA15 by Angelo Graziosi
!    which can be found at
! 
!    http://www.webalice.it/angelo.graziosi/Everhart-Integrator.pdf
!
! References:
! 
! [1] E. Everhart, An Efficient Integrator That Uses Gauss-Radau Spacings,
!     in A. Carusi and G. B. Valsecchi - Dynamics of Comets: Their Origin and
!     Evolution, 185-202. 1985 by D. Reidel Publishing Company.
! [2] E. Everhart, Implicit Single-Sequence Methods For Integrating Orbits,
!     Celestial Mechanics, 10, pp. 35-55, 1974.
! [3] D. Amato, G. Baù, C. Bombardelli, Orbit propagation of close encounters
!     with regularized formulations (provisional title), to be submitted to
!     MNRAS.
! [4] H. Rein, D.S. Spiegel, ... TODO
! 
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
! 
! ==============================================================================

use KINDS,   only: dp=>dk,qp=>qk
use SUNDMAN, only: DSUND,DDSUND,TIMETOUT
implicit none
private

! COMMON PARAMETERS
integer,parameter   ::  nor = 15
integer,parameter   ::  nsteps = 7, ncoef = (nsteps*(nsteps-1))/2
integer,parameter   ::  max_steps = 1000
integer,parameter   ::  max_iter_root = 100
real(dp),parameter  ::  zero = 0._dp, one = 1._dp, half = one/2._dp
real(dp),parameter  ::  eps_tf_match = 1._dp*epsilon(0._dp)
real(dp),parameter  ::  eps_root = 5._dp*epsilon(0._dp)
logical,parameter   ::  err_global = .false.
! Gauss-Radau spacings:
! These hs(:) values are the Gauss-Radau spacings, scaled to the
! range 0 to 1, for integrating to order 15. hs(0) = 0. always.
! The sum of these H-values should be 3.7(3) = 3.7333333... = 56/15
! (Viete formulas for the polynomial of degree 7 whose roots are
! hs(1:NSTEPS)-values).
! Values are taken from Ref. [2] up to 25 significant digits. Value for h2
! is corrected according to Ref. [1].
! 
real(qp),parameter  ::  hs(0:nsteps) = [0._qp,&
&0.0562625605369221464656522_qp, 0.1802406917368923649875799_qp,&
&0.3526247171131696373739078_qp, 0.5471536263305553830014486_qp,&
&0.7342101772154105315232106_qp, 0.8853209468390957680903598_qp,&
&0.9775206135612875018911745_qp]

! COMMON VARIABLES
logical :: npq, ncl, nes
integer :: iy

public :: XRA15

contains




subroutine XRA15(nclass,neq,eqs,itask,FORCE,ss,mip,xin,vin,h0in,tain,tzin,tstar&
&,tip,yt,istate,idiag,ht)
! Description:
!     Driver subroutine for the XRA15 numerical solver, which uses the numerical
!     scheme by Everhart [1,2].
! 
! INPUTS:
!   neq: number of scalar ODEs to be integrated.
!   ss: local truncation error tolerance
!   nclass: class of the ODE system. Check [1,2] for explanation on "nclass"
!           values.
!   itask: integer flag for root-finding. If itask = 1 the integration stops
!          at the first root found.
!   eqs:   integer flag for the type of regularized equations of motion being
!          integrated. It is used to choose the appropriate Sundman
!          transformation for output at prescribed times.
!          eqs == 1: Cowell (No root-finding is performed)
!          eqs == 2: K-S, eqs == 3: EDromo, eqs == 4: GDromo.
!   tstar: for regularized formulations, prescribed value of physical time at
!          which the integration is to be stopped through event location.
!   
! LOCALS:
!   ncl: flag which says if the equations are of 1st order (.TRUE.) or
!        of 2nd order (.FALSE.) :
!        y' = F(t,y),    NCL == .TRUE.
!        y" = F(t,y),    NCL == .FALSE.
!        y" = F(t,y,y'), NCL == .FALSE.
!   npq: flag which says if the equations are of 2nd order general
!        (.FALSE.) or NOT 2nd order general (.TRUE.), i.e. of 1st order or
!        2nd order Special (without y'):
!        NCLASS == 1,  NPQ == .TRUE.
!        NCLASS == −2, NPQ == .TRUE.
!        NCLASS == 2,  NPQ == .FALSE.
!   nes: is .TRUE. only if ss > 0. Then the sequence size is fixed to h0.
! 
! OUTPUTS:
! yt: array of size [1:steps,1:(1+NV)] containing the propagated trajectory.
!     For each row i, yt(i,1) is the value of the independent variable; the
!     remaining NV elements are the values of the state vector. NV = neq for
!     Type I, NV = 2*neq for Type II.
! ht: optional array of size [1:steps,1:2] containing the size of the past
!     step for each value of the independent variable.
! istate: integer flag for integration diagnostics. Negative values indicate
!         that the integration was stopped due to an exception. The meaning is
!         the following:
!         NORMAL OUTPUT
!         1: The integration was successfully concluded in the normal mode
!            (itask = 1).
!         2: The integration was successfully concluded in the interpolation
!            mode (itask = 2).
!         3: The integration was successfully concluded in the root-finding 
!            mode, i.e. a root was found. (itask = 3).
!         EXCEPTIONS:
!         -1: No initial step size was found (ncount > max_ncount).
!         -2: Too many integration or interpolation steps (iy = max_steps).
!         -3: NaNs detected in the state vector.
!         -4: A change in sign of g was detected, but the iterative methods
!             did not converge.
!         -5: Failure in allocation/deallocation of variables.
! idiag: integer array containing diagnostic information.
!   idiag(1): Total number of calls to FORCE, idiag(1) = idiag(2) + idiag(3).
!   idiag(2): Evaluations due to normal integration steps performed by XRA15_STEP
!             or outside when the value of f0 is obtained.
!   idiag(3): Evaluations due to event location.
!   idiag(4): Total number of integration (not interpolation) steps taken.
!   idiag(5): If = 1, BISECTION was used to find the initial value for Newton's
!             method.
!   idiag(6): Total number of root-finding steps.
!   idiag(7): If = 1, the B-value iterations did not converge at some time
!             during the integration.
!   idiag(8): Maximum number of B-value iterations during the integration.
! 
! ==============================================================================

! INPUTS
integer,intent(in)  ::  neq,nclass,eqs,mip,itask
real(dp),intent(in) ::  ss,xin(1:neq),vin(1:neq),tip(1:mip),h0in,tain,tzin,tstar
procedure()         ::  FORCE
! OUTPUTS
integer,intent(out) :: istate,idiag(1:10)
real(dp),allocatable,intent(out)  ::  yt(:,:)
real(dp),allocatable,optional,intent(out)  ::  ht(:,:)

! LOCALS
real(dp) :: b(1:nsteps,1:neq), g(1:nsteps,1:neq), e(1:nsteps,1:neq)
real(dp),allocatable  ::  auxy(:,:), auxh(:,:)
real(dp) ::  x(1:neq),v(1:neq),h,hp,t,ta,tz,tf
real(dp) ::  tprev,hprev,xprev(1:neq),vprev(1:neq)
real(dp) ::  f0(1:neq)
real(dp) ::  f0prev(1:neq),bprev(1:nsteps,1:neq),eprev(1:nsteps,1:neq)
real(dp) ::  wc(1:nsteps),uc(1:nsteps),wc0,c(1:ncoef),d(1:ncoef),r(1:ncoef)
real(dp) ::  bd(1:nsteps)
logical  ::  nsf,nper,itrp,evloc
integer  ::  nf,ni,ns,ncount,jip
! Root-finding
real(dp) ::  gprev,gcur,dg2prev,dg2cur
real(dp) ::  troot,xroot(1:neq),vroot(1:neq),groot,dgroot
integer  ::  nin,nfn
integer  ::  bisec
logical  ::  rconv

! PARAMETERS
integer :: ierr

! DIAGNOSTICS
logical :: exc

! ==============================================================================


! ==============================================================================
! 01. GENERAL INITIALIZATIONS.
! ==============================================================================

! INPUT DATA INITIALIZATION.
ncl = (nclass == 1)
npq = (nclass < 2)
nes = (ss >= ONE)

ierr = 0
! WORK SPACE ALLOCATION.
if (ncl) then
  allocate(auxy(1:max_steps,1:1+neq),stat=ierr); if (EXCEPTION(ierr,1,istate)) return
else
  allocate(auxy(1:max_steps,1:1+2*neq),stat=ierr); if (EXCEPTION(ierr,1,istate)) return
end if
if (present(ht)) then
  allocate(auxh(1:max_steps,1:2),stat=ierr)
  if (EXCEPTION(ierr,1,istate)) return
end if

! COEFFICIENTS SETUP.
call RA15_COEFS(ncl,wc,uc,wc0,c,d,r)

! NSF, NPER setup:
! NSF is .FALSE. on starting sequence, otherwise .TRUE.
! NPER is .TRUE. only on the last sequence of integration.
nsf = .false.
nper = .false.
itrp  = (itask == 2)
evloc = (itask == 3)

! Copy inputs into x,v,h
x = xin
v = vin
h = h0in
! If mip > 0, we are in INTERPOLATION MODE. Take first and last times from the
! tip vector.
if (itrp) then
  ta = tip(1)
  tz = tip(mip)
else
  ta = tain
  tz = tzin
end if
t = ta

! Set workspace to zero.
if (ncl) v(:) = zero
auxy = zero
if (present(ht)) auxh = zero
if (nclass == 1) then
  auxy(1,1:1+neq) = [ta,x]
else if (nclass == 2 .or. nclass == -2) then
  auxy(1,1:1+2*neq) = [ta,x,v]
end if
iy = 2

! INITIAL STEP SIZE.
tf = tz - ta
h  = sign(abs(h),tf)
! Set in an estimate to HP based on experience. Same sign as TF. Save into auxh.
hp = sign(0.1_dp,tf)
if (h /= zero) hp = h
if (hp/tf > half) hp = half*tf
if(present(ht)) auxh(1,:) = [ta,0._dp]

! Initialize counters and diagnostics: steps taken, number of B-iterations, 
! tip index, ODE field evaluations, diagnostic array.
ns  = 0
!ni  = 6
jip = 2
nf  = 0
idiag = 0
rconv = .false.

! Compute f0 at the first step
if (npq) then
  call FORCE(neq,t,x,f0)
else
  call FORCE(neq,t,x,v,f0)
end if
idiag(2) = idiag(2) + 1

! NCOUNT is the number of attempts to find the optimal sequence size.
! If NCOUNT > MAX_NCOUNT XRA15 returns to the caller.
ncount = 0

hprev = h; h = hp
b = 0._dp; bprev = 0._dp
e = 0._dp; eprev = 0._dp
bd = 0._dp

istate = 0

! ==============================================================================
! 02. MAIN LOOP
! ==============================================================================

main_loop: do
  ! Check on maximum number of steps
  if (iy >= max_steps) then
    istate = -2
    call CLEAN()
    return
  end if
  
  step_control: do
    ! B-VALUES PREDICTION
    bprev = b; eprev = e
    if (ns /= 0) call BPREDICT(neq,hp,h,ns,b,e)
    
    ! Save values of x, v, t for interpolation and root-finding.
    tprev = t; xprev = x; vprev = v; h = hp
    ! RA15 STEP
    call XRA15_STEP(npq,ncl,FORCE,neq,x,v,t,h,f0,nf,ni,b,d,wc,wc0,uc,r,c,idiag(7))
    if (ni > idiag(8)) idiag(8) = ni
    idiag(2) = idiag(2) + nf
    ns = ns + 1
    ! Check on NaNs in the state vector
    if (any(x /= x) .or. any(v /= v)) then
      istate = -3
      ! Quit without deallocating the outputs (save the trajectory up to t)
      exit main_loop
    end if
    
    ! STEP CONTROL
    hprev = h
    call XRA15_SCONTROL(neq,ns,nes,nsf,nper,ncount,ss,(t-ta),tf,b,f0,wc,h,hp,exc,istate)
    
    if (istate < 0) then
      ! Problems with finding the stepsize. Is there a problem in FORCE?
      call CLEAN()
      return
    end if
    
    if (exc .or. nsf .or. nper) then
      exit step_control
      
    else if (.not.(nsf)) then
      ! If the step must be repeated, reset.
      x = xprev; v = vprev; t = tprev; b = bprev
      h = hp
      ns = ns - 1
    
    end if
        
  end do step_control
  if (exc) exit
  
  ! Compute f0 for the next step.
  f0prev = f0
  if (npq) then
    call FORCE(neq,t,x,f0)
  else
    call FORCE(neq,t,x,v,f0)
  end if
  idiag(2) = idiag(2) + 1
  
  if (itrp) then
    ! Interpolate to output station. TODO: Check b-values.
    call XRA15_INTRP(neq,ncl,jip,mip,tip,tprev,xprev,vprev,h,f0prev,bprev,&
    &wc0,wc,uc,auxy,auxh)
    
    ! Event location is not active in interpolation mode.
    
  else
    ! START OF EVENT LOCATION
    ! We only search for a root in non-interpolation mode. Check for a sign
    ! change of TIMETOUT in this step.
    if (evloc) then
      gprev = TIMETOUT(eqs,neq,tstar,tprev,xprev)
      gcur  = TIMETOUT(eqs,neq,tstar,t,x)
      if (gcur*gprev <= eps_root) then
        ! If sign change, start root-finding. If any of gcur <= eps_root,
        ! gprev <= eps_root, output root (we were very lucky).
        if (abs(gprev) <= eps_root) then
          troot = tprev; xroot = xprev; vroot = vprev
        else if (abs(gcur) <= eps_root) then
          troot = t; xroot = x; vroot = v
        else
          ! Compute second derivatives at the extrema.
          dg2prev = DDSUND(eqs,neq,tprev,xprev,f0prev)
          dg2cur  = DDSUND(eqs,neq,t,x,f0)
          bisec = 0
          
          ! Assign initial value for Newton's method.
          ! Cases: 1) Inflection point at the extrema 2) Monotonous first
          ! derivative 3) Nonmonotonous first derivative
          if (abs(dg2prev) <= eps_root .or. abs(dg2cur) <= eps_root) then
            if (dg2cur > eps_root) then
              troot = t; groot = gcur; xroot = x; vroot = v
            else
              troot = tprev; groot = gprev; xroot = xprev; vroot = vprev
            end if
            
          else if (dg2prev*dg2cur > eps_root) then
            if (dg2prev > eps_root) then
              troot = t; groot = gcur; xroot = x; vroot = v
            else
              troot = tprev; groot = gprev; xroot = xprev; vroot = vprev
            end if
            
          else
            ! Start bisection method to find the initial value for Newton's method.
            idiag(5) = 1
            call BISECTION(eqs,neq,ns,ni,idiag(7),FORCE,tstar,tprev,xprev,vprev,&
            &gprev,dg2prev,hprev,f0prev,bprev,eprev,t,x,v,gcur,dg2cur,nfn,nin,&
            &bisec,troot,groot,xroot,vroot,d,wc,wc0,uc,r,c)
            if (bisec /= 1) then
              ! Bisection did not converge. Exit while saving the trajectory.
              istate = -4
              exit
            end if
            if (ni > idiag(8)) idiag(8) = ni
            idiag(3) = idiag(3) + nfn
            idiag(6) = idiag(6) + nin
            
          end if
          if (bisec /= 2) then
            ! NEWTON'S METHOD.
            dgroot = DSUND(eqs,neq,troot,xroot,vroot)
            call NEWTON(npq,ncl,neq,eqs,FORCE,tstar,troot,groot,dgroot,xroot,vroot,&
            &rconv,nin,nfn,tprev,xprev,vprev,hprev,f0prev,bprev,eprev,ns,ni,&
            &idiag(7),d,wc,wc0,uc,r,c)
            if (ni > idiag(8)) idiag(8) = ni
            idiag(3) = idiag(3) + nfn
            idiag(6) = idiag(6) + nin
          end if
          
          if (rconv .or. bisec == 2) then
            t = troot; x = xroot; v = vroot
            nper = .true.
          else
            ! Newton's method did not converge. Exit while saving the trajectory
            ! up to t.
            istate = -4
            exit
          end if
  
        end if
        
      end if
      ! END OF EVENT LOCATION
    end if
    
    ! Save data to auxy
    if (ncl) then
      auxy(iy,1:1+neq) = [t,x]
    else
      auxy(iy,1:1+2*neq) = [t,x,v]
    end if
    if (present(ht)) auxh(iy,1:2) = [t,h]
    iy = iy + 1
    
  end if
  
  ! EXIT CONDITION
  if (nper) then
    if (istate /= -4) istate = 1
    if (itrp) then
      ! In interpolation mode, save data of the last step.
      if (ncl) then
        auxy(iy,1:1+neq) = [t,x]
      else
        auxy(iy,1:1+2*neq) = [t,x,v]
      end if
      if (present(ht)) auxh(iy,1:2) = [t,h]
      iy = iy + 1
      istate = 2
    end if
    if (rconv) istate = 3
    exit
        
  end if
  
  ! ADJUSTMENT OF LAST STEP
  if (abs((t-ta) + hp) >= abs(tf) - eps_tf_match) then
    ! Check on the final time: if the last step overshoots the final time, adjust
    ! it to cover exactly the integration span.
    hp = tf-(t-ta)
    nper = .true.
    
  end if
  
end do main_loop

! == POST-PROCESSING ==
idiag(1) = idiag(2) + idiag(3)
idiag(4) = ns
yt = auxy(1:iy-1,:)
if (present(ht)) ht = auxh(1:iy-1,:)

! WORK SPACE DEALLOCATION.
if (allocated(auxh)) deallocate(auxh,stat=ierr); if (EXCEPTION(ierr,2,istate)) return
if (allocated(auxy)) deallocate(auxy,stat=ierr); if (EXCEPTION(ierr,2,istate)) return

contains


  function EXCEPTION(ierr,msg,istate)
  integer,intent(in)    ::  ierr,msg
  integer,intent(inout) ::  istate
  logical  ::  EXCEPTION
  
  if (ierr /= 0) then
    select case (msg)
      case(1)
        write(*,*) 'XRA15: Workspace allocation has failed with ierr = ',ierr,'.'
      case(2)
        write(*,*) 'XRA15: Workspace deallocation has failed with ierr = ',ierr,'.'
    end select
    write(*,*) 'XRA15: Returning control to the caller.'
    EXCEPTION = .true.
    istate = -5
  end if
  
  end function EXCEPTION
  
  subroutine CLEAN()
    if (allocated(auxy)) deallocate(auxy)
    if (allocated(auxh)) deallocate(auxh)
  
  end subroutine CLEAN

end subroutine XRA15




subroutine XRA15_STEP(npq,ncl,FORCE,neq,x,v,t,h,f0,nf,ni,b,d,wc,wc0,uc,r,c,warn)
! Description:
!    Performs one integration step with Everhart's Gauss-Radau numerical scheme of
!    order 15. This version implements the predictor-corrector iteration criteria
!    as in [4].
!
! ==============================================================================

! INPUTS
integer,intent(in)   ::  neq
real(dp),intent(in)  ::  f0(1:neq)
real(dp),intent(in)  ::  wc0,wc(1:nsteps),uc(1:nsteps)
real(dp),intent(in)  ::  d(1:ncoef),c(1:ncoef),r(1:ncoef)
logical,intent(in)   ::  npq,ncl
real(dp),intent(in)  ::  h
procedure()          ::  FORCE

! OUTPUTS
integer,intent(out)  ::  nf,warn,ni

! INPUTS/OUTPUTS
real(dp),intent(inout) ::  x(1:neq),v(1:neq),t
real(dp),intent(inout) ::  b(1:nsteps,1:neq)

! LOCALS AND PARAMETERS
real(dp)  ::  h2,s,q
real(dp)  ::  fj(1:neq),tempv(1:neq),y(1:neq),yp(1:neq),gk(1:neq)
real(dp)  ::  g(1:nsteps,1:neq)
integer   ::  j
real(dp)  ::  pc_err,pc_last


! ==============================================================================


nf = 0
warn = 0
pc_err = huge(0._dp)
pc_last = huge(0._dp)
ni = 0

! Find new G-values from the predicted B-values, following Eqs. (7) in [1].
g(1,:) = b(1,:)+D(1)*b(2,:)+D(2)*b(3,:)+D(4)*b(4,:)+D(7)*b(5,:)+D(11)*b(6,:)+&
        &D(16)*b(7,:)
g(2,:) = b(2,:)+D(3)*b(3,:)+D(5)*b(4,:)+D(8)*b(5,:)+D(12)*b(6,:)+D(17)*b(7,:)
g(3,:) = b(3,:)+D(6)*b(4,:)+D(9)*b(5,:)+D(13)*b(6,:)+D(18)*b(7,:)
g(4,:) = b(4,:)+D(10)*b(5,:)+D(14)*b(6,:)+D(19)*b(7,:)
g(5,:) = b(5,:)+D(15)*b(6,:)+D(20)*b(7,:)
g(6,:) = b(6,:)+D(21)*b(7,:)
g(7,:) = b(7,:)

! H    is the sequence size
! HP   is the guessed sequence size
! HVAL is the absolute value of sequence size
! TM   is the current time relative to TA
! TMF  is the current time (time to be passed to the FORCE/FCN)
h2 = h*h
if (ncl) h2 = h

! better_B_loop has 6 iterations on first sequence and 2 iterations therafter
better_B_loop: do
  
  ! EXIT CONDITIONS
  if (pc_err <= 1.E-16) then
    exit
    
  else if (ni > 2 .and. pc_last <= pc_err) then
!    write(*,*) 'Predictor-corrector oscillating.'
    warn = 1
    exit
    
  else if (ni >= 12) then
    warn = 1
    exit
    
  end if
  
  ! This loop is for each substep within a sequence.
  substep_loop: do j = 1,nsteps
    s = hs(j)
    q = s
    if (ncl) q = one
    
    ! ==== PREDICTOR FORMULAS ====
    ! Here Y is used for the value of y at substep n. We use Eq. (9) from [1].
    tempv = wc(3)*b(3,:)+s*(wc(4)*b(4,:)+s*(wc(5)*b(5,:) &
              &+s*(wc(6)*b(6,:)+s*wc(7)*b(7,:))))
    y(:) = x(:)+q*(h*v(:)+h2*s*(f0(:)*wc0+s*(wc(1)*b(1,:) & 
              &+s*(wc(2)*b(2,:)+s*tempv))))
    
    ! If the equations are 1st order or 2nd order-special, we don't need the
    ! velocity predictor.
    
    if (.not.(npq)) then
      ! The velocity predictors are calculated next, if needed for
      ! general Class II. Here YP is used as the value of y' at
      ! substep n (Eq. (10)).
      tempv = uc(3)*b(3,:)+s*(uc(4)*b(4,:)+s*(uc(5)*b(5,:) &
                &+s*(uc(6)*b(6,:)+s*uc(7)*b(7,:))))
      yp(:) = v(:)+s*h*(f0(:)+s*(uc(1)*b(1,:) &
                &+s*(uc(2)*b(2,:)+s*tempv)))
    end if
    
    ! ==== FORCES ====
    if (npq) then
      call FORCE(neq,t+s*h,y,fj)
    else
      call FORCE(neq,t+s*h,y,yp,fj)
    end if
    nf = nf + 1

    ! ==== UPDATE OF B-COEFFICIENTS ====
    ! 
    ! (A)
    ! Find G-values from the force FJ found at current substep.
    ! This section uses Eqs. (4) of [1].
    ! Before save in TEMP the current value.
    ! 
    ! (B)
    ! TEMP is now the improvement on G(J,K) over its former
    ! value. Now we upgrade the B-value using this difference
    ! in the one term.
    ! This section is based on Eqs. (5) of [1].
    ! 
    gk = (fj-f0)/s
    select case (j)
    case (1)
      ! See comment (A) above...
      tempv = g(1,:)
      g(1,:) = gk
      ! See comment (B) above...
      tempv = g(1,:) - tempv
      b(1,:) = b(1,:) + tempv
   
   case (2)
      ! See comment (A) above...
      tempv = g(2,:)
      g(2,:) = (gk-g(1,:))*r(1)
      
      ! See comment (B) above...
      tempv = g(2,:)-tempv
      b(1,:) = b(1,:)+c(1)*tempv
      b(2,:) = b(2,:)+tempv
   
   case(3)
      ! See comment (A) above...
      tempv = g(3,:)
      g(3,:) = ((gk-g(1,:))*r(2)-g(2,:))*r(3)
      
      ! See comment (B) above...
      tempv = g(3,:)-tempv
      b(1,:) = b(1,:)+c(2)*tempv
      b(2,:) = b(2,:)+c(3)*tempv
      b(3,:) = b(3,:)+tempv
   
   case(4)
      ! See comment (A) above...
      tempv = g(4,:)
      g(4,:) = (((gk-g(1,:))*r(4)-g(2,:))*r(5)-g(3,:))*r(6)
      
      ! See comment (B) above...
      tempv = g(4,:)-tempv
      b(1,:) = b(1,:)+c(4)*tempv
      b(2,:) = b(2,:)+c(5)*tempv
      b(3,:) = b(3,:)+c(6)*tempv
      b(4,:) = b(4,:)+tempv
   
   case(5)
      ! See comment (A) above...
      tempv = g(5,:)
      g(5,:) = ((((gk-g(1,:))*r(7)-g(2,:))*r(8)-g(3,:))*r(9)-g(4,:))*r(10)
      
      ! See comment (B) above...
      tempv = g(5,:)-tempv
      b(1,:) = b(1,:)+c(7)*tempv
      b(2,:) = b(2,:)+c(8)*tempv
      b(3,:) = b(3,:)+c(9)*tempv
      b(4,:) = b(4,:)+c(10)*tempv
      b(5,:) = b(5,:)+tempv
    
    case(6)
      ! See comment (A) above...
      tempv = g(6,:)
      g(6,:) = (((((gk-g(1,:))*r(11)-g(2,:))*r(12)-g(3,:))*r(13)-g(4,:))*&
      &r(14)-g(5,:))*r(15)
      
      ! See comment (B) above...
      tempv = g(6,:)-tempv
      b(1,:) = b(1,:)+c(11)*tempv
      b(2,:) = b(2,:)+c(12)*tempv
      b(3,:) = b(3,:)+c(13)*tempv
      b(4,:) = b(4,:)+c(14)*tempv
      b(5,:) = b(5,:)+c(15)*tempv
      b(6,:) = b(6,:)+tempv
    
    case(7)
      ! See comment (A) above...
      tempv = g(7,:)
      g(7,:) = ((((((gk-g(1,:))*r(16)-g(2,:))*r(17)-g(3,:))*r(18)-&
      &g(4,:))*r(19)-g(5,:))*r(20)-g(6,:))*R(21)
      
      ! See comment (B) above...
      tempv = g(7,:)-tempv
      b(1,:) = b(1,:)+c(16)*tempv
      b(2,:) = b(2,:)+c(17)*tempv
      b(3,:) = b(3,:)+c(18)*tempv
      b(4,:) = b(4,:)+c(19)*tempv
      b(5,:) = b(5,:)+c(20)*tempv
      b(6,:) = b(6,:)+c(21)*tempv
      b(7,:) = b(7,:)+tempv
      
      ! Estimate predictor corrector error
      pc_last = pc_err
      if (err_global) then
        pc_err = maxval(abs(tempv))/maxval(abs(f0))
      else
        pc_err = maxval(abs(tempv)/abs(f0))
      end if
      
    end select
  end do substep_loop
  
  ni = ni + 1
  
end do better_B_loop
  
! === CORRECTOR FORMULAS ===
x(:) = x(:) + v(:)*h + h2*(f0(:)*wc0+b(1,:)*wc(1)+b(2,:)*wc(2)+b(3,:)*wc(3)+&
&b(4,:)*wc(4)+b(5,:)*wc(5)+b(6,:)*wc(6)+b(7,:)*wc(7))
if (.not.(ncl)) then
  v(:) = v(:)+h*(f0(:)+b(1,:)*uc(1)+b(2,:)*uc(2)+b(3,:)*uc(3)+b(4,:)*uc(4)+&
  &b(5,:)*uc(5)+b(6,:)*uc(6)+b(7,:)*uc(7))
end if
! We have done a sequence and can update current time and sequence counter.
t = t + h

end subroutine XRA15_STEP




subroutine XRA15_INTRP(neq,ncl,jip,mip,tip,tprev,xprev,vprev,h,f0,b,wc0,wc,uc,&
&auxy,auxh)
! Description:
!    Interpolates the numerical solution inside the last step to the prescribed
!    times vector "tip".
! 
! ==============================================================================

! INPUTS
logical,intent(in)   ::  ncl
integer,intent(in)   ::  mip,neq
real(dp),intent(in)  ::  tprev,h,wc0
real(dp),intent(in)  ::  xprev(1:neq),vprev(1:neq),f0(1:neq),tip(1:mip)
real(dp),intent(in)  ::  b(1:nsteps,1:neq),wc(1:nsteps),uc(1:nsteps)

! INPUTS/OUTPUTS
integer,intent(inout)   ::  jip
real(dp),allocatable,intent(inout)  ::  auxy(:,:),auxh(:,:)

! LOCALS
real(dp)  ::  hout,hout2,sout,tcur
real(dp)  ::  xout(1:neq),vout(1:neq)


! ==============================================================================

tcur = tprev + h
do
  ! Check for values of tip included in the last step.
  ! We also sanity-check if jip > mip.
  if (jip > mip) exit
  ! Remember that tip(mip) == tf. If tip(jip) == tcur at step i, we compute the
  ! values of x, v at step i from the B-values at step (i + 1). This is
  ! mathematically irrelevant since at the beginning of step (i+ 1) hout = hout2 =
  ! sout = 0 anyway, but it is important for the implementation. In fact we
  ! always have tip(mip) == tf, i.e. the final values of x,v are never computed
  ! here, but in the exit "if" block in the XRA15 driver.
  if ((abs(tip(jip)/tcur) - one) >= zero) exit
  hout = tip(jip) - tprev
  hout2 = hout*hout
  sout = hout/h
  
  ! These are just the predictor/corrector equations. The sum was expanded and
  ! rewritten in a more human-readable way.
  if (ncl) then
    xout = xprev + hout*(f0 + b(1,:)*sout*wc(1) + b(2,:)*sout**2*wc(2) + &
    & b(3,:)*sout**3*wc(3) + b(4,:)*sout**4*wc(4) + b(5,:)*sout**5*wc(5) + &
    & b(6,:)*sout**6*wc(6) + b(7,:)*sout**7*wc(7))
    
    ! If the equations are Type 1, there is no velocity predictor.
    vout = 0._dp
    
  else
    xout = xprev + vprev*hout + hout2*(f0*wc0 + b(1,:)*sout*wc(1) + &
    & b(2,:)*sout**2*wc(2) + b(3,:)*sout**3*wc(3) + b(4,:)*sout**4*wc(4)&
    & + b(5,:)*sout**5*wc(5) + b(6,:)*sout**6*wc(6) + b(7,:)*sout**7*wc(7))
    
    vout = vprev + hout*(f0 + b(1,:)*sout*uc(1) + b(2,:)*sout**2*uc(2) +&
    & b(3,:)*sout**3*uc(3) + b(4,:)*sout**4*uc(4) + b(5,:)*sout**5*uc(5) &
    & + b(6,:)*sout**6*uc(6) + b(7,:)*sout**7*uc(7))
    
  end if
  
  ! Save into auxy
  if (ncl) then
    auxy(iy,1:1+neq) = [tip(jip),xout]
  else
    auxy(iy,1:1+2*neq) = [tip(jip),xout,vout]
  end if
  if (allocated(auxh)) auxh(iy,1:2) = [tip(jip),h]
  iy = iy + 1
  
  jip = jip + 1
  
end do


end subroutine XRA15_INTRP




subroutine XRA15_SCONTROL(neq,ns,nes,accept,nper,ncount,ss,tm,tf,b,f,wc,h,hp,exc,istate)
! Description:
!    Computes the step size for the next step, "hp". If "hp" is "sfact" times
!    smaller than the current step size "h", signals to the caller to reject the
!    step. If "hp" is smaller than "h" and this is the first step, it does the same.
!    Uses the logic for step size selection from the IAS15 integrator [4].
! 
! ==============================================================================

! INPUTS
logical,intent(in)      ::  nes
integer,intent(in)      ::  neq,ns
real(dp),intent(in)     ::  wc(1:nsteps)
real(dp),intent(in)     ::  ss,tm,tf
real(dp),intent(in)     ::  b(1:nsteps,1:neq)
real(dp),intent(in)     ::  f(1:neq)
! OUTPUTS
real(dp),intent(out)    ::  hp
logical,intent(out)     ::  exc

! INPUT/OUTPUTS
integer,intent(inout)   ::  istate

! INPUTS/OUTPUTS
logical,intent(inout)   ::  accept,nper
real(dp),intent(inout)  ::  h
integer,intent(inout)   ::  ncount

! LOCALS AND PARAMETERS
integer,parameter  :: max_ncount = 10
real(dp),parameter :: pw = one/7._dp, sr = 0.25_dp
real(dp) :: b7max,fmax,err_est

! ==============================================================================


exc = .false.

if (nes) then
  ! If the step size is fixed, just set accept = .true. if this is the first
  ! sequence.
  if (.not.(accept)) accept = .true.
  hp = h
  
else
  ! Predict required step size using global error estimate, Eq. (9) of [4]
  if (err_global) then
    b7max = maxval(abs(b(7,:)),1)
    fmax = maxval(abs(f),1)
    err_est = b7max/fmax
  else
    err_est = maxval(abs(b(7,:))/abs(f), 1)
    
  end if
  hp = ((ss/err_est)**pw)*h
  
  if (abs(hp) < 2._dp*eps_tf_match) hp = sign(2._dp*eps_tf_match,hp)
  
  accept = .true.
  if (hp/h < sr) then
    ! Reject last step and retry with hp if hp/h < 0.25.
    accept = .false.
    
  else if (hp/h > 1._dp/sr) then
    ! Limit the growth ratio between subsequent steps to 1/0.25 = 4.
    hp = h/sr
    
  end if
  
end if

end subroutine XRA15_SCONTROL




subroutine BPREDICT(neq,hp,h,ns,b,e)
! Description:
!    Applies the correction BD for the B-values as specified in Sec. 2.5 of [1].
! 
! ==============================================================================

! INPUTS
integer,intent(in)   ::  neq,ns
real(dp),intent(in)  ::  hp,h

! INPUTS/OUTPUTS
real(dp),intent(inout)  ::  e(1:nsteps,1:neq)
real(dp),intent(inout)  ::  b(1:nsteps,1:neq)

! LOCALS AND PARAMETERS
real(dp),parameter  ::  Z2 = 2, Z3 = 3, Z4 = 4, Z5 = 5, Z6 = 6, Z7 = 7, &
                        Z10 = 10, Z15 = 15, Z20 = 20, Z21 = 21, Z35 = 35
real(dp)  ::  q,q2,q3,q4,q5,q6,q7
real(dp)  ::  bd(1:nsteps,1:neq)
integer   ::  j

! ==============================================================================

q = hp/h
q2 = q*q
q3 = q*q2
q4 = q2*q2
q5 = q2*q3
q6 = q3*q3
q7 = q3*q4

bd = 0._dp
if (ns /= 1) then
  do j=1,nsteps
    bd(j,:) = b(j,:) - e(j,:)
  
  end do
  
end if

e(1,:) = q*(b(1,:)+Z2*b(2,:)+Z3*b(3,:)+Z4*b(4,:)+Z5*b(5,:)+Z6*b(6,:)+Z7*b(7,:))
e(2,:) = q2*(b(2,:)+Z3*b(3,:)+Z6*b(4,:)+Z10*b(5,:)+Z15*b(6,:)+Z21*b(7,:))
e(3,:) = q3*(b(3,:)+Z4*b(4,:)+Z10*b(5,:)+Z20*b(6,:)+Z35*b(7,:))
e(4,:) = q4*(b(4,:)+Z5*b(5,:)+Z15*b(6,:)+Z35*b(7,:))
e(5,:) = q5*(b(5,:)+Z6*b(6,:)+Z21*b(7,:))
e(6,:) = q6*(b(6,:)+Z7*b(7,:))
e(7,:) = q7*b(7,:)

! Apply the correction. Notice that when we have done ONLY
! one sequence (NS == 1), BD == 0 from its initialization, i.e.
! we are doing B = E. It is only when NS > 1 that we are applying
! the correction BD.
do j = 1,nsteps
   b(j,:) = e(j,:)+bd(j,:)
end do

end subroutine BPREDICT




subroutine RA15_COEFS(ncl,wc,uc,wc0,c,d,r)
! Description:
!    Computes the coefficients for Everhart's numerical scheme of 15th order [1,2].
! 
! ==============================================================================

! INPUTS
logical,intent(in)  ::  ncl
! OUTPUTS
real(dp),intent(out) ::  wc(1:nsteps),uc(1:nsteps),wc0
real(dp),intent(out) ::  c(1:ncoef),d(1:ncoef),r(1:ncoef)
! LOCALS AND PARAMETERS
real(dp) ::  temp
integer  ::  l,k,la,lb,lc,ld,le
integer,parameter :: nw(0:nsteps)= [0,0,1,3,6,10,15,21]


! CONSTANT coefficients setup
k = 2
do  l = 1,nsteps
    temp = k+k*k
    if (ncl) temp = k
    wc(l) = one/temp
    temp = k
    uc(l) = one/temp
    k = k+1
end do

wc0 = half
if (ncl) wc0 = one

c(1) = -hs(1)
d(1) = hs(1)
r(1) = one/(hs(2)-hs(1))
la = 1
lc = 1
do k = 3,nsteps
   lb = la
   la = lc+1
   lc = nw(k)
   c(la) = -hs(k-1)*c(lb)
   c(lc) = c(la-1)-hs(k-1)
   
   d(la) = hs(1)*d(lb)
   d(lc) = -c(lc)
   
   r(la) = one/(hs(k)-hs(1))
   r(lc) = one/(hs(k)-hs(k-1))
   
   if (k == 3) cycle
   
   do l = 4, k
      ld = la+l-3
      le = lb+l-4
      c(ld) = c(le)-hs(k-1)*c(le+1)
      d(ld) = d(le)+hs(l-2)*d(le+1)
      r(ld) = one/(hs(k)-hs(l-2))
   end do
end do

end subroutine RA15_COEFS




subroutine NEWTON(npq,ncl,neq,eqs,FORCE,tstar,sroot,groot,dgroot,xroot,vroot,&
&conv,niter,neval,s0,x0,v0,h0,f0,b0,e0,ns,ni,warn,d,wc,wc0,uc,r,c)
! Description:
!    Finds the value of the root 'sroot' of the scalar function G, G(sroot) = 0,
!    starting from the initial guess 's0'.
!    The user has to provide the function DG evaluating the derivative of G.
!    The input variables with the '0' suffix refer to their value at the beginning
!    of the current integration step.
!
! ==============================================================================

! INPUTS
logical,intent(in)   ::  npq,ncl
integer,intent(in)   ::  neq,eqs,ns
real(dp),intent(in)  ::  tstar
real(dp),intent(in)  ::  s0,x0(1:neq),v0(1:neq)
real(dp),intent(in)  ::  h0,f0(1:neq),b0(1:nsteps,1:neq),e0(1:nsteps,1:neq)
real(dp),intent(in)  ::  d(1:ncoef),r(1:ncoef),c(1:ncoef)
real(dp),intent(in)  ::  wc(1:nsteps),wc0,uc(1:nsteps)
procedure()  ::  FORCE

! INPUTS/OUTPUTS
real(dp),intent(inout)  ::  sroot,groot,dgroot

! OUTPUTS
real(dp),intent(out) :: xroot(1:neq),vroot(1:neq)
logical,intent(out)  :: conv
integer,intent(out)  :: niter,neval,ni,warn

! LOCALS
integer :: nf
real(dp)  ::  b(1:nsteps,1:neq),e(1:nsteps,1:neq)
real(dp)  ::  s
real(qp)  ::  qsroot,qgroot,qdgroot,grootprev

! ==============================================================================

! Initialize number of iterations and evaluations, convergence flag.
niter = 0
nf    = 0
neval = 0
conv = .false.

qsroot = real(sroot,qp); qgroot = real(groot,qp); qdgroot = real(dgroot,qp)

do
  niter = niter + 1
  ! Check on maximum number of iterations
  if (niter > max_iter_root) then
    conv = .false.
    exit
  end if
  
  ! Newton iteration (performed in quad precision).
  qsroot = real(sroot,qp); qgroot = real(groot,qp); qdgroot = real(dgroot,qp)
  qsroot = qsroot - qgroot/qdgroot
  sroot = real(qsroot,dp)
  
  ! COMPUTATION OF G, DG.
  ! Perform XRA15 substeps to get the value of the state vector for each
  ! iteration. The starting point is always at the beginning of the integration
  ! step which was just performed.
  ! Before performing the step, extrapolate the B-values.
  b = b0; e = e0
  call BPREDICT(neq,(sroot-s0),h0,ns-1,b,e)
  xroot = x0; vroot = v0; s = s0;
  call XRA15_STEP(npq,ncl,FORCE,neq,xroot,vroot,s,(sroot - s0)&
  &,f0,nf,ni,b,d,wc,wc0,uc,r,c,warn)
  neval = neval + nf
  
  grootprev = groot
  groot  = TIMETOUT(eqs,neq,tstar,sroot,xroot)
  dgroot = DSUND(eqs,neq,sroot,xroot,vroot)
  
  ! Nominal exit condition: exit if groot underflows or if it oscillates
  ! due to roundoff.
  if (abs(groot) <= eps_root .or. (abs(groot) >= abs(grootprev))) then
    conv = .true.
    exit
  end if

end do

end subroutine NEWTON




subroutine BISECTION(eqs,neq,ns,ni,warn,FORCE,tstar,sprev,xprev,vprev,gprev,&
&dg2prev,hprev,f0prev,bprev,eprev,scur,xcur,vcur,gcur,dg2cur,neval,niter,cflag,&
&s0_guess,g0_guess,x0_guess,v0_guess,d,wc,wc0,uc,r,c)
! Description:
!    Initial guess for Newton's method when dg2prev*dg2cur < 0.
! 
! ==============================================================================

! INPUTS
integer,intent(in)  :: eqs,neq,ns
real(dp),intent(in) :: sprev,hprev,xprev(1:neq),vprev(1:neq)
real(dp),intent(in) :: scur,xcur(1:neq),vcur(1:neq)
real(dp),intent(in) :: dg2prev,dg2cur
real(dp),intent(in) :: f0prev(1:neq),bprev(1:nsteps,1:neq),eprev(1:nsteps,1:neq)
real(dp),intent(in) :: tstar
real(dp),intent(in) :: d(1:ncoef),r(1:ncoef),c(1:ncoef)
real(dp),intent(in) :: wc(1:nsteps),wc0,uc(1:nsteps)
real(dp),intent(in) ::  gprev,gcur
procedure()  ::  FORCE

! OUTPUTS
real(dp),intent(out) ::  s0_guess,g0_guess,x0_guess(1:neq),v0_guess(1:neq)
integer,intent(out)  ::  neval,niter,ni,warn
integer,intent(out)  ::  cflag

! LOCALS
integer   ::  nf
real(dp)  ::  sa,salast,sb,sblast,smid,s
real(dp)  ::  xa(1:neq),va(1:neq),xb(1:neq),vb(1:neq),xmid(1:neq),vmid(1:neq)
real(dp)  ::  fa(1:neq),fb(1:neq)
real(dp)  ::  b(1:nsteps,1:neq),e(1:nsteps,1:neq)
real(dp)  ::  gmid,dg2a,dg2b

! ==============================================================================

cflag = 0
! Set the extrema for the first iteration to be the ones of the
! current integration step.
sa   = sprev;     salast = sa
sb   = scur;      sblast = sb
xa   = xprev;     xb     = xcur
dg2a = dg2prev;   dg2b   = dg2cur
smid = 0.5_dp*(sa + sb)

neval = 0
nf    = 0
niter = 0

do
  niter = niter + 1
  ! Check on max. number of iterations
  if (niter >= 2*max_iter_root) then
    cflag = 0
    return
  end if
  
  ! COMPUTE G IN THE MIDPOINT.
  b = bprev; e = eprev
  call BPREDICT(neq,(smid-sprev),hprev,ns-1,b,e);
  xmid = xprev; vmid = vprev; s = sprev;
  call XRA15_STEP(npq,ncl,FORCE,neq,xmid,vmid,s,(smid-s),f0prev,&
  &nf,ni,b,d,wc,wc0,uc,r,c,warn)
  neval = neval + nf
  
  gmid = TIMETOUT(eqs,neq,tstar,smid,xmid)
  
  if (gmid < -eps_root) then
    ! Bisect to the right.
    sa = smid; xa = xmid; va = vmid
    sb = sblast
    if (npq) then
      call FORCE(neq,sa,xa,fa)
    else
      call FORCE(neq,sa,xa,va,fa)
    end if
    neval = neval + 1
    dg2a = DDSUND(eqs,neq,sa,xa,fa)
    
  else if (gmid > eps_root) then
    ! Bisect to the left.
    sb = smid; xb = xmid; vb = vmid
    sa = salast
    if (npq) then
      call FORCE(neq,sa,xb,fb)
    else
      call FORCE(neq,sa,xb,vb,fb)
    end if
    neval = neval + 1
    dg2b = DDSUND(eqs,neq,sb,xb,fb)
  
  else
    ! We have found a root by bisection
    s0_guess = smid; x0_guess = xmid; v0_guess = vmid
    cflag = 2
    return
  
  end if
  
  ! If the root is close to the inflection point, find it through bisection
  if (dg2a*dg2b >= eps_root .and. (abs(gmid) > 10._dp*eps_root)) then
    cflag = 1
    exit
  end if
  ! Bisect interval
  smid = 0.5_dp*(sa + sb)
  salast = sa; sblast = sb
  
end do

! Assign the initial guess for Newton's method. 'smid' = s_{i,M} in Ref. [3].
if (gmid > eps_root) then
  ! Case t* < tF (Ref. [3])
  if (dg2prev < -eps_root) then
    s0_guess = sprev; g0_guess = gprev; x0_guess = xprev; v0_guess = vprev
  else if (dg2prev >= eps_root) then
    s0_guess = smid; g0_guess = gmid; x0_guess = xmid; v0_guess = vmid
  end if
  
else if (gmid < -eps_root) then
  ! Case t* > tF (Ref. [3])
  if (dg2cur <= -eps_root) then
    s0_guess = smid; g0_guess = gmid; x0_guess = xmid; v0_guess = vmid
  else if (dg2cur > eps_root) then
    s0_guess = scur; g0_guess = gcur; x0_guess = xcur; v0_guess = vcur
  end if
  
end if

end subroutine BISECTION

end module XEVERHART
