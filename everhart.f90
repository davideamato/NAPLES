!
! Fortran Interface to the Everhart Integrator Library
! by Angelo Graziosi   (firstname.lastnameATalice.it)
! Copyright Angelo Graziosi
! 
! It is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! 
! This is the 'everhart' module.
! 
! 
! A simple module which tries to re-implement in modern Fortran the
! Everhart's RADAU integrator.
! 
! Ref. :
! 
! [1] E. Everhart, An Efficient Integrator That Uses Gauss-Radau Spacings,
!     in A. Carusi and G. B. Valsecchi - Dynamics of Comets: Their Origin and
!     Evolution, 185-202. 1985 by D. Reidel Publishing Company.
! [2] E. Everhart, Implicit Single-Sequence Methods For Integrating Orbits,
!     Celestial Mechanics, 10, pp. 35-55, 1974.
!
! 
! Revision: Davide Amato
!           Technical University of Madrid
!           d.amato@upm.es
! 
! The following modifications to the original code by A. Graziosi have been
! introduced:
! Version 1:
!    - Suppression of debug_flag, save_data_flag. The function calls are now
!      always saved, and the data is output through a new allocatable array YT.
!    - Output of diagnostics through arrays idiag, rdiag.
!    - Each call to "force" in RA15 has been modified to accomodate
!      LSODAR-style procedures.
!    - In RADAU_ON, input "ss0" directly instead of passing through "ll0". Also,
!      the fixed stepsize is set for ss0 > 1. If a fixed stepsize is wanted, set
!      ss0 >> 1. + epsilon(0._dp), e.g. ss0 = 2 is fine.
!    - Minor typographical and style changes.
! Version 2:
!    - Low-order interpolation using the predictor equations.
!    - Output of the step size history.
!
! The order of this integrator has been checked on April 26th, 2016. (D. Amato)
! It has been tested on Keplerian orbits on May 4th, 2016. (D. Amato)
!

module EVERHART

use KINDS, only: dp=>dk,qp=>qk
implicit none
private

integer,parameter  ::  NOR = 15
integer,parameter  ::  NSTEPS = 7, NCOEF = (NSTEPS*(NSTEPS-1))/2
integer,parameter  ::  CHUNK_SIZE = 100000
real(DP), parameter :: ZERO = 0._dp, ONE = 1._dp, HALF = ONE/2._dp
!
! These hs(:) values are the Gauss-Radau spacings, scaled to the
! range 0 to 1, for integrating to order 15. hs(0) = 0. always.
! The sum of these H-values should be 3.7(3) = 3.7333333... = 56/15
! (Viete formulas for the polynomial of degree 7 whose roots are
! hs(1:NSTEPS)-values).
! Values are taken from Ref. [2] up to 25 significant digits. Value for h2
! is corrected according to Ref. [1].
! 

! GAUSS-RADAU SPACINGS:
real(qp),parameter  ::  hs(0:NSTEPS) = [0._qp,&
&0.0562625605369221464656522_qp, 0.1802406917368923649875799_qp,&
&0.3526247171131696373739078_qp, 0.5471536263305553830014486_qp,&
&0.7342101772154105315232106_qp, 0.8853209468390957680903598_qp,&
&0.9775206135612875018911745_qp]

integer :: nv, nclass, log_unit, data_unit
logical :: npq, ncl, nes
! WC, UC, WC0, SS, C, D, R are, really, CONSTANTS
real(dp) :: WC(NSTEPS), UC(NSTEPS), WC0, SS
real(dp) :: C(NCOEF), D(NCOEF), R(NCOEF)
!
! The workspace would be NV x 3*NSTEPS+4 = NV x 3*7+4 −−> w(NV,25)
! (BSG uses w(NEQ,36), being NEQ the number of equations of 1st order)
!
real(dp),allocatable :: f0(:), fj(:), y(:), yp(:)
real(dp),allocatable :: b(:,:), g(:,:), e(:,:), AUXY(:,:), AUXH(:,:)

public  ::  radau_on,radau_off,ra15

contains


subroutine radau_on(nv0,ss0,nclass0)
integer,  intent(in) :: nv0, nclass0
real(dp), intent(in) :: ss0
integer, parameter :: NW(0:NSTEPS)= [ 0, 0, 1, 3, 6, 10, 15, 21 ]
real(DP) :: temp
integer :: l, la, lb, lc, ld, le, k, ierr

! Work space allocation
allocate(b(NSTEPS,nv0),stat=ierr)
if (ierr /= 0) then
    write(*,*) '*** FATAL ERROR ***'
    write(*,*) 'Allocation failure for B(:,:). Exiting...'
    stop
end if
allocate(g(NSTEPS,nv0),stat=ierr)
if (ierr /= 0) then
    write(*,*) '*** FATAL ERROR ***'
    write(*,*) 'Allocation failure for G(:,:). Exiting...'
    stop
end if
allocate(e(NSTEPS,nv0),stat=ierr)
if (ierr /= 0) then
    write(*,*) '*** FATAL ERROR ***'
    write(*,*) 'Allocation failure for E(:,:). Exiting...'
    stop
end if
allocate(f0(nv0),stat=ierr)
if (ierr /= 0) then
    write(*,*) '*** FATAL ERROR ***'
    write(*,*) 'Allocation failure for F0(:). Exiting...'
    stop
end if
allocate(fj(nv0),stat=ierr)
if (ierr /= 0) then
    write(*,*) '*** FATAL ERROR ***'
    write(*,*) 'Allocation failure for FJ(:). Exiting...'
    stop
end if
allocate(y(nv0),stat=ierr)
if (ierr /= 0) then
    write(*,*) '*** FATAL ERROR ***'
    write(*,*) 'Allocation failure for Y(:). Exiting...'
    stop
end if
allocate(yp(nv0),stat=ierr)
if (ierr /= 0) then
    write(*,*) '*** FATAL ERROR ***'
    write(*,*) 'Allocation failure for YP(:). Exiting...'
    stop
end if
ierr = 1
if (nclass0 == 1) then
  allocate(AUXY(1:CHUNK_SIZE,1:1+nv0),stat=ierr)
else if (nclass0 == 2 .or. nclass0 == -2) then
  allocate(AUXY(1:CHUNK_SIZE,1:1+2*nv0),stat=ierr)
end if
if (ierr /= 0) then
  write(*,*) '*** FATAL ERROR ***'
  write(*,*) 'Allocation failure for AUXY(:,:). Exiting...'
  stop
end if
allocate(AUXH(1:CHUNK_SIZE,1:2),stat=ierr)
if (ierr /= 0) then
  write(*,*) '*** FATAL ERROR ***'
  write(*,*) 'Allocation failure for AUXH(:,:). Exiting...'
  stop
end if


! Input data initialization
nv = nv0
nclass = nclass0
!
! Global logical data initialization
!
! NCL is a flag which says if the equations are of 1st order (.TRUE.) or
! of 2nd order (.FALSE.) :
!
! y' = F(t,y),    NCL == .TRUE.
! y" = F(t,y),    NCL == .FALSE.
! y" = F(t,y,y'), NCL == .FALSE.
!
! NPQ is a flag which says if the equations are of 2nd order general
! (.FALSE.) or NOT 2nd order general (.TRUE.), i.e. of 1st order or
! 2nd order Special (without y') :
!
! NCLASS == 1,  NPQ == .TRUE.
! NCLASS == −2, NPQ == .TRUE.
! NCLASS == 2,  NPQ == .FALSE.
! NES is .TRUE. only if LL is negative. Then the sequence size is H0.
!
ncl = (nclass == 1)
npq = (nclass < 2)
nes = (ss0 >= ONE)

! CONSTANT coefficients setup
k = 2
do  l = 1, NSTEPS
    temp = k+k*k
    if (ncl) temp = k
    WC(l) = ONE/temp
    temp = k
    UC(l) = ONE/temp
    k = k+1
end do

WC0 = HALF
if (ncl) WC0 = ONE

C(1) = -HS(1)
D(1) = HS(1)
R(1) = ONE/(HS(2)-HS(1))
la = 1
lc = 1
do k = 3, NSTEPS
   lb = la
   la = lc+1
   lc = NW(k)
   C(la) = -HS(k-1)*C(lb)
   C(lc) = C(la-1)-HS(k-1)
   
   D(la) = HS(1)*D(lb)
   D(lc) = -C(lc)
   
   R(la) = ONE/(HS(k)-HS(1))
   R(lc) = ONE/(HS(k)-HS(k-1))
   
   if (k == 3) cycle
   
   do l = 4, k
      ld = la+l-3
      le = lb+l-4
      C(ld) = C(le)-HS(k-1)*C(le+1)
      D(ld) = D(le)+HS(l-2)*D(le+1)
      R(ld) = ONE/(HS(k)-HS(l-2))
   end do
end do

! SS is, really, a CONSTANT (like WC, UC, and WC0)
!SS = 10.0_DP ** (-ll)
SS = SS0
!
! The statements above are used only once in an integration to set up
! the constants. They use less than a second of execution time.
!
end subroutine radau_on


subroutine radau_off()
integer :: ierr

! Freeing work space
if (allocated(AUXH)) deallocate(AUXH,stat=ierr)
if (ierr /= 0) then
   write(*,*) '*** FATAL ERROR ***'
   write(*,*) 'Deallocation failure for AUXH(:,:). Exiting...'
   stop
end if
if (allocated(AUXY)) deallocate(AUXY,stat=ierr)
if (ierr /= 0) then
   write(*,*) '*** FATAL ERROR ***'
   write(*,*) 'Deallocation failure for AUXY(:,:). Exiting...'
   stop
end if
if (allocated(yp)) deallocate(yp,stat=ierr)
if (ierr /= 0) then
   write(*,*) '*** FATAL ERROR ***'
   write(*,*) 'Deallocation failure for YP(:). Exiting...'
   stop
end if
if (allocated(y)) deallocate(y,stat=ierr)
if (ierr /= 0) then
   write(*,*) '*** FATAL ERROR ***'
   write(*,*) 'Deallocation failure for Y(:). Exiting...'
   stop
end if
if (allocated(fj)) deallocate(fj,stat=ierr)
if (ierr /= 0) then
   write(*,*) '*** FATAL ERROR ***'
   write(*,*) 'Deallocation failure for FJ(:). Exiting...'
   stop
end if
if (allocated(f0)) deallocate(f0,stat=ierr)
if (ierr /= 0) then
   write(*,*) '*** FATAL ERROR ***'
   write(*,*) 'Deallocation failure for F0(:). Exiting...'
   stop
end if
if (allocated(e)) deallocate(e,stat=ierr)
if (ierr /= 0) then
   write(*,*) '*** FATAL ERROR ***'
   write(*,*) 'Deallocation failure for E(:,:). Exiting...'
   stop
end if
if (allocated(g)) deallocate(g,stat=ierr)
if (ierr /= 0) then
   write(*,*) '*** FATAL ERROR ***'
   write(*,*) 'Deallocation failure for G(:,:). Exiting...'
   stop
end if
if (allocated(b)) deallocate(b,stat=ierr)
if (ierr /= 0) then
   write(*,*) '*** FATAL ERROR ***'
   write(*,*) 'Deallocation failure for B(:,:). Exiting...'
   stop
end if
end subroutine radau_off


subroutine ra15(tain,tzin,xin,vin,YT,h0in,mip,tip,force,idiag,rdiag,HT)
integer,  intent(in) :: mip
real(dp), intent(in) :: tain, tzin
real(dp), intent(in) :: xin(1:nv), vin(1:nv), tip(1:mip)
real(dp), intent(in) :: h0in
real(dp), allocatable, intent(out) :: YT(:,:)
integer, intent(out)  ::  idiag(:)
real(dp), intent(out) ::  rdiag(:)
real(dp), allocatable, optional, intent(out)  ::  HT(:,:)
procedure() ::  force
!
! Integrator by E. Everhart, Physics Department, University of Denver
! Revision : Angelo Graziosi (Sept. 12, 2014)
!
! This 15th-order version is called RA15. Order NOR is 15.
!
! y' = F(t,y)    is NCLASS == 1,   y" = F(t,y) is NCLASS == -2,
! y" = F(t,y,y') is NCLASS == 2,
!
! TF is t(final)-t(initial). (Negative when integrating backward.)
! NV = the number of simultaneous differential equations.
! 
! 
! LL controls accuracy. Thus SS = 10.**(-LL) controls the size of
! the last term in a series. Try LL = 8 and work up or down from there.
! However, if LL < 0, then H0 is the constant sequence size used.
! A non-zero H0 sets the size of the first sequence regardless of
! LL sign. Zero's and Oh's could look alike here. Use care!
! 
! Xin and Vin are the starting position and velocity vectors (values of y
! and y' at t = ta).
! 
! YT contains the integrated state vector. The rows of YT contain the time and
! the state vector at the end of each integration sequence (i.e. at each time
! step of the integration).
! 
! h0in is the initial stepsize given by the user.
! 
! tip is a vector of times at which output is requested by interpolation.
! In interpolation mode, the first and last element of tip MUST BE the start and
! end times of integration respectively.
! 
! mip is the size of tip, and controls the INTERPOLATION MODE. If mip > 0, then
! itrp = .true. and the integrator is in interpolation mode: the output is given
! at the times specified by the tip vector.
! 
! idiag is an integer array containing integration diagnostics. The contents
! are the following:
! idiag(1): total number of ODE field evaluations
! idiag(2): total number of steps taken
!
! rdiag is a real array containing integration diagnostics. The contents are the
! following:
! rdiag(1): stepsize used in the last step.
!
! HT is a real array of length NS if not in interpolation mode, MIP otherwise.
! It contains information on the stepping. The contents are the following:
! HT(i,1):    previous internal time step at step "i". If itrp = .FALSE., it
!             coincides with the times vector, YT(:,1).
! HT(i,2):    previous internal step size at step "i". If itrp = .FALSE., then
!             HT(i+1,1) = HT(i,1) + HT(i,2).
!
integer, parameter :: MAX_NCOUNT = 10
real(DP), parameter :: SR = 1.4_DP, PW = ONE/9., EPS_TF_MATCH = 1.0E-08_DP, &
     Z2 = 2, Z3 = 3, Z4 = 4, Z5 = 5, Z6 = 6, Z7 = 7, &
     Z10 = 10, Z15 = 15, Z20 = 20, Z21 = 21, Z35 = 35
integer, save :: i, k, ncount, ns, nf, ni, j
logical, save :: nsf, nper
logical  ::  itrp
real(dp), save :: tf, hval, hp, tm, tmf, h, h2, s, q, temp, hv, &
     bd(NSTEPS), q2, q3, q4, q5, q6, q7
real(dp) :: x(1:nv), v(1:nv), xprev(1:nv), vprev(1:nv), xout(1:nv), vout(1:nv)
real(dp) :: tprev, hout, hout2, h0, ta, tz, sout, qout, tempv(1:nv)
integer  :: jip


! 
! NSF is .FALSE. on starting sequence, otherwise .TRUE.
! NPER is .TRUE. only on last sequence of integration.
! 
nsf = .false.
nper = .false.
itrp = (mip > 0)

! Copy inputs into x,v,h0
x = xin
v = vin
h0 = h0in
! If mip > 0, we are in INTERPOLATION MODE. Take first and last times from the
! tip vector.
if (itrp) then
  ta = tip(1)
  tz = tip(mip)
else
  ta = tain
  tz = tzin
end if

! Initialize the working space. We need to initialize B, BD and the auxiliary
! array AUXY. 
if (ncl) v(:) = ZERO
b(:,:) = ZERO
bd(:) = ZERO
AUXY = ZERO
AUXH = ZERO
if (nclass == 1) then
  AUXY(1,1:1+nv) = [ta,x]
else if (nclass == 2 .or. nclass == -2) then
  AUXY(1,1:1+2*nv) = [ta,x,v]
end if


tf = tz - ta
h0 = sign(abs(h0),tf)
!
! Now set in an estimate to HP based on experience. Same sign as TF.
! 
hp = sign(0.1_DP,tf)
if (h0 /= ZERO) hp = h0
if (hp/tf > HALF) hp = HALF*tf

! Save into AUXH
if (present(HT)) AUXH(1,:) = [ta,hp]

! NCOUNT is the number of attempts to find the optimal sequence size.
! If NCOUNT > MAX_NCOUNT it returns to the caller: integration failed.
ncount = 0

! 
! Now the loop regarding the first sequence, aka THE MAIN LOOP, or
! if you prefer, the main sequence loop.
! 
! NS is the number of sequences done
! NF is the number of calls to FORCE subroutine
! NI is the number of iterations to predict the B-values. NI is 6 for
! the first sequence, 2 after it.
!
main_loop: do
   ns = 0
   jip = 2
   nf = 0
   ni = 6
   tm = ZERO
   tmf = ta
!   call force(tmf,x,v,f0)
   select case (nclass)
   case(1,-2)
     call force(nv,tmf,x,f0)
   case(2)
     call force(nv,tmf,x,v,f0)
   end select
   nf = nf + 1
   
   ! Now begins every sequence after the first. First find new
   ! G-values from the predicted B-values, following Eqs. (7) in text.
   every_sequence_loop: do
      do k = 1, nv
         g(1,k) = b(1,k)+D(1)*b(2,k)+D(2)*b(3,k)+D(4)*b(4,k)+D(7)*b(5,k) &
              +D(11)*b(6,k)+D(16)*b(7,k)
         g(2,k) = b(2,k)+D(3)*b(3,k)+D(5)*b(4,k)+D(8)*b(5,k)+D(12)*b(6,k) &
              +D(17)*b(7,k)
         g(3,k) = b(3,k)+D(6)*b(4,k)+D(9)*b(5,k)+D(13)*b(6,k)+D(18)*b(7,k)
         g(4,k) = b(4,k)+D(10)*b(5,k)+D(14)*b(6,k)+D(19)*b(7,k)
         g(5,k) = b(5,k)+D(15)*b(6,k)+D(20)*b(7,k)
         g(6,k) = b(6,k)+D(21)*b(7,k)
         g(7,k) = b(7,k)
      end do
      ! H    is the sequence size
      ! HP   is the guessed sequence size
      ! HVAL is the absolute value of sequence size
      ! TM   is the current time relative to TA
      ! TMF  is the current time (time to be passed to the force/FCN)
      h = hp
      h2 = h*h
      if (ncl) h2 = h
      hval = abs(h)
      
      ! better_B_loop is 6 iterations on first sequence and
      ! 2 iterations therafter
      better_B_loop: do i = 1, ni
         ! This loop is for each substep within a sequence.
         substep_loop: do j = 1, NSTEPS
            s = HS(j)
            q = s
            if (ncl) q = ONE
            ! Here Y is used for the value of y at substep n.
            ! We use Eq. (9). The collapsed series are broken in two parts
            ! because an otherwise excellent compiler could not handle the
            ! complicated expression.
            do k = 1, nv
               temp = WC(3)*b(3,k)+s*(WC(4)*b(4,k)+s*(WC(5)*b(5,k) &
                    +s*(WC(6)*b(6,k)+s*WC(7)*b(7,k))))
               y(k) = x(k)+q*(h*v(k)+h2*s*(f0(k)*WC0+s*(WC(1)*b(1,k) & 
                    +s*(WC(2)*b(2,k)+s*temp))))
               
               ! If equations are 1st order or 2nd order special (i.e.
               ! without y', continue oops.. cycle..
               if (npq) cycle
               
               ! The velocity predictors are calculated next, if needed for
               ! general Class II. Here YP is used as the value of y' at
               ! substep n (Eq. (10)).
               temp = UC(3)*b(3,k)+s*(UC(4)*b(4,k)+s*(UC(5)*b(5,k) &
                    +s*(UC(6)*b(6,k)+s*UC(7)*b(7,k))))
               yp(k) = v(k)+s*h*(f0(k)+s*(UC(1)*b(1,k) &
                    +s*(UC(2)*b(2,k)+s*temp)))
            end do
            
            !!! DEBUG
!            write(*,*) b
!            write(*,*) y
!            write(*,*) yp
            !!! DEBUG
            
            ! Find forces at each substep.
!            call force(tmf+s*h,y,yp,fj)
            select case (nclass)
            case(1,-2)
              call force(nv,tmf+s*h,y,fj)
            case(2)
              call force(nv,tmf+s*h,y,yp,fj)
            end select
            nf = nf + 1
            ! 
            ! (A)
            ! Find G-values from the force FJ found at current substep.
            ! This section uses Eqs. (4) of text.
            ! Before save in TEMP the current value.
            ! 
            ! (B)
            ! TEMP is now the improvement on G(J,K) over its former
            ! value. Now we upgrade the B-value using this difference
            ! in the one term.
            ! This section is based on Eqs. (5).
            ! 
            select case (j)
            case (1)
               do k = 1, nv
                  q = (fj(k)-f0(k))/s
                  
                  ! See comment (A) above...
                  temp = g(1,k)
                  g(1,k) = q
                  
                  ! See comment (B) above...
                  temp = g(1,k)-temp
                  b(1,k) = b(1,k)+temp
               end do
            case (2)
               do k = 1, nv
                  q = (fj(k)-f0(k))/s
                  
                  ! See comment (A) above...
                  temp = g(2,k)
                  g(2,k) = (q-g(1,k))*R(1)
                  
                  ! See comment (B) above...
                  temp = g(2,k)-temp
                  b(1,k) = b(1,k)+C(1)*temp
                  b(2,k) = b(2,k)+temp
               end do
            case (3)
               do k = 1, nv
                  q = (fj(k)-f0(k))/s
                  
                  ! See comment (A) above...
                  temp = g(3,k)
                  g(3,k) = ((q-g(1,k))*R(2)-g(2,k))*R(3)
                  
                  ! See comment (B) above...
                  temp = g(3,k)-temp
                  b(1,k) = b(1,k)+C(2)*temp
                  b(2,k) = b(2,k)+C(3)*temp
                  b(3,k) = b(3,k)+temp
               end do
            case (4)
               do k = 1, nv
                  q = (fj(k)-f0(k))/s
                  
                  ! See comment (A) above...
                  temp = g(4,k)
                  g(4,k) = (((q-g(1,k))*R(4)-g(2,k))*R(5)-g(3,k))*R(6)
            
                  ! See comment (B) above...
                  temp = g(4,k)-temp
                  b(1,k) = b(1,k)+C(4)*temp
                  b(2,k) = b(2,k)+C(5)*temp
                  b(3,k) = b(3,k)+C(6)*temp
                  b(4,k) = b(4,k)+temp
               end do
            case (5)
               do k = 1, nv
                  q = (fj(k)-f0(k))/s
                  
                  ! See comment (A) above...
                  temp = g(5,k)
                  g(5,k) = ((((q-g(1,k))*R(7)-g(2,k))*R(8)-g(3,k))*R(9) &
                       -g(4,k))*R(10)
                  
                  ! See comment (B) above...
                  temp = g(5,k)-temp
                  b(1,k) = b(1,k)+C(7)*temp
                  b(2,k) = b(2,k)+C(8)*temp
                  b(3,k) = b(3,k)+C(9)*temp
                  b(4,k) = b(4,k)+C(10)*temp
                  b(5,k) = b(5,k)+temp
               end do
            case (6)
               do k = 1, nv
                  q = (fj(k)-f0(k))/s
                  
                  ! See comment (A) above...
                  temp = g(6,k)
                  g(6,k) = (((((q-g(1,k))*R(11)-g(2,k))*R(12) &
                       -g(3,k))*R(13)-g(4,k))*R(14)-g(5,k))*R(15)
                  
                  ! See comment (B) above...
                  temp = g(6,k)-temp
                  b(1,k) = b(1,k)+C(11)*temp
                  b(2,k) = b(2,k)+C(12)*temp
                  b(3,k) = b(3,k)+C(13)*temp
                  b(4,k) = b(4,k)+C(14)*temp
                  b(5,k) = b(5,k)+C(15)*temp
                  b(6,k) = b(6,k)+temp
               end do
            case (7)
               do k = 1, nv
                  q = (fj(k)-f0(k))/s
                  
                  ! See comment (A) above...
                  temp = g(7,k)
                  g(7,k) = ((((((q-g(1,k))*R(16)-g(2,k))*R(17) &
                       -g(3,k))*R(18)-g(4,k))*R(19)-g(5,k))*R(20) &
                       -g(6,k))*R(21)
                  
                  ! See comment (B) above...
                  temp = g(7,k)-temp
                  b(1,k) = b(1,k)+C(16)*temp
                  b(2,k) = b(2,k)+C(17)*temp
                  b(3,k) = b(3,k)+C(18)*temp
                  b(4,k) = b(4,k)+C(19)*temp
                  b(5,k) = b(5,k)+C(20)*temp
                  b(6,k) = b(6,k)+C(21)*temp
                  b(7,k) = b(7,k)+temp
               end do
            end select
         end do substep_loop
         
         if (nes .or. i < ni) cycle better_B_loop
         
         ! Integration of sequence is over. Next is sequence size control.
         hv = ZERO
         do k = 1, nv
            hv = max(hv,abs(b(7,k)))
         end do
         hv = hv*WC(7)/hval**7
      end do better_B_loop
      
   ! If this is the 1st sequence... we still have to adjust the
   ! time step
   if (.not. nsf) then
     if (.not. nes) hp = sign((SS/hv)**PW,tf)
     if (nes .or. hp/h > ONE) then
        if (nes) hp = h0
        nsf = .true.
     else
        hp = 0.8_DP*hp
        ncount = ncount+1
        if (ncount > MAX_NCOUNT) then
           write(*,*)
           write(*,*) '*************************************'
           write(*,*) 'NCOUNT > ', MAX_NCOUNT
           write(*,*) 'Cannot find an optimal sequence size.'
           write(*,*) 'RA15 returns to the caller.'
           write(*,*) '*************************************'
           write(*,*)
           ! Exiting the main loop should be the same as RETURN.
           exit main_loop
           !return
        end if
        
        ! Restart with HP = 0.8*H if new HP is smaller than original
        ! H on 1st sequence.
        cycle main_loop
     end if
   end if
   
   ! This loop finds new X and V values at end of sequence.
   ! Eqs. (11), (12).
   xprev = x
   vprev = v
   do k = 1, nv
      x(k) = x(k)+v(k)*h+h2*(f0(k)*WC0+b(1,k)*WC(1)+b(2,k)*WC(2) &
           +b(3,k)*WC(3)+b(4,k)*WC(4)+b(5,k)*WC(5)+b(6,k)*WC(6) &
           +b(7,k)*WC(7))
      
      ! If equations are 1st order, skip to compute y' (aka V)
      ! at end of sequence, i.e. cycle..
      if (ncl) cycle
      
      v(k) = v(k)+h*(f0(k)+b(1,k)*UC(1)+b(2,k)*UC(2)+b(3,k)*UC(3) &
           +b(4,k)*UC(4)+b(5,k)*UC(5)+b(6,k)*UC(6)+b(7,k)*UC(7))
   end do
   ! We have done a sequence and can update current time and
   ! sequence counter.
   tprev = tm
   tm = tm + h
   tmf = tmf + h
   ns = ns + 1
!   write(*,*) 'RA15: tprev = ',tprev
!   write(*,*) 'RA15: tmf = ',tmf
!   write(*,*) 'RA15: nf = ',nf
!   write(*,*) 'RA15: x = ',x
!   write(*,*) 'RA15: v = ',v
!   read(*,*)
   
   if (itrp) then
     ! INTERPOLATION MODE.
     ! Check if output is requested in the last performed sequence. If this
     ! is the final sequence of the integration, don't interpolate (we have
     ! already calculated x,v in the previous loop).
     do
       if ((abs(tip(jip)/tmf) - ONE) > ZERO .or. jip >= mip) exit
       hout = tip(jip) - tprev
       hout2 = hout*hout
       sout = hout/h

       ! These are just the predictor equations. The sum was expanded and
       ! rewritten in a more human-readable way (we're in 2016 after all).
       if (ncl) then
         xout = xprev + hout*(f0 + b(1,:)*sout*WC(1) + b(2,:)*sout**2*WC(2) + &
         & b(3,:)*sout**3*WC(3) + b(4,:)*sout**4*WC(4) + b(5,:)*sout**5*WC(5) +&
         & b(6,:)*sout**6*WC(6) + b(7,:)*sout**7*WC(7))
         
         ! If the equations are Type 1, there is no velocity predictor.
         vout = 0._dp
         
       else
         xout = xprev + vprev*hout + hout2*(f0*WC0 + b(1,:)*sout*WC(1) + &
         & b(2,:)*sout**2*WC(2) + b(3,:)*sout**3*WC(3) + b(4,:)*sout**4*WC(4)&
         & + b(5,:)*sout**5*WC(5) + b(6,:)*sout**6*WC(6) + b(7,:)*sout**7*WC(7))
         
         vout = vprev + hout*(f0 + b(1,:)*sout*UC(1) + b(2,:)*sout**2*UC(2) +&
         & b(3,:)*sout**3*UC(3) + b(4,:)*sout**4*UC(4) + b(5,:)*sout**5*UC(5) &
         & + b(6,:)*sout**6*UC(6) + b(7,:)*sout**7*UC(7))
         
       end if
              
       ! SAVE INTO AUXY
       select case (nclass)
       case (1)
         AUXY(jip,1:1+nv) = [tip(jip),xout]
       case (2,-2)
         AUXY(jip,1:1+2*nv) = [tip(jip),xout,vout]
       end select
       
       ! Save past time step into output
       if (present(HT)) AUXH(jip,:) = [tprev,h]
       
       jip = jip + 1
       
     end do
   
   else
     ! Save current position and velocity in AUXY.
     select case (nclass)
     case (1)
       AUXY(1+ns,1:1+nv) = [tmf,x]
     case (2,-2)
       AUXY(1+ns,1:1+2*nv) = [tmf,x,v]
     end select
     
     ! Save past time step into output
     if (present(HT)) AUXH(1+ns,:) = [tprev,h]
     
   end if
   
   ! Return if done.
   if (nper) then
      ! Save trajectory into output array YT. In Fortran 2003, the following
      ! statements (re)allocates automatically YT depending on the dimension
      ! of AUXY.
      if (itrp) then
       select case (nclass)
       case (1)
         AUXY(mip,1:1+nv) = [tmf,x]
       case (2,-2)
         AUXY(mip,1:1+2*nv) = [tmf,x,v]
       end select
        YT = AUXY(1:mip,:)
      else
        YT = AUXY(1:1+ns,:)
      end if
      
      ! On exit, rdiag(1) contains the last computed (signed) sequence size
      rdiag(1) = h
      
      ! Save integration diagnostics
      idiag(1) = nf
      idiag(2) = ns
      
      ! Save step size history
      if (present(HT)) then
        if (itrp) then
          AUXH(mip,:) = [tprev,h]
          HT = AUXH(1:mip,:)
        else
          AUXH(1+ns,:) = [tprev,h]
          HT = AUXH(1:1+ns,:)
        end if
      end if
      
      ! Exiting the main loop should be the same as RETURN.
      exit main_loop
      !return
   end if
   
   ! Control on size of next sequence and adjust last sequence to exactly
   ! cover the integration span. NPER = .TRUE. set on last sequence.
!   call force(tmf,x,v,f0)
   select case (nclass)
   case(1,-2)
     call force(nv,tmf,x,f0)
   case(2)
     call force(nv,tmf,x,v,f0)
   end select
   nf = nf+1
   
   if (nes) then
      hp = h0
   else
      hp = sign((SS/hv)**PW,tf)
      if (hp/h > SR) hp = h*SR
   end if
   if (abs(tm+hp) >= abs(tf)-EPS_TF_MATCH) then
      hp = tf-tm
      nper = .true.
   end if
   
   ! Now predict B−values for next step using Eqs. (13). Values from the
   ! preceding sequence were saved in the E−matrix. The correction BD
   ! is applied in the following loop as described in Sec. 2.5.
   q = hp/h
   
   ! To avoid re−computation of the same expession (q**2, q**3,...)
   ! for each K...
   q2 = q*q!q**2
   q3 = q*q2!q**3
   q4 = q2*q2!q**4
   q5 = q2*q3!q**5
   q6 = q3*q3!q**6
   q7 = q3*q4!q**7
   do k = 1, nv
      ! If we have done at least TWO sequences..
      if (ns /= 1) then
         do j = 1, NSTEPS
            bd(j) = b(j,k)-e(j,k)
         end do
      end if
      
      e(1,k) = q*(b(1,k)+Z2*b(2,k)+Z3*b(3,k)+Z4*b(4,k)+Z5*b(5,k) &
           +Z6*b(6,k)+Z7*b(7,k))
      e(2,k) = q2*(b(2,k)+Z3*b(3,k)+Z6*b(4,k)+Z10*b(5,k)+Z15*b(6,k) &
           +Z21*b(7,k))
      e(3,k) = q3*(b(3,k)+Z4*b(4,k)+Z10*b(5,k)+Z20*b(6,k)+Z35*b(7,k))
      e(4,k) = q4*(b(4,k)+Z5*b(5,k)+Z15*b(6,k)+Z35*b(7,k))
      e(5,k) = q5*(b(5,k)+Z6*b(6,k)+Z21*b(7,k))
      e(6,k) = q6*(b(6,k)+Z7*b(7,k))
      e(7,k) = q7*b(7,k)
      
      ! Apply the correction.. Notice that when we have done ONLY
      ! one sequence (NS == 1), BD == 0 from its initialization, i.e.
      ! we are doing B = E. It is only when NS > 1 that we are applying
      ! the correction BD.
      do j = 1, NSTEPS
         b(j,k) = e(j,k)+bd(j)
      end do
   end do
   
   ! Two iterations for every sequence. (Use 3 for 23rd and 27th order.)
   ni = 2
   end do every_sequence_loop
end do main_loop
end subroutine ra15
    
end module EVERHART

!do
!  if ((abs(tip(jip)/tmf) - ONE) > ZERO .or. jip >= mip) exit
!  hout = abs(tip(jip) - tprev)
!  hout2 = hout*hout
!  if (ncl) hout2 = hout
!  xout = xprev + vprev*hout + hout2*(f0*WC0+b(1,:)*WC(1)+b(2,:)*WC(2)&
!  &+b(3,:)*WC(3)+b(4,:)*WC(4)+b(5,:)*WC(5)+b(6,:)*WC(6)+b(7,:)*WC(7))
!  
!  if (.not.(ncl)) then
!    vout = vprev + hout*(f0+b(1,:)*UC(1)+b(2,:)*UC(2)+b(3,:)*UC(3)&
!    &+b(4,:)*UC(4)+b(5,:)*UC(5)+b(6,:)*UC(6)+b(7,:)*UC(7))
!  end if
!  
!  ! SAVE INTO AUXY
!  select case (nclass)
!  case (1)
!    AUXY(jip,1:1+nv) = [tip(jip),xout]
!  case (2,-2)
!    AUXY(jip,1:1+2*nv) = [tip(jip),xout,vout]
!  end select
!  
!  jip = jip + 1
!  
!end do

!       qout = sout
!       if (ncl) qout = ONE
!       if (ncl) hout2 = hout
!       tempv = WC(3)*b(3,:) + sout*(WC(4)*b(4,:) + sout*(WC(5)*b(5,:) &
!               + sout*(WC(6)*b(6,:) + sout*WC(7)*b(7,:))))
!       xout = xprev + qout*(hout*vprev + hout2*sout*(f0*WC0 + sout*(WC(1)*b(1,:) & 
!            + sout*(WC(2)*b(2,:) + sout*tempv))))
!       
!       if (.not.(npq)) then
!         tempv = UC(3)*b(3,:) + sout*(UC(4)*b(4,:) + sout*(UC(5)*b(5,:) &
!               + sout*(UC(6)*b(6,:) + sout*UC(7)*b(7,:))))
!         vout = vprev + sout*hout*(f0 + sout*(UC(1)*b(1,:) + sout*(UC(2)*b(2,:)&
!               + sout*tempv)))
!       end if
!         xout = xprev + hout*(f0 + b(1,:)*sout*WC(2,:) + b(2,:)*sout**2/3._dp + &
!         & b(3,:)*sout**3/4._dp + b(4,:)*sout**4/5._dp + b(5,:)*sout**5/6._dp +&
!         & b(6,:)*sout**6/7._dp + b(7,:)*sout**7/8._dp)
