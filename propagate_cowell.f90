module PROP_COWELL_MOD

use KINDS, only: ik,dk,qk

implicit none

contains

subroutine QLEG_COW_LSODAR(R_i,V_i,JD_i,JD_f,dt_D,tol,Yt,istate,jroot,ind_sw,calls)
! Description:
!    Propagates from JD_i to JD_f starting from (R_i,V_i), using Cowell's formulation with timestep
!    dt_D (dimensional).
!
! ==============================================================================

! MODULES
use AUXILIARIES,  only: inSoI
use CONSTANTS,    only: DU,TU,muEarth,smaEarth,wEarth,secsPerDay
use EQS_COWELL,   only: QCOWELL_RHS
use EVTS_COWELL,  only: QCOWELL_EVT
use SETTINGS,     only: mxsteps
use THIRDBODIES,  only: QPOS_VEL_CIRC

! VARIABLES
implicit none
! Arguments IN
real(qk),intent(in)  ::  R_i(1:3),V_i(1:3)     ! Initial state [km,km/s]
real(qk),intent(in)  ::  JD_i,JD_f             ! Initial and final times [day]
real(qk),intent(in)  ::  dt_D                  ! Time step [day]
real(qk),intent(in)  ::  tol                   ! Integration tolerance
! Arguments OUT
real(qk),allocatable,intent(out)  ::  Yt(:,:)  ! Trajectory array
integer(ik),intent(out)  ::  istate            ! Integration diagnostics
integer(ik),intent(out)  ::  jroot(1)          ! Roots
integer(ik),intent(out)  ::  ind_sw            ! Index of switch time
integer(ik),intent(out)  ::  calls             ! Function calls for this leg

! States and times
real(qk)  ::  y(1:6)
real(qk)  ::  t,dt,t_f
real(qk)  ::  lon_Earth,r_Earth(1:3),v_Earth(1:3)
real(qk)  ::  r_helio(1:3),v_helio(1:3),t_D
real(qk)  ::  r_geo(1:3),v_geo(1:3)

! Memory
real(qk),allocatable  ::  CHUNK(:,:)

! Diagnostics
real(qk)  ::  sma,ecc

! LSODAR variables
integer(ik),parameter  ::  neq = 6_ik
integer(ik),parameter  ::  itol = 1_ik
integer(ik),parameter  ::  itask = 1_ik
integer(ik),parameter  ::  iopt = 1_ik
integer(ik),parameter  ::  lrw = 200_ik, liw = 200_ik
integer(ik),parameter  ::  jt = 2
integer(ik)            ::  nevts
integer                ::  steps
real(qk)     ::  rwork(1:lrw)
integer(ik)  ::  iwork(1:liw)
external DLSODAR

! Misc.
integer(ik)    ::  i_int
real(qk),parameter  ::  zero = epsilon(0._qk)

! ==============================================================================

! ==============================================================================
! 01. INITIALIZATIONS
! ==============================================================================

! Non-dimensionalizations
t = JD_i*secsPerDay*TU
t_f = JD_f*secsPerDay*TU
dt = dt_D*secsPerDay*TU
steps = (t_f - t)/dt
y(1:3) = R_i/DU
y(4:6) = V_i/(DU*TU)

! LSODAR variables
istate = 1_ik
iwork = 0_ik
rwork = 0._qk
iwork(6) = 20000_ik
if (inSoI) then  ! Event detection is necessary only in phase CE
  nevts = 1
else
  nevts = 0
end if

! Misc.
i_int = 1_ik

! If Yt is allocated, deallocate it to avoid problems
if (allocated(Yt)) deallocate(Yt)

! Allocate temporary array for trajectory data
if (allocated(CHUNK)) deallocate(CHUNK)
allocate(CHUNK(1:2*mxsteps,1:7))
CHUNK(1,1:7) = STATE_HELIO(t,y,inSoI)

! ==============================================================================
! 02. INTEGRATION
! ==============================================================================

do
    call DLSODAR(QCOWELL_RHS,neq,y,t,t+dt,itol,tol,tol,itask,istate,iopt,rwork,&
    &lrw,iwork,liw,FAKE,jt,QCOWELL_EVT,nevts,jroot)
    if (istate < 0) exit
    
    ! ==========================================================================
    ! 02a. Save current time and heliocentric state.
    ! ==========================================================================
    
    CHUNK(i_int+1,1:7) = STATE_HELIO(t,y,inSoI)

    ! ==========================================================================
    ! 02b. Exit conditions.
    ! ==========================================================================
    
    if (inSoI) then
        if (istate == 3) then
            ind_sw = i_int + 1
            Yt = CHUNK(1:ind_sw,:)
            exit
        end if
    else
        if (i_int == steps - 1) then
            dt = t_f - t
        else if (i_int == steps) then
            Yt = CHUNK(1:steps+1,:)
            exit
        end if
    end if
    
    ! Reaching of maximum number of steps
    if (i_int == mxsteps) then
        istate = -11_ik
        exit
    end if
    
    i_int = i_int + 1
    
end do

if (allocated(CHUNK)) deallocate(CHUNK)


contains

  subroutine FAKE()
  implicit none
  end subroutine FAKE
  
  function STATE_HELIO(t,y,inSoI)
    real(qk),intent(in) ::  t,y(1:6)
    logical,intent(in)  ::  inSoI
    real(qk)  ::  STATE_HELIO(1:7)
    real(qk)  ::  t_D,r_D(1:3),v_D(1:3)
    real(qk)  ::  yD_Earth(1:6)
    
    t_D = t/TU
    STATE_HELIO(1) = t_D
    
    if (inSoI) then
      yD_Earth = QPOS_VEL_CIRC(t_D,wEarth,smaEarth)
      STATE_HELIO(2:4) = yD_Earth(1:3) + y(1:3)*DU 
      STATE_HELIO(5:7) = yD_Earth(4:6) + y(4:6)*DU*TU
    
    else
      STATE_HELIO(2:4) = y(1:3)*DU
      STATE_HELIO(5:7) = y(4:6)*DU*TU
    
    end if
    
  end function STATE_HELIO
  
end subroutine QLEG_COW_LSODAR



end module PROP_COWELL_MOD

!if (iwork(19) == 2_ik .and. .not.(meth_switch)) then
!  if (dt > 0._qk) then
!      write(*,'(a,l)') 'BDF method used in the reference propagation (fwd)'
!      meth_flag(2) = 1_ik
!  else
!      write(*,'(a,l)') 'BDF method used in the reference propagation (bwd)'
!      meth_flag(1) = 2_ik
!  end if
!  meth_switch = .true.
!end if
