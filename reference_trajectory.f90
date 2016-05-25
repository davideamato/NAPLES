module REFERENCE_TRAJECTORY_MOD

use KINDS, only: ik,dk,qk
implicit none

contains

subroutine REFERENCE_TRAJECTORY_COWELL(R_i,V_i,JD_i,JD_f,tol,Yt,JD_sw,ind_sw,istate)

! MODULES
use AUXILIARIES,     only: inSoI,RSwitch
use CONSTANTS,       only: DU,TU,muSun,muEarth,smaEarth,wEarth,secsPerDay,R_SoI,pi
use CONSTANTS,       only: SET_UNITS
use SETTINGS,        only: dt_H,dt_CE
use PROP_COWELL_MOD, only: QLEG_COW_LSODAR

! VARIABLES
implicit none
! Arguments IN
real(qk),intent(in)    ::  R_i(1:3),V_i(1:3)    ! Initial position and velocity [km,km/s]
real(qk),intent(in)    ::  JD_i,JD_f            ! Initial and final times [days]
real(qk),intent(in)    ::  tol                  ! Integration tolerance
! Arguments OUT
real(qk),allocatable,intent(out)   ::  Yt(:,:)    ! Trajectory array [km,s]
real(qk),intent(out)     ::  JD_sw                      ! Switch time [days]
integer(ik),intent(out)  ::  ind_sw                     ! Index of switch time
integer(ik),intent(out)  ::  istate                     ! Integration diagnostics
! Times
!real(qk)  ::  dt_helio,dt_geo
! States and trajectories
real(qk)              ::  R_helio(1:3),V_helio(1:3)
real(qk)              ::  R_geo(1:3),V_geo(1:3)
real(qk),allocatable  ::  traj_H(:,:),traj_CE(:,:)
integer               ::  length_CE,length_H,length,ii
! Integration
real(qk)            ::  dir
integer(ik)         ::  jroot(1),calls
! Misc
integer :: temp

! ==============================================================================

! ==============================================================================
! 01. INITIALIZATIONS
! ==============================================================================

! In this subroutine, we start propagating from t_CA = 0, i.e. at the time of
! minimum geocentric distance.

! Reference quantities for non-dimensionalization
call SET_UNITS(R_i,muEarth)

dir = sign(1._qk,JD_f - JD_i)

! Initialize geocentric state
R_geo = R_i
V_geo = V_i

! Initialize istate, jroot, switch time index
istate = 1
jroot  = 0
ind_sw = -1

! ===============================================================================
! 02. PROPAGATION WITH ONLINE TRAJECTORY MATCHING (needed to get the switch time)
! ===============================================================================

inSoI = .true.

! Propagate IN
call QLEG_COW_LSODAR(R_geo,V_geo,JD_i,JD_f,dir*dt_CE,tol,traj_CE,istate,jroot,ind_sw,calls)
if (istate < 0) return

! Switch IN-> OUT:
inSoI = .false.
length_CE = size(traj_CE,1)
JD_sw = traj_CE(length_CE,1)/secsPerDay
R_helio = traj_CE(length_CE,2:4)
V_helio = traj_CE(length_CE,5:7)

call SET_UNITS(R_helio,muSun)

! Propagate OUT:
call QLEG_COW_LSODAR(R_helio,V_helio,JD_sw,JD_f,dir*dt_H,tol,traj_H,istate,jroot,temp,calls)
if (istate < 0) return

length_H = size(traj_H,1)

! ==============================================================================
! 03. MERGE TRAJECTORIES AND QUIT
! ==============================================================================

length = length_CE + length_H - 1

if(allocated(Yt)) deallocate(Yt)
allocate(Yt(1:length,1:7))
Yt(1:length_CE-1,:) = traj_CE(1:length_CE-1,:)
Yt(length_CE:length,:) = traj_H(1:length_H,:)

deallocate(traj_CE)
deallocate(traj_H)


!contains
!  
!  function TIMESTEP_HYPER(R_i,V_i)
!  ! Computes the time step to use in the close encounter (CE) phase.
!  ! The duration of the close encounter is computed as the time of flight
!  ! necessary to reach RSwitch in a Keplerian hyperbolic trajectory.
!  ! Reference: Battin R., An Introduction to the Mathematics and Methods of
!  !            Astrodynamics, 1999.
!    real(qk),intent(in)  ::  R_i(1:3), V_i(1:3)
!    real(qk)  ::  TIMESTEP_HYPER
!    real(qk)  ::  Rmag,sma_i,ecc_vec(1:3),ecc_i,GA_switch,hyper_ToF
!    
!    Rmag = sqrt(dot_product(R_i,R_i))
!    sma_i = Rmag*muEarth/(2._qk*muEarth - Rmag*dot_product(V_i,V_i))
!    ecc_vec = (R_i*dot_product(V_i,V_i) - V_i*dot_product(R_i,V_i))/muEarth - R_i/Rmag
!    ecc_i = sqrt(dot_product(ecc_vec,ecc_vec))
!    GA_switch  = acos((ecc_i)/(1._qk - RSwitch/sma_i))
!    hyper_ToF = 2._qk*sqrt(-sma_i**3_ik/muEarth)*(ecc_i*tan(GA_switch) - &
!    &log(tan(.5_qk*GA_switch + pi/4._qk)))
!    TIMESTEP_HYPER = hyper_ToF/steps * ((JD_f - JD_i)/abs((JD_f - JD_i)))
!    
!  end function

end subroutine REFERENCE_TRAJECTORY_COWELL

!dt_helio = ((JD_f - JD_sw)*secsPerDay)/steps

end module REFERENCE_TRAJECTORY_MOD
