module REFERENCE_TRAJECTORY_MOD
! Description:
!    Contains the REFERENCE_TRAJECTORY subroutine.
! 
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
! 
! ==============================================================================

use KINDS, only: ik,dk,qk
implicit none

contains

subroutine REFERENCE_TRAJECTORY_COWELL(R_i,V_i,JD_i,JD_f,tol,Yt,JD_sw,ind_sw,istate)
! Description:
!    Computes the reference trajectory in quadruple precision, using trajectory
!    splitting with Cowell equations integrated by DLSODAR.
!    Provides the reference trajectory and times of switch IN and OUT.
!
! ==============================================================================

! MODULES
use AUXILIARIES,     only: inSoI
use CONSTANTS,       only: muSun,muEarth,smaEarth,secsPerDay,pi
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
! States and trajectories
real(qk)              ::  R_helio(1:3),V_helio(1:3)
real(qk)              ::  R_geo(1:3),V_geo(1:3)
real(qk),allocatable  ::  traj_H(:,:),traj_CE(:,:)
integer               ::  length_CE,length_H,length
! Integration
real(qk)            ::  dir
integer(ik)         ::  jroot(1)
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
call QLEG_COW_LSODAR(R_geo,V_geo,JD_i,JD_f,dir*dt_CE,tol,traj_CE,istate,jroot,ind_sw)
if (istate < 0) return

! Switch IN-> OUT:
inSoI = .false.
length_CE = size(traj_CE,1)
JD_sw = traj_CE(length_CE,1)/secsPerDay
R_helio = traj_CE(length_CE,2:4)
V_helio = traj_CE(length_CE,5:7)

call SET_UNITS(R_helio,muSun)

! Propagate OUT:
call QLEG_COW_LSODAR(R_helio,V_helio,JD_sw,JD_f,dir*dt_H,tol,traj_H,istate,jroot,temp)
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

end subroutine REFERENCE_TRAJECTORY_COWELL

end module REFERENCE_TRAJECTORY_MOD
