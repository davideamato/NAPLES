module EQS_EDROMO

use KINDS, only: ik,dk,qk
implicit none

contains

subroutine DEDROMO_RHS(neq,phi,z,zdot)
! Description:
!    Computes the value of the right-hand side of the equations of motion of the EDromo
!    formulation.
! 
! Versions:
!   17/09/2015: v1.
!
! Author:
!    Davide Amato
!    Space Dynamics Group
!    Technical University of Madrid
!    d.amato@upm.es
!
! ==============================================================================================
!
!                                  VARIABLES AND DECLARATIONS
!
! ==============================================================================================

! MODULES
use AUXILIARIES, only: inSoI
use SETTINGS,    only: flag_time_EDr
use TRANSFORM,   only: DINERT2ORB_EDROMO
use THIRDBODIES, only: DACC3B_ND
use CONSTANTS,   only: DU,TU,wEarth,smaEarth,muEarth,muSun

! VARIABLES
implicit none

! Arguments
integer(ik),intent(in)     ::  neq             ! Number of equations
real(dk),intent(in)        ::  phi             ! EDromo independent variable
real(dk),intent(in)        ::  z(1:neq)        ! EDromo state vector, ND
real(dk),intent(out)       ::  zdot(1:neq)    	! RHS of EoM's, ND

! Auxiliary quantities
real(dk)    ::  sph,cph
real(dk)    ::  rho,rmag,zeta,emme
real(dk)    ::  cnu,snu
real(dk)    ::  enne,L3,wz
real(dk)    ::  aux0,aux1,aux2,aux3,aux4

! State
real(dk)    ::  rV(1:3),t
real(dk)    ::  x_vec(1:3),y_vec(1:3)

! Potential
real(dk)    ::  A,Upot,dUdt,dUdr(1:3)

! Perturbations
real(dk)    ::  lon_Earth,r2(1:3)
real(dk)    ::  f(1:3),p(1:3)

! Debug
integer,save  ::  called = 1

! ==============================================================================================
!
!                                            EXECUTION
!
! ==============================================================================================

! INDEX OF STATE VECTORS

! z(1) = zeta1
! z(2) = zeta2
! z(3) = zeta3
! z(4) = zeta4
! z(5) = zeta5
! z(6) = zeta6
! z(7) = zeta7
! z(8) = zeta8 (physical time,
!                      constant time element,
!                      linear time element,
!                      depending on flag_time)

! ==============================================================================================
! 01. COMPUTE AUXILIARY QUANTITIES (1)
! ==============================================================================================

! Store trig functions
sph = sin(phi)
cph = cos(phi)

rho   = 1._dk - z(1)*cph - z(2)*sph
rmag  = z(3)*rho
zeta  = z(1)*sph - z(2)*cph
emme  = sqrt(1._dk - z(1)**2_ik - z(2)**2_ik)

cnu = (cph - z(1) + (zeta*z(2))/(emme + 1._dk))/rho
snu = (sph - z(2) - (zeta*z(1))/(emme + 1._dk))/rho

! ==============================================================================================
! 02. COMPUTE POSITION IN INERTIAL FRAME AND POTENTIAL
! ==============================================================================================

! Intermediate frame unit vectors
x_vec = 2._dk*[ .5_dk - z(5)**2_ik - z(6)**2_ik,  &
           &  z(4)*z(5) + z(6)*z(7),  &
           &  z(4)*z(6) - z(5)*z(7)]
y_vec = 2._dk*[  z(4)*z(5) - z(6)*z(7),  &
           & .5_dk - z(4)**2_ik - z(6)**2_ik,  &
           &  z(5)*z(6) + z(4)*z(7)]

! Position in inertial frame
rV = rmag*(x_vec*cnu + y_vec*snu)

! Potential
Upot = 0._dk; dUdt = 0._dk; dUdr = 0._dk

!! Velocity in inertial frame
!i_vec = x_vec*cnu + y_vec*snu
!j_vec = y_vec*cnu - x_vec*snu
!v_rad = zeta/(sqrt(z(3))*rho)
!v_tan = sqrt((1. - z(1)**2 - z(2)**2)/(z(3)*rho**2) - 2.*Upot)
!vV = v_rad*i_vec + v_tan*j_vec

! ==============================================================================================
! 03. COMPUTE PERTURBING ACCELERATIONS
! ==============================================================================================

p = 0._dk; f = 0._dk

! Get time
if ( flag_time_EDr == 0_ik ) then
    ! Physical time	
    t = z(8)
elseif  ( flag_time_EDr == 1_ik ) then
    ! Constant time element
    t = z(8) - z(3)**1.5_dk*(zeta - phi)
elseif  ( flag_time_EDr == 2_ik ) then
    ! Linear time element
    t = z(8) - z(3)**1.5_dk*zeta
end if

lon_Earth = wEarth*(t/TU)
if (inSoI) then
    ! Third body: Sun
    r2 = -(smaEarth/DU)*[cos(lon_Earth),sin(lon_Earth),0._dk]
    p  = DACC3B_ND(rV,r2,muEarth,muSun)
else
    ! Third body: Earth
    r2 = (smaEarth/DU)*[cos(lon_Earth),sin(lon_Earth),0._dk]
    p  = DACC3B_ND(rV,r2,muSun,muEarth)
end if

! Inertial -> Orbital frame
p    = DINERT2ORB_EDROMO(p,z,cnu,snu)
!dUdr = DINERT2ORB_EDROMO(dUdr,z,cnu,snu)

f = p + dUdr

! ==============================================================================================
! 04. COMPUTE AUXILIARY QUANTITIES (2)
! ==============================================================================================

enne = sqrt(emme**2_ik - 2._dk*z(3)*rho**2_ik*Upot)

aux0 = (p(1)*zeta + p(2)*enne + dUdt*sqrt(z(3))*rho)
zdot(3) = 2._dk*z(3)**3_ik*aux0
L3 = zdot(3)/(2._dk*z(3))

aux1 = ((2._dk*Upot - f(1)*rmag)*(2._dk - rho + emme)*rmag)/(emme*(emme+1._dk))
aux2 = (L3*zeta*(rho-emme))/(emme*(emme+1._dk))
wz   = (enne - emme)/rho + aux1 + aux2

! ==============================================================================================
! 05. COMPUTE RIGHT-HAND SIDE
! ==============================================================================================

! In-plane
aux3 = (f(1)*rmag - 2._dk*Upot)*rmag

zdot(1) =  aux3*sph + L3*((1._dk+rho)*cph - z(1))
zdot(2) = -aux3*cph + L3*((1._dk+rho)*sph - z(2))

! Out-of-plane
aux4 = (f(3)*rmag**2_ik)/(2._dk*enne)

zdot(4) =  aux4*(z(7)*cnu - z(6)*snu) + 0.5_dk*wz*z(5)
zdot(5) =  aux4*(z(6)*cnu + z(7)*snu) - 0.5_dk*wz*z(4)
zdot(6) = -aux4*(z(5)*cnu - z(4)*snu) + 0.5_dk*wz*z(7)
zdot(7) = -aux4*(z(4)*cnu + z(5)*snu) - 0.5_dk*wz*z(6)

! Time / Time Element
if (flag_time_EDr == 0_ik) then
    zdot(8) = sqrt(z(3))*rmag
else if (flag_time_EDr == 1_ik) then   ! Constant Time Element
    zdot(8) = z(3)**1.5_dk * ( aux3 + (zeta - 1.5_dk*phi)*zdot(3)/z(3) )
else if (flag_time_EDr == 2_ik) then   ! Linear Time Element
    zdot(8) = z(3)**1.5_dk * ( 1._dk + aux3 + 2._dk*L3*zeta )
end if

end subroutine DEDROMO_RHS
    
end module EQS_EDROMO
