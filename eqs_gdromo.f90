module EQS_GDROMO
! Description:
!    Contains the subroutine for computing the RHS of the GDromo equations in
!    double precision.
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

subroutine DGDROMO_RHS(neq,phi,z,zdot)
! Description:
!    Computes the value of the right-hand side of the equations of motion of the
!    GDromo formulation.
!
! ==============================================================================

! MODULES
use AUXILIARIES, only: inSoI
use SETTINGS,    only: flag_time_GDr
use TRANSFORM,   only: DGDROMO2CART_CORE,DINERT2ORB_GDROMO
use THIRDBODIES, only: DACC3B_ND
use CONSTANTS,   only: DU,TU,wEarth,smaEarth,muEarth,muSun

! VARIABLES
implicit none

! Arguments
integer(ik),intent(in)     ::  neq             ! Number of equations
real(dk),intent(in)        ::  phi             ! EDromo independent variable
real(dk),intent(in)        ::  z(1:neq)        ! EDromo state vector, ND
real(dk),intent(out)       ::  zdot(1:neq)     ! RHS of EoM's, ND

! Auxiliary quantities
real(dk)    ::  sph,cph
real(dk)    ::  rho,zeta,g,emme,gamma
real(dk)    ::  rmag,cnu,snu
real(dk)    ::  enne,L3,wz
real(dk)    ::  aux0,aux1,aux2,aux3,aux4
real(dk)    ::  s1,s2

! State
real(dk)    ::  rV(1:3),vV(1:3),t

! Potential
real(dk)    ::  Upot,dUdt,dUdr(1:3)

! Perturbations
real(dk)    ::  p(1:3),f(1:3)
real(dk)    ::  lon_Earth,r2(1:3)

! ==============================================================================

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

! ==============================================================================
! 01. COMPUTE AUXILIARY QUANTITIES (1)
! ==============================================================================

! Store trig functions
sph = sinh(phi)
cph = cosh(phi)

rho   = z(1)*cph + z(2)*sph - 1._dk
zeta  = z(1)*sph + z(2)*cph
gamma = sqrt(z(1)**2_ik + z(2)**2_ik)
emme  = sqrt(z(1)**2_ik - z(2)**2_ik - 1._dk)
g     = sqrt(z(1)**2_ik - z(2)**2_ik)

rmag  = z(3)*rho
s1 = (emme*zeta/rho)/(g*gamma)
s2 = (emme**2_ik/rho - 1._dk)/(g*gamma)
cnu = z(1)*s2 - z(2)*s1
snu = z(1)*s1 + z(2)*s2

! ==============================================================================
! 02. COMPUTE POTENTIAL, CARTESIAN COORDINATES
! ==============================================================================

Upot = 0._dk; dUdt = 0._dk; dUdr = 0._dk
call DGDROMO2CART_CORE(z,rmag,cnu,snu,zeta,rho,emme,rV,vV,Upot)

! ==============================================================================
! 03. COMPUTE PERTURBATIONS
! ==============================================================================

p = 0._dk; f = 0._dk

! Get time
if ( flag_time_GDr == 0_ik ) then
	! Physical time	
    t = z(8)
elseif  ( flag_time_GDr == 1_ik ) then
    ! Constant time element
    t = z(8) + z(3)**1.5_dk*(zeta - phi)
elseif  ( flag_time_GDr == 2_ik ) then
    ! Linear time element
    t = z(8) + z(3)**1.5_dk*zeta
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
p    = DINERT2ORB_GDROMO(p,z,cnu,snu)
dUdr = DINERT2ORB_GDROMO(dUdr,z,cnu,snu)

f = p + dUdr

! ==============================================================================
! 04. COMPUTE AUXILIARY QUANTITIES (2)
! ==============================================================================

enne = sqrt(emme**2_ik - 2._dk*z(3)*rho**2_ik*Upot)

aux0 = (p(1)*zeta + p(2)*enne + dUdt*sqrt(z(3))*rho)
zdot(3) = -2._dk*z(3)**3_ik * aux0
L3 = zdot(3)/(2._dk*z(3))

aux1 = (rmag/emme) * (f(1)*rmag - 2._dk*Upot) * (1._dk + (rho + 1._dk)*(1._dk/g**2_ik + emme/gamma**2_ik))
aux2 = zeta*L3 * ((rho - 1._dk)/gamma**2_ik + (emme + rho/emme)/g**2_ik)
wz   = (enne - emme)/rho - aux1 - aux2

! ==============================================================================
! 05. COMPUTE RIGHT-HAND SIDE
! ==============================================================================

! In-plane
aux3 = (f(1)*rmag - 2._dk*Upot)*rmag

zdot(1) = -aux3*sph + L3*((1._dk-rho)*cph - z(1))
zdot(2) =  aux3*cph + L3*((rho-1._dk)*sph - z(2))

! Out-of-plane
aux4 = (f(3)*rmag**2_ik)/(2._dk*enne)

zdot(4) =  aux4*(z(7)*cnu - z(6)*snu) + 0.5_dk*wz*z(5)
zdot(5) =  aux4*(z(6)*cnu + z(7)*snu) - 0.5_dk*wz*z(4)
zdot(6) = -aux4*(z(5)*cnu - z(4)*snu) + 0.5_dk*wz*z(7)
zdot(7) = -aux4*(z(4)*cnu + z(5)*snu) - 0.5_dk*wz*z(6)

! Time / Time Element
if (flag_time_GDr == 0_ik) then
    zdot(8) = sqrt(z(3))*rmag
else if (flag_time_GDr == 1_ik) then   ! Constant Time Element
    zdot(8) = z(3)**2.5_dk*(z(3)*aux0*(2._dk*zeta - 3._dk*phi) - aux3/z(3))
else if (flag_time_GDr == 2_ik) then   ! Linear Time Element
    zdot(8) = z(3)**2.5_dk*(2._dk*z(3)*zeta*aux0 - aux3/z(3)) - z(3)**1.5_dk
end if

end subroutine DGDROMO_RHS
    
end module EQS_GDROMO
