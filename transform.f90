module TRANSFORM

use KINDS, only: ik,dk,qk
implicit none

contains

function DKSMAT(u)
! Description:
!    Computes the KS-matrix from the KS extended state vector.
! 
! ==============================================================================================

! VARIABLES
implicit none
! Arguments
real(dk),intent(in)  ::  u(:)                ! KS extended state vector
real(dk)             ::  DKSMAT(1:4,1:4)     ! KS-matrix

! ==============================================================================================
!                                            EXECUTION
! ==============================================================================================

DKSMAT = RESHAPE([ u(1),  u(2),  u(3),  u(4)&
                &,-u(2),  u(1),  u(4), -u(3)&
                &,-u(3), -u(4),  u(1),  u(2)&
                &, u(4), -u(3),  u(2), -u(1)],[4,4])

end function DKSMAT 

subroutine DKS2CART(u,x,xdot)
! Description:
!    Computes the Cartesian state (ND) from the KS state vector.
! 
! ==============================================================================================
!
!                                  VARIABLES AND DECLARATIONS
!
! ==============================================================================================

! VARIABLES
implicit none
! Arguments
real(dk),intent(in)   ::  u(:)                  ! KS (extended) state vector
real(dk),intent(out)  ::  x(:),xdot(:)          ! Position and velocity, ND
! Locals
real(dk)              ::  r                     ! Radius magnitude, ND

! ==============================================================================================

! Position
x(1) = u(1)**2_ik - u(2)**2_ik - u(3)**2_ik + u(4)**2_ik
x(2) = 2._dk*(u(1)*u(2) - u(3)*u(4))
x(3) = 2._dk*(u(1)*u(3) + u(2)*u(4))
r    = u(1)**2_ik + u(2)**2_ik + u(3)**2_ik + u(4)**2_ik

! Velocity
xdot(1) = 2._dk*(u(1)*u(5) - u(2)*u(6) - u(3)*u(7) + u(4)*u(8))/r
xdot(2) = 2._dk*(u(2)*u(5) + u(1)*u(6) - u(4)*u(7) - u(3)*u(8))/r
xdot(3) = 2._dk*(u(3)*u(5) + u(4)*u(6) + u(1)*u(7) + u(2)*u(8))/r

! If x, xdot are in R^4 their fourth component is null.
if (size(x,1) > 3)    x(4:) = 0._dk
if (size(xdot,1) > 3) xdot(4:) = 0._dk

end subroutine DKS2CART

function DTO_ENERGY(z)
! Description:
!    Transforms the third element of the Dromo state vector from the inverse of the specific
!    angular momentum to the total energy.
!
! ==============================================================================================
!                                  VARIABLES AND DECLARATIONS
! ==============================================================================================

! VARIABLES
implicit none
! Function definition
real(dk)              ::  DTO_ENERGY(1:8)
! Arguments
real(dk),intent(in)   ::  z(1:8)

! ==============================================================================================
!                                            EXECUTION
! ==============================================================================================

DTO_ENERGY(:) = z(:)
DTO_ENERGY(3) = .5_dk*(z(1)**2_ik + z(2)**2_ik - z(3)**2_ik)
end function DTO_ENERGY

function DTO_MOMENTUM(z)
! Description:
!    Transforms the third element of the Dromo state vector from the total energy
!    to the inverse of the specific angular momentum.
!
! ==============================================================================================
!                                  VARIABLES AND DECLARATIONS
! ==============================================================================================

! VARIABLES
implicit none
! Function definition
real(dk)             ::  DTO_MOMENTUM(1:8)
! Arguments
real(dk),intent(in)  ::  z(1:8)

! ==============================================================================================
!                                            EXECUTION
! ==============================================================================================

DTO_MOMENTUM(:) = z(:)
DTO_MOMENTUM(3) = sqrt(z(1)**2_ik + z(2)**2_ik - 2._dk*z(3))
end function DTO_MOMENTUM

subroutine DDROMO2CART(z,phi,phi0,RV,VV,pot,propagateEnergy)
! Description:
!    Transforms from Dromo state vector to Cartesian coordinates. Optionally, it takes into
!    account the potential deriving from Earth's oblateness.
!
! ==============================================================================================
!
!                                  VARIABLES AND DECLARATIONS
!
! ==============================================================================================

! VARIABLES
implicit none
! Arguments IN
real(dk),intent(in)   ::  z(1:8)            ! Dromo(P) state vector
real(dk),intent(in)   ::  phi               ! Dromo independent variable
real(dk),intent(in)   ::  phi0              ! Dromo independent variable - initial value
logical,intent(in)    ::  propagateEnergy   ! .TRUE. if z(3) is the total energy.
real(dk),intent(in)   ::  pot               ! Potential
! Arguments OUT
real(dk),intent(out)  ::  RV(1:3), VV(1:3)  ! Cartesian position and velocity.
! Locals
real(dk)              ::  r, u, s, lam, sdp, cdp
real(dk)              ::  z3               ! Inverse of specific angular momentum.

! ==============================================================================================

! ==============================================================================================
! 01. TOTAL ENERGY / SPECIFIC ANGULAR MOMENTUM
! ==============================================================================================

if (propagateEnergy) then
    z3 = sqrt(z(1)**2_ik + z(2)**2_ik - 2._dk*z(3))
else
    z3 = z(3)
end if

! ==============================================================================================
! 02. COMPUTATION OF AUXILIARIES
! ==============================================================================================

sdp = sin(phi-phi0); cdp = cos(phi-phi0)
s = z3 + z(1)*cos(phi) + z(2)*sin(phi)
r = 1._dk/(z3*s)
u = z(1)*sin(phi) - z(2)*cos(phi)

! ==============================================================================================
! 03. COMPUTATION OF STATE
! ==============================================================================================

lam = sqrt(s**2_ik - 2._dk*pot)

RV(1) =    r*( ( 1._dk - 2._dk*z(5)**2_ik - 2._dk*z(6)**2_ik )*cdp &
&+ 2._dk*( z(4)*z(5) - z(6)*z(7) )*sdp ) 
RV(2) =    r*( ( 1._dk - 2._dk*z(4)**2_ik - 2._dk*z(6)**2_ik )*sdp &
&+ 2._dk*( z(4)*z(5) + z(6)*z(7) )*cdp )
RV(3) = 2._dk*r*( ( z(5)*z(6) + z(4)*z(7) )*sdp + ( z(4)*z(6) - z(5)*z(7) )*cdp )

VV(1) = ( 2._dk*lam*( z(4)*z(5) - z(6)*z(7) ) + u*( 1._dk - 2._dk*z(5)**2_ik - 2._dk*z(6)**2_ik ) )*cdp + &
      & ( 2._dk*u*( z(4)*z(5) - z(6)*z(7) ) -&
      & lam*( 1._dk - 2._dk*z(5)**2_ik - 2._dk*z(6)**2_ik ) )*sdp
VV(2) = ( lam*( 1._dk - 2._dk*z(4)**2_ik - 2._dk*z(6)**2_ik ) + 2._dk*u*( z(4)*z(5) + z(6)*z(7) ) )*cdp + &
      & ( u*( 1._dk - 2._dk*z(4)**2_ik - 2._dk*z(6)**2_ik ) -&
      & 2._dk*lam*( z(4)*z(5) + z(6)*z(7) ) )*sdp
VV(3) = 2._dk*( lam*( z(5)*z(6) + z(4)*z(7) ) + u*( z(4)*z(6) - z(5)*z(7) ) )*cdp + &
      & 2._dk*( u*( z(5)*z(6) + z(4)*z(7) ) - lam*( z(4)*z(6) - z(5)*z(7) ) )*sdp
end subroutine DDROMO2CART

function DINERT2ORB(vI,z,phi,phi0)
! Description:
!    Transforms a vector vI from orbital to inertial axes through a rotation matrix obtained
!    fromo Dromo(P) elements.
!
! ==============================================================================================
        
! VARIABLES
implicit none

! Arguments IN
real(dk),intent(in)          ::  vI(1:3)
real(dk),intent(in)          ::  z(:)
real(dk),intent(in)          ::  phi, phi0
! Function definition
real(dk)                     ::  DINERT2ORB(1:3)
! Locals
real(dk)                     ::  Q0T(1:3,1:3), MphiT(1:3,1:3), QIR(1:3,1:3)
real(dk)                     ::  dphi

! ==============================================================================================

dphi  = phi - phi0

! Compute rotation matrix QIR (Inertial -> Orbital)
MphiT = reshape( (/ cos(dphi), -sin(dphi), 0._dk, sin(dphi), cos(dphi), 0._dk,&
& 0._dk, 0._dk, 1._dk /), (/ 3, 3 /) )
Q0T   = reshape( (/ 1._dk - 2._dk*z(5)**2_ik - 2._dk*z(6)**2_ik,&
& 2._dk*z(4)*z(5) - 2._dk*z(6)*z(7), 2._dk*z(4)*z(6) + 2._dk*z(5)*z(7),&
& 2._dk*z(4)*z(5) + 2._dk*z(6)*z(7),&
& 1._dk - 2._dk*z(4)**2_ik - 2._dk*z(6)**2_ik, 2._dk*z(5)*z(6) - 2._dk*z(4)*z(7),&
& 2._dk*z(4)*z(6) - 2._dk*z(5)*z(7), 2._dk*z(5)*z(6) + 2._dk*z(4)*z(7), 1._dk -&
& 2._dk*z(4)**2_ik - 2._dk*z(5)**2_ik /),&
& (/ 3, 3 /) )
QIR   = matmul(MphiT, Q0T)

! Rotate VI into orbital axes
DINERT2ORB = matmul(QIR,vI)
end function DINERT2ORB


function DINERT2ORB_EDROMO(vI,z,cnu,snu)
! Description:
!    Transforms a vector vI from inertial to orbital axes through a rotation matrix obtained
!    fromo EDromo(P) elements.
!
! ==============================================================================================
!
!                                  VARIABLES AND DECLARATIONS
!
! ==============================================================================================

! MODULES
        
! VARIABLES

implicit none

! Arguments
real(dk),intent(in)  ::  vI(1:3)
real(dk),intent(in)  ::  z(:)
real(dk),intent(in)  ::  cnu,snu
real(dk)             ::  DINERT2ORB_EDROMO(1:3)
! Locals
real(dk)         ::  R1(1:3,1:3)    ! Inertial -> Intermediate rotation matrix
real(dk)         ::  R2(1:3,1:3)    ! Intermediate -> Orbital rotation matrix
real(dk)         ::  RTOT(1:3,1:3)  ! Inertial -> Orbital rotation matrix (= R2*R1)
! Auxiliary variables
real(dk)  ::  z42,z52,z62,z4z5,z6z7,z4z6,z5z7,z5z6,z4z7

! =============EXECUTION============

z42 = z(4)**2
z52 = z(5)**2
z62 = z(6)**2
z4z5 = z(4)*z(5)
z6z7 = z(6)*z(7)
z4z6 = z(4)*z(6)
z5z7 = z(5)*z(7)
z5z6 = z(5)*z(6)
z4z7 = z(4)*z(7)

R1 = 2._dk*reshape([ .5_dk - z52 - z62, z4z5 + z6z7, z4z6 - z5z7,&
            &  z4z5 - z6z7, .5_dk - z42 - z62, z5z6 + z4z7, &
            &  z4z6 + z5z7, z5z6 - z4z7, .5_dk - z42 - z52 ], [3_ik,3_ik])
R1 = transpose(R1)
R2    = reshape([ cnu, -snu, 0._dk,&
              &   snu,  cnu, 0._dk,&
              &   0._dk,    0._dk, 1._dk ],[3_ik,3_ik])
RTOT = matmul(R2,R1)

DINERT2ORB_EDROMO = matmul(RTOT,vI)

end function DINERT2ORB_EDROMO


subroutine DEDROMO2CART(z,phi,Upot,r_vec,v_vec)

! VARIABLES
implicit none
! Arguments IN
real(dk),intent(in)   ::  z(1:8),phi,Upot
! Arguments OUT
real(dk),intent(out)  ::  r_vec(1:3),v_vec(1:3)

! Auxiliaries
real(dk)    ::  sph,cph
real(dk)    ::  rho,rmag,zeta,emme
real(dk)    ::  cnu,snu
real(dk)    ::  x_vec(1:3),y_vec(1:3)
real(dk)    ::  i_vec(1:3),j_vec(1:3)
real(dk)    ::  v_rad,v_tan

! ==============================================================================================

! ==============================================================================================
! 01. COMPUTE AUXILIARY QUANTITIES
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
! 02. COMPUTE POSITION IN INERTIAL FRAME
! ==============================================================================================

! Intermediate frame unit vectors
x_vec = 2._dk*[ .5_dk - z(5)**2_ik - z(6)**2_ik,  &
           &  z(4)*z(5) + z(6)*z(7),  &
           &  z(4)*z(6) - z(5)*z(7)]
y_vec = 2._dk*[  z(4)*z(5) - z(6)*z(7),  &
           & .5_dk - z(4)**2_ik - z(6)**2_ik,  &
           &  z(5)*z(6) + z(4)*z(7)]

! Position in inertial frame
r_vec = rmag*(x_vec*cnu + y_vec*snu)
    
! ==============================================================================================
! 03. COMPUTE VELOCITY IN INERTIAL FRAME
! ==============================================================================================

! Radial and tangential unit vectors in the inertial frame
i_vec = x_vec*cnu + y_vec*snu
j_vec = y_vec*cnu - x_vec*snu

! Radial and tangential components
v_rad = zeta/(sqrt(z(3))*rho)
v_tan = sqrt((1._dk - z(1)**2_ik - z(2)**2_ik)/(z(3)*rho**2_ik) - 2._dk*Upot)

! Velocity in the inertial frame
v_vec = v_rad*i_vec + v_tan*j_vec

end subroutine DEDROMO2CART


function DEDROMO_TE2TIME(z,phi,flag_time)

! VARIABLES
implicit none
! Arguments IN
real(dk),intent(in)     ::  z(1:8),phi
integer(ik),intent(in)  ::  flag_time
! Function definition
real(dk)                ::  DEDROMO_TE2TIME

! ==============================================================================================

if ( flag_time == 0_ik ) then
    ! Physical time	
    DEDROMO_TE2TIME = z(8)
else if  ( flag_time == 1_ik ) then
    ! Constant time element
    DEDROMO_TE2TIME = z(8) - z(3)**1.5_dk*(z(1)*sin(phi) - z(2)*cos(phi) - phi)
else if  ( flag_time == 2_ik ) then
    ! Linear time element
    DEDROMO_TE2TIME = z(8) - z(3)**1.5_dk*(z(1)*sin(phi) - z(2)*cos(phi))
end if

end function DEDROMO_TE2TIME


function DGDROMO_TE2TIME(z,phi,flag_time)
! Description:
!    Computes value of physical time from the GDromo state vector.
!
! ==============================================================================

! VARIABLES
implicit none
! Arguments IN
real(dk),intent(in)     ::  z(1:8),phi
integer(ik),intent(in)  ::  flag_time       ! Flag indicating the type of time element.
                                            ! = 0: physical time
                                            ! = 1: constant time element
                                            ! = 2: linear time element
! Locals
real(dk)  ::  zeta                          ! Auxiliary variable
! Function definition
real(dk)  ::  DGDROMO_TE2TIME

! ==============================================================================

zeta = z(1)*sinh(phi) + z(2)*cosh(phi)

! Time/Time Element    
if ( flag_time == 0_ik ) then
	! Physical time	
    DGDROMO_TE2TIME = z(8)
else if  ( flag_time == 1_ik ) then
    ! Constant time element
    DGDROMO_TE2TIME = z(8) + z(3)**1.5_dk*(zeta - phi)
else if  ( flag_time == 2_ik ) then
    ! Linear time element
    DGDROMO_TE2TIME = z(8) + z(3)**1.5_dk*zeta
end if

end function DGDROMO_TE2TIME


subroutine DGDROMO2CART(z,phi,pot,rV,vV)
! Description:
!    Computes Cartesian position and velocity from the GDromo state vector "z", 
!    independent variable "phi" and potential "pot".
! 
! ==============================================================================

! VARIABLES
implicit none
! Arguments IN
real(dk),intent(in)   ::  z(1:8),phi,pot
! Arguments OUT
real(dk),intent(out)  ::  rV(1:3),vV(1:3)
! Locals
real(dk)              ::  sph,cph
real(dk)              ::  ecg
real(dk)              ::  gamma,aux1,aux2
real(dk)              ::  rmag,cnu,snu,zeta,rho,emme

! ==============================================================================

! Store hyperbolic function
sph = sinh(phi)
cph = cosh(phi)

! Orbital radius
rho = z(1)*cph + z(2)*sph - 1._dk
rmag = z(3)*rho

! Generalized eccentricity
ecg = sqrt(z(1)**2_ik - z(2)**2_ik)

! Compute auxiliary variables
emme  = sqrt(ecg**2_ik-1._dk)
zeta  = z(1)*sph + z(2)*cph
gamma = sqrt(z(1)**2_ik + z(2)**2_ik)
aux1  = emme*zeta/(rho*ecg*gamma)
aux2  = (emme**2_ik/rho - 1._dk)/(ecg*gamma)

! cosine and sine of the angle of rotation nu between the LVLH and the
! intermediate frames
cnu = z(1)*aux2 - z(2)*aux1
snu = z(1)*aux1 + z(2)*aux2

call DGDROMO2CART_CORE(z,rmag,cnu,snu,zeta,rho,emme,rV,vV,pot)

end subroutine DGDROMO2CART

function DINERT2ORB_GDROMO(vI,z,cnu,snu)
! Description:
!    Transforms a vector vI from orbital to inertial axes through a rotation matrix obtained
!    fromo GDromo(P) elements.
!
! ==============================================================================================
        
! VARIABLES
implicit none

! Arguments
real(dk),intent(in)  ::  vI(1:3)
real(dk),intent(in)  ::  z(:)
real(dk),intent(in)  ::  cnu,snu
real(dk)  ::  DINERT2ORB_GDROMO(1:3)
! Locals
real(dk)  ::  Q0T(1:3,1:3), MnuT(1:3,1:3), QIR(1:3,1:3)

! =============EXECUTION============

! Compute rotation matrix QIR (Inertial -> Orbital)
MnuT = reshape( (/ cnu, -snu, 0._dk,&
                 & snu,  cnu, 0._dk,&
                 & 0._dk, 0._dk, 1._dk /), (/ 3_ik, 3_ik /) )
Q0T   = reshape( (/ 1._dk - 2._dk*z(5)**2_ik - 2._dk*z(6)**2_ik,&
                &   2._dk*z(4)*z(5) - 2._dk*z(6)*z(7),&
                &   2._dk*z(4)*z(6) + 2._dk*z(5)*z(7),&
                &   2._dk*z(4)*z(5) + 2._dk*z(6)*z(7),&
                &   1._dk - 2._dk*z(4)**2_ik - 2._dk*z(6)**2_ik,&
                &   2._dk*z(5)*z(6) - 2._dk*z(4)*z(7),&
                &   2._dk*z(4)*z(6) - 2._dk*z(5)*z(7),&
                &   2._dk*z(5)*z(6) + 2._dk*z(4)*z(7),&
                &   1._dk - 2._dk*z(4)**2_ik - 2._dk*z(5)**2_ik /), &
                & (/ 3_ik, 3_ik /) )
QIR   = matmul(MnuT, Q0T)

! Rotate vI into orbital axes
DINERT2ORB_GDROMO = matmul(QIR,vI)

end function DINERT2ORB_GDROMO

subroutine DGDROMO2CART_CORE(z,rmag,cnu,snu,zeta,rho,emme,RV,VV,pot)
! Description:
!    Transforms from HDromo state vector to Cartesian coordinates. Optionally, it takes into
!    account the potential deriving from Earth's oblateness.
!
! ==============================================================================================
!
!                                  VARIABLES AND DECLARATIONS
!
! ==============================================================================================

implicit none

! Arguments
real(dk),intent(in)   ::  z(1:8)                ! HDromo state vector
real(dk),intent(in)   ::  rmag                  ! Radius vector magnitude
real(dk),intent(in)   ::  pot                   ! Potential
real(dk),intent(in)   ::  zeta,rho,emme         ! Auxiliary quantities
real(dk),intent(in)   ::  cnu,snu               ! cosine and sine of the angle of rotation nu between the LVLH and the
                                                ! intermediate frames
real(dk),intent(out)  ::  RV(1:3), VV(1:3)  	! Cartesian position and velocity
! Locals
real(dk)  ::  u,lam

! ==============================================================================================

! Position vector
RV(1) =  rmag*( ( 1._dk - 2._dk*z(5)**2_ik - 2._dk*z(6)**2_ik )*cnu +&
& 2._dk*( z(4)*z(5) - z(6)*z(7) )*snu ) 
RV(2) =  rmag*( ( 1._dk - 2._dk*z(4)**2_ik - 2._dk*z(6)**2_ik )*snu +&
& 2._dk*( z(4)*z(5) + z(6)*z(7) )*cnu )
RV(3) = 2._dk*rmag*( ( z(5)*z(6) + z(4)*z(7) )*snu +&
& ( z(4)*z(6) - z(5)*z(7) )*cnu )

! Velocity vector
! Radial velocity
u = zeta/(rho*sqrt(z(3)))
! Transverse velocity
lam = sqrt(emme**2_ik/(z(3)*rho**2_ik) - 2._dk*pot)

VV(1) = ( 2._dk*lam*( z(4)*z(5) - z(6)*z(7) ) + u*( 1._dk - 2._dk*z(5)**2_ik - 2._dk*z(6)**2_ik ) )*cnu + &
          	& ( 2._dk*u*( z(4)*z(5) - z(6)*z(7) ) - lam*( 1._dk - 2._dk*z(5)**2_ik - 2._dk*z(6)**2_ik ) )*snu
VV(2) = ( lam*( 1._dk - 2._dk*z(4)**2_ik - 2._dk*z(6)**2_ik ) + 2._dk*u*( z(4)*z(5) + z(6)*z(7) ) )*cnu + &
          	& ( u*( 1._dk - 2._dk*z(4)**2_ik - 2._dk*z(6)**2_ik ) - 2._dk*lam*( z(4)*z(5) + z(6)*z(7) ) )*snu
VV(3) = 2._dk*( lam*( z(5)*z(6) + z(4)*z(7) ) + u*( z(4)*z(6) - z(5)*z(7) ) )*cnu + &
          	& 2._dk*( u*( z(5)*z(6) + z(4)*z(7) ) - lam*( z(4)*z(6) - z(5)*z(7) ) )*snu

end subroutine DGDROMO2CART_CORE


end module TRANSFORM
