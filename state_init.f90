module STATE_INIT

use KINDS, only: ik,dk,qk
implicit none

contains

subroutine DINIT_KS(R0_d,V0_d,t0_d,mu,DU,TU,u0,Vpot_d)
! Description:
!     Starting from dimensional Cartesian coordinates and time, initializes the extended state
!     vector of the Kustaanheimo-Stiefel formulation (from now on "KS"). We try to follow the
!     naming convention used in "Stiefel E. L., Scheifele G. - Linear and Regular Celestial
!     Mechanics, 1974".
!
! ==============================================================================================
!
!                                            PREAMBLE
!
! ==============================================================================================

! VARIABLES
implicit none
! Arguments
real(dk),intent(in)     ::  R0_d(1:3),V0_d(1:3)  ! Dimensional position and velocity [km,km/s]
real(dk),intent(in)     ::  t0_d                 ! Dimensional time [s]
real(qk),intent(in)     ::  mu                   ! Grav parameter of the main body [km^3/s^2]
real(qk),intent(in)     ::  DU,TU                ! Ref. quantities for ND [km,1/s]
real(dk),intent(in)     ::  Vpot_d               ! Dimensional potential [km^2/s^2]
real(dk),intent(out)    ::  u0(:)                ! Extended KS state vector

! Locals
! -- State vectors and associated quantities
real(dk)            ::  x0(1:4)              ! Radius vector in R^4, ND
real(dk)            ::  xdot0(1:4)           ! Velocity vector in R^4, ND
real(dk)            ::  t0                   ! Initial time, ND
real(dk)            ::  r                    ! Radius magnitude, ND
! -- Energies and potentials
real(dk)            ::  h                    ! Total energy, ND
real(dk)            ::  Ksq                  ! Grav parameter of the main body, ND
real(dk)            ::  Kin                  ! Kinetic energy, ND
real(dk)            ::  Vpot                 ! Perturbing potential, ND

! ==============================================================================================

! STATE VECTOR DICTIONARY
! u(1:4)        u1,...,u4; KS-position, in R^4
! u(5:8)        u1',...,u4'; KS-velocity, in R^4
! u(9)          h; (-total energy) = (-Keplerian energy) + (-potential)
! u(10)         t; non-dimensional physical time

! ==============================================================================================
! 01. INITIALIZE POSITION AND VELOCITY IN R^4, NON-DIMENSIONALIZATIONS
! ==============================================================================================

x0(1:3) = R0_d/DU
x0(4)   = 0._dk
r       = sqrt(dot_product(x0,x0))
t0      = t0_d*TU

xdot0(1:3) = V0_d/(DU*TU)
xdot0(4)   = 0._dk

! Also, non-dimensionalize potential and grav parameter
Vpot = Vpot_d/(DU*TU)**2_ik
Ksq  = mu/(DU**3_ik*TU**2_ik)

! ==============================================================================================
! 02. INITIALIZE KS-POSITION AND KS-VELOCITY (state vector u)
! ==============================================================================================

! KS-POSITION
if ( x0(1) >= 0.) then
    u0(1) = 0._dk
    u0(4) = sqrt(.5_dk*(r + x0(1)) - u0(1)**2_ik)
    u0(2) = (x0(2)*u0(1) + x0(3)*u0(4))/(r + x0(1))
    u0(3) = (x0(3)*u0(1) - x0(2)*u0(4))/(r + x0(1))
else
    u0(2) = 0._dk
    u0(3) = sqrt(.5_dk*(r - x0(1)) - u0(2)**2_ik)
    u0(1) = (x0(2)*u0(2) + x0(3)*u0(3))/(r - x0(1))
    u0(4) = (x0(3)*u0(2) - x0(2)*u0(3))/(r - x0(1))
end if

! KS-VELOCITY
u0(5) = .5_dk*dot_product(u0(1:3),xdot0(1:3))
u0(6) = .5_dk*dot_product([-u0(2),u0(1),u0(4)],xdot0(1:3))
u0(7) = .5_dk*dot_product([-u0(3),-u0(4),u0(1)],xdot0(1:3))
u0(8) = .5_dk*dot_product([u0(4),-u0(3),u0(2)],xdot0(1:3))

! Total energy
Kin   = .5_dk*dot_product(xdot0,xdot0)
h     = Ksq/r - Kin - Vpot
u0(9) = h

! Initial time
u0(10) = t0

end subroutine DINIT_KS

function DFTIME_2D(R,V,mu)
! Computes the initial value of the fictitious time as the eccentric
! anomaly (for elliptic orbits) or Gudermannian anomaly (for hyperbolic orbits)
! for the K-S formulation.
! Only works for 2D.

! VARIABLES
implicit none
! Inputs
real(qk)  ::  R(1:3),V(1:3),mu
! Function definition
real(dk)  ::  DFTIME_2D
! Locals
real(dk)  ::  vsq,rv
real(dk)  ::  rN
real(dk)  ::  e_vec(1:3),ecc
real(dk)  ::  lon_per
real(dk)  ::  R_ep(1:2)
real(dk)  ::  f0

! ==============================================================================

! Compute eccentricity and eccentricity vector
vsq = dot_product(V,V)
rv  = dot_product(R,V)
rN = sqrt(dot_product(R,R))
e_vec = ((vsq - mu/rN)*R - rv*V)/mu
ecc = sqrt(dot_product(e_vec,e_vec))

! Compute true anomaly
lon_per = atan2(e_vec(2),e_vec(1))
R_ep = [R(1)*cos(lon_per) + R(2)*sin(lon_per),&
-R(1)*sin(lon_per) + R(2)*cos(lon_per)]/rN
f0   = atan2(R_ep(2),R_ep(1))

! Compute fictitious time
if (ecc <= 1._dk) then
    ! s0 = Eccentric anomaly
    DFTIME_2D = atan2(sqrt(1._dk - ecc**2)*sin(f0),(ecc + cos(f0)))
else
    ! s0 = Gudermannian anomaly
    DFTIME_2D = 2._dk*atan(sqrt((ecc-1._dk)/(ecc+1._dk))*tan(.5_dk*f0))
end if

end function DFTIME_2D

function DHYPAN(R,V)
! Computes hyperbolic anomaly from DIMENSIONLESS R,V; where the reference
! quantities are "DU" for the length and "TU = sqrt(DU**3/mu)" for time.
real(qk)  ::  R(1:3), V(1:3)
real(qk)  ::  vsq,rmag,En,rv
real(dk)  ::  DHYPAN

vsq = dot_product(V,V)
rmag = sqrt(dot_product(R,R))
En  = .5_qk*vsq - 1._qk/rmag
rv  = dot_product(R,V)
DHYPAN = atanh((sqrt(2._qk*En)*rv)/(1._qk + 2._qk*En*rmag))

end function DHYPAN

subroutine DINIT_EDROMO(R0,V0,t0,DU,TU,z0,phi0,W_d,flag_time)
! Description:
!    Initializes the state vector for the EDromo formulation, starting from Cartesian
!    coordinates.
!
! Author:
!    Davide Amato
!    Space Dynamics Group
!    Technical University of Madrid
!    d.amato@upm.es
!    Last revision: 05-09-2015
! ==============================================================================================
!
!                                            PREAMBLE
!
! ==============================================================================================
    
! VARIABLES
implicit none

! Arguments
integer(ik),parameter   ::  neq = 8		        ! Number of elements of EDromo state vector
real(dk),intent(in)     ::  R0(1:3),V0(1:3)		! Dimensional position and velocity [km,km/s]
real(qk),intent(in)     ::  DU,TU               ! Ref. quantities for non-dimensionalization
real(dk),intent(in)     ::  phi0                ! Initial phi value
real(dk),intent(in)     ::  t0               	! Initial time value [s]
real(dk),intent(in)     ::  W_d                 ! Potential [km^2/s^2]
integer(ik),intent(in)  ::  flag_time		    ! = 0 physical time
         								        ! = 1 constant time element
								                ! = 2 linear time element 
real(dk),intent(out)    ::  z0(1:neq)           ! EDromo state vector

! Local variables
real(dk)  ::  y(1:6)           		            ! Cartesian state vector, ND
real(dk)  ::  r0mag,v0mag      	                ! Radius and velocity magnitudes, ND
real(dk)  ::  rv0                               ! r dot v, ND
real(dk)  ::  h0(1:3),hmag0    	                ! Initial angular momentum and magnitude
real(dk)  ::  pot0           		            ! Potential, ND
real(dk)  ::  totEn0           		            ! Initial total energy, ND
real(dk)  ::  c0                                ! Generalized angular momentum
real(dk)  ::  cph,sph
real(dk)  ::  nu0,cnu,snu          	            ! Auxiliary variables
real(dk)  ::  i_vec(1:3),j_vec(1:3),k_vec(1:3)  ! Orbital frame unit vectors in inertial RF
real(dk)  ::  x_vec(1:3),y_vec(1:3)             ! Intermediate frame unit vectors in inertial RF
real(dk)  ::  aux                	            ! Auxiliary variables
real(dk)  ::  zero             		            ! Reference machine zero

! ==============================================================================================
!
!                                            EXECUTION
!
! ==============================================================================================

! ==============================================================================================
! 01. NON-DIMENSIONALIZATION
! ==============================================================================================
y(1:3) = R0/DU
y(4:6) = V0/(DU*TU)
pot0   = W_d/(DU*TU)**2_ik

! Compute machine zero for comparisons
zero = 10._dk*epsilon(0._dk)

! ==============================================================================================
! 02. IN-PLANE ELEMENTS
! ==============================================================================================
r0mag = sqrt(dot_product(y(1:3),y(1:3)))
v0mag = sqrt(dot_product(y(4:6),y(4:6)))

cph = cos(phi0)
sph = sin(phi0)

! Total energy
totEn0  = .5_dk*v0mag**2_ik - 1._dk/r0mag + pot0

! Angular momentum
h0 = DCROSS_PRODUCT(y(1:3),y(4:6))
hmag0 = sqrt(dot_product(h0,h0))

! Generalized angular momentum
c0 = sqrt(hmag0**2_ik + 2._dk*r0mag**2_ik*pot0)

! Dot product
rv0 = dot_product(y(1:3),y(4:6))

! IN-PLANE ELEMENTS
z0(1) = (1._dk + 2._dk*totEn0*r0mag)*cph + rv0*sqrt(-2._dk*totEn0)*sph
z0(2) = (1._dk + 2._dk*totEn0*r0mag)*sph - rv0*sqrt(-2._dk*totEn0)*cph
z0(3) = -1._dk/(2._dk*totEn0)

! ==============================================================================================
! 03. QUATERNION ELEMENTS
! ==============================================================================================

nu0 = phi0 + 2._dk*atan2(rv0,(c0 + r0mag*sqrt(-2._dk*totEn0)))
cnu = cos(nu0); snu = sin(nu0)

! Orbital frame unit vectors in IRF
i_vec = y(1:3)/r0mag
k_vec = h0/hmag0
j_vec = DCROSS_PRODUCT(k_vec,i_vec)

! Intermediate frame unit vectors in IRF
x_vec = i_vec*cnu - j_vec*snu
y_vec = j_vec*cnu + i_vec*snu

! SAFE INITIALIZATION OF THE QUATERNION
! The arguments of the roots have to be always positive, therefore we can safely use ABS()
aux = abs(1._dk + x_vec(1) + y_vec(2) + k_vec(3))
z0(7) = .5_dk*sqrt(aux)
! Check for singularities and NaNs
if ( aux <= zero ) then
    aux    = abs((.5_dk*(k_vec(3) + 1._dk)))
    z0(6)  = sqrt(aux)
    if ( aux <= zero ) then
        aux   = abs(.5_dk*(1._dk - y_vec(2)))
        z0(4) = sqrt(aux)
        if ( aux <= zero ) then
            z0(5) = 1._dk
        else
            z0(5) = y_vec(1)/(2._dk*z0(4))
        end if   
    else
        z0(4) = k_vec(1)/(2._dk*z0(6))
        z0(5) = k_vec(2)/(2._dk*z0(6))
    end if
else
    z0(4) = (y_vec(3) - k_vec(2))/(4._dk*z0(7))
    z0(5) = (k_vec(1) - x_vec(3))/(4._dk*z0(7))
    z0(6) = (x_vec(2) - y_vec(1))/(4._dk*z0(7))
end if

! ==============================================================================
! 04. TIME / TIME ELEMENT
! ==============================================================================

if ( flag_time == 0_ik ) then
    ! Physical time
    z0(8) = t0*TU
elseif  ( flag_time == 1_ik ) then
    ! Constant time element
    z0(8) = t0*TU + z0(3)**1.5_dk*(z0(1)*sph - z0(2)*cph - phi0)
elseif  ( flag_time == 2_ik ) then
    ! Linear time element
    z0(8) = t0*TU + z0(3)**1.5_dk*(z0(1)*sph - z0(2)*cph)
end if

end subroutine DINIT_EDROMO

subroutine DINIT_GDROMO(R0,V0,t0,DU,TU,z0,phi0,W_d,flag_time)
! Description:
!    Initializes the state vector for the EDromo formulation, starting from Cartesian
!    coordinates.
!
! Author:
!    Davide Amato
!    Space Dynamics Group
!    Technical University of Madrid
!    d.amato@upm.es
!    Last revision: 05-09-2015
! ==============================================================================================
!
!                                            PREAMBLE
!
! ==============================================================================================
    
! VARIABLES

implicit none

! Arguments
integer(ik),parameter   ::  neq = 8		        ! Number of elements of HDromo state vector
real(dk),intent(in)     ::  R0(1:3),V0(1:3)		! Dimensional position and velocity [km,km/s]
real(qk),intent(in)     ::  DU,TU               ! Ref. quantities for non-dimensionalization
real(dk),intent(in)     ::  phi0                ! Initial phi value
real(dk),intent(in)     ::  t0               	! Initial time value [s]
real(dk),intent(in)     ::  W_d                 ! Potential [km^2/s^2]
integer(ik),intent(in)  ::  flag_time		    ! = 0 physical time
         								        ! = 1 constant time element
								                ! = 2 linear time element 
real(dk),intent(out)    ::  z0(:)               ! HDromo state vector

! Local variables
real(dk)  ::  y(1:6)           		            ! Cartesian state vector, ND
real(dk)  ::  r0mag,v0mag      	                ! Radius and velocity magnitudes, ND
real(dk)  ::  h0(1:3),hmag0    	                ! Initial angular momentum and magnitude
real(dk)  ::  uu0,u0,lam0,s0      	            ! Initial radial and tangential velocities
real(dk)  ::  pot0           		            ! Potential, ND
real(dk)  ::  totEn0           		            ! Initial total energy, ND
real(dk)  ::  hg0,ecg0                          ! Generalized angular momentum and eccentricity
real(dk)  ::  cnu_0,snu_0          	            ! Auxiliary variables
real(dk)  ::  zeta,gamma          	            ! Auxiliary variables
real(dk)  ::  aux,aux1,aux2      	            ! Auxiliary variables
real(dk)  ::  zero             		            ! Reference machine zero
real(dk)  ::  xO_0(1:3),yO_0(1:3),zO_0(1:3)     ! Unit vectors of the LVLH reference frame (orbital frame)
real(dk)  ::  xI_0(1:3),yI_0(1:3)				! Unit vectors of the intermediate reference frame

! ==============================================================================================
!
!                                            EXECUTION
!
! ==============================================================================================

! ==============================================================================================
! 01. NON-DIMENSIONALIZATION
! ==============================================================================================
y(1:3) = R0/DU
y(4:6) = V0/(DU*TU)
pot0   = W_d/(DU*TU)**2_ik

! Compute machine zero for comparison
zero = 10._dk*epsilon(0._dk)

! ==============================================================================================
! 02. CONVERSION TO HDROMO ELEMENTS
! ==============================================================================================

r0mag = sqrt(dot_product(y(1:3),y(1:3)))
v0mag = sqrt(dot_product(y(4:6),y(4:6)))

! Radial velocity
uu0 = dot_product(y(1:3),y(4:6))
u0 = uu0/r0mag
! Angular momentum
h0 = DCROSS_PRODUCT(y(1:3),y(4:6))
hmag0 = sqrt(dot_product(h0,h0))
! Tangential velocity
lam0 = hmag0/r0mag
    
! Generalized transverse velocity
s0 = sqrt(lam0**2_ik + 2._dk*pot0)
    
! Total energy
totEn0  = v0mag**2_ik/2._dk - 1._dk/r0mag + pot0

z0(1) = (1._dk + 2._dk*totEn0*r0mag)*cosh(phi0) - uu0*sqrt(2._dk*totEn0)*sinh(phi0)
z0(2) = -(1._dk + 2._dk*totEn0*r0mag)*sinh(phi0) + uu0*sqrt(2._dk*totEn0)*cosh(phi0)
z0(3) = 1._dk/(2._dk*totEn0)
    
! Generalized angular momentum
hg0 = r0mag*s0

! Generalized eccentricity
ecg0 = sqrt(z0(1)**2_ik - z0(2)**2_ik)
    
! SAFE INITIALIZATION OF THE QUATERNION
! Auxiliary variables 
zeta = z0(1)*sinh(phi0) + z0(2)*cosh(phi0)
gamma = sqrt(z0(1)**2_ik + z0(2)**2_ik)
aux1 = u0*hg0/(ecg0*gamma)
aux2 = (hg0**2_ik/r0mag - 1._dk)/(ecg0*gamma)
    
! Cosine and sine of the angle of rotation nu between the LVLH and the intermediate frames
cnu_0 = z0(1)*aux2 - z0(2)*aux1
snu_0 = z0(1)*aux1 + z0(2)*aux2
    
! x-y-z unit vectors of the LVLH reference frame (orbital frame)
xO_0 = y(1:3)/r0mag
zO_0 = h0/hmag0
yO_0 = DCROSS_PRODUCT(zO_0,xO_0)
! x-y-z unit vectors of the intermediate reference frame
xI_0 = xO_0*cnu_0 - yO_0*snu_0
yI_0 = xO_0*snu_0 + yO_0*cnu_0
! zI_0 = zO_0
    
! The arguments of the roots have to be always positive, therefore we can safely use ABS()
aux = abs(1._dk + xI_0(1) + yI_0(2) + zO_0(3))
z0(7) = .5_dk*sqrt(aux)
! Check for singularities and nans
if ( aux <= zero ) then
    aux    = abs((.5_dk*(zO_0(3) + 1._dk)))
    z0(6)  = sqrt(aux)
    if ( aux <= zero ) then
        aux   = abs((.5_dk*(1._dk-yI_0(2))))
        z0(4) = SQRT(aux)
        if ( aux <= zero ) then
            z0(5) = 1._dk
        else
            z0(5) = 1._dk/(2._dk*z0(4)) * yI_0(1)
        end if   
    else
        z0(4) = 1._dk/(2._dk*z0(6)) * zO_0(1)
        z0(5) = 1._dk/(2._dk*z0(6)) * zO_0(2)
    end if
else
    z0(4) = 1._dk/(4._dk*z0(7))*( yI_0(3) - zO_0(2) )
    z0(5) = 1._dk/(4._dk*z0(7))*( zO_0(1) - xI_0(3) )
    z0(6) = 1._dk/(4._dk*z0(7))*( xI_0(2) - yI_0(1) )
end if

! Time/Time Element    
if ( flag_time == 0_ik ) then
	! Physical time	
    z0(8) = t0*TU
elseif  ( flag_time == 1_ik ) then
    ! Constant time element
    z0(8) = t0*TU - z0(3)**1.5_dk*(zeta - phi0)
elseif  ( flag_time == 2 ) then
    ! Linear time element
    z0(8) = t0*TU - z0(3)**1.5_dk*zeta
end if

end subroutine DINIT_GDROMO


function DCROSS_PRODUCT(a,b)
! Description:
!    Computes cross product of a*b.
!
! ==============================================================================

! VARIABLES
implicit none
! Function definition
real(dk)        :: DCROSS_PRODUCT(1:3)
! Arguments IN
real(dk)        :: a(1:3),b(1:3)

! ==============================================================================

DCROSS_PRODUCT(1) = a(2)*b(3) - a(3)*b(2)
DCROSS_PRODUCT(2) = a(3)*b(1) - a(1)*b(3)
DCROSS_PRODUCT(3) = a(1)*b(2) - a(2)*b(1)

end function DCROSS_PRODUCT


end module
