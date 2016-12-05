module STEPSIZES
! Description:
!    Contains the subroutine DSTEPSIZE used to compute the initial step size for
!    the distinct formulations according to the reference.
! 
! Reference:
!    Amato D., Baù G., and Bombardelli C., "Accurate numerical orbit propagation
!    of planetary close encounters". Submitted to MNRAS. 2016.
! 
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
! 
! ==============================================================================

use KINDS, only: dk,qk
implicit none


contains


function DSTEPSIZE(R,V,mu,DU,eqs,inSoI,tol)
! Description:
!    Choose the initial stepsize according to the formulation and the sign of the
!    orbital energy.
! 
! Reference:
! [1] Battin R.H., "An Introduction to the Mathematics and Methods of
!     Astrodynamics", Revised Ed., AIAA, Reston, VA, USA. 1999.
! [2] Amato D., Baù G., and Bombardelli C., "Accurate numerical orbit
!    propagation of planetary close encounters". Submitted to MNRAS. 2016.
! 
! ==============================================================================

use CONSTANTS, only: twopi,pi,smaEarth

! VARIABLES
! Arguments
real(qk),intent(in)  ::  R(1:3),V(1:3),mu,DU
real(dk),intent(in)  ::  tol
integer,intent(in)   ::  eqs
logical,intent(in)   ::  inSoI
! Function definition
real(dk)  ::  DSTEPSIZE
! Angular momentum, energy, sma, eccentricity, semi-parameter
real(dk)  ::  rn,vsq
real(dk)  ::  h(1:3),En,sma,ecc,p
! Intervals in true, Gudermannian, hyperbolic and hyperbolic mean anomalies
real(dk)  ::  fplus,zplus,Hplus,Nplus
real(dk)  ::  df,dH,dz,dN
! Elliptic and hyperbolic intervals
real(dk)  ::  rlim
real(dk)  ::  ell,hyp
integer,parameter :: nell = 100, nhyp = 600

! ==============================================================================

! Compute energy and sma
rn = sqrt(dot_product(R,R))
vsq = dot_product(V,V)
En = 0.5_dk*vsq - mu/rn
sma = -0.5_dk*mu/En

if (En < 0._dk) then
  select case (eqs)
    case (-1,1)
      ell = twopi*sqrt(sma**3/mu)
    
    case (2)
      ell = twopi
      
    case (3)
      ell = twopi*sqrt(sma/mu)
    
    case (4)
      ! GDromo cannot integrate closed orbits.
      ell = 0._dk
      
  end select
  
  DSTEPSIZE = ell/nell
  
else
  if (inSoI) then
    ! Set rlim = 0.1 AU
    rlim = 0.1_dk*(smaEarth/DU)
  else
    rlim = 5._dk*(smaEarth/DU)
  end if
  
  ! Limiting true anomaly "fplus" corresponding to the point on the
  ! hyperbola at a distance "rlim" away from the primary, 0 <= fplus <= pi.
  h = [R(2)*V(3) - R(3)*V(2), R(3)*V(1) - R(1)*V(3), R(1)*V(2) - R(2)*V(1)]
  p = dot_product(h,h)/mu
  ecc = sqrt(1._dk - p/sma)
  fplus = acos((p/rlim - 1._dk)/ecc)
  df = 2._dk*fplus
  
  ! Limiting Gudermannian anomaly and interval of Gudermannian anomaly.
  zplus = 2._dk*atan(sqrt((ecc - 1._dk)/(ecc + 1._dk))*tan(0.5_dk*fplus))
  dz = 2._dk*zplus
  
  select case (eqs)
    case (-1,1)
      ! Limiting hyperbolic mean anomaly
      Nplus = ecc*tan(zplus) - log(tan(0.5_dk*zplus + 0.25_dk*pi))
      dN = 2._dk*Nplus
      
      hyp = sqrt((-sma)**3/mu)*dN
    
    case (2)
      hyp = dz
    
    case (3)
      ! EDromo cannot integrate open orbits.
      hyp = 0._dk
    
    case (4)
      ! Limiting hyperbolic anomaly
      Hplus = 2._dk*atanh(tan(0.5_dk*zplus))
      dH = 2._dk*Hplus
      
      hyp = dH
      
  end select
  
  DSTEPSIZE = hyp/nhyp
  
end if

DSTEPSIZE = -DSTEPSIZE/log10(tol)

end function DSTEPSIZE

end module STEPSIZES
