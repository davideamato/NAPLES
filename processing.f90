module PROCESSING

use KINDS, only: dk,qk
implicit none

contains


subroutine RESIDUALS(npts,ref,test,AU_km,mu,RSw,dR_abs,dR_rel,dV_abs,dV_rel,dEn_rel,&
&dSMA_abs,RSw_eff,RSw_diff_max)

! VARIABLES
implicit none
! Arguments IN
integer,intent(in)    ::  npts              ! Number of points to evaluate
real(qk),intent(in)   ::  ref(1:npts,1:7)   ! Heliocentric reference state [s,km,km/s]
real(qk),intent(in)   ::  test(1:npts,1:7)  ! Heliocentric test state      [s,km,km/s]
real(qk),intent(in)   ::  AU_km,mu,RSw
real(dk),intent(in)   ::  RSw_eff(1:2)
! Arguments OUT
real(dk),intent(out)  ::  dR_abs(1:npts),dR_rel(1:npts)     ! Heliocentric position residuals [AU,-]
real(dk),intent(out)  ::  dV_abs(1:npts),dV_rel(1:npts)     ! Heliocentric velocity residuals [km/s,-]
real(dk),intent(out)  ::  dEn_rel(1:npts)                   ! Energy residuals [(km/s)^2,-]
real(dk),intent(out)  ::  dSMA_abs(1:npts)                  ! SMA residuals [AU]
real(dk),intent(out)  ::  RSw_diff_max

! Locals
! QUAD
real(qk)     ::  kms_AUPerDay
real(qk)     ::  dR_vec(1:3),dV_vec(1:3)
real(qk)     ::  r_ref,v_ref
! DOUBLE
real(dk)     ::  r_test,v_test
real(dk)     ::  En_ref,En_test
real(dk)     ::  sma_ref,sma_test
real(dk)     ::  dEn_abs
real(dk)     ::  RSw_diff(1:2)

integer  ::  i

! ==============================================================================================

do i=1,npts
  r_ref = sqrt(dot_product(ref(i,2:4),ref(i,2:4)))
  v_ref = sqrt(dot_product(ref(i,5:7),ref(i,5:7)))
  r_test = sqrt(dot_product(test(i,2:4),test(i,2:4)))
  v_test = sqrt(dot_product(test(i,5:7),test(i,5:7)))
  
  ! 1. POSITION RESIDUALS
  dR_vec = ref(i,2:4)/AU_km - test(i,2:4)/AU_km
  dR_abs(i) = sqrt(dot_product(dR_vec,dR_vec))
  dR_rel(i) = dR_abs(i)/(r_ref/AU_km)
  
  ! 2. VELOCITY RESIDUALS
  dV_vec = ref(i,5:7) - test(i,5:7)
  dV_abs(i) = sqrt(dot_product(dV_vec,dV_vec))
  dV_rel(i) = dV_abs(i)/v_ref
  
  ! 3. ENERGY RESIDUALS
  En_ref = 0.5_qk*v_ref**2 - mu/r_ref
  En_test = 0.5_dk*v_test**2 - mu/r_test
  dEn_abs = abs(En_test - En_ref)
  dEn_rel(i) = dEn_abs/En_ref
  
  ! 4. SMA RESIDUALS
  sma_ref = -mu/(v_ref**2 - 2._qk*mu/r_ref)  ! Note: we compute it in quad
  sma_test = -0.5_dk*mu/En_test
  dSMA_abs(i) = (sma_test - sma_ref)/AU_km
  
end do

! Check on the effective switch radius
RSw_diff = abs((RSw_eff - RSw)/RSw)
RSw_diff_max = maxval(RSw_diff)

end subroutine RESIDUALS

end module PROCESSING
