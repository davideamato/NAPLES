module SETTINGS
    
    use KINDS, only: qk,dk
    implicit none
    
    ! SIMULATION SETTINGS
    real(qk)     ::  d_int(1:2)     ! Interval in d [R_Earth]
    real(qk)     ::  e_int(1:2)     ! Interval in e [-]
    real(qk)     ::  th_int(1:2)    ! Interval in theta0 [deg]
    real(qk)     ::  RSw_int(1:2)   ! Interval in switch radii [AU]
    integer  ::  n_d                ! Number of simulations in minimum approach distance
    integer  ::  n_e                ! Number of simulations in hyperbolic eccentricity
    integer  ::  n_th               ! Number of simulations in initial perigee angle
    integer  ::  n_rs               ! Number of simulations in the switch radius
    integer  ::  espacing           ! Spacing of the hyperbolic eccentricity points
    
    ! INTEGRATOR SETTINGS
    integer  ::  integ              ! 1: LSODAR, 2: Radau
    real(qk)    ::  tolref          ! Tolerance for the accurate integration
    real(dk)    ::  tol             ! Tolerance for the test integration
!    integer  ::  steps              ! Number of integration steps
    real(qk) ::  dt_H,dt_CE         ! dt for the H-/H+ and CE phases (days)
    integer  ::  mxsteps            ! Max number of integration steps
    
    ! EQUATIONS OF MOTION
    integer  ::  eqs_H              ! 1: Cowell, 2: K-S, 3: EDromo
    integer  ::  eqs_CE             ! 1: Cowell, 2: K-S, 3: GDromo
    
    ! REGULARIZED FORMULATIONS SETTINGS
    integer  ::  flag_time_EDr      ! 0: physical time, 1: Constant Time Element, 2: Linear Time Element
    integer  ::  flag_time_GDr      ! 0: physical time, 1: Constant Time Element, 2: Linear Time Element
    
    ! OUTPUT PATH
    character(len=4096)  ::  path
    
end module SETTINGS
