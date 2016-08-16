program NAPLES
! Description: program NAPLES (Numerical Analysis of Planetary EncounterS).
! Propagates orbits in the Earth-Sun circular, restricted three-body problem 
! starting from initial conditions assigned at the minimum approach distance. 
! The initial conditions are parametrized in (d,e,theta), and a change of 
! primaries is performed when passing the geocentric distance RSw, according to 
! the Online Trajectory Matching algorithm.
! 
! Davide Amato
! d.amato@upm.es

! MODULES
use KINDS,       only: ik,dk,qk
use SETTINGS,    only: d_int,e_int,th_int,RSW_int,n_d,n_e,n_th,n_rs,espacing,&
&tolref,tol,eqs_H,eqs_CE,integ,path
use CONSTANTS,   only: REarth,wEarth,muEarth,muSun,R_SoI,smaEarth,d2r,pi,&
&secsPerDay,twopi
use AUXILIARIES, only: JD_i,JD_stop,JD_CA,time_steps,meth_switch,RSwitch,&
&ind_time,inSoI
use INTEGRATION, only: er_avg,eval_ratio,called
use IO,          only: id_set,id_ares,id_rres,id_enres,id_stats
use IO,          only: READ_SETTINGS,WRITE_SETTINGS,CREATE_ENRES,CREATE_ABSRES,&
&CREATE_RELRES,CREATE_SMARES,CREATE_STATS,WRITE_RESULTS,WRITE_EXCEPTION
use REFERENCE_TRAJECTORY_MOD, only: REFERENCE_TRAJECTORY_COWELL
use TEST_PROP,   only: DP_TRAJECTORY
use PROCESSING,  only: RESIDUALS
use THIRDBODIES, only: QPOS_VEL_CIRC


! VARIABLES
implicit none
! Parameter meshing
integer               ::  i_d,j_e,k_t,l_r,i_sim ! Counters
real(qk),allocatable  ::  d_arr(:)              ! [R_Earth]
real(qk),allocatable  ::  e_arr(:)              ! [-]
real(qk),allocatable  ::  th_arr(:)             ! [deg]
real(qk),allocatable  ::  RSw_arr(:)            ! [AU]
real(qk)  ::  dmin,ecc,theta                    ! [km,-,rad]
real(qk)  ::  delta
! Initial conditions
real(qk)  ::  R_CA_pl(1:3),V_CA_pl(1:3)  ! Geo position and velocity @ close approach (t = 0) [km,km/s]
real(qk)  ::  R_CA(1:3),V_CA(1:3)        ! Helio position and velocity @ close approach (t = 0) [km,km/s]
real(qk)  ::  TEarth
! Trajectory data
real(qk),allocatable  ::  traj_bwd(:,:),traj_fwd(:,:),ref_traj(:,:)
real(dk),allocatable  ::  dp_H1(:,:), dp_CE(:,:), dp_H2(:,:), dp_full(:,:)
real(dk)  ::  state_end(1:7)
real(qk)  ::  state_cp(1:4,1:7)          ! State at checkpoints (test)
real(qk)  ::  state_cp_ref(1:4,1:7)      ! State at checkpoints (reference)
real(qk)  ::  yD_Earth(1:6)
real(qk)  ::  R_i(1:3), V_i(1:3)
integer   ::  lenref, len_H1, len_CE, len_H2
! Times
integer   ::  ind_in_bwd,ind_out_fwd
integer   ::  ind_in,ind_out
integer   ::  ind_in_ref,ind_out_ref
real(qk)  ::  JD_in,JD_out,JD_f
! Results
real(dk)  ::  dR_abs(1:3),dR_rel(1:3),dV_abs(1:3),dV_rel(1:3)
real(dk)  ::  dEn_rel(1:3),dSMA_abs(1:3)
! Sanity checks, settings and diagnostics
integer   ::  istate_bwd,istate_fwd,istate_curr
integer   ::  idiag(1:10)
real(dk)  ::  rdiag(1:3)
integer   ::  totsim,step_mess,i_fail
integer   ::  eqs
real(dk)  ::  RSw_eff(1:2),RSw_diff
integer   ::  XC,MF(1:3)
integer   ::  callsArr(1:2),calls_end,istepsArr(1:2),isteps_end
real(dk)  ::  mordArr(1:2),mord_end
integer   ::  ncheck
logical   ::  sanity
! CPU time
real(dk)  ::  cputime_start,cputime_end
! Auxiliary variables
real(qk)  ::  yt_quad(1:7)

! Output
character(len=8)             ::  date,date_end
character(len=10)            ::  time,time_end
character(len=14)            ::  runID
character(len=4096)          ::  ppath

! Debug
integer  ::  ii

write(*,'(80(''*''))')
write(*,'(''*'',35x,a7,36x,''*'')') 'NAPLES'
write(*,'(''*'',17x,a43,18x,''*'')') 'Numerical Analysis of PLanetary EncounterS.'
write(*,'(80(''*''))')

! ==============================================================================
! 01. READ INPUTS FROM FILE
! ==============================================================================

call READ_SETTINGS()
write(*,'(/,a)') 'Press ENTER to continue... '
read(*,*)

call DATE_AND_TIME(date,time)
runID = date//time(1:6)
write(*,'(4a,/)') 'Simulation ',date,time(1:6),' starting.'

! ==============================================================================
! 02. GENERATE (d,e,th) MESH
! ==============================================================================

allocate(d_arr(1:n_d))
allocate(e_arr(1:n_e))
allocate(th_arr(1:n_th))
allocate(RSw_arr(1:maxval([n_rs,1])))

! d mesh (linear)
delta = (d_int(2) - d_int(1))/(n_d - 1)
d_arr(1) = d_int(1)
do i_d = 2,n_d
    d_arr(i_d) = d_arr(1) + (i_d - 1)*delta
end do

! e mesh (linear or log)
select case (espacing)
case(1)
    delta = (e_int(2) - e_int(1))/(n_e - 1)
    e_arr(1) = e_int(1)
    do j_e = 2,n_e
        e_arr(j_e) = e_arr(1) + (j_e - 1)*delta
    end do
case(2)
    delta = (log10(e_int(2)) - log10(e_int(1)))/(n_e - 1)
    e_arr(1) = e_int(1)
    do j_e = 2,n_e
        e_arr(j_e) = 10._qk**(log10(e_arr(1)) + (j_e - 1)*delta)
    end do
end select

! theta mesh (linear)
delta = (th_int(2) - th_int(1))/(n_th - 1)
th_arr(1) = th_int(1)
do k_t = 2,n_th
    th_arr(k_t) = th_arr(1) + (k_t - 1)*delta
end do

! RSw mesh (linear)
delta = (RSw_int(2) - RSw_int(1))/(n_rs - 1)
RSw_arr(1) = RSW_int(1)
do l_r = 2,n_rs
    RSw_arr(l_r) = RSw_arr(1) + (l_r - 1)*delta
end do

! ==============================================================================
! 03. OPEN OUTPUT FILES
! ==============================================================================

ppath = trim(path)//runID//'/'
call SYSTEM('mkdir '//trim(ppath))
call WRITE_SETTINGS(trim(ppath)//'settings_'//runID//'.txt',runID)
call CREATE_ENRES(trim(ppath)//'enres_'//runID//'.dat',runID)
call CREATE_ABSRES(trim(ppath)//'absres_'//runID//'.dat',runID)
call CREATE_RELRES(trim(ppath)//'relres_'//runID//'.dat',runID)
call CREATE_SMARES(trim(ppath)//'smares_'//runID//'.dat',runID)
call CREATE_STATS(trim(ppath)//'stats_'//runID//'.dat',runID)

! ==============================================================================
! 04. LOOP OVER THE PARAMETERS
! ==============================================================================

! Further misc initializations
wEarth = sqrt((muEarth+muSun)/smaEarth**3)
R_SoI  = (muEarth/muSun)**(2._qk/5._qk)*smaEarth
d2r = pi/180._qk
i_sim = 0; i_fail = 0
totsim = n_d*n_e*n_th*maxval([n_rs,1])
er_avg = 0; eval_ratio = 0
if (totsim <= 10000) then
  step_mess = 10
else if (totsim > 10000 .and. totsim < 1000000) then
  step_mess = 100
else
  step_mess = 1000
end if
lenref = 0;
len_H1 = 0; len_CE = 0; len_H2 = 0 

d_loop: do i_d = 1,n_d
  
  dmin = d_arr(i_d)*REarth
  
  e_loop: do j_e = 1,n_e
    
    ecc = e_arr(j_e)
    
    th_loop: do k_t = 1,n_th
      
      theta = th_arr(k_t)*d2r
      
      RSw_loop: do l_r = 1,maxval([n_rs,1])
        
        RSwitch = Rsw_arr(l_r)*smaEarth
        
        i_sim = i_sim + 1
        
        ! OUTPUT TO USER
        if (mod(i_sim,step_mess)==0) then
          write(*,'(a,f6.2,3(a,i6),2(a,f6.2))') 'Progress: ',real(i_sim)/real(totsim)*100.,&
          &'%; completed ',i_sim,' out of ',totsim,'; failed ',i_fail,' ( ',&
          &real(i_fail)/real(i_sim)*100.,'%)'
        end if
        
        ! ======================================================================
        ! 05. INITIALIZE HELIOCENTRIC STATE
        ! ======================================================================
        
        ! Geocentric state
        R_CA_pl = dmin*[cos(theta),sin(theta),0._qk]
        V_CA_pl = sqrt(muEarth/sqrt(dot_product(R_CA_Pl,R_CA_Pl))*(1._qk + ecc))*&
        [-sin(theta),cos(theta),0._qk]
        
        ! Heliocentric state
        R_CA = [smaEarth,0._qk,0._qk] + R_CA_pl
        V_CA = [0._qk,sqrt((muSun+muEarth)/smaEarth),0._qk] + V_CA_pl
        
        ! Propagation interval: 6 months before and after close approach.
        JD_CA = 0._qk/secsPerDay
        TEarth = twopi*sqrt(smaEarth**3/(muSun+muEarth))
        JD_i  = (JD_CA - 0.5_qk*TEarth)/secsPerDay
        JD_f  = (JD_CA + 0.5_qk*TEarth)/secsPerDay
        
        ! Integration diagnostics
        istate_bwd = 1_ik; istate_fwd = 1_ik
        meth_switch = .false.
        MF = 0
        XC = 0
        calls_end = 0; callsArr = 0
        isteps_end = 0; istepsArr = 0
        mord_end = 0._dk; mordArr = 0
        RSw_eff = 0._dk
        
        ! ======================================================================
        ! 05a. REFERENCE TRAJECTORY - BWD
        ! ======================================================================
        
        call CPU_TIME(cputime_start)
        
        call REFERENCE_TRAJECTORY_COWELL(R_CA_pl,V_CA_pl,JD_CA,JD_i,tolref,traj_bwd,&
        &JD_in,ind_in_bwd,istate_bwd)
        if (EXCEPTION_CHECK(istate_bwd,phase=0)) cycle RSw_loop
        
        ! ======================================================================
        ! 05b. REFERENCE TRAJECTORY - FWD
        ! ======================================================================
        
        call REFERENCE_TRAJECTORY_COWELL(R_CA_pl,V_CA_pl,JD_CA,JD_f,tolref,&
        &traj_fwd,JD_out,ind_out_fwd,istate_fwd)
        if (EXCEPTION_CHECK(istate_fwd,phase=0)) cycle RSw_loop
        
        call CPU_TIME(cputime_end)
        
        ! Save reference trajectory
        lenref = size(traj_bwd,1) - 1_ik + size(traj_fwd,1)
        allocate(ref_traj(lenref,1:7)); ref_traj = 0._qk
        ref_traj(1:size(traj_bwd,1) - 1,:) = traj_bwd(size(traj_bwd,1):2:-1,:)
        ref_traj(size(traj_bwd,1):lenref,:) = traj_fwd
        
        ! Save indices of times corresponding to switches (useful in the case in
        ! wich the reference and test trajectories are output at different times)
        ind_in_ref = size(traj_bwd,1) - ind_in_bwd + 1
        ind_out_ref = size(traj_bwd,1) + ind_out_fwd - 1 
        
        ! Save time steps vector (excluding the switch times).
        ! This is for output by LSODAR in the test propagation.
        if (allocated(time_steps)) deallocate(time_steps)
        allocate(time_steps(1:lenref),source=ref_traj(:,1))

        ! ======================================================================
        ! 06. TEST PROPAGATION
        ! ======================================================================
        
        ! === START OF PHASE H1 (heliocentric) ===
        
        ! Initialize sphere of influence flag and set the first output time
        inSoI = .false.
        ind_time = 2
        
        state_cp(1,:) = ref_traj(1,:)
        if (n_rs /= 0) then
          JD_stop = JD_in
        else if (n_rs == 0) then
          JD_stop = JD_f
        end if
        eqs = eqs_H
        
        call DP_TRAJECTORY(state_cp(1,2:4),state_cp(1,5:7),JD_i,JD_stop,eqs,&
        &integ,tol,dp_H1,idiag,rdiag)
        if (EXCEPTION_CHECK(idiag(2),phase=1)) cycle RSw_loop
        
        len_H1 = size(dp_H1,1)
        if (n_rs == 0) then
          state_cp(4,:) = dp_H1(len_H1,:)
          
          ! Save diagnostics
          calls_end = calls_end + idiag(1)
          isteps_end = isteps_end + idiag(3)
          mord_end = mord_end + rdiag(2)
        
        else
          state_cp(2,:) = dp_H1(len_H1,:)
          
          ! Save diagnostics
          calls_end = calls_end + idiag(1);   callsArr(1) = idiag(1)
          isteps_end = isteps_end + idiag(3); istepsArr(1) = idiag(3)
          mord_end = mord_end + rdiag(2);     mordArr(1) = rdiag(2)
        
          ! === START OF PHASE CE (planetocentric) ===
          
          inSoI = .true.
          
          ! Switch to the planetocentric frame in quad precision.
          yD_Earth = QPOS_VEL_CIRC(real(dp_H1(len_H1,1),qk),wEarth,smaEarth)
          R_i = state_cp(2,2:4) - yD_Earth(1:3)
          V_i = state_cp(2,5:7) - yD_Earth(4:6)
          
          JD_stop = JD_out
          eqs = eqs_CE
          
          call DP_TRAJECTORY(R_i,V_i,state_cp(2,1)/secsPerDay,JD_stop,eqs,&
          &integ,tol,dp_CE,idiag,rdiag)
          if (EXCEPTION_CHECK(idiag(2),phase=2)) cycle RSw_loop

          ! Save diagnostics
          calls_end = calls_end + idiag(1);   callsArr(2) = idiag(1)
          isteps_end = isteps_end + idiag(3); istepsArr(2) = idiag(3)
          mord_end = mord_end + rdiag(2);     mordArr(2) = rdiag(2)
          
          len_CE = size(dp_CE,1)
          ! Save effective switch radii
          RSw_eff(1) = sqrt(dot_product(R_i,R_i))
          RSw_eff(2) = sqrt(dot_product(dp_CE(len_CE,2:4),dp_CE(len_CE,2:4)))
          RSw_eff = RSw_eff/smaEarth
          
          ! dp_CE is planetocentric, convert to heliocentric in quad precision.
          helio_conv: do ii = 1,len_CE
            yt_quad = dp_CE(ii,1:7)
            dp_CE(ii,2:7) = yt_quad(2:7) + QPOS_VEL_CIRC(yt_quad(1),wEarth,smaEarth)
          end do helio_conv
          
          state_cp(3,:) = dp_CE(len_CE,:)
          
          ! === START OF PHASE H2 (heliocentric) ===
          
          inSoI = .false.
          R_i = state_cp(3,2:4)
          V_i = state_cp(3,5:7)
          
          JD_stop = JD_f
          eqs = eqs_H
          
          call DP_TRAJECTORY(R_i,V_i,state_cp(3,1)/secsPerDay,JD_stop,eqs,&
          &integ,tol,dp_H2,idiag,rdiag)
          if (EXCEPTION_CHECK(idiag(2),phase=3)) cycle RSw_loop

          ! Save diagnostics
          calls_end = calls_end + idiag(1)
          isteps_end = isteps_end + idiag(3)
          mord_end = mord_end + rdiag(2); mord_end = mord_end/3._dk
          
          len_H2 = size(dp_H2,1)
          state_cp(4,:) = dp_H2(len_H2,:)
          
          ! === END OF THE PROPAGATION ===
          
        end if
        
        ! ======================================================================
        ! 07. PROCESS AND SAVE OUTPUT
        ! ======================================================================
        
        ! If n_rs == 0 (HELIO mode) then ind_in = len_H1 = size(dp_H1,1).
        ind_in = len_H1
        ind_out = len_H1 + len_CE - 1
        
        ! Merge the trajectories together in dp_full(:,:)
        ! Temporary fix: only if integ == 1
        if (allocated(dp_full)) deallocate(dp_full)
        allocate(dp_full(1:lenref,1:7))
        ! Sanity check on the length of the trajectories for LSODAR
        if (integ == 1) then
          if ( (len_H1 + len_CE + len_H2 - 2 /= lenref) .and. n_rs /= 0) then
            sanity = EXCEPTION_CHECK(-13,3)
            cycle RSw_loop
          end if
          dp_full(1:ind_in,:) = dp_H1(1:len_H1,:)
          if (n_rs /= 0) then
            dp_full(ind_in+1:ind_out,:) = dp_CE(2:len_CE,:)
            dp_full(ind_out+1:lenref,:) = dp_H2(2:len_H2,:)
          end if
          
        end if
        
        ! Save checkpoint values of the reference trajectories
        state_cp_ref(1,:) = ref_traj(1,:)
        state_cp_ref(2,:) = ref_traj(ind_in_ref,:)
        state_cp_ref(3,:) = ref_traj(ind_out_ref,:)
        state_cp_ref(4,:) = ref_traj(lenref,:)
        ncheck = 3
        
        call RESIDUALS(ncheck,state_cp_ref(2:4,:),state_cp(2:4,:),smaEarth,muSun,&
        &RSw_arr(l_r),dR_abs,dR_rel,dV_abs,dV_rel,dEn_rel,dSMA_abs,RSw_eff,RSw_diff)
        
        ! OUTPUT RESIDUALS
        call WRITE_RESULTS(i_d,j_e,k_t,l_r,d_arr(i_d),e_arr(j_e),th_arr(k_t),&
        &RSw_arr(l_r),dR_abs,dR_rel,dV_abs,dV_rel,dEn_rel,dSMA_abs,RSw_diff,MF,&
        &XC,callsArr,istepsArr,mordArr,calls_end,mord_end,isteps_end)
        
        ! Flush output to disk every 100 propagations
        call FLUSHO()
        
        ! Deallocation for next iteration
        call CLEAN()
        
      end do RSw_loop
    end do th_loop
  end do e_loop
end do d_loop

! ==============================================================================
! 08. FINALIZATION
! ==============================================================================

! Deallocations
deallocate(d_arr)
deallocate(e_arr)
deallocate(th_arr)
deallocate(RSw_arr)
call CLEAN()
close(id_set); close(id_ares); close(id_rres); close(id_enres); close(id_stats)

call DATE_AND_TIME(date_end,time_end)
write(*,'(/,a)') 'Simulation '//runID//&
&' completed on '//date_end(1:4)//'/'//date_end(5:6)//'/'//date_end(7:8)//&
&' at '//time_end(1:2)//':'//time_end(3:4)//':'//time_end(5:6)//'.'

if (integ == 2) then
  ! Display the average ratio of function evaluations
  er_avg = er_avg/called
  write(*,'(a,f7.1,a)') 'Root-finding evaluations / stepping evaluations = ',&
  &er_avg*100._dk,'% (avg.)'

end if

contains


subroutine CLEAN()
    
  implicit none
  if (allocated(ref_traj)) deallocate(ref_traj)
  if (allocated(traj_bwd)) deallocate(traj_bwd)
  if (allocated(traj_fwd)) deallocate(traj_fwd)
  if (allocated(time_steps)) deallocate(time_steps)
  if (allocated(dp_H1)) deallocate(dp_H1)
  if (allocated(dp_CE)) deallocate(dp_CE)
  if (allocated(dp_H2)) deallocate(dp_H2)
  if (allocated(dp_full)) deallocate(dp_full)

end subroutine


function EXCEPTION_CHECK(istate,phase)

integer,intent(in)  ::  istate,phase
logical  ::  EXCEPTION_CHECK

  EXCEPTION_CHECK = .false.
  if (meth_switch) MF(phase) = 1
  if (istate < 0) then
    XC = phase*100 - istate
    call WRITE_EXCEPTION(i_d,j_e,k_t,l_r,&
    &d_arr(i_d),e_arr(j_e),th_arr(k_t),Rsw_arr(l_r),MF,XC)
    call CLEAN()
    call FLUSHO()
    i_fail = i_fail + 1
    EXCEPTION_CHECK = .true.
  end if

end function EXCEPTION_CHECK


subroutine FLUSHO()
! Close and reopen the output files every 100th propagation
! to flush buffer to disk.
  
  if (mod(i_sim,100)==0) then
    close(id_set); close(id_ares); close(id_rres); close(id_enres)
    close(id_stats)
    
    open(unit=id_set,file=trim(ppath)//'settings_'//runID//'.txt',&
    &status='old',action='write',position='append')
    open(unit=id_enres,file=trim(ppath)//'enres_'//runID//'.dat',&
    &status='old',action='write',position='append')
    open(unit=id_ares,file=trim(ppath)//'absres_'//runID//'.dat',&
    &status='old',action='write',position='append')
    open(unit=id_rres,file=trim(ppath)//'relres_'//runID//'.dat',&
    &status='old',action='write',position='append')
    open(unit=id_stats,file=trim(ppath)//'stats_'//runID//'.dat',&
    &status='old',action='write',position='append')
  end if

end subroutine FLUSHO

end program NAPLES

! ========================================
! ========== DEBUGGING GARBAGE ===========
! ========================================

!do ii=1,lenref
!  write(52,'(7(es22.15,'',''))') dp_traj(ii,1:7)
!end do

!!subroutine progress(i,nsim)

!!implicit none
!!integer,intent(in)  ::  i,nsim
!!character(len=80)   ::  bar
!!integer  ::  barlen

!!bar = 'Progress:    % '//'[===============================================================]'
!!if (real(i/nsim,dk) <= 1._dk) then
!!  barlen = (i/nsim)*60
!!else
!!  barlen = 60
!!end if

!end subroutine progress
!      write(*,*) 'X0 = ',ref_traj(1,1)
!      write(*,*) 'Y0 = ',ref_traj(1,2)
!      write(*,*) 'Z0 = ',ref_traj(1,3)
!      write(*,*) 'Xf = ',ref_traj(lenref,1)
!      write(*,*) 'Yf = ',ref_traj(lenref,2)
!      write(*,*) 'Zf = ',ref_traj(lenref,3)
!      
!!! DEBUG
!write(*,*) size(traj_bwd,1), size(traj_fwd,1)
!write(*,*) lenref
!do ii=1,size(traj_bwd,1)
!  write(48,'(7(es22.15,'',''))') traj_bwd(ii,:)
!end do
!do ii=1,size(traj_fwd,1)
!  write(49,'(7(es22.15,'',''))') traj_fwd(ii,:)
!end do
!do ii=1,lenref
!  write(50,'(7(es22.15,'',''))') ref_traj(ii,:)
!end do
!!!! DEBUG
! Save time step sizes used in the reference integration.
        ! dt_H1 = t2 - t1
        ! dt_CE = t_(CA + 1) - t_CA
        ! dt_H2 = tf - t(f-1)
!        dt_H1 = time_steps(2) - time_steps(1)
!        dt_CE = time_steps(size(traj_bwd) + 1) - time_steps(size(traj_bwd))
!        dt_H2 = time_steps(lenref) - time_steps(lenref-1)

!        !! DEBUG
!        do ii=1,lenref
!          write(20,'(7(es22.15,'',''))') ref_traj(ii,:)
!          write(30,'(7(es22.15,'',''))') dp_full(ii,:)
!        end do
!        do ii=1,len_H1
!          write(21,'(7(es22.15,'',''))') dp_H1(ii,:)
!        end do
!        do ii=1,len_CE
!          write(22,'(7(es22.15,'',''))') dp_CE(ii,:)
!        end do
!        do ii=1,len_H2
!          write(23,'(7(es22.15,'',''))') dp_H2(ii,:)
!        end do
!        stop
!        !! DEBUG
!        write(*,'(7(es22.15))') ( state_cp(ii,:) , ii = 1,4 )
!        write(*,'(7(es22.15))') ( state_cp_ref(ii,:) , ii = 1,4 )
!        write(*,*) d_arr(i_d)
!        write(*,*) e_arr(j_e)
!        write(*,*) th_arr(k_t)
!        write(*,*) RSw_arr(l_r)
!        do ii=1,lenref
!          write(20,'(7(es22.15,'',''))') ref_traj(ii,:)
!        end do
!        write(30,'(7(es22.15,1x))') (state_cp(ii,:), ii=1,4)
