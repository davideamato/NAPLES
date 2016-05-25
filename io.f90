module IO

! MODULE
use KINDS, only: dk,qk

! VARIABLES
implicit none
! Input file unit
integer,parameter  ::  iin = 10, id_set = 11
integer,parameter  ::  id_ares = 12, id_rres = 13, id_enres = 14
integer,parameter  ::  id_stats = 15


contains

subroutine READ_SETTINGS()

! MODULES
use SETTINGS

! VARIABLES
implicit none
character(len=13)   ::  title

! ==============================================================================

write(*,'(/,a)') 'Reading input file...'

open(unit=iin,file='input.txt',status='old',action='read')

write(*,'(a,/)') 'Done.'

write(*,'(a)') 'SIMULATION SETTINGS:'

call SKIP(10)

! Read parameter space intervals
read(iin,'(a13,2(e23.15))') title, d_int
read(iin,'(a13,2(e23.15))') title, e_int
read(iin,'(a13,2(e23.15))') title, th_int
read(iin,'(a13,2(e23.15))') title, RSw_int
write(*,'(a16,g14.7,'', '',g14.7,a)') 'd interval: ', d_int, ' Earth radii'
write(*,'(a16,g14.7,'', '',g14.7)')   'e interval: ', e_int
write(*,'(a16,g14.7,'', '',g14.7,a)') 'theta interval: ', th_int, ' deg'
write(*,'(a16,g14.7,'', '',g14.7,a)') 'R_sw interval: ', RSW_int, ' AU'

call SKIP(7)

! Read mesh settings
read(iin,'(a13,i4)') title, n_d
read(iin,'(a13,i4)') title, n_e
read(iin,'(a13,i4)') title, n_th
read(iin,'(a13,i4)') title, n_rs
read(iin,'(a13,i4)') title, espacing

write(*,'(a16,i4)') 'd points: ',n_d
write(*,'(a16,i4)') 'e points: ', n_e
write(*,'(a16,i4)') 'theta points: ', n_th
write(*,'(a16,i4)') 'R_Sw points: ', n_rs
write(*,'(a16,i10,/)') 'Total: ',n_d*n_e*n_th*n_rs
select case(espacing)
case(1)
    write(*,'(a16,a)') 'e spacing: ','  Linear'
case(2)
    write(*,'(a16,a)') 'e spacing: ','  Log'
end select

call SKIP(8)

! Read integrator settings
read(iin,'(a13,i2)') title, integ
read(iin,'(a13,e22.15)') title, tolref
read(iin,'(a13,e22.15)') title, tol
read(iin,'(a13,e22.15)') title, dt_H
read(iin,'(a13,e22.15)') title, dt_CE
read(iin,'(a13,i9)')     title, mxsteps

select case (integ)
case(1)
  write(*,'(a16,a)') 'Integrator: ','  LSODAR'
case(2)
  write(*,'(a16,a)') 'Integrator ','  Radau'
end select
write(*,'(a16,g14.7)') 'Ref. tol: ', tolref
write(*,'(a16,g14.7)') 'Test tol: ', tol
write(*,'(a16,g14.7,a)') 'dt_H: ',dt_H,' days'
write(*,'(a16,g14.7,a)') 'dt_CE: ',dt_CE,' days'
write(*,'(a16,i9)')    'Max steps: ',mxsteps

call SKIP(5)

! Read choice of equations of motion (heliocentric)
read(iin,'(a13,i2)') title, eqs_H
select case(eqs_H)
case default
   write(*,'(a,i2)') 'eqs_H = ',eqs_H
   write(*,'(a)') 'Invalid choice of the equations of motion for the heliocentric phase.'
   write(*,'(a)') 'Stopping...'
   stop
case(-1)
   write(*,'(a16,a)') 'Equations (H): ', '  Cowell (II)'
case(1)
   write(*,'(a16,a)') 'Equations (H): ', '  Cowell (I)'
case(2)
   write(*,'(a16,a)') 'Equations (H): ', '  K-S'
case(3)
   write(*,'(a16,a)') 'Equations (H): ', '  EDromo'
end select
! Read choice of equations of motion (close encounter)
read(iin,'(a13,i2)') title, eqs_CE
select case(eqs_CE)
case default
   write(*,'(a,i2)') 'eqs_CE = ',eqs_CE
   write(*,'(a)') 'Invalid choice of the equations of motion for the close encounter phase.'
   write(*,'(a)') 'Stopping...'
   stop
case(-1)
   write(*,'(a16,a)') 'Equations (CE): ', '  Cowell (II)'
case(1)
   write(*,'(a16,a)') 'Equations (CE): ', '  Cowell (I)'
case(2)
   write(*,'(a16,a)') 'Equations (CE): ', '  K-S'
case(4)
   write(*,'(a16,a)') 'Equations (CE): ', '  GDromo'
end select

call SKIP(6)

! Read regularized formulation settings
read(iin,'(a13,i2)') title, flag_time_EDr
select case(flag_time_EDr)
case(0)
    write(*,'(a16,a)') 'EDromo TE: ', '  None'
case(1)
    write(*,'(a16,a)') 'EDromo TE: ', '  Constant'
case(2)
    write(*,'(a16,a)') 'EDromo TE: ', '  Linear'
end select
read(iin,'(a13,i2)') title, flag_time_GDr
select case(flag_time_GDr)
case(0)
    write(*,'(a16,a)') 'GDromo TE: ', '  None'
case(1)
    write(*,'(a16,a)') 'GDromo TE: ', '  Constant'
case(2)
    write(*,'(a16,a)') 'GDromo TE: ', '  Linear'
end select

call SKIP(2)

! Read output path
read(iin,'(a13,a)') title, path
write(*,'(a16,a)') 'Out path: ', '  '//trim(path)

close(iin)

contains

subroutine SKIP(n)
integer,intent(in)  ::  n
integer  ::  i
character(len=120)  ::  dummy
do i=1,n
    read(iin,'(a)') dummy
end do
end subroutine SKIP

end subroutine READ_SETTINGS


subroutine WRITE_SETTINGS(pathname,runID)

use SETTINGS
implicit none
character(len=*),intent(in)  ::  pathname,runID

open(id_set,file=pathname,status='replace',action='write')

write(id_set,'(120(''#''))')
write(id_set,'(2(''#''),55x,a7,54x,2(''#''))') 'NACE'
write(id_set,'(2(''#''),49x,a19,48x,2(''#''))') 'SIMULATION SETTINGS'
write(id_set,'(120(''#''),/)')

write(id_set,'(a16,a)') 'Run ID: ', '  '//runID
write(id_set,'(a16,g14.7,'', '',g14.7,a)') 'd interval: ', d_int, ' Earth radii'
write(id_set,'(a16,g14.7,'', '',g14.7)')   'e interval: ', e_int
write(id_set,'(a16,g14.7,'', '',g14.7,a)') 'theta interval: ', th_int, ' deg'
write(id_set,'(a16,g14.7,'', '',g14.7,a)') 'RSw interval: ', RSw_int, ' AU'
write(id_set,'(a16,i4)') 'd points: ',n_d
write(id_set,'(a16,i4)') 'e points: ', n_e
write(id_set,'(a16,i4)') 'theta points: ', n_th
write(id_set,'(a16,i4)') 'RSw points: ', n_rs
select case(espacing)
case(1)
    write(id_set,'(a16,a)') 'e spacing: ','  Linear'
case(2)
    write(id_set,'(a16,a)') 'e spacing: ','  Log'
end select
select case (integ)
case(1)
  write(id_set,'(a16,a)') 'Integrator: ','  LSODAR'
case(2)
  write(id_set,'(a16,a)') 'Integrator ','  Radau'
end select
write(id_set,'(a16,g14.7)') 'Ref. tol: ', tolref
write(id_set,'(a16,g14.7)') 'Test tol: ', tol
write(id_set,'(a16,g14.7,a)') 'dt_H: ', dt_H,' days'
write(id_set,'(a16,g14.7,a)') 'dt_CE: ', dt_CE,' days'
write(id_set,'(a16,i9)')    'Max steps: ',mxsteps
select case(eqs_H)
case(1)
   write(id_set,'(a16,a)') 'Equations (H): ', '  Cowell'
case(2)
   write(id_set,'(a16,a)') 'Equations (H): ', '  K-S'
case(3)
   write(id_set,'(a16,a)') 'Equations (H): ', '  EDromo'
end select
select case(eqs_CE)
case(1)
   write(id_set,'(a16,a)') 'Equations (CE): ', '  Cowell'
case(2)
   write(id_set,'(a16,a)') 'Equations (CE): ', '  K-S'
case(4)
   write(id_set,'(a16,a)') 'Equations (CE): ', '  GDromo'
end select
select case(flag_time_EDr)
case(0)
    write(id_set,'(a16,a)') 'EDromo TE: ', '  None'
case(1)
    write(id_set,'(a16,a)') 'EDromo TE: ', '  Constant'
case(2)
    write(id_set,'(a16,a)') 'EDromo TE: ', '  Linear'
end select
select case(flag_time_GDr)
case(0)
    write(id_set,'(a16,a)') 'GDromo TE: ', '  None'
case(1)
    write(id_set,'(a16,a)') 'GDromo TE: ', '  Constant'
case(2)
    write(id_set,'(a16,a)') 'GDromo TE: ', '  Linear'
end select
write(id_set,'(a16,a)') 'Out path: ', '  '//trim(path)

close(id_set)

end subroutine WRITE_SETTINGS


subroutine CREATE_ENRES(pathname,runID)

implicit none
! Arguments IN
character(len=*),intent(in)    ::  pathname,runID

! ==============================================================================================

open(unit=id_enres,file=pathname,status='replace',action='write')
write(id_enres,'(120(''#''))')
write(id_enres,'(2(''#''),55x,a7,54x,2(''#''))') 'NACE'
write(id_enres,'(2(''#''),50x,a16,50x,2(''#''))') 'ENERGY RESIDUALS'
write(id_enres,'(120(''#''),/,''#'')')

write(id_enres,'(''#'',a2,1x,3(a3,1x),7(a14,1x),a2,1x,a3)') &
&'i','j','k','l','dmin [R_E]','ecc [-]','theta [deg]', 'Rsw [AU]',&
&'dEN_IN [-]','dEn_OUT [-]','dEn_END [-]','MF','XC'
write(id_enres,'(127(''#''))')

end subroutine CREATE_ENRES


subroutine CREATE_ABSRES(pathname,runID)

implicit none
! Arguments IN
character(len=*),intent(in)    ::  pathname,runID

! ==============================================================================================

open(unit=id_ares,file=pathname,status='replace',action='write')
write(id_ares,'(120(''#''))')
write(id_ares,'(2(''#''),55x,a7,54x,2(''#''))') 'NACE'
write(id_ares,'(2(''#''),38x,a40,38x,2(''#''))') 'ABSOLUTE POSITION AND VELOCITY RESIDUALS'
write(id_ares,'(120(''#''),/,''#'')')

write(id_ares,'(''#'',a2,1x,3(a3,1x),10(a14,1x),a2,1x,a3)') &
&'i','j','k','l','dmin [R_E]','ecc [-]','theta [deg]','RSw [AU]',&
&'dR_IN [AU]','dR_OUT [AU]','dR_END [AU]','dV_IN [km/s]', 'dV_OUT [km/s]',&
&'dV_END [km/s]','MF','XC'
write(id_ares,'(172(''#''))')

end subroutine CREATE_ABSRES


subroutine CREATE_RELRES(pathname,runID)

implicit none
! Arguments IN
character(len=*),intent(in)    ::  pathname,runID

! ==============================================================================

open(unit=id_rres,file=pathname,status='replace',action='write')
write(id_rres,'(120(''#''))')
write(id_rres,'(2(''#''),55x,a7,54x,2(''#''))') 'NACE'
write(id_rres,'(2(''#''),38x,a40,38x,2(''#''))') 'RELATIVE POSITION AND VELOCITY RESIDUALS'
write(id_rres,'(120(''#''),/,''#'')')

write(id_rres,'(''#'',a2,1x,3(a3,1x),10(a14,1x),a2,1x,a3)') &
&'i','j','k','l','dmin [R_E]','ecc [-]','theta [deg]','RSw [AU]',&
&'dR_IN [-]','dR_OUT [-]','dR_END [-]','dV_IN [-]','dV_OUT [-]','dV_END [-]',&
&'MF','XC'
write(id_rres,'(172(''#''))')

end subroutine CREATE_RELRES


subroutine CREATE_STATS(pathname,runID)

implicit none
! Arguments IN
character(len=*),intent(in)    ::  pathname,runID

! ==============================================================================

open(unit=id_stats,file=pathname,status='replace',action='write')
write(id_stats,'(120(''#''))')
write(id_stats,'(2(''#''),55x,a7,54x,2(''#''))') 'NACE'
write(id_stats,'(2(''#''),47x,a22,47x,2(''#''))') 'INTEGRATION STATISTICS'
write(id_stats,'(120(''#''),/,''#'')')

write(id_stats,'(''#'',a2,1x,3(a3,1x),4(a14,1x),3(a9,1x),3(a14,1x),3(a10,1x),a14,1x,a3,1x,a3)')&
&'i','j','k','l','dmin [R_E]','ecc [-]','theta [deg]','RSw [AU]',&
&'FEVAL_H1','FEVAL_CE','FEVAL_END','MORD_H1','MORD_CE','MORD_END',&
&'ISTEPS_H1','ISTEPS_CE','ISTEPS_END','RSw_diff','MF','XC'
write(id_stats,'(206(''#''))')

end subroutine CREATE_STATS


subroutine WRITE_RESULTS(i,j,k,l,dmin,ecc,theta,RSw,&
&dR_abs,dR_rel,dV_abs,dV_rel,dEn_rel,RSw_diff,MF,XC,callsArr,istepsArr,mordArr,&
&calls_end,mord_end,isteps_end)

implicit none
! Arguments IN
integer,intent(in)   ::  i,j,k,l
integer,intent(in)   ::  MF(1:3),XC,callsArr(1:2),calls_end,istepsArr(1:2),isteps_end
real(qk),intent(in)  ::  dmin,ecc,theta,RSw
real(dk),intent(in)  ::  RSw_diff,dR_abs(1:3),dR_rel(1:3),dV_abs(1:3),dV_rel(1:3)
real(dk),intent(in)  ::  dEn_rel(1:3),mordArr(1:2),mord_end

write(id_enres,'(4(i3,'',''),4(g14.7,'',''),3(es14.7,'',''),3i1,'','',i3.3)')&
& i,j,k,l,dmin,ecc,theta,RSw,dEn_rel,MF,XC
write(id_ares,'(4(i3,'',''),4(g14.7,'',''),6(es14.7,'',''),3i1,'','',i3.3)')&
& i,j,k,l,dmin,ecc,theta,RSw,dR_abs,dV_abs,MF,XC
write(id_rres,'(4(i3,'',''),4(g14.7,'',''),6(es14.7,'',''),3i1,'','',i3.3)')&
& i,j,k,l,dmin,ecc,theta,RSw,dR_rel,dV_rel,MF,XC
write(id_stats,&
&'(4(i3,'',''),4(g14.7,'',''),3(i9,'',''),3(g14.7,'',''),3(i10,'',''),g14.7,'','',3i1,'','',i3.3)')&
&i,j,k,l,dmin,ecc,theta,RSw,callsArr,calls_end,mordArr,mord_end,istepsArr,isteps_end,RSw_diff,MF,XC


end subroutine WRITE_RESULTS


subroutine WRITE_EXCEPTION(i,j,k,l,dmin,ecc,theta,RSw,MF,XC)

implicit none
! Arguments IN
integer,intent(in)  ::  i,j,k,l,MF(1:3),XC
real(qk),intent(in) ::  dmin,ecc,theta,RSw
! Locals

write(id_enres,'(4(i3,'',''),4(g14.7,'',''),3(a14,'',''),3i1,'','',i3.3)')&
& i,j,k,l,dmin,ecc,theta,RSw,'NaN','NaN','NaN',MF,XC
write(id_ares,'(4(i3,'',''),4(g14.7,'',''),6(a14,'',''),3i1,'','',i3.3)')&
& i,j,k,l,dmin,ecc,theta,RSw,'NaN','NaN','NaN','NaN','NaN','NaN',MF,XC
write(id_rres,'(4(i3,'',''),4(g14.7,'',''),6(a14,'',''),3i1,'','',i3.3)')&
& i,j,k,l,dmin,ecc,theta,RSw,'NaN','NaN','NaN','NaN','NaN','NaN',MF,XC
write&
&(id_stats,'(4(i3,'',''),4(g14.7,'',''),3(a9,'',''),3(a14,'',''),3(a10,'',''),a14,'','',3i1,'','',i3.3)')&
&i,j,k,l,dmin,ecc,theta,RSw,'NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN',MF,XC

end subroutine WRITE_EXCEPTION



end module IO
