module KINDS
! Description:
!    Contains values for the kind parameters used in the program. To work in
!    quadruple precision, set dk = qk = 33. Remember that in this case only the
!    DLSODAR subroutines can be used.
! 
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
! 
! ==============================================================================

implicit none

integer,parameter  ::  ik = selected_int_kind(9)
integer,parameter  ::  dk = selected_real_kind(15)
integer,parameter  ::  qk = selected_real_kind(33)

end module KINDS
