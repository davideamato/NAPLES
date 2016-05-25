module TYPEDEF
    
    ! Description:
    !    Contains derived type definitions for the Benchmark_Flyby program.
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
    
    use KINDS, only: dk,qk
    
    ! VARIABLES AND DECLARATIONS
    
    implicit none
    
    ! Cartesian trajectory type definition (quad)
    type  ::  qtraj
        real(qk)             ::  rad(1:3)
        real(qk)             ::  vel(1:3)
        real(qk)             ::  t
        real(qk)             ::  rmag
        real(qk)             ::  vmag
        logical              ::  inSOI
        type(qtraj),pointer  ::  next
        type(qtraj),pointer  ::  prev
    contains
        procedure,pass  ::  QSET_TRACK
    end type qtraj
    
    ! Cartesian trajectory type definition (double)
    type  ::  dtraj
        real(dk)             ::  rad(1:3)
        real(dk)             ::  vel(1:3)
        real(dk)             ::  t
        real(dk)             ::  rmag
        real(dk)             ::  vmag
        logical              ::  inSOI
        type(dtraj),pointer  ::  next
        type(dtraj),pointer  ::  prev
    contains
        procedure,pass  ::  DSET_TRACK
    end type dtraj
    
    contains
    
    ! ==============================================================================================
    !
    !                                          PROCEDURES
    !
    ! ==============================================================================================
    
    subroutine QTRAJ_LINK(node)
    
    ! Description:
    !    Allocates and links a new node of the trajectory linked list and fills in the appropriate
    !    fields.
    ! 
    ! ==============================================================================================
    !                                            PREAMBLE
    ! ==============================================================================================
    !
    ! MODULES
    
    ! VARIABLES
    
    implicit none
    
    ! Arguments
    type(qtraj),pointer,intent(inout)   ::  node                     ! List node
    
    ! Locals
    integer             ::  all_err   ! Allocation error code
    character(len=80)   ::  all_msg   ! Allocation error message
    
    ! ==============================================================================================
    !                                            EXECUTION
    ! ==============================================================================================
    
    ! Build new node
    allocate(node%next,stat=all_err,errmsg=all_msg); if (all_err/=0) write(*,*),all_msg
    node%next%prev => node
    node           => node%next
    nullify(node%next)
    
    end subroutine QTRAJ_LINK
    
    subroutine DTRAJ_LINK(node)
    
    ! Description:
    !    Allocates and links a new node of the trajectory linked list and fills in the appropriate
    !    fields.
    ! 
    ! ==============================================================================================
    !                                            PREAMBLE
    ! ==============================================================================================
    !
    ! MODULES
    
    ! VARIABLES
    
    implicit none
    
    ! Arguments
    type(dtraj),pointer,intent(inout)   ::  node                     ! List node
    
    ! Locals
    integer             ::  all_err   ! Allocation error code
    character(len=80)   ::  all_msg   ! Allocation error message
    
    ! ==============================================================================================
    !                                            EXECUTION
    ! ==============================================================================================
    
    ! Build new node
    allocate(node%next,stat=all_err,errmsg=all_msg); if (all_err/=0) write(*,*),all_msg
    node%next%prev => node
    node           => node%next
    nullify(node%next)
    
    end subroutine DTRAJ_LINK
    
    subroutine QTRAJ_INIT(node_h,node_t)
    
    ! Description:
    !    Initializes the trajectory linked list.
    ! 
    ! ==============================================================================================
    !                                            PREAMBLE
    ! ==============================================================================================
    !
    ! MODULES
    
    ! VARIABLES
    
    implicit none
    
    ! Arguments
    ! State vectors
    type(qtraj),pointer,intent(inout)   ::  node_h,node_t            ! Tail node (points to %next)
    
    ! Locals
    integer             ::  all_err   ! Allocation error code
    character(len=80)   ::  all_msg   ! Allocation error message
    
    ! ==============================================================================================
    !                                            EXECUTION
    ! ==============================================================================================
    
    ! Deallocate nodes if already allocated
    if (associated(node_h)) nullify(node_h)
    if (associated(node_t)) nullify(node_t)
    
    ! Build head of the list
    allocate(node_h,stat=all_err,errmsg=all_msg); if (all_err/=0) write(*,*),all_msg
    nullify(node_h%prev)
    node_t => node_h
    nullify(node_t%next)
    node_t%prev => node_h

    end subroutine QTRAJ_INIT
    
    subroutine DTRAJ_INIT(node_h,node_t)
    
    ! Description:
    !    Initializes the trajectory linked list.
    ! 
    ! ==============================================================================================
    !                                            PREAMBLE
    ! ==============================================================================================
    !
    ! MODULES
    
    ! VARIABLES
    
    implicit none
    
    ! Arguments
    ! State vectors
    type(dtraj),pointer,intent(inout)   ::  node_h,node_t            ! Tail node (points to %next)
    
    ! Locals
    integer             ::  all_err   ! Allocation error code
    character(len=80)   ::  all_msg   ! Allocation error message
    
    ! ==============================================================================================
    !                                            EXECUTION
    ! ==============================================================================================
    
    ! Deallocate nodes if already allocated
    if (associated(node_h)) nullify(node_h)
    if (associated(node_t)) nullify(node_t)
    
    ! Build head of the list
    allocate(node_h,stat=all_err,errmsg=all_msg); if (all_err/=0) write(*,*),all_msg
    nullify(node_h%prev)
    node_t => node_h
    nullify(node_t%next)
    node_t%prev => node_h

    end subroutine DTRAJ_INIT
    
    subroutine QSET_TRACK(this,rad,vel,t,inSOI)
    
    ! VARIABLES
    implicit none
    
    class(qtraj),intent(out)  ::  this
    real(qk),intent(in)       ::  rad(1:3),vel(1:3),t
    logical,intent(in)        ::  insoi
    real(qk)                  ::  radmag,velmag
    
    ! STATEMENTS
    radmag = sqrt(dot_product(rad,rad))
    velmag = sqrt(dot_product(vel,vel))
    
    ! Store position and velocity and compute their magnitudes
    this%rad   = rad
    this%vel   = vel
    this%t     = t
    this%rmag  = radmag
    this%vmag  = velmag
    this%inSOI = inSOI
    
    end subroutine QSET_TRACK
    
    subroutine DSET_TRACK(this,rad,vel,t,inSOI)
    
    ! VARIABLES
    implicit none
    
    class(dtraj),intent(out)  ::  this
    real(dk),intent(in)       ::  rad(1:3),vel(1:3),t
    logical,intent(in)        ::  insoi
    real(dk)                  ::  radmag,velmag
    
    ! STATEMENTS
    radmag = sqrt(dot_product(rad,rad))
    velmag = sqrt(dot_product(vel,vel))
    
    ! Store position and velocity and compute their magnitudes
    this%rad   = rad
    this%vel   = vel
    this%t     = t
    this%rmag  = radmag
    this%vmag  = velmag
    this%inSOI = inSOI
    
    end subroutine DSET_TRACK
    
    subroutine QLIST2VEC(TH,y)
    ! Description:
    !     Converts a trajectory list into an array of state vectors.
    !
    ! ==============================================================================================
    !                                            PREAMBLE
    ! ==============================================================================================
    
    ! VARIABLES
    implicit none
    ! Arguments IN
    type(qtraj),pointer,intent(in)  ::  TH               ! Head of trajectory list
    ! Arguments OUT
    real(qk),allocatable            ::  y(:,:)           ! Trajectory state vector array
    
    ! Locals
    type(qtraj),pointer  ::  ptr  ! Tjectory pointer
    integer              ::  i
    integer              ::  n    ! Number of elements in TH
    character(len=80)    ::  amsg
    integer              ::  astat
    
    ! ==============================================================================================
    !                                            EXECUTION
    ! ==============================================================================================
    
    ! ==============================================================================================
    ! 01. GET NUMBER OF ELEMENTS OF THE LIST
    ! ==============================================================================================
    
    ! Initialize ptr, n
    ptr => TH
    n = 0
    
    do
        if(.not.(associated(ptr%next))) exit
        n = n + 1
        ptr => ptr%next
    end do
    
    ! ==============================================================================================
    ! 02. COPY TRAJECTORY INTO VECTOR
    ! ==============================================================================================
    
    ! Allocate the output vector
    allocate(y(1:n,1:7),stat=astat,errmsg=amsg); if (astat /=0) write(*,*),amsg
    
    ! Copy values into the output vector
    ptr => TH
    do i = 1,n
        y(i,1)   = ptr%t
        y(i,2:4) = ptr%rad
        y(i,5:7) = ptr%vel
        ptr => ptr%next
    end do
    
    end subroutine QLIST2VEC
    
    subroutine DLIST2VEC(TH,y)
    ! Description:
    !     Converts a trajectory list into an array of state vectors.
    !
    ! ==============================================================================================
    !                                            PREAMBLE
    ! ==============================================================================================
    
    ! VARIABLES
    implicit none
    ! Arguments IN
    type(dtraj),pointer,intent(in)  ::  TH               ! Head of trajectory list
    ! Arguments OUT
    real(dk),allocatable            ::  y(:,:)           ! Trajectory state vector array
    
    ! Locals
    type(dtraj),pointer  ::  ptr  ! Trajectory pointer
    integer              ::  i
    integer              ::  n    ! Number of elements in TH
    character(len=80)    ::  amsg
    integer              ::  astat
    
    ! ==============================================================================================
    !                                            EXECUTION
    ! ==============================================================================================
    
    ! ==============================================================================================
    ! 01. GET NUMBER OF ELEMENTS OF THE LIST
    ! ==============================================================================================
    
    ! Initialize ptr, n
    ptr => TH
    n = 0
    
    do
        if(.not.(associated(ptr%next))) exit
        n = n + 1
        ptr => ptr%next
    end do
    
    ! ==============================================================================================
    ! 02. COPY TRAJECTORY INTO VECTOR
    ! ==============================================================================================
    
    ! Allocate the output vector
    allocate(y(1:n,1:7),stat=astat,errmsg=amsg); if (astat /=0) write(*,*),amsg
    
    ! Copy values into the output vector
    ptr => TH
    do i = 1,n
        y(i,1)   = ptr%t
        y(i,2:4) = ptr%rad
        y(i,5:7) = ptr%vel
        ptr => ptr%next
    end do
    
    end subroutine DLIST2VEC
    
    subroutine QTRAJ_DEALLOC(head)
    ! Description:
    !    Deallocates a trajectory list.
    ! 
    ! ==============================================================================================
    !                                            PREAMBLE
    ! ==============================================================================================
        
    ! VARIABLES
    implicit none
    
    ! Arguments
    type(qtraj),pointer,intent(in)     ::  head        ! Head of list
    type(qtraj),pointer                ::  curr,next   ! Pointers
    
    ! ==============================================================================================
    !                                            EXECUTION
    ! ==============================================================================================
    
    curr => head
    next => head%next
    do
        if (.not.(associated(curr))) exit
        next => curr%next
        deallocate(curr)
        curr => next
    end do
    
    end subroutine QTRAJ_DEALLOC
    
    subroutine DTRAJ_DEALLOC(head)
    ! Description:
    !    Deallocates a trajectory list.
    ! 
    ! ==============================================================================================
    !                                            PREAMBLE
    ! ==============================================================================================
        
    ! VARIABLES
    implicit none
    
    ! Arguments
    type(dtraj),pointer,intent(in)     ::  head        ! Head of list
    type(dtraj),pointer                ::  curr,next   ! Pointers
    
    ! ==============================================================================================
    !                                            EXECUTION
    ! ==============================================================================================
    
    curr => head
    next => head%next
    do
        if (.not.(associated(curr))) exit
        next => curr%next
        deallocate(curr)
        curr => next
    end do
    
    end subroutine DTRAJ_DEALLOC
    
end module TYPEDEF
