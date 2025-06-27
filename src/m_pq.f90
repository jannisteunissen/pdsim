!> Implementation of a min priority queue using a binary min-heap.
!>
!> Author: Jannis Teunissen
module m_pq

  implicit none
  private

  !> Floating point type for priorities
  integer, parameter :: dp = kind(0.0d0)

  !> Node with integer priority
  type node_ip_t
     integer :: priority
     integer :: val
  end type node_ip_t

  !> Node with real(dp) priority
  type node_rp_t
     real(dp) :: priority
     integer  :: val
  end type node_rp_t

  !> Priority queue with integer priority
  type pqi_t
     type(node_ip_t), allocatable :: nodes(:)
     integer                      :: n_stored
  end type pqi_t

  !> Priority queue with real(dp) priority
  type pqr_t
     type(node_rp_t), allocatable :: nodes(:)
     integer                      :: n_stored
  end type pqr_t

  public :: pqi_t
  public :: pqi_create
  public :: pqi_destroy
  public :: pqi_push
  public :: pqi_pop
  public :: pqi_peep

  public :: pqr_t
  public :: pqr_create
  public :: pqr_destroy
  public :: pqr_push
  public :: pqr_pop
  public :: pqr_peep

contains

  !> Create an empty priority queue with a given initial size
  subroutine pqi_create(pq, initial_size)
    type(pqi_t), intent(out)      :: pq
    integer, intent(in), optional :: initial_size

    !> Default size for new priority queues (if not given)
    integer, parameter :: default_initial_size = 256

    if (present(initial_size)) then
       allocate(pq%nodes(initial_size))
    else
       allocate(pq%nodes(default_initial_size))
    end if

    pq%n_stored = 0
  end subroutine pqi_create

  !> Clear all storage associated with the priority queue
  subroutine pqi_destroy(pq)
    type(pqi_t), intent(inout) :: pq

    if (allocated(pq%nodes)) deallocate(pq%nodes)
    pq%n_stored = 0
  end subroutine pqi_destroy

  !> Add a new element to the queue
  subroutine pqi_push(pq, val, priority)
    type(pqi_t), intent(inout)   :: pq
    integer, intent(in)          :: val
    integer, intent(in)          :: priority
    integer                      :: i, i_parent, new_size
    type(node_ip_t), allocatable :: tmp(:)

    ! Resize if array is full
    if (pq%n_stored == size(pq%nodes)) then
       new_size = 2 * size(pq%nodes)
       allocate(tmp(new_size))
       tmp(1:pq%n_stored) = pq%nodes
       call move_alloc(tmp, pq%nodes)
    end if

    i = pq%n_stored + 1
    i_parent = i / 2

    ! Bubble up. The new element will repeatedly be swapped with its parent at
    ! i/2 until its parent has a lower or the same priority
    do while (i > 1)
       if (pq%nodes(i_parent)%priority <= priority) exit
       pq%nodes(i) = pq%nodes(i_parent)
       i = i_parent
       i_parent = i_parent / 2
    end do

    pq%nodes(i)%priority = priority
    pq%nodes(i)%val = val
    pq%n_stored = pq%n_stored + 1
  end subroutine pqi_push

  !> Get the element with the lowest priority
  subroutine pqi_peep(pq, val, priority)
    type(pqi_t), intent(in) :: pq
    integer, intent(out)    :: val
    integer, intent(out)    :: priority
    val      = pq%nodes(1)%val
    priority = pq%nodes(1)%priority
  end subroutine pqi_peep

  !> Get the item with the lowest priority and remove it from the queue
  subroutine pqi_pop(pq, val, priority)
    type(pqi_t), intent(inout) :: pq
    integer, intent(out)       :: val
    integer, intent(out)       :: priority
    integer                    :: i, i_child

    if (pq%n_stored <= 0) error stop "Cannot pop from empty queue"

    val      = pq%nodes(1)%val
    priority = pq%nodes(1)%priority

    ! Bubble down. 'Virtually' place last element at index 1, and repeatedly
    ! swap with the lowest-priority child (at 2*i or 2*i+1) as long as it has
    ! a lower priority.
    pq%n_stored = pq%n_stored - 1
    i = 1
    do
       i_child = 2 * i

       if (i_child < pq%n_stored) then
          ! Select child with lowest priority
          if (pq%nodes(i_child+1)%priority < pq%nodes(i_child)%priority) then
             i_child = i_child + 1
          end if
       else if (i_child > pq%n_stored) then
          ! No children
          exit
       end if

       if (pq%nodes(i_child)%priority < pq%nodes(pq%n_stored+1)%priority) then
          pq%nodes(i) = pq%nodes(i_child)
          i = i_child
       else
          exit
       end if
    end do

    pq%nodes(i) = pq%nodes(pq%n_stored+1)
  end subroutine pqi_pop

  !> Create an empty priority queue with a given initial size
  subroutine pqr_create(pq, initial_size)
    type(pqr_t), intent(out)      :: pq
    integer, intent(in), optional :: initial_size

    !> Default size for new priority queues (if not given)
    integer, parameter :: default_initial_size = 256

    if (present(initial_size)) then
       allocate(pq%nodes(initial_size))
    else
       allocate(pq%nodes(default_initial_size))
    end if

    pq%n_stored = 0
  end subroutine pqr_create

  !> Clear all storage associated with the priority queue
  subroutine pqr_destroy(pq)
    type(pqr_t), intent(inout) :: pq

    if (allocated(pq%nodes)) deallocate(pq%nodes)
    pq%n_stored = 0
  end subroutine pqr_destroy

  !> Add a new element to the queue
  subroutine pqr_push(pq, val, priority)
    type(pqr_t), intent(inout)   :: pq
    integer, intent(in)          :: val
    real(dp), intent(in)         :: priority
    integer                      :: i, i_parent, new_size
    type(node_rp_t), allocatable :: tmp(:)

    ! Resize if array is full
    if (pq%n_stored == size(pq%nodes)) then
       new_size = 2 * size(pq%nodes)
       allocate(tmp(new_size))
       tmp(1:pq%n_stored) = pq%nodes
       call move_alloc(tmp, pq%nodes)
    end if

    i = pq%n_stored + 1
    i_parent = i / 2

    ! Bubble up. The new element will repeatedly be swapped with its parent at
    ! i/2 until its parent has a lower or the same priority
    do while (i > 1)
       if (pq%nodes(i_parent)%priority <= priority) exit
       pq%nodes(i) = pq%nodes(i_parent)
       i = i_parent
       i_parent = i_parent / 2
    end do

    pq%nodes(i)%priority = priority
    pq%nodes(i)%val = val
    pq%n_stored = pq%n_stored + 1
  end subroutine pqr_push

  !> Get the element with the lowest priority
  subroutine pqr_peep(pq, val, priority)
    type(pqr_t), intent(in) :: pq
    integer, intent(out)    :: val
    real(dp), intent(out)   :: priority
    val      = pq%nodes(1)%val
    priority = pq%nodes(1)%priority
  end subroutine pqr_peep

  !> Get the item with the lowest priority and remove it from the queue
  subroutine pqr_pop(pq, val, priority)
    type(pqr_t), intent(inout) :: pq
    integer, intent(out)       :: val
    real(dp), intent(out)      :: priority
    integer                    :: i, i_child

    if (pq%n_stored <= 0) error stop "Cannot pop from empty queue"

    val      = pq%nodes(1)%val
    priority = pq%nodes(1)%priority

    ! Bubble down. 'Virtually' place last element at index 1, and repeatedly
    ! swap with the lowest-priority child (at 2*i or 2*i+1) as long as it has
    ! a lower priority.
    pq%n_stored = pq%n_stored - 1
    i = 1
    do
       i_child = 2 * i

       if (i_child < pq%n_stored) then
          ! Select child with lowest priority
          if (pq%nodes(i_child+1)%priority < pq%nodes(i_child)%priority) then
             i_child = i_child + 1
          end if
       else if (i_child > pq%n_stored) then
          ! No children
          exit
       end if

       if (pq%nodes(i_child)%priority < pq%nodes(pq%n_stored+1)%priority) then
          pq%nodes(i) = pq%nodes(i_child)
          i = i_child
       else
          exit
       end if
    end do

    pq%nodes(i) = pq%nodes(pq%n_stored+1)
  end subroutine pqr_pop

end module m_pq
