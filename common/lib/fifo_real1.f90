!> Implements fifo for rank 1 real (double precision) arrays.
module fifo_real1_module
  use fifobase_module
  implicit none
  private

  public :: fifo_real1, size

  integer, parameter :: dp = kind(1.0d0)

  !> Extended data type.
  type :: mydata
    real(dp), allocatable :: data(:)
  end type mydata

  !> Extendid fifo.
  type, extends(fifobase) :: fifo_real1
  contains
    procedure :: push => fifo_real1_push
    procedure :: pop => fifo_real1_pop
    procedure :: get => fifo_real1_get
    procedure :: push_alloc => fifo_real1_push_alloc
    procedure :: pop_alloc => fifo_real1_pop_alloc
    procedure :: popall => fifo_real1_popall
    procedure :: popall_concat => fifo_real1_popall_concat
    ! Workaround: should be private, but NAG fails to override private routines.
    procedure :: datatofile => fifo_real1_datatofile
    procedure :: datafromfile => fifo_real1_datafromfile
  end type fifo_real1


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! FIFO_REAL1 Routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Makes a copy of item and stores it in the collection.
  !! \param self  Instance.
  !! \param item  Item to store.
  subroutine fifo_real1_push(self, item)
    class(fifo_real1), intent(inout) :: self
    real(dp), intent(in) :: item(:)

    class(*), pointer :: wrapper

    allocate(mydata :: wrapper)
    select type(wrapper)
    type is (mydata)
      wrapper%data = item        ! Automatic allocation
    end select
    call self%pushptr(wrapper)

  end subroutine fifo_real1_push


  !> Retrieves the next item (fifo) and removes it from the collection.
  !! \param self  Instance.
  !! \param item  Item storing the result.
  subroutine fifo_real1_pop(self, item)
    class(fifo_real1), intent(inout) :: self
    real(dp), intent(out) :: item(:)

    class(*), pointer :: wrapper

    call self%popptr(wrapper)
    select type (wrapper)
    type is (mydata)
      item(:) = wrapper%data
    end select
    deallocate(wrapper)
    
  end subroutine fifo_real1_pop


  !> Retrieves the next item without removing it from the collection.
  !!
  !! \details At first call the first element of the fifo is retrieved. At
  !! subsequent calls the elements are returned following the fifo principle. If
  !! the last element in the fifo had been returned, the first will be returned
  !! again.
  !!
  !! \param self  Instance.
  !! \param item  Item storing the result.
  subroutine fifo_real1_get(self, item)
    class(fifo_real1), intent(inout) :: self
    real(dp), intent(out) :: item(:)

    class(*), pointer :: wrapper

    call self%getptr(wrapper)
    select type (wrapper)
    type is (mydata)
      item(:) = wrapper%data
    end select

  end subroutine fifo_real1_get


  !> Moves an allocatable item into the collection.
  !!
  !! \details Similar to push but for allocatable elements. The allocation
  !! status of the item is moved to the collection, so that the original item is
  !! automatically deallocated. No temporary copy of the item is created.
  !!
  !! \param self  Instance.
  !! \param item  Item to store. Deallocated on return.
  subroutine fifo_real1_push_alloc(self, item)
    class(fifo_real1), intent(inout) :: self
    real(dp), allocatable, intent(inout) :: item(:)

    class(*), pointer :: wrapper

    allocate(mydata :: wrapper)
    select type (wrapper)
    type is (mydata)
      call move_alloc(item, wrapper%data)
    end select
    call self%pushptr(wrapper)

  end subroutine fifo_real1_push_alloc


  !> Retrieves the next item (fifo) and removes it from the collection.
  !!
  !! \details Similar to pop but for allocatable elements. The allocation status
  !! is moved from the collection to the item, so that the item will be
  !! automatically allocated. No temporary copy of the item is created.
  !!
  !! \param self  Instance.
  !! \param item  Item storing the result.
  subroutine fifo_real1_pop_alloc(self, item)
    class(fifo_real1), intent(inout) :: self
    real(dp), allocatable, intent(out) :: item(:)

    class(*), pointer :: wrapper

    call self%popptr(wrapper)
    select type (wrapper)
    type is (mydata)
      call move_alloc(wrapper%data, item)
    end select
    deallocate(wrapper)
    
  end subroutine fifo_real1_pop_alloc


  !> Retrieves all items from the collection as an allocatable array and deletes
  !! them.
  !!
  !! \details The routine allocates an array with the given shape and an
  !! additional dimension with the size of the collectoin.
  !!
  !! \param self  Instance.
  !! \param itemshape  Shape of the items in the collection.
  !! \param items  Array containing the items.
  !!
  !! \warning It is the responsibility of the caller to invoke this method
  !! only on collections containing elements with the same shape. No checking
  !! of shape conformance is done.
  subroutine fifo_real1_popall(self, items)
    class(fifo_real1), intent(inout) :: self
    real(dp), allocatable, intent(out) :: items(:,:)

    class(*), pointer :: wrapper
    integer :: itemshape(1)
    integer :: ii

    call self%getptr(wrapper)
    select type (wrapper)
    type is (mydata)
      itemshape(:) = shape(wrapper%data)
    end select
    allocate(items(itemshape(1), size(self)))
    do ii = 1, size(self)
      call self%pop(items(:,ii))
    end do

  end subroutine fifo_real1_popall


  !> Retrieves all items from the collection as an allocatable array by
  !! concatenating them and deletes them.
  !!
  !! \details The routine allocates an array with the given shape times
  !! the size of the collection.
  !!
  !! \param self  Instance.
  !! \param items  Array containing the items.
  subroutine fifo_real1_popall_concat(self, items)
    class(fifo_real1), intent(inout) :: self
    real(dp), allocatable, intent(out) :: items(:)

    integer :: ii, ind, total, nn
    class(*), pointer :: wrapper

    total = 0
    do ii = 1, size(self)
      call self%getptr(wrapper)
      select type (wrapper)
      type is (mydata)
        total = total + size(wrapper%data)
      end select
    end do
    allocate(items(total))
    ind = 1
    do ii = 1, size(self)
      call self%popptr(wrapper)
      select type (wrapper)
      type is (mydata)
        nn = size(wrapper%data)
        items(ind:ind+nn-1) = wrapper%data(:)
        ind = ind + nn
      end select
      deallocate(wrapper)
    end do
    
  end subroutine fifo_real1_popall_concat


  !> Overrides the datatofile method of the base class.
  !! \param self  Instance.
  !! \param fileid  Id of the file in which data should be written.
  !! \param filepos  Position in the file, to which data should be written.
  !! \param data  Data node to save to file.
  subroutine fifo_real1_datatofile(self, fileid, filepos, data)
    class(fifo_real1), intent(inout) :: self
    integer, intent(in) :: fileid, filepos
    class(*), pointer, intent(inout) :: data

    select type (data)
    type is (mydata)
      write(fileid, pos=filepos) shape(data%data)
      write(fileid) data%data
    end select
    deallocate(data)

  end subroutine fifo_real1_datatofile


  !> Overides the datafromfile method of the base class.
  !! \param self  Instance.
  !! \param fileid  Id of the file from which data should be read.
  !! \param filepos  Position in the file, from which data should be read.
  !! \param data  Data node to create from file.
  subroutine fifo_real1_datafromfile(self, fileid, filepos, data)
    class(fifo_real1), intent(inout) :: self
    integer, intent(in) :: fileid, filepos
    class(*), pointer, intent(out) :: data

    integer :: itemshape(1)

    allocate(mydata :: data)
    select type (data)
    type is (mydata)
      read(fileid, pos=filepos) itemshape
      allocate(data%data(itemshape(1)))
      read(fileid) data%data
    end select

  end subroutine fifo_real1_datafromfile


end module fifo_real1_module
