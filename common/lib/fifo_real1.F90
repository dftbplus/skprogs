!> Implements fifo for rank 1 real (double precision) arrays.
module common_fifo_real1

  use common_fifobase, only : TFiFoBase, size

  implicit none
  private

  public :: TFiFoReal1

  integer, parameter :: dp = kind(1.0d0)


  !> Extended data type.
  type :: TMyData

    real(dp), allocatable :: data(:)

  end type TMyData


  !> Extended fifo.
  type, extends(TFiFoBase) :: TFiFoReal1

  contains

    procedure :: push => TFiFoReal1_push
    procedure :: pop => TFiFoReal1_pop
    procedure :: get => TFiFoReal1_get
    procedure :: push_alloc => TFiFoReal1_push_alloc
    procedure :: pop_alloc => TFiFoReal1_pop_alloc
    procedure :: popall => TFiFoReal1_popall
    procedure :: popall_concat => TFiFoReal1_popall_concat

    ! Workaround: should be private, but NAG fails to override private routines.
    procedure :: datatofile => TFiFoReal1_datatofile
    procedure :: datafromfile => TFiFoReal1_datafromfile

  end type TFiFoReal1


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! TFIFOREAL1 Routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Makes a copy of item and stores it in the collection.
  !! \param this  Instance.
  !! \param item  Item to store.
  subroutine TFiFoReal1_push(this, item)

    class(TFiFoReal1), intent(inout) :: this

    real(dp), intent(in) :: item(:)

    class(*), pointer :: wrapper

    allocate(TMyData :: wrapper)
    select type(wrapper)
    type is (TMyData)
      wrapper%data = item
    end select
    call this%pushptr(wrapper)

  end subroutine TFiFoReal1_push


  !> Retrieves the next item (fifo) and removes it from the collection.
  !! \param this  Instance.
  !! \param item  Item storing the result.
  subroutine TFiFoReal1_pop(this, item)

    class(TFiFoReal1), intent(inout) :: this

    real(dp), intent(out) :: item(:)

    class(*), pointer :: wrapper

    call this%popptr(wrapper)
    select type (wrapper)
    type is (TMyData)
      item(:) = wrapper%data
    end select
    deallocate(wrapper)
    
  end subroutine TFiFoReal1_pop


  !> Retrieves the next item without removing it from the collection.
  !!
  !! \details At first call the first element of the fifo is retrieved. At
  !! subsequent calls the elements are returned following the fifo principle. If
  !! the last element in the fifo had been returned, the first will be returned
  !! again.
  !!
  !! \param this  Instance.
  !! \param item  Item storing the result.
  subroutine TFiFoReal1_get(this, item)

    class(TFiFoReal1), intent(inout) :: this

    real(dp), intent(out) :: item(:)

    class(*), pointer :: wrapper

    call this%getptr(wrapper)
    select type (wrapper)
    type is (TMyData)
      item(:) = wrapper%data
    end select

  end subroutine TFiFoReal1_get


  !> Moves an allocatable item into the collection.
  !!
  !! \details Similar to push but for allocatable elements. The allocation
  !! status of the item is moved to the collection, so that the original item is
  !! automatically deallocated. No temporary copy of the item is created.
  !!
  !! \param this  Instance.
  !! \param item  Item to store. Deallocated on return.
  subroutine TFiFoReal1_push_alloc(this, item)

    class(TFiFoReal1), intent(inout) :: this

    real(dp), allocatable, intent(inout) :: item(:)

    class(*), pointer :: wrapper

    allocate(TMyData :: wrapper)
    select type (wrapper)
    type is (TMyData)
      call move_alloc(item, wrapper%data)
    end select
    call this%pushptr(wrapper)

  end subroutine TFiFoReal1_push_alloc


  !> Retrieves the next item (fifo) and removes it from the collection.
  !!
  !! \details Similar to pop but for allocatable elements. The allocation status
  !! is moved from the collection to the item, so that the item will be
  !! automatically allocated. No temporary copy of the item is created.
  !!
  !! \param this  Instance.
  !! \param item  Item storing the result.
  subroutine TFiFoReal1_pop_alloc(this, item)

    class(TFiFoReal1), intent(inout) :: this

    real(dp), allocatable, intent(out) :: item(:)

    class(*), pointer :: wrapper

    call this%popptr(wrapper)
    select type (wrapper)
    type is (TMyData)
      call move_alloc(wrapper%data, item)
    end select
    deallocate(wrapper)
    
  end subroutine TFiFoReal1_pop_alloc


  !> Retrieves all items from the collection as an allocatable array and deletes
  !! them.
  !!
  !! \details The routine allocates an array with the given shape and an
  !! additional dimension with the size of the collectoin.
  !!
  !! \param this  Instance.
  !! \param itemshape  Shape of the items in the collection.
  !! \param items  Array containing the items.
  !!
  !! \warning It is the responsibility of the caller to invoke this method
  !! only on collections containing elements with the same shape. No checking
  !! of shape conformance is done.
  subroutine TFiFoReal1_popall(this, items)

    class(TFiFoReal1), intent(inout) :: this

    real(dp), allocatable, intent(out) :: items(:,:)

    class(*), pointer :: wrapper

    integer :: itemshape(1)

    !> Auxiliary variable
    integer :: ii

    call this%getptr(wrapper)
    select type (wrapper)
    type is (TMyData)
      itemshape(:) = shape(wrapper%data)
    end select
    allocate(items(itemshape(1), size(this)))
    do ii = 1, size(this)
      call this%pop(items(:, ii))
    end do

  end subroutine TFiFoReal1_popall


  !> Retrieves all items from the collection as an allocatable array by
  !! concatenating them and deletes them.
  !!
  !! \details The routine allocates an array with the given shape times
  !! the size of the collection.
  !!
  !! \param this  Instance.
  !! \param items  Array containing the items.
  subroutine TFiFoReal1_popall_concat(this, items)

    class(TFiFoReal1), intent(inout) :: this

    real(dp), allocatable, intent(out) :: items(:)

    integer :: ii, ind, total, nn

    class(*), pointer :: wrapper

    total = 0
    do ii = 1, size(this)
      call this%getptr(wrapper)
      select type (wrapper)
      type is (TMyData)
        total = total + size(wrapper%data)
      end select
    end do
    allocate(items(total))
    ind = 1
    do ii = 1, size(this)
      call this%popptr(wrapper)
      select type (wrapper)
      type is (TMyData)
        nn = size(wrapper%data)
        items(ind:ind+nn-1) = wrapper%data(:)
        ind = ind + nn
      end select
      deallocate(wrapper)
    end do
    
  end subroutine TFiFoReal1_popall_concat


  !> Overrides the datatofile method of the base class.
  !! \param this  Instance.
  !! \param fileid  Id of the file in which data should be written.
  !! \param filepos  Position in the file, to which data should be written.
  !! \param data  Data node to save to file.
  subroutine TFiFoReal1_datatofile(this, fileid, filepos, data)

    class(TFiFoReal1), intent(inout) :: this

    integer, intent(in) :: fileid, filepos

    class(*), pointer, intent(inout) :: data

    select type (data)
    type is (TMyData)
      write(fileid, pos=filepos) shape(data%data)
      write(fileid) data%data
    end select
    deallocate(data)

  end subroutine TFiFoReal1_datatofile


  !> Overides the datafromfile method of the base class.
  !! \param this  Instance.
  !! \param fileid  Id of the file from which data should be read.
  !! \param filepos  Position in the file, from which data should be read.
  !! \param data  Data node to create from file.
  subroutine TFiFoReal1_datafromfile(this, fileid, filepos, data)

    class(TFiFoReal1), intent(inout) :: this

    integer, intent(in) :: fileid, filepos

    class(*), pointer, intent(out) :: data

    integer :: itemshape(1)

    allocate(TMyData :: data)
    select type (data)
    type is (TMyData)
      read(fileid, pos=filepos) itemshape
      allocate(data%data(itemshape(1)))
      read(fileid) data%data
    end select

  end subroutine TFiFoReal1_datafromfile


end module common_fifo_real1
