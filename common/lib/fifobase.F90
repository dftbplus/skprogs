!> Contains the base fifo class.
module common_fifobase

  implicit none
  private

  public :: TFiFoBase, size


  !> Returns the size of the collection.
  interface size
    module procedure TFiFoBase_size
  end interface size


  !> Base fifo implementation managing pointers.
  type :: TFiFoBase
    private

    integer :: nitem = 0
    integer :: inmemory = 0
    integer :: memorylimit = -1

    integer :: fileid
    character(len=:), allocatable :: filename

    class(TFiFoNode), pointer :: head => null()
    class(TFiFoNode), pointer :: tail => null()
    class(TFiFoNode), pointer :: current => null()
    class(TFiFoNode), pointer :: previous => null()

  contains

    procedure :: initswap => TFiFoBase_initswap
    procedure :: pushptr => TFiFoBase_pushptr
    procedure :: popptr => TFiFoBase_popptr
    procedure :: getptr => TFiFoBase_getptr
    procedure :: getsize => TFiFoBase_size
    procedure :: reset => TFiFoBase_reset

    final :: TFiFoBase_destruct

    procedure, private :: writenodedata => TFiFoBase_writenodedata
    procedure, private :: readnodedata => TFiFoBase_readnodedata
    procedure, private :: freeresources => TFiFoBase_freeresources

    ! Workaround: should be private, but NAG fails to override private routines.
    procedure :: datafromfile => TFiFoBase_datafromfile
    procedure :: datatofile => TFiFoBase_datatofile

  end type TFiFoBase


  !> Represents one node in the fifo.
  type TFiFoNode

    class(*), pointer :: data => null()
    class(TFiFoNode), pointer :: next => null()
    integer :: filepos = -1

  end type TFiFoNode


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! FIFO Routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of items in the collection.
  !! \param obj  Collection instance.
  !! \return Number of items.
  pure function TFiFoBase_size(obj) result(res)

    class(TFiFoBase), intent(in) :: obj
    integer :: res

    res = obj%nitem

  end function TFiFoBase_size


  !> Initializes a swap for the collection.
  !!
  !! \details If swap is initialized for the collection, all entries above
  !! a given number are written to a file, instead of keeping them in memory.
  !! When the entries are read from the collection, a read buffer must be
  !! allocated, so the total number of elements kept in the memory will be
  !! increased by one.
  !!
  !! \param memorylimit  Maximal number of entries to keep in memory (-1: all
  !!   or 0: none or any positive number).
  !! \param filename  Name of the swap file.
  !! \param fileid  File id to use for handling the swap file.
  subroutine TFiFoBase_initswap(this, memorylimit, filename, fileid)
    class(TFiFoBase), intent(inout) :: this
    integer, intent(in) :: memorylimit
    character(len=*), intent(in) :: filename
    integer, intent(in) :: fileid

    if (this%memorylimit /= -1) then
      stop "FIFO swap can be initialized only once"
    end if
    this%memorylimit = memorylimit
    this%filename = filename
    this%fileid = fileid
  
  end subroutine TFiFoBase_initswap


  !> Pushes a pointer to the collection.
  !! \param this  Instance.
  !! \param data  Pointer to the data object.
  subroutine TFiFoBase_pushptr(this, data)
    class(TFiFoBase), intent(inout) :: this
    class(*), pointer, intent(in) :: data

    class(TFiFoNode), pointer :: node

    allocate(node)
    node%data => data
    if (.not. associated(this%head)) then
      this%head => node
      this%current => node
    else
      this%tail%next => node
    end if
    this%tail => node
    this%nitem = this%nitem + 1
    this%inmemory = this%inmemory + 1
    if (this%memorylimit /= -1 .and. this%inmemory > this%memorylimit) then
      call this%writenodedata(node)
    end if

  end subroutine TFiFoBase_pushptr


  !> Pops a pointer from the collection.
  !! \param this  Instance.
  !! \param data  Pointer to the data object on return.
  subroutine TFiFoBase_popptr(this, data)
    class(TFiFoBase), intent(inout) :: this
    class(*), pointer, intent(out) :: data

    class(TFiFoNode), pointer :: node

    if (.not. associated(this%head)) then
      data => null()
      return
    end if

    node => this%head
    this%head => node%next
    if (associated(node, this%current)) then
      this%current => node%next
    end if
    if (associated(node, this%previous)) then
      nullify(this%previous)
    end if
    if (.not. associated(node%data)) then
      call this%readnodedata(node)
    end if
    data => node%data
    deallocate(node)
    this%nitem = this%nitem - 1
    this%inmemory = this%inmemory - 1

  end subroutine TFiFoBase_popptr


  !> Gets a copy of a pointer from the collection.
  !! \param this  Instance.
  !! \param data  Pointer to the data object on return.
  subroutine TFiFoBase_getptr(this, data)
    class(TFiFoBase), intent(inout) :: this
    class(*), pointer, intent(out) :: data

    if (.not. associated(this%current)) then
      data => null()
      return
    end if

    ! If previous get read something from file, clear the buffer.
    if (associated(this%previous)) then
      if (this%previous%filepos /= -1 .and. associated(this%previous%data)) then
        deallocate(this%previous%data)
        this%inmemory = this%inmemory - 1
      end if
    end if
    
    if (.not. associated(this%current%data)) then
      call this%readnodedata(this%current)
    end if
    data => this%current%data

    this%previous => this%current
    if (associated(this%current%next)) then
      this%current => this%current%next
    else
      this%current => this%head
    end if

  end subroutine TFiFoBase_getptr


  !> Restets the collection to it initial (empty) state.
  !! \param this  Instance.
  subroutine TFiFoBase_reset(this)
    class(TFiFoBase), intent(inout) :: this

    call this%freeresources()
    this%nitem = 0
    this%inmemory = 0
    this%memorylimit = -1
    nullify(this%head, this%tail, this%current, this%previous)

  end subroutine TFiFoBase_reset


  !> Destructor for the class.
  !! \param this Instance.
  subroutine TFiFoBase_destruct(this)
    type(TFiFoBase), intent(inout) :: this

    call this%freeresources()

  end subroutine TFiFoBase_destruct


  !> Destroys the nodes in the collections and closes open files.
  !! \param this  Instance variable.
  subroutine TFiFoBase_freeresources(this)
    class(TFiFoBase), intent(inout) :: this

    class(TFiFoNode), pointer :: node
    logical :: opened

    node => this%head
    do while (associated(node))
      deallocate(node%data)
      this%head => node%next
      deallocate(node)
      node => this%head
    end do
    
    if (this%memorylimit /= -1) then
      inquire(this%fileid, opened=opened)
      if (opened) then
        close(this%fileid, status="delete")
      end if
    end if
    
  end subroutine TFiFoBase_freeresources


  !> Writes the data of a node to the disc and deallocates the data object.
  !! \param this  Instance.
  !! \param node  Node with the data that should be stored in a file.
  !! \note This routine invokes the data types write method instead of
  !!   writing the data directly.
  subroutine TFiFoBase_writenodedata(this, node)
    class(TFiFoBase), intent(inout) :: this
    class(TFiFoNode), pointer, intent(inout) :: node

    character(len=10) :: action

    inquire(this%fileid, action=action)
    if (action == "UNDEFINED") then
      ! No saved entries, create new swap file
      open(this%fileid, file=this%filename, access="stream", status="replace",&
          & action="write", form="unformatted", position="rewind")
    elseif (action == "READ") then
      ! Last commmand was pop/get, close file and and reopen in append mode.
      close(this%fileid)
      open(this%fileid, file=this%filename, access="stream", status="old",&
          & action="write", form="unformatted", position="append")
    end if

    inquire(this%fileid, pos=node%filepos)
    call this%datatofile(this%fileid, node%filepos, node%data)
    this%inmemory = this%inmemory - 1

  end subroutine TFiFoBase_writenodedata


  !> Reads the data of a node from file and allocates the data object.
  !! \param this  Instance.
  !! \param node  Node with the data that should be read from a file.
  !! \note This routine invokes the data types read method instead of
  !!   reading the data directly.
  subroutine TFiFoBase_readnodedata(this, node)
    class(TFiFoBase), intent(inout) :: this
    class(TFiFoNode), pointer, intent(inout) :: node

    character(len=10) :: action

    inquire(this%fileid, action=action)
    if (action == "WRITE") then
      close(this%fileid)
      open(this%fileid, file=this%filename, access="stream", status="old",&
          & action="read", form="unformatted")
    end if

    call this%datafromfile(this%fileid, node%filepos, node%data)
    this%inmemory = this%inmemory + 1

  end subroutine TFiFoBase_readnodedata


  !> Writes the content of a data node to a file.
  !!
  !! \details Extensions of the data object should rewrite it according to
  !! the data they contain.
  !!
  !! \param this  Instance.
  !! \param data  Pointer to a data node, will be deallocated at exit.
  !! \param fileid  File in which the data should be written.
  !! \param filepos  Position in the file, where the data must be written.
  subroutine TFiFoBase_datatofile(this, fileid, filepos, data)
    class(TFiFoBase), intent(inout) :: this
    integer, intent(in) :: fileid, filepos
    class(*), intent(inout), pointer :: data
    
    stop "Collection does not support swapping to file."
    
  end subroutine TFiFoBase_datatofile
  

  !> Reads the content of a data node from a file.
  !!
  !! \details Extensions of the data object should rewrite it according to
  !! the data they contain.
  !!
  !! \param this  Instance.
  !! \param fileid  File from which the data should be read.
  !! \param filepos  Position in the file, where the data should be read from.
  subroutine TFiFoBase_datafromfile(this, fileid, filepos, data)
    class(TFiFoBase), intent(inout) :: this
    integer, intent(in) :: fileid, filepos
    class(*), intent(out), pointer :: data

    stop "Collection does not support swapping to file."
    
  end subroutine TFiFoBase_datafromfile


end module common_fifobase
