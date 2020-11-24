!> Contains the base fifo class.
module fifobase_module
  implicit none
  private

  public :: fifobase, size


  !> Returns the size of the collection.
  interface size
    module procedure fifo_size
  end interface size

  !> Base fifo implementation managing pointers.
  type :: fifobase
    private
    integer :: nitem = 0
    integer :: inmemory = 0
    integer :: memorylimit = -1
    class(fifonode), pointer :: head => null()
    class(fifonode), pointer :: tail => null()
    class(fifonode), pointer :: current => null()
    class(fifonode), pointer :: previous => null()
    integer :: fileid
    character(:), allocatable :: filename
  contains
    procedure :: initswap => fifo_initswap
    procedure :: pushptr => fifo_pushptr
    procedure :: popptr => fifo_popptr
    procedure :: getptr => fifo_getptr
    procedure :: getsize => fifo_size
    procedure :: reset => fifo_reset
    final :: fifo_destruct
    procedure, private :: writenodedata => fifo_writenodedata
    procedure, private :: readnodedata => fifo_readnodedata
    procedure, private :: freeresources => fifo_freeresources
    ! Workaround: should be private, but NAG fails to override private routines.
    procedure :: datafromfile => fifo_datafromfile
    procedure :: datatofile => fifo_datatofile
  end type fifobase

  !> Represents one node in the fifo.
  type fifonode
    class(*), pointer :: data => null()
    class(fifonode), pointer :: next => null()
    integer :: filepos = -1
  end type fifonode


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! FIFO Routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of items in the collection.
  !! \param obj  Collection instance.
  !! \return Number of items.
  pure function fifo_size(obj) result(res)
    class(fifobase), intent(in) :: obj
    integer :: res

    res = obj%nitem

  end function fifo_size


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
  subroutine fifo_initswap(self, memorylimit, filename, fileid)
    class(fifobase), intent(inout) :: self
    integer, intent(in) :: memorylimit
    character(*), intent(in) :: filename
    integer, intent(in) :: fileid

    if (self%memorylimit /= -1) then
      stop "FIFO swap can be initialized only once"
    end if
    self%memorylimit = memorylimit
    self%filename = filename
    self%fileid = fileid
  
  end subroutine fifo_initswap


  !> Pushes a pointer to the collection.
  !! \param self  Instance.
  !! \param data  Pointer to the data object.
  subroutine fifo_pushptr(self, data)
    class(fifobase), intent(inout) :: self
    class(*), pointer, intent(in) :: data

    class(fifonode), pointer :: node

    allocate(node)
    node%data => data
    if (.not. associated(self%head)) then
      self%head => node
      self%current => node
    else
      self%tail%next => node
    end if
    self%tail => node
    self%nitem = self%nitem + 1
    self%inmemory = self%inmemory + 1
    if (self%memorylimit /= -1 .and. self%inmemory > self%memorylimit) then
      call self%writenodedata(node)
    end if

  end subroutine fifo_pushptr


  !> Pops a pointer from the collection.
  !! \param self  Instance.
  !! \param data  Pointer to the data object on return.
  subroutine fifo_popptr(self, data)
    class(fifobase), intent(inout) :: self
    class(*), pointer, intent(out) :: data

    class(fifonode), pointer :: node

    if (.not. associated(self%head)) then
      data => null()
      return
    end if

    node => self%head
    self%head => node%next
    if (associated(node, self%current)) then
      self%current => node%next
    end if
    if (associated(node, self%previous)) then
      nullify(self%previous)
    end if
    if (.not. associated(node%data)) then
      call self%readnodedata(node)
    end if
    data => node%data
    deallocate(node)
    self%nitem = self%nitem - 1
    self%inmemory = self%inmemory - 1

  end subroutine fifo_popptr


  !> Gets a copy of a pointer from the collection.
  !! \param self  Instance.
  !! \param data  Pointer to the data object on return.
  subroutine fifo_getptr(self, data)
    class(fifobase), intent(inout) :: self
    class(*), pointer, intent(out) :: data

    if (.not. associated(self%current)) then
      data => null()
      return
    end if

    ! If previous get read something from file, clear the buffer.
    if (associated(self%previous)) then
      if (self%previous%filepos /= -1 .and. associated(self%previous%data)) then
        deallocate(self%previous%data)
        self%inmemory = self%inmemory - 1
      end if
    end if
    
    if (.not. associated(self%current%data)) then
      call self%readnodedata(self%current)
    end if
    data => self%current%data

    self%previous => self%current
    if (associated(self%current%next)) then
      self%current => self%current%next
    else
      self%current => self%head
    end if

  end subroutine fifo_getptr


  !> Restets the collection to it initial (empty) state.
  !! \param self  Instance.
  subroutine fifo_reset(self)
    class(fifobase), intent(inout) :: self

    call self%freeresources()
    self%nitem = 0
    self%inmemory = 0
    self%memorylimit = -1
    nullify(self%head, self%tail, self%current, self%previous)

  end subroutine fifo_reset


  !> Destructor for the class.
  !! \param self Instance.
  subroutine fifo_destruct(self)
    type(fifobase), intent(inout) :: self

    call self%freeresources()

  end subroutine fifo_destruct


  !> Destroys the nodes in the collections and closes open files.
  !! \param self  Instance variable.
  subroutine fifo_freeresources(self)
    class(fifobase), intent(inout) :: self

    class(fifonode), pointer :: node
    logical :: opened

    node => self%head
    do while (associated(node))
      deallocate(node%data)
      self%head => node%next
      deallocate(node)
      node => self%head
    end do
    
    if (self%memorylimit /= -1) then
      inquire(self%fileid, opened=opened)
      if (opened) then
        close(self%fileid, status="delete")
      end if
    end if
    
  end subroutine fifo_freeresources


  !> Writes the data of a node to the disc and deallocates the data object.
  !! \param self  Instance.
  !! \param node  Node with the data that should be stored in a file.
  !! \note This routine invokes the data types write method instead of
  !!   writing the data directly.
  subroutine fifo_writenodedata(self, node)
    class(fifobase), intent(inout) :: self
    class(fifonode), pointer, intent(inout) :: node

    character(10) :: action

    inquire(self%fileid, action=action)
    if (action == "UNDEFINED") then
      ! No saved entries, create new swap file
      open(self%fileid, file=self%filename, access="stream", status="replace",&
          & action="write", form="unformatted", position="rewind")
    elseif (action == "READ") then
      ! Last commmand was pop/get, close file and and reopen in append mode.
      close(self%fileid)
      open(self%fileid, file=self%filename, access="stream", status="old",&
          & action="write", form="unformatted", position="append")
    end if

    inquire(self%fileid, pos=node%filepos)
    call self%datatofile(self%fileid, node%filepos, node%data)
    self%inmemory = self%inmemory - 1

  end subroutine fifo_writenodedata


  !> Reads the data of a node from file and allocates the data object.
  !! \param self  Instance.
  !! \param node  Node with the data that should be read from a file.
  !! \note This routine invokes the data types read method instead of
  !!   reading the data directly.
  subroutine fifo_readnodedata(self, node)
    class(fifobase), intent(inout) :: self
    class(fifonode), pointer, intent(inout) :: node

    character(10) :: action

    inquire(self%fileid, action=action)
    if (action == "WRITE") then
      close(self%fileid)
      open(self%fileid, file=self%filename, access="stream", status="old",&
          & action="read", form="unformatted")
    end if

    call self%datafromfile(self%fileid, node%filepos, node%data)
    self%inmemory = self%inmemory + 1

  end subroutine fifo_readnodedata


  !> Writes the content of a data node to a file.
  !!
  !! \details Extensions of the data object should rewrite it according to
  !! the data they contain.
  !!
  !! \param self  Instance.
  !! \param data  Pointer to a data node, will be deallocated at exit.
  !! \param fileid  File in which the data should be written.
  !! \param filepos  Position in the file, where the data must be written.
  subroutine fifo_datatofile(self, fileid, filepos, data)
    class(fifobase), intent(inout) :: self
    integer, intent(in) :: fileid, filepos
    class(*), intent(inout), pointer :: data
    
    stop "Collection does not support swapping to file."
    
  end subroutine fifo_datatofile
  

  !> Reads the content of a data node from a file.
  !!
  !! \details Extensions of the data object should rewrite it according to
  !! the data they contain.
  !!
  !! \param self  Instance.
  !! \param fileid  File from which the data should be read.
  !! \param filepos  Position in the file, where the data should be read from.
  subroutine fifo_datafromfile(self, fileid, filepos, data)
    class(fifobase), intent(inout) :: self
    integer, intent(in) :: fileid, filepos
    class(*), intent(out), pointer :: data

    stop "Collection does not support swapping to file."
    
  end subroutine fifo_datafromfile


end module fifobase_module
