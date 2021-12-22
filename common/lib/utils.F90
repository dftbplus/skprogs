!> Module that contains various general utilities, based on a code by John E. Pask, LLNL.
module common_utils

  use common_accuracy, only: dp

  implicit none
  private

  public :: upcase, lowcase, whitechar, blank, numstrings, getstring, stop_error, arange, loadtxt,&
      & savetxt, assert, str


  interface str
    module procedure str_int, str_real, str_real_n
  end interface str


contains

  !> Returns string 'strIn' in uppercase.
  function upcase(strIn) result(strOut)

    !> string to convert
    character(len=*), intent(in) :: strIn

    !> uppercase string
    character(len=len(strIn)) :: strOut

    !! auxiliary variables
    integer :: ii, diff

    strOut = strIn
    diff = ichar('A') - ichar('a')

    do ii = 1, len(strOut)
      if ((ichar(strOut(ii:ii)) >= ichar('a')) .and. (ichar(strOut(ii:ii)) <= ichar('z'))) then
        ! if lowercase, make uppercase
        strOut(ii:ii) = char(ichar(strOut(ii:ii)) + diff)
      end if
    end do

  end function upcase


  !> Returns string 'strIn' in lowercase.
  function lowcase(strIn) result(strOut)

    !>
    character(len=*), intent(in) :: strIn

    !>
    character(len=len(strIn)) :: strOut

    !! auxiliary variables
    integer :: ii, diff

    strOut = strIn
    diff = ichar('A') - ichar('a')

    do ii = 1, len(strOut)
      if ((ichar(strOut(ii:ii)) >= ichar('A')) .and. (ichar(strOut(ii:ii)) <= ichar('Z'))) then
        ! if uppercase, make lowercase
        strOut(ii:ii) = char(ichar(strOut(ii:ii)) - diff)
      end if
    end do

  end function lowcase


  !> Identifies 'white' characters.
  !! Returns .true. if char is space (32) or tab (9), .false. otherwise.
  pure function whitechar(char) result(tWhitechar)

    !> character to check
    character, intent(in) :: char

    !> .true. if char is space (32) or tab (9), .false. otherwise
    logical :: tWhitechar

    if ((iachar(char) == 32) .or. (iachar(char) == 9)) then
      tWhitechar = .true.
    else
      tWhitechar = .false.
    end if

  end function whitechar


  !> Returns true, if string contains only white characters.
  pure function blank(str) result(tBlank)

    !> string to check
    character(len=*), intent(in) :: str

    !> true, if string contains only white characters
    logical :: tBlank

    !! iterates over all characters of the string
    integer :: ii

    do ii = 1, len(str)
      if (.not. whitechar(str(ii:ii))) exit
    end do

    tBlank = (ii > len(str))

  end function blank


  !> Returns number of substrings contained in input 'string' delimited by white space.
  pure function numstrings(string) result(nn)

    !> input string to search
    character(len=*), intent(in) :: string

    !> temporary string to facilitate analysis
    character(len=len(string)+2) :: tmp

    !> number of substrings contained in input string
    integer :: nn

    !! iterates over temporary string
    integer :: ii

    tmp = " " // string // " "

    nn = 0
    do ii = 1, len(tmp) - 1
      if (whitechar(tmp(ii:ii)) .and. (.not. whitechar(tmp(ii+1:ii+1)))) nn = nn + 1
    end do

  end function numstrings


  !> Returns first substring 'sub' in string 'str', delimited by white space, starting at index 'is'
  !! in 'str'. If 'sub' is found, 'is' is set to (index of last character of 'sub' in 'str') + 1;
  !! else 'is' is set to 0. If 'is' is out of range on input, routine terminates with 'is' = -1.
  subroutine getstring(str, is, sub)

    !> input string
    character(len=*), intent(in) :: str

    !> on input: starting index for search for 'sub' in 'str'
    !! on output: (index of last character of 'sub' in 'str') + 1
    integer, intent(inout) :: is

    !> first substring in 'str', starting from index 'is'
    character(len=*), intent(out) :: sub

    !> temporary string to facilitate search
    character(len=len(str)+1) :: tmp

    !! true, if previous/current character is 'white'
    logical tPrevwhite, tCurwhite

    !! auxiliary variables
    integer ii, i1, i2

    if (is <= 0 .or. is > len(str)) then
      sub = ""
      is = -1
      return
    end if
    tmp = str // " "
    if (is == 1) then
      tPrevwhite = .true.
    else
      tPrevwhite = whitechar(tmp(is-1:is-1))
    end if
    i1 = 0
    i2 = 0
    do ii = is, len(tmp)
      tCurwhite = whitechar(tmp(ii:ii))
      if (tPrevwhite .and. (.not. tCurwhite)) i1 = ii   ! beginning of substring
      if ((i1 > 0) .and. tCurwhite) then                ! end of substring
        i2 = ii - 1
        exit
      end if
      tPrevwhite = tCurwhite
    end do
    if (i2 > 0) then
      sub = tmp(i1:i2)
      is = i2 + 1
    else
      sub = ""
      is = 0
    end if

  end subroutine getstring


  !> Aborts the program with nonzero exit code.
  !!
  !! The statement "stop msg" will return 0 exit code when compiled using gfortran. stop_error()
  !! uses the statement "stop 1" which returns an exit code 1 and a statement to print the message.
  !!
  !! Example
  !! -------
  !!
  !! call stop_error("Invalid argument")
  subroutine stop_error(msg)

    !> message to print on stdout
    character(len=*) :: msg

    print *, msg
    stop 1

  end subroutine stop_error


  !> Loads a 2D array from a text file.
  !!
  !! Example
  !! -------
  !!
  !! real(dp), allocatable :: data(:,:)
  !! call loadtxt("log.txt", data)  ! 'data' will be automatically allocated
  !!
  !! Where 'log.txt' contains for example::
  !!
  !!      1  2  3
  !!      2  4  6
  !!      8  9 10
  !!     11 12 13
  !!     ...
  !!
  subroutine loadtxt(filename, array)

    !> filename to load the array from
    character(len=*), intent(in) :: filename

    !> data array, will be automatically allocated with the correct dimensions
    real(dp), intent(out), allocatable :: array(:,:)

    !! file identifier
    integer :: fd

    !! dummy character, used to iterate over columns of file
    character :: cc

    !! number of rows and columns of file
    integer :: ncol, nrow

    !! error status
    integer :: iErr

    !! iterates over rows of file
    integer :: ii

    !! true, if last character was 'white'
    logical :: tLastwhite

    !! dummy real variable, used to iterate over rows of file
    real(dp) :: rr

    open(newunit=fd, file=filename, status="old")

    ! determine number of columns
    ncol = 0
    tLastwhite = .true.
    do
      read(fd, '(a)', advance='no', iostat=iErr) cc
      if (iErr /= 0) exit
      if (tLastwhite .and. .not. whitechar(cc)) ncol = ncol + 1
      tLastwhite = whitechar(cc)
    end do

    rewind(fd)

    ! determine number or rows
    nrow = 0
    do
      read(fd, *, iostat=iErr) rr
      if (iErr /= 0) exit
      nrow = nrow + 1
    end do

    rewind(fd)

    allocate(array(nrow, ncol))
    do ii = 1, nrow
      read(fd, *) array(ii, :)
    end do

    close(fd)

  end subroutine loadtxt


  !> Saves a 2D array into a textfile.
  !!
  !! Example
  !! -------
  !!
  !! real(dp) :: data(3, 2)
  !! call savetxt("log.txt", data)
  subroutine savetxt(fname, array)

    !> file to save the array to
    character(len=*), intent(in) :: fname

    !> two-dimensional array to save
    real(dp), intent(in) :: array(:,:)

    !! file identifier
    integer :: fd

    !! iterates over lines of array
    integer :: ii

    open(newunit=fd, file=fname, status="replace")

    do ii = 1, size(array, dim=1)
      write(fd, *) array(ii, :)
    end do

    close(fd)

  end subroutine savetxt


  !> Returns an array = [start, start+dx, start+2*dx, ..., end-dx].
  !!
  !! Example
  !! -------
  !!
  !! real(dp), allocatable :: array(:)
  !! call arange(1, 5, 1, array)   ! array = [1, 2, 3, 4]
  subroutine arange(start, end, step, array)

    !> start, end, step
    real(dp), intent(in) :: start, end, step

    !> generated array
    real(dp), intent(out), allocatable :: array(:)

    !! number of entries in array
    integer :: nn

    !! iterates over entries in array
    integer :: ii

    nn = int((end - start) / step)
    allocate(array(nn))

    do ii = 1, nn
      array(ii) = start + (ii - 1) * step
    end do

  end subroutine arange


  !> Aborts the program, if (condition == .false.).
  !!
  !! Example
  !! -------
  !!
  !! call assert(a == 5)
  subroutine assert(tCondition)

    !> condition
    logical, intent(in) :: tCondition

    if (.not. tCondition) call stop_error("Assert failed.")

  end subroutine assert


  !> Returns the length of the string representation of 'ii'.
  pure function str_int_len(ii) result(sz)

    !> integer
    integer, intent(in) :: ii

    !> length of string representation of 'ii'
    integer :: sz

    !! maximum length of string representation of 'ii'
    integer, parameter :: MAX_STR = 100

    !! string representation of 'ii'
    character(len=MAX_STR) :: str

    ! If 'str' is too short (MAX_STR too small), Fortan will abort with:
    ! "Fortran runtime error: End of record"

    write(str, '(i0)') ii
    sz = len_trim(str)

  end function str_int_len


  !> Converts integer 'ii' to string.
  pure function str_int(ii) result(str)

    !> integer
    integer, intent(in) :: ii

    !> string representation of 'ii'
    character(len=str_int_len(ii)) :: str

    write(str, '(i0)') ii

  end function str_int


  !> Returns the length of the string representation of 'rr'.
  pure function str_real_len(rr, fmt) result(sz)

    !> real
    real(dp), intent(in) :: rr

    !> desired format of string representation of 'rr'
    character(len=*), intent(in) :: fmt

    !> length of the string representation of 'rr'
    integer :: sz

    !! maximum length of string representation of 'rr'
    integer, parameter :: MAX_STR = 100

    !! string representation of 'rr'
    character(len=MAX_STR) :: str

    ! If 's' is too short (MAX_STR too small), Fortan will abort with:
    ! "Fortran runtime error: End of record"

    write(str, fmt) rr
    sz = len_trim(str)

  end function str_real_len


  !> Converts the real number 'rr' to string with 7 decimal digits.
  pure function str_real(rr) result(str)

    !> real
    real(dp), intent(in) :: rr

    !> format string
    character(len=*), parameter :: fmt = "(f0.6)"

    !! string with 7 decimal digits
    character(len=str_real_len(rr, fmt)) :: str

    write(str, fmt) rr

  end function str_real


  !> Converts the real number 'rr' to string with 'nn' decimal digits.
  pure function str_real_n(rr, nn) result(str)

    !> real
    real(dp), intent(in) :: rr

    !> desired number of decimal digits
    integer, intent(in) :: nn

    !> string with 'nn' decimal digits
    character(len=str_real_len(rr, "(f0." // str_int(nn) // ")")) :: str

    write(str, "(f0." // str_int(nn) // ")") rr

  end function str_real_n

end module common_utils
