module input
  use accuracy
  use gridorbital
  use twocnt, only: twocnt_in, atomdata
  implicit none
  private

  public :: readinput

  integer, parameter :: maxlen = 1024
  character(len=*), parameter :: lineformat = "(A1024)"
  character, parameter :: comment = "#"

contains

  subroutine readinput(inp, inputfile)
    type(twocnt_in), intent(out) :: inp
    character(*), intent(in) :: inputfile

    integer :: fp, iline
    character(maxlen) :: line, buffer1, buffer2
    integer :: iostat
    integer, allocatable :: potcomps(:)
    logical :: readradderivs

    fp = 14
    open(fp, file=inputfile, form="formatted", action="read")
    !! General part
    iline = 0
    call nextline_(fp, iline, line)
    read(line, *, iostat=iostat) buffer1, buffer2
    call checkerror_(inputfile, line, iline, iostat)
    if (buffer1 /= "hetero" .and. buffer1 /= "homo") then
      call error_("Wrong interaction (must be hetero or homo)", inputfile, &
          &line, iline)
    end if
    inp%hetero = (buffer1 == "hetero")
    select case (buffer2)
    case("potential")
      inp%density = .false.
      inp%ixc = 0
    case("density_lda")
      inp%density = .true.
      inp%ixc = 1
    case("density_pbe")
      inp%density = .true.
      inp%ixc = 2
    case default
      call error_("Wrong superposition mode (must be potential, density_lda &
          &or density_pbe", inputfile, line, iline)
    end select
    
    call nextline_(fp, iline, line)
    read(line, *, iostat=iostat) inp%r0, inp%dr, inp%epsilon, inp%maxdist
    call checkerror_(inputfile, line, iline, iostat)

    call nextline_(fp, iline, line)
    read(line, *, iostat=iostat) inp%ninteg1, inp%ninteg2
    call checkerror_(inputfile, line, iline, iostat)

    if (inp%density) then
      allocate(potcomps(2))
      potcomps = [ 2, 3 ]
    else
      allocate(potcomps(3))
      potcomps = [ 2, 3, 4 ]
    end if
    readradderivs = .not. inp%hetero
    call readatom_(inputfile, fp, iline, potcomps, inp%density, readradderivs, &
        &inp%atom1)
    if (inp%hetero) then
      call readatom_(inputfile, fp, iline, potcomps, inp%density, .true., &
          &inp%atom2)
    end if
    
    close(fp)
    
  end subroutine readinput



  subroutine readatom_(fname, fp, iline, potcomps, density, radderivs, atom)
    character(*), intent(in) :: fname
    integer, intent(in) :: fp
    integer, intent(inout) :: iline
    integer, intent(in) :: potcomps(:)
    logical, intent(in) :: density, radderivs
    type(atomdata), intent(out) :: atom

    character(maxlen) :: line, buffer
    real(dp), allocatable :: data(:,:), potval(:)
    real(dp) :: vals(1)
    integer :: ii, iostat, imax


    call nextline_(fp, iline, line)
    read(line, *, iostat=iostat) atom%nbasis
    call checkerror_(fname, line, iline, iostat)
    
    allocate(atom%angmoms(atom%nbasis))
    allocate(atom%rad(atom%nbasis))
    if (radderivs) then
      allocate(atom%drad(atom%nbasis))
      allocate(atom%ddrad(atom%nbasis))
    end if
    do ii = 1, atom%nbasis
      call nextline_(fp, iline, line)
      read(line, *, iostat=iostat) buffer, atom%angmoms(ii)
      call checkerror_(fname, line, iline, iostat)
      if (radderivs) then
        call readdata_(buffer, [ 1, 3, 4, 5 ], data)
        call init(atom%rad(ii), data(:,1), data(:,2))
        call init(atom%drad(ii), data(:,1), data(:,3))
        call init(atom%ddrad(ii), data(:,1), data(:,4))
      else
        call readdata_(buffer, [ 1, 3 ], data)
        call init(atom%rad(ii), data(:,1), data(:,2))
      end if
      ! Check if wave function follows the sign convention
      ! (positive where abs(r * R(r)) has its maximum)
      imax = maxloc(abs(data(:,1) * data(:,2)), dim=1)
      if (data(imax,2) < 0.0_dp) then
          write(*, "(A,F5.2,A)") "Wave function negative at the maximum of&
              & radial probability (r =", data(imax,1), " Bohr)"
          write(*, "(A)") "Please change the sign of the wave function (and of&
              & its derivatives)!"
          write(*, "(A,A,A)") "File: '", trim(buffer), "'"
          stop
        end if
    end do
    call checkangmoms_(atom%angmoms)

    call nextline_(fp, iline, line)
    read(line, *, iostat=iostat) buffer
    call checkerror_(fname, line, iline, iostat)
    call readdata_(buffer, [ 1, 3, 4, 5 ], data)
    allocate(potval(size(data, dim=1)))
    potval = 0.0_dp
    do ii = 1, size(potcomps)
      potval = potval + data(:,potcomps(ii))
    end do
    call init(atom%pot, data(:,1), potval)
    
    call nextline_(fp, iline, line)
    read(line, *, iostat=iostat) buffer
    call checkerror_(fname, line, iline, iostat)
    if (density) then
      call readdata_(buffer, [ 1, 3, 4, 5 ], data)
      call init(atom%rho, data(:,1), data(:,2))
      call init(atom%drho, data(:,1), data(:,3))
      call init(atom%ddrho, data(:,1), data(:,4))
    else
      if (trim(line) /= "noread") then
        write(*,"(A,I0,A)") "Line ", iline, &
            &" ignored since density is not needed."
      end if
    end if

  end subroutine readatom_
  
  
  subroutine readdata_(fname, cols, data)
    character(*), intent(in) :: fname
    integer, intent(in) :: cols(:)
    real(dp), allocatable, intent(out) :: data(:,:)

    real(dp), allocatable :: tmp(:)
    character(maxlen) :: line
    integer :: ngrid, ii, fp, iline, iostat

    fp = 12
    allocate(tmp(maxval(cols)))
    iline = 1
    open(fp, file=fname, action="read", form="formatted")
    call nextline_(fp, iline, line)
    read(line, *, iostat=iostat) ngrid
    call checkerror_(fname, line, iline, iostat)
    allocate(data(ngrid, size(cols)))
    do ii = 1, ngrid
      call nextline_(fp, iline, line)
      read(line, *, iostat=iostat) tmp(:)
      call checkerror_(fname, line, iline, iostat)
      data(ii,:) = tmp(cols)
    end do
    close(fp)
    deallocate(tmp)
    
  end subroutine readdata_

  

  subroutine nextline_(fp, iline, line)
    integer, intent(in) :: fp
    integer, intent(inout) :: iline
    character(maxlen), intent(out) :: line

    integer :: ii
    character(maxlen) :: buffer

    do while (.true.)
      iline = iline + 1
      read(fp, lineformat) buffer
      ii = index(buffer, comment)
      if (ii == 0) then
        line = adjustl(buffer)
      else
        line = adjustl(buffer(1:ii-1))
      end if
      if (len_trim(line) > 0) then
        exit
      end if
    end do

  end subroutine nextline_


  
  subroutine checkangmoms_(angmoms)
    integer, intent(in) :: angmoms(:)

    integer :: ii

    if (maxval(angmoms) > 4) then
      write(*,*) "Only angular momentum up to 'f' is allowed."
      stop
    end if

  end subroutine checkangmoms_


  subroutine checkerror_(fname, line, iline, iostat)
    character(*), intent(in) :: fname, line
    integer, intent(in) :: iline, iostat

    if (iostat /= 0) then
      call error_("Bad syntax", fname, line, iline)
    end if
    
  end subroutine checkerror_


  subroutine error_(txt, fname, line, iline)
    character(*), intent(in) :: txt, fname, line
    integer, intent(in) :: iline

    write(*,"(A,A)") "!!! Parsing error: ", txt
    write(*,"(2X,A,A)") "File: ", trim(fname)
    write(*,"(2X,A,I0)") "Line number: ", iline
    write(*,"(2X,A,A,A)") "Line: '", trim(line), "'"
    stop
  end subroutine error_
  
    
    
end module input
