!> Module that handles input parsing of configuration and raw data.
module input

  use common_accuracy, only : dp
  use gridorbital, only : TGridorb2_init
  use twocnt, only : TTwocntInp, TAtomdata
  use xcfunctionals, only : xcFunctional

  implicit none
  private

  public :: readInput

  !> maximum line length of sktwocnt.in file
  integer, parameter :: maxlen = 1024

  !> expected line format when reading sktwocnt.in file
  character(len=*), parameter :: lineformat = "(A1024)"

  !> comment string
  character, parameter :: comment = "#"


contains

  !> Reads and extracts relevant information from 'sktwocnt.in' file.
  subroutine readInput(inp, fname)

    !> instance of parsed input for twocnt
    type(TTwocntInp), intent(out) :: inp

    !> filename
    character(len=*), intent(in) :: fname

    !! file identifier
    integer :: fp

    !! current line index
    integer :: iLine

    !! character buffer
    character(len=maxlen) :: line, buffer1, buffer2

    !! error status
    integer :: iErr

    !! xc-functional type
    !! (1: LDA-PW91, 2: GGA-PBE96, 3: GGA-BLYP, 4: LCY-PBE96, 5: LCY-BNL, 6: PBE0, 7: B3LYP,
    !! 8: CAMY-B3LYP, 9: CAMY-PBEh, 10: TPSS, 11: SCAN, 12: r2SCAN, 13: r4SCAN, 14: TASK)
    integer :: iXC

    !! potential data columns, summed up in order to receive the total atomic potential
    integer, allocatable :: potcomps(:)

    !! true, if radial grid-orbital 1st/2nd derivative shall be read
    logical :: tReadRadDerivs

    inp%tMGGA = .false.
    inp%tLC = .false.
    inp%tCam = .false.
    inp%tGlobalHybrid = .false.

    open(newunit=fp, file=fname, form="formatted", action="read")
    ! general part
    iLine = 0

    call nextline_(fp, iLine, line)
    read(line, *, iostat=iErr) buffer1, buffer2, iXC
    call checkerror_(fname, line, iLine, iErr)
    if ((buffer1 /= "hetero") .and. (buffer1 /= "homo")) then
      call error_("Wrong interaction (must be 'hetero' or 'homo')!", fname, line, iLine)
    end if
    inp%tHetero = (buffer1 == "hetero")

    select case (buffer2)
    case("potential")
      inp%tDensitySuperpos = .false.
    case("density")
      inp%tDensitySuperpos = .true.
    case default
      call error_("Wrong superposition mode (must be 'potential' or 'density')!", fname, line,&
          & iline)
    end select

    select case (iXC)
    case(xcFunctional%LDA_PW91)
      ! LDA-PW91
    case(xcFunctional%GGA_PBE96)
      ! GGA-PBE96
    case(xcFunctional%GGA_BLYP)
      ! GGA-BLYP
    case(xcFunctional%LCY_PBE96)
      ! LCY-PBE96 (purely long-range corrected)
      inp%tLC = .true.
    case(xcFunctional%LCY_BNL)
      ! LCY-BNL (purely long-range corrected)
      inp%tLC = .true.
    case(xcFunctional%HYB_PBE0)
      ! PBE0 (global hybrid)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%HYB_B3LYP)
      ! B3LYP (global hybrid)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%CAMY_B3LYP)
      ! CAMY-B3LYP (general CAM form)
      inp%tCam = .true.
    case(xcFunctional%CAMY_PBEh)
      ! CAMY-PBEh (general CAM form)
      inp%tCam = .true.
    case(xcFunctional%MGGA_TPSS, xcFunctional%MGGA_SCAN, xcFunctional%MGGA_r2SCAN,&
        & xcFunctional%MGGA_r4SCAN, xcFunctional%MGGA_TASK, xcFunctional%MGGA_TASK_CC)
      inp%tMGGA = .true.
    case default
      call error_("Unknown exchange-correlation functional!", fname, line, iline)
    end select
    inp%iXC = iXC

    if (inp%iXC == xcFunctional%HYB_B3LYP) then
      ! 20% fraction of HFX hard-coded at the moment
      inp%camAlpha = 0.2_dp
      inp%camBeta = 0.0_dp
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) inp%nRadial, inp%nAngular, inp%ll_max, inp%rm
      call checkerror_(fname, line, iLine, iErr)
    elseif (inp%iXC == xcFunctional%HYB_PBE0) then
      inp%camBeta = 0.0_dp
      call nextline_(fp, iLine, line)
      ! currently only HYB-PBE0 does support arbitrary HFX portions (HYB-B3LYP does not)
      read(line, *, iostat=iErr) inp%camAlpha
      call checkerror_(fname, line, iLine, iErr)
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) inp%nRadial, inp%nAngular, inp%ll_max, inp%rm
      call checkerror_(fname, line, iLine, iErr)
    elseif (inp%tLC) then
      inp%camAlpha = 0.0_dp
      inp%camBeta = 1.0_dp
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) inp%omega
      if (inp%omega < 1.0e-08_dp) then
        write(*,'(a)') 'Chosen omega too small!'
        stop
      end if
      call checkerror_(fname, line, iLine, iErr)
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) inp%nRadial, inp%nAngular, inp%ll_max, inp%rm
      call checkerror_(fname, line, iLine, iErr)
    elseif (inp%tCam) then
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) inp%omega, inp%camAlpha, inp%camBeta
      if (inp%omega < 1.0e-08_dp) then
        write(*,'(a)') 'Chosen omega too small!'
        stop
      end if
      call checkerror_(fname, line, iLine, iErr)
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) inp%nRadial, inp%nAngular, inp%ll_max, inp%rm
      call checkerror_(fname, line, iLine, iErr)
    end if

    call nextline_(fp, iLine, line)
    read(line, *, iostat=iErr) inp%r0, inp%dr, inp%epsilon, inp%maxdist
    call checkerror_(fname, line, iLine, iErr)

    call nextline_(fp, iLine, line)
    read(line, *, iostat=iErr) inp%ninteg1, inp%ninteg2
    call checkerror_(fname, line, iLine, iErr)

    if (inp%tDensitySuperpos) then
      allocate(potcomps(2))
      potcomps = [2, 3]
    else
      allocate(potcomps(3))
      potcomps = [2, 3, 4]
    end if
    tReadRadDerivs = .not. inp%tHetero

    call readatom_(fname, fp, iLine, potcomps, inp%tDensitySuperpos, tReadRadDerivs,&
        & (inp%tGlobalHybrid .or. inp%tLC .or. inp%tCam), inp%tMGGA, inp%atom1)
    if (inp%tHetero) then
      call readatom_(fname, fp, iLine, potcomps, inp%tDensitySuperpos, .true., (inp%tGlobalHybrid&
          & .or. inp%tLC .or. inp%tCam), inp%tMGGA, inp%atom2)
    end if

    close(fp)

  end subroutine readInput


  !> Fills TAtomdata instance based on slateratom's output.
  subroutine readatom_(fname, fp, iLine, potcomps, tDensitySuperpos, tReadRadDerivs, tNonLocal,&
      & tMGGA, atom)

    !> filename
    character(len=*), intent(in) :: fname

    !> file identifier
    integer, intent(in) :: fp

    !> current line index
    integer, intent(inout) :: iLine

    !> potential data columns, summed up in order to receive the total atomic potential
    integer, intent(in) :: potcomps(:)

    !> true, if density superposition is requested, otherwise potential superposition is applied
    logical, intent(in) :: tDensitySuperpos

    !> true, if radial grid-orbital 1st/2nd derivative shall be read
    logical, intent(in) :: tReadRadDerivs

    !! true, there are non-local exchange contributions to calculate
    logical, intent(in) :: tNonLocal

    !! true, if kinetic energy density (tau) needs to be read from density file
    logical, intent(in) :: tMGGA

    !> atomic properties instance
    type(TAtomdata), intent(out) :: atom

    !! character buffer
    character(maxlen) :: line, buffer

    !! temporary storage for checking radial wavefunction sign
    real(dp) :: vals(2)

    !! temporarily stores atomic wavefunction and potential
    real(dp), allocatable :: data(:,:), potval(:)

    !! error status
    integer :: iErr

    !! auxiliary variables
    integer :: ii, imax

    call nextline_(fp, iLine, line)
    if (tNonLocal) then
      read(line, *, iostat=iErr) atom%nBasis, atom%nCore
    else
      read(line, *, iostat=iErr) atom%nBasis
    end if
    call checkerror_(fname, line, iLine, iErr)

    allocate(atom%angmoms(atom%nBasis))
    allocate(atom%rad(atom%nBasis))
    if (tReadRadDerivs) then
      allocate(atom%drad(atom%nBasis))
      allocate(atom%ddrad(atom%nBasis))
    end if

    do ii = 1, atom%nBasis
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) buffer, atom%angmoms(ii)
      call checkerror_(fname, line, iLine, iErr)
      if (tReadRadDerivs) then
        call readdata_(buffer, [1, 3, 4, 5], data)
        call TGridorb2_init(atom%rad(ii), data(:, 1), data(:, 2))
        call TGridorb2_init(atom%drad(ii), data(:, 1), data(:, 3))
        call TGridorb2_init(atom%ddrad(ii), data(:, 1), data(:, 4))
      else
        call readdata_(buffer, [1, 3], data)
        call TGridorb2_init(atom%rad(ii), data(:, 1), data(:, 2))
      end if
      ! check if wave function follows the sign convention
      ! (positive where abs(r * R(r)) has its maximum)
      imax = maxloc(abs(data(:, 1) * data(:, 2)), dim=1)
      if (data(imax, 2) < 0.0_dp) then
        write(*, "(A,F5.2,A)") "Wave function negative at the maximum of radial probability&
            & (r =", data(imax, 1), " Bohr)"
        write(*, "(A)") "Please change the sign of the wave function (and of its derivatives)!"
        write(*, "(A,A,A)") "File: '", trim(buffer), "'"
        stop
      end if
    end do
    call checkangmoms_(atom%angmoms)

    ! read core orbitals, LC functionals only
    if (tNonLocal) then
      allocate(atom%coreAngmoms(atom%nCore))
      allocate(atom%coreOcc(atom%nCore))
      allocate(atom%coreRad(atom%nCore))

      do ii = 1, atom%nCore
        call nextline_(fp, iLine, line)
        read(line, *, iostat=iErr) buffer, atom%coreAngmoms(ii), atom%coreOcc(ii)
        call checkerror_(fname, line, iLine, iErr)
        call readdata_(buffer, [1, 3], data)
        call TGridorb2_init(atom%coreRad(ii), data(:, 1), data(:, 2))
        vals = atom%coreRad(ii)%getValue([0.01_dp, 0.02_dp])
        if ((vals(2) - vals(1)) < 0.0_dp) then
          call atom%coreRad(ii)%rescale(-1.0_dp)
        end if
      end do

    end if

    call nextline_(fp, iLine, line)
    read(line, *, iostat=iErr) buffer
    call checkerror_(fname, line, iLine, iErr)
    call readdata_(buffer, [1, 3, 4, 5], data)
    allocate(potval(size(data, dim=1)))
    potval(:) = 0.0_dp
    do ii = 1, size(potcomps)
      potval(:) = potval + data(:, potcomps(ii))
    end do
    call TGridorb2_init(atom%pot, data(:, 1), potval)

    call nextline_(fp, iLine, line)
    read(line, *, iostat=iErr) buffer
    call checkerror_(fname, line, iLine, iErr)
    if (tDensitySuperpos) then
      if (tMGGA) then
        call readdata_(buffer, [1, 3, 4, 5, 6], data)
        call TGridorb2_init(atom%tau, data(:, 1), data(:, 5))
      else
        call readdata_(buffer, [1, 3, 4, 5], data)
      end if
      call TGridorb2_init(atom%rho, data(:, 1), data(:, 2))
      call TGridorb2_init(atom%drho, data(:, 1), data(:, 3))
      call TGridorb2_init(atom%ddrho, data(:, 1), data(:, 4))
    else
      if (trim(line) /= "noread") then
        write(*, "(A,I0,A)") "Line ", iLine, " ignored since density is not needed."
      end if
    end if

  end subroutine readatom_


  !> Reads desired colums of a data file.
  subroutine readdata_(fname, cols, data)

    !> filename
    character(len=*), intent(in) :: fname

    !> desired columns to read from file
    integer, intent(in) :: cols(:)

    !> obtained data on grid with nGrid entries
    real(dp), intent(out), allocatable :: data(:,:)

    !! temporarily stores all columns of a single line in file
    real(dp), allocatable :: tmp(:)

    !! character buffer for current line of file
    character(maxlen) :: line

    !! number of grid points stored in file
    integer :: nGrid

    !! error status
    integer :: iErr

    !! current line
    integer :: iLine

    !! file identifier
    integer :: fp

    !! auxiliary variable
    integer :: ii

    iLine = 1

    allocate(tmp(maxval(cols)))

    open(newunit=fp, file=fname, action="read", form="formatted")

    call nextline_(fp, iLine, line)
    read(line, *, iostat=iErr) nGrid
    call checkerror_(fname, line, iLine, iErr)

    allocate(data(nGrid, size(cols)))
    do ii = 1, nGrid
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) tmp(:)
      call checkerror_(fname, line, iLine, iErr)
      data(ii, :) = tmp(cols)
    end do

    close(fp)

  end subroutine readdata_


  !> Iterates through lines of a file, while respecting an user-def. comment string and empty lines.
  subroutine nextline_(fp, iLine, line)

    !> file identifier
    integer, intent(in) :: fp

    !> current line of the file
    integer, intent(inout) :: iLine

    !> line buffer
    character(maxlen), intent(out) :: line

    !! position of comment string in line if present, otherwise zero
    integer :: ii

    !! temporarily stores an entire line
    character(maxlen) :: buffer

    do while (.true.)
      iLine = iLine + 1
      read(fp, lineformat) buffer
      ii = index(buffer, comment)
      if (ii == 0) then
        line = adjustl(buffer)
      else
        line = adjustl(buffer(1:ii - 1))
      end if
      if (len_trim(line) > 0) exit
    end do

  end subroutine nextline_


  !> Checks range of angular momenta w.r.t. program compatibility.
  subroutine checkangmoms_(angmoms)

    !> angular momenta
    integer, intent(in) :: angmoms(:)

    if (maxval(angmoms) > 4) then
      write(*,*) "Only angular momentum up to 'f' is allowed."
      stop
    end if

  end subroutine checkangmoms_


  !> Error handling.
  subroutine checkerror_(fname, line, iLine, iErr)

    !> filename
    character(len=*), intent(in) :: fname

    !> content of current line
    character(len=*), intent(in) :: line

    !> current line of parsed file
    integer, intent(in) :: iLine

    !> error status
    integer, intent(in) :: iErr

    if (iErr /= 0) then
      call error_("Bad syntax", fname, line, iLine)
    end if

  end subroutine checkerror_


  !> Throws error message.
  subroutine error_(txt, fname, line, iLine)

    !> user-specified error message
    character(len=*), intent(in) :: txt

    !> filename
    character(len=*), intent(in) :: fname

    !> content of erroneous line
    character(len=*), intent(in) :: line

    !> index of erroneous line
    integer, intent(in) :: iLine

    write(*, "(A,A)") "!!! Parsing error: ", txt
    write(*, "(2X,A,A)") "File: ", trim(fname)
    write(*, "(2X,A,I0)") "Line number: ", iLine
    write(*, "(2X,A,A,A)") "Line: '", trim(line), "'"

    stop

  end subroutine error_

end module input
