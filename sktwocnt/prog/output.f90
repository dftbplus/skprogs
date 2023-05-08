!> Output routines for the sktwocnt code.
module output

  use common_accuracy, only : dp

  implicit none
  private

  public :: write_sktables


contains

  !> Writes tabulated Hamiltonian and overlap matrix to file.
  subroutine write_sktables(skham, skover)

    !> Hamiltonian and overlap matrix
    real(dp), intent(in) :: skham(:,:), skover(:,:)

    call write_sktable_("at1-at2.ham.dat", skham)
    call write_sktable_("at1-at2.over.dat", skover)

  end subroutine write_sktables


  !> Helper routine writing the SK files.
  subroutine write_sktable_(fname, sktable)

    !> file name
    character(len=*), intent(in) :: fname

    !> Slater-Koster type integrals (Hamiltonian or overlap)
    real(dp), intent(in) :: sktable(:,:)

    !! file identifier
    integer :: fp

    !! number of all nonzero two-center integrals
    integer :: ninteg

    !! number of dimer distances, i.e. lines of written file
    integer :: nline

    !! formatting string
    character(len=11) :: formstr

    ninteg = size(sktable, dim=1)
    nline = size(sktable, dim=2)
    write(formstr, "(A,I0,A)") "(", ninteg, "ES21.12)"

    open(newunit=fp, file=fname, status="replace", action="write")
    write(fp, "(I0)") nline
    write(fp, formstr) sktable
    close(fp)

  end subroutine write_sktable_

end module output
