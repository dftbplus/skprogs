!> Output routines for the sktwocnt code.
module output
  use accuracy
  implicit none
  private

  public :: write_sktables

  ! Maximal angular momentum in the old and the extended old SK file.
  integer, parameter :: LMAX_OLD = 2
  integer, parameter :: LMAX_EXTENDED = 3
  

contains

  subroutine write_sktables(skham, skover)
    real(dp), intent(in) :: skham(:,:), skover(:,:)

    call write_sktable_("at1-at2.ham.dat", skham)
    call write_sktable_("at1-at2.over.dat", skover)

  end subroutine write_sktables
  

  !> Helper routine writing the SK files.
  !! \param fname  File name.
  !! \param sktable  Slater-Koster type integrals (Hamiltonian or overlap).
  subroutine write_sktable_(fname, sktable)
    character(*), intent(in) :: fname
    real(dp), intent(in) :: sktable(:,:)

    integer :: fp, ninteg, nline
    character(11) :: formstr

    ninteg = size(sktable, dim=1)
    print *, "NINTEG:", ninteg
    nline = size(sktable, dim=2)
    write(formstr, "(A,I0,A)") "(", ninteg, "ES21.12)"
    fp = 14
    open(fp, file=fname, status="replace", action="write")
    write(fp, "(I0)") nline
    write(fp, formstr) sktable
    close(fp)
    
  end subroutine write_sktable_

  
end module output
