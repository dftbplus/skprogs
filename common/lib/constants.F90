!> Contains a list of physical constants for the code.
module common_constants

  use common_accuracy, only : dp

  implicit none
  private

  public :: pi, Bohr__AA, AA__Bohr, Hartree__eV, eV__Hartree, cc

  !> pi
  real(dp), parameter :: pi = 3.14159265358979323846_dp

  !> Bohr->Angstrom
  real(dp), parameter :: Bohr__AA = 0.529177249_dp

  !> Angstrom->Bohr
  real(dp), parameter :: AA__Bohr = 1.0_dp / Bohr__AA

  !> Hartre -> eV
  real(dp), parameter :: Hartree__eV = 27.2113845_dp

  !> eV->Hartree
  real(dp), parameter :: eV__Hartree = 1.0_dp / Hartree__eV

  !> speed of light
  real(dp), parameter :: cc = 137.0359997_dp

end module common_constants
