module constants
  use accuracy
  implicit none

  real(dp), parameter :: pi = 3.14159265358979323846_dp
  real(dp), parameter :: r_Bohr   = 0.529177249_dp      !< Bohr radius (&#197;)
  real(dp), parameter :: Bohr__AA = r_Bohr              !< Bohr->Angstrom
  real(dp), parameter :: AA__Bohr = 1.0_dp / Bohr__AA   !< Angstrom->Bohr
  real(dp), parameter :: Hartree = 27.2113845_dp        !< H -> eV (CODATA)
  real(dp), parameter :: sol = 137.0359997_dp           !< Speed of Light a.u.  
end module constants
