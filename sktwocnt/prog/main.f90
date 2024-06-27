!> Program to calculate two-center integrals of Slater-Koster tables.
program main

  use common_accuracy, only : dp
  use input, only : readInput
  use twocnt, only : TTwocntInp, TIntegMap, get_twocenter_integrals
  use output, only : write_sktables
  use cmdargs, only : parse_command_arguments

  implicit none

  !> representation of parsed input for sktwocnt.
  type(TTwocntInp) :: inp

  !> specifies type for mapping integrals.
  type(TIntegMap) :: imap

  !> resulting Hamiltonian and overlap matrices
  real(dp), allocatable :: skham(:,:), skover(:,:)

  call parse_command_arguments()

  call readInput(inp, "sktwocnt.in")
  write(*, "(A)") "Input done."
  
  call get_twocenter_integrals(inp, imap, skham, skover)
  write(*, "(A)") "Twocnt done."

  call write_sktables(skham, skover)

end program main
