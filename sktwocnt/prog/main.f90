program main

  use common_accuracy, only: dp
  use input
  use twocnt
  use output
  use cmdargs
  implicit none

  type(twocnt_in) :: inp
  type(integmap) :: imap
  real(dp), allocatable :: skham(:, :), skover(:, :)

  call parse_command_arguments()
  call readinput(inp, "sktwocnt.in")
  write(*, "(A)") "Input done."
  call get_twocenter_integrals(inp, imap, skham, skover)
  write(*, "(A)") "Twocnt done."
  call write_sktables(skham, skover)

end program main
