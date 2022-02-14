!> Contains a list of constants for the control of precision of the calculation, both for the
!! fortran numerical model and defaults for the various algorithms in the code.
!! Not all routines use the string length specifications to set their character string lengths.
module common_accuracy

  use, intrinsic :: iso_fortran_env, only : real64

  implicit none
  private

  public :: dp, cp, r8, r4, rdp, rsp, sc, mc, lc

  !> precision of the real data type
  integer, parameter :: dp = real64

  !> precision of the complex data type
  integer, parameter :: cp = dp

  !> parameter for real with 8 byte
  !! (use those only for interfacing F77 routines)
  integer, parameter :: r8 = kind(1.0d0)

  !> parameter for real with 4 byte
  !! (use those only for interfacing F77 routines)
  integer, parameter :: r4 = kind(1.0)

  !> parameter for real with 8 byte
  !! (use those only for interfacing F77 routines)
  integer, parameter :: rdp = kind(1.0d0)

  !> parameter for real with 4 byte
  !! (use those only for interfacing F77 routines)
  integer, parameter :: rsp = kind(1.0)

  !> length of a short string
  integer, parameter :: sc = 10

  !> length of a medium length string
  integer, parameter :: mc = 50

  !> length of a long string
  integer, parameter :: lc = 200

end module common_accuracy

