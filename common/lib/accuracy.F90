!> Contains a list of constants for the control of precision of the calculation, both for the
!! fortran numerical model and defaults for the various algorithms in the code.
!! Not all routines use the string length specifications to set their character string lengths.
module common_accuracy

  implicit none
  private

  public :: dp, cp, sc, mc, lc

  !> precision of the real data type
  integer, parameter :: dp = 8

  !> precision of the complex data type
  integer, parameter :: cp = dp

  !> length of a short string
  integer, parameter :: sc = 10

  !> length of a medium length string
  integer, parameter :: mc = 50

  !> length of a long string
  integer, parameter :: lc = 200

end module common_accuracy

