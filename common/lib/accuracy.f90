!> Some global accuracy settings.
module accuracy
  implicit none
  public

  integer, parameter :: dp = 8
  integer, parameter :: cp = dp  !* precision of the complex data type
  integer, parameter :: sc = 10  !* length of a short string
  integer, parameter :: mc = 50  !* length of a medium length string
  integer, parameter :: lc = 200 !* length of a long string

end module accuracy

