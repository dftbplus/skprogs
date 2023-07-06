!> Wrapper module for the old F77 Lebedev-Laikov routines.
module common_lebedev_laikov

  use common_accuracy, only : r8

  implicit none
  private

  public :: ld0006, ld0014, ld0026, ld0038, ld0050, ld0074, ld0086, ld0110, ld0146, ld0170, ld0194,&
      & ld0230, ld0266, ld0302, ld0350


  interface

    subroutine ld0006(x, y, z, w, n)
      import :: r8
      real(r8), intent(out) :: x(6), y(6), z(6), w(6)
      integer, intent(out) :: n
    end subroutine ld0006

    subroutine ld0014(x, y, z, w, n)
      import :: r8
      real(r8), intent(out) :: x(14), y(14), z(14), w(14)
      integer, intent(out) :: n
    end subroutine ld0014

    subroutine ld0026(x, y, z, w, n)
      import :: r8
      real(r8), intent(out) :: x(26), y(26), z(26), w(26)
      integer, intent(out) :: n
    end subroutine ld0026

    subroutine ld0038(x, y, z, w, n)
      import :: r8
      real(r8), intent(out) :: x(38), y(38), z(38), w(38)
      integer, intent(out) :: n
    end subroutine ld0038

    subroutine ld0050(x, y, z, w, n)
      import :: r8
      real(r8), intent(out) :: x(50), y(50), z(50), w(50)
      integer, intent(out) :: n
    end subroutine ld0050

    subroutine ld0074(x, y, z, w, n)
      import :: r8
      real(r8), intent(out) :: x(74), y(74), z(74), w(74)
      integer, intent(out) :: n
    end subroutine ld0074

    subroutine ld0086(x, y, z, w, n)
      import :: r8
      real(r8), intent(out) :: x(86), y(86), z(86), w(86)
      integer, intent(out) :: n
    end subroutine ld0086

    subroutine ld0110(x, y, z, w, n)
      import :: r8
      real(r8), intent(out) :: x(110), y(110), z(110), w(110)
      integer, intent(out) :: n
    end subroutine ld0110

    subroutine ld0146(x, y, z, w, n)
      import :: r8
      real(r8), intent(out) :: x(146), y(146), z(146), w(146)
      integer, intent(out) :: n
    end subroutine ld0146

    subroutine ld0170(x, y, z, w, n)
      import :: r8
      real(r8), intent(out) :: x(170), y(170), z(170), w(170)
      integer, intent(out) :: n
    end subroutine ld0170

    subroutine ld0194(x, y, z, w, n)
      import :: r8
      real(r8), intent(out) :: x(194), y(194), z(194), w(194)
      integer, intent(out) :: n
    end subroutine ld0194

    subroutine ld0230(x, y, z, w, n)
      import :: r8
      real(r8), intent(out) :: x(230), y(230), z(230), w(230)
      integer, intent(out) :: n
    end subroutine ld0230

    subroutine ld0266(x, y, z, w, n)
      import :: r8
      real(r8), intent(out) :: x(266), y(266), z(266), w(266)
      integer, intent(out) :: n
    end subroutine ld0266

    subroutine ld0302(x, y, z, w, n)
      import :: r8
      real(r8), intent(out) :: x(302), y(302), z(302), w(302)
      integer, intent(out) :: n
    end subroutine ld0302

    subroutine ld0350(x, y, z, w, n)
      import :: r8
      real(r8), intent(out) :: x(350), y(350), z(350), w(350)
      integer, intent(out) :: n
    end subroutine ld0350

  end interface

end module common_lebedev_laikov
