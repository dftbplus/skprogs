!> Module that builds up various supervectors.
module confinement

  use common_accuracy, only : dp
  use utilities, only : fak
  use core_overlap, only : v

  implicit none
  private

  public :: TConf, confType


  !> Confinement potential structure.
  type :: TConf

    !> type of confinement potential (0: none, 1: power, 2: Woods-Saxon)
    !! at the moment different shells are compressed by the same type of potential
    integer :: type = 0

    !> confinement radii (power compression)
    real(dp) :: r0(0:4)

    !> power of confinement (power compression)
    real(dp) :: power(0:4)

    !> onset radii (Woods-Saxon compression)
    real(dp) :: onset(0:4)

    !> cutoff of confinement (Woods-Saxon compression)
    real(dp) :: cutoff(0:4)

    !> potential well height of confinement (Woods-Saxon compression)
    real(dp) :: vmax(0:4)

  contains

    procedure :: getConfPower => TConf_getConfPower

  end type TConf


  !> Enumerator for type of confinement potential.
  type :: TConfEnum

    !> no compression
    integer :: none = 0

    !> power compression
    integer :: power = 1

    !> Woods-Saxon compression
    integer :: ws = 2

  end type TConfEnum

  !> Container for enumerated types of confinement potentials.
  type(TConfEnum), parameter :: confType = TConfEnum()


contains

  !> Calculates analytic matrix elements of confining potential.
  !! No checking for power, e.g. power==0 or power<0 etc. !
  subroutine TConf_getConfPower(this, max_l, num_alpha, alpha, poly_order, vconf)

    !> instance
    class(TConf), intent(in) :: this

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> confinement supervector
    real(dp), intent(out) :: vconf(0:,:,:)

    !! temporary storage
    real(dp) :: alpha1

    !! auxiliary variables
    integer :: ii, jj, kk, ll, mm, nn, oo, nlp, nlq

    vconf(:,:,:) = 0.0_dp

    do ii = 0, max_l
      if (this%power(ii) > 1.0e-06_dp) then
        nn = 0
        do jj = 1, num_alpha(ii)
          do ll = 1, poly_order(ii)
            nn = nn + 1
            oo = 0
            nlp = ll + ii
            do kk = 1, num_alpha(ii)
              alpha1 = 0.5_dp * (alpha(ii, jj) + alpha(ii, kk))
              do mm = 1, poly_order(ii)
                oo = oo + 1
                nlq = mm + ii
                vconf(ii, nn, oo) = 1.0_dp / sqrt(v(alpha(ii, jj), 2 * nlp)&
                    & * v(alpha(ii, kk), 2 * nlq)) / (this%r0(ii) * 2.0_dp)**this%power(ii)&
                    & * v(alpha1, nlp + nlq + this%power(ii))
              end do
            end do
          end do
        end do
      end if
    end do

  end subroutine TConf_getConfPower

end module confinement
