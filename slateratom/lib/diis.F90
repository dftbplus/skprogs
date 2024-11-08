#:include 'common.fypp'

module diis
  use common_accuracy, only : dp

  implicit none
    
  private
  ! public :: TDIISMixer, TDIISMixer_init, TDIISMixer_reset, TDIISMixer_mix

  !> DIIS mixer class
  type, public :: TDIISMixer
    private
    
    contains

end type TDIISMixer

contains    

end module diis