!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Contains computer environment settings
module common_environment
  use common_globalenv, only : shutdown, stdOut
#:if WITH_MPI
  use common_globalenv, only : globalMpiComm
  use common_mpienv, only : TMpiEnv, TMpiEnv_init, TMpiEnv_final
#:endif
  implicit none

  private
  public :: TEnvironment, TEnvironment_init


  !> Contains environment settings.
  type :: TEnvironment
    private

    !> Whether this process is the lead?
    logical, public :: tGlobalLead = .true.

    !> Nr. of groups in the system
    integer, public :: nGroup = 1

    !> Id of current group (starts with 0)
    integer, public :: myGroup = 0

  #:if WITH_MPI
    !> Global mpi settings
    type(TMpiEnv), public :: mpi

    !> Whether MPI environment had been initialised
    logical :: mpiInitialised = .false.
  #:endif

  contains
    procedure :: destruct => TEnvironment_destruct
    procedure :: shutdown => TEnvironment_shutdown

  #:if WITH_MPI
    procedure :: initMpi => TEnvironment_initMpi
  #:endif

  end type TEnvironment


contains

  !> Returns an initialized instance.
  subroutine TEnvironment_init(this)

    !> Instance
    type(TEnvironment), intent(out) :: this

    continue

  end subroutine TEnvironment_init


  !> Finalizes the environment.
  subroutine TEnvironment_destruct(this)

    !> Instance
    class(TEnvironment), intent(inout) :: this

    #:if WITH_MPI
      if (this%mpiInitialised) then
        call TMpiEnv_final(this%mpi)
        this%mpiInitialised = .false.
      end if
    #:endif

    flush(stdOut)

  end subroutine TEnvironment_destruct


  !> Gracefully cleans up and shuts down.
  !!
  !! Note: This routine must be collectively called by all processes.
  subroutine TEnvironment_shutdown(this)

    !> Instance
    class(TEnvironment), intent(inout) :: this

    call this%destruct()
    call shutdown()

  end subroutine TEnvironment_shutdown


#:if WITH_MPI
  !> Initializes MPI environment.
  subroutine TEnvironment_initMpi(this, nGroup)

    !> Instance
    class(TEnvironment), intent(inout) :: this

    !> Number of process groups to create
    integer, intent(in) :: nGroup

    ! MPI settings
    call TMpiEnv_init(this%mpi, globalMpiComm=globalMpiComm, nGroup=nGroup)
    this%tGlobalLead = this%mpi%tGlobalLead
    this%nGroup = this%mpi%nGroup
    this%myGroup = this%mpi%myGroup
    this%mpiInitialised = .true.

  end subroutine TEnvironment_initMpi
#:endif

end module common_environment
