!> Contains routines to write out various data structures in a comprehensive tagged format.
module common_taggedout

  use common_accuracy, only : dp

  implicit none
  private

  public :: TTaggedwriter, TTaggedwriter_init, writetag, lenLabel


  !> Length of permissible tag labels. Tag names should be shorter than lenLabel!
  integer, parameter :: lenLabel = 20

  !> Max length of the format strings for individual items
  integer, parameter :: lenFormStr = 20


  !> Tag format writer type.
  type :: TTaggedwriter
    character(lenFormStr) :: formReal
    character(lenFormStr) :: formCmplx
    character(lenFormStr) :: formInt
    character(lenFormStr) :: formLogical
  end type TTaggedwriter


  !> Writes objects in a standardized tagged form to a given file.
  interface writetag
    module procedure TTaggedwriter_real0
    module procedure TTaggedwriter_real1
    module procedure TTaggedwriter_real2
    module procedure TTaggedwriter_real3
    module procedure TTaggedwriter_real4
    module procedure TTaggedwriter_cplx0
    module procedure TTaggedwriter_cplx1
    module procedure TTaggedwriter_cplx2
    module procedure TTaggedwriter_cplx3
    module procedure TTaggedwriter_cplx4
    module procedure TTaggedwriter_int0
    module procedure TTaggedwriter_int1
    module procedure TTaggedwriter_int2
    module procedure TTaggedwriter_int3
    module procedure TTaggedwriter_int4
    module procedure TTaggedwriter_logical0
    module procedure TTaggedwriter_logical1
    module procedure TTaggedwriter_logical2
    module procedure TTaggedwriter_logical3
    module procedure TTaggedwriter_logical4
  end interface


contains

  !> Initializes the tagged writer.
  subroutine TTaggedwriter_init(this)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(out) :: this

    !> Number of decimal, exponent, and character places
    integer :: nDec, nExp, nChar

    !> Number of resulting field entries
    integer :: nField

    ! example: "-3.1234567E-123 " would correspond to nDec = 7, nExp = 3, nChar = 16
    nexp = ceiling(log(maxexponent(1.0_dp) / log(10.0)) / log(10.0))
    ndec = precision(1.0_dp)

    nchar = ndec + nexp + 6
    nfield = 80 / nchar

    if (nfield == 0) then
      nfield = 1
    end if

    write(this%formReal, "('(', I2.2, 'ES', I2.2, '.', I2.2, 'E', I3.3, ')')") nField, nChar,&
        & nDec, nExp

    write(this%formCmplx, "('(', I2.2, '(2ES', I2.2, '.', I2.2, 'E', I3.3, '))')") nfield / 2,&
        & nchar, ndec, nexp

    !! "-12345 "
    nchar = digits(1) + 2
    nfield = 80 / nchar

    if (nfield == 0) then
      nfield = 1
    end if

    write(this%formInt, "('(', I2.2, 'I', I2.2, ')')") nfield, nchar

    write(this%formLogical, "('(40L2)')")

  end subroutine TTaggedwriter_init


  subroutine TTaggedwriter_real0(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    real(dp), intent(in) :: data

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formReal
    end if

    write(file, "('@', A, ':real:0:')") trim(tag)
    write(file, form) data

  end subroutine TTaggedwriter_real0


  subroutine TTaggedwriter_real1(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    real(dp), intent(in) :: data(:)

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formReal
    end if

    write(file, "('@', A, ':real:1:', I0)") trim(tag), shape(data)
    write(file, form) data

  end subroutine TTaggedwriter_real1


  subroutine TTaggedwriter_real2(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    real(dp), intent(in) :: data(:,:)

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formReal
    end if

    write(file, "('@', A, ':real:2:', I0, ',', I0)") trim(tag), shape(data)
    write(file, form) data
    
  end subroutine TTaggedwriter_real2


  subroutine TTaggedwriter_real3(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    real(dp), intent(in) :: data(:,:,:)

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formReal
    end if

    write(file, "('@', A, ':real:3:', I0, ',', I0, ',', I0)") trim(tag), shape(data)
    write(file, form) data
    
  end subroutine TTaggedwriter_real3


  subroutine TTaggedwriter_real4(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    real(dp), intent(in) :: data(:,:,:,:)

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formReal
    end if

    write(file, "('@', A, ':real:4:', I0, ',', I0, ',', I0, ',', I0)") trim(tag), shape(data)
    write(file, form) data

  end subroutine TTaggedwriter_real4


  subroutine TTaggedwriter_cplx0(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    complex(dp), intent(in) :: data

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formCmplx
    end if

    write(file, "('@', A, ':complex:0:')") trim(tag)
    write(file, form) data

  end subroutine TTaggedwriter_cplx0


  subroutine TTaggedwriter_cplx1(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    complex(dp), intent(in) :: data(:)

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formCmplx
    end if

    write(file, "('@', A, ':complex:1:', I0)") trim(tag), shape(data)
    write(file, form) data
    
  end subroutine TTaggedwriter_cplx1


  subroutine TTaggedwriter_cplx2(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    complex(dp), intent(in) :: data(:,:)

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formCmplx
    end if

    write(file, "('@', A, ':complex:2:', I0, ',', I0)") trim(tag), shape(data)
    write(file, form) data
    
  end subroutine TTaggedwriter_cplx2


  subroutine TTaggedwriter_cplx3(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    complex(dp), intent(in) :: data(:,:,:)

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formCmplx
    end if

    write(file, "('@', A, ':complex:3:', I0, ',', I0, ',', I0)") trim(tag), shape(data)
    write(file, form) data
    
  end subroutine TTaggedwriter_cplx3


  subroutine TTaggedwriter_cplx4(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    complex(dp), intent(in) :: data(:,:,:,:)

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formCmplx
    end if

    write(file, "('@', A, ':complex:4:', I0, ',', I0, ',', I0, ',', I0)") trim(tag), shape(data)
    write(file, form) data
    
  end subroutine TTaggedwriter_cplx4


  subroutine TTaggedwriter_int0(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    integer, intent(in) :: data

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formInt
    end if

    write(file, "('@', A, ':integer:0:')") trim(tag)
    write(file, form) data

  end subroutine TTaggedwriter_int0


  subroutine TTaggedwriter_int1(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    integer, intent(in) :: data(:)

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formInt
    end if

    write(file, "('@', A, ':integer:1:', I0)") trim(tag), shape(data)
    write(file, form) data
    
  end subroutine TTaggedwriter_int1


  subroutine TTaggedwriter_int2(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    integer, intent(in) :: data(:,:)

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formInt
    end if

    write(file, "('@', A, ':integer:2:', I0, ',', I0)") trim(tag), shape(data)
    write(file, form) data
    
  end subroutine TTaggedwriter_int2


  subroutine TTaggedwriter_int3(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    integer, intent(in) :: data(:,:,:)

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formInt
    end if

    write(file, "('@', A, ':integer:3:', I0, ',', I0, ',', I0)") trim(tag), shape(data)
    write(file, form) data
    
  end subroutine TTaggedwriter_int3


  subroutine TTaggedwriter_int4(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    integer, intent(in) :: data(:,:,:,:)

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formInt
    end if

    write(file, "('@', A, ':integer:4:', I0, ',', I0, ',', I0, ',', I0)") trim(tag), shape(data)
    write(file, form) data

  end subroutine TTaggedwriter_int4


  subroutine TTaggedwriter_logical0(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    logical, intent(in) :: data

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formLogical
    end if

    write(file, "('@', A, ':logical:0:')") trim(tag)
    write(file, form) data

  end subroutine TTaggedwriter_logical0


  subroutine TTaggedwriter_logical1(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    logical, intent(in) :: data(:)

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formLogical
    end if

    write(file, "('@', A, ':logical:1:', I0)") trim(tag), shape(data)
    write(file, form) data
    
  end subroutine TTaggedwriter_logical1


  subroutine TTaggedwriter_logical2(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    logical, intent(in) :: data(:,:)

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formLogical
    end if

    write(file, "('@', A, ':logical:2:', I0, ',', I0)") trim(tag), shape(data)
    write(file, form) data
    
  end subroutine TTaggedwriter_logical2


  subroutine TTaggedwriter_logical3(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    logical, intent(in) :: data(:,:,:)

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formLogical
    end if

    write(file, "('@', A, ':logical:3:', I0, ',', I0, ',', I0)") trim(tag), shape(data)
    write(file, form) data
    
  end subroutine TTaggedwriter_logical3


  subroutine TTaggedwriter_logical4(this, file, tag, data, optform)

    !> Instance of a tag format writer
    type(TTaggedwriter), intent(in) :: this

    !> File ID
    integer, intent(in) :: file

    !> Tag label
    character(len=*), intent(in) :: tag

    !> Data to print
    logical, intent(in) :: data(:,:,:,:)

    !> Optional formatting string
    character(lenFormStr), optional, intent(in) :: optform

    !> Actual formatting string
    character(lenFormStr) :: form

    if (present(optform)) then
      form = optform
    else
      form = this%formLogical
    end if

    write(file, "('@', A, ':logical:4:', I0, ',', I0, ',', I0, ',', I0)") trim(tag), shape(data)
    write(file, form) data
    
  end subroutine TTaggedwriter_logical4

end module common_taggedout
