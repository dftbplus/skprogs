!> Contains routines to write out various data structures in a comprehensive
!! tagged format.
module taggedout
  use accuracy, only : dp
  implicit none
  private

  public :: taggedwriter, init, writetag, taglen

  integer, parameter :: taglen = 40
  integer, parameter :: formlen = 20

  type :: taggedwriter
    character(formlen) :: form_real
    character(formlen) :: form_cmplx
    character(formlen) :: form_int
    character(formlen) :: form_logical
  end type taggedwriter

  interface init
    module procedure taggedwriter_init
  end interface

  !> Writes objects in a standardized tagged form to a given file.
  interface writetag
    module procedure taggedwriter_real0
    module procedure taggedwriter_real1
    module procedure taggedwriter_real2
    module procedure taggedwriter_real3
    module procedure taggedwriter_real4
    module procedure taggedwriter_cplx0
    module procedure taggedwriter_cplx1
    module procedure taggedwriter_cplx2
    module procedure taggedwriter_cplx3
    module procedure taggedwriter_cplx4
    module procedure taggedwriter_int0
    module procedure taggedwriter_int1
    module procedure taggedwriter_int2
    module procedure taggedwriter_int3
    module procedure taggedwriter_int4
    module procedure taggedwriter_logical0
    module procedure taggedwriter_logical1
    module procedure taggedwriter_logical2
    module procedure taggedwriter_logical3
    module procedure taggedwriter_logical4
  end interface


contains

  !> Initializes the tagged writer.
  subroutine taggedwriter_init(self)
    type(taggedwriter), intent(out) :: self

    integer :: ndec, nexp, nchar, nfield

    !! "-3.1234567E-123 ": nDec = 7, nexp = 3, nchar = 16
    nexp = ceiling(log(maxexponent(1.0_dp)/log(10.0))/log(10.0))
    ndec = precision(1.0_dp)
    nchar = ndec + nexp + 6
    nfield = 80 / nchar
    if (nfield == 0) then
      nfield = 1
    end if

99000 format('(', I2.2, 'ES', I2.2, '.', I2.2, 'E', I3.3, ')')
    write(self%form_real, 99000) nfield, nchar, ndec, nexp

99010 format('(', I2.2, '(2ES', I2.2, '.', I2.2, 'E', I3.3, '))')
    write(self%form_cmplx, 99010) nfield/2, nchar, ndec, nexp

    !! "-12345 "
    nchar = digits(1) + 2
    nfield = 80 / nchar
    if (nfield == 0) then
      nfield = 1
    end if

99020 format('(', I2.2, 'I', I2.2, ')')
    write (self%form_int, 99020) nfield, nchar

99030 format('(40L2)')
    write(self%form_logical, 99030)

  end subroutine taggedwriter_init



  subroutine taggedwriter_real0(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    real(dp), intent(in) :: val
    character(formlen), optional, intent(in) :: optform

    character(formlen) :: form    

    if (present(optform)) then
      form = optform
    else
      form = self%form_real
    end if

99040 format('@', A, ':real:0:')
    write(file, 99040) trim(tag)
    write(file, form) val

  end subroutine taggedwriter_real0



  subroutine taggedwriter_real1(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(len=*), intent(in) :: tag
    real(dp), intent(in) :: val(:)
    character(formlen), optional, intent(in) :: optform

    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_real
    end if

99050 format('@', A, ':real:1:', I0)
    write(file, 99050) trim(tag), shape(val)
    write(file, form) val

  end subroutine taggedwriter_real1



  subroutine taggedwriter_real2(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    real(dp), intent(in) :: val(:,:)
    character(formlen), optional, intent(in) :: optform

    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_real
    end if

99060 format('@', A, ':real:2:', I0, ',', I0)
    write(file, 99060) trim(tag), shape(val)
    write(file, form) val
    
  end subroutine taggedwriter_real2



  subroutine taggedwriter_real3(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    real(dp), intent(in) :: val(:,:,:)
    character(formlen), optional, intent(in) :: optform

    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_real
    end if

99070 format('@', A, ':real:3:', I0, ',', I0, ',', I0)
    write(file, 99070) trim(tag), shape(val)
    write(file, form) val
    
  end subroutine taggedwriter_real3



  subroutine taggedwriter_real4(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    real(dp), intent(in) :: val(:,:,:,:)
    character(formlen), optional, intent(in) :: optform

    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_real
    end if

99080 format('@', A, ':real:4:', I0, ',', I0, ',', I0, ',', I0)
    write(file, 99080) trim(tag), shape(val)
    write(file, form) val

  end subroutine taggedwriter_real4



  subroutine taggedwriter_cplx0(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    complex(dp), intent(in) :: val
    character(formlen), optional, intent(in) :: optform 

    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_cmplx
    end if

99090 format('@', A, ':complex:0:')
    write(file, 99090) trim(tag)
    write(file, form) val

  end subroutine taggedwriter_cplx0



  subroutine taggedwriter_cplx1(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    complex(dp), intent(in) :: val(:)
    character(formlen), optional, intent(in) :: optform

    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_cmplx
    end if

99100 format('@', A, ':complex:1:', I0)
    write(file, 99100) trim(tag), shape(val)
    write(file, form) val
    
  end subroutine taggedwriter_cplx1



  subroutine taggedwriter_cplx2(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    complex(dp), intent(in) :: val(:,:)
    character(formlen), optional, intent(in) :: optform
    
    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_cmplx
    end if

99110 format('@', A, ':complex:2:', I0, ',', I0)
    write(file, 99110) trim(tag), shape(val)
    write(file, form) val
    
  end subroutine taggedwriter_cplx2



  subroutine taggedwriter_cplx3(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    complex(dp), intent(in) :: val(:,:,:)
    character(formlen), optional, intent(in) :: optform
    
    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_cmplx
    end if

99120 format('@', A, ':complex:3:', I0, ',', I0, ',', I0)
    write(file, 99120) trim(tag), shape(val)
    write(file, form) val
    
  end subroutine taggedwriter_cplx3



  subroutine taggedwriter_cplx4(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    complex(dp), intent(in) :: val(:,:,:,:)
    character(formlen), optional, intent(in) :: optform

    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_cmplx
    end if

99130 format('@', A, ':complex:4:', I0, ',', I0, ',', I0, ',', I0)
    write(file, 99130) trim(tag), shape(val)
    write(file, form) val
    
  end subroutine taggedwriter_cplx4



  subroutine taggedwriter_int0(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    integer, intent(in) :: val
    character(formlen), optional, intent(in) :: optform

    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_int
    end if

99140 format('@', A, ':integer:0:')
    write(file, 99140) trim(tag)
    write(file, form) val

  end subroutine taggedwriter_int0



  subroutine taggedwriter_int1(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    integer, intent(in) :: val(:)
    character(formlen), optional, intent(in) :: optform
    
    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_int
    end if

99150 format('@', A, ':integer:1:', I0)
    write(file, 99150) trim(tag), shape(val)
    write(file, form) val
    
  end subroutine taggedwriter_int1



  subroutine taggedwriter_int2(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    integer, intent(in) :: val(:,:)
    character(formlen), optional, intent(in) :: optform

    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_int
    end if

99160 format('@', A, ':integer:2:', I0, ',', I0)
    write(file, 99160) trim(tag), shape(val)
    write(file, form) val
    
  end subroutine taggedwriter_int2



  subroutine taggedwriter_int3(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    integer, intent(in) :: val(:,:,:)
    character(formlen), optional, intent(in) :: optform
    
    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_int
    end if

99170 format('@', A, ':integer:3:', I0, ',', I0, ',', I0)
    write(file, 99170) trim(tag), shape(val)
    write(file, form) val
    
  end subroutine taggedwriter_int3



  subroutine taggedwriter_int4(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    integer, intent(in) :: val(:,:,:,:)
    character(formlen), optional, intent(in) :: optform
    
    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_int
    end if

99180 format('@', A, ':integer:4:', I0, ',', I0, ',', I0, ',', I0)
    write(file, 99180) trim(tag), shape(val)
    write(file, form) val

  end subroutine taggedwriter_int4



  subroutine taggedwriter_logical0(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    logical, intent(in) :: val
    character(formlen), optional, intent(in) :: optform

    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_logical
    end if

99190 format('@', A, ':logical:0:')
    write(file, 99190) trim(tag)
    write(file, form) val

  end subroutine taggedwriter_logical0



  subroutine taggedwriter_logical1(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    logical, intent(in) :: val(:)
    character(formlen), optional, intent(in) :: optform
    
    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_logical
    end if

99200 format('@', A, ':logical:1:', I0)
    write(file, 99200) trim(tag), shape(val)
    write(file, form) val
    
  end subroutine taggedwriter_logical1



  subroutine taggedwriter_logical2(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    logical, intent(in) :: val(:,:)
    character(formlen), optional, intent(in) :: optform
    
    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_logical
    end if

99210 format('@', A, ':logical:2:', I0, ',', I0)
    write(file, 99210) trim(tag), shape(val)
    write(file, form) val
    
  end subroutine taggedwriter_logical2



  subroutine taggedwriter_logical3(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    logical, intent(in) :: val(:,:,:)
    character(formlen), optional, intent(in) :: optform
    
    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_logical
    end if

99220 format('@', A, ':logical:3:', I0, ',', I0, ',', I0)
    write(file, 99220) trim(tag), shape(val)
    write(file, form) val
    
  end subroutine taggedwriter_logical3



  subroutine taggedwriter_logical4(self, file, tag, val, optform)
    type(taggedwriter), intent(in) :: self
    integer, intent(in) :: file
    character(*), intent(in) :: tag
    logical, intent(in) :: val(:,:,:,:)
    character(formlen), optional, intent(in) :: optform
    
    character(formlen) :: form

    if (present(optform)) then
      form = optform
    else
      form = self%form_logical
    end if

99230 format('@', A, ':logical:4:', I0, ',', I0, ',', I0, ',', I0)
    write(file, 99230) trim(tag), shape(val)
    write(file, form) val
    
  end subroutine taggedwriter_logical4


end module taggedout
