!!* Read input from stdin
module input

  use common_accuracy, only : dp

  implicit none
  private 
  
  public :: read_input_1, read_input_2, echo_input
  
contains

  subroutine read_input_1(nuc,max_l,occ_shells,maxiter,poly_order,&
      &min_alpha,max_alpha,num_alpha,generate_alpha,alpha,&
      &conf_r0,conf_power,num_occ,num_power,num_alphas,xcnr,&
      &eigprint,zora,broyden,mixing_factor,xalpha_const)

    ! Read in everything except occupation numbers.

    integer :: ii,jj
    integer, intent(out)  :: nuc,max_l,maxiter,num_occ,num_power
    real(dp), intent(out) :: conf_power(0:)
    integer, intent(out) :: num_alphas,xcnr
    logical, intent(out) :: generate_alpha,eigprint,zora,broyden
    real(dp), intent(out) :: conf_r0(0:),min_alpha,max_alpha,mixing_factor
    real(dp), intent(out) :: alpha(0:,:),xalpha_const
    integer, intent(out)  :: occ_shells(0:),num_alpha(0:),poly_order(0:)


    write(*,'(A)') 'Enter nuclear charge, maximal angular momentum (s=0), &
        &max. SCF, ZORA'
    read(*,*) nuc,max_l,maxiter,zora

    write(*,'(A)') 'Enter XC functional, 0=HF, 1=X-Alpha, 2=PW-LDA, 3=PBE'
    read(*,*) xcnr
    if (xcnr==0) write(*,'(A)') 'WARNING: ONLY CORRECT FOR CLOSED SHELL 1S !'
    if ((xcnr==0).and.zora) then
      write(*,'(A)') 'ZORA only available for DFT !'
      STOP
    end if
    if (xcnr==1) then
      write(*,'(A)') 'Enter empirical parameter for X-Alpha exchange'
      read(*,*) xalpha_const
    end if

    if (max_l>4) then
      write(*,'(A)') 'Sorry, l=',max_l,' is a bit too large. No nuclear weapons&
          &allowed.'
      STOP
    end if

    write(*,'(A)') 'Enter Confinement: r_0 and integer power, power=0 -> off'
    do ii=0,max_l
      write(*,'(A,I3)') 'l=',ii
      read(*,*) conf_r0(ii),conf_power(ii)
    end do

    write(*,'(A)') 'Enter number of occupied shells for each angular momentum'
    do ii=0,max_l
      write(*,'(A,I3)') 'l=',ii
      read(*,*) occ_shells(ii)
    end do

    write(*,'(A)') 'Enter number of exponents and polynomial coefficients for each angular momentum'
    do ii=0,max_l
      write(*,'(A,I3)') 'l=',ii
      read(*,*) num_alpha(ii),poly_order(ii)
      if (num_alpha(ii)>10) then
        write(*,'(A)') ' Sorry, num_alpha must be smaller than 11.'
        STOP
      end if
    end do

    !  write(*,'(A)') 'Enter number of exponents for each angular momentum'
    !  do ii=0,max_l
    !    write(*,'(A,I3)') 'l=',ii
    !    read(*,*) num_alpha(ii)
    !    if (num_alpha(ii)>10) then
    !      write(*,'(A)') ' Sorry, num_alpha must be smaller than 11.'
    !      STOP
    !    end if
    !  end do

    write(*,'(A)') 'Do you want to generate the exponents ? .true./.false.'
    read(*,*) generate_alpha

    if (generate_alpha) then
      ! generate alphas
      !
      do ii=0,max_l
        write(*,'(A)') 'Enter smallest exponent and largest exponent.'
        read(*,*) min_alpha,max_alpha
        !
        call gen_alphas(min_alpha,max_alpha,num_alpha(ii),alpha(ii,:))
      end do
    else
      do ii=0,max_l
        write(*,'(A,I3,A,I3,A)') 'Enter ',num_alpha(ii),'exponents for l=',&
            &ii,' one per line'
        do jj=1,num_alpha(ii)
          read(*,*) alpha(ii,jj)
        end do
      end do
    end if

    num_occ=0
    do ii=0,max_l
      num_occ=max(num_occ,occ_shells(ii))
    end do

    num_power=0
    do ii=0,max_l
      num_power=max(num_power,poly_order(ii))
    end do

    num_alphas=0
    do ii=0,max_l
      num_alphas=max(num_alphas,num_alpha(ii))
    end do

    write(*,'(A)') 'Print Eigenvectors ? .true./.false.'
    read(*,*) eigprint

    write(*,'(A)') ' Use Broyden mixer ? .true./.false. and mixing parameter <1'
    read(*,*) broyden,mixing_factor

  end subroutine read_input_1

  subroutine read_input_2(occ,max_l,occ_shells, qnvalorbs)

    ! Read in occupation numbers.

    real(dp), intent(out) :: occ(:,0:,:)
    integer, intent(in) :: max_l,occ_shells(0:)
    integer, intent(out) :: qnvalorbs(:,0:)
    integer :: ii,jj

    write(*,'(A)') 'Enter the occupation numbers for each angular momentum&
        & and shell, up and down in one row'

    occ=0.0d0             

    write(*,'(A)') ' '
    write(*,'(A)') 'UP Electrons DOWN Electrons'            
    do ii=0,max_l
      do jj=1,occ_shells(ii)
        write(*,'(A,I3,A,I3)') 'l= ',ii,' and shell ',jj
        read(*,*) occ(1,ii,jj),occ(2,ii,jj)
      end do
    end do

    write(*,"(A)") "Quantum numbers of wavefunctions to be written:"
    do ii = 0, max_l
      write(*, "(A,I0,A)") "l= ", ii, ": from to"
      read(*,*) qnvalorbs(:, ii)
      qnvalorbs(:,ii) = (/ minval(qnvalorbs(:,ii)), maxval(qnvalorbs(:,ii)) /)
      qnvalorbs(:,ii) = qnvalorbs(:,ii) - ii
    end do

  end subroutine read_input_2

  subroutine echo_input(nuc,max_l,occ_shells,maxiter,poly_order,num_alpha,&
      &alpha,conf_r0,conf_power,occ,num_occ,num_power,&
      &num_alphas,xcnr,zora,num_mesh_points,xalpha_const)

    ! Echo completed input to stdout.

    integer :: ii,jj
    integer, intent(in)  :: nuc,max_l,maxiter,num_occ,num_power
    real(dp), intent(in) :: conf_power(0:)
    integer, intent(in)  :: num_alphas,xcnr,num_mesh_points
    real(dp), intent(in) :: conf_r0(0:),occ(:,0:,:)
    real(dp), intent(in) :: alpha(0:,:),xalpha_const
    integer, intent(in)  :: occ_shells(0:),num_alpha(0:),poly_order(0:)
    logical, intent(in) :: zora

    write(*,'(A)') ' '
    write(*,'(A)') '--------------'
    write(*,'(A)') 'INPUT SUMMARY '
    write(*,'(A)') '--------------'

    if (zora) write(*,'(A)') 'SCALAR RELATIVISTIC ZORA CALCULATION'
    if (.not.zora) write(*,'(A)') 'NON-RELATIVISTIC CALCULATION'
    write(*,'(A)') ' '
    write(*,'(A,I3)') 'Nuclear Charge: ',nuc
    if (xcnr==0) write(*,'(A,I3)') 'HF Exchange, only correct for closed shell !'
    if (xcnr==1) write(*,'(A,F12.8)') 'X-Alpha, alpha= ',xalpha_const
    if (xcnr==2) write(*,'(A,I3)') 'LDA, Perdew-Wang Parametrization'
    if (xcnr==3) write(*,'(A,I3)') 'PBE'
    write(*,'(A,I1)') 'Max. angular momentum: ',max_l
    write(*,'(A,I5)') 'Number of points for numerical radial integration: ',&
        &num_mesh_points

    write(*,'(A)') ' '
    do ii=0,max_l
      write(*,'(A,I1,A,I2)') 'Occupied Shells for l=',ii,': ',occ_shells(ii)
    end do

    write(*,'(A)') ' '
    do ii=0,max_l
      write(*,'(A,I1,A,I2)') 'Number of Polynomial Coeff. for l=',ii,': ',poly_order(ii)
    end do

    write(*,'(A)') ' '
    do ii=0,max_l
      write(*,'(A,I1)') 'Exponents for l=',ii
      do jj=1,num_alpha(ii)
        write(*,'(F12.8)') alpha(ii,jj)
      end do
    end do

    write(*,'(A)') ' '
    write(*,'(A)') 'Occupation Numbers UP/DWN'
    do ii=0,max_l
      do jj=1,occ_shells(ii)
        write(*,'(A,I1,A,I2,A,2F12.8)') 'Angular Momentum ',ii,' Shell ',jj,&
            &': ',occ(1,ii,jj),occ(2,ii,jj)
      end do
    end do
    !
    !  write(*,'(A)') ' '
    !  write(*,'(A)') 'Occupation Numbers DWN'
    !  do ii=0,max_l
    !    do jj=1,occ_shells(ii)
    !      write(*,'(A,I1,A,I2,A,F12.8)') 'Angular Momentum ',ii,' Shell ',jj,&
    !                                     &': ',occ(2,ii,jj)
    !    end do
    !  end do

    write(*,'(A)') ' '
    !  write(*,'(A,F12.8,A,I1)') 'Confining Radius is ',conf_r0,' a.u. with power of ',conf_power
    do ii=0,max_l
      if (conf_power(ii) > 1d-6) then
        write(*,'(A,I3,A,E15.7,A,E15.7)') 'l= ',ii,', r0= ',conf_r0(ii),' power= ',&
            conf_power(ii)
      else
        write(*,'(A,I3,A)') 'l= ',ii,' no confinement '
      end if
    end do

    write(*,'(A)') ' '
    write(*,'(A,I2,A)') 'There are at maximum ',num_occ,' occ. shells for one l'
    write(*,'(A,I2,A)') 'There are at maximum ',num_power,' coefficients for one&
        & exponent'
    write(*,'(A,I2,A)') 'There are at maximum ',num_alphas,' exponents'

    write(*,'(A)') ' '
    write(*,'(A)') '------------------'
    write(*,'(A)') 'END INPUT SUMMARY '
    write(*,'(A)') '------------------'
    write(*,'(A)') ' '

  end subroutine echo_input

  subroutine gen_alphas(min_alpha,max_alpha,num_alpha,alpha)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Generate alpha coefficients for Slater expansion
    ! 
    ! min_alpha               : smallest alpha
    ! max_alpha               : largest alpha
    ! num_alpha               : number of alphas
    ! alpha                   : output, generated alphas
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none

    integer :: ii
    real(dp), intent(in) :: min_alpha,max_alpha
    real(dp) :: alpha(10)
    integer :: num_alpha
    real(dp) :: beta(10),f
    do ii=1,10
      alpha(ii)=0.0_dp
    end do
    alpha(1)=min_alpha
    if (num_alpha==1) return
    f=(max_alpha/alpha(1))**(1.0d0/FLOAT((num_alpha-1)))
    do ii=1,(num_alpha-1)
      alpha(1+ii)=alpha(ii)*f
    end do
    do ii=1,num_alpha
      beta(num_alpha+1-ii)=alpha(ii)
    end do
    do ii=1,num_alpha
      alpha(ii)=beta(ii)
    end do
    return
  end subroutine gen_alphas

end module input
