!
!   test_mlfv.f90
!
!   A fortran f90 code to test mlfv_mod (Mittag-Leffler f90 implementation)
!   
!   
!   Created by Davide Verotta on 3/12/10.
!   Copyright 2010 UCSF. All rights reserved.
!   
!   Modified by Eduardo Mendes (with the permission of Davide Verotta) on 5/12/15 
!              to fix bugs and to deal with complex numbers.
!
! Obs: Please nota that the way the Mittag-Leffler function is called.
!      In case z is not complex, use real to get real part of the result.

program test_mlfv

use mlfv_mod

implicit none

integer, parameter :: n=10
real(8) :: alpha,beta,y
real(8), dimension(n) :: x
real(8), dimension(n) :: t
complex(8) :: z,aux,r
integer :: fi
integer ::  i
character(len=3) :: imag_unit = '+i*'

	
x(2)=-0.453750000000000   
x(1)=-1.15500000000000
x(3)=-1.61250000000000
x(4)=-1.72500000000000
x(5)=0
	
print*,"-------------------------------------------------------------------------------"
print*, " x                   alpha               beta                y"
print*,"-------------------------------------------------------------------------------"


do i=1,5
    fi=4
    alpha=.2
    beta=.2
    z=dcmplx(x(i),0)
    y=real(mlfv(alpha,beta,z,fi))
    write(*,"(4E20.12)"),x(i),alpha,beta,y
end do

print *," "

print *,'*******************************************'
print *,'Mittag-Leffer function for a complex number'
print *,'*******************************************'
print *, ' '


print*,"-----------------------------------------------------------------------------------------------------------------------"
print*, " x                                      alpha               beta                y"
print*,"-----------------------------------------------------------------------------------------------------------------------"

r=dcmplx(0d0,5.1962d0)
alpha=0.6d0
beta=1.0d0
fi=6

aux=mlfv(alpha,beta,r,fi)
if (aimag(aux)<0.) then 
   imag_unit = '-i*'
end if

write(*,"(1E18.12,a3,1E18.12,2x,1E18.12,2x,1E18.12,2x,1E18.12,a3,1E18.12)") real(r),'+i*',aimag(r),alpha,beta,&
        &real(aux),imag_unit,aimag(aux)

	
print *," "
	
!  Added by Eduardo Mendes on 5/12/15

alpha=0.5d0
beta=1.0d0
fi=3


print *,'***********************************'
print *,'Mittag-Leffer function for t**alpha'
print *,'***********************************'
print *, ' '

print*,"-------------------------------------------------------------------------------"
print*, " t                   alpha               beta                y                "
print*,"-------------------------------------------------------------------------------"

	
do i=1,n
     t(i)=(dble(i-1)*0.5d0)
     z=dcmplx(t(i)**alpha,0)
     y=real(mlfv(alpha,beta,z,fi))
     write(*,"(4E20.12)") t(i),alpha,beta,y
end do

alpha=1.2d0
beta=1.0d0
fi=5
	
print *,' '
print *,'************************************'
print *,'Mittag-Leffer function for -t**alpha'
print *,'************************************'
print *, ' '

print*,"-------------------------------------------------------------------------------"
print*, " t                   alpha               beta                y                "
print*,"-------------------------------------------------------------------------------"

	
do i=1,n
     t(i)=(dble(i-1)*0.5d0)
     z=dcmplx(-t(i)**alpha,0)
     y=real(mlfv(alpha,beta,z,fi))
     write(*,"(4E20.12)") t(i),alpha,beta,y
end do
	
end program test_mlfv