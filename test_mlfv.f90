! A fortran f90 code to test mlfv_mod (Mittag-Leffler f90 implementation)
!   
!   
! Created by Davide Verotta on 3/12/10.
! Copyright 2010 UCSF. All rights reserved.
!   
! Completely modified by Eduardo Mendes (with the permission of Davide Verotta) on 5/12/15 
!    to fix bugs, to deal with complex numbers and to include the derivative
!    of the Mittag-Leffler function
!
! Examples: I have tried to reproduce Figres 5 to 15 of the paper Computation of the Mittag-Leffler
!           function and its derivatives. Fract. Calc. Appl. Anal. 5(2002), 491-518 
!           by R.Gorenflo, J.Loutchko, and Yu. Luchko
!
! Obs.:  The derivative shown here is not the time derivative shown on the paper but the derivative
!        shown on Section 4 of the aforementioned paper.
!
!        Please nota that the way the Mittag-Leffler function is called.
!        In case where z is not complex, use real to get the real part of the result.

program test_mlfv

use mlfv_mod

implicit none

integer, parameter :: npoints = 4

real(8) :: alpha,beta,y,x
real(8), dimension(npoints+1) :: t
complex(8) :: z,aux
integer :: fi
integer ::  i
character(len=3) :: imag_unity = '+i*'
character(len=3) :: imag_unitz = '+i*'
real(8), parameter :: pi = 3.1415926535897932384626434D0

fi=8

! Figure 5

alpha=0.25d0
beta=1.0d0


print *,' '
print *,'********************************************************************************************************'
write(*,"(a36,f4.2,a9,f4.2,a20)") " Mittag-Leffer function for alpha = ",alpha,", beta = ", beta," E_{alpha,beta}(-t)"
print *,' '
print *,'Figure 5 of the paper Fract. Calc. Appl. Anal. 5(2002), 491-518 by R.Gorenflo, J.Loutchko and Yu. Luchko'
print *,'********************************************************************************************************'
print *, ' '

print*,"---------------------------------------------------------------------------------------------------"
print*, " t                   alpha               beta                y                   Derivative of E "
print*,"---------------------------------------------------------------------------------------------------"

	
do i=1,(npoints+1)
     t(i)=(dble(i-1)*(10d0-0d0)/dble(npoints)) ! 10
     z=dcmplx(-t(i),0)
     y=real(mlfv(alpha,beta,z,fi))
     x=real(mlfvderiv(alpha,beta,z,fi))
     write(*,"(5E20.12)") t(i),alpha,beta,y,x
end do

print *,' '

! Figure 6

alpha=1.75d0
beta=1.0d0


print *,' '
print *,'********************************************************************************************************'
write(*,"(a36,f4.2,a9,f4.2,a20)") " Mittag-Leffer function for alpha = ",alpha,", beta = ", beta," E_{alpha,beta}(-t)"
print *,' '
print *,'Figure 6 of the paper Fract. Calc. Appl. Anal. 5(2002), 491-518 by R.Gorenflo, J.Loutchko and Yu. Luchko'
print *,'********************************************************************************************************'
print *, ' '

print*,"---------------------------------------------------------------------------------------------------"
print*, " t                   alpha               beta                y                   Derivative of E "
print*,"---------------------------------------------------------------------------------------------------"

	
do i=1,(npoints+1)
     t(i)=(dble(i-1)*(50d0-0d0)/dble(npoints)) ! 50
     z=dcmplx(-t(i),0)
     y=real(mlfv(alpha,beta,z,fi))
     x=real(mlfvderiv(alpha,beta,z,fi))
     write(*,"(5E20.12)") t(i),alpha,beta,y,x
end do

print *,' '

! Figure 7

alpha=2.25d0
beta=1.0d0


print *,' '
print *,'********************************************************************************************************'
write(*,"(a36,f4.2,a9,f4.2,a20)") " Mittag-Leffer function for alpha = ",alpha,", beta = ", beta," E_{alpha,beta}(-t)"
print *,' '
print *,'Figure 7 of the paper Fract. Calc. Appl. Anal. 5(2002), 491-518 by R.Gorenflo, J.Loutchko and Yu. Luchko'
print *,'********************************************************************************************************'
print *, ' '

print*,"---------------------------------------------------------------------------------------------------"
print*, " t                   alpha               beta                y                   Derivative of E "
print*,"---------------------------------------------------------------------------------------------------"

	
do i=1,(npoints+1)
     t(i)=(dble(i-1)*(100d0-0d0)/dble(npoints)) ! 100
     z=dcmplx(-t(i),0)
     y=real(mlfv(alpha,beta,z,fi))
     x=real(mlfvderiv(alpha,beta,z,fi))
     write(*,"(5E20.12)") t(i),alpha,beta,y,x
end do

print *,' '

!  Figure 8

alpha=0.75d0
beta=1.0d0

print *,'********************************************************************************************************'
write(*,"(a36,f4.2,a9,f4.2,a20)") " Mittag-Leffer function for alpha = ",alpha,", beta = ", beta,"arg(z) = alpha*pi/4"
print *,' '
print *,'Figure 8 of the paper Fract. Calc. Appl. Anal. 5(2002), 491-518 by R.Gorenflo, J.Loutchko and Yu. Luchko'
print *,'********************************************************************************************************'
print *, ' '

print*,"-----------------------------------------------------------------------------------------------------------------------"
print*, " t                  z               alpha   beta                y                  abs(y) "
print*,"-----------------------------------------------------------------------------------------------------------------------"

do i=1,(npoints+1)
     t(i)=(dble(i-1)*(5d0-0d0)/dble(npoints))
     z=dcmplx(t(i),0)*exp(dcmplx(0,1)*alpha*pi/4d0)
     if (aimag(z)<0.) then 
       imag_unitz = '-i*'
     else
       imag_unitz = '+i*'
     end if
     aux=mlfv(alpha,beta,z,fi)
     if (aimag(aux)<0.) then 
       imag_unity = '-i*'
     else
       imag_unity = '+i*'
     end if
!     write(*,"(4E20.12)") t(i),alpha,beta,y
     write(*,"(f6.2,3x,e11.5,a3,e11.5,3x,f4.2,4x,f4.2,3x,e11.5,a3,e11.5,4x,e11.5)") t(i),real(z),&
                         &imag_unitz,abs(aimag(z)),alpha,beta,real(aux),imag_unity,abs(aimag(aux)),abs(aux)
end do

print *,' '

!  Figure 9

alpha=0.75d0
beta=1.0d0

print *,'********************************************************************************************************'
write(*,"(a36,f4.2,a9,f4.2,a20)") " Mittag-Leffer function for alpha = ",alpha,", beta = ", beta,"arg(z) = alpha*pi/2"
print *,' '
print *,'Figure 9 of the paper Fract. Calc. Appl. Anal. 5(2002), 491-518 by R.Gorenflo, J.Loutchko and Yu. Luchko'
print *,'********************************************************************************************************'
print *, ' '

print*,"-----------------------------------------------------------------------------------------------------------------------"
print*, " t                  z               alpha   beta                y                  abs(y) "
print*,"-----------------------------------------------------------------------------------------------------------------------"

do i=1,(npoints+1)
     t(i)=(dble(i-1)*(50d0-0d0)/dble(npoints))
     z=dcmplx(t(i),0)*exp(dcmplx(0,1)*alpha*pi/2d0)
     if (aimag(z)<0.) then 
       imag_unitz = '-i*'
     else
       imag_unitz = '+i*'
     end if
     aux=mlfv(alpha,beta,z,fi)
     if (aimag(aux)<0.) then 
       imag_unity = '-i*'
     else
       imag_unity = '+i*'
     end if
!     write(*,"(4E20.12)") t(i),alpha,beta,y
     write(*,"(f6.2,3x,e11.5,a3,e11.5,3x,f4.2,4x,f4.2,3x,e11.5,a3,e11.5,4x,e11.5)") t(i),real(z),&
                         &imag_unitz,abs(aimag(z)),alpha,beta,real(aux),imag_unity,abs(aimag(aux)),abs(aux)
end do

print *,' '

!  Figure 10

alpha=0.75d0
beta=1.0d0

print *,'**********************************************************************************************************'
write(*,"(a36,f4.2,a9,f4.2,a22)") " Mittag-Leffer function for alpha = ",alpha,", beta = ", beta,"arg(z) = alpha*3*pi/4"
print *,' '
print *,'Figure 10 of the paper Fract. Calc. Appl. Anal. 5(2002), 491-518 by R.Gorenflo, J.Loutchko and Yu. Luchko'
print *,'**********************************************************************************************************'
print *, ' '

print*,"-----------------------------------------------------------------------------------------------------------------------"
print*, " t                  z               alpha   beta                y                  abs(y) "
print*,"-----------------------------------------------------------------------------------------------------------------------"

do i=1,(npoints+1)
     t(i)=(dble(i-1)*(50d0-0d0)/dble(npoints))
     z=dcmplx(t(i),0)*exp(dcmplx(0,1)*alpha*3d0*pi/4d0)
     if (aimag(z)<0.) then 
       imag_unitz = '-i*'
     else
       imag_unitz = '+i*'
     end if
     aux=mlfv(alpha,beta,z,fi)
     if (aimag(aux)<0.) then 
       imag_unity = '-i*'
     else
       imag_unity = '+i*'
     end if
!     write(*,"(4E20.12)") t(i),alpha,beta,y
     write(*,"(f6.2,3x,e11.5,a3,e11.5,3x,f4.2,4x,f4.2,3x,e11.5,a3,e11.5,4x,e11.5)") t(i),real(z),&
                         &imag_unitz,abs(aimag(z)),alpha,beta,real(aux),imag_unity,abs(aimag(aux)),abs(aux)
end do

print *,' '

!  Figure 11

alpha=0.75d0
beta=1.0d0

print *,'**********************************************************************************************************'
write(*,"(a36,f4.2,a9,f4.2,a22)") " Mittag-Leffer function for alpha = ",alpha,", beta = ", beta,"arg(z) = pi"
print *,' '
print *,'Figure 11 of the paper Fract. Calc. Appl. Anal. 5(2002), 491-518 by R.Gorenflo, J.Loutchko and Yu. Luchko'
print *,'**********************************************************************************************************'
print *, ' '

print*,"-----------------------------------------------------------------------------------------------------------------------"
print*, " t                  z               alpha   beta                y                  real(y) "
print*,"-----------------------------------------------------------------------------------------------------------------------"

do i=1,(npoints+1)
     t(i)=(dble(i-1)*(100d0-0d0)/dble(npoints))
     z=dcmplx(t(i),0)*exp(dcmplx(0,1)*pi)
     if (aimag(z)<0.) then 
       imag_unitz = '-i*'
     else
       imag_unitz = '+i*'
     end if
     aux=mlfv(alpha,beta,z,fi)
     if (aimag(aux)<0.) then 
       imag_unity = '-i*'
     else
       imag_unity = '+i*'
     end if
!     write(*,"(4E20.12)") t(i),alpha,beta,y
     write(*,"(f6.2,3x,e11.5,a3,e11.5,3x,f4.2,4x,f4.2,3x,e11.5,a3,e11.5,4x,e11.5)") t(i),real(z),&
                         &imag_unitz,abs(aimag(z)),alpha,beta,real(aux),imag_unity,abs(aimag(aux)),real(aux)
end do

print *,' '

!  Figure 12

alpha=1.25d0
beta=1.0d0

print *,'********************************************************************************************************'
write(*,"(a36,f4.2,a9,f4.2,a20)") " Mittag-Leffer function for alpha = ",alpha,", beta = ", beta,"arg(z) = alpha*pi/4"
print *,' '
print *,'Figure 12 of the paper Fract. Calc. Appl. Anal. 5(2002), 491-518 by R.Gorenflo, J.Loutchko and Yu. Luchko'
print *,'********************************************************************************************************'
print *, ' '

print*,"-----------------------------------------------------------------------------------------------------------------------"
print*, " t                  z               alpha   beta                y                  abs(y) "
print*,"-----------------------------------------------------------------------------------------------------------------------"

do i=1,(npoints+1)
     t(i)=(dble(i-1)*(5d0-0d0)/dble(npoints))
     z=dcmplx(t(i),0)*exp(dcmplx(0,1)*alpha*pi/4d0)
     if (aimag(z)<0.) then 
       imag_unitz = '-i*'
     else
       imag_unitz = '+i*'
     end if
     aux=mlfv(alpha,beta,z,fi)
     if (aimag(aux)<0.) then 
       imag_unity = '-i*'
     else
       imag_unity = '+i*'
     end if
!     write(*,"(4E20.12)") t(i),alpha,beta,y
     write(*,"(f6.2,3x,e11.5,a3,e11.5,3x,f4.2,4x,f4.2,3x,e11.5,a3,e11.5,4x,e11.5)") t(i),real(z),&
                         &imag_unitz,abs(aimag(z)),alpha,beta,real(aux),imag_unity,abs(aimag(aux)),abs(aux)
end do

print *,' '

!  Figure 13

alpha=1.25d0
beta=1.0d0

print *,'********************************************************************************************************'
write(*,"(a36,f4.2,a9,f4.2,a20)") " Mittag-Leffer function for alpha = ",alpha,", beta = ", beta,"arg(z) = alpha*pi/2"
print *,' '
print *,'Figure 13 of the paper Fract. Calc. Appl. Anal. 5(2002), 491-518 by R.Gorenflo, J.Loutchko and Yu. Luchko'
print *,'********************************************************************************************************'
print *, ' '

print*,"-----------------------------------------------------------------------------------------------------------------------"
print*, " t                  z               alpha   beta                y                  abs(y) "
print*,"-----------------------------------------------------------------------------------------------------------------------"

do i=1,(npoints+1)
     t(i)=(dble(i-1)*(50d0-0d0)/dble(npoints))
     z=dcmplx(t(i),0)*exp(dcmplx(0,1)*alpha*pi/2d0)
     if (aimag(z)<0.) then 
       imag_unitz = '-i*'
     else
       imag_unitz = '+i*'
     end if
     aux=mlfv(alpha,beta,z,fi)
     if (aimag(aux)<0.) then 
       imag_unity = '-i*'
     else
       imag_unity = '+i*'
     end if
!     write(*,"(4E20.12)") t(i),alpha,beta,y
     write(*,"(f6.2,3x,e11.5,a3,e11.5,3x,f4.2,4x,f4.2,3x,e11.5,a3,e11.5,4x,e11.5)") t(i),real(z),&
                         &imag_unitz,abs(aimag(z)),alpha,beta,real(aux),imag_unity,abs(aimag(aux)),abs(aux)
end do

print *,' '

!  Figure 14

alpha=1.25d0
beta=1.0d0

print *,'**********************************************************************************************************'
write(*,"(a36,f4.2,a9,f4.2,a22)") " Mittag-Leffer function for alpha = ",alpha,", beta = ", beta,"arg(z) = alpha*3*pi/4"
print *,' '
print *,'Figure 14 of the paper Fract. Calc. Appl. Anal. 5(2002), 491-518 by R.Gorenflo, J.Loutchko and Yu. Luchko'
print *,'**********************************************************************************************************'
print *, ' '

print*,"-----------------------------------------------------------------------------------------------------------------------"
print*, " t                  z               alpha   beta                y                  abs(y) "
print*,"-----------------------------------------------------------------------------------------------------------------------"

do i=1,(npoints+1)
     t(i)=(dble(i-1)*(50d0-0d0)/dble(npoints))
     z=dcmplx(t(i),0)*exp(dcmplx(0,1)*alpha*3d0*pi/4d0)
     if (aimag(z)<0.) then 
       imag_unitz = '-i*'
     else
       imag_unitz = '+i*'
     end if
     aux=mlfv(alpha,beta,z,fi)
     if (aimag(aux)<0.) then 
       imag_unity = '-i*'
     else
       imag_unity = '+i*'
     end if
!     write(*,"(4E20.12)") t(i),alpha,beta,y
     write(*,"(f6.2,3x,e11.5,a3,e11.5,3x,f4.2,4x,f4.2,3x,e11.5,a3,e11.5,4x,e11.5)") t(i),real(z),&
                         &imag_unitz,abs(aimag(z)),alpha,beta,real(aux),imag_unity,abs(aimag(aux)),abs(aux)
end do

print *,' '

!  Figure 15

alpha=1.25d0
beta=1.0d0

print *,'**********************************************************************************************************'
write(*,"(a36,f4.2,a9,f4.2,a22)") " Mittag-Leffer function for alpha = ",alpha,", beta = ", beta,"arg(z) = pi"
print *,' '
print *,'Figure 15 of the paper Fract. Calc. Appl. Anal. 5(2002), 491-518 by R.Gorenflo, J.Loutchko and Yu. Luchko'
print *,'**********************************************************************************************************'
print *, ' '

print*,"-----------------------------------------------------------------------------------------------------------------------"
print*, " t                  z               alpha   beta                y                  real(y) "
print*,"-----------------------------------------------------------------------------------------------------------------------"

do i=1,(npoints+1)
     t(i)=(dble(i-1)*(100d0)/dble(npoints))
     z=dcmplx(t(i),0)*exp(dcmplx(0,1)*pi)
     if (aimag(z)<0.) then 
       imag_unitz = '-i*'
     else
       imag_unitz = '+i*'
     end if
     aux=mlfv(alpha,beta,z,fi)
     if (aimag(aux)<0.) then 
       imag_unity = '-i*'
     else
       imag_unity = '+i*'
     end if
!     write(*,"(4E20.12)") t(i),alpha,beta,y
     write(*,"(f6.2,3x,e11.5,a3,e11.5,3x,f4.2,4x,f4.2,3x,e11.5,a3,e11.5,4x,e11.5)") t(i),real(z),&
                         &imag_unitz,abs(aimag(z)),alpha,beta,real(aux),imag_unity,abs(aimag(aux)),real(aux)
end do

print *,' '
	
end program test_mlfv