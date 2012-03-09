! File Jameson.f90
      subroutine jameson (xin,yin,sxx,sxy,syx,syy,uin,b)
      real(kind=8), dimension(64,128), intent (in) :: xin,yin,sxx&
&,sxy,syx,syy
      real(kind=8), dimension(4,64,128), intent (in) :: uin
      real(kind=8), dimension(128,64), intent (out) :: b
      integer :: n,m
      real(kind=8) :: testVal
      real(kind=8), dimension(64,128) :: x,y
      real(kind=8), dimension(4,64,128) :: u, unew
      real(kind=8), dimension(4) :: ff
      real(kind=8), dimension(2) :: f, testSol
!      external :: f, g


      print*, "***************************"
      print*, "** Begining Jameson.f90! **"
      print*, "***************************"

      print*, "** Setting Constants ect **"
      print*, "***************************"

      n = 65
      m = 130
      
      print*, "** Adding in ghost cells **"
      print*, "***************************"
      
      !Move uin into a larger u vector to include ghost cells
      do i = 2,n-1
         do j = 2,m-2
            do k = 1,4
               u(k,i,j) = uin(k,i-1,j-1)
            enddo
         enddo
      enddo
       print*,"** Adding in more cells  **"
      print*, "***************************"
 
     !Include ghost cells for x and y vectors
      do i = 2,n
         do j = 2,m-1
         x(i,j) = xin(i-1,j-1)
         y(i,j) = yin(i-1,j-1)
      enddo
   enddo

   testVal = 10.0
    print*, "The value of testVal is :",testVal
    testVal = f(testVal)
    print*, "The value of testVal is :",testSol

   !ff = f(u(:,5,5))
!   print*, "the value of u is: ",u(:,5,5)
!   ff = f(u(:,5,5))
!   print*, "the value of ff is now: ",ff

      return 
    end subroutine jameson


function f(x)
    implicit none
    real(kind=8), intent(in) :: x
    real(kind=8), dimension(2) :: f
    f(1) = x
    f(2) = x**2
    print*, "The vector f is: ",f
end function f


!    function f(uu)
!       implicit none
!        real(kind=8), dimension(4), intent(in) :: uu
!        real(kind=8), dimension(4) :: f
!        real(kind=8) :: g, p
!        print*, "the value of u is: ",uu
!        g = 1.4d0
!!        p = (g-1)*(u(1)-(u(2)**2+u(3)**2)/(2*u(1)))

!!        print*, "The p value is: ",p
!        print*, "The g value is: ",g
!        f(:) = uu
!!        f(1) = u(2)
!!        f(2) = u(2)**2/u(1) + p 
!!        f(3) = u(2)*u(3)/u(1)
!!        f(4) = (u(4)+p)*(u(2)/u(1))
!!        print*, "The u vector is: ",u(1)
!!        print*, "The f vector is: ",f
!      end function f
