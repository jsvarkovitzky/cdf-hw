! File Jameson.f90
      subroutine jameson (xin,yin,sxx,sxy,syx,syy,uin,b)
      real(kind=8), dimension(64,128), intent (in) :: xin,yin,sxx&
&,sxy,syx,syy
      real(kind=8), dimension(4,64,128), intent (in) :: uin
      real(kind=8), dimension(128,64), intent (out) :: b
      integer :: n,m
      real(kind=8), dimension(64,128) :: x,y
      real(kind=8), dimension(4,64,128) :: u
      external :: f
      real(kind=8) :: f
      print*, "***************************"
      print*, "** Begining Jameson.f90! **"
      print*, "***************************"

      print*, "** Setting Constants ect **"
      print*, "***************************"

      n = 65
      m = 130
      
      print*, "** Adding in ghost cells **"
      print*, "***************************"
      

      do i = 2,n
         do j = 2,m-1
        do k = 1,4
           u(k,i,j) = uin(k,i-1,j-1)
        enddo
     enddo
  enddo
      do i = 2,n
         do j = 2,m-1
         x(i,j) = xin(i-1,j-1)
         y(i,j) = yin(i-1,j-1)
      enddo
   enddo

      return 
    end subroutine jameson

     real(kind=8) function f(x)
        implicit none
        real(kind=8), intent(in) :: x
        f = x**2
        print *, "The value of f is: ", f
      end function f

