C File Jameson.f90
      subroutine jameson (x,y,sxx,sxy,syx,syy,b)
      real(kind=8), dimension(64,128), intent (in) :: x,y,sxx,sxy,syx,syy
      real(kind=8), dimension(128,64), intent (out) :: b
      print*, "Hello from Fortran!"

      b = transpose(x)
      
      end
