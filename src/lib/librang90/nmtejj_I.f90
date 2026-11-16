      module nmtejj_i
      implicit none (type, external)
      interface
!
      integer function nmtejj(i2q,i2j,j,nk,nd)
      implicit none (type, external)
      integer, intent(in) :: i2q, i2j, j, nk, nd
      end function
      end interface
      end module
