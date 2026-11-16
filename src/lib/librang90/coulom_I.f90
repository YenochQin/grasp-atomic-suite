      module coulom_i
      implicit none (type, external)
      interface
!
      subroutine coulom(j1,j2,j3,j4,l1,l2,l3,l4,k,aa)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: j1, j2, j3, j4, l1, l2, l3, l4, k
      real(real64),intent(out)  :: aa
      end subroutine
      end interface
      end module
