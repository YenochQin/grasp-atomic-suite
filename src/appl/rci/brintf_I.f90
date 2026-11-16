      module brintf_i
      implicit none (type, external)
      interface
      real(real64) function brintf (itype, ia, ib, ic, id, k)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: itype
      integer, intent(in) :: ia
      integer, intent(in) :: ib
      integer, intent(in) :: ic
      integer, intent(in) :: id
      integer, intent(in) :: k
      end function brintf
      end interface
      end module brintf_i
