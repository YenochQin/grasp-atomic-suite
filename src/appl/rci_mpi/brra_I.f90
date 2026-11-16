      module brra_i
      implicit none (type, external)
      interface
      real(real64) function brra (itype, ia, ic, ib, id, k)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: itype
      integer, intent(in) :: ia
      integer, intent(in) :: ic
      integer, intent(in) :: ib
      integer, intent(in) :: id
      integer :: k
      end function
      end interface
      end module
