      module talk_i
      implicit none (type, external)
      interface
      subroutine talk (ja, jb, nu, ia, ib, ic, id, itype, coef)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ja
      integer, intent(in) :: jb
      integer, intent(in) :: nu
      integer, intent(in) :: ia
      integer, intent(in) :: ib
      integer, intent(in) :: ic
      integer, intent(in) :: id
      integer, intent(in) :: itype
      real(real64), intent(in) :: coef
      end subroutine talk
      end interface
      end module talk_i
