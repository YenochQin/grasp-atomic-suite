      module oneparticlejj_i
      implicit none (type, external)
      interface
!
      subroutine oneparticlejj(ka,iopar,ja,jb,ia1,ia2,vshell)
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only:  max_orbitals
      implicit none (type, external)
!GG      INTEGER NNNW
!GG      PARAMETER (NNNW = 214)
      integer, intent(in)  :: ka,iopar,ja,jb
      integer, intent(out) :: ia1,ia2
      real(real64), intent(out) :: vshell(max_orbitals)
      end subroutine
      end interface
      end module
