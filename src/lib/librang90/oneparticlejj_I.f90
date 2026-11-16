      module oneparticlejj_i
      interface
!
      subroutine oneparticlejj(ka,iopar,ja,jb,ia1,ia2,vshell)
      use vast_kind_param, only:  double
      use parameter_def,   only:  nnnw
!GG      INTEGER NNNW
!GG      PARAMETER (NNNW = 214)
      integer, intent(in)  :: ka,iopar,ja,jb
      integer, intent(out) :: ia1,ia2
      real(double), dimension(nnnw), intent(out) :: vshell
      end subroutine
      end interface
      end module
