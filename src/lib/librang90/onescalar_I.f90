      module onescalar_i
      interface
!
      subroutine onescalar(ja,jb,ia1,ia2,vshell)
      use vast_kind_param, only:  double
      use parameter_def,   only:  max_orbitals
!GG      INTEGER NNNW
!GG      PARAMETER (NNNW = 214)
      integer, intent(in)  :: ja,jb
      integer, intent(out) :: ia1,ia2
      real(double), dimension(max_orbitals), intent(out) :: vshell
      end subroutine
      end interface
      end module
