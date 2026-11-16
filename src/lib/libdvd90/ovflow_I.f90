      module ovflow_i
      implicit none (type, external)
      interface
      subroutine ovflow (n, nume, kpass, scra1, basis, ab, s, svec, eigval)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: n
      integer, intent(in) :: nume
      integer, intent(inout) :: kpass
!CFF   .. added dimension to scra1
      real(real64) :: scra1(nume)
      real(real64) :: svec(kpass*nume)
      real(real64) :: basis(n*kpass)
      real(real64) :: ab(n*kpass)
      real(real64), intent(out) :: s((kpass*(kpass + 1))/2)
      real(real64), intent(in) :: eigval(kpass)
!VAST...Calls: MULTBC, DINIT
      end subroutine
      end interface
      end module
