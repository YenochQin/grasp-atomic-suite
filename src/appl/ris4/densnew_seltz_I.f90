      module densnew_seltz_i
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine densnew_seltz (doit,dint1,dint2,dint3,dint4,dint5,dint6,    &
                          dint7,dint1vec,dens1vec,nrnuc)
      use vast_kind_param, only: double
      use parameter_def,   only: nnnw, nnnp
      use prnt_c,          only : nvec
      real(double), dimension(nnnw,nnnw), intent(in) :: dint1
      real(double), dimension(nnnw,nnnw), intent(in) :: dint2
      real(double), dimension(nnnw,nnnw), intent(in) :: dint3
      real(double), dimension(nnnw,nnnw), intent(in) :: dint4
      real(double), dimension(nnnw,nnnw), intent(in) :: dint5
      real(double), dimension(nnnw,nnnw), intent(in) :: dint6
      real(double), dimension(nnnw,nnnw), intent(in) :: dint7
      integer, intent(in) :: nrnuc
      real(double), dimension(nvec,nrnuc), intent(out) :: dens1vec
      real(double), dimension(nnnw,nnnw,nrnuc), intent(in) :: dint1vec
      integer, intent(in) :: doit
      end subroutine
      end interface
      end module
