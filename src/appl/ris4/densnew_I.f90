      module densnew_i
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine densnew (doit,dint1,dint2,dint3,dint4,dint5,dint6,    &
                          dint7)
      use vast_kind_param, only: double
      use parameter_def,   only: max_orbitals, max_radial_points
      use prnt_c,          only : nvec
      real(double), dimension(max_orbitals,max_orbitals), intent(in) :: dint1
      real(double), dimension(max_orbitals,max_orbitals), intent(in) :: dint2
      real(double), dimension(max_orbitals,max_orbitals), intent(in) :: dint3
      real(double), dimension(max_orbitals,max_orbitals), intent(in) :: dint4
      real(double), dimension(max_orbitals,max_orbitals), intent(in) :: dint5
      real(double), dimension(max_orbitals,max_orbitals), intent(in) :: dint6
      real(double), dimension(max_orbitals,max_orbitals), intent(in) :: dint7
      integer, intent(in) :: doit
      end subroutine
      end interface
      end module
