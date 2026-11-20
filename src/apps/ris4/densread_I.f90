      module densread_i
      implicit none (type, external)
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine densread (dint1,dint2,dint3,                          &
                           dint4,dint5,dint6,                          &
                           dint7)
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: max_orbitals, max_radial_points
      use selected_level_indices,           only: nvec
      implicit none (type, external)
      real(real64), intent(in) :: dint1(max_orbitals, max_orbitals), dint2(max_orbitals, max_orbitals),  &
                                          dint3(max_orbitals, max_orbitals), dint4(max_orbitals, max_orbitals), &
                                             dint5(max_orbitals, max_orbitals), dint6(max_orbitals, max_orbitals),  &
                                          dint7(max_orbitals, max_orbitals)
      end subroutine densread
      end interface
      end module densread_i
