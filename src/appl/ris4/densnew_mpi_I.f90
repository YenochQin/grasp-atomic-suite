      module densnew_mpi_i
      implicit none (type, external)
      interface
      subroutine densnew_mpi (doit,dint1,dint2,dint3,dint4,dint5,dint6,dint7)
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def, only: max_orbitals
      implicit none (type, external)
      integer, intent(in) :: doit
      real(real64), intent(in) :: dint1(max_orbitals, max_orbitals), dint2(max_orbitals, max_orbitals), &
                                                        dint3(max_orbitals, max_orbitals), dint4(max_orbitals, max_orbitals), &
                                                        dint5(max_orbitals, max_orbitals), dint6(max_orbitals, max_orbitals), &
                                                        dint7(max_orbitals, max_orbitals)
      end subroutine
      end interface
      end module
