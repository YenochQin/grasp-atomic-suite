      module densread_seltz_mpi_i
      implicit none (type, external)
      interface
      subroutine densread_seltz_mpi (dint1,dint2,dint3,dint4,dint5,dint6, &
         dint7,dint1vec,dens1vec,nrnuc)
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def, only: max_orbitals
      use prnt_c, only: nvec
      implicit none (type, external)
      integer, intent(in) :: nrnuc
      real(real64), intent(in) :: dint1(max_orbitals, max_orbitals), dint2(max_orbitals, max_orbitals),  &
                                          dint3(max_orbitals, max_orbitals), dint4(max_orbitals, max_orbitals), &
                                             dint5(max_orbitals, max_orbitals), dint6(max_orbitals, max_orbitals),  &
                                          dint7(max_orbitals, max_orbitals)
      real(real64), intent(out) :: dens1vec(nvec, nrnuc)
      real(real64), intent(in) :: dint1vec(max_orbitals, max_orbitals, nrnuc)
      end subroutine densread_seltz_mpi
      end interface
      end module densread_seltz_mpi_i
