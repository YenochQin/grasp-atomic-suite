      module smsread_mpi_i
      implicit none (type, external)
      interface
      subroutine smsread_mpi (vint,vint2)
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters, only: max_orbitals
      implicit none (type, external)
      real(real64), intent(in) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
      end subroutine smsread_mpi
      end interface
      end module smsread_mpi_i
