!
!***********************************************************************
!                                                                      *
      module eigv_c
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
      real(real64) :: hamiltonian_energy_shift, eavff, eavii
      real(real64), pointer :: shifted_eigenvalues(:), evalff(:), evalii(:)
      real(real64), pointer :: csf_mixing_coefficients(:), evecff(:), evecii(:)
      end module eigv_c
