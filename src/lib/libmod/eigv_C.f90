!
!***********************************************************************
!                                                                      *
      module eigv_c
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      real(real64) :: hamiltonian_energy_shift, eavff, eavii
      real(real64), pointer :: shifted_eigenvalues(:), evalff(:), evalii(:)
      real(real64), pointer :: csf_mixing_coefficients(:), evecff(:), evecii(:)
      end module eigv_c
