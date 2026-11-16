!
!***********************************************************************
!                                                                      *
      module eigv_c
!                                                                      *
!***********************************************************************
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      real(double) :: hamiltonian_energy_shift, eavff, eavii
      real(double), dimension(:), pointer :: shifted_eigenvalues, evalff, evalii
      real(double), dimension(:), pointer :: csf_mixing_coefficients, evecff, evecii
      end module eigv_c
