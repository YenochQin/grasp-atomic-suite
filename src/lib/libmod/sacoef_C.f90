!
!***********************************************************************
!                                                                      *
      module sacoef_c
!                                                                      *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2017   *
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!
      real(real64), pointer :: channel_coefficients(:, :)
      integer, pointer :: channel_csf_columns(:, :), channel_matrix_positions(:, :), channel_packed_orbital_labels(:, :)
!
      end module sacoef_c
