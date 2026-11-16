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
      use vast_kind_param, only:  double
!
      real(double), dimension(:,:), pointer :: channel_coefficients
      integer, dimension(:,:), pointer :: channel_csf_columns, channel_matrix_positions, channel_packed_orbital_labels
!
      end module sacoef_c
