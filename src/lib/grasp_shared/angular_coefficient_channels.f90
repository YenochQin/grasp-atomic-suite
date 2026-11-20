! Stores angular coefficients with their CSF columns, matrix positions, and orbital labels by channel.
      module angular_coefficient_channels
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      private :: real64
      real(real64), pointer :: channel_coefficients(:, :)
      integer, pointer :: channel_csf_columns(:, :), channel_matrix_positions(:, :), channel_packed_orbital_labels(:, :)
      end module angular_coefficient_channels
