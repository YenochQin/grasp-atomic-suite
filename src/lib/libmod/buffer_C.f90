      module buffer_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
      integer :: angular_buffer_capacity, angular_coefficient_count
      integer, pointer :: angular_coefficient_labels(:, :)
      real(real64), pointer :: angular_coefficients(:)
      end module buffer_c
