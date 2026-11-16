      module buffer_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  13:03:28   1/25/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: angular_buffer_capacity, angular_coefficient_count
      integer, pointer :: angular_coefficient_labels(:, :)
      real(real64), pointer :: angular_coefficients(:)
      end module buffer_c
