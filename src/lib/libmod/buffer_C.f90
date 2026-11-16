      module buffer_c
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  13:03:28   1/25/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: angular_buffer_capacity, angular_coefficient_count
      integer, dimension(:,:), pointer :: angular_coefficient_labels
      real(double), dimension(:), pointer :: angular_coefficients
      end module buffer_c
