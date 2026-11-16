      module setsda_i
      implicit none (type, external)
      interface
      subroutine setsda (write_coefficient_callback, nnonz, lprint, nb, myid, nprocs, fhead)
      use outsdampi_i, only: output_callback => outsdampi
      implicit none (type, external)
      procedure(output_callback) :: write_coefficient_callback
      integer, intent(in) :: nnonz
      logical             :: lprint
      integer, intent(in) :: nb
      integer, intent(in) :: myid
      integer, intent(in) :: nprocs
      character (len = *), intent(in) :: fhead
      end subroutine setsda
      end interface
      end module setsda_i
