      module ris_cal_mpi_i
      implicit none (type, external)
      interface
      subroutine ris_cal_mpi (name)
      implicit none (type, external)
      character(len=24), intent(in) :: name
      end subroutine
      end interface
      end module
