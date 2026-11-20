      module setmcpmpi_i
      use csf_descriptor_reader, only: descriptor_source
      implicit none (type, external)
      interface
!cjb  myid, nprocs = NOT args
      subroutine setmcpmpi (core_subshell_count, block_labels, filehead, descriptors)
      import descriptor_source
      implicit none (type, external)
!     INTEGER , INTENT(IN) :: MYID
!     INTEGER , INTENT(IN) :: NPROCS
      integer , intent(in) :: core_subshell_count
      character (len = 8), intent(in) :: block_labels(*)
      character (len = *), intent(in) :: filehead
      type(descriptor_source), intent(in) :: descriptors
!...This routine performs I/O.
      end subroutine setmcpmpi
      end interface
      end module setmcpmpi_i
