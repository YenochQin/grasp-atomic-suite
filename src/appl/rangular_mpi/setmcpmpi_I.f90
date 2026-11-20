      module setmcpmpi_i
      use csf_descriptor_reader, only: descriptor_source
      implicit none (type, external)
      interface
!cjb  myid, nprocs = NOT args
      subroutine setmcpmpi (ncore, idblk, filehead, descriptors)
      import descriptor_source
      implicit none (type, external)
!     INTEGER , INTENT(IN) :: MYID
!     INTEGER , INTENT(IN) :: NPROCS
      integer , intent(in) :: ncore
      character (len = 8), intent(in) :: idblk(*)
      character (len = *), intent(in) :: filehead
      type(descriptor_source), intent(in) :: descriptors
!...This routine performs I/O.
      end subroutine setmcpmpi
      end interface
      end module setmcpmpi_i
