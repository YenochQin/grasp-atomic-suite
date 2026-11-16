      module hmoutmpi_i
      implicit none (type, external)
      interface
      subroutine hmoutmpi (myid, nprocs, ncf)
      implicit none (type, external)
      integer :: myid
!VAST...Dummy argument MYID is not referenced in this routine.
      integer :: nprocs
!VAST...Dummy argument NPROCS is not referenced in this routine.
      integer, intent(in) :: ncf
      end subroutine hmoutmpi
      end interface
      end module hmoutmpi_i
