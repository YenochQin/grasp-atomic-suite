      module hmout_i
      implicit none (type, external)
      interface
      subroutine hmout (myid, nprocs, ncf)
      implicit none (type, external)
      integer :: myid
!VAST...Dummy argument MYID is not referenced in this routine.
      integer :: nprocs
!VAST...Dummy argument NPROCS is not referenced in this routine.
      integer, intent(in) :: ncf
      end subroutine
      end interface
      end module
