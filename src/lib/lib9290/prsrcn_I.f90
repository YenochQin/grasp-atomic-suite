      module prsrcn_i
      implicit none (type, external)
      interface
      subroutine prsrcn (record, ncore, ioccs, ierr)
      use parameter_def, only: max_orbitals
      implicit none (type, external)
      character (len = 256), intent(in) :: record
      integer, intent(in) :: ncore
      integer, intent(out) :: ioccs(max_orbitals)
      integer, intent(out) :: ierr
!VAST.../ORB2/ NW(IN)
!VAST.../ORB4/ NP(IN)
!VAST.../ORB5/ NKJ(IN)
!VAST.../ORB10/ NH(IN)
!...This routine performs I/O.
      end subroutine
      end interface
      end module
