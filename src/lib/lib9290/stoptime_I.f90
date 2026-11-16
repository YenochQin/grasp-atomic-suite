      module stoptime_i
      implicit none (type, external)
      interface
      subroutine stoptime (ncount1, progname)
      implicit none (type, external)
      integer, intent(in) :: ncount1
      character (len = *), intent(in) :: progname
!VAST...Calls: SYSTEM_CLOCK, DATE_AND_TIME
!...This routine performs I/O.
      end subroutine
      end interface
      end module
