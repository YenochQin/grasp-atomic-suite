      module starttime_i
      implicit none (type, external)
      interface
      subroutine starttime (ncount1, progname)
      implicit none (type, external)
      integer, intent(out) :: ncount1
      character (len = *), intent(in) :: progname
!VAST...Calls: DATE_AND_TIME, SYSTEM_CLOCK
!...This routine performs I/O.
      end subroutine
      end interface
      end module
