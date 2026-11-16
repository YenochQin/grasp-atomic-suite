      module qqsort_i
      implicit none (type, external)
      interface
      subroutine qqsort (nfile, number, kstart, name, kamax)
      implicit none (type, external)
      integer, intent(in) :: nfile
      integer, intent(in) :: number
      integer :: kstart
!...Dummy argument KSTART is not referenced in this routine.
      character (len = 24) :: name
!...Dummy argument NAME is not referenced in this routine.
      integer, intent(in) :: kamax
      end subroutine qqsort
      end interface
      end module qqsort_i
