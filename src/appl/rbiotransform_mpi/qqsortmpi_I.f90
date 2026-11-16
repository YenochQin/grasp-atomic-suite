      module qqsort_i
      implicit none (type, external)
      interface
      subroutine qqsort (nfile, number, kstart, name, kamax)
      implicit none (type, external)
      integer, intent(in) :: nfile
      integer, intent(in) :: number
      integer :: kstart
      character (len = 24) :: name
      integer, intent(in) :: kamax
      end subroutine
      end interface
      end module
