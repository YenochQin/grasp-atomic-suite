      module posfile_i
      implicit none (type, external)
      interface
      subroutine posfile (mode, nunit, nrec)
      implicit none (type, external)
      integer, intent(in) :: mode
      integer, intent(in) :: nunit
      integer, intent(in) :: nrec
!...This routine performs I/O.
      end subroutine
      end interface
      end module
