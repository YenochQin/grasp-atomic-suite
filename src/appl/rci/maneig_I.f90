      module maneig_i
      implicit none (type, external)
      interface
      subroutine maneig (iatjpo, iaspar)
      implicit none (type, external)
      integer, intent(out) :: iatjpo
      integer, intent(out) :: iaspar
      end subroutine
      end interface
      end module
