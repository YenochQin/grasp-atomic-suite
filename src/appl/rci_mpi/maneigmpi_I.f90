      module maneig_i
      implicit none (type, external)
      interface
      subroutine maneig (iatjpo, iaspar, nelmnt_a)
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none (type, external)
      integer, intent(out) :: iatjpo
      integer, intent(out) :: iaspar
      integer(int64)        :: nelmnt_a
      end subroutine
      end interface
      end module
