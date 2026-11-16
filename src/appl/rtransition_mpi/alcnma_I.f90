      module alcnma_i
      implicit none (type, external)
      interface
     subroutine alcnma(iptr, isldr, isldr1, xsldr, nmdim, imode)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, pointer :: iptr(:), isldr(:), isldr1(:)
      real(real64), pointer :: xsldr(:)
      integer, intent(inout) :: nmdim
      integer, intent(in) :: imode
!...This routine performs I/O.
      end subroutine
      end interface
      end module
