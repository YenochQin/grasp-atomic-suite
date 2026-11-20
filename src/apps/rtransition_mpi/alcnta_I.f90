      module alcnta_i
      implicit none (type, external)
      interface
      subroutine alcnta(isldr, isldr1, xsldr, ntdim, imode)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, pointer :: isldr(:), isldr1(:)
      real(real64), pointer :: xsldr(:)
      integer, intent(inout) :: ntdim
      integer, intent(in) :: imode
!...This routine performs I/O.
      end subroutine alcnta
      end interface
      end module alcnta_i
