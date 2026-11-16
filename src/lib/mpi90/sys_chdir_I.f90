      module sys_chdir_i
      implicit none (type, external)
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine sys_chdir (dir, lendir, ierr)
      implicit none (type, external)
      character(len=*), intent(in):: dir
      integer, intent(in):: lendir
      integer, intent(out):: ierr
      end subroutine
      end interface
      end module
