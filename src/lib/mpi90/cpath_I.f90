       module cpath_i
       implicit none (type, external)
       interface
!
       subroutine cpath(startdir, permdir, tmpdir)
       implicit none (type, external)
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
       character(len=*), intent(in) :: startdir, permdir, tmpdir
       end subroutine
       end interface
       end module
