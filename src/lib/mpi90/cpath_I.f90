       module cpath_i
       implicit none (type, external)
       interface
!
       subroutine cpath(startdir, permdir, tmpdir)
       implicit none (type, external)
       character(len=*), intent(in) :: startdir, permdir, tmpdir
       end subroutine
       end interface
       end module
