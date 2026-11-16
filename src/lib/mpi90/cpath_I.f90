       module cpath_i
       interface
!
       subroutine cpath(startdir, permdir, tmpdir)
!      USE vast_kind_param, ONLY:  DOUBLE
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
       character(len=*), intent(in) :: startdir, permdir, tmpdir
       end subroutine
       end interface
       end module
