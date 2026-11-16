      module oscl_i
      implicit none (type, external)
      interface
      subroutine oscl (name,fullname,tmpdir,startdir,idstring)
      implicit none (type, external)
      character (len = 24), intent(inout) :: name(2)
      character(len=3)  ::  idstring
      character (len = 128), intent(inout) :: fullname(2)
      character(len=128) :: startdir,tmpdir
      end subroutine
      end interface
      end module
