      module oscl_i
      implicit none (type, external)
      interface
      subroutine oscl (name,fullname,tmpdir,startup_directory,block_label_text)
      implicit none (type, external)
      character (len = 24), intent(inout) :: name(2)
      character(len=3)  ::  block_label_text
      character (len = 128), intent(inout) :: fullname(2)
      character(len=128) :: startup_directory,tmpdir
      end subroutine oscl
      end interface
      end module oscl_i
