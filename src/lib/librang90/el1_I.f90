      module el1_i
      implicit none (type, external)
      interface
!
      subroutine el1(jja,jjb,ja,jb,iire,icolbrei)
      implicit none (type, external)
      integer, intent(in) :: jja,jjb,ja,jb,iire,icolbrei
      end subroutine
      end interface
      end module
