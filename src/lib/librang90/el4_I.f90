      module el4_i
      implicit none (type, external)
      interface
!
      subroutine el4(jja,jjb,ja,jb,jc,jd,icolbrei)
      implicit none (type, external)
      integer, intent(in) :: jja,jjb,ja,jb,jc,jd,icolbrei
      end subroutine
      end interface
      end module
