      module el33_i
      implicit none (type, external)
      interface
!
      subroutine el33(jjja,jjjb,ja,jb,jc,irez,jja,jjb,jjc,jjd   &
                                                        ,icolbrei)
      implicit none (type, external)
      integer, intent(in) :: jjja,jjjb,ja,jb,jc,irez,jja,jjb,jjc,jjd &
                                                        ,icolbrei
      end subroutine
      end interface
      end module
