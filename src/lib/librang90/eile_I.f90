      module eile_i
      implicit none (type, external)
      interface
!                                                                  *
      subroutine eile(ja,jb,jc,jaa,jbb,jcc)
      implicit none (type, external)
      integer, intent(in)  :: ja, jb, jc
      integer, intent(out) :: jaa, jbb, jcc
      end subroutine
      end interface
      end module
