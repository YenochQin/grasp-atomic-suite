      module rumtjj_i
      implicit none (type, external)
      interface
!
      subroutine rumtjj(knt,jj,lq,lv,l)
      implicit none (type, external)
      integer, intent(in)  :: jj, knt
      integer, intent(out) :: lq, lv, l
      end subroutine
      end interface
      end module
