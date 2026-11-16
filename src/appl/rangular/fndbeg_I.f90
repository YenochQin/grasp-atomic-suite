      module fndbeg_i
      implicit none (type, external)
      interface
      subroutine fndbeg (jastrt, jbstrt, index, llistt, llistv)
      implicit none (type, external)
      integer, intent(inout) :: jastrt
      integer, intent(inout) :: jbstrt
      integer, intent(inout) :: index
      integer, intent(out) :: llistt
      integer, pointer, intent(in) :: llistv(:)
      end subroutine
      end interface
      end module
