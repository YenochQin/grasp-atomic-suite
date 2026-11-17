      module lodcslbio_i
      implicit none (type, external)
      interface
      subroutine lodcslbio (ncore,igg)
      implicit none (type, external)
      integer, intent(in) ::  igg
      integer, intent(out) :: ncore
      end subroutine lodcslbio
      end interface
      end module lodcslbio_i
