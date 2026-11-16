      module sortmem_i
      interface
      subroutine sortmem (nfile,iggmax_k,ncoeff,ntgral,lprint,nb,fhead)
      integer, intent(in)  :: nfile
      integer, intent(in)  :: ncoeff
      integer, intent(out) :: ntgral
      logical, intent(in)  :: lprint
      integer, intent(in)  :: nb,iggmax_k
      character (len = *), intent(in) :: fhead
      end subroutine
      end interface
      end module
