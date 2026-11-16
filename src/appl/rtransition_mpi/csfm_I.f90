      module csfm_i
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine csfm (asfa,asfb,lev1,lev2)
      use vast_kind_param, only:  double
      real(double), intent(out) :: asfa, asfb
      integer, intent(in) :: lev1, lev2
      end subroutine
      end interface
      end module
