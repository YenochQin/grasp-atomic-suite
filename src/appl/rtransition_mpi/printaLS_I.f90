      module printals_i
      interface
      subroutine printals (inum_ii,inum_ff,asfa,asfb,i,j,omega,factor)
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      use vast_kind_param,only: double
      integer  :: i, j, inum_ii, inum_ff
      real(double), intent(in) :: asfa, asfb
      real(double), intent(in) :: omega, factor
      end subroutine
      end interface
      end module
