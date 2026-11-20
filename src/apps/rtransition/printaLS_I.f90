      module printals_i
      implicit none (type, external)
      interface
      subroutine printals (inum_ii,inum_ff,asfa,asfb,i,j,omega,factor)
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer  :: i, j, inum_ii, inum_ff
      real(real64), intent(in) :: asfa, asfb
      real(real64), intent(in) :: omega, factor
      end subroutine printals
      end interface
      end module printals_i
