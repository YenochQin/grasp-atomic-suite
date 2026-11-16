!***********************************************************************
!                                                                      *
      subroutine cxk (s,is,kaps,nu,k,ibr,iex)
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      implicit none
!      DIMENSION IS(4),KAPS(4),S(12)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nu
      integer :: k
      integer, intent(in) :: ibr
      integer, intent(in) :: iex
      integer, intent(in) :: is(4)
      integer, intent(in) :: kaps(4)
      real(double) , intent(inout) :: s(12)
!-----------------------------------------------
!
      stop 'CXK: Error '
      end
