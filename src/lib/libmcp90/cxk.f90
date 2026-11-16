module cxk_i
    implicit none (type, external)
    private
    public :: cxk
contains
!***********************************************************************
!                                                                      *
      subroutine cxk (s,is,kaps,nu,k,ibr,iex)
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!      DIMENSION IS(4),KAPS(4),S(12)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nu
      integer, intent(in) :: k
      integer, intent(in) :: ibr
      integer, intent(in) :: iex
      integer, intent(in) :: is(4)
      integer, intent(in) :: kaps(4)
      real(real64) , intent(inout) :: s(12)
!-----------------------------------------------
!
      stop 'CXK: Error '
      end
end module cxk_i
