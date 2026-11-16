!***********************************************************************
!                                                                      *
      subroutine alloccheck(n, irez)
!                                                                      *
!   Written by G. Gaigalas            Last revision: 27 October 2017   *
!                                                                      *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use memory_man
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)  :: n
      integer, intent(out) :: irez
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double), pointer, dimension(:) :: pr
      integer, pointer, dimension(:) :: pi
      integer :: error
!-----------------------------------------------
      irez = 1
      allocate(pr(n),stat=error)
      if (error /= 0) then
         irez = 0
         return
      end if
      allocate(pi(3*n),stat=error)
      if (error /= 0) then
         irez = 0
         call dalloc (pr,'PR','allocCheck')
         return
      end if
      call dalloc (pr,'PR','allocCheck')
      call dalloc (pi,'PI','allocCheck')
      return
      end subroutine alloccheck
