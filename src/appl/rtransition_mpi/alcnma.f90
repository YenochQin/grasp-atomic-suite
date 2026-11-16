!***********************************************************************
!                                                                      *
      subroutine alcnma(iptr, isldr, isldr1, xsldr, nmdim, imode)
!                                                                      *
!   This  subprogram allocates (IMODE = 1), reallocates (IMODE = 2),   *
!   and deallocates (IMODE = 3) storage for  certain arrays that are   *
!   local to SUBROUTINE TRSORT.                                        *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use memory_man
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: nmdim
      integer , intent(in) :: imode
      integer, pointer :: iptr(:), isldr(:), isldr1(:)
      real(real64), pointer :: xsldr(:)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: newsiz
!-----------------------------------------------
!
      select case (imode)
      case (1)
!
!   Initial array dimension
!
         nmdim = 64
!
!   Allocate storage for arrays
!
         call alloc (iptr, nmdim, 'IPTR', 'ALCNMA')
         call alloc (isldr, nmdim, 'ISLDR', 'ALCNMA')
         call alloc (isldr1, nmdim, 'ISLDR1', 'ALCNMA')
         call alloc (xsldr, nmdim, 'XSLDR', 'ALCNMA')
!
      case (2)
!
!   Double the allocation of storage for the arrays
!
         newsiz = 2*nmdim
!
         call ralloc (iptr, newsiz, 'IPTR', 'ALCNMA')
         call ralloc (isldr, newsiz, 'ISLDR', 'ALCNMA')
         call ralloc (isldr1, newsiz, 'ISLDR1', 'ALCNMA')
         call ralloc (xsldr, newsiz, 'XSLDR', 'ALCNMA')
!
         nmdim = newsiz
!
      case (3)
!
!   Deallocate the storage for the arrays
!
         call dalloc (iptr,  'IPTR', 'ALCNMA')
         call dalloc (isldr, 'ISLDR', 'ALCNMA')
         call dalloc (isldr1, 'ISLDR1', 'ALCNMA')
         call dalloc (xsldr,  'XSLDR', 'ALCNMA')
!
      case default
!
         write (6, *) 'ALCNMA: Invalid argument IMODE = ', imode
         stop
!
      end select
!
      return
      end subroutine alcnma
