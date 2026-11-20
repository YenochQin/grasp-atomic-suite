!***********************************************************************
!                                                                      *
      subroutine alcnta(isldr, isldr1, xsldr, ntdim, imode)
!                                                                      *
!   This  subprogram allocates (IMODE = 1), reallocates (IMODE = 2),   *
!   and deallocates (IMODE = 3) storage for  certain arrays that are   *
!   used in OSCL.                                                      *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
!                                                                      *
!   Farid A. Parpia.                      Last revision: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use array_memory_management, only: alloc, dalloc, ralloc
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ntdim
      integer, intent(in) :: imode
      integer, pointer :: isldr(:), isldr1(:)
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
         ntdim = 1
!
!   Allocate storage for arrays
!
         call alloc (isldr, ntdim, "ISLDR", "ALCNTA" )
         call alloc (isldr1, ntdim,"ISLDR1", "ALCNTA" )
         call alloc (xsldr, ntdim, "XSLDR", "ALCNTA" )
!
      case (2)
!
!   Double the allocation of storage for the arrays
!
         newsiz = 2*ntdim
!
         call ralloc (isldr, newsiz, "ISLDR", "ALCNTA" )
         call ralloc (isldr1, newsiz, "ISLDR1", "ALCNTA" )
         call ralloc (xsldr, newsiz, "XSLDR", "ALCNTA" )
!
         ntdim = newsiz
!
      case (3)
!
!   Deallocate the storage for the arrays
!
         call dalloc (isldr, "ISLDR", "ALCNTA" )
         call dalloc (isldr1, "ISLDR1", "ALCNTA" )
         call dalloc (xsldr,  "XSLDR", "ALCNTA" )
!
      case default
!
         write (6, *) "ALCNTA: Invalid argument IMODE = ", imode
         stop
!
      end select
!
      return
      end subroutine alcnta
