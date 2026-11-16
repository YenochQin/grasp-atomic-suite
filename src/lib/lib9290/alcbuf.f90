!***********************************************************************
!                                                                      *
      subroutine alcbuf(mode)
!                                                                      *
!   The arrays in  MODULE BUFFER_C are allocated (MODE = 1), realloca-   *
!   ted (MODE = 2), and deallocated (MODE = 3) in this routine.        *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
!                                                                      *
!   Written by Farid A. Parpia.             Last update: 06 Oct 1992   *
!   Editted after translation by C. Froese Fischer       29 Apr 2007   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:02:30   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!====================================================================

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use memory_man
      use buffer_c,       only: label, coeff, nbdim
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: mode
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: newsiz, ierr
!
!
      select case (mode)
      case (1)
!
!   Allocation
!
         nbdim = 10
         call alloc (label, 6,nbdim, 'LABEL', 'ALCBUF')
         call alloc (coeff, nbdim, 'COEFF', 'ALCBUF')
!
      case (2)
!
!   Reallocation to double storage
!
         newsiz = nbdim + nbdim
         call ralloc(label,6,newsiz, 'LABEL', 'ALCBUF')
         call ralloc(coeff,newsiz, 'COEFF', 'ALCBUF')
         nbdim = newsiz
!
      case (3)
!
!   Deallocation
!
          call dalloc(label, 'LABEL', 'ALCBUF')
          call dalloc(coeff, 'COEFF', 'ALCBUF')
!
      case default
!
!   Argument error
!
         write (6, *) 'ALCBUF: Invalid argument MODE = ', mode
         stop
!
      end select
!
      return
!
      end subroutine alcbuf
