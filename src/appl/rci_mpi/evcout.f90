!***********************************************************************
!                                                                      *
      subroutine evcout
!                                                                      *
!   Routine for printing eigenvectors.                                 *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use memory_man
      use eigv_c
      use orb_c,           only: csf_count, orbital_count, csf_subshell_occupations
      use prnt_c
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V A R I A B L E S
!-----------------------------------------------
      integer :: i, j1, j2, nvex, jp, l, m
      integer, pointer :: iloc(:)
!-----------------------------------------------
!
!   Header
!
      write (24,300)
!
!   Determine which eigenvectors are to be printed out
!
      nvex = nvec
      call alloc (iloc,nvec,'ILOC', 'EVCOUT')
      do 1 i = 1,nvec
         iloc(i) = ivec(i)
    1 continue
!
!   There are eight columns across the width of a page
!
      jp = 8
!
!   Set up for the first set of columns
!
      j1 = 1
      j2 = jp
!
!   Loop over sets of columns
!
    3 if (j2 .gt. nvex) j2 = nvex
!
!   Write out the column numbers; skip a line
!
      write (24,301) (iloc(l),l = j1,j2)
      write (24,302)
!
!   Write out the rows
!
      do 4 m = 1,csf_count
!
         write (24,303) (csf_mixing_coefficients(m+(l-1)*csf_count),l = j1,j2)
!
!   Skip a line after every tenth row
!
         if (mod(m,10) .eq. 0) write (24,302)
!
    4 continue
!
!   Set up for the next set of columns
!
      write (24,303)
      if (j2 .lt. nvex) then
         j1 = j1+jp
         j2 = j2+jp
         goto 3
      endif
!
!   Deallocate storage for ILOC
!
      call dalloc (iloc, 'ILOC', 'EVCOUT')
!
      return
!
  300 format (//' Eigenvectors:'/)
  301 format (1x,8(i8,7x))
  302 format (1x)
  303 format (1x,1p,8d15.7)
!
      end
