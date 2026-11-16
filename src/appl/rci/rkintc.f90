!***********************************************************************
!                                                                      *
      subroutine rkintc(ia, ib, ic, id, k, tegral)
!                                                                      *
!                         k                                            *
!   This routine returns R (abcd) integrals.                           *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use cteilsrk_c
      use orb_c
      use kkstart_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(inout) :: ia, ib, ic, id
      integer, intent(in) :: k
      real(real64), intent(out) :: tegral
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: key, iswap, index, jl, ju, jm, loc
      logical :: found, first
!-----------------------------------------------
!
      key = orbital_count + 1
!
!   Ensure that the indices are in `canonical' order
!   Compute the composite (packed) index
!


      if (ia > ic) then
         iswap = ic
         ic = ia
         ia = iswap
      endif
      if (ib > id) then
         iswap = id
         id = ib
         ib = iswap
      endif
      if (ia > ib) then
         iswap = ib
         ib = ia
         ia = iswap
         iswap = id
         id = ic
         ic = iswap
      endif

      index = ((ia*key + ib)*key + ic)*key + id
!
      jl = kstart(k)
      ju = kstart(k+1) - 1

      if (index<indteirk(jl) .or. index>indteirk(ju)) then
         write (*, *) 'Something wrong in rkintc'
         stop
      endif
!
!   The index is within the range of the indices stored; search
!   for it in the list of indices
!
    1 continue
      if (ju - jl > 1) then
         jm = (ju + jl)/2
         if (indteirk(jm) > index) then
            ju = jm
         else
            jl = jm
         endif
         go to 1
      endif
!
!   The range is bracketed to the extent possible
!
      if (index == indteirk(ju)) then
         loc = ju
      else if (index == indteirk(jl)) then
         loc = jl
      else
         write (*, *) k, ia, ib, ic, id, index
         write (*, *) 'Rkintc Integral not found'
         stop
      endif
!
!   Return the value of the integral
!   from storage

      tegral = valteirk(loc)
!
      return
      end subroutine rkintc
