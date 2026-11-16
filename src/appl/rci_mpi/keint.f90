!***********************************************************************
!                                                                      *
      subroutine keint (ia,ib,tegral)
!                                                                      *
!   This routine returns kinetic energy integrals.                     *
!      A list of previously computed integrals is maintained. If the   *
!   required  integral is not already  available it is computed  and   *
!   stored in the list.                                                *
!      Note that indices IA and IB may be interchanged on output.      *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, RALLOC, RINTI.                         *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: orbital_index_radix
      use memory_man
      use keilst_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rinti_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(inout) :: ia, ib
      real(real64), intent(out) :: tegral
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical ::  found
      integer ::  i, index, loc, ju, jl, jm, key, newsiz, iswap
!-----------------------------------------------
!
!   Ensure that the indices are in `canonical' order
!
      if (ia .gt. ib) then
         iswap = ib
         ib = ia
         ia = iswap
      endif
!
!   Compute the composite (packed) index
!
      index = ia*orbital_index_radix+ib
!
      if (.not. frstki) then
!
!   This branch is executed on all entries except the first
!
         if (index .gt. indkei(nkei)) then
!
!   The index is greater than the largest stored
!
            found = .false.
            loc = nkei
!
         elseif (index .lt. indkei(1)) then
!
!   The index is less than the smallest stored
!
            found = .false.
            loc = 0
!
         else
!
!   The index is within the range of the indices stored; search
!   for it in the list of indices
!
            ju = nkei
            jl = 1
!
    1       if (ju-jl .gt. 1) then
!
               jm = (ju+jl)/2
               if (indkei(jm) .gt. index) then
                  ju = jm
               else
                  jl = jm
               endif
               goto 1
!
            else
!
!   The range is bracketed to the extent possible
!
               if (index .eq. indkei(ju)) then
!
                  found = .true.
                  loc = ju
!
               elseif (index .eq. indkei(jl)) then
!
                  found = .true.
                  loc = jl
!
               else
!
                  found = .false.
                  loc = jl
!
               endif
!
            endif
!
         endif
!
         if (found) then
!
!   Found the index in the list; return the value of the integral
!   from storage
!
            tegral = valkei(loc)
!
         else
!
!   Index not found; compute the integral
!
            tegral = rinti (ia,ib,1)
!
!   Increment the integral counter
!
            nkei = nkei+1
!
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
            if (nkei .gt. ndkea) then
               newsiz = ndkea+ndkea/2
               call ralloc (indkei,newsiz,'INDKEI', 'KEINT')
               call ralloc (valkei,newsiz,'VALKEI', 'KEINT')
               ndkea = newsiz
            endif
!
            do 2 i = nkei,loc+2,-1
               indkei(i) = indkei(i-1)
               valkei(i) = valkei(i-1)
    2       continue
!
!   Put the new index and value into storage
!
            indkei(loc+1) = index
            valkei(loc+1) = tegral
!
         endif
!
      else
!
!   This branch is executed only once
!
         frstki = .false.
!
!   Allocate the initial storage for arrays INDOEI and VALOEI;
!   NDCOEA stores the array dimension
!
         ndkea = 2
         call alloc (indkei,ndkea,'INDKEI','KEINT' )
         call alloc (valkei,ndkea,'VALKEI','KEINT')
!
!   Compute the integral's value
!
         tegral = rinti (ia,ib,1)
!
!   Initialise the integral counter
!
         nkei = 1
!
!   Store the integral and its value
!
         indkei(1) = index
         valkei(1) = tegral
!
      endif
!
      return
      end subroutine keint
