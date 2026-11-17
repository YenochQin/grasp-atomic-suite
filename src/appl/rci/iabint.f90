!***********************************************************************
!                                                                      *
      subroutine iabint (ia,ib,tegral)
!                                                                      *
!   This routine returns I(ab) integrals.                              *
!      A list of previously computed integrals is maintained. If the   *
!   required  integral is not already  available it is computed  and   *
!   stored in the list.                                                *
!      Note that indices IA and IB may be interchanged on output.      *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, RALLOC, RINTI                          *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: orbital_index_radix
      use memory_man, only: alloc, ralloc, double
      use coeils_c, only: double, frstco, indoei, ncoei, ndcoea, valoei
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rinti_i, only: rinti
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
      if (ia > ib) then
         iswap = ib
         ib = ia
         ia = iswap
      end if
!
!   Compute the composite (packed) index
!
      index = ia*orbital_index_radix+ib
!
      if (.not. frstco) then
!
!   This branch is executed on all entries except the first
!
         if (index > indoei(ncoei)) then
!
!   The index is greater than the largest stored
!
            found = .false.
            loc = ncoei
!
         else if (index < indoei(1)) then
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
            ju = ncoei
            jl = 1
!
    1       if (ju-jl > 1) then
!
               jm = (ju+jl)/2
               if (indoei(jm) > index) then
                  ju = jm
               else
                  jl = jm
               end if
               goto 1
!
            else
!
!   The range is bracketed to the extent possible
!
               if (index == indoei(ju)) then
!
                  found = .true.
                  loc = ju
!
               else if (index == indoei(jl)) then
!
                  found = .true.
                  loc = jl
!
               else
!
                  found = .false.
                  loc = jl
!
               end if
!
            end if
!
         end if
!
         if (found) then
!
!   Found the index in the list; return the value of the integral
!   from storage
!
            tegral = valoei(loc)
!
         else
!
!   Index not found; compute the integral
!
            tegral = rinti (ia,ib,0)
!
!   Increment the integral counter
!
            ncoei = ncoei+1
!
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
            if (ncoei > ndcoea) then
               newsiz = ndcoea+ndcoea/2
               call ralloc (indoei,newsiz,"INDOEI","IABINT" )
               call ralloc (valoei,newsiz,"VALOEI","IABINT")
               ndcoea = newsiz
            end if
!
            do 2 i = ncoei,loc+2,-1
               indoei(i) = indoei(i-1)
               valoei(i) = valoei(i-1)
    2       continue
!
!   Put the new index and value into storage
!
            indoei(loc+1) = index
            valoei(loc+1) = tegral
!
         end if
!
      else
!
!   This branch is executed only once
!
         frstco = .false.
!
!   Allocate the initial storage for arrays INDOEI and VALOEI;
!   NDCOEA stores the array dimension
!
         ndcoea = 2
         call alloc (indoei,ndcoea,"INDOEI","IABINT" )
         call alloc (valoei,ndcoea,"VALOEI","IABINT")
!
!   Compute the integral's value
!
         tegral = rinti (ia,ib,0)
!
!   Initialise the integral counter
!
         ncoei = 1
!
!   Store the integral and its value
!
         indoei(1) = index
         valoei(1) = tegral
!
      end if
!
      return
      end subroutine iabint
