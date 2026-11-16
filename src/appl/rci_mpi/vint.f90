!***********************************************************************
!                                                                      *
      subroutine vint (ia,ib,tegral)
!                                                                      *
!   This routine returns Vinti integrals.                              *
!      A list of previously computed integrals is maintained. If the   *
!   required  integral is not already  available it is computed  and   *
!   stored in the list.                                                *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, RALLOC, VINTI.                         *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 06 Jun 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: orbital_index_radix
      use memory_man
      use vinlst_c, first=>frstvi
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use vinti_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ia, ib
      real(double), intent(out) :: tegral
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical :: found
      integer :: index, loc, ju, jl, jm, newsiz, i
!------------------------------------------------
!
!   Compute the composite (packed) index
!
      index = ia*orbital_index_radix+ib
!
      if (.not. first) then
!
!   This branch is executed on all entries except the first
!
         if (index .gt. indtei(nvinti)) then
!
!   The index is greater than the largest stored
!
            found = .false.
            loc = nvinti
!
         elseif (index .lt. indtei(1)) then
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
            ju = nvinti
            jl = 1
!
    1       if (ju-jl .gt. 1) then
!
               jm = (ju+jl)/2
               if (indtei(jm) .gt. index) then
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
               if (index .eq. indtei(ju)) then
!
                  found = .true.
                  loc = ju
!
               elseif (index .eq. indtei(jl)) then
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
            tegral = valtei(loc)
!
         else
!
!   Index not found; compute the integral
!
            tegral = vinti (ia,ib)
!
!   Increment the integral counter
!
            nvinti = nvinti+1
!
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
            if (nvinti .gt. ndvin) then
               newsiz = ndvin+ndvin/2
               call ralloc (indtei,newsiz, 'INDTEI', 'VINT')
               call ralloc (valtei,newsiz, 'VALTEI', 'VINT')
               ndvin = newsiz
            endif
!
            do 2 i = nvinti,loc+2,-1
               indtei(i) = indtei(i-1)
               valtei(i) = valtei(i-1)
    2       continue
!
!   Put the new index and value into storage
!
            indtei(loc+1) = index
            valtei(loc+1) = tegral
!
         endif
!
      else
!
!   This branch is executed only once
!
         first = .false.
!
!   Allocate the initial storage for arrays INDTEI and VALTEI;
!   NDVIN stores the array dimension
!
         ndvin = 2
         call alloc (indtei,ndvin, 'INDTEI', 'VINT')
         call alloc (valtei,ndvin, 'VALTEI', 'VINT')
!
!   Compute the integral's value
!
         tegral = vinti (ia,ib)
!
!   Initialise the integral counter
!
         nvinti = 1
!
!   Store the integral and its value
!
         indtei(1) = index
         valtei(1) = tegral
!
      endif
!
      return
      end subroutine vint
