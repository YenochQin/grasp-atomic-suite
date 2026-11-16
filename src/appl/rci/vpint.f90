!***********************************************************************
!                                                                      *
      subroutine vpint (ia,ib,tegral)
!                                                                      *
!   This  routine computes nuclear vacuum polarization integrals.      *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, RALLOC.                                *
!               [RCI92]: VPINTF.                                       *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
     use vast_kind_param, only: double
     use parameter_def,   only: orbital_index_radix
     use memory_man
     use grid_c
     use vpilst_c, indoei=>indvpi, valoei=>valvpi, first=>frstvp
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
     use vpintf_i
     implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(inout) :: ia, ib
      real(double), intent(out) :: tegral
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical :: found
      integer :: iswap, index, loc, ju, jl, jm, newsiz, i
!------------------------------------------------
!
!   Ensure that the indices are in `canonical' order
!
      if (ia > ib) then
         iswap = ib
         ib = ia
         ia = iswap
      endif
!
!   Compute the composite (packed) index
!
      index = ia*key+ib
!
      if (.not. first) then
!
!   This branch is executed on all entries except the first
!
         if (index > indoei(nvpi)) then
!
!   The index is greater than the largest stored
!
            found = .false.
            loc = nvpi
!
         elseif (index < indoei(1)) then
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
            ju = nvpi
            jl = 1
!
    1       if (ju-jl .gt. 1) then
!
               jm = (ju+jl)/2
               if (indoei(jm) .gt. index) then
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
               if (index .eq. indoei(ju)) then
!
                  found = .true.
                  loc = ju
!
               elseif (index .eq. indoei(jl)) then
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
            tegral = valoei(loc)
!
         else
!
!   Index not found; compute the integral
!
            tegral = vpintf (ia,ib)
!
!   Increment the integral counter
!
            nvpi = nvpi+1
!
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
            if (nvpi .gt. ndvpa) then
               newsiz = ndvpa+ndvpa/2
               call ralloc (indoei,newsiz,'INDOEI', 'VPINT')
               call ralloc (valoei,newsiz,'VALOEI', 'VPINT')
               ndvpa = newsiz
            endif
!
            do 2 i = nvpi,loc+2,-1
               indoei(i) = indoei(i-1)
               valoei(i) = valoei(i-1)
    2       continue
!
!   Put the new index and value into storage
!
            indoei(loc+1) = index
            valoei(loc+1) = tegral
!
         endif
!
      else
!
!   This branch is executed only once
!
         first = .false.
!
!   Allocate the initial storage for arrays INDOEI and VALOEI;
!   NDVPA stores the array dimension
!
         ndvpa = 2
         call alloc (indoei,ndvpa,'INDOEI', 'VPINT')
         call alloc (valoei,ndvpa,'VALOEI', 'VPINT')
!
!   Compute the integral's value
!
         tegral = vpintf (ia,ib)
!
!   Initialise the integral counter
!
         nvpi = 1
!
!   Store the integral and its value
!
         indoei(1) = index
         valoei(1) = tegral
!
      endif
!
      return
      end subroutine vpint
