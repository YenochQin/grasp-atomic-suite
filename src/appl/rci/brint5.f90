!***********************************************************************
!                                                                      *
      subroutine brint5 (ia,ib,ic,id,nu,tegral)
!                                                                      *
!   Returns integrals for the transverse photon interaction.           *
!      Integrals are stored in ordered lists. If the integral cannot   *
!   be read from a list, it is computed by calling BRINTF.             *
!                                                                      *
!   Observ that it is not possible to use more than 100 orbitals when  *
!   the Breit interaction is included. If 100 orbitals are used then   *
!   NU <= 19.                                                          *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, RALLOC.                                *
!               [RCI92]: BRINTF                                        *
!                                                                      *
!   Written by Farid A Parpia               Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use memory_man
      use bilst_c
      use orb_c,           only: csf_count, orbital_count, csf_subshell_occupations
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use brintf_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ia, ib, ic, id, nu
      real(double), intent(out) :: tegral
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      logical :: found
      integer :: key, index, loc, ju, jl, jm, newsiz, i
!-----------------------------------------------
!
      key = orbital_count + 1
!
!   Compute the integral label
!
      index = (((nu*key+ia)*key+ib)*key+ic)*key+id
!
      if (.not. first(5)) then
!
!   This branch is executed on all entries except the first
!
         if (index .gt. indtp5(ntpi(5))) then
!
!   The index is greater than the largest stored
!
            found = .false.
            loc = ntpi(5)
!
         elseif (index .lt. indtp5(1)) then
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
            ju = ntpi(5)
            jl = 1
    1       if (ju-jl .gt. 1) then
               jm = (ju+jl)/2
               if (indtp5(jm) .gt. index) then
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
               if (index .eq. indtp5(ju)) then
!
                  found = .true.
                  loc = ju
!
               elseif (index .eq. indtp5(jl)) then
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
            tegral = valtp5(loc)
!
         else
!
!   Index not found; compute the integral
!
            tegral = brintf (5,ia,ib,ic,id,nu)
!
!   Increment the integral counter
!
            ntpi(5) = ntpi(5)+1
!
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
!
            if (ntpi(5) .gt. ndtpa(5)) then
               newsiz = ndtpa(5)+ndtpa(5)/2
               call ralloc (indtp5,newsiz,'INDTP5', 'BRINT5')
               call ralloc (valtp5,newsiz,'VALTP5', 'BRINT5')
               ndtpa(5)= newsiz
            endif
            do 4 i = ntpi(5),loc+2,-1
               indtp5(i) = indtp5(i-1)
               valtp5(i) = valtp5(i-1)
    4       continue
!
!   Put the new index and value into storage
!
            indtp5(loc+1) = index
            valtp5(loc+1) = tegral
!
         endif
!
      else
!
!   This branch is executed only once per type of integral
!
         first(5) = .false.
!
!   Designate the initial storage for arrays INDTPx and VALTPx;
!   Array NDTPA stores the array dimensions
!
         ndtpa(5) = 10000
!
!   Compute the integral's value
!
         tegral = brintf (5,ia,ib,ic,id,nu)
!
!   Initialise the integral counter
!
         ntpi(5) = 1
!
!   Store the integral and its value
!
         call alloc (indtp5,ndtpa(5),'INDTP5', 'BRINT5')
         call alloc (valtp5,ndtpa(5),'VALTP5', 'BRINT5')
         indtp5(1) = index
         valtp5(1) = tegral
!
      endif
!
      return
      end
