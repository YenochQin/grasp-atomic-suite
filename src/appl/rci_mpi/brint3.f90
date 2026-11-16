!***********************************************************************
!                                                                      *
      subroutine brint3 (ia,ib,ic,id,nu,tegral)
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use memory_man
      use bilst_c
      use orb_c,           only: csf_count, orbital_count, csf_subshell_occupations
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use brintf_i, only: brintf
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ia, ib, ic, id, nu
      real(real64), intent(out) :: tegral
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
      if (.not. first(3)) then
!
!   This branch is executed on all entries except the first
!
         if (index > indtp3(ntpi(3))) then
!
!   The index is greater than the largest stored
!
            found = .false.
            loc = ntpi(3)
!
         else if (index < indtp3(1)) then
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
            ju = ntpi(3)
            jl = 1
    1       if (ju-jl > 1) then
               jm = (ju+jl)/2
               if (indtp3(jm) > index) then
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
               if (index == indtp3(ju)) then
!
                  found = .true.
                  loc = ju
!
               else if (index == indtp3(jl)) then
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
            tegral = valtp3(loc)
!
         else
!
!   Index not found; compute the integral
!
            tegral = brintf (3,ia,ib,ic,id,nu)
!
!   Increment the integral counter
!
            ntpi(3) = ntpi(3)+1
!
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
!
            if (ntpi(3) > ndtpa(3)) then
               newsiz = ndtpa(3)+ndtpa(3)/2
               call ralloc (indtp3,newsiz,"INDTP3", "BRINT3")
               call ralloc (valtp3,newsiz,"VALTP3", "BRINT3")
               ndtpa(3)= newsiz
            end if
            do 4 i = ntpi(3),loc+2,-1
               indtp3(i) = indtp3(i-1)
               valtp3(i) = valtp3(i-1)
    4       continue
!
!   Put the new index and value into storage
!
            indtp3(loc+1) = index
            valtp3(loc+1) = tegral
!
         end if
!
      else
!
!   This branch is executed only once per type of integral
!
         first(3) = .false.
!
!   Designate the initial storage for arrays INDTPx and VALTPx;
!   Array NDTPA stores the array dimensions
!
         ndtpa(3) = 10000
!
!   Compute the integral's value
!
         tegral = brintf (3,ia,ib,ic,id,nu)
!
!   Initialise the integral counter
!
         ntpi(3) = 1
!
!   Store the integral and its value
!
         call alloc (indtp3,ndtpa(3),"INDTP3", "BRINT3")
         call alloc (valtp3,ndtpa(3),"VALTP3", "BRINT3")
         indtp3(1) = index
         valtp3(1) = tegral
!
      end if
!
      return
      end subroutine brint3
