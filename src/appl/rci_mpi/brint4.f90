!***********************************************************************
!                                                                      *
      subroutine brint4 (ia,ib,ic,id,nu,tegral)
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
      use brintf_i
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
      if (.not. first(4)) then
!
!   This branch is executed on all entries except the first
!
         if (index .gt. indtp4(ntpi(4))) then
!
!   The index is greater than the largest stored
!
            found = .false.
            loc = ntpi(4)
!
         elseif (index .lt. indtp4(1)) then
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
            ju = ntpi(4)
            jl = 1
    1       if (ju-jl .gt. 1) then
               jm = (ju+jl)/2
               if (indtp4(jm) .gt. index) then
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
               if (index .eq. indtp4(ju)) then
!
                  found = .true.
                  loc = ju
!
               elseif (index .eq. indtp4(jl)) then
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
            tegral = valtp4(loc)
!
         else
!
!   Index not found; compute the integral
!
            tegral = brintf (4,ia,ib,ic,id,nu)
!
!   Increment the integral counter
!
            ntpi(4) = ntpi(4)+1
!
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
!
            if (ntpi(4) .gt. ndtpa(4)) then
               newsiz = ndtpa(4)+ndtpa(4)/2
               call ralloc (indtp4,newsiz,'INDTP4', 'BRINT4')
               call ralloc (valtp4,newsiz,'VALTP4', 'BRINT4')
               ndtpa(4)= newsiz
            endif
            do 4 i = ntpi(4),loc+2,-1
               indtp4(i) = indtp4(i-1)
               valtp4(i) = valtp4(i-1)
    4       continue
!
!   Put the new index and value into storage
!
            indtp4(loc+1) = index
            valtp4(loc+1) = tegral
!
         endif
!
      else
!
!   This branch is executed only once per type of integral
!
         first(4) = .false.
!
!   Designate the initial storage for arrays INDTPx and VALTPx;
!   Array NDTPA stores the array dimensions
!
         ndtpa(4) = 10000
!
!   Compute the integral's value
!
         tegral = brintf (4,ia,ib,ic,id,nu)
!
!   Initialise the integral counter
!
         ntpi(4) = 1
!
!   Store the integral and its value
!
         call alloc (indtp4,ndtpa(4),'INDTP4', 'BRINT4')
         call alloc (valtp4,ndtpa(4),'VALTP4', 'BRINT4')
         indtp4(1) = index
         valtp4(1) = tegral
!
      endif
!
      return
      end
