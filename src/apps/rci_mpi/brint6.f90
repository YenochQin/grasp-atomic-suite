!***********************************************************************
!                                                                      *
      subroutine brint6 (orbital_a,orbital_b,orbital_c,orbital_d,radial_integral_order,radial_integral_value)
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
      use array_memory_management, only: alloc, ralloc
      use breit_integral_tables, only: ndtpa, ntpi, indtp1, indtp2, indtp3, indtp4, indtp5, indtp6, &
          valtp1, valtp2, valtp3, valtp4, valtp5, valtp6, first
      use orbital_csf_basis,           only: csf_count, orbital_count, csf_subshell_occupations
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use brintf_i, only: brintf
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: orbital_a, orbital_b, orbital_c, orbital_d, radial_integral_order
      real(real64), intent(out) :: radial_integral_value
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
      index = (((radial_integral_order*key+orbital_a)*key+orbital_b)*key+orbital_c)*key+orbital_d
!
      if (.not. first(6)) then
!
!   This branch is executed on all entries except the first
!
         if (index > indtp6(ntpi(6))) then
!
!   The index is greater than the largest stored
!
            found = .false.
            loc = ntpi(6)
!
         else if (index < indtp6(1)) then
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
            ju = ntpi(6)
            jl = 1
    1       if (ju-jl > 1) then
               jm = (ju+jl)/2
               if (indtp6(jm) > index) then
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
               if (index == indtp6(ju)) then
!
                  found = .true.
                  loc = ju
!
               else if (index == indtp6(jl)) then
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
            radial_integral_value = valtp6(loc)
!
         else
!
!   Index not found; compute the integral
!
            radial_integral_value = brintf (6,orbital_a,orbital_b,orbital_c,orbital_d,radial_integral_order)
!
!   Increment the integral counter
!
            ntpi(6) = ntpi(6)+1
!
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
!
            if (ntpi(6) > ndtpa(6)) then
               newsiz = ndtpa(6)+ndtpa(6)/2
               call ralloc (indtp6,newsiz,"INDTP6", "BRINT6")
               call ralloc (valtp6,newsiz,"VALTP6", "BRINT6")
               ndtpa(6)= newsiz
            end if
            do 4 i = ntpi(6),loc+2,-1
               indtp6(i) = indtp6(i-1)
               valtp6(i) = valtp6(i-1)
    4       continue
!
!   Put the new index and value into storage
!
            indtp6(loc+1) = index
            valtp6(loc+1) = radial_integral_value
!
         end if
!
      else
!
!   This branch is executed only once per type of integral
!
         first(6) = .false.
!
!   Designate the initial storage for arrays INDTPx and VALTPx;
!   Array NDTPA stores the array dimensions
!
         ndtpa(6) = 10000
!
!   Compute the integral's value
!
         radial_integral_value = brintf (6,orbital_a,orbital_b,orbital_c,orbital_d,radial_integral_order)
!
!   Initialise the integral counter
!
         ntpi(6) = 1
!
!   Store the integral and its value
!
         call alloc (indtp6,ndtpa(6),"INDTP6", "BRINT6")
         call alloc (valtp6,ndtpa(6),"VALTP6", "BRINT6")

         indtp6(1) = index
         valtp6(1) = radial_integral_value
!
      end if
!
      return
      end subroutine brint6
