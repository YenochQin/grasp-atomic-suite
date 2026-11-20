!***********************************************************************
!                                                                      *
      subroutine vint (first_orbital,second_orbital,radial_integral_value)
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
!   Converted by CFF
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: orbital_index_radix
      use array_memory_management, only: alloc, ralloc
      use vinti_integral_table, only: first_vinti_lookup => frstvi, indtei, ndvin, nvinti, valtei
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use vinti_i, only: vinti
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: first_orbital, second_orbital
      real(real64), intent(out) :: radial_integral_value
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: orbital_label_radix = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical :: found
      integer :: index, loc, ju, jl, jm, newsiz, i
!------------------------------------------------
!
!   Compute the composite (packed) index
!
      index = first_orbital*orbital_index_radix+second_orbital
!
      if (.not. first_vinti_lookup) then
!
!   This branch is executed on all entries except the first
!
         if (index > indtei(nvinti)) then
!
!   The index is greater than the largest stored
!
            found = .false.
            loc = nvinti
!
         else if (index < indtei(1)) then
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
    1       if (ju-jl > 1) then
!
               jm = (ju+jl)/2
               if (indtei(jm) > index) then
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
               if (index == indtei(ju)) then
!
                  found = .true.
                  loc = ju
!
               else if (index == indtei(jl)) then
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
            radial_integral_value = valtei(loc)
!
         else
!
!   Index not found; compute the integral
!
            radial_integral_value = vinti (first_orbital,second_orbital)
!
!   Increment the integral counter
!
            nvinti = nvinti+1
!
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
            if (nvinti > ndvin) then
               newsiz = ndvin+ndvin/2
               call ralloc (indtei,newsiz, "INDTEI", "VINT")
               call ralloc (valtei,newsiz, "VALTEI", "VINT")
               ndvin = newsiz
            end if
!
            do 2 i = nvinti,loc+2,-1
               indtei(i) = indtei(i-1)
               valtei(i) = valtei(i-1)
    2       continue
!
!   Put the new index and value into storage
!
            indtei(loc+1) = index
            valtei(loc+1) = radial_integral_value
!
         end if
!
      else
!
!   This branch is executed only once
!
         first_vinti_lookup = .false.
!
!   Allocate the initial storage for arrays INDTEI and VALTEI;
!   NDVIN stores the array dimension
!
         ndvin = 2
         call alloc (indtei,ndvin, "INDTEI", "VINT")
         call alloc (valtei,ndvin, "VALTEI", "VINT")
!
!   Compute the integral's value
!
         radial_integral_value = vinti (first_orbital,second_orbital)
!
!   Initialise the integral counter
!
         nvinti = 1
!
!   Store the integral and its value
!
         indtei(1) = index
         valtei(1) = radial_integral_value
!
      end if
!
      return
      end subroutine vint
