!***********************************************************************
!                                                                      *
      subroutine vpint (first_orbital,second_orbital,radial_integral_value)
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
     use, intrinsic :: iso_fortran_env, only: real64
     use grasp_parameters,   only: orbital_index_radix
     use array_memory_management, only: alloc, ralloc
     use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
         radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
     use vacuum_polarization_integral_table, only: vacuum_polarization_integral_indices => indvpi, &
        & vacuum_polarization_integral_values => valvpi, first_vacuum_polarization_lookup => frstvp, ndvpa, nvpi
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
     use vpintf_i, only: vpintf
     implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(inout) :: first_orbital, second_orbital
      real(real64), intent(out) :: radial_integral_value
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: orbital_label_radix = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical :: found
      integer :: iswap, index, loc, ju, jl, jm, newsiz, i
!------------------------------------------------
!
!   Ensure that the indices are in `canonical' order
!
      if (first_orbital > second_orbital) then
         iswap = second_orbital
         second_orbital = first_orbital
         first_orbital = iswap
      end if
!
!   Compute the composite (packed) index
!
      index = first_orbital*orbital_label_radix+second_orbital
!
      if (.not. first_vacuum_polarization_lookup) then
!
!   This branch is executed on all entries except the first
!
         if (index > vacuum_polarization_integral_indices(nvpi)) then
!
!   The index is greater than the largest stored
!
            found = .false.
            loc = nvpi
!
         else if (index < vacuum_polarization_integral_indices(1)) then
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
    1       if (ju-jl > 1) then
!
               jm = (ju+jl)/2
               if (vacuum_polarization_integral_indices(jm) > index) then
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
               if (index == vacuum_polarization_integral_indices(ju)) then
!
                  found = .true.
                  loc = ju
!
               else if (index == vacuum_polarization_integral_indices(jl)) then
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
            radial_integral_value = vacuum_polarization_integral_values(loc)
!
         else
!
!   Index not found; compute the integral
!
            radial_integral_value = vpintf (first_orbital,second_orbital)
!
!   Increment the integral counter
!
            nvpi = nvpi+1
!
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
            if (nvpi > ndvpa) then
               newsiz = ndvpa+ndvpa/2
               call ralloc (vacuum_polarization_integral_indices,newsiz,"INDOEI", "VPINT")
               call ralloc (vacuum_polarization_integral_values,newsiz,"VALOEI", "VPINT")
               ndvpa = newsiz
            end if
!
            do 2 i = nvpi,loc+2,-1
               vacuum_polarization_integral_indices(i) = vacuum_polarization_integral_indices(i-1)
               vacuum_polarization_integral_values(i) = vacuum_polarization_integral_values(i-1)
    2       continue
!
!   Put the new index and value into storage
!
            vacuum_polarization_integral_indices(loc+1) = index
            vacuum_polarization_integral_values(loc+1) = radial_integral_value
!
         end if
!
      else
!
!   This branch is executed only once
!
         first_vacuum_polarization_lookup = .false.
!
!   Allocate the initial storage for arrays INDOEI and VALOEI;
!   NDVPA stores the array dimension
!
         ndvpa = 2
         call alloc (vacuum_polarization_integral_indices,ndvpa,"INDOEI", "VPINT")
         call alloc (vacuum_polarization_integral_values,ndvpa,"VALOEI", "VPINT")
!
!   Compute the integral's value
!
         radial_integral_value = vpintf (first_orbital,second_orbital)
!
!   Initialise the integral counter
!
         nvpi = 1
!
!   Store the integral and its value
!
         vacuum_polarization_integral_indices(1) = index
         vacuum_polarization_integral_values(1) = radial_integral_value
!
      end if
!
      return
      end subroutine vpint
