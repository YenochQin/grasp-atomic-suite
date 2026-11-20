!***********************************************************************
!                                                                      *
      subroutine evcout
      use calculation_file_units, only: summary_output_unit
!                                                                      *
!   Routine for printing eigenvectors.                                 *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use array_memory_management, only: alloc, dalloc
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use orbital_csf_basis,           only: csf_count, orbital_count, csf_subshell_occupations
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V A R I A B L E S
!-----------------------------------------------
      integer :: i, j1, j2, nvex, jp, l, csf_index
      integer, pointer :: iloc(:)
!-----------------------------------------------
!
!   Header
!
      write (summary_output_unit,'(//" Eigenvectors:"/)')
!
!   Determine which eigenvectors are to be printed out
!
      nvex = nvec
      call alloc (iloc,nvec,"ILOC", "EVCOUT")
      do 1 i = 1,nvec
         iloc(i) = ivec(i)
    1 continue
!
!   There are eight columns across the width of a page
!
      jp = 8
!
!   Set up for the first set of columns
!
      j1 = 1
      j2 = jp
!
!   Loop over sets of columns
!
    3 if (j2 > nvex) j2 = nvex
!
!   Write out the column numbers; skip a line
!
      write (summary_output_unit,'(1x,8(i8,7x))') (iloc(l),l = j1,j2)
      write (summary_output_unit,'(1x)')
!
!   Write out the rows
!
      do 4 csf_index = 1,csf_count
!
         write (summary_output_unit,'(1x,1p,8d15.7)') (csf_mixing_coefficients(csf_index+(l-1)*csf_count),l = j1,j2)
!
!   Skip a line after every tenth row
!
         if (mod(csf_index,10) == 0) write (summary_output_unit,'(1x)')
!
    4 continue
!
!   Set up for the next set of columns
!
      write (summary_output_unit,'(1x,1p,8d15.7)')
      if (j2 < nvex) then
         j1 = j1+jp
         j2 = j2+jp
         goto 3
      end if
!
!   Deallocate storage for ILOC
!
      call dalloc (iloc, "ILOC", "EVCOUT")
!
      return
!
!
      end subroutine evcout
