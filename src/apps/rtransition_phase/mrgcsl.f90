!***********************************************************************
!                                                                      *
      subroutine mrgcsl(state_name)
!                                                                      *
!   Entry routine for merging two csl lists                            *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use atomic_calculation_state, only:  emn, ioncty, electron_count, nuclear_charge
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ldcsl1_i, only: ldcsl1
      use ldcsl2_i, only: ldcsl2
      use merg12_i, only: merg12
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character  :: state_name(2)*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: initial_core_subshell_count, final_core_subshell_count
!-----------------------------------------------
!

      write (6, *)
      write (6, *) "MRGCSL: Execution begins ..."
!
!   Load the first  .csl  file
!
      call ldcsl1 (initial_core_subshell_count, state_name(1))
!
!   Load the second  .csl  file
!
      call ldcsl2 (final_core_subshell_count, state_name(2))
!
!   Merge the two  .csl  lists, observe that there may be doublets
!   among the CSF's
!
      call merg12 (state_name, initial_core_subshell_count, final_core_subshell_count)

      return
      end subroutine mrgcsl
