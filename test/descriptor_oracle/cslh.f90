module cslh_i
    implicit none (type, external)
    private
    public :: cslh
contains
!***********************************************************************
!
      subroutine cslh(name, ncore, nblkin, idblk, csf_input_unit)
!  A container which calls setcsll to open, read <name> file to get
!     nblock, ncfblk(), idblk(), ncftot.
!  It then calls lib92/lodcsh to get
!     ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()
!  The file pointer points to the first CSL record after this routine.
!
!  Called by rcivu, mcpvu
!
!  Xinghong He 98-06-23
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use array_memory_management, only: alloc, ralloc
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, &
          highest_selected_level_by_block, symmetry_block_count
      use orbital_csf_basis, only: csf_count
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setcsll_i, only: setcsll
      use lodcsh_i, only: lodcsh
      implicit none (type, external)
      integer, intent(out) :: csf_input_unit
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(out) :: ncore
      integer, intent(in) :: nblkin
      character(len=*), intent(inout) :: name
      character(len=8), intent(out) :: idblk(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iqadum
      integer :: ierr
!-----------------------------------------------
!

! node-0 does exactly the same as the serial code does

      call alloc (csf_counts_by_block, nblkin + 1, "NCFBLK", "CSLH")

      call setcsll (csf_input_unit, name, nblkin, symmetry_block_count, csf_counts_by_block, csf_count, idblk)
      call ralloc (csf_counts_by_block, symmetry_block_count + 1, "NCFBLK", "CSLH")
      rewind (csf_input_unit)
      read (csf_input_unit, *)
      !..Load header of <name> file
      call lodcsh (csf_input_unit, ncore)

      return
      end subroutine cslh
end module cslh_i
