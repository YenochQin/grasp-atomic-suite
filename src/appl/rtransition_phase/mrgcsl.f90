!***********************************************************************
!                                                                      *
      subroutine mrgcsl(name)
!                                                                      *
!   Entry routine for merging two csl lists                            *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use debug_c
      use def_c, only:  emn, ioncty, electron_count, nuclear_charge
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ldcsl1_i
      use ldcsl2_i
      use merg12_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character  :: name(2)*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ncorer, ncore
!-----------------------------------------------
!

      write (6, *)
      write (6, *) 'MRGCSL: Execution begins ...'
!
!   Load the first  .csl  file
!
      call ldcsl1 (ncorer, name(1))
!
!   Load the second  .csl  file
!
      call ldcsl2 (ncore, name(2))
!
!   Merge the two  .csl  lists, observe that there may be doublets
!   among the CSF's
!
      call merg12 (name, ncorer, ncore)

      return
      end subroutine mrgcsl
