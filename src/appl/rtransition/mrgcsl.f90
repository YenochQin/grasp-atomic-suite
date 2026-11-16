!***********************************************************************
!                                                                      *
      subroutine mrgcsl(name)
!                                                                      *
!   Entry routine for merging two csl lists                            *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use debug_c
      use def_c, only:  emn, ioncty, electron_count, nuclear_charge
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ldcsl1_i
      use ldcsl2_i
      use merg12_i
      implicit none
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
