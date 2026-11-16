!***********************************************************************
!                                                                      *
      subroutine radfile(name)
!                                                                      *
!  This subroutine outputs the transformed radial orbitals             *
!                                                                      *
!  Written by Per Jonsson                                              *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use biorb_c
      use grid_c
      use orb_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character, intent(in) :: name(2)*128
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, k, i
!-----------------------------------------------
!
      j = index(name(1),' ')
      open(unit=30, file=name(1)(1:j-1)//'.bw', form='UNFORMATTED', status=&
         'UNKNOWN', position='asis')

      write (30) 'G92RWF'
      write (*, *) 'NWII', nwii
      do k = 1, nwii
         write (30) npii(k), nakii(k), eii(k), mfii(k)
         write (30) pzii(k), (pfii(i,k),i=1,mfii(k)), (qfii(i,k),i=1,mfii(k))
         write (30) (radial_grid(i),i=1,mfii(k))
      end do

      close(30)

      j = index(name(2),' ')
      open(unit=30, file=name(2)(1:j-1)//'.bw', form='UNFORMATTED', status=&
         'UNKNOWN', position='asis')

      write (30) 'G92RWF'
      write (*, *) 'NWFF', nwff
      do k = 1, nwff
         write (30) npff(k), nakff(k), eff(k), mfff(k)
         write (30) pzff(k),(pfff(i,k),i=1,mfff(k)),(qfff(i,k),i=1,mfff(k))
         write (30) (radial_grid(i),i=1,mfff(k))
      end do

      close(30)

      return
      end subroutine radfile
