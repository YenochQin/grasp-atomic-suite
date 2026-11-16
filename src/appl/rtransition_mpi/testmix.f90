!***********************************************************************
!                                                                      *
      subroutine testmix
!                                                                      *
!   This routine checks the mixing coefficients                        *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use biorb_c
      use def_c
      use eigv_c
      use prnt_c
      use syma_c

      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
!-----------------------------------------------
!
!

      write (*, *) ' ****************'
      write (*, *) ' Entering testmix'
      write (*, *) ' ****************'
      write (*, *)
      write (*, *) 'Initial state'
      write (*, *) 'EVALII', (eavii + evalii(i),i=1,nvecii)
      write (*, *) 'NELECII,NCFII,NWII,NVECMXII', nelecii, ncfii, nwii, nvecmxii
      write (*, *) nvecii
      write (*, *) (ivecii(i),i=1,nvecii)
      write (*, *) (iatjpoii(i),iasparii(i),i=1,nvecii)
      write (*, *) ((evecii(i + (j - 1)*ncfii),i=1,ncfii),j=1,nvecii)

      write (*, *) 'Final state'
      write (*, *) 'EVALFF', (eavff + evalff(i),i=1,nvecff)
      write (*, *) 'NELECFF,NCFFF,NWFF,NVECMXFF', nelecff, ncfff, nwff, nvecmxff
      write (*, *) nvecff
      write (*, *) (ivecff(i),i=1,nvecff)
      write (*, *) (iatjpoff(i),iasparff(i),i=1,nvecff)
      write (*, *) ((evecff(i + (j - 1)*ncfff),i=1,ncfff),j=1,nvecff)
      write (*, *)
      write (*, *) ' ***************'
      write (*, *) ' Leaving testmix'
      write (*, *) ' ***************'

      return
      end subroutine testmix
