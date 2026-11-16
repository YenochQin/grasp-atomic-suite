!***********************************************************************
!                                                                      *
      subroutine readmix(name, inpci, init)
!                                                                      *
!   Open and read the mixing coefficent files                          *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:38   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use memory_man
      use biorb_c
      use def_c
      use orb_c
      use eigv_c
      use prnt_c
      use syma_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: inpci
      integer , intent(in) :: init
      character  :: name(2)*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ib, iatjp, iaspa, i, j
      character :: g92mix*6
!-----------------------------------------------
!
      if (init == 1) then
!
!   Read the initial state mixing file
!
         read (68) ib, ncfii, nvecii, iatjp, iaspa
         call alloc (evalii, nvecii, 'EVALII', 'READMIX')
         call alloc (evecii, ncfii*nvecii, 'EVECII', 'READMIX')
         call alloc (ivecii, nvecii,'IVECII', 'READMIX' )
         call alloc (iatjpoii, nvecii, 'IATJPOII', 'READMIX')
         call alloc (iasparii, nvecii, 'ISPARII', 'READMIX')
         read (68) (ivecii(i),i=1,nvecii)
         iatjpoii(:nvecii) = iatjp
         iasparii(:nvecii) = iaspa
         read (68) eavii, (evalii(i),i=1,nvecii)

         read (68) ((evecii(i + (j - 1)*ncfii),i=1,ncfii),j=1,nvecii)
!
!     CLOSE(68)
!
      else
!
!   Read the final state mixing file
!

         read (78) ib, ncfff, nvecff, iatjp, iaspa
         call alloc (evalff, nvecff, 'EVALFF', 'READMIX')
         call alloc (evecff, ncfff*nvecff, 'EVECFF', 'READMIX')
         call alloc (ivecff, nvecff,'IVECFF', 'READMIX' )
         call alloc (iatjpoff, nvecff, 'IATJPOFF', 'READMIX')
         call alloc (iasparff, nvecff, 'ISPARFF', 'READMIX')
         read (78) (ivecff(i),i=1,nvecff)
         iatjpoff(:nvecff) = iatjp
         iasparff(:nvecff) = iaspa
         read (78) eavff, (evalff(i),i=1,nvecff)

         read (78) ((evecff(i + (j - 1)*ncfff),i=1,ncfff),j=1,nvecff)
!
!   Close the initial state mixing  file
!
!      CLOSE(78)
      endif

      return
      end subroutine readmix
