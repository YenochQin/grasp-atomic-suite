!***********************************************************************
!                                                                      *
      subroutine getmixblock(name, nci)
!                                                                      *
!   Reads mixing coefficient file from block-structured format         *
!                                                                      *
!   Note:                                                              *
!   eav is not compatible with the non-block version if some blocks    *
!   were not diagonalized                                              *
!                                                                      *
!   This is a modified version of cvtmix.f                             *
!                                                                      *
!   Written by Per Jonsson,                           September 2003   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!   Modified by C. Cychen                                       2021   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64, int64
      use memory_man
      use def_c
      use eigv_c
      use orb_c
      use prnt_c, only: byte, double, ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use syma_c, only: double, iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
      use blk_c,            only: nevinblk, ncfinblk, two_j
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i, only: openfl
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nci
      character(len=24), intent(in) :: name
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, ierr, ios, ncftot, nvectot, nvecsiz, nblock, i, nvecpat, &
         ncfpat, neavsum, jb, nb, ncfblk, nevblk, iatjp, iaspa, j
      real(real64) :: eavsum
      integer(int64) :: nvecsizpat, ncftot_i8, nvecsiz_i8
      character(len=3)   :: status
      character(len=6)   :: g92mix
      character(len=11)  :: form
      character(len=256) :: filnam
!-----------------------------------------------
!
!   The  .mix  file is UNFORMATTED; it must exist
!
      k = index(name," ")
      if (nci == 0) then
         filnam = name(1:k-1)//".cm"
      else
         filnam = name(1:k-1)//".m"
      end if
      form = "UNFORMATTED"
      status = "OLD"
!
      call openfl (25, filnam, form, status, ierr)
      if (ierr == 1) then
         write (stderr_unit, *) "Error when opening", filnam
         stop
      end if
!
!   Check the header of the file; if not as expected, try again
!
      read (25, iostat=ios) g92mix
      if (ios/=0 .or. g92mix/="G92MIX") then
         write (stderr_unit, *) "Not a GRASP92 MIXing Coefficients File;"
         close(25)
         stop
      end if

      read (25) electron_count, ncftot, orbital_count, nvectot, nvecsiz, nblock
      write (*, *) "   nelec  = ", electron_count
      write (*, *) "   ncftot = ", ncftot
      write (*, *) "   nw     = ", orbital_count
      write (*, *) "   nblock = ", nblock
      write (*, *)
      ncftot_i8 = ncftot

!***********************************************************************
! Allocate memory for old format data
!***********************************************************************

!GG      CALL ALLOC (EVAL, NVECTOT, 'EVAL', 'GETMIXBLOCK')
!GG      CALL ALLOC (EVEC, NCFTOT*NVECTOT, 'EVEC', 'GETMIXBLOCK')
!GG      CALL ALLOC (IVEC, NVECTOT, 'IVEC', 'GETMIXBLOCK')
!GG      CALL ALLOC (IATJPO, NVECTOT, 'IATJPO', 'GETMIXBLOCK')
!GG      CALL ALLOC (IASPAR, NVECTOT, 'IASPAR', 'GETMIXBLOCK')
      allocate (shifted_eigenvalues(nvectot))
      allocate (csf_mixing_coefficients(ncftot_i8*nvectot))
      allocate (ivec(nvectot))
      allocate (iatjpo(nvectot))
      allocate (iaspar(nvectot))

!***********************************************************************
! Initialize mixing coefficients to zero; others are fine
!***********************************************************************
      csf_mixing_coefficients(:nvectot*ncftot_i8) = 0.0_real64

!***********************************************************************
! Initialize counters and sum registers
!
!    nvecpat:    total number of eigenstates of the previous blocks
!    ncfpat:     total number of CSF of the previous blocks
!    nvecsizpat: vector size of the previous blocks
!    eavsum:     sum of diagonal elements of the previous blocks where
!                at least one eigenstate is calculated
!    neavsum:    total number CSF contributing to eavsum
!***********************************************************************

      nvecpat = 0
      ncfpat = 0
      nvecsizpat = 0
      neavsum = 0
      eavsum = 0.0_real64

      write (*, *) "  block     ncf     nev    2j+1  parity"
      do jb = 1, nblock

         read (25) nb, ncfblk, nevblk, iatjp, iaspa
         write (*, "(5I8)") nb, ncfblk, nevblk, iatjp, iaspa
         nevinblk(jb) = nevblk
         ncfinblk(jb) = ncfblk
         two_j(jb) = iatjp - 1
         if (jb /= nb) stop "jb .NE. nb"

         if (nevblk > 0) then

            read (25) (ivec(nvecpat + i),i=1,nevblk)
               ! ivec(i)   = ivec(i) + ncfpat ! serial # of the state
            iatjpo(nvecpat+1:nevblk+nvecpat) = iatjp
            iaspar(nvecpat+1:nevblk+nvecpat) = iaspa

            read (25) hamiltonian_energy_shift, (shifted_eigenvalues(nvecpat+i),i=1,nevblk)

!           ...Construct the true energy by adding up the average
            shifted_eigenvalues(nvecpat+1:nevblk+nvecpat) = shifted_eigenvalues(nvecpat+1:nevblk+nvecpat) + &
               hamiltonian_energy_shift
!           ...For overal (all blocks) average energy
            eavsum = eavsum + hamiltonian_energy_shift*ncfblk
            neavsum = neavsum + ncfblk

            read (25) ((csf_mixing_coefficients(nvecsizpat+ncfpat+i+(j-1)*ncftot_i8),i=1,ncfblk),j=1,&
               nevblk)
         end if
!
         nvecpat = nvecpat + nevblk
         ncfpat = ncfpat + ncfblk
         nvecsizpat = nvecsizpat + nevblk*ncftot_i8
!
      end do

!     ...Here eav is the average energy of the blocks where at least
!        one eigenstate is calculated. It is not the averge of the
!        total Hamiltonian.

      hamiltonian_energy_shift = eavsum/neavsum

      if (ncftot /= neavsum) write (6, *) &
         "Not all blocks are diagonalized --- Average E ", "not correct"

!     ...Substrct the overal average energy
      shifted_eigenvalues(:nvectot) = shifted_eigenvalues(:nvectot) - hamiltonian_energy_shift
!
      close(25)
!
      csf_count = ncftot
      nvec = nvectot
!
      return
      end subroutine getmixblock
