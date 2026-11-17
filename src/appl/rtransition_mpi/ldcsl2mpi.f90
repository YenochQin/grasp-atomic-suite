!***********************************************************************
!                                                                      *
      subroutine ldcsl2(ncore, name)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!                                                                      *
!   Open, check, load data from and close the  .csl  file. This file   *
!   is always attached to stream 21.                                   *
!                                                                      *
!   Call(s) to: [LIB92]: IQ, LENGTH, LODCSL, OPENFL.                   *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use def_c, only: double, electron_count, nelecff
      use orb_c, only: csf_count, double, max_orbitals, orbital_count, orbital_kappa, orbital_l, orbital_principal_n, &
          orbital_symmetry_labels, orbital_two_j
      use biorb_c, only: double, max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
          nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use csf_basis_mpi, only: load_csf_basis_mpi
      use blk_c, only: nblock, ncfblk, nblocki, nblockf, ncfi, ncff
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(out) :: ncore
      character , intent(in) :: name*128
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, ios, ncorer, i
      character :: record*15
!-----------------------------------------------
!
!
!     Common relevant for the radial wave functions of the final state.
!
!
!   The  .csl  file is FORMATTED; it must exist
!
      call load_csf_basis_mpi(trim(name),ncorer)
      ios=0
      if (ios /= 0) error stop "Transition descriptor load failed"
      if (nblock > size(ncff)) error stop "Transition block capacity exceeded"
      nblockf=nblock
      ncff=0
      ncff(:nblock)=ncfblk(:nblock)
!
!   Check if the core should be redefined
!
      ncore = ncorer
!      DO 3 I = NCORER+1,NW
!         IFULLI = NKJ(I)+1
!         DO 2 J = 1,NCF
!            IF (IQ (I,J) .NE. IFULLI) GOTO 4
!    2    CONTINUE
!         CALL CONVRT (NP(I),RECORD,LENTH)
!         PRINT *, 'Subshell '//RECORD(1:LENTH)//NH(I)//' is full'
!         PRINT *, ' in all CSFs; including this'
!         PRINT *, ' subshell in the core;'
!         NCORE = NCORE+1
!    3 CONTINUE
!
      nelecff = electron_count
      nwff = orbital_count
      ncfff = csf_count
      nhff(:orbital_count) = orbital_symmetry_labels(:orbital_count)
      npff(:orbital_count) = orbital_principal_n(:orbital_count)
      nakff(:orbital_count) = orbital_kappa(:orbital_count)
      nklff(:orbital_count) = orbital_l(:orbital_count)
      nkjff(:orbital_count) = orbital_two_j(:orbital_count)
      return
      end subroutine ldcsl2
