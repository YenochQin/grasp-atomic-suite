!***********************************************************************
!                                                                      *
      subroutine ldcsl1 (ncorer,name)
!                                                                      *
!   Open, check, load data from and close the  .csl  file. This file   *
!   is always attached to stream 21.                                   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, IQ, LENGTH, LODCSL, OPENFL.            *
!                                                                      *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def,   only: max_orbitals
      use memory_man
      use def_c
      use biorb_c
      use orb_c
      use stat_c, only: byte, csf_coupling_quantum_codes, csf_subshell_quantum_codes, double, jcupar, jqsar
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use iq_i, only: iq
      use lodcsl_i, only: lodcsl
      use openfl_i, only: openfl
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: ncorer
      character(len=128) :: name
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       integer :: i, j, ios, ncore
       character(len=15) :: record
!
!   Check the first record of the file; if not as expected, try again
!
      j = index(name," ")
      open (unit = 21,file=name(1:j-1)//".c",form="FORMATTED",          &
     &     status="OLD")

      read (21,"(1A15)",iostat = ios) record
      if ((ios /= 0) .or.                                             &
     &    (record(1:15) /= "Core subshells:")) then
         print *, "Not a Configuration Symmetry List File;"
         close (21)
         stop
      end if
!
!   Load data from the  .csl  file
!
      call lodcsl (ncore)
!
!   Close the  .csl  file
!
      close (21)
!
!   Check if the core should be redefined
!
      ncorer = ncore
!      DO 3 I = NCORE+1,NW
!         IFULLI = NKJ(I)+1
!         DO 2 J = 1,NCF
!            IF (IQ (I,J) .NE. IFULLI) GOTO 4
!    2    CONTINUE
!         CALL CONVRT (NP(I),RECORD,LENTH)
!         PRINT *, 'Subshell '//RECORD(1:LENTH)//NH(I)//' is full'
!         PRINT *, ' in all CSFs; including this'
!         PRINT *, ' subshell in the core;'
!         NCORER = NCORER+1
!    3 CONTINUE
!
!   Copy data to reference arrays
!
    4 nelecr = electron_count
      nelecii = electron_count
      ncfr = csf_count
      nwr = orbital_count
!
      call alloc (iqar, max_orbitals, csf_count, "IQAR", "LDCSL1" )
      call alloc (jqsar,max_orbitals,3,csf_count, "JQSAR", "LDCSL1")
      call alloc (jcupar,max_orbitals, csf_count, "JCUPAR", "LDCSL1")
!
      do 6 i = 1,csf_count
         do 5 j = 1,max_orbitals
            iqar(j,i) = csf_subshell_occupations(j,i)
            jqsar(j,1,i) = csf_subshell_quantum_codes(j,1,i)
            jqsar(j,2,i) = csf_subshell_quantum_codes(j,2,i)
            jqsar(j,3,i) = csf_subshell_quantum_codes(j,3,i)
            jcupar(j,i) = csf_coupling_quantum_codes(j,i)
    5    continue
    6 continue
!
      nwii=orbital_count
      ncfii=csf_count
      do 7 i = 1,orbital_count
         nhii(i) = orbital_symmetry_labels(i)
         npii(i) = orbital_principal_n(i)
         nakii(i) = orbital_kappa(i)
         nklii(i) = orbital_l(i)
         nkjii(i) = orbital_two_j(i)
         npr(i) = orbital_principal_n(i)
         nakr(i) = orbital_kappa(i)
         nhr(i) = orbital_symmetry_labels(i)
    7 continue
!
!   Deallocate storage so that LODCSL can restart correctly
!
      call dalloc (csf_subshell_occupations, "IQA", "LDCSL1")
      call dalloc (csf_subshell_quantum_codes, "JQSA", "LDCSL1")
      call dalloc (csf_coupling_quantum_codes, "JCUPA", "LDCSL1")
!
      close (21)

      return
      end subroutine ldcsl1
