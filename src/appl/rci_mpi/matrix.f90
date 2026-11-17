!***********************************************************************
!                                                                      *
      subroutine matrix (ncore, j2max)
!                                                                      *
!   This SUBROUTINE calls routines to  form the  Hamiltonian  matrix   *
!   and to diagonalise it.   The total angular momenta and parity of   *
!   the ASF are found;  the eigenvectors are chosen so that the sign   *
!   of the largest element is positive.                                *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, HMOUT, IQ              *
!                        WGHTD5.                                       *
!               [RCI92]: MANEIG, QED, SETHAM.                          *
!                                                                      *
!                                         Last revision: 28 Dec 1992   *
!   Modified by Xinghong He               Last revision: 31 Dec 1997   *
!   Block version Xinghong He             Last revision: 12 Jun 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64, int64
      use parameter_def,   only: max_orbitals
      use memory_man, only: alloc, dalloc, byte, double, long
      use decide_c, only: double, lfordr, lnms, lse, lsms, ltrans, lvp
      use debug_c, only: cutoff, double, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use def_c, only: double, hartree_to_ev, hartree_to_wavenumber_cm, iccmin
      use eigv_c, only: double, hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use hblock_c, only: csf_counts_by_block, double, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use hmat_c, only: double, long, stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
          matrix_partner_indices, ntpitmp, nelmnttmp, ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp, ncoretmp, &
          nvpitmp, nkeitmp, nvintitmp, ncftmp, cutofftmp
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
      use jlabl_c, only: labj => jlbr, labp => jlbp
      use orb_c, only: csf_count, orbital_count, csf_subshell_occupations
      use prnt_c, only: byte, double, ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use stat_c, only: byte, csf_coupling_quantum_codes, csf_subshell_quantum_codes, double, jcupar, jqsar
      use wave_c, only: double, max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
      use where_c, only: double, imcdf
      use blim_c, only: coeffcut1, coeffcut2, double, iprerun, ncsfpre
      use eigvec1_c, only: double, evec1
      use iccu_c, only: iccut, iccutblk
      use cteilsrk_c, only: double, indteirk, valteirk, rk_local_count
      use coeils_c, only: double, frstco, indoei, ncoei, ndcoea, valoei
      use bilst_c, only: double, ndtpa, ntpi, indtp1, indtp2, indtp3, indtp4, indtp5, indtp6, &
          valtp1, valtp2, valtp3, valtp4, valtp5, valtp6, first
      use vpilst_c, only: double, frstvp, indvpi, ndvpa, nvpi, valvpi
      use keilst_c, only: double, frstki, indkei, ndkea, nkei, valkei
      use vinlst_c, only: double, frstvi, indtei, ndvin, nvinti, valtei
      use fposition_c, only: nposition
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use auxblk_i, only: auxblk
      use lodcslmpi_i, only: lodcslmpi
      use qed_slfen_i, only: qed_slfen
      use genmat_i, only: genmat
      use genmat2_i, only: genmat2
      use hmout_i, only: hmout
      use iq_i, only: iq
      use maneig_i, only: maneig
      use engout_i, only: engout
      use wghtd5_i, only: wghtd5
      use qed_i, only: qed
      use mpi_c, only: double, mpi_process_count, mpi_rank
      use mpiu_i, only: stopmpi
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ncore
      integer,  intent(in):: j2max
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64) :: slfint(max_orbitals)
      character(len=8) :: cnum
      real(real64) :: atwinv, elsto, eau, ecm, eev, elemnt
      real(real64), pointer :: slf_en(:), ucf(:), etot(:)
      integer(int64) :: nelmnt_a
      integer :: iiatjpo, iiaspar
      integer :: i, j, irestart, ncminpas, jblock, ic, ip, mode
!-----------------------------------------------
!     ...Common to all blocks - place here to save CPU time
      call auxblk (j2max, atwinv)

!***************************************************************
!      Loop over blocks
!***************************************************************
      ncminpas = 0

      do 100 jblock = 1, symmetry_block_count
         csf_count    =   csf_counts_by_block(jblock)
         nvec   =   eigenvector_counts_by_block(jblock)
         nvecmx = ncmaxblk(jblock)
         iccut(1)  = iccutblk(jblock)
         !.. Determine position of the previous block in the .res file

         nposition = 7 + orbital_count + orbital_count  ! File position of the previous block
                                  ! in the .res file
         do i = 1, jblock - 1
            j = (csf_counts_by_block(i) - mpi_rank - 1 + mpi_process_count) / mpi_process_count
            if (csf_counts_by_block(i) < mpi_process_count) j = csf_counts_by_block(i) / (mpi_rank+1)
            nposition = nposition + j + 1
         end do

         !.. SETHAM does not handle this extrem case
         if (mpi_process_count > csf_count)                                          &
                         call stopmpi ("matrix: too many nodes", mpi_rank)

!        ...Obtain ivec() from iccmin()
         if (nvec > 0) then
            call alloc (ivec, nvec, "IVEC", "MATRIX")
            do i = 1, nvec
               ivec(i) = iccmin(i+ncminpas)
            end do
            ncminpas = ncminpas + nvec
         end if

!        ...These 3 were allocated in lodcsh2 and deallocated at the end
!        ... of this routine and in the setham. In this block version,
!        ... both allocation and deallocation are placed here. See the
!        ... following goto 80 for reason.
         call alloc (csf_subshell_occupations, max_orbitals, csf_count, "IQA", "MATRIX")
         call alloc (csf_subshell_quantum_codes, max_orbitals, 3, csf_count, "JQSA", "MATRIX")
         call alloc (csf_coupling_quantum_codes, max_orbitals, csf_count, "JCUPA", "MATRIX")
         call alloc (slf_en,csf_count,"SLF_EN", "MATRIX")
         call alloc (ucf, orbital_count,"UCF", "MATRIX")
          do ic=1,csf_count
               slf_en(ic) = 0.0
          end do

!      ...Load CSF list of the current block
        call lodcslmpi (21, ncore, jblock)
! zou

         if (lse) then
            if (mpi_rank == 0) then
              print *, "Entering QED ..."

              write (24,*)
              write (24,*) " Self Energy Corrections: "
              write (24,*)
            end if
            call qed_slfen (slfint)
            do ic = 1, csf_count
               elemnt = 0.0_real64
               do i = 1,orbital_count
                  elemnt = elemnt+iq (i,ic)*slfint(i)
               end do
               slf_en(ic) = elemnt
            end do

            if (mpi_rank == 0) then
               write (24,*)
               write (24,*) "Self-energy corrections estimated"       &
                       //" --- these will influence the data"
               write (24,*) " in the RCI92 MIXing coefficients File."
               end if
            end if

! zou
         if (nvec <= 0) then
            if (mpi_rank == 0) write (25) jblock, csf_count, nvec, 999, 999
!           ...Generate H matrix of the current block and write to disk
!           ...eav, nelmnt are also obtained from genmat
            call genmat (atwinv, jblock, mpi_rank, mpi_process_count, elsto, irestart,&
           slf_en)
             ! This call is optional, meaning the matrix of this block
               ! does not need to be computed. But don't comment it out
               ! since other programs may assume thet existence of them.
            call genmat2 (irestart, nelmnt_a, elsto)
            goto 80 ! need to clear memory
         end if
!        ------------------------
         call genmat (atwinv, jblock, mpi_rank, mpi_process_count, elsto, irestart,   &
           slf_en)
         call genmat2 (irestart, nelmnt_a, elsto)
!
!   Allocate and deallocate memory for the mixing coefficients
!   from the prerun
!
      if (iprerun==1) call alloc (evec1,csf_count*nvec, "EVEC1", "MATRIX")
      if (iprerun==2) call dalloc (evec1, "EVEC1", "MATRIX")

! Since maneig needs both nelmnt and nelmnt_a, # of non-zero
! matrix elements computed by current node and the total #
! from all nodes, a new var nelmnt_a is created rathan than
! using the one in the common block.

!        write (*,*) 'Matrix61,myid=',myid
        call maneig (iiatjpo, iiaspar, nelmnt_a)
!GG        CALL MANEIG (iiatjpo, iiaspar, nelmnt_a, imethod)
!GG         CALL MANEIG (iiatjpo, iiaspar)
!
!   Write out eigenvalues (ENGOUT), dominant components of the
!   eigenvectors (WGHTD5) to stream 24; write out ASF symmetries,
!   eigenvalues  and eigenvectors to RCI92 mixing coefficients file.
!   EAV and ELSTO are added back to energy here
!
        if (mpi_rank == 0) then
! ELSTO has never been in Hamiltonian matrix, yet it was
! added to EAV which was later substracted from H. Thus at
! this point, EAV is correct (it has ELSTO added), EVAL()
! need ELSTO and the correct EAV.
           if (csf_count > 1) then
              do i = 1, nvec
                 shifted_eigenvalues(i) = shifted_eigenvalues(i) + elsto
              end do
           end if

           call engout (hamiltonian_energy_shift,shifted_eigenvalues,iiatjpo,iiaspar,ivec,nvec,3)
           call wghtd5 (iiatjpo, iiaspar)

!          ...Write ASF symmetries, eigenvalues, and eigenvectors to RCI92
!          ...MIXing coefficients File; close the file; print a report
           write (25) jblock, csf_count, nvec, iiatjpo, iiaspar
           write (25) (ivec(i), i = 1, nvec)
           write (25) hamiltonian_energy_shift,(shifted_eigenvalues(i),i = 1,nvec)
           write (25) ((csf_mixing_coefficients(i+(j-1)*csf_count),i = 1,csf_count),j = 1,nvec)

           print *, "RCI90 MIXing coefficients File generated."
        end if
!
!   Save the mixing coefficients from the prerun
!
      if (iprerun == 1) then
         do j = 1, nvec
            do i = 1, csf_count
               evec1(i+(j-1)*csf_count) = csf_mixing_coefficients(i+(j-1)*csf_count)
            end do
         end do
      end if
!
!   Estimate diagonal self-energy corrections; the stored
!   eigenvalues and eigenvectors are not modified by these
!   estimates
!
      if (.not.lse) then
         if (mpi_rank == 0)            print *, "Entering QED ..."
            call alloc (etot,nvec,"ETOT", "MATRIX")
         if (mpi_rank == 0) then
            write (24,*)
            write (24,*) " Self Energy Corrections: "
            write (24,*)
            write (24,301)
            write (24,*)
         end if
  301 format (" Level  J Parity",7x,"Hartrees",14x,"Kaysers",          &
               16x,"eV" )
  302 format (1i3,2x,2a4,1p,3d22.14)
         do j = 1, nvec
            call qed (j,slfint,ucf)
            elemnt = 0.0_real64
            ic = ivec(j)
            do i = 1,orbital_count
               elemnt = elemnt+ucf(i)*slfint(i)
!              ELEMNT = ELEMNT+IQ (I,IC)*SLFINT(I)
            end do
            etot(j) = shifted_eigenvalues(j)+elemnt
!
            eau = elemnt
            ecm = eau*hartree_to_wavenumber_cm
            eev = eau*hartree_to_ev
            ip = (iiaspar+3)/2
            if (mpi_rank == 0)                                           &
               write (24,302) j,labj(iiatjpo),labp(ip),eau,ecm,eev
!
         end do
         if (mpi_rank == 0) then
            write (24,*)
            write (24,*) "Self-energy corrections estimated"           &
                    //" --- these do not influence the data"
            write (24,*) " in the RCI92 MIXing coefficients File."
            call engout (hamiltonian_energy_shift,etot,iiatjpo,iiaspar,ivec,nvec,mode)
! zou       CALL ENGOUT (EAV+elsto,ETOT,IiATJPO,iIASPAR,IVEC,NVEC,MODE)
         end if
         call dalloc (etot, "ETOT", "MATRIX")
      end if

!        ...Locals
         call dalloc (ivec, "IVEC", "MATRIX")
!        ...Allocated in maneig
         call dalloc (shifted_eigenvalues, "EVAL", "MATRIX")
         call dalloc (csf_mixing_coefficients, "EVEC", "MATRIX")

   80    continue

!        ...Locals
         call dalloc (csf_subshell_occupations, "IQA", "MATRIX")
         call dalloc (csf_subshell_quantum_codes, "JQSA", "MATRIX")
         call dalloc (csf_coupling_quantum_codes, "JCUPA", "MATRIX")
         call dalloc (slf_en, "SLF_EN", "MATRIX")
         call dalloc (ucf, "UCF", "MATRIX")

  100 continue
!
!   Close the restart files; nothing will be added to them now
!
      close (imcdf)
!     ...Clean up
      call dalloc (csf_counts_by_block, "NCFBLK", "MATRIX")
      call dalloc (eigenvector_counts_by_block, "NEVBLK", "MATRIX")
      call dalloc (ncmaxblk, "NCMAXBLK", "MATRIX")
      call dalloc (iccutblk, "ICUTTBLK", "MATRIX")
      call dalloc (iccmin, "ICCMIN", "MATRIX") ! allocated in items as pnccmin

      call dalloc (valteirk, "VALTEIRK", "MATRIX") ! allocated in genintrk
      call dalloc (indteirk, "INDTEIRK", "MATRIX") ! allocated in genintrk
!
!   Deallocate storage for the integral lists from the
!   Dirac-Coulomb operator; the storage was allocated
!   in IABINT and RKINTC
!
      if (ncoei > 0) then
         call dalloc (indoei, "INDOEI", "MATRIX")
         call dalloc (valoei, "VALOEI", "MATRIX")
      end if
!
!   Deallocate storage for the integral lists from the
!   transverse photon interaction operator; this storage
!   was allocated in BRINT1, brint2,...brint6
!
      if (ltrans) then
         if (ntpi(1) > 0) then
            call dalloc (indtp1, "INDTP1", "MATRIX")
            call dalloc (valtp1, "VALTP1", "MATRIX")
         end if
         if (ntpi(2) > 0) then
            call dalloc (indtp2, "INDTP2", "MATRIX")
            call dalloc (valtp2, "VALTP2", "MATRIX")
         end if
         if (ntpi(3) > 0) then
            call dalloc (indtp3, "INDTP3", "MATRIX")
            call dalloc (valtp3, "VALTP3", "MATRIX")
         end if
         if (ntpi(4) > 0) then
            call dalloc (indtp4, "INDTP4", "MATRIX")
            call dalloc (valtp4, "VALTP4", "MATRIX")
         end if
         if (ntpi(5) > 0) then
            call dalloc (indtp5, "INDTP5", "MATRIX")
            call dalloc (valtp5, "VALTP5", "MATRIX")
         end if
         if (ntpi(6) > 0) then
            call dalloc (indtp6, "INDTP6", "MATRIX")
            call dalloc (valtp6, "VALTP6", "MATRIX")
         end if
      end if
!
!   Deallocate storage for the nuclear motional energy integral
!   lists; this was allocated in KEINT and VINT
!
      if (lnms) then
         if (nkei > 0) then
            call dalloc (indkei, "INDKEI", "MATRIX")
            call dalloc (valkei, "VALKEI", "MATRIX")
         end if
      end if
      if (lsms) then
         if (nvinti > 0) then
            call dalloc (indtei, "INDTEI", "MATRIX")
            call dalloc (valtei, "VALTIE", "MATRIX")
         end if
      end if
!
!   Deallocate storage for the vacuum-polarisation integral list;
!   this was allocated in VPINT
!
      if (lvp) then
         if (nvpi > 0) then
            call dalloc (indvpi, "INDVPI", "MATRIX")
            call dalloc (valvpi, "VALVPI", "MATRIX")
         end if
      end if

      call dalloc (radial_large_components, "PF", "MATRIX") ! lodrwf or lodres
      call dalloc (radial_small_components, "QF", "MATRIX") ! lodrwf or lodres

      return

      end subroutine matrix
