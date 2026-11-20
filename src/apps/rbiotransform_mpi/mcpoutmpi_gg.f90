!***********************************************************************
!                                                                      *
      subroutine mcpout(name, startup_directory, ik, ntestg, inpci)
!                                                                      *
!   This routine controls the computation  and storage of the values   *
!   and all indices of the angular coefficients                        *
!                                                                      *
!                                                                      *
!                   T  (ab)                                            *
!                    rs                                                *
!                                                                      *
!   k is the multipolarity of a two-particle Coulomb integral. a, b,   *
!   are orbital sequence numbers.  r and s are configuration           *
!   state function indices.                                            *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, TNSRJJ                         *
!               [GENMCP]: QSORT.                                       *
!                                                                      *
!   Written by Per Jonsson                Last revision: JUne 1996     *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_orbitals, orbital_index_radix
      use array_memory_management, only: dalloc
!      USE def_C
!      USE EIGV_C
      use first_order_csf_cutoffs, only: iccut
      use mcp_generation_options, only: diagonal_only, max_coulomb_multipole_rank, use_first_order_restriction
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
      use csf_coupling_data, only: csf_coupling_quantum_codes, csf_subshell_quantum_codes, jcupar, jqsar
      use transition_basis_blocks, only: ncfblk, nevinblk, ncfinblk, two_j, nblock, idxblk, nelectot, ncftot, nwtot, nvectot, &
          nvecsiztot, nblock1, nblocki, nblockf, ncfi, ncff
      use orbital_csf_basis,           only: csf_count, orbital_count,orbital_kappa, csf_subshell_occupations
      use biorthogonal_transition_selection, only: itjqj1, itjqj2, jqj1, jqj2, ntrans
      use biorthogonal_orbital_order, only: nordff, nordii
      use mcp_coefficient_storage, only: cnn, intgrl, intptr, jann, jbnn, ncoeff, nintg
      use biorthogonal_subshell_layout,         only: kamax,nlmax,nakinvii,nakinvff
      use mpi_state, only: mpi_process_count, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use one_electron_scalar_coefficients, only: onescalar
      use closed_shell_mcp, only: cord
      use csf_quantum_queries, only: itjpo
      use angdata_i, only: angdata
      use qqsort_i, only: qqsort
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ik
      integer, intent(in) :: ntestg
      integer  :: inpci
      character  :: name*24
      character  :: startup_directory*128
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nf = 200
      integer, parameter :: nvmax = 100
      integer, parameter :: orbital_label_radix = orbital_index_radix
      integer, parameter :: keysq = orbital_label_radix*orbital_label_radix
      real(real64), parameter :: cutoff = 1.0e-10_real64
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: llistt(nlmax)
      integer :: nakinv(max_orbitals)
      integer :: nshl(nlmax), ninl(nlmax)
      integer :: i, ntestl, ntest, ka, iopar, j, ncf0, iblk, k, &
         bra_csf_index, ket_csf_index, ibb, orbital_index, ket_orbital_index, lab, jan, jbn, l
      real(real64) :: shell_coefficients(max_orbitals)
      real(real64) :: sc(nvmax)
      real(real64) :: cirot(20*nlmax*nlmax)
      real(real64) :: evsc(10000)
      logical   :: f0int, lincr, is_restart, comp, avail
      character :: number_text*20, ck*2
!-----------------------------------------------
!
! Locals
!     POINTER (pscr, scr(1))
!     POINTER (pciout, ciout(1))
!

      write (6, *) "NBLOCK,(NCFBLK(i),i=1,NBLOCK)"
      write (6, *) nblock, (ncfblk(i),i=1,nblock)

      ntestl = 0
      ntest = max(ntestg,ntestl)
      ntest = 0
!
!   Set the rank (zero) and parity (even) for the one-particle
!   coefficients
!
      ka = 0
      iopar = 1
!
!   Check if angular data is available. If available read this data.
!   If not available calculate the data
!
      call angdata (name, avail, kamax)
      write (6, *) "AVAIL=", avail
      if (avail) return

      if (ik == 1) then
         nakinv(:max_orbitals) = nakinvii(:max_orbitals)
      else
         nakinv(:max_orbitals) = nakinvff(:max_orbitals)
      end if



      write (6, *) " open sorted ang. file .TB(NF)", nf
      j = index(name," ")
      open(unit=nf, file=name(1:j-1)//".TB", status="UNKNOWN", form=&
         "UNFORMATTED", position="asis")

      ncf0 = 1
      do iblk = 1, nblock
!
!   Open scratchfiles to dump the T coefficients for each kappa
!
         do k = 1, kamax
            open(unit=80 + k, status="UNKNOWN", form="UNFORMATTED", position=&
               "asis")
         end do
!
!   Initialize the counters for the total number of T coefficients
!
         llistt(:nlmax) = 0
!
!   JA and JB respectively refer to the initial and final states
!   in the list of NCF configurations
!
         do bra_csf_index = mpi_rank+ncf0,ncfblk(iblk),mpi_process_count
            if (mod(bra_csf_index,100)==0 .and. ik==1) then
               write (*, *) " JA1 =", bra_csf_index, bra_csf_index - ncf0 + 1
            else if (mod(bra_csf_index,100)==0 .and. ik==2) then
               write (*, *) " JA2 =", bra_csf_index, bra_csf_index - ncf0 + 1
            end if
!
            do ket_csf_index = ncf0, ncfblk(iblk)
!
!   Call the MCT package to compute T coefficients
!
               if (ntrans == 1) then
                  comp = .false.
                  if (ik == 1) then
                     do ibb = 1, jqj1
                        if (itjpo(bra_csf_index) /= itjqj1(ibb)) cycle
                        comp = .true.
                     end do
                  else
                     do ibb = 1, jqj2
                        if (itjpo(bra_csf_index) /= itjqj2(ibb)) cycle
                        comp = .true.
                     end do
                  end if
               else
                  comp = .true.
               end if


               if (.not.comp) cycle
!            write(*,*) JA,JB
!GG               CALL TNSRJJ (KA, IOPAR, JA, JB, IA, IB, TSHELL)
              call onescalar(bra_csf_index,ket_csf_index,orbital_index,ket_orbital_index,shell_coefficients)
               if (orbital_index == 0) cycle
               if (orbital_index == ket_orbital_index) then
                  do orbital_index = 1, orbital_count
!
!   If T coefficient is greater than zero and the kappa quantum numbers
!   of the two orbitals are the same dump to file
!   In a later version use a buffer with a reasonable record length
!
                     if (dabs(shell_coefficients(orbital_index)) <= cutoff) cycle
                     llistt(nakinv(orbital_index)) = llistt(nakinv(orbital_index)) + 1
                     lab = orbital_index*orbital_label_radix + orbital_index
                     write (80 + nakinv(orbital_index)) bra_csf_index - ncf0 + 1, ket_csf_index - ncf0 + 1, lab&
                        , shell_coefficients(orbital_index)
                  end do
               else
                  if (dabs(shell_coefficients(1))>cutoff .and. orbital_kappa(orbital_index)==orbital_kappa(ket_orbital_index)) then
                     llistt(nakinv(orbital_index)) = llistt(nakinv(orbital_index)) + 1
                     if (nordii==0 .and. nordff==0) then
!
!   Experssion for normal orbital ordering
!
                        lab = orbital_index*orbital_label_radix + ket_orbital_index
                        jan = bra_csf_index - ncf0 + 1
                        jbn = ket_csf_index - ncf0 + 1
                     else if (nordii==1 .and. nordff==1) then
!
!   Experssion for reversed orbital ordering
!
                        lab = ket_orbital_index*orbital_label_radix + orbital_index
                        jan = ket_csf_index - ncf0 + 1
                        jbn = bra_csf_index - ncf0 + 1
                     else
                        write (*, *) "SOMETHING WRONG"
                        stop
                     end if
                     write (80 + nakinv(orbital_index)) jan, jbn, lab, shell_coefficients(1)
                  end if
               end if
!
            end do
         end do


!
!   sort the MCP data into inegral based lists.
!
         do l = 1, kamax
            if (llistt(l) > 0) then
               call qqsort (l, llistt(l), ik, name, kamax)
            else
               if (l == 1) write (nf) csf_count, orbital_count, kamax
               write (nf) llistt(l), llistt(l)
            end if
         end do
!
!  Close the angular files
!
         do l = 1, kamax
            close(l + 80, status="DELETE")
         end do
!
         ncf0 = ncfblk(iblk) + 1
      end do
!
!   Deallocate storage that is no longer required. This was
!   allocated in lodcsl.
!
      call dalloc (csf_subshell_occupations, "IQA", "MCPOUT")
      call dalloc (csf_subshell_quantum_codes, "JQSA", "MCPOUT")
      call dalloc (csf_coupling_quantum_codes, "JCUPA", "MCPOUT")
      return
      end subroutine mcpout
