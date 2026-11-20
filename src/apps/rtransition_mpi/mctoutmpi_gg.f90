!***********************************************************************
!                                                                      *
      subroutine mctout(operator_parity, transition_type_index, name)
!                                                                      *
!   This routine loads  coefficients with parity and  rank specified   *
!   by KP(JKP) into the arrays ISLDR and XSLDR.  IOPAR is the parity   *
!   (+/- 1) and is determined from the sign of  KP(JKP).               *
!                                                                      *
!                                         Last revision: 28 Dec 1992   *
!   Updated by Jacek Bieron               Last revision: 10 Mar 1994   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: orbital_index_radix, max_orbitals
      use array_memory_management, only: alloc, dalloc, ralloc
      use transition_basis_blocks, only: ncfblk, nevinblk, ncfinblk, two_j, nblock, idxblk, nelectot, ncftot, nwtot, nvectot, &
          nvecsiztot, nblock1, nblocki, nblockf, ncfi, ncff
      use debug_output_settings,         only: ldbpa
      use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
      use first_order_csf_cutoffs, only: iccut
      use orbital_csf_basis, only: max_orbitals, orbital_count
      use transition_off_diagonal_flags, only: noffd1, noffd2
      use transition_matrix_workspace, only: nsdim, hb1, hb2, hc1, hc2, hm1, hm2, jja, jjb, lk, kk, ntdim, &
          xsldr, totc, totb, isldr, isldr1, nint, ninteg, nptr, lab, nkp, kp, ltc
      use mpi_state, only: mpi_process_count, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use angdata_i, only: angdata
      use csf_quantum_queries, only: itjpo
      use one_electron_tensor_coefficients, only: oneparticlejj
      use trsort_i, only: trsort
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: operator_parity
      integer, intent(in) :: transition_type_index
      character  :: name(2)*24
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(real64), parameter :: cutoff = 1.0e-10_real64
      integer, parameter :: nfile = 93
      integer, parameter :: nfile1 = 237
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nfile2, iblki, iblkf, nmct, nlabel, ncfi0, ncff0, bra_csf_index, ket_csf_index, ncr&
         , orbital_index, ket_orbital_index, newsiz, i
      logical :: cache_exists
      integer, pointer :: label(:)
      real(real64), pointer:: coeff(:)
      real(real64) :: shell_coefficients(max_orbitals)
!-----------------------------------------------
!
!
!
! NCFI(I): the end position of the Ith block for the initial states in the globle CSF list
! NCFF(I): the end position of the Ith block for the final states in the globle CSF list
!
!   Check if angular data is available on file
!
      nfile2 = nfile1 + transition_type_index
      call angdata (name, cache_exists, transition_type_index, nfile2)

!
!   If angular data is not available open the scratch file to store the
!   coefficients; position file
!   to beginning
!
!
      lk = abs(kp(transition_type_index))
      operator_parity = isign(1,kp(transition_type_index))
      if (cache_exists) return
!      WRITE (6, *) 'LK,IOPAR,from MCTOUT'
!      WRITE (6, *) LK, IOPAR
!
! Start of the block loops
      do iblki = 1, nblocki
         do iblkf = 1, nblockf
!        OPEN (NFILE,STATUS = 'new', FORM = 'UNFORMATTED')
!        Sometimes, when there has been an error, status need to be "unknown"
!CFF
            open(nfile,status="unknown",form="UNFORMATTED")
            nmct = 0
!
!
!   If angular data is not available
!   Generate MCT coefficients for given rank/parity combination and
!   store them by CSF on NFILE
!
!
!   Allocate storage to buffer arrays
!
            nlabel = 32
            call alloc (label, nlabel, "LABEL", "MCTOUT")
            call alloc (coeff, nlabel, "COEFF", "MCTOUT")

            if (iblki == 1) then
               ncfi0 = 1
            else
               ncfi0 = ncfi(iblki-1) + 1
            end if
!
            if (iblkf == 1) then
               ncff0 = ncfi(nblocki) + 1
            else
               ncff0 = ncfi(nblocki) + ncff(iblkf-1) + 1
            end if
            do bra_csf_index = mpi_rank+ncfi0, ncfi(iblki), mpi_process_count
               do ket_csf_index = ncff0, ncfi(nblocki) + ncff(iblkf)
!
!           IR = IC
!
                  ncr = 0

!
!   In many case one is interested only in M1 and E2 transitions between
!   levels with different J values. If this is the case then the do check
!   on the J quantum numbers of the CSFs before calling TNSRJJ.
!
                  if (kp(transition_type_index)==1 .and. noffd1==1) then
                     if (itjpo(bra_csf_index) == itjpo(ket_csf_index)) cycle
                  end if
                  if (kp(transition_type_index)==2 .and. noffd2==1) then
                     if (itjpo(bra_csf_index) == itjpo(ket_csf_index)) cycle
                  end if
!          if(ispar(ic)*ispar(ir)*iopar.ne.1.
!    &        or.itrig(itjpo(ic),itjpo(ir),2*lk+1).ne.1) go to 13
!          if(ichkq1(IC,IR).eq.0) go to 13
                  call oneparticlejj(lk,operator_parity,bra_csf_index,ket_csf_index,orbital_index,ket_orbital_index, &
                     & shell_coefficients)
                  if (orbital_index /= 0) then
                     if (orbital_index == ket_orbital_index) then
                        do orbital_index = 1, orbital_count
                           if (abs(shell_coefficients(orbital_index)) <= cutoff) cycle
                           ncr = ncr + 1
                           if (ncr > nlabel) then
                              newsiz = 2*nlabel
                              call ralloc (label,  newsiz, "LABEL", "MCTOUT")
                              call ralloc (coeff,  newsiz, "COEFF", "MCTOUT")
                              nlabel = newsiz
                           end if
                           label(ncr) = orbital_index*orbital_index_radix + orbital_index
                           coeff(ncr) = shell_coefficients(orbital_index)
                        end do
                     else
                        if (abs(shell_coefficients(1)) > cutoff) then
                           ncr = ncr + 1
                           if (ncr > nlabel) then
                              newsiz = 2*nlabel
                              call ralloc (label,  newsiz, "LABEL", "MCTOUT")
                              call ralloc (coeff,  newsiz, "COEFF", "MCTOUT")
                              nlabel = newsiz
                           end if
                           label(ncr) = orbital_index*orbital_index_radix + ket_orbital_index
                           coeff(ncr) = shell_coefficients(1)
                        end if
                     end if
                  end if
                  if (ncr <= 0) cycle

                  write (nfile) bra_csf_index - ncfi0 + 1, ket_csf_index - ncff0 + 1, ncr
                  write (nfile) (label(i),coeff(i),i=1,ncr)
                  nmct = nmct + ncr

!
               end do
            end do
!
!   Deallocate storage for buffer arrays
!
            call dalloc (label,  "LABEL", "MCTOUT")
            call dalloc (coeff,  "COEFF", "MCTOUT")
!
!            WRITE (*, 301) NMCT, LK, IOPAR
!
!   Sort the MCT coefficients by integral labels
!
            call trsort (name, nfile, nfile2, ldbpa(2), transition_type_index, iblki, iblkf)
            close(nfile, status="delete")
! end of the loops for blocks
         end do
      end do
!
!   Read the data back as required by OSCL conventions
!
      rewind (nfile2)
      return
!
      return
!
      end subroutine mctout
