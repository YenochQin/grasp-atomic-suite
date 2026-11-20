!***********************************************************************
!                                                                      *
      subroutine lodcslbio(ncore,igg)
!                                                                      *
!   Loads the data from the  .csl  file. A number of checks are made   *
!   to ensure correctness and consistency.                             *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, IQ, ISPAR, ITJPO, JCUP, JQS,   *
!                        PACK, PARSJL, PRSRCN, PRSRSL                  *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Updated by Xinghong He                               29 Oct 1997   *
!       To accept both block and non-block formats                     *
!                                                                      *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use grasp_parameters,   only: max_orbitals
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use atomic_calculation_state, only: electron_count
      use orbital_csf_basis, only: csf_count, csf_subshell_occupations, max_orbitals, orbital_count, orbital_kappa, &
          orbital_l, orbital_principal_n, orbital_symmetry_labels, orbital_two_j
      use csf_coupling_data, only: csf_coupling_quantum_codes, csf_subshell_quantum_codes, jcupar, jqsar
      use subshell_term_table,          only: jtab, ntab
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use transition_basis_blocks,            only: nblock,ncfblk
      use array_memory_management, only: alloc, dalloc, ralloc
      use, intrinsic :: iso_fortran_env, only: iostat_end
      use mpi_io_errors, only: fail_input_io => fail_input_io_mpi
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use prsrsl_i, only: prsrsl
      use integer_format, only: convrt
      use prsrcn_i, only: prsrcn
      use parsjl_i, only: parsjl
      use csf_byte_packing, only: pack
      use csf_quantum_queries, only: iq, jqs, jcup, itjpo, ispar
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,  intent(in) ::  igg
      integer,  intent(out) :: ncore
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nw2 = 2*max_orbitals
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer  :: iocc(max_orbitals)
      integer   :: iqsub(nw2)
      integer  :: jx(max_orbitals)
      integer :: i
      integer :: ncorp1, npeel, npeel2, j, npj, nakj, lenth, ncfd, nrec &
         , ios, ierr, loc, nqs, newsiz, isparc, njx, ioc, ipty, nqsn    &
         , njxn, npeeln, nopen, jlast, ilast, iocci, nkji, ifulli, nu   &
         , jsub, iqt, nbeg, nend, jxn, jpi, ii, itemp, ncorel
      logical :: empty, full
      character          :: recl
      character(len=256) :: record
!-----------------------------------------------
!
!
!   Entry message
!
      write (6, *) "Loading Configuration Symmetry List File ..."
!
!   Get the list of subshells
!
      orbital_count = 0
!
!   Read the list of core subshells; set up the arrays NP, NAK,
!   NKL, NKJ, NH for these subshells
!
      call prsrsl (21, 1)
      ncore = orbital_count
      ncorp1 = orbital_count + 1
!
!   Skip the peel subshell identification header; read the list of
!   peel subshells; set up the arrays NP, NAK, NKL, NKJ, NH for
!   these subshells
!
      read (21, *)
      call prsrsl (21, 2)
      npeel = orbital_count - ncore
      npeel2 = npeel*2
!
!   Ensure that the sets of core and peel subshell are disjoint
!
      do j = ncore + 1, orbital_count
         npj = orbital_principal_n(j)
         nakj = orbital_kappa(j)
         do i = 1, ncore
            if (orbital_principal_n(i)/=npj .or. orbital_kappa(i)/=nakj) cycle
            write (stderr_unit, *) "LODCSL: The lists of core and", &
               " peel subshells must form disjoint sets."
            stop
         end do
      end do
!
!   Print the number of relativistic subshells
!
      if (orbital_count > 1) then
         call convrt (orbital_count, record, lenth)
         write (6, *) "There are "//record(1:lenth)// &
                      " relativistic subshells;"
      else
         write (6, *) "There is 1 relativistic subshell;"
      end if
!
!   Initial allocation for arrays with a dimension dependent
!   on the number of CSFs; the initial allocation must be
!   greater than 1
!
      if(igg == 1) then
         ncfd = 1000
         call alloc (csf_subshell_occupations,  max_orbitals,  ncfd, "IQA",   "LODCSL")
         call alloc (csf_subshell_quantum_codes, max_orbitals,3,ncfd, "JQSA",  "LODCSL")
         call alloc (csf_coupling_quantum_codes,max_orbitals,  ncfd, "JCUPA", "LODCSL")
      else if(igg == 2) then
         ncfd = 1000 + csf_count
         call ralloc (csf_subshell_occupations,  max_orbitals,  ncfd, "IQA",   "LODCSL")
         call ralloc (csf_subshell_quantum_codes, max_orbitals,3,ncfd, "JQSA",  "LODCSL")
         call ralloc (csf_coupling_quantum_codes,max_orbitals,  ncfd, "JCUPA", "LODCSL")
      else if(igg == 3) then
         call dalloc (csf_subshell_occupations,   "IQA",   "LODCSL")
         call dalloc (csf_subshell_quantum_codes,  "JQSA",  "LODCSL")
         call dalloc (csf_coupling_quantum_codes, "JCUPA", "LODCSL")
         ncfd = 1000
         call alloc (csf_subshell_occupations,  max_orbitals,  ncfd, "IQA",   "LODCSL")
         call alloc (csf_subshell_quantum_codes, max_orbitals,3,ncfd, "JQSA",  "LODCSL")
         call alloc (csf_coupling_quantum_codes,max_orbitals,  ncfd, "JCUPA", "LODCSL")
      end if

!
!   Skip the header for the list of CSFs
!
      read (21, *)
!
!   NREC is the sequence number of the last record read in the
!   Configuration Symmetry List File
!
      nrec = 5
!
!   There must be three records for each CSF: For instance,
!
!    4s ( 2) 4p-( 2) 4p ( 3) 4d-( 2) 4d ( 5) 4f-( 6) 4f ( 4)
!                        3/2       0     5/2             2;4
!                                           1               3-
!   Each record is as follows:
!      (1) Peel subshell occupation number (q) specification.
!      (2) Subshell total angular momentum quantum number
!          (v_, J_sub) specifications. Subshell total angular
!          momentum quantum numbers (J_sub) must be specified
!          for all open subshells, even if this can be
!          deduced from the subshell occupation. Seniority
!          quantum numbers (v) must be specified for
!          subshells with j = 7/2, q = 4, J_sub = 2 or 4.
!          The seniority quantum number must precede the
!          total angular momentum quantum number: in the
!          example above, for the 4f_7/2 subshell, q = 4
!          and J_sub = 4, whence it is necessary to specify
!          v --- 2 in this case.
!      (3) The minimum number of intermediate angular
!          momentum quantum numbers (X) as well as the final
!          angular momentum quantum number (J) immediately
!          followed by the sign of the parity (P) must be
!          specified on this record. In the example above,
!          the J_sub = 0 for the 4d_3/2 subshell, whence
!          it is unnecessary to specify its coupling to all
!          preceding subshells.
!   These conventions have been chosen so as to render the CSF
!   specifications easily interpreted by the user
!
      csf_count = 0
      nblock = 0
    3 continue
      csf_count = csf_count + 1
!
      read (21, "(A)", iostat=ios) record
!**********************************************************************
!blk*
!   To skip the border line added to mark the end of a block
!
      if (ios == 0) then
         if (record(1:2) == " *") then
            nblock = nblock + 1
            ncfblk(nblock) = csf_count -1
            read (21, "(A)", iostat=ios) record
         end if
      end if
      if (ios /= 0 .and. ios /= iostat_end) &
         call fail_input_io(ios, 21, "CSF lodcslbio occupation record")
!**********************************************************************

      if (ios == 0) then
!
!   Read in the occupations (q) of the peel shells; stop with a
!   message if an error occurs
!
         call prsrcn (record, ncore, iocc, ierr)
         if (ierr /= 0) goto 26
!
!   Read the J_sub and v quantum numbers
!
         read (21, "(A)", iostat=ios) record
         if (ios /= 0) then
            write (stderr_unit, *) "LODCSL: Expecting subshell quantum", &
               " number specification;"
            call fail_input_io(ios, 21, "CSF lodcslbio subshell quantum numbers")
         end if
         loc = len_trim(record)
         call parsjl (1, ncore, record, loc, iqsub, nqs, ierr)
         if (ierr /= 0) goto 26
!
!   Read the X, J, and (sign of) P quantum numbers
!
         read (21, "(A)", iostat=ios) record
         if (ios /= 0) then
            write (stderr_unit, *) "LODCSL: Expecting intermediate ", &
               "and final angular momentum"
            write (stderr_unit, *) "quantum number and final parity ", &
               "specification;"
            call fail_input_io(ios, 21, "CSF lodcslbio coupling quantum numbers")
         end if
!
!   Allocate additional storage if necessary
!
!CFF     It is possible that this should be moved to "3 Continue"
!        where NCF is incremented
         if (csf_count > ncfd) then
            newsiz = ncfd + ncfd/2
            call ralloc (csf_subshell_occupations,  max_orbitals,  newsiz, "IQA",   "LODCSL")
            call ralloc (csf_subshell_quantum_codes, max_orbitals,3,newsiz, "JQSA",  "LODCSL")
            call ralloc (csf_coupling_quantum_codes,max_orbitals,  newsiz, "JCUPA", "LODCSL")
            ncfd = newsiz
         end if
!
!   Zero out the arrays that store packed integers
!
         do i = 1,max_orbitals
            csf_subshell_occupations(i,csf_count)    = 0
            csf_subshell_quantum_codes(i,1,csf_count) = 0
            csf_subshell_quantum_codes(i,2,csf_count) = 0
            csf_subshell_quantum_codes(i,3,csf_count) = 0
            csf_coupling_quantum_codes(i,csf_count)  = 0
         end do
!
!   Determine the parity and all intermediate and the final
!   angular momentum quantum numbers
!
         do i = 256, 1, -1
            if (record(i:i) == " ") cycle
            loc = i
            exit
         end do
         recl = record(loc:loc)
         if (recl == "+") then
            isparc = 1
         else if (recl == "-") then
            isparc = -1
         else
            write (stderr_unit, *) "LODCSL: Incorrect parity ", &
                             "specification;"
            goto 26
         end if
         loc = loc - 1
!
         call parsjl (2, ncore, record, loc, jx, njx, ierr)
         if (ierr /= 0) goto 26
!
!   Set the occupation and subshell quantum number array elements
!   in IQ, JQS for the core subshells
!
         do i = 1, ncore
            call pack (orbital_two_j(i) + 1, i, csf_subshell_occupations(1:max_orbitals,csf_count))
            call pack (0, i, csf_subshell_quantum_codes(1:max_orbitals,1,csf_count))
            call pack (0, i, csf_subshell_quantum_codes(1:max_orbitals,2,csf_count))
            call pack (1, i, csf_subshell_quantum_codes(1:max_orbitals,3,csf_count))
         end do
!
!   Check all subshell, intermediate and final angular momentum
!   quantum numbers; set the array elements in IQ, JQS for the peel
!   subshells; set the coupling array element in JCUP and the total
!   angular momentum array element in ITJPO
!
         ioc = 0
         ipty = 0
         nqsn = 0
         njxn = 0
         npeeln = 0
         nopen = 0
         jlast = 0
         ilast = 0
         do i = ncorp1, orbital_count
            iocci = iocc(i)
            npeeln = npeeln + iocci
            nkji = orbital_two_j(i)
            ifulli = nkji + 1
            empty = iocci == 0
            if (.not.empty) ioc = ioc + 1
            full = iocci == ifulli
            if (empty .or. full) then
               nu = 0
               jsub = 0
            else
               ipty = ipty + orbital_l(i)*iocci
               if (nkji /= 7) then
                  nqsn = nqsn + 1
                  if (nqsn > nqs) then
                     write (stderr_unit, *) "LODCSL: Too few subshell quantum", &
                        " numbers specified;"
                     goto 26
                  end if
                  nu = 0
                  jsub = iqsub(nqsn)
               else
                  if (iocci /= 4) then
                     nqsn = nqsn + 1
                     if (nqsn > nqs) then
                        write (stderr_unit, *) "LODCSL: Too few subshell ", &
                           "quantum numbers specified;"
                        goto 26
                     end if
                     nu = 0
                     jsub = iqsub(nqsn)
                  else
                     nqsn = nqsn + 1
                     if (nqsn > nqs) then
                        write (stderr_unit, *) "LODCSL: Too few subshell ", &
                           "quantum numbers specified;"
                        goto 26
                     end if
                     jsub = iqsub(nqsn)
                     if (jsub==4 .or. jsub==8) then
                        nu = jsub/2
                        nqsn = nqsn + 1
                        if (nqsn > nqs) then
                           write (stderr_unit, *) "LODCSL: Too few subshell", &
                              " quantum numbers specified;"
                           goto 26
                        end if
                        jsub = iqsub(nqsn)
                     else
                        nu = 0
                     end if
                  end if
               end if
               iqt = min(iocci,ifulli - iocci)
               loc = (ifulli - 2)/2
               loc = (loc*(loc + 1))/2 + iqt
               nbeg = jtab(loc+1) + 1
               nend = jtab(loc+2)
               do j = nbeg, nend, 3
                  if (ntab(j+2) /= jsub + 1) cycle
                  if (nu == 0) then
                     nu = ntab(j)
                     goto 9
                  else
                     if (ntab(j) == nu) goto 9
                  end if
               end do
               call convrt (orbital_principal_n(i), record, lenth)
               write (stderr_unit, *) "LODCSL: Subshell quantum numbers ", &
                  "specified incorrectly for "//record(1:lenth)//orbital_symmetry_labels(i)//&
                  " subshell."
               goto 26
            end if
    9       continue
            if (.not.empty .and. .not.full) then
               nopen = nopen + 1
               if (nopen > 1) then
                  if (jsub == 0) then
                     jxn = jlast
                  else
                     ilast = ioc
                     njxn = njxn + 1
                     if (njxn > njx) then
                        write (stderr_unit, *) "LODCSL: Too few intermediate", &
                           " and final angular momentum", &
                           " quantum numbers specified;"
                        goto 26
                     end if
                     jxn = jx(njxn)
                     do j = abs(jlast - jsub), jlast + jsub, 2
                        if (jxn == j) goto 11
                     end do
                     call convrt (orbital_principal_n(i), record, lenth)
                     write (stderr_unit, *) &
                        "LODCSL: coupling of "//record(1:lenth)//orbital_symmetry_labels(i),&
                        " subshell to previous subshells is incorrect."
                     goto 26
                  end if
   11             continue
                  call pack (jxn + 1, nopen - 1, csf_coupling_quantum_codes(1:max_orbitals,csf_count))
                  jlast = jxn
               else
                  jlast = jsub
               end if
            end if
            call pack (iocci, i, csf_subshell_occupations(1:max_orbitals,csf_count))
            call pack (nu, i, csf_subshell_quantum_codes(1:max_orbitals,1,csf_count))
            call pack (0, i, csf_subshell_quantum_codes(1:max_orbitals,2,csf_count))
            call pack (jsub + 1, i, csf_subshell_quantum_codes(1:max_orbitals,3,csf_count))
         end do
!
         do i = max(1,nopen), orbital_count
            call pack (0, i, csf_coupling_quantum_codes(1:max_orbitals,csf_count))
         end do
!
         if (nqsn /= nqs) then
            write (stderr_unit, *) "LODCSL: Too many subshell", &
               " quantum numbers specified;"
            goto 26
         end if
!
         if (ilast /= ioc) njxn = njxn + 1
         if (njxn /= njx) then
            write (stderr_unit, *) "LODCSL: Too many intermediate", &
            " and final angular momentum", " quantum numbers specified;"
            goto 26
         end if
!
         if (jx(njxn) /= jlast) then
            write (stderr_unit, *) "LODCSL: Final angular momentum", &
               " incorrectly specified;"
            goto 26
         end if
!
         ipty = (-1)**ipty
         if (ipty /= isparc) then
            write (stderr_unit, *) "LODCSL: Parity specified incorrectly;"
            goto 26
         end if
!
         jpi = (jlast + 1)*ipty
         call pack (jpi, max_orbitals, csf_coupling_quantum_codes(1:max_orbitals,csf_count))
!
         if (csf_count > 1) then
            if (npeeln /= npeel) then
               write (stderr_unit, *) "LODCSL: Inconsistency in the number", &
                  " of electrons."
               goto 26
            end if
         else
            npeel = npeeln
         end if
!
!   Check if this CSF was already in the list; stop with a
!   message if this is the case
!
         if (csf_count > 1) then
            do j = 1, csf_count - 1
               do i = ncorp1, orbital_count
                  if (iq(i,j) /= iq(i,csf_count)) goto 17
                  if (jqs(1,i,j) /= jqs(1,i,csf_count)) goto 17
                  if (jqs(2,i,j) /= jqs(2,i,csf_count)) goto 17
                  if (jqs(3,i,j) /= jqs(3,i,csf_count)) goto 17
               end do
               do i = 1, nopen - 1
                  if (jcup(i,j) /= jcup(i,csf_count)) goto 17
               end do
            end do
            write (stderr_unit, *) "LODCSL: Repeated CSF;"
            goto 26
         end if
!
!   Successfully read a CSF; update NREC and read another CSF
!
   17    continue
         nrec = nrec + 3
         goto 3
!
      else
!
!   There is always at least one CSF
!
         if (csf_count == 1) then
            do i = 1, ncore
               call pack (orbital_two_j(i) + 1, i, csf_subshell_occupations(1:max_orbitals,1))
               call pack (0, i, csf_subshell_quantum_codes(1:max_orbitals,1,1))
               call pack (0, i, csf_subshell_quantum_codes(1:max_orbitals,2,1))
               call pack (1, i, csf_subshell_quantum_codes(1:max_orbitals,3,1))
            end do
            call pack (0, 1, csf_coupling_quantum_codes(1:max_orbitals,1))
            call pack (1, max_orbitals, csf_coupling_quantum_codes(1:max_orbitals,1))
         else
            csf_count = csf_count - 1
         end if
!
      end if
!
!   Check if any subshell is empty; eliminate it from the
!   list if this is the case; issue a message
!
      i = ncorp1
   19 continue
      if (i <= orbital_count) then
         do j = 1, csf_count
            if (iq(i,j) /= 0) goto 23
         end do
         call convrt (orbital_principal_n(i), record, lenth)
         write (6, *) "Subshell "//record(1:lenth)//orbital_symmetry_labels(i)//" is empty", &
            " in all CSFs; eliminating this", " subshell from the list;"
         orbital_count = orbital_count - 1
         do ii = i, orbital_count
            orbital_principal_n(ii) = orbital_principal_n(ii+1)
            orbital_kappa(ii) = orbital_kappa(ii+1)
            orbital_l(ii) = orbital_l(ii+1)
            orbital_two_j(ii) = orbital_two_j(ii+1)
            orbital_symmetry_labels(ii) = orbital_symmetry_labels(ii+1)
            do j = 1, csf_count
               itemp = iq(ii + 1,j)
               call pack (itemp, ii, csf_subshell_occupations(1:max_orbitals,j))
               itemp = jqs(1,ii + 1,j)
               call pack (itemp, ii, csf_subshell_quantum_codes(ii:max_orbitals,1,j))
               itemp = jqs(2,ii + 1,j)
               call pack (itemp, ii, csf_subshell_quantum_codes(ii:max_orbitals,2,j))
               itemp = jqs(3,ii + 1,j)
               call pack (itemp, ii, csf_subshell_quantum_codes(ii:max_orbitals,3,j))
            end do
         end do
   23    continue
         i = i + 1
         goto 19
      end if
!
!   Store the number of electrons in the COMMON variable
!
      ncorel = 0
      ncorel = sum(orbital_two_j(:ncore)+1)
      electron_count = ncorel + npeel
!
!   All done; report
!
      call convrt (csf_count, record, lenth)
      write (6, *) "There are "//record(1:lenth)//" relativistic CSFs;"
      write (6, *) " ... load complete;"
!
!   Debug printout
!
      if (ldbpa(1)) then
         write (99, *) "From LODCSL:"
         do i = 1, csf_count
            write (99, *) "CSF ", i
            write (99, *) "ITJPO: ", itjpo(i)
            write (99, *) "ISPAR: ", ispar(i)
            write (99, *) "IQ: ", (iq(j,i),j=1,orbital_count)
            write (99, *) "JQS(1): ", (jqs(1,j,i),j=1,orbital_count)
            write (99, *) "JQS(2): ", (jqs(2,j,i),j=1,orbital_count)
            write (99, *) "JQS(3): ", (jqs(3,j,i),j=1,orbital_count)
            write (99, *) "JCUP: ", (jcup(j,i),j=1,orbital_count - 1)
         end do
      end if


      nblock = nblock + 1
      ncfblk(nblock) = csf_count
!
      return
!
   26 continue
      call convrt (csf_count, record, lenth)
      write (stderr_unit, *) " CSF sequence number: "//record(1:lenth)//":"
      rewind (21)
      do i = 1, nrec
         read (21, *)
      end do
      do i = 1, 3
         read (21,"(A)",err = 29,end = 29) record
         lenth = len_trim(record)
         write (stderr_unit, *) record(1:lenth)
      end do
   29 close(21)
      stop
!
      end subroutine lodcslbio
