!***********************************************************************
!                                                                      *
      subroutine lodcslm(ncore)
!                                                                      *
!   Loads the data from the  .csl  file. A number of checks are made   *
!   to ensure correctness and consistency.                             *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Updated by Xinghong He                               31 Oct 1997
!       To accept both block and non-block formats
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:21:55   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: max_orbitals
      use memory_man
      use debug_c
      use def_c
      use orb_c
      use stat_c
      use terms_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use prsrsl_i
      use convrt_i
      use prsrcn_i
      use parsjl_i
      use pack_i
      use iq_i
      use jqs_i
      use jcup_i
      use itjpo_i
      use ispar_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ncore
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nw2 = 2*max_orbitals
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(max_orbitals) :: iocc
      integer , dimension(nw2) :: iqsub
      integer , dimension(max_orbitals) :: jx
      integer :: ncorp1, npeel, npeel2, j, npj, nakj, i, lenth, ncfd, nrec, ios&
         , ierr, loc, nqs, newsiz, isparc, njx, ioc, ipty, nqsn, njxn, npeeln, &
         nopen, jlast, ilast, iocci, nkji, ifulli, nu, jsub, iqt, nbeg, nend, &
         jxn, jpi, ii, itemp, ncorel
      logical :: empty, full
      character :: record*256, recl
!-----------------------------------------------
!
!
!   Entry message
!
!      WRITE (6, *) 'Loading Configuration Symmetry List File ...'
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
            write (6, *) 'LODCSL: The lists of core and'
            write (6, *) ' peel subshells must form'
            write (6, *) ' disjoint sets.'
            stop
         end do
      end do
!
!   Print the number of relativistic subshells
!
      if (orbital_count > 1) then
         call convrt (orbital_count, record, lenth)
!         WRITE (6, *) ' there are '//RECORD(1:LENTH)//&
!            ' relativistic subshells;'
      else
         write (6, *) ' there is 1 relativistic subshell;'
      endif
!
!   Initial allocation for arrays with a dimension dependent
!   on the number of CSFs; the initial allocation must be
!   greater than 1
!
      call dalloc (csf_subshell_occupations,   'IQA',   'LODCSLM')
      call dalloc (csf_subshell_quantum_codes,  'JQSA',  'LODCSLM')
      call dalloc (csf_coupling_quantum_codes, 'JCUPA', 'LODCSLM')
      ncfd = 2
      call alloc (csf_subshell_occupations,   max_orbitals,    ncfd, 'IQA',   'LODCSLM')
      call alloc (csf_subshell_quantum_codes,  max_orbitals, 3, ncfd, 'JQSA',  'LODCSLM')
      call alloc (csf_coupling_quantum_codes, max_orbitals,    ncfd, 'JCUPA', 'LODCSLM')
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
    3 continue
      csf_count = csf_count + 1
!
      read (21, '(A)', iostat=ios) record
!**********************************************************************
!blk*
!   To skip the border line added to mark the end of a block
!
      if (record(1:2) == ' *') read (21, '(A)', iostat=ios) record
!**********************************************************************

      if (ios == 0) then
!
!   Read in the occupations (q) of the peel shells; stop with a
!   message if an error occurs
!
         call prsrcn (record, ncore, iocc, ierr)
         if (ierr /= 0) go to 26
!
!   Read the J_sub and v quantum numbers
!
         read (21, '(A)', iostat=ios) record
         if (ios /= 0) then
            write (6, *) 'LODCSL: Expecting subshell quantum'
            write (6, *) ' number specification;'
            go to 26
         endif
         loc = len_trim(record)
         call parsjl (1, ncore, record, loc, iqsub, nqs, ierr)
         if (ierr /= 0) go to 26
!
!   Read the X, J, and (sign of) P quantum numbers
!
         read (21, '(A)', iostat=ios) record
         if (ios /= 0) then
            write (6, *) 'LODCSL: Expecting intermediate'
            write (6, *) ' and final angular momentum'
            write (6, *) ' quantum number and final parity'
            write (6, *) ' specification;'
            go to 26
         endif
!
!   Allocate additional storage if necessary
!
         if (csf_count > ncfd) then
            newsiz = ncfd + ncfd/2
            call ralloc (csf_subshell_occupations, max_orbitals, newsiz, 'IQA', 'LODCSLM')
            call ralloc (csf_subshell_quantum_codes, max_orbitals, 3, newsiz, 'JQSA', 'LODCSLM')
            call ralloc (csf_coupling_quantum_codes, max_orbitals, newsiz, 'JCUPA', 'LODCSLM')
            ncfd = newsiz
         endif
!
!   Zero out the arrays that store packed integers
!
         csf_subshell_occupations(:max_orbitals,csf_count) = 0
         csf_subshell_quantum_codes(:max_orbitals,1,csf_count) = 0
         csf_subshell_quantum_codes(:max_orbitals,2,csf_count) = 0
         csf_subshell_quantum_codes(:max_orbitals,3,csf_count) = 0
         csf_coupling_quantum_codes(:max_orbitals,csf_count) = 0
!
!   Determine the parity and all intermediate and the final
!   angular momentum quantum numbers
!
         do i = 256, 1, -1
            if (record(i:i) == ' ') cycle
            loc = i
            exit
         end do
         recl = record(loc:loc)
         if (recl == '+') then
            isparc = 1
         else if (recl == '-') then
            isparc = -1
         else
            write (6, *) 'LODCSL: Incorrect parity'
            write (6, *) ' specification;'
            go to 26
         endif
         loc = loc - 1
!
         call parsjl (2, ncore, record, loc, jx, njx, ierr)
         if (ierr /= 0) go to 26
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
                     write (6, *) 'LODCSL: Too few subshell'
                     write (6, *) ' quantum numbers specified;'
                     go to 26
                  endif
                  nu = 0
                  jsub = iqsub(nqsn)
               else
                  if (iocci /= 4) then
                     nqsn = nqsn + 1
                     if (nqsn > nqs) then
                        write (6, *) 'LODCSL: Too few subshell'
                        write (6, *) ' quantum numbers specified;'
                        go to 26
                     endif
                     nu = 0
                     jsub = iqsub(nqsn)
                  else
                     nqsn = nqsn + 1
                     if (nqsn > nqs) then
                        write (6, *) 'LODCSL: Too few subshell'
                        write (6, *) ' quantum numbers specified;'
                        go to 26
                     endif
                     jsub = iqsub(nqsn)
                     if (jsub==4 .or. jsub==8) then
                        nu = jsub/2
                        nqsn = nqsn + 1
                        if (nqsn > nqs) then
                           write (6, *) 'LODCSL: Too few subshell'
                           write (6, *) ' quantum numbers specified;'
                           go to 26
                        endif
                        jsub = iqsub(nqsn)
                     else
                        nu = 0
                     endif
                  endif
               endif
               iqt = min(iocci,ifulli - iocci)
               loc = (ifulli - 2)/2
               loc = (loc*(loc + 1))/2 + iqt
               nbeg = jtab(loc+1) + 1
               nend = jtab(loc+2)
               do j = nbeg, nend, 3
                  if (ntab(j+2) /= jsub + 1) cycle
                  if (nu == 0) then
                     nu = ntab(j)
                     go to 9
                  else
                     if (ntab(j) == nu) go to 9
                  endif
               end do
               call convrt (orbital_principal_n(i), record, lenth)
               write (6, *) 'LODCSL: Subshell quantum numbers'
               write (6, *) ' specified incorrectly for'
               write (6, *) ' '//record(1:lenth)//orbital_symmetry_labels(i)//' subshell.'
               go to 26
            endif
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
                        write (6, *) 'LODCSL: Too few intermediate'
                        write (6, *) ' and final angular momentum'
                        write (6, *) ' quantum numbers specified;'
                        go to 26
                     endif
                     jxn = jx(njxn)
                     do j = abs(jlast - jsub), jlast + jsub, 2
                        if (jxn == j) go to 11
                     end do
                     call convrt (orbital_principal_n(i), record, lenth)
                     write (6, *) 'LODCSL: coupling of '//record(1:lenth)//orbital_symmetry_labels(i&
                        )
                     write (6, *) ' subshell to previous subshells'
                     write (6, *) ' is incorrect.'
                     go to 26
                  endif
   11             continue
                  call pack (jxn + 1, nopen - 1, csf_coupling_quantum_codes(1:max_orbitals,csf_count))
                  jlast = jxn
               else
                  jlast = jsub
               endif
            endif
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
            write (6, *) 'LODCSL: Too many subshell'
            write (6, *) ' quantum numbers specified;'
            go to 26
         endif
!
         if (ilast /= ioc) njxn = njxn + 1
         if (njxn /= njx) then
            write (6, *) 'LODCSL: Too many intermediate'
            write (6, *) ' and final angular momentum'
            write (6, *) ' quantum numbers specified;'
            go to 26
         endif
!
         if (jx(njxn) /= jlast) then
            write (6, *) 'LODCSL: Final angular momentum'
            write (6, *) ' incorrectly specified;'
            go to 26
         endif
!
         ipty = (-1)**ipty
         if (ipty /= isparc) then
            write (6, *) 'LODCSL: Parity specified incorrectly;'
            go to 26
         endif
!
         jpi = (jlast + 1)*ipty
         call pack (jpi, max_orbitals, csf_coupling_quantum_codes(1:max_orbitals,csf_count))
!
         if (csf_count > 1) then
            if (npeeln /= npeel) then
               write (6, *) 'LODCSL: Inconsistency in the number'
               write (6, *) ' of electrons.'
               go to 26
            endif
         else
            npeel = npeeln
         endif
!
!   Check if this CSF was already in the list; stop with a
!   message if this is the case
!
!         print *, 'Check duplicated CSFs'
         if (csf_count > 1) then
            do j = 1, csf_count - 1
!                  print *,'j= ',j,ncf
               do i = ncorp1, orbital_count
!                  print *,i
!                  print *, IQ(I,J), JQS(1,I,J), JQS(2,I,J),JQS(3,I,J)
!          print *, IQ(I,ncf), JQS(1,I,ncf), JQS(2,I,ncf),JQS(3,I,ncf)
                  if (iq(i,j) /= iq(i,csf_count)) go to 17
                  if (jqs(1,i,j) /= jqs(1,i,csf_count)) go to 17
                  if (jqs(2,i,j) /= jqs(2,i,csf_count)) go to 17
                  if (jqs(3,i,j) /= jqs(3,i,csf_count)) go to 17
               end do
               do i = 1, nopen - 1
!                  WRITE (6, *) I
!                  WRITE (6, *) JCUP(I,J), JCUP(I,NCF)
                  if (jcup(i,j) /= jcup(i,csf_count)) go to 17
               end do
            end do
            write (6, *) 'LODCSL: Repeated CSF;'
            go to 26
         endif
!
!   Successfully read a CSF; update NREC and read another CSF
!
   17    continue
         nrec = nrec + 3
         go to 3
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
         endif
!
      endif
!
!   Check if any subshell is empty; eliminate it from the
!   list if this is the case; issue a message
!
      i = ncorp1
   19 continue
      if (i <= orbital_count) then
         do j = 1, csf_count
            if (iq(i,j) /= 0) go to 23
         end do
         call convrt (orbital_principal_n(i), record, lenth)
         write (6, *) 'Subshell '//record(1:lenth)//orbital_symmetry_labels(i)//' is empty'
         write (6, *) ' in all CSFs; eliminating this'
         write (6, *) ' subshell from the list;'
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
         go to 19
      endif
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
!      WRITE (6, *) ' there are '//RECORD(1:LENTH)//' relativistic CSFs;'
!      WRITE (6, *) ' ... load complete;'
!
!   Debug printout
!
      if (ldbpa(1)) then
         write (99, *) 'From LODCSL:'
         do i = 1, csf_count
            write (99, *) 'CSF ', i
            write (99, *) 'ITJPO: ', itjpo(i)
            write (99, *) 'ISPAR: ', ispar(i)
            write (99, *) 'IQ: ', (iq(j,i),j=1,orbital_count)
            write (99, *) 'JQS(1): ', (jqs(1,j,i),j=1,orbital_count)
            write (99, *) 'JQS(2): ', (jqs(2,j,i),j=1,orbital_count)
            write (99, *) 'JQS(3): ', (jqs(3,j,i),j=1,orbital_count)
            write (99, *) 'JCUP: ', (jcup(j,i),j=1,orbital_count - 1)
         end do
      endif
!
      return
!
   26 continue
      call convrt (csf_count, record, lenth)
      write (6, *) ' CSF sequence number: '//record(1:lenth)//':'
      rewind (21)
      do i = 1, nrec
         read (21, *)
      end do
      do i = 1, 3
         read (21, '(A)', err=29, end=29) record
         lenth = len_trim(record)
         write (6, *) record(1:lenth)
      end do
   29 continue
      close(21)
      stop
!
      end subroutine lodcslm
