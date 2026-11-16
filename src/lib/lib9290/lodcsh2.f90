module lodcsh2_i
    implicit none (type, external)
    private
    public :: lodcsh2
contains
!***********************************************************************
!                                                                      *
      subroutine lodcsh2(nfile, ncore, jb)
!
! IMPORTANT:
! ==========
!   If jb == LOADALL, then it loads ALL blocks. In this case, the
!   NCFblock in the common should be the total ncf and allocations
!   should be changed to this ncf. Both should be done before calling
!   this routine - like the genuine block case.
!
!   Loads the data from the  .csl  file for the current block          *
!   All parameters are inputs:
!
!     nfile:    The unit number, usually 21
!     jb:       The block number (see above)
!     ncore:    The number of the core (sub-)shells.
!   Since NCFblock is known in the block version, allocation is done
!   once outside this routine. NCFblock is checked against the one (NCF)
!   obtained here.
!
!   Call(s) to: [LIB92]: CONVRT, IQ, ISPAR, ITJPO, JCUP, JQS,
!                        PACK, PARSJL, PRSRCN, PRSRSL
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Modified by C. F. Fischer for block calculation      22 May 1997   *
!   Updated by Xinghong He                               08 Jul 1998   *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def,   only: max_orbitals
      use debug_c
      use def_c
      use orb_c, ncfblock => csf_count
      use stat_c
      use terms_c,         only: jtab, ntab
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use prsrcn_i
      use parsjl_i
      use pack_i
      use convrt_i
      use iq_i
      use jqs_i
      use jcup_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nfile
      integer, intent(in) :: ncore
      integer, intent(in) :: jb
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: loadall = -119
      character(len=7), parameter :: myname = 'LODCSH2'
      integer, parameter :: nw2 = 2*max_orbitals
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iocc(max_orbitals)
      integer :: iqsub(nw2)
      integer :: jx(max_orbitals)
      integer :: ncorp1, nrec, ncf, npeel, i, j
      integer :: ios, ierr, loc, nqs, isparc, njx, ioc, ipty
      integer :: nqsn, njxn, npeeln, nopen, jlast, ilast, iocci
      integer :: nkji, ifulli, nu, jsub, iqt, nbeg, nend
      integer :: lenth, jxn, jpi, ncorel, iqgg
      logical :: empty, full
      character :: str*256, recl
!-----------------------------------------------
!
      if (jb /= loadall) then
         write (6, *) 'Loading CSF File for block ', jb
      else
         write (6, *) 'Loading CSF File for ALL blocks '
      endif

      ncorp1 = ncore + 1
      npeel = orbital_count - ncore
!
! NPEEL is used as 1) number of peel orbitals (here) and
!                  2) number of peel electrons (later in this routine)
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
!
!   Zero out the arrays that store packed integers - only when ncfblock>0
!
      csf_subshell_occupations(:max_orbitals,:ncfblock) = 0
      csf_subshell_quantum_codes(:max_orbitals,1,:ncfblock) = 0
      csf_subshell_quantum_codes(:max_orbitals,2,:ncfblock) = 0
      csf_subshell_quantum_codes(:max_orbitals,3,:ncfblock) = 0
      csf_coupling_quantum_codes(:max_orbitals,:ncfblock) = 0

      ncf = 0
    3 continue
      ncf = ncf + 1
!
      read (nfile, '(A)', iostat=ios) str

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This IF...READ makes the routine load the entire file (all blocks)
! by ignoring the end-of-block mark

      if (ios .eq. 0 .and. str(1:2) .eq. ' *' .and. jb .eq. loadall) &
        read (nfile, '(A)', iostat = ios) str
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (ios==0 .and. str(1:2)/=' *') then
!
!   Read in the occupations (q) of the peel shells; stop with a
!   message if an error occurs
!
         call prsrcn (str, ncore, iocc, ierr)
         if (ierr /= 0) go to 28
!
!   Read the J_sub and v quantum numbers
!
         read (nfile,'(A)',iostat = ios) str
         if (ios /= 0) then
            write (stderr_unit, *) myname//': Expecting subshell quantum', &
               ' number specification;'
            go to 27
         endif
         loc = len_trim(str)
         call parsjl (1, ncore, str, loc, iqsub, nqs, ierr)
         if (ierr /= 0) go to 27
!
!   Read the X, J, and (sign of) P quantum numbers
!
         read (nfile,'(A)',iostat = ios) str
         if (ios /= 0) then
            write (stderr_unit, *) myname//': Expecting intermediate ', &
               'and final angular momentum'
            write (stderr_unit, *) 'quantum number and final parity ', &
               'specification;'
            go to 26
         endif
!
!   Zero out the arrays that store packed integers
!
         csf_subshell_occupations(:max_orbitals,ncf) = 0
         csf_subshell_quantum_codes(:max_orbitals,1,ncf) = 0
         csf_subshell_quantum_codes(:max_orbitals,2,ncf) = 0
         csf_subshell_quantum_codes(:max_orbitals,3,ncf) = 0
         csf_coupling_quantum_codes(:max_orbitals,ncf) = 0

!   Determine the parity and all intermediate and the final
!   angular momentum quantum numbers
!
         loc = len_trim(str)
         recl = str(loc:loc)
         if (recl == '+') then
            isparc = 1
         else if (recl == '-') then
            isparc = -1
         else
            write (stderr_unit, *) myname//': Incorrect parity ', &
                            'specification;'
            go to 26
         endif
         loc = loc - 1
!
         call parsjl (2, ncore, str, loc, jx, njx, ierr)
         if (ierr /= 0) go to 26
!
!   Set the occupation and subshell quantum number array elements
!   in IQ, JQS for the core subshells
!
         do i = 1, ncore
            call pack (orbital_two_j(i) + 1, i, csf_subshell_occupations(1:max_orbitals,ncf))
            call pack (0, i, csf_subshell_quantum_codes(1:max_orbitals,1,ncf))
            call pack (0, i, csf_subshell_quantum_codes(1:max_orbitals,2,ncf))
            call pack (1, i, csf_subshell_quantum_codes(1:max_orbitals,3,ncf))
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
                     write (stderr_unit, *) myname//': Too few subshell ', &
                        'quantum numbers specified;'
                     go to 26
                  endif
                  nu = 0
                  jsub = iqsub(nqsn)
               else
                  if (iocci /= 4) then
                     nqsn = nqsn + 1
                     if (nqsn > nqs) then
                        write (stderr_unit, *) myname//': Too few subshell ', &
                           'quantum numbers specified;'
                        go to 26
                     endif
                     nu = 0
                     jsub = iqsub(nqsn)
                  else
                     nqsn = nqsn + 1
                     if (nqsn > nqs) then
                        write (stderr_unit, *) myname//': Too few subshell ', &
                           'quantum numbers specified;'
                        go to 26
                     endif
                     jsub = iqsub(nqsn)
                     if (jsub==4 .or. jsub==8) then
                        nu = jsub/2
                        nqsn = nqsn + 1
                        if (nqsn > nqs) then
                           write (stderr_unit, *) myname//': Too few subshell', &
                              ' quantum numbers specified;'
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
               call convrt (orbital_principal_n(i), str, lenth)
               write (stderr_unit, *) myname//': Subshell quantum numbers ', &
                  'specified incorrectly for '//str(1:lenth)//orbital_symmetry_labels(i)//&
                  ' subshell.'
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
                        write (stderr_unit, *) myname//': Too few intermediate', &
                           ' and final angular momentum', &
                           ' quantum numbers specified;'
                        go to 26
                     endif
                     jxn = jx(njxn)
                     do j = abs(jlast - jsub), jlast + jsub, 2
                        if (jxn == j) go to 11
                     end do
                     call convrt (orbital_principal_n(i), str, lenth)
                     write (stderr_unit, *) myname//': coupling of '//str(1:lenth)//&
                        orbital_symmetry_labels(i), ' subshell to previous subshells is incorrect.'
                     go to 26
                  endif
   11             continue
                  call pack (jxn + 1, nopen - 1, csf_coupling_quantum_codes(1:max_orbitals,ncf))
                  jlast = jxn
               else
                  jlast = jsub
               endif
            endif
            call pack(iocci, i, csf_subshell_occupations(1:max_orbitals,ncf))
            call pack (nu, i, csf_subshell_quantum_codes(1:max_orbitals,1,ncf))
            call pack (0, i, csf_subshell_quantum_codes(1:max_orbitals,2,ncf))
            call pack (jsub + 1, i, csf_subshell_quantum_codes(1:max_orbitals,3,ncf))
         end do
!
         do i = max(1,nopen), orbital_count
            call pack (0, i, csf_coupling_quantum_codes(1:max_orbitals,ncf))
         end do
!
         if (nqsn /= nqs) then
            write (stderr_unit, *) myname//': Too many subshell', &
               ' quantum numbers specified;'
            go to 26
         endif
!
         if (ilast /= ioc) njxn = njxn + 1
         if (njxn /= njx) then
            write (stderr_unit, *) myname//': Too many intermediate', &
               ' and final angular momentum', ' quantum numbers specified;'
            go to 26
         endif
!
         if (jx(njxn) /= jlast) then
            write (stderr_unit, *) myname//': Final angular momentum', &
               ' incorrectly specified;'
            go to 26
         endif
!
         ipty = (-1)**ipty
         if (ipty /= isparc) then
            write (stderr_unit, *) myname//': Parity specified incorrectly;'
            go to 26
         endif
!
         jpi = (jlast + 1)*ipty
         call pack (jpi, max_orbitals, csf_coupling_quantum_codes(1:max_orbitals,ncf))
!
         if (ncf > 1) then
            if (npeeln /= npeel) then
               write (stderr_unit, *) myname//': Inconsistency in the number', &
                  ' of electrons.'
               go to 26
            endif
         else
            npeel = npeeln
         endif
!
!   Check if this CSF was already in the list; stop with a
!   message if this is the case
!
         if (ncf > 1) then
            do j = 1, ncf - 1
               do i = ncorp1, orbital_count
                  if (iq(i,j) /= iq(i,ncf)) go to 17
                  if (jqs(1,i,j) /= jqs(1,i,ncf)) go to 17
                  if (jqs(2,i,j) /= jqs(2,i,ncf)) go to 17
                  if (jqs(3,i,j) /= jqs(3,i,ncf)) go to 17
               end do
               do i = 1, nopen - 1
                  if (jcup(i,j) /= jcup(i,ncf)) go to 17
               end do
            end do
            write (stderr_unit, *) myname//': Repeated CSF;'
            go to 26
         endif
!
!   Successfully read a CSF; update NREC and read another CSF
!
   17    continue
         nrec = nrec + 3

         go to 3
!
      else                                       ! the record just read is either ' *' or EOF, marking
            ! the end of a block or end of the file
!
!   There is always at least one CSF
!
         if (ncf == 1) then
            do i = 1, ncore
               call pack (orbital_two_j(i) + 1, i, csf_subshell_occupations(1:max_orbitals,1))
               call pack (0, i, csf_subshell_quantum_codes(1:max_orbitals,1,1))
               call pack (0, i, csf_subshell_quantum_codes(1:max_orbitals,2,1))
               call pack (1, i, csf_subshell_quantum_codes(1:max_orbitals,3,1))
            end do
            call pack (0, 1, csf_coupling_quantum_codes(1:max_orbitals,1))
            call pack (1, max_orbitals, csf_coupling_quantum_codes(1:max_orbitals,1))
         else
            ncf = ncf - 1
         endif
!
      endif

      if (ncf /= ncfblock) then
         write (stderr_unit, *) myname//': ncf=', ncf, 'ncfblock=', ncfblock
         stop
      endif
!
!   Check if any subshell is empty; eliminate it from the
!   list if this is the case; issue a message
!
      i = ncorp1
   19 continue
      if (i <= orbital_count) then
         do j = 1, ncf
            if (iq(i,j) /= 0) go to 23
         end do
         call convrt (orbital_principal_n(i), str, lenth)
         write (6, *) 'Subshell '//str(1:lenth)//orbital_symmetry_labels(i)//' is empty', &
            ' in all CSFs'
   23    continue
         i = i + 1
         go to 19
      endif
!
!   Store the number of electrons in the COMMON variable
!   This will act as a check now - it's been determined in lodcsh
!
      ncorel = 0
      ncorel = sum(orbital_two_j(:ncore)+1)
!      NELEC = NCOREL+NPEEL
      if (ncorel + npeel /= electron_count) then
         write (stderr_unit, *) myname//': nelec not equal to that in lodcsh'
         stop
      endif
      write (6,*)'There are ',ncf,' relativistic CSFs... load complete;'
      return
!
   26 continue
      backspace (nfile)
   27 continue
      backspace (nfile)
   28 continue
      backspace (nfile)
      write (stderr_unit, *) ' CSF sequence number: ', ncf
      do i = 1, 3
         read  (nfile,'(A)',err = 29,end = 29) str
         write (stderr_unit, *) str(1:len_trim(str))
      end do
   29 continue
      close(nfile)

      stop
      end subroutine lodcsh2
end module lodcsh2_i
