module lodcsh2_i
    implicit none (type, external)
    private
    public :: lodcsh2
contains
!***********************************************************************
!                                                                      *
      subroutine lodcsh2(nfile, ncore, jb)
      use, intrinsic :: iso_fortran_env, only: output_unit
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
      use grasp_parameters,   only: max_orbitals
      use debug_output_settings, only: cutoff => matrix_element_screening_cutoff, angular_coefficient_debug_flag, ibug2, ibug3, &
          & ibug4, ibug5, ibug6, ldbpa, &
          ldbpg, rmcdhf_debug_flags
      use atomic_calculation_state, only: electron_count
      use orbital_csf_basis, only: csf_count, csf_subshell_occupations, max_orbitals, orbital_count, orbital_l, &
          orbital_principal_n, orbital_symmetry_labels, orbital_two_j
      use csf_coupling_data, only: csf_coupling_quantum_codes, csf_subshell_quantum_codes, left_coupling_codes, &
          left_subshell_quantum_codes
      use subshell_term_table,         only: state_offsets, state_values
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use prsrcn_i, only: prsrcn
      use parsjl_i, only: parsjl
      use csf_byte_packing, only: pack_csf_subshell_byte
      use integer_format, only: format_integer
      use csf_quantum_queries, only: csf_occupation, csf_subshell_quantum_code, csf_intermediate_coupling
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
      character(len=7), parameter :: myname = "LODCSH2"
      integer, parameter :: nw2 = 2*max_orbitals
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iocc(max_orbitals)
      integer :: iqsub(nw2)
      integer :: jx(max_orbitals)
      integer :: ncorp1, nrec, loaded_csf_count, npeel, i, j
      integer :: ios, ierr, loc, nqs, isparc, njx, ioc, ipty
      integer :: nqsn, njxn, npeeln, nopen, jlast, ilast, iocci
      integer :: nkji, ifulli, nu, jsub, iqt, nbeg, nend
      integer :: lenth, jxn, jpi, ncorel, iqgg
      logical :: empty, full
      character :: str*256, recl
!-----------------------------------------------
!
      if (jb /= loadall) then
         write (output_unit, *) "Loading CSF File for block ", jb
      else
         write (output_unit, *) "Loading CSF File for ALL blocks "
      end if

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
      csf_subshell_occupations(:max_orbitals,:csf_count) = 0
      csf_subshell_quantum_codes(:max_orbitals,1,:csf_count) = 0
      csf_subshell_quantum_codes(:max_orbitals,2,:csf_count) = 0
      csf_subshell_quantum_codes(:max_orbitals,3,:csf_count) = 0
      csf_coupling_quantum_codes(:max_orbitals,:csf_count) = 0

      loaded_csf_count = 0
    3 continue
      loaded_csf_count = loaded_csf_count + 1
!
      read (nfile, "(A)", iostat=ios) str

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This IF...READ makes the routine load the entire file (all blocks)
! by ignoring the end-of-block mark

      if (ios == 0 .and. str(1:2) == " *" .and. jb == loadall) &
        read (nfile, "(A)", iostat = ios) str
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (ios==0 .and. str(1:2)/=" *") then
!
!   Read in the occupations (q) of the peel shells; stop with a
!   message if an error occurs
!
         call prsrcn (str, ncore, iocc, ierr)
         if (ierr /= 0) goto 28
!
!   Read the J_sub and v quantum numbers
!
         read (nfile,"(A)",iostat = ios) str
         if (ios /= 0) then
            write (stderr_unit, *) myname//": Expecting subshell quantum", &
               " number specification;"
            goto 27
         end if
         loc = len_trim(str)
         call parsjl (1, ncore, str, loc, iqsub, nqs, ierr)
         if (ierr /= 0) goto 27
!
!   Read the X, J, and (sign of) P quantum numbers
!
         read (nfile,"(A)",iostat = ios) str
         if (ios /= 0) then
            write (stderr_unit, *) myname//": Expecting intermediate ", &
               "and final angular momentum"
            write (stderr_unit, *) "quantum number and final parity ", &
               "specification;"
            goto 26
         end if
!
!   Zero out the arrays that store packed integers
!
         csf_subshell_occupations(:max_orbitals,loaded_csf_count) = 0
         csf_subshell_quantum_codes(:max_orbitals,1,loaded_csf_count) = 0
         csf_subshell_quantum_codes(:max_orbitals,2,loaded_csf_count) = 0
         csf_subshell_quantum_codes(:max_orbitals,3,loaded_csf_count) = 0
         csf_coupling_quantum_codes(:max_orbitals,loaded_csf_count) = 0

!   Determine the parity and all intermediate and the final
!   angular momentum quantum numbers
!
         loc = len_trim(str)
         recl = str(loc:loc)
         if (recl == "+") then
            isparc = 1
         else if (recl == "-") then
            isparc = -1
         else
            write (stderr_unit, *) myname//": Incorrect parity ", &
                            "specification;"
            goto 26
         end if
         loc = loc - 1
!
         call parsjl (2, ncore, str, loc, jx, njx, ierr)
         if (ierr /= 0) goto 26
!
!   Set the occupation and subshell quantum number array elements
!   in IQ, JQS for the core subshells
!
         do i = 1, ncore
            call pack_csf_subshell_byte (orbital_two_j(i) + 1, i, &
                       csf_subshell_occupations(1:max_orbitals,loaded_csf_count))
            call pack_csf_subshell_byte (0, i, csf_subshell_quantum_codes(1:max_orbitals,1,loaded_csf_count))
            call pack_csf_subshell_byte (0, i, csf_subshell_quantum_codes(1:max_orbitals,2,loaded_csf_count))
            call pack_csf_subshell_byte (1, i, csf_subshell_quantum_codes(1:max_orbitals,3,loaded_csf_count))
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
                     write (stderr_unit, *) myname//": Too few subshell ", &
                        "quantum numbers specified;"
                     goto 26
                  end if
                  nu = 0
                  jsub = iqsub(nqsn)
               else
                  if (iocci /= 4) then
                     nqsn = nqsn + 1
                     if (nqsn > nqs) then
                        write (stderr_unit, *) myname//": Too few subshell ", &
                           "quantum numbers specified;"
                        goto 26
                     end if
                     nu = 0
                     jsub = iqsub(nqsn)
                  else
                     nqsn = nqsn + 1
                     if (nqsn > nqs) then
                        write (stderr_unit, *) myname//": Too few subshell ", &
                           "quantum numbers specified;"
                        goto 26
                     end if
                     jsub = iqsub(nqsn)
                     if (jsub==4 .or. jsub==8) then
                        nu = jsub/2
                        nqsn = nqsn + 1
                        if (nqsn > nqs) then
                           write (stderr_unit, *) myname//": Too few subshell", &
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
               nbeg = state_offsets(loc+1) + 1
               nend = state_offsets(loc+2)
               do j = nbeg, nend, 3
                  if (state_values(j+2) /= jsub + 1) cycle
                  if (nu == 0) then
                     nu = state_values(j)
                     goto 9
                  else
                     if (state_values(j) == nu) goto 9
                  end if
               end do
               call format_integer (orbital_principal_n(i), str, lenth)
               write (stderr_unit, *) myname//": Subshell quantum numbers ", &
                  "specified incorrectly for "//str(1:lenth)//orbital_symmetry_labels(i)//&
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
                        write (stderr_unit, *) myname//": Too few intermediate", &
                           " and final angular momentum", &
                           " quantum numbers specified;"
                        goto 26
                     end if
                     jxn = jx(njxn)
                     do j = abs(jlast - jsub), jlast + jsub, 2
                        if (jxn == j) goto 11
                     end do
                     call format_integer (orbital_principal_n(i), str, lenth)
                     write (stderr_unit, *) myname//": coupling of "//str(1:lenth)//&
                        orbital_symmetry_labels(i), " subshell to previous subshells is incorrect."
                     goto 26
                  end if
   11             continue
                  call pack_csf_subshell_byte (jxn + 1, nopen - 1, csf_coupling_quantum_codes(1:max_orbitals,loaded_csf_count))
                  jlast = jxn
               else
                  jlast = jsub
               end if
            end if
            call pack_csf_subshell_byte(iocci, i, csf_subshell_occupations(1:max_orbitals,loaded_csf_count))
            call pack_csf_subshell_byte (nu, i, csf_subshell_quantum_codes(1:max_orbitals,1,loaded_csf_count))
            call pack_csf_subshell_byte (0, i, csf_subshell_quantum_codes(1:max_orbitals,2,loaded_csf_count))
            call pack_csf_subshell_byte (jsub + 1, i, csf_subshell_quantum_codes(1:max_orbitals,3,loaded_csf_count))
         end do
!
         do i = max(1,nopen), orbital_count
            call pack_csf_subshell_byte (0, i, csf_coupling_quantum_codes(1:max_orbitals,loaded_csf_count))
         end do
!
         if (nqsn /= nqs) then
            write (stderr_unit, *) myname//": Too many subshell", &
               " quantum numbers specified;"
            goto 26
         end if
!
         if (ilast /= ioc) njxn = njxn + 1
         if (njxn /= njx) then
            write (stderr_unit, *) myname//": Too many intermediate", &
               " and final angular momentum", " quantum numbers specified;"
            goto 26
         end if
!
         if (jx(njxn) /= jlast) then
            write (stderr_unit, *) myname//": Final angular momentum", &
               " incorrectly specified;"
            goto 26
         end if
!
         ipty = (-1)**ipty
         if (ipty /= isparc) then
            write (stderr_unit, *) myname//": Parity specified incorrectly;"
            goto 26
         end if
!
         jpi = (jlast + 1)*ipty
         call pack_csf_subshell_byte (jpi, max_orbitals, csf_coupling_quantum_codes(1:max_orbitals,loaded_csf_count))
!
         if (loaded_csf_count > 1) then
            if (npeeln /= npeel) then
               write (stderr_unit, *) myname//": Inconsistency in the number", &
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
         if (loaded_csf_count > 1) then
            do j = 1, loaded_csf_count - 1
               do i = ncorp1, orbital_count
                  if (csf_occupation(i,j) /= csf_occupation(i,loaded_csf_count)) goto 17
                  if (csf_subshell_quantum_code(1,i,j) /= csf_subshell_quantum_code(1,i,loaded_csf_count)) goto 17
                  if (csf_subshell_quantum_code(2,i,j) /= csf_subshell_quantum_code(2,i,loaded_csf_count)) goto 17
                  if (csf_subshell_quantum_code(3,i,j) /= csf_subshell_quantum_code(3,i,loaded_csf_count)) goto 17
               end do
               do i = 1, nopen - 1
                  if (csf_intermediate_coupling(i,j) /= csf_intermediate_coupling(i,loaded_csf_count)) goto 17
               end do
            end do
            write (stderr_unit, *) myname//": Repeated CSF;"
            goto 26
         end if
!
!   Successfully read a CSF; update NREC and read another CSF
!
   17    continue
         nrec = nrec + 3

         goto 3
!
      else                                       ! the record just read is either ' *' or EOF, marking
            ! the end of a block or end of the file
!
!   There is always at least one CSF
!
         if (loaded_csf_count == 1) then
            do i = 1, ncore
               call pack_csf_subshell_byte (orbital_two_j(i) + 1, i, csf_subshell_occupations(1:max_orbitals,1))
               call pack_csf_subshell_byte (0, i, csf_subshell_quantum_codes(1:max_orbitals,1,1))
               call pack_csf_subshell_byte (0, i, csf_subshell_quantum_codes(1:max_orbitals,2,1))
               call pack_csf_subshell_byte (1, i, csf_subshell_quantum_codes(1:max_orbitals,3,1))
            end do
            call pack_csf_subshell_byte (0, 1, csf_coupling_quantum_codes(1:max_orbitals,1))
            call pack_csf_subshell_byte (1, max_orbitals, csf_coupling_quantum_codes(1:max_orbitals,1))
         else
            loaded_csf_count = loaded_csf_count - 1
         end if
!
      end if

      if (loaded_csf_count /= csf_count) then
         write (stderr_unit, *) myname//": ncf=", loaded_csf_count, "ncfblock=", csf_count
         stop
      end if
!
!   Check if any subshell is empty; eliminate it from the
!   list if this is the case; issue a message
!
      i = ncorp1
   19 continue
      if (i <= orbital_count) then
         do j = 1, loaded_csf_count
            if (csf_occupation(i,j) /= 0) goto 23
         end do
         call format_integer (orbital_principal_n(i), str, lenth)
         write (output_unit, *) "Subshell "//str(1:lenth)//orbital_symmetry_labels(i)//" is empty", &
            " in all CSFs"
   23    continue
         i = i + 1
         goto 19
      end if
!
!   Store the number of electrons in the COMMON variable
!   This will act as a check now - it's been determined in lodcsh
!
      ncorel = 0
      ncorel = sum(orbital_two_j(:ncore)+1)
!      NELEC = NCOREL+NPEEL
      if (ncorel + npeel /= electron_count) then
         write (stderr_unit, *) myname//": nelec not equal to that in lodcsh"
         stop
      end if
      write (output_unit,*)"There are ",loaded_csf_count," relativistic CSFs... load complete;"
      return
!
   26 continue
      backspace (nfile)
   27 continue
      backspace (nfile)
   28 continue
      backspace (nfile)
      write (stderr_unit, *) " CSF sequence number: ", loaded_csf_count
      do i = 1, 3
         read  (nfile,"(A)",err = 29,end = 29) str
         write (stderr_unit, *) str(1:len_trim(str))
      end do
   29 continue
      close(nfile)

      stop
      end subroutine lodcsh2
end module lodcsh2_i
