!***********************************************************************
!                                                                      *
      subroutine sortmem(nfile,iggmax_k,ncoeff,ntgral,lprint,nb,fhead)
!                                                                      *
!   This routine sorts lists                                           *
!                                                                      *
!                     (ICLMN,INDEX,LABEL,COEFF)                        *
!   into lists                                                         *
!                     (LABEL,ICLMN,INDEX,COEFF)                        *
!                                                                      *
!   using Heapsort. File NFILE is closed by this routine.  NCOEFF is   *
!   the number of triads (INDEX, ...) and is an input. NTGRAL is the   *
!   number of different values of LABEL, and is an output. If LPRINT   *
!   is  .TRUE. , the contents of the sorted file are interpreted and   *
!   printed to unit 99.                                                *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 21 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:24:58   1/ 5/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: orbital_index_radix
      use memory_man
      use orb_c,           only: orbital_principal_n, csf_count, orbital_symmetry_labels
      use iounit_c
      use sacoef_c, coeff=>channel_coefficients,label=>channel_packed_orbital_labels,iclmn=>channel_csf_columns, &
         & index=>channel_matrix_positions
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nfile
      integer, intent(in) :: ncoeff
      integer, intent(out)  :: ntgral
      integer, intent(in)   :: nb,iggmax_k
      logical, intent(in)   :: lprint
      character, intent(in) :: fhead*(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!      INTEGER, DIMENSION(:), pointer :: ICLMN, INDEX, LABEL, NSWAP
      integer, pointer :: nswap(:)
!      REAL(DOUBLE), DIMENSION(:), pointer :: COEFF
      integer :: lcnum, lck, ierr, mb, i, l, ir, icl, ind, lab, nsw, j,&
         last, ibeg, iend, ncontr, ia, ib, k, id, ic, igg
      real(real64) :: cof
      character :: cnum*20, srtlab*8, mcplab*3, ck*2
      character (len = len (fhead) + 3):: fullname
!-----------------------------------------------
!
      if(nfile .eq. 31) then
         igg = 1
      else
         igg = nfile - 30
      end if
!      IGG = IGGMAX_K*(IGG-1)
!   Message
!
      call convrt (ncoeff, cnum, lcnum)
      if (nfile > 31) then
         call convrt (nfile - 32, ck, lck)
         write (6, *) 'Sorting '//cnum(1:lcnum)//' V(k='//ck(1:lck)//&
            ') coefficients ...', nfile
      else
         write (6, *) 'Sorting '//cnum(1:lcnum)//' T coefficients ...', nfile
      endif
!
      call convrt (nfile, ck, lck)
      if (lck > 2) then
         write (stderr_unit, *) 'sort: nfile > 99; check fullname'
         stop
      endif

      fullname = fhead//'.'//ck(1:2)

      open(29,file=fullname,status='OLD',form='UNFORMATTED', &
         iostat=ierr, position='APPEND')
      if (ierr /= 0) then
         write (stderr_unit, *) ' Error when opening the file ', fullname
         stop
      endif
!
!   Sort the list
!
      if (ncoeff > 0) then
!
!   Allocate storage for all required arrays
!
         if (nfile == 33) call alloc (nswap, ncoeff,'NSWAP', 'SORTMEM')
      endif
!
!   Sort LABEL into ascending order using the heapsort algorithm;
!   move the associated members of COEFF and INDEX in the same
!   manner; the code below is adapted from Press et al.
!
      if (nfile==33 .and. ncoeff>1) then

         l = ncoeff/2 + 1
         ir = ncoeff
  234    continue
         if (l > 1) then
            l = l - 1
            cof = coeff(l,igg)
            icl = iclmn(l,igg)
            ind = index(l,igg)
            lab = label(l,igg)
            nsw = nswap(l)
         else
            cof = coeff(ir,igg)
            icl = iclmn(ir,igg)
            ind = index(ir,igg)
            lab = label(ir,igg)
            nsw = nswap(ir)
            coeff(ir,igg) = coeff(1,igg)
            iclmn(ir,igg) = iclmn(1,igg)
            index(ir,igg) = index(1,igg)
            label(ir,igg) = label(1,igg)
            nswap(ir) = nswap(1)
            ir = ir - 1
            if (ir == 1) then
               coeff(1,igg) = cof
               iclmn(1,igg) = icl
               index(1,igg) = ind
               label(1,igg) = lab
               nswap(1) = nsw
               go to 456
            endif
         endif
         i = l
         j = l + l
  345    continue
         if (j <= ir) then
            if (j < ir) then
               if (label(j,igg) < label(j+1,igg)) j = j + 1
            endif
            if (lab < label(j,igg)) then
               coeff(i,igg) = coeff(j,igg)
               iclmn(i,igg) = iclmn(j,igg)
               index(i,igg) = index(j,igg)
               label(i,igg) = label(j,igg)
               nswap(i) = nswap(j)
               i = j
               j = j + j
            else
               j = ir + 1
            endif
            go to 345
         endif
         coeff(i,igg) = cof
         iclmn(i,igg) = icl
         index(i,igg) = ind
         label(i,igg) = lab
         nswap(i) = nsw
         go to 234

      else if (nfile/=33 .and. ncoeff>1) then
!
!   Sort LABEL into ascending order using the heapsort algorithm;
!   move the associated members of COEFF and INDEX in the same
!   manner; the code below is adapted from Press et al.
!
         l = ncoeff/2 + 1
         ir = ncoeff
   92    continue
         if (l > 1) then
            l = l - 1
            cof = coeff(l,igg)
            icl = iclmn(l,igg)
            ind = index(l,igg)
            lab = label(l,igg)
         else
            cof = coeff(ir,igg)
            icl = iclmn(ir,igg)
            ind = index(ir,igg)
            lab = label(ir,igg)
            coeff(ir,igg) = coeff(1,igg)
            iclmn(ir,igg) = iclmn(1,igg)
            index(ir,igg) = index(1,igg)
            label(ir,igg) = label(1,igg)
            ir = ir - 1
            if (ir == 1) then
               coeff(1,igg) = cof
               iclmn(1,igg) = icl
               index(1,igg) = ind
               label(1,igg) = lab
               go to 456
            endif
         endif
         i = l
         j = l + l
   93    continue
         if (j <= ir) then
            if (j < ir) then
               if (label(j,igg) < label(j+1,igg)) j = j + 1
            endif
            if (lab < label(j,igg)) then
               coeff(i,igg) = coeff(j,igg)
               iclmn(i,igg) = iclmn(j,igg)
               index(i,igg) = index(j,igg)
               label(i,igg) = label(j,igg)
               i = j
               j = j + j
            else
               j = ir + 1
            endif
            go to 93
         endif
         coeff(i,igg) = cof
         iclmn(i,igg) = icl
         index(i,igg) = ind
         label(i,igg) = lab
         go to 92

      endif
!
!   Sorting complete; rewrite the file header
!
!
  456 continue
      write (29) 'MCP', nb, csf_count, ncoeff
!GG      WRITE (9999,*) 'MCP', NB, NCF, NCOEFF, NFILE
!
!   Write the sorted list to mcp.xx
!
      if (ncoeff > 0) then
!
         last = label(1,igg)
         ibeg = 1
         iend = 1
         ntgral = 1
!
         do i = 2, ncoeff
            if (label(i,igg) == last) then
               iend = iend + 1
            else
               write (29) last, iend - ibeg + 1
!GG               WRITE (9999,*) LAST, IEND - IBEG + 1
               write (29) (iclmn(j,igg),index(j,igg),coeff(j,igg),     &
                                                          j=ibeg,iend)
!cjb-GG              DO J = IBEG,IEND
!GG               WRITE (9999,'(2I12,E25.15)') ICLMN(J),INDEX(J),COEFF(J)
!cjb-GG              END DO

!              IF (NFILE.EQ.33) WRITE (20) (NSWAP(J),J = IBEG,IEND)
               ntgral = ntgral + 1
               last = label(i,igg)
               ibeg = iend + 1
               iend = ibeg
            endif
         end do
!
         if (ibeg <= ncoeff) then
            write (29) last, ncoeff - ibeg + 1
!GG            WRITE (9999,*) LAST, NCOEFF - IBEG + 1
            write (29) (iclmn(j,igg),index(j,igg),coeff(j,igg),       &
                                                       j=ibeg,ncoeff)
!GG          DO J=IBEG,NCOEFF
!GG            WRITE (9999,'(2I12,E25.15)') ICLMN(J),INDEX(J),COEFF(J)
!GG          END DO

!           IF (NFILE.EQ.33) WRITE (20) (NSWAP(J),J = IBEG,NCOEFF)
         endif
!
      else
!
         ntgral = 0
!
      endif
!
!     write the terminator record for this block
!
      write (29) 0, 0
!GG      WRITE (9999,*) 0, 0
      close(29)
!
!   Completion message
!
      !PRINT *, ' ... sort complete; ', ntgral, ' integrals;'
!
!   Debug printout
!
      if (lprint) then
         write (99, 300)
         write (6, 300)
         if (ncoeff > 0) then
!
            last = label(1,igg)
            ibeg = 1
            iend = 1
!
            do i = 2, ncoeff
               if (label(i,igg) == last) iend = iend + 1
  567          continue
               if (label(i,igg)==last .and. i/=ncoeff) cycle
               lab = last
               ncontr = iend - ibeg + 1
               if (nfile == 31) then
                  ia = mod(lab,key)
                  ib = lab/key
                  write (99, 301) orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ib), &
                     & orbital_symmetry_labels(ib)
                  do j = ibeg, iend
                     write (99, 302) iclmn(j,igg),index(j,igg),       &
                                                        coeff(j,igg)
                  end do
               else
                  k = nfile - 32
                  id = mod(lab,key)
                  lab = lab/key
                  ib = mod(lab,key)
                  lab = lab/key
                  ic = mod(lab,key)
                  ia = lab/key
                  write (99, 304) k, orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ib), &
                     & orbital_symmetry_labels(ib), orbital_principal_n(ic), orbital_symmetry_labels&
                     (ic), orbital_principal_n(id), orbital_symmetry_labels(id)
                  do j = ibeg, iend
                     write (99, 305) k,iclmn(j,igg),index(j,igg),     &
                                                          coeff(j,igg)
                  end do
               endif
               last = label(i,igg)
               ibeg = iend + 1
               iend = ibeg
               if (iend == ncoeff) go to 567
            end do
         endif
         write (99, 303) ntgral
      endif
!
!   Deallocate storage
!
      if (ncoeff > 0) then
         if (nfile == 33) call dalloc (nswap, 'NSWAP', 'SORTMEM')
      endif

  300 format(/,'From SORT:')
  301 format(' I(',1i2,1a2,',',1i2,1a2,'):')
  302 format('  T_[',1i2,',',1i4,'] = ',1p,d19.12)
  303 format('  Number of integrals is ',1i4)
  304 format(' R^[(',1i2,')] (',1i2,1a2,',',1i2,1a2,';',1i2,1a2,',',1i2,1a2,&
         '):')
  305 format('  V^[(',1i2,')]_[',1i8,',',1i8,'] = ',1p,d19.12)

      return
      end subroutine sortmem
