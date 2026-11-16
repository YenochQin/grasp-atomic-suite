!***********************************************************************
!                                                                      *
      subroutine sort(nfile, ncoeff, ntgral, lprint, nb, fhead)
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
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: orbital_index_radix
      use memory_man
      use orb_c,           only: orbital_principal_n, csf_count, orbital_symmetry_labels
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: nfile
      integer  :: ncoeff
      integer, intent(out)  :: ntgral
      integer, intent(in)   :: nb
      logical, intent(in)   :: lprint
      character, intent(in) :: fhead*(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, dimension(:), pointer :: iclmn, index, label, nswap
      real(double), dimension(:), pointer :: coeff
      integer :: lcnum, lck, ierr, mb, i, l, ir, icl, ind, lab, nsw, j,&
         last, ibeg, iend, ncontr, ia, ib, k, id, ic
      real(double) :: cof
      character :: cnum*20, srtlab*8, mcplab*3, ck*2
      character (len = len (fhead) + 3):: fullname
!-----------------------------------------------
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
!
      fullname = fhead//'.'//ck(1:2)
      open(29,file=fullname,status='OLD',form='UNFORMATTED', &
         iostat=ierr, position='APPEND')
      if (ierr /= 0) then
         write (stderr_unit, *) ' Error when opening the file ', fullname
         stop
      endif
!
      read (nfile) mcplab, mb
      if (nb /= mb) then
         write (stderr_unit, *) 'sort: nb = ', nb, '.NE. mb (=', mb, ')'
         stop
      endif
!
!   Sort the list
!
      if (ncoeff > 0) then
!
!   Allocate storage for all required arrays
!
         call alloc (coeff, ncoeff, 'COEFF', 'SORT')
         call alloc (iclmn, ncoeff, 'ICLMN', 'SORT')
         call alloc (index, ncoeff, 'INDEX', 'SORT')
         call alloc (label, ncoeff, 'LABEL', 'SORT')
         if (nfile == 33) call alloc (nswap, ncoeff,'NSWAP', 'SORT'  )
!
!   Read arrays into memory from NFILE
!
         do i = 1, ncoeff
            read (nfile) iclmn(i), index(i), label(i), coeff(i)
         end do
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
            cof = coeff(l)
            icl = iclmn(l)
            ind = index(l)
            lab = label(l)
            nsw = nswap(l)
         else
            cof = coeff(ir)
            icl = iclmn(ir)
            ind = index(ir)
            lab = label(ir)
            nsw = nswap(ir)
            coeff(ir) = coeff(1)
            iclmn(ir) = iclmn(1)
            index(ir) = index(1)
            label(ir) = label(1)
            nswap(ir) = nswap(1)
            ir = ir - 1
            if (ir == 1) then
               coeff(1) = cof
               iclmn(1) = icl
               index(1) = ind
               label(1) = lab
               nswap(1) = nsw
               go to 456
            endif
         endif
         i = l
         j = l + l
  345    continue
         if (j <= ir) then
            if (j < ir) then
               if (label(j) < label(j+1)) j = j + 1
            endif
            if (lab < label(j)) then
               coeff(i) = coeff(j)
               iclmn(i) = iclmn(j)
               index(i) = index(j)
               label(i) = label(j)
               nswap(i) = nswap(j)
               i = j
               j = j + j
            else
               j = ir + 1
            endif
            go to 345
         endif
         coeff(i) = cof
         iclmn(i) = icl
         index(i) = ind
         label(i) = lab
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
            cof = coeff(l)
            icl = iclmn(l)
            ind = index(l)
            lab = label(l)
         else
            cof = coeff(ir)
            icl = iclmn(ir)
            ind = index(ir)
            lab = label(ir)
            coeff(ir) = coeff(1)
            iclmn(ir) = iclmn(1)
            index(ir) = index(1)
            label(ir) = label(1)
            ir = ir - 1
            if (ir == 1) then
               coeff(1) = cof
               iclmn(1) = icl
               index(1) = ind
               label(1) = lab
               go to 456
            endif
         endif
         i = l
         j = l + l
   93    continue
         if (j <= ir) then
            if (j < ir) then
               if (label(j) < label(j+1)) j = j + 1
            endif
            if (lab < label(j)) then
               coeff(i) = coeff(j)
               iclmn(i) = iclmn(j)
               index(i) = index(j)
               label(i) = label(j)
               i = j
               j = j + j
            else
               j = ir + 1
            endif
            go to 93
         endif
         coeff(i) = cof
         iclmn(i) = icl
         index(i) = ind
         label(i) = lab
         go to 92

      endif
!
!   Sorting complete; rewrite the file header
!
  456 continue
      write (29) 'MCP', nb, csf_count, ncoeff
!GG      WRITE (9999,*) 'MCP', NB, NCF, NCOEFF, NFILE
!
!   Write the sorted list to mcp.xx
!
      if (ncoeff > 0) then
!
         last = label(1)
         ibeg = 1
         iend = 1
         ntgral = 1
!
         do i = 2, ncoeff
            if (label(i) == last) then
               iend = iend + 1
            else
               write (29) last, iend - ibeg + 1
!GG               WRITE (9999,*) LAST, IEND - IBEG + 1
               write (29) (iclmn(j),index(j),coeff(j),j=ibeg,iend)
!cjb-GG       DO J = IBEG,IEND
!GG               WRITE (9999,'(2I12,E25.15)') ICLMN(J),INDEX(J),COEFF(J)
!cjb-GG        END DO

               ntgral = ntgral + 1
               last = label(i)
               ibeg = iend + 1
               iend = ibeg
            endif
         end do
!
         if (ibeg <= ncoeff) then
            write (29) last, ncoeff - ibeg + 1
!GG            WRITE (9999,*) LAST, NCOEFF - IBEG + 1
            write (29) (iclmn(j),index(j),coeff(j),j=ibeg,ncoeff)
          do j=ibeg,ncoeff
!GG            WRITE (9999,'(2I12,E25.15)') ICLMN(J),INDEX(J),COEFF(J)
          end do

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
!      PRINT *, ' ... sort complete; ', ntgral, ' integrals;'
!
!   Debug printout
!
      if (lprint) then
         write (99, 300)
         write (6, 300)
         if (ncoeff > 0) then
!
            last = label(1)
            ibeg = 1
            iend = 1
!
            do i = 2, ncoeff
               if (label(i) == last) iend = iend + 1
  567          continue
               if (label(i)==last .and. i/=ncoeff) cycle
               lab = last
               ncontr = iend - ibeg + 1
               if (nfile == 31) then
                  ia = mod(lab,key)
                  ib = lab/key
                  write (99, 301) orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ib), &
                     & orbital_symmetry_labels(ib)
                  do j = ibeg, iend
                     write (99, 302) iclmn(j), index(j), coeff(j)
                  end do
               else
                  k = nfile - 32
                  id = mod(lab,key)
                  lab = lab/key
                  ib = mod(lab,key)
                  lab = lab/key
                  ic = mod(lab,key)
                  ia = lab/key
                  write (99, 304)k,orbital_principal_n(ia),orbital_symmetry_labels(ia),orbital_principal_n(ib), &
                     & orbital_symmetry_labels(ib),orbital_principal_n(ic), &
                                 orbital_symmetry_labels(ic),orbital_principal_n(id),orbital_symmetry_labels(id)
                  do j = ibeg, iend
                     write (99, 305)k,iclmn(j),index(j),coeff(j)
                  end do
               endif
               last = label(i)
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
         call dalloc (coeff,'COEFF', 'SORT')
         call dalloc (iclmn, 'ICLMN', 'SORT')
         call dalloc (index, 'INDEX', 'SORT')
         call dalloc (label, 'LABEL', 'SORT')
         if (nfile == 33) call dalloc (nswap, 'NSWAP', 'SORT')
      endif

  300 format(/,'From SORT:')
  301 format(' I(',1i2,1a2,',',1i2,1a2,'):')
  302 format('  T_[',1i2,',',1i4,'] = ',1p,d19.12)
  303 format('  Number of integrals is ',1i4)
  304 format(' R^[(',1i2,')] (',1i2,1a2,',',1i2,1a2,';',1i2,1a2,',',1i2,1a2,&
         '):')
  305 format('  V^[(',1i2,')]_[',1i8,',',1i8,'] = ',1p,d19.12)

      return
      end subroutine sort
