!***********************************************************************
!                                                                      *
      subroutine setham (myid, nprocs, jblock, elsto,icstrt, nelmntt, &
                         atwinv,slf_en)
!                                                                      *
!   Sets up the Hamiltonian matrix and determines the average energy.  *
!
!   Serial I/O moved out; able to run on single/multiple processors
!   For this purpose a new common /setham_to_genmat2/ is created
!                                                                      *
!   Call(s) to: [LIB92]: ALCBUF, CONVRT, DALLOC, ICHOP, RKCO, TNSRJJ.  *
!               [RCI92]: BRINT, IABINT, KEINT, RKINTC, VINT, VPINT.    *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 30 Oct 1992   *
!   Block version by Xinghong He          Last revision: 15 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double, long
      use parameter_def,   only: nnnp, keyorb
      use memory_man
!cjb iso_fortran_env :: output_unit, error_unit
      use, intrinsic :: iso_fortran_env
!-----------------------------------------------
!   C O M M O N    B l o c k s
!-----------------------------------------------
      use bcore_c, only: icore
      use bilst_c
      use buffer_c
      use coeils_c
      use decide_c
      use debug_c
      use def_c
      use eigv_c
      use iccu_c
      use grid_c
      use hmat_c
      use keilst_c
      use ncdist_c
      use orb_c, only: ncf,  nw
      use prnt_c
      use stat_c
      use stor_c
      use tatb_c
      use vinlst_c
      use vpilst_c
      use wave_c
      use cteilsrk_c
      use blim_c
      use where_c
      use eigvec1_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use alcbuf_i
      use iabint_i
      use brint1_i
      use brint2_i
      use brint3_i
      use brint4_i
      use brint5_i
      use brint6_i
      use convrt_i
      use ichop_i
      use keint_i
      use rkintc_i
      use vint_i
      use vpint_i
      use onescalar_i, only: onescalar
      use rkco_gg_i, only: rkco_gg
      implicit none
      external breid,cord
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer(long)       :: nelmntt
      integer             :: jblock, icstrt
      integer, intent(in) :: myid, nprocs
      real(double)        :: elsto, atwinv
      real(double), dimension(*) :: slf_en
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      real(double), dimension(nnnw) :: tshell
      real(double) :: tgrl1, tgrl2, tegral

      integer, parameter :: key = keyorb
!
!     Matrix elements smaller than CUTOFF are not accumulated
!
!cjb cutoff is use associated and cannot be redeclared
!     REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-20
!cjb CUTOFF = 1.0D-12 below
!cjb
      integer :: ipi, ipj, inc1, inc2, kt, ipt, incor, ncoec, nctec, &
                 i, j, nmcbp, ncore, ic, nelc, irstart, ir, ia, ib,  &
                 itype, iia
       real(double) :: elemnt, precoeff, tcoeff, vcoeff, contr
!-----------------------------------------------------------------------
      nelmnt = nelmntt
!cjb
!cjb  CUTOFF = 1.0D-12
      cutoff = 1.0d-12

      atwinv = 1.d0/emn

      if (iprerun .eq. 2) then
         do ipi = 1,nvec
            do ipj = 1,ncf
               write (*,*) ipi,ipj,evec1(ipj+(ipi-1)*ncf)
            enddo
         enddo
      endif
!
!   Allocate storage to arrays in COMMON/BUFFER/; these are
!   used for the Coulomb and transverse two-electron integrals
!
      call alcbuf (1)

!     ...Locals
      call alloc (emt, ncf, 'EMT','SETHAM' )
      call alloc (irow, ncf, 'IROW', 'SETHAM')
!
      inc1 = 1
      inc2 = 1
!
!   Initialisations for contributions from the Dirac-Coulomb
!   operator
!
      kt  = 0
      ipt = 1
!
      incor = 1

      ncoec = 0
!
      nctec   = 0

      if (ltrans) then

!        ...Initialisations for transverse interaction correction
         do 2 i = 1, nw
            icore(i) = 0
            do j = 1, ncf
               if (ichop (i,j) .le. 0) goto 2
            enddo
            icore(i) = 1
    2    continue

         nmcbp = 0
         ncore = 0
      endif

! Loop over rows of the Hamiltonian matrix - distributed

      do 10 ic = icstrt, ncf, nprocs

         nelc = 0    ! counter - Number of non-zeros of this row

! Loop over columns of the current row

         irstart = 1
         do 85 ir = irstart, ic

! PER
            if (lfordr .and. (ir .gt. iccut(1))) then
               if (ir.ne.ic) cycle
            end if
! PER

            elemnt = 0.d0     ! accumulates various contributions to H

!
!   Generate the integral list for the matrix element of the
!   one-body operators
!
            if (iprerun .eq. 1) then
               inc1 = 0
               inc2 = 0
               if (ic.le.ncsfpre .or. ic.eq.ir) then
                  inc1 = 1
               endif
            endif

            if (iprerun .eq. 2) then
!
!   Diagonal elements are always included
!   Off diagonal elements are included only if the
!   products of the weights from the prerun are larger
!   than the cutoff.
!
               if (ic .eq. ir) then
                  inc1 = 1
                  inc2 = 1
               else
                  inc1 = 0
                  inc2 = 0
               endif
               do ipi = 1,nvec
                  precoeff =                                           &
                   dabs(evec1(ic+(ipi-1)*ncf)*evec1(ir+(ipi-1)*ncf))
                  if (precoeff .gt. coeffcut1) inc1 = 1
                  if (precoeff .gt. coeffcut2) inc2 = 1
               enddo
            endif

!            ...INC1.EQ.1 ------------>
         if (inc1 .eq. 1) then   !inc1 is always 1 without PRE-RUN
           call onescalar(ic,ir,ia,ib,tshell)
!
!   Accumulate the contribution from the one-body operators:
!   kinetic energy, electron-nucleus interaction; update the
!   angular integral counter
!
         if (ia .ne. 0) then
            if (ia .eq. ib) then
               do ia = 1,nw
                  iia = ia
                  tcoeff = tshell(ia)
                  if (abs (tcoeff) .gt. cutoff) then
                     ncoec = ncoec + 1
                     call iabint (iia, iia, tegral)
                        !------------------------
                     elemnt = elemnt + tegral*tcoeff
                     if (lnms) then
                        call keint (iia,iia,tegral)
                        !------------------------
                        elemnt = elemnt + tegral*atwinv*tcoeff
                     endif
                     if (lvp) then
                        call vpint (iia, iia, tegral)
                        !------------------------
                        elemnt = elemnt + tegral*tcoeff
                     endif
                  endif
               enddo
            else
               tcoeff = dble(tshell(1))
               if (dabs (tcoeff) .gt. cutoff) then
                  ncoec = ncoec + 1
                  call iabint (ia, ib, tegral)
                        !------------------------
                  elemnt = elemnt + tegral*tcoeff
                  if (lnms) then
                     call keint (ia, ib, tegral)
                        !------------------------
                     elemnt = elemnt + tegral*atwinv*tcoeff
                  endif
                  if (lvp) then
                     call vpint (ia, ib, tegral)
                        !------------------------
                     elemnt = elemnt + tegral*tcoeff
                  endif
               endif
            endif
         endif
!
         ibug1 = 0
!
!   Accumulate the contributions from the two-electron
!   Coulomb operator and the mass polarisation; the latter
!   is computed first because the orbital indices may be
!   permuted by RKINTC
!
         nvcoef = 0
!
         call rkco_gg (ic, ir, cord, incor, 1)
!
         do 7 i = 1, nvcoef
            vcoeff = dble(coeff(i))
            if (dabs (vcoeff) .gt. cutoff) then
               nctec = nctec + 1
               if (lsms) then
                  if (label(5,i) .eq. 1) then
                     call vint (label(1,i), label(3,i), tgrl1)
                     call vint (label(2,i), label(4,i), tgrl2)
                     elemnt = elemnt - tgrl1*tgrl2*atwinv*vcoeff
                  endif
               endif
               call rkintc (label(1,i), label(2,i),                  &
                            label(3,i), label(4,i),                  &
                            label(5,i), tegral)
               elemnt = elemnt + tegral*vcoeff
            endif
    7    continue
!
         ibug1 = 0

         endif  !inc1 is always 1 without PRE-RUN
!            ...INC1.EQ.1 <------------
!***********************************************************************
!            ...LTRANS .AND. (INC2.EQ.1) ------------>
         if (ltrans .and. (inc2.eq.1)) then
            !IF (INC2 .EQ. 1) THEN  !inc2 is always 1 without PRE-RUN
!
!   Accumulate the contribution from the two-electron
!   transverse interaction operator
!
            nvcoef = 0
!
            call rkco_gg (ic, ir, breid, 1, 2)
!
            do 8 i = 1, nvcoef
               if (dabs (coeff(i)) > cutoff) then
                  nmcbp = nmcbp + 1
                  itype = abs (label(6,i))
                  if (itype == 1) then
                     call brint1 (label(1,i), label(2,i),               &
                                  label(3,i), label(4,i),               &
                                  label(5,i), tegral)
                  elseif (itype == 2) then
                     call brint2 (label(1,i), label(2,i),               &
                                  label(3,i), label(4,i),               &
                                  label(5,i), tegral)
                  elseif (itype == 3) then
                     call brint3 (label(1,i), label(2,i),               &
                                  label(3,i), label(4,i),               &
                                  label(5,i), tegral)
                  elseif (itype == 4) then
                     call brint4 (label(1,i), label(2,i),               &
                                  label(3,i), label(4,i),               &
                                  label(5,i), tegral)
                  elseif (itype == 5) then
                     call brint5 (label(1,i), label(2,i),               &
                                  label(3,i), label(4,i),               &
                                  label(5,i), tegral)
                  elseif (itype == 6) then
                     call brint6 (label(1,i), label(2,i),               &
                                  label(3,i), label(4,i),               &
                                  label(5,i), tegral)
                  endif
                  contr = coeff(i)*tegral
                  if (label(6,i) > 0) then
                     elemnt = elemnt + contr
                  else
!                        ...It comes here only when ic=ir=1
!                           clue: rkco<-breid<-talk<-label(6,i)
                     ncore = ncore + 1
                     elsto = elsto + contr
                  endif
               endif
    8       continue
!
            ibug1 = 0
!
!               ...ELSTO is a constant over all diagonals, thus its
!                  contribution to the total energy can be added later
!            IF (IR .EQ. IC) ELEMNT = ELEMNT + ELSTO
!
            !ENDIF   !inc2 is always 1 without PRE-RUN
         endif
!            ...LTRANS .AND. (INC2.EQ.1) <------------
!***********************************************************************
!
! Store this element if it is diagonal or its value is greater than
! CUTOFF
!
         if ( (ir == ic) .or. (dabs (elemnt) > cutoff) ) then
            nelc       = nelc + 1
            emt(nelc)  = elemnt
            irow(nelc) = ir
         endif
!
   85    continue
! zou
!        print *, ic,SLF_EN(IC)
         if(lse) emt(nelc) = emt(nelc) + slf_en(ic)
! zou
!
!   This column is done; write it to disk
!
         write (imcdf) nelc, elsto, (emt(ir), ir = 1, nelc),       &
                                   (irow(ir), ir = 1, nelc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This EAV (and the above EMT) does not have ELSTO.
         eav = eav + emt(nelc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Output IC on screen to show how far the calculation has proceeded
!
!-----------------------------------------------------------------------
!cjb      if (mod(ic-1,100).eq.0) then
!cjb         PRINT '(A4,I10,A9,I10,A9,I5,A8,I5)', 'Row ', ic,       &
!cjb                 ' nnonz = ', nelc,                             &
!cjb                 ' block = ', jblock,' myid = ', myid
!cjb      ENDIF
!
!cjb          CALL CONVRT (IC,CNUM,LCNUM)
!cjb          if (mod(IC,100).eq.0) then
!cjb             PRINT *, 'Column '//CNUM(1:LCNUM)//' complete;'
!cjb          end if
!cjb
!
!cjb MPI progress begin --------------------------------------------------
!
!cjb just started
       if (ic .le. nprocs) then
            print '(A6,I10,A9,I10,A9,I5,A8,I5)', 'Start ', ic,       &
                    ' nnonz = ', nelc,                             &
                    ' block = ', jblock,' myid = ', myid
        flush output_unit
        flush error_unit
       endif
!
!cjb progress
       if (ic .gt. nprocs .and. ic .lt. ncf-nprocs .and. &
           mod (ic-1,100*nprocs) .eq. myid) then
         if (myid .eq. 0) then
            print '(A5,I10,A9,I10,A9,I5,A8,I5)', 'Done ', ic,       &
                    ' nnonz = ', nelc,                             &
                    ' block = ', jblock,' myid = ', myid
        flush output_unit
        flush error_unit
         else
            print '(A5,I10,A9,I10,A9,I5,A8,I5)', 'Row  ', ic,       &
                    ' nnonz = ', nelc,                             &
                    ' block = ', jblock,' myid = ', myid
        flush output_unit
        flush error_unit
         endif
       endif
!
!cjb almost finished
       if (ic .gt. ncf-nprocs) then

            print '(A7,I10,A9,I10,A9,I5,A8,I5)', 'Finish ', ic,       &
                    ' nnonz = ', nelc,                             &
                    ' block = ', jblock,' myid = ', myid
        flush output_unit
        flush error_unit
       endif
!cjb MPI progress end ----------------------------------------------------
!
!
!   Update the counter for the total number of elements
!
         nelmnt = nelmnt + nelc
!
   10 continue
!***********************************************************************
!
!   Deallocate storage for the arrays in /BUFFER/
!
      call alcbuf (3)

!     ...Locals
      call dalloc (emt, 'EMT', 'SETHAM')
      call dalloc (irow, 'IROW', 'SETHAM')

!  Fill the common block /setham_to_genmat2/ for use in genmat2

      cutofftmp = cutoff
      ncoeitmp = ncoei
      ncoectmp = ncoec
      ncteitmp = rk_local_count(size(valteirk), myid, nprocs)
      nctectmp = nctec
      ntpitmp = ntpi
      nmcbptmp = nmcbp
      ncoretmp = ncore
      nvpitmp = nvpi
      nkeitmp = nkei
      nvintitmp = nvinti
      nelmnttmp = nelmnt
      ncftmp = ncf

      return
      end subroutine setham
