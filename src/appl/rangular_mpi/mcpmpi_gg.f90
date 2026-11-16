!***********************************************************************
!                                                                      *
!cjb  myid, nprocs = NOT args
      subroutine mcpmpi (nb, restrt, fhead)
!                                                                      *
!   This routine controls the computation  and storage of the values   *
!   and all indices of the angular coefficients                        *
!                                                                      *
!                                       k                              *
!                   T  (ab)            V  (abcd)                       *
!                    rs                 rs                             *
!                                                                      *
!   k is the multipolarity of a two-particle Coulomb integral. a, b,   *
!   c and d are orbital sequence numbers.  r and s are configuration   *
!   state function indices.                                            *
!                                                                      *
!   Call(s) to: [LIB92]: ALCBUF, ALLOC, CONVRT, DALLOC, RKCO_GG,       *
!                        TNSRJJ.                                       *
!               [GENMCP]: FNDBEG, SETSDA, SORT.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 28 Sep 1993   *
!   Modified by C. Froese Fischer for block computation.               *
!   Modified by G. Gaigalas and J. Bieron for new spin-angular         *
!   integration                                        01 April 2012   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use iso_fortran_env, only: int64
      use rangular_profile, only: stage_clock, record_stage
      use parameter_def,   only:  nnnw, keyorb
      use memory_man
      use mpi_c,           only:  myid, nprocs
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use outsdampi_i
      use fndbeg_i
      use alcbuf_i
      use rkco_gg_i
      use setsda_i
      use sort_i
!-----------------------------------------------
!   C O M M O N  B L O C K S
!-----------------------------------------------
      use buffer_c,  only: nvcoef, label, coeff
      use debug_c,   only: ldbpa
      use default_c, only: ndef
      use iccu_c,    only: iccut
      use mcp_c,     only: kmax, diag, lfordr
      use orb_c,     only: ncf, iqa
      use stat_c,    only: jqsa, jcupa
      use onescalar_i, only: onescalar
      implicit none
      external cord
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: nb
!     INTEGER  :: MYID
!     INTEGER  :: NPROCS
      logical , intent(in) :: restrt
      character(len=*), intent(in) :: fhead
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = keyorb
      integer, parameter :: keysq = keyorb*keyorb
!GG      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-10
!GG      REAL(DOUBLE), PARAMETER :: CUTOFF = 0.0D0
!cjb  REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-08
      real(double), parameter :: cutoff = 1.0d-10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, dimension(:), pointer :: llistv
      integer :: jastrt,jbstrt,npos,llistt,k,ja,jb,ia,ib,lab,i,ic, &
         id, nswap, iswap, lac, lbd, ntgi
      real(double), dimension(nnnw) :: tshell
      real(double) :: vcoeff
      logical :: f0int, lincr
      integer(int64) :: stage_start_count
!-----------------------------------------------
!
!cjb           PRINT *, ' mcpmpi_gg.f90 FHEAD '
!cjb           PRINT *, ' mcpmpi_gg.f90 FHEAD = ', FHEAD
!
!   Allocate storage to the array that stores the lengths
!   of the lists of V coefficients of each multipolarity
!
      stage_start_count = stage_clock()
      call alloc(llistv, 0, kmax, 'LLISTV', 'MCPMPI')
!
!   If this is a restart, determine the pair of CSFs with which
!   the computation should begin.
!   If not restart, initialize the counters for the total number of
!   T coefficients and V coefficients of each multipolarity.
!
      if (restrt) then
         print *, 'Not ready for RESTART'
         stop
         call fndbeg (jastrt, jbstrt, npos, llistt, llistv)
      else
         jastrt = 1 + myid
         npos = 0
         llistt = 0
         do k = 0, kmax
            llistv(k) = 0
         enddo
      endif
!
!   Allocate storage for the arrays in BUFFER
!
      call alcbuf (1)
!
!   Write header to  .dbg  file if appropriate
!
      if (ldbpa(2) .or. ldbpa(3)) write (99,300)
!
!   JA and JB respectively refer to the initial and final states
!   in the list of NCF configurations
!
      call record_stage('coefficient_buffer_allocation',nb,stage_start_count)
      stage_start_count = stage_clock()
      do 5 ja = jastrt, ncf, nprocs
         jbstrt = 1
         do 4 jb = jbstrt, ja
            if (diag .or. (lfordr .and. (jb .gt. iccut(nb))) ) then
               if (jb /= ja) cycle
            end if
!
!   LINCR is .TRUE. if npos is to be incremented by 1; there
!   is always a diagonal element in each column
!
            if (jb /= ja) then
               lincr = .true.
            else
               npos = npos + 1
               lincr = .false.
            endif
            if (jb /= ja) then
!   Compute T coefficients
               call onescalar(ja,jb,ia,ib,tshell)
               if (ia /= 0 .and. ia /= ib .and.       &
                      abs (tshell(1)) > cutoff) then
                  npos = npos + 1
                  lincr = .false.
                  llistt = llistt + 1
                  lab = min (ia,ib) * key + max (ia,ib)
                  write (31) ja, npos, lab, tshell(1)
               endif
            endif
!
!   Call the MCP package to generate V coefficients; ac and bd
!   are the density pairs. NVCOEF is initialized here but changed
!   in /rkco/cor[d] via COMMON.
!
            nvcoef = 0
!            CALL RKCO_GG (JA, JB, 0, 1)
            call rkco_gg (ja, jb, cord, 0, 1)
            do 2 i = 1, nvcoef
               vcoeff = coeff(i)
               if (abs (vcoeff) > cutoff) then
                  ia = label(1,i)
                  ib = label(2,i)
                  ic = label(3,i)
                  id = label(4,i)
                  k  = label(5,i)
                  f0int=(k.eq.0).and.(ia.eq.ic).and.(ib.eq.id)
                  if (.not. f0int) then
! Swap index to make sure IA <= IC, IB <= ID and record the number
! of swaps
                     nswap = 0
                     if (ia > ic) then
                        iswap = ic
                        ic = ia
                        ia = iswap
                        nswap = nswap + 1
                     endif
                     if (ib > id) then
                        iswap = id
                        id = ib
                        ib = iswap
                        nswap = nswap + 1
                     endif
                     if (lincr) then
                        npos = npos + 1
                        lincr = .false.
                     endif
                     llistv(k) = llistv(k) + 1
                     lac = ia * key + ic
                     lbd = ib * key + id
                     if (lac .lt. lbd) then
                        lab = lac * keysq + lbd
                     else
                        lab = lbd * keysq + lac
                     endif
                     write (32+k) ja, npos, lab, vcoeff
                  endif
               endif
    2       continue
!
!   All angular coefficients for this pair of CSFs have been
!   generated; update file 30
!
            if ((jb .eq. ja) .or. (.not. lincr))   &
               write (30) ja, jb, npos
    4    continue
         if (mod(ja-1,100) == 0) print '(A4,I7,A9,I10,A9,I3,A8,I3)',  &
            'Row ',ja,' nnonz = ',npos,' block = ',nb,' myid = ',myid
    5 continue
      call record_stage('pair_generation_and_buffering',nb,stage_start_count, &
                        int(llistt,int64)+sum(int(llistv,int64)),storage='disk')
      stage_start_count = stage_clock()
!
!   Deallocate storage that is no longer required
!
      call dalloc (iqa, 'IQA', 'MCPMPI')
      call dalloc (jqsa, 'JQSA', 'MCPMPI')
      call dalloc (jcupa, 'JCUPA', 'MCPMPI')
      call alcbuf (3)
      call record_stage('csf_and_pair_buffer_release',nb,stage_start_count)
!
!   Write out a report for this run
!
      if (ndef .ne. 0 .and. myid .eq. 0) then
         write (24,*)
         write (24,*) llistt, ' T coefficients generated;'
         do k = 0, kmax
            write (24,*) llistv(k),' V(k=',k,') coefficients generated;'
         enddo
      endif
!
!   Set up sparse structure definition arrays in file 30
!

      stage_start_count = stage_clock()
      rewind (30)
      call setsda (outsdampi, npos, ldbpa(4), nb, myid, nprocs, fhead)
      close (30, status = 'DELETE')
      call record_stage('sparse_structure_and_output',nb,stage_start_count,int(npos,int64),storage='disk')
!
!   Sort MCP coefficients into integral-based lists
!
!   1) the T coefficient
      stage_start_count = stage_clock()
      rewind (31)
      call sort (31, (llistt), ntgi, ldbpa(2), nb, fhead)
      if (ndef .ne. 0 .and. myid .eq. 0)      &
         write (24,*) ntgi, ' I(ab) integrals;'
      close (31, status = 'DELETE')
      call record_stage('sort_and_output',nb,stage_start_count,int(llistt,int64),31,'disk')
!   2) the V coefficient
      do k = 0, kmax
         stage_start_count = stage_clock()
         rewind (k+32)
         call sort (k+32, (llistv(k)), ntgi, ldbpa(3), nb, fhead)
         if (ndef .ne. 0 .and. myid .eq. 0)       &
            write (24,*) ' k = ', k, ': ', ntgi, ' Slater integrals;'
         close (k+32, status = 'DELETE')
         call record_stage('sort_and_output',nb,stage_start_count,int(llistv(k),int64),k+32,'disk')
      enddo
       call dalloc(llistv, 'LLISTV', 'MCPMPI')
  300 format (/'From MCPMPI:')
      return
      end subroutine mcpmpi
