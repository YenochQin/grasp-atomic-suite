!***********************************************************************
!                                                                      *
      subroutine mcp (nb, restrt, myid, nprocs, fhead)
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
!   Modified by Gediminas Gaigalas:                         Feb 2017   *
!                               1) for new spin-angular integration,   *
!                               2) for sorting in the memory.          *
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
      use rangular_parallel, only: parallel_pairs_available, begin_pair_block, end_pair_block, load_parallel_pair
      use rangular_reuse, only: reuse_requested, begin_reuse_block, end_reuse_block, load_reused_pair
      use parameter_def,   only:  nnnw, keyorb
      use memory_man
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use outsda_i
      use fndbeg_i
      use alcbuf_i
      use rkco_gg_i
      use setsda_i
      use sort_i
      use sortmem_i
      use alloccheck_i
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
      use sacoef_c
      use onescalar_i, only: onescalar
      use settmpgg_i, only: settmpgg
      implicit none
      external cord
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: nb
      integer  :: myid
      integer  :: nprocs
      logical, intent(in) :: restrt
      character(len=*), intent(in) :: fhead
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = keyorb
      integer, parameter :: keysq = keyorb*keyorb
      real(double), parameter :: cutoff = 1.0d-10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! GG begin
!     2^(23)
      integer, parameter :: nggmax = 8388608
!     2^(24)
!        PARAMETER (NGGMAX = 16777216)
!     2^(27)
!        PARAMETER (NGGMAX   = 134217728)
!        PARAMETER (NGGMAX = 536870912)
! GG end
      integer, dimension(:), pointer :: llistv, mark
      integer :: jastrt,jbstrt,npos,llistt,k,ja,jb,ia,ib,lab,i,ic, &
         id, nswap, iswap, lac, lbd, ntgi
      real(double), dimension(nnnw) :: tshell
      real(double) :: vcoeff
      logical :: f0int, lincr
      integer :: irez, iggmax, iggmax_k, markt, newsiz, jj
      integer(int64) :: stage_start_count
!-----------------------------------------------
!
!   Allocate storage to the array that stores the lengths
!   of the lists of V coefficients of each multipolarity
!
! GG begin
      stage_start_count = stage_clock()
      irez = 1
      iggmax = nggmax
      iggmax_k = iggmax/(2+kmax)
      iggmax = iggmax_k*(2+kmax)
!GG      PRINT*, "The program alloc Memory for sorting:",IGGMAX
      call alloc(llistv, 0, kmax, 'LLISTV', 'MCP')
      call alloc (mark, kmax + 1,'MARK', 'MCP')
      call alloc (tcoeff,iggmax_k,kmax+2,'TCOEFF','MCP')
      call alloc (iiclmn,iggmax_k,kmax+2,'TICLMN','MCP')
      call alloc (iindex,iggmax_k,kmax+2,'TINDEX','MCP')
      call alloc (ilabel,iggmax_k,kmax+2,'TLABEL','MCP')
      markt = 0
      mark = 0
! GG end
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
      if (parallel_pairs_available) call begin_pair_block(nb)
      if (reuse_requested) call begin_reuse_block()
      call record_stage('coefficient_buffer_allocation', nb, stage_start_count)
      stage_start_count = stage_clock()
!
!   Write header to  .dbg  file if appropriate
!
      if (ldbpa(2) .or. ldbpa(3)) write (99,300)
!
!   JA and JB respectively refer to the initial and final states
!   in the list of NCF configurations
!
      do 5 ja = jastrt, ncf, nprocs
         jbstrt = 1
         do 4 jb = jbstrt, ja
            if (diag .or. (lfordr .and. (jb .gt. iccut(nb))) ) then
               if (jb.ne.ja) cycle
            end if
!   LINCR is .TRUE. if npos is to be incremented by 1; there
!   is always a diagonal element in each column
            if (jb .ne. ja) then
               lincr = .true.
            else
               npos = npos + 1
               lincr = .false.
            endif
            if (parallel_pairs_available) call load_parallel_pair(ja,jb,ia,ib,tshell(1))
            if (reuse_requested) call load_reused_pair(ja,jb,ia,ib,tshell(1))
            if (jb .ne. ja) then
               if (.not. parallel_pairs_available .and. .not. reuse_requested) call onescalar(ja,jb,ia,ib,tshell)
               if (ia .ne. 0 .and. ia .ne. ib .and.       &
                      abs (tshell(1)) .gt. cutoff) then
                  npos = npos + 1
                  lincr = .false.
                  llistt = llistt + 1
                  lab = min (ia,ib) * key + max (ia,ib)
! GG begin
                  if(markt .eq. 0) then
                    if(llistt .gt. iggmax_k) then
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
                      newsiz = (iggmax_k+iggmax_k/2)*(2+kmax)
                      call alloccheck (newsiz,irez)
                      if(irez .eq. 1) then
                        iggmax = newsiz
                        iggmax_k = iggmax/(2+kmax)
                        call ralloc                                    &
                             (tcoeff,iggmax_k,kmax+2,'TCOEFF','MCP')
                        call ralloc                                    &
                             (iiclmn,iggmax_k,kmax+2,'TICLMN','MCP')
                        call ralloc                                    &
                             (iindex,iggmax_k,kmax+2,'TINDEX','MCP')
                        call ralloc                                    &
                             (ilabel,iggmax_k,kmax+2,'TLABEL','MCP')
                        iiclmn(llistt,1) = ja
                        iindex(llistt,1) = npos
                        ilabel(llistt,1) = lab
                        tcoeff(llistt,1) = tshell(1)
                      else
!GG                        print*,"LLISTT =",LLISTT
!GG                        print*,
!GG     :             "The program switches to the disk version of
!sorting"
                        call settmpgg (nb, 31, 'tmp')
                        markt = 1
                        do jj = 1, llistt-1
                          write (31) iiclmn(jj,1), iindex(jj,1),       &
                                     ilabel(jj,1), tcoeff(jj,1)
                        end do
                        write (31) ja, npos, lab, tshell(1)
                      end if
                    else
                      iiclmn(llistt,1) = ja
                      iindex(llistt,1) = npos
                      ilabel(llistt,1) = lab
                      tcoeff(llistt,1) = tshell(1)
                    end if
                  else
                    write (31) ja, npos, lab, tshell(1)
                  end if
! GG end
               endif
            endif
!
!   Call the MCP package to generate V coefficients; ac and bd
!   are the density pairs. NVCOEF is initialized here but changed
!   in /rkco/cor[d] via COMMON.
!
            if (.not. parallel_pairs_available .and. .not. reuse_requested) then
               nvcoef = 0
               call rkco_gg (ja, jb, cord, 0, 1)
            end if
            do 2 i = 1, nvcoef
               vcoeff = coeff(i)
               if (abs (vcoeff) .gt. cutoff) then
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
                     if (ia .gt. ic) then
                        iswap = ic
                        ic = ia
                        ia = iswap
                        nswap = nswap + 1
                     endif
                     if (ib .gt. id) then
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
! GG begin
                  if(mark(k+1) .eq. 0) then
                    if(llistv(k) .gt. iggmax_k) then
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
                      newsiz = (iggmax_k+iggmax_k/2)*(2+kmax)
                      call alloccheck (newsiz,irez)
                      if(irez .eq. 1) then
                        iggmax = newsiz
                        iggmax_k = iggmax/(2+kmax)
                        call ralloc                                    &
                             (tcoeff,iggmax_k,kmax+2,'TCOEFF','MCP')
                        call ralloc                                    &
                             (iiclmn,iggmax_k,kmax+2,'TICLMN','MCP')
                        call ralloc                                    &
                             (iindex,iggmax_k,kmax+2,'TINDEX','MCP')
                        call ralloc                                    &
                             (ilabel,iggmax_k,kmax+2,'TLABEL','MCP')
                        iiclmn(llistv(k),k+2) = ja
                        iindex(llistv(k),k+2) = npos
                        ilabel(llistv(k),k+2) = lab
                        tcoeff(llistv(k),k+2) = vcoeff
                      else
!GG                        print*,"K=",K," LLISTV =",LLISTV(K)
!GG                        print*,
!GG     :             "The program switches to the disk version of
!sorting"
                        call settmpgg (nb, 32+k, 'tmp')
                        mark(k+1) = 1
                        do jj = 1,iggmax_k
                          write (32+k) iiclmn(jj,k+2),iindex(jj,k+2),  &
                                       ilabel(jj,k+2),tcoeff(jj,k+2)
                        end do
                        write (32+k) ja, npos, lab, vcoeff
                      end if
                    else
                      iiclmn(llistv(k),k+2) = ja
                      iindex(llistv(k),k+2) = npos
                      ilabel(llistv(k),k+2) = lab
                      tcoeff(llistv(k),k+2) = vcoeff
                    end if
                  else
                    write (32+k) ja, npos, lab, vcoeff
                  end if
! GG end


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
      if (mod(ja,100) == 0) then
         if (ja .eq. 1 .or. ja .eq. ncf .or. mod (ja,100) .eq. 0) then
            print *, 'Row ', ja, '; nnonz = ', npos,';  block = ', nb
         endif
      end if
    5 continue
      call record_stage('pair_generation_and_buffering', nb, stage_start_count, &
                        int(llistt,int64) + sum(int(llistv,int64)), storage='memory_or_spill')
      if (parallel_pairs_available) call end_pair_block(nb)
      if (reuse_requested) call end_reuse_block(nb)
      stage_start_count = stage_clock()
!
!   Deallocate storage that is no longer required
!
      call dalloc (iqa, 'IQA', 'MCP')
      call dalloc (jqsa, 'JQSA', 'MCP')
      call dalloc (jcupa, 'JCUPA', 'MCP')
      call alcbuf (3)
      call record_stage('csf_and_pair_buffer_release', nb, stage_start_count)
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
      call setsda (outsda, npos, ldbpa(4), nb, myid, nprocs, fhead)
      close (30, status = 'DELETE')
      call record_stage('sparse_structure_and_output', nb, stage_start_count, int(npos,int64), storage='disk')
!
!   Sort MCP coefficients into integral-based lists
!
! 1) the T coefficient
      stage_start_count = stage_clock()
      if(markt .eq. 0) then
         call sortmem(31,iggmax_k,(llistt), ntgi, ldbpa(2), nb, fhead)
         call record_stage('sort_and_output', nb, stage_start_count, int(llistt,int64), 31, 'memory')
      else
         rewind (31)
         call sort (31, (llistt), ntgi, ldbpa(2), nb, fhead)
         close (31, status = 'DELETE')
         call record_stage('sort_and_output', nb, stage_start_count, int(llistt,int64), 31, 'disk')
      end if
      if (ndef .ne. 0 .and. myid .eq. 0)                   &
                       write (24,*) ntgi, ' I(ab) integrals;'
! 2) the V coefficient
      do k = 0, kmax
        stage_start_count = stage_clock()
        if(mark(k+1) .eq. 0) then
         call sortmem(k+32,iggmax_k,(llistv(k)),ntgi,ldbpa(3),nb,fhead)
         call record_stage('sort_and_output', nb, stage_start_count, int(llistv(k),int64), k+32, 'memory')
        else
          rewind (k+32)
          call sort (k+32, (llistv(k)), ntgi, ldbpa(3), nb, fhead)
          close (k+32, status = 'DELETE')
          call record_stage('sort_and_output', nb, stage_start_count, int(llistv(k),int64), k+32, 'disk')
        end if
          if (ndef .ne. 0 .and. myid .eq. 0)               &
            write (24,*) ' k = ', k, ': ', ntgi, ' Slater integrals;'
      enddo
      call dalloc (mark, 'MARK','MCP')
      call dalloc (tcoeff,'TCOEFF','MCP')
      call dalloc (iiclmn,'TICLMN','MCP')
      call dalloc (iindex,'TINDEX','MCP')
      call dalloc (ilabel,'TLABLE','MCP')
      call dalloc(llistv, 'LLISTV', 'MCP')
  300 format (/'From MCP:')
      return
      end subroutine mcp
