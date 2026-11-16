!***********************************************************************
!                                                                      *
      subroutine matrix (ncore, j2max)
!                                                                      *
!   This SUBROUTINE calls routines to  form the  Hamiltonian  matrix   *
!   and to diagonalise it.   The total angular momenta and parity of   *
!   the ASF are found;  the eigenvectors are chosen so that the sign   *
!   of the largest element is positive.                                *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, HMOUT, IQ              *
!                        WGHTD5.                                       *
!               [RCI92]: MANEIG, QED, SETHAM.                          *
!                                                                      *
!                                         Last revision: 28 Dec 1992   *
!   Modified by Xinghong He               Last revision: 31 Dec 1997   *
!   Block version Xinghong He             Last revision: 12 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double, long
      use parameter_def,   only: nnnw
      use memory_man
      use decide_c
      use debug_c
      use def_c
      use eigv_c
      use hblock_c
      use hmat_c
      use iounit_c
      use jlabl_c, labj=>jlbr, labp=>jlbp
      use orb_c, only: ncf, nw, iqa
      use prnt_c
      use stat_c
      use wave_c
      use where_c
      use blim_c
      use eigvec1_c
      use iccu_c
      use cteilsrk_c
      use coeils_c
      use bilst_c
      use vpilst_c
      use keilst_c
      use vinlst_c
      use fposition_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use auxblk_i
      use lodcslmpi_i
      use qed_slfen_i
      use genmat_i
      use genmat2_i
      use hmout_i
      use iq_i
      use maneig_i
      use engout_i
      use wghtd5_i
      use qed_i
      use mpi_c
      use mpiu_i, only: stopmpi
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ncore
      integer,  intent(in):: j2max
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double), dimension(nnnw) :: slfint
      character(len=8) :: cnum
      real(double) :: atwinv, elsto, eau, ecm, eev, elemnt
      real(double), dimension(:), pointer :: slf_en, ucf, etot
      integer(long) :: nelmnt_a
      integer :: iiatjpo, iiaspar
      integer :: i, j, irestart, ncminpas, jblock, ic, ip, mode
!-----------------------------------------------
!     ...Common to all blocks - place here to save CPU time
      call auxblk (j2max, atwinv)

!***************************************************************
!      Loop over blocks
!***************************************************************
      ncminpas = 0

      do 100 jblock = 1, nblock
         ncf    =   ncfblk(jblock)
         nvec   =   nevblk(jblock)
         nvecmx = ncmaxblk(jblock)
         iccut(1)  = iccutblk(jblock)
         !.. Determine position of the previous block in the .res file

         nposition = 7 + nw + nw  ! File position of the previous block
                                  ! in the .res file
         do i = 1, jblock - 1
            j = (ncfblk(i) - myid - 1 + nprocs) / nprocs
            if (ncfblk(i) .lt. nprocs) j = ncfblk(i) / (myid+1)
            nposition = nposition + j + 1
         enddo

         !.. SETHAM does not handle this extrem case
         if (nprocs .gt. ncf)                                          &
                         call stopmpi ('matrix: too many nodes', myid)

!        ...Obtain ivec() from iccmin()
         if (nvec .gt. 0) then
            call alloc (ivec, nvec, 'IVEC', 'MATRIX')
            do i = 1, nvec
               ivec(i) = iccmin(i+ncminpas)
            enddo
            ncminpas = ncminpas + nvec
         endif

!        ...These 3 were allocated in lodcsh2 and deallocated at the end
!        ... of this routine and in the setham. In this block version,
!        ... both allocation and deallocation are placed here. See the
!        ... following goto 80 for reason.
         call alloc (iqa, nnnw, ncf, 'IQA', 'MATRIX')
         call alloc (jqsa, nnnw, 3, ncf, 'JQSA', 'MATRIX')
         call alloc (jcupa, nnnw, ncf, 'JCUPA', 'MATRIX')
         call alloc (slf_en,ncf,'SLF_EN', 'MATRIX')
         call alloc (ucf, nw,'UCF', 'MATRIX')
          do ic=1,ncf
               slf_en(ic) = 0.0
          enddo

!      ...Load CSF list of the current block
        call lodcslmpi (21, ncore, jblock)
! zou

         if (lse) then
            if (myid .eq. 0) then
              print *, 'Entering QED ...'

              write (24,*)
              write (24,*) ' Self Energy Corrections: '
              write (24,*)
            endif
            call qed_slfen (slfint)
            do ic = 1, ncf
               elemnt = 0.0d00
               do i = 1,nw
                  elemnt = elemnt+iq (i,ic)*slfint(i)
               enddo
               slf_en(ic) = elemnt
            enddo

            if (myid .eq. 0) then
               write (24,*)
               write (24,*) 'Self-energy corrections estimated'       &
                       //' --- these will influence the data'
               write (24,*) ' in the RCI92 MIXing coefficients File.'
               endif
            endif

! zou
         if (nvec <= 0) then
            if (myid .eq. 0) write (25) jblock, ncf, nvec, 999, 999
!           ...Generate H matrix of the current block and write to disk
!           ...eav, nelmnt are also obtained from genmat
            call genmat (atwinv, jblock, myid, nprocs, elsto, irestart,&
           slf_en)
             ! This call is optional, meaning the matrix of this block
               ! does not need to be computed. But don't comment it out
               ! since other programs may assume thet existence of them.
            call genmat2 (irestart, nelmnt_a, elsto)
            goto 80 ! need to clear memory
         endif
!        ------------------------
         call genmat (atwinv, jblock, myid, nprocs, elsto, irestart,   &
           slf_en)
         call genmat2 (irestart, nelmnt_a, elsto)
!
!   Allocate and deallocate memory for the mixing coefficients
!   from the prerun
!
      if (iprerun.eq.1) call alloc (evec1,ncf*nvec, 'EVEC1', 'MATRIX')
      if (iprerun.eq.2) call dalloc (evec1, 'EVEC1', 'MATRIX')

! Since maneig needs both nelmnt and nelmnt_a, # of non-zero
! matrix elements computed by current node and the total #
! from all nodes, a new var nelmnt_a is created rathan than
! using the one in the common block.

!        write (*,*) 'Matrix61,myid=',myid
        call maneig (iiatjpo, iiaspar, nelmnt_a)
!GG        CALL MANEIG (iiatjpo, iiaspar, nelmnt_a, imethod)
!GG         CALL MANEIG (iiatjpo, iiaspar)
!
!   Write out eigenvalues (ENGOUT), dominant components of the
!   eigenvectors (WGHTD5) to stream 24; write out ASF symmetries,
!   eigenvalues  and eigenvectors to RCI92 mixing coefficients file.
!   EAV and ELSTO are added back to energy here
!
        if (myid .eq. 0) then
! ELSTO has never been in Hamiltonian matrix, yet it was
! added to EAV which was later substracted from H. Thus at
! this point, EAV is correct (it has ELSTO added), EVAL()
! need ELSTO and the correct EAV.
           if (ncf > 1) then
              do i = 1, nvec
                 eval(i) = eval(i) + elsto
              enddo
           end if

           call engout (eav,eval,iiatjpo,iiaspar,ivec,nvec,3)
           call wghtd5 (iiatjpo, iiaspar)

!          ...Write ASF symmetries, eigenvalues, and eigenvectors to RCI92
!          ...MIXing coefficients File; close the file; print a report
           write (25) jblock, ncf, nvec, iiatjpo, iiaspar
           write (25) (ivec(i), i = 1, nvec)
           write (25) eav,(eval(i),i = 1,nvec)
           write (25) ((evec(i+(j-1)*ncf),i = 1,ncf),j = 1,nvec)

           print *, 'RCI90 MIXing coefficients File generated.'
        endif
!
!   Save the mixing coefficients from the prerun
!
      if (iprerun .eq. 1) then
         do j = 1, nvec
            do i = 1, ncf
               evec1(i+(j-1)*ncf) = evec(i+(j-1)*ncf)
            enddo
         enddo
      endif
!
!   Estimate diagonal self-energy corrections; the stored
!   eigenvalues and eigenvectors are not modified by these
!   estimates
!
      if (.not.lse) then
         if (myid .eq. 0)            print *, 'Entering QED ...'
            call alloc (etot,nvec,'ETOT', 'MATRIX')
         if (myid .eq. 0) then
            write (24,*)
            write (24,*) ' Self Energy Corrections: '
            write (24,*)
            write (24,301)
            write (24,*)
         end if
  301 format (' Level  J Parity',7x,'Hartrees',14x,'Kaysers',          &
               16x,'eV' )
  302 format (1i3,2x,2a4,1p,3d22.14)
         do j = 1, nvec
            call qed (j,slfint,ucf)
            elemnt = 0.0d00
            ic = ivec(j)
            do i = 1,nw
               elemnt = elemnt+ucf(i)*slfint(i)
!              ELEMNT = ELEMNT+IQ (I,IC)*SLFINT(I)
            enddo
            etot(j) = eval(j)+elemnt
!
            eau = elemnt
            ecm = eau*aucm
            eev = eau*auev
            ip = (iiaspar+3)/2
            if (myid .eq. 0)                                           &
               write (24,302) j,labj(iiatjpo),labp(ip),eau,ecm,eev
!
         enddo
         if (myid .eq. 0) then
            write (24,*)
            write (24,*) 'Self-energy corrections estimated'           &
                    //' --- these do not influence the data'
            write (24,*) ' in the RCI92 MIXing coefficients File.'
            call engout (eav,etot,iiatjpo,iiaspar,ivec,nvec,mode)
! zou       CALL ENGOUT (EAV+elsto,ETOT,IiATJPO,iIASPAR,IVEC,NVEC,MODE)
         endif
         call dalloc (etot, 'ETOT', 'MATRIX')
      endif

!        ...Locals
         call dalloc (ivec, 'IVEC', 'MATRIX')
!        ...Allocated in maneig
         call dalloc (eval, 'EVAL', 'MATRIX')
         call dalloc (evec, 'EVEC', 'MATRIX')

   80    continue

!        ...Locals
         call dalloc (iqa, 'IQA', 'MATRIX')
         call dalloc (jqsa, 'JQSA', 'MATRIX')
         call dalloc (jcupa, 'JCUPA', 'MATRIX')
         call dalloc (slf_en, 'SLF_EN', 'MATRIX')
         call dalloc (ucf, 'UCF', 'MATRIX')

  100 continue
!
!   Close the restart files; nothing will be added to them now
!
      close (imcdf)
!     ...Clean up
      call dalloc (ncfblk, 'NCFBLK', 'MATRIX')
      call dalloc (nevblk, 'NEVBLK', 'MATRIX')
      call dalloc (ncmaxblk, 'NCMAXBLK', 'MATRIX')
      call dalloc (iccutblk, 'ICUTTBLK', 'MATRIX')
      call dalloc (iccmin, 'ICCMIN', 'MATRIX') ! allocated in items as pnccmin

      call dalloc (valteirk, 'VALTEIRK', 'MATRIX') ! allocated in genintrk
      call dalloc (indteirk, 'INDTEIRK', 'MATRIX') ! allocated in genintrk
!
!   Deallocate storage for the integral lists from the
!   Dirac-Coulomb operator; the storage was allocated
!   in IABINT and RKINTC
!
      if (ncoei .gt. 0) then
         call dalloc (indoei, 'INDOEI', 'MATRIX')
         call dalloc (valoei, 'VALOEI', 'MATRIX')
      endif
!
!   Deallocate storage for the integral lists from the
!   transverse photon interaction operator; this storage
!   was allocated in BRINT1, brint2,...brint6
!
      if (ltrans) then
         if (ntpi(1) .gt. 0) then
            call dalloc (indtp1, 'INDTP1', 'MATRIX')
            call dalloc (valtp1, 'VALTP1', 'MATRIX')
         endif
         if (ntpi(2) .gt. 0) then
            call dalloc (indtp2, 'INDTP2', 'MATRIX')
            call dalloc (valtp2, 'VALTP2', 'MATRIX')
         endif
         if (ntpi(3) .gt. 0) then
            call dalloc (indtp3, 'INDTP3', 'MATRIX')
            call dalloc (valtp3, 'VALTP3', 'MATRIX')
         endif
         if (ntpi(4) .gt. 0) then
            call dalloc (indtp4, 'INDTP4', 'MATRIX')
            call dalloc (valtp4, 'VALTP4', 'MATRIX')
         endif
         if (ntpi(5) .gt. 0) then
            call dalloc (indtp5, 'INDTP5', 'MATRIX')
            call dalloc (valtp5, 'VALTP5', 'MATRIX')
         endif
         if (ntpi(6) .gt. 0) then
            call dalloc (indtp6, 'INDTP6', 'MATRIX')
            call dalloc (valtp6, 'VALTP6', 'MATRIX')
         endif
      endif
!
!   Deallocate storage for the nuclear motional energy integral
!   lists; this was allocated in KEINT and VINT
!
      if (lnms) then
         if (nkei .gt. 0) then
            call dalloc (indkei, 'INDKEI', 'MATRIX')
            call dalloc (valkei, 'VALKEI', 'MATRIX')
         endif
      endif
      if (lsms) then
         if (nvinti .gt. 0) then
            call dalloc (indtei, 'INDTEI', 'MATRIX')
            call dalloc (valtei, 'VALTIE', 'MATRIX')
         endif
      endif
!
!   Deallocate storage for the vacuum-polarisation integral list;
!   this was allocated in VPINT
!
      if (lvp) then
         if (nvpi .gt. 0) then
            call dalloc (indvpi, 'INDVPI', 'MATRIX')
            call dalloc (valvpi, 'VALVPI', 'MATRIX')
         endif
      endif

      call dalloc (pf, 'PF', 'MATRIX') ! lodrwf or lodres
      call dalloc (qf, 'QF', 'MATRIX') ! lodrwf or lodres

      return

      end
