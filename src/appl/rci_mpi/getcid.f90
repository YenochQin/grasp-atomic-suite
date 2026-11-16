!***********************************************************************
!                                                                      *
      subroutine getcid (isofile, rwffile, idblk)
!                                                                      *
!   Interactively determines the data governing the CI problem.        *
!   iccut is replaced by an array iccutblk(1:nblock)
!                                                                      *
!   Call(s) to: [LIB92]: GETYN, NUCPOT, RADGRD, SETQIC.                *
!               [RCI92]: SETISO, SETRWF.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!   Block version by Xinghong He          Last revision: 15 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: nnnp
      use decide_c
      use def_c,           only: cvac, c, z, accy, nelec, emn
      use default_c
      use grid_c
      use hblock_c
      use iccu_c
      use iounit_c
      use npar_c
      use npot_c,          only: zz, nnuc
      use nsmdat_c
      use orb_c
      use wave_c
      use wfac_c
      use blim_c
      use where_c
      use qedcut_c
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setisompi_i
      use setrwfmpi_i
      use getyn_i, only: getyn
      use nucpot_i, only: nucpot
      use radgrd_i, only: radgrd
      use setqic_i, only: setqic
      implicit none
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      character(len=*) :: isofile, rwffile
      character(len=8), dimension(*) ::  idblk
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       logical :: yes
       integer :: i, jblock, ntmp, np10, j
!-----------------------------------------------------------------------
!
! Open, check, load data from, and close the  .iso  file
! Data loaded are: EMN,Z,/NPAR/,/NSMDAT/ where /.../ means whole
! common block.
!
      print *, 'Calling SETISO ...'
      call setisompi (isofile)
!
! The speed of light, if non-default then spread from node-0
! Quantities to be obtained: C
!
      if (ndef .ne. 0) then
!cjb node0 input
        if (myid .eq. 0) then
         write (istde,*) 'Revise the physical speed of light (',cvac,  &
                           ' in a.u.) ?'
         yes = getyn ()
         if (yes) then
            write (istde,*) 'Enter the revised value:'
            read *,c
         else
            c = cvac
         endif
        endif
        call mpi_bcast (c, 1, mpi_double_precision, 0, &
                               mpi_comm_world, ierr)
      else
         c = cvac
      endif
!
! Treat some CSF's as perturbation ? Broadcasting is used only in
! non-default mode, as the above case for speed of light.
! Quantities to be obtained: LFORDR
!
      if (ndef .ne. 0) then
         if (myid .eq. 0) then
            write (istde,*) 'Treat contributions of some CSFs',       &
                            ' as first-order perturbations?'
            lfordr = getyn ()
         endif
         call mpi_bcast (lfordr,1,mpi_logical,0,mpi_comm_world,ierr)
      else
         lfordr = .false.
      endif

! Get iccutblk() from the user-input

      if (.not. lfordr) then
         !...Default first
         do i = 1, nblock
           iccutblk(i) = ncfblk(i)
         enddo
      else

      ! Let master do the i/o, then broadcast
         if (myid .eq. 0) then
            write (istde,*) 'There are ', nblock, 'blocks. They are:'
            write (istde,*) '  block     J Parity     No of CSFs'
            do i = 1, nblock
               write (istde,*) i, idblk(i)(1:5), ncfblk(i)
            enddo

            write (istde,*)
            write (istde,*) 'Enter iccut for each block'
            do jblock = 1, nblock
               write (istde,*) 'Block ', jblock, '   ncf = ',        &
                              ncfblk(jblock)                         &
                      , ' id = ', idblk(jblock)(1:5)
  123          read (istdi,*) ntmp
               if (ntmp .ge. 0 .and. ntmp .le. ncfblk(jblock)) then
                  iccutblk(jblock) = ntmp
               else
                  write (istde,*) 'ICCUT out of range, re-enter:'
                  goto 123
               endif
               write(734,*) ntmp,'! ICCUT for block',jblock
            enddo
         endif

         call mpi_bcast (iccutblk,nblock,mpi_integer,0,            &
                                               mpi_comm_world,ierr)
      endif ! .NOT. LFORDR

!*****************************************************************
!
! Pre-run ?
!
!     IF (IPRERUN .EQ. 0) THEN

!        WRITE (istde,*) ' Prerun with limited interaction?'
!        YES = GETYN ()

!        IF (YES) THEN
!           IPRERUN = 1
!           LTRANS = .FALSE.
!           LVP = .FALSE.
!           LNMS = .FALSE.
!           LSMS = .FALSE.
!           LSE = .FALSE.

!           WRITE (istde,*)  ' Give CSL cut'
!           READ *, NCSFPRE
!           WRITE (istde,*)  ' Give coefficient cut for H_0'
!           READ *, COEFFCUT1
!           WRITE (istde,*) ' Give coefficient cut for the transvers'
!    &,                  ' interaction'
!           GOTO 99
!        ENDIF
!     ENDIF
!*****************************************************************
!
! Include transverse ?
! Quantities to be obtained: LTRANS, WFACT
!
      if (myid .eq. 0) then
         write (istde,*) 'Include contribution of H (Transverse)?'
         ltrans = getyn ()
         write (istde,*) 'Modify all transverse photon frequencies?'
         yes = getyn ()
         if (yes) then
            write (istde,*) 'Enter the scale factor:'
            read *, wfact
         else
            wfact = 1.0d00
         endif
      endif ! myid .EQ. 0
      call mpi_bcast (ltrans, 1, mpi_logical, 0, mpi_comm_world, ierr)
      call mpi_bcast (wfact, 1, mpi_double_precision, 0,            &
                           mpi_comm_world, ierr)
!
! Other interactions ? One logical for each case. Done altogether
!
      if (myid .eq. 0) then
         write (istde,*) 'Include H (Vacuum Polarisation)?'
         lvp = getyn ()

         write (istde,*) 'Include H (Normal Mass Shift)?'
         lnms = getyn ()

         write (istde,*) 'Include H (Specific Mass Shift)?'
         lsms = getyn ()

         write (istde,*) 'Estimate self-energy?'
         lse = getyn ()
         if (lse.eqv..true.) then
             nqedcut = 1
            write (istde,*)                                           &
        'Largest n quantum number for including self-energy for orbital'
            write (istde,*) 'n should be less or equal 8'
            read *, nqedmax
         else
            nqedcut = 0
         end if

         if (ltrans) then
          write(734,'(a)')'y            ! Contribution of H Transverse?'
         else
          write(734,'(a)')'n            ! Contribution of H Transverse?'
         end if
         if (yes) then
           write(734,'(a)') 'y            ! Modify photon frequencies?'
           write(734,*) wfact,'! Scale factor'
         else
           write(734,'(a)') 'n            ! Modify photon frequencies?'
         end if
         if (lvp) then
           write(734,'(a)') 'y            ! Vacuum polarization?'
         else
           write(734,'(a)') 'n            ! Vacuum polarization?'
         end if
         if (lnms) then
           write(734,'(a)') 'y            ! Normal mass shift?'
         else
           write(734,'(a)') 'n            ! Normal mass shift?'
         end if
         if (lsms) then
           write(734,'(a)') 'y            ! Specific mass shift?'
         else
           write(734,'(a)') 'n            ! Specific mass shift?'
         end if
         if (lse) then
           write(734,'(a)') 'y            ! Self energy?'
           write(734,*) nqedmax, '! Max n for including self energy'
         else
           write(734,'(a)') 'n            ! Self energy?'
         end if
      endif ! myid .EQ. 0

      call mpi_bcast (lvp,  1, mpi_logical, 0, mpi_comm_world, ierr)
      call mpi_bcast (lnms, 1, mpi_logical, 0, mpi_comm_world, ierr)
      call mpi_bcast (lsms, 1, mpi_logical, 0, mpi_comm_world, ierr)
      call mpi_bcast (lse,  1, mpi_logical, 0, mpi_comm_world, ierr)
!
! Parameters controlling the radial grid
!
!cff-cjb grid
!cff use default grid
!
! 99  IF (NPARM .EQ. 0) THEN
!        RNT = EXP (-65.D0/16.D0) / Z
!        H = 0.5D0**4
!        N = MIN (220,NNNP)
!     ELSE
!CFF     .. should be Z-dependent
         rnt = 2.0d-06/z
         h = 5.d-2
         n = nnnp
!     ENDIF
      hp = 0.d0
!     IF ( NDEF.NE.0) THEN
!        IF (myid .EQ. 0) THEN
!           WRITE (istde,*) 'The default radial grid parameters',      &
!                          ' for this case are:'
!           WRITE (istde,*) ' RNT = ',RNT,';'
!           WRITE (istde,*) ' H = ',H,';'
!           WRITE (istde,*) ' HP = ',HP,';'
!           WRITE (istde,*) ' N = ',N,';'
!           WRITE (istde,*) ' revise these values?'
!           YES = GETYN ()
!        ENDIF
!
!         ...To prevent subsequent BCAST when YES is false
!        CALL MPI_Bcast (YES, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

!        IF (YES) THEN
!           IF (myid .EQ. 0) THEN
!              WRITE (istde,*) 'Enter RNT:'
!              READ *, RNT
!              WRITE (istde,*) 'Enter H:'
!              READ *, H
!              WRITE (istde,*) 'Enter HP:'
!              READ *, HP
!              WRITE (istde,*) 'Enter N:'
!              READ *, N
!           ENDIF
!           CALL MPI_Bcast (RNT, 1, MPI_DOUBLE_PRECISION, 0,    &
!                               MPI_COMM_WORLD, ierr)
!           CALL MPI_Bcast (H,   1, MPI_DOUBLE_PRECISION, 0,    &
!                               MPI_COMM_WORLD, ierr)
!           CALL MPI_Bcast (HP,  1, MPI_DOUBLE_PRECISION, 0,    &
!                               MPI_COMM_WORLD, ierr)
!           CALL MPI_Bcast (N,   1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!        ENDIF
!     ENDIF ! NDEF.NE.0
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      accy = h**6
!
!   Set up the coefficients for the numerical procedures
!
      call setqic
!
!   Generate the radial grid and all associated arrays
!
      call radgrd
!
!   Generate $- r \times V_ (r)$
!
      call nucpot
!
!   Load the radial wavefunctions
!
      call setrwfmpi (rwffile)
!
!   Write the basic parameters of the model electron cloud to the
!   .res  file; this is the second record on the file --- the
!   first is the header (see SUBROUTINE SETRES)
!
      write (imcdf) nelec, ncf, nw, nblock ! ncf is ncftot, from setcsll
!
!   Write the nuclear parameters and $- r \times V_ (r)$
!   to the  .res  file; these are the third, fourth, and fifth
!   records on the file
!
      write (imcdf) z, emn
      write (imcdf) nparm,(parm(i), i = 1, nparm)
      write (imcdf) n, (zz(i), i = 1, n), nnuc
!
!   Write the physical effects specification to the  .res  file.
!   iccutblk() is now an array of length nblock.
!   This is the sixth record on the file
!
      write (imcdf) c, lfordr, (iccutblk(i), i = 1, nblock),           &
                    ltrans, wfact, lvp, lnms, lsms
!
!   Write the grid data to the  .res  file; this is the seventh
!   record on the file
!
      np10 = n + 10
      write (imcdf) rnt, h, hp, (r(i), i = 1, np10),                   &
                 (rp(i), i = 1, np10), (rpor(i), i = 1, np10) ! (imcdf)
!
!   Write out the interpolated radial wavefunctions; there are
!   2*NW such records; thus the total number of records written
!   at this point is 7+2*NW
!
      do j = 1, nw
         write (imcdf) e(j), gama(j), pz(j), mf(j)
         write (imcdf) (pf(i,j), i = 1, mf(j)), (qf(i,j), i = 1, mf(j))
      enddo
!
      return
      end
