

!***********************************************************************
!                                                                      *
      subroutine genrwf
!                                                                      *
!   Controls the computation of the subshell radial wavefunctions.     *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, GETRSL, GETYN.                         *
!               [ERWF]: FRMHYD, FRMRWF, FRMTFP, PRTREM, SUMMRY,        *
!                       TFPOT.                                         *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
!   Updated by Xinghong He                Last revision: 11 Nov 1997   *
!    for menu-driven and less questions
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use memory_man
      use coun_c
      use def_c, only: z, c, nnnp
      use default_c
      use grid_c
      use left_c, only: set
      use npar_c
      use orb_c
      use wave_c, only: pf, qf
      use whfrom_c, only: source
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use tfpot_i
      use prtrem_i
      use getrsl_i
      use frmrwf_i
      use frmtfp_i
      use frmhyd_i
      use summry_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(nnnw) :: index
      integer :: j, i, k, nradial, nsubs, loc
      real(double) :: con, fkk
      logical :: all, yes, modify=.false., existed
      character :: infile*128
!-----------------------------------------------

!
!   Set the threshold for node counting
!
      thresh = 0.05
!
!   Set up the Thomas-Fermi potential
!
      call tfpot
!
!   Allocate storage to the arrays that store the subshell
!   radial wavefunction arrays; initialise these and all
!   associated arrays
!
      call alloc (pf,nnnp,nw,'PF', 'GENRWF')
      call alloc (qf,nnnp,nw,'QF', 'GENRWF')

      con = z/c
      con = con*con

      do j = 1, nw
         set(j) = .false.
         source(j) = ' '
         pf(:n,j) = 0.d0
         qf(:n,j) = 0.d0

         k = abs(nak(j))
         if (nparm > 0) then
            gama(j) = dble(k)
         else if (nparm == 0) then
            fkk = dble(k*k)
            if (fkk >= con) then
               gama(j) = sqrt(fkk - con)
            else
               write (istde, *) 'LODRWF: Imaginary gamma parameter ', 'for ', &
                  np(j), nh(j), ' orbital;'
               write (istde, *) 'the point model for the nucleus ', &
                  'is inappropriate for Z > ', c, '.'
               stop
            endif
         endif
      end do
!
!   Write out the complete list of subshell radial wave functions
!
      call prtrem (all)
!
! Direct to read radial functions till finish
!
  123 continue
      if (.not.all) then
  234    continue
         write (istde, *)
         write (istde, *) 'Read subshell radial wavefunctions. ', &
            'Choose one below'
         write (istde, *) '    1 -- GRASP92 File'
         write (istde, *) '    2 -- Thomas-Fermi'
         write (istde, *) '    3 -- Screened Hydrogenic'
         write (istde, *) '    4 -- Screened Hydrogenic [custom Z]'

         read (istdi, *) nradial
         if (nradial<1 .or. nradial>4) then
            write (istde, *) nradial, 'is not a valid choice, redo'
            go to 234
         endif

         if (nradial == 1) then
  345       continue
            write (istde, *) 'Enter the file name (Null then "rwfn.out")'
            read (istdi, '(A)') infile
            if (len_trim(infile) == 0) infile = 'rwfn.out'

            inquire(file=infile, exist=existed)
            if (.not.existed) then
               write (istde, *) ' File "', infile(1:len_trim(infile)), &
                  '" does not exist, redo'
               go to 345
            endif
         endif

         write (istde, *) 'Enter the list of relativistic subshells:'
         open(unit=734, file='tmp_734', status='NEW')
         call getrsl (index, nsubs)
         close(734, status='DELETE')

         if (nradial == 1) then
            call frmrwf (index, nsubs, infile)
         else if (nradial == 2) then
            call frmtfp (index, nsubs)
         else if (nradial == 3) then
            call frmhyd (index, nsubs, modify,nradial)
         else if (nradial == 4) then
            call frmhyd (index, nsubs, modify,nradial)
         endif

         call prtrem (all)
         if (.not.all) go to 234

      endif
!
! All read. Let know, and allow modifying if non default
!
      write (istde, *) 'All required subshell radial wavefunctions ', &
         ' have been estimated:'
      call summry (istde)

      if (ndef == 0) then
         modify = .false.
      else
         write (istde, *) 'Revise any of these estimates?'
         modify = getyn()
      endif

      if (modify) then
  456    continue
         write (istde, *) 'Enter the list of subshells whose radial ', &
            'wavefunctions are to be revised:'
         open(unit=734, file='tmp_734', status='NEW')
         call getrsl (index, nsubs)
         close(734, status='DELETE')
         if (nsubs == 0) go to 456
         do j = 1, nsubs
            loc = index(j)
            set(loc) = .false.
            pf(:n,loc) = 0.d0
            qf(:n,loc) = 0.d0
            source(loc) = ' '
         end do
         all = .false.
         go to 123
      endif

      if (ndef /= 0) call summry (24)

      return
      end subroutine genrwf
