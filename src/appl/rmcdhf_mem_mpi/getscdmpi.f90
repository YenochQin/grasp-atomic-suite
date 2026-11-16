!***********************************************************************
!                                                                      *
      subroutine getscdmpi(eol, idblk, isofile, rwffile, npfile)
!                                                                      *
!   Interactively determines the data governing the SCF problem.       *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT, GETRSL, GETYN, NUCPOT, RADGRD,        *
!                        SETISO, SETQIC, SETRWF.                       *
!               [RSCF92]: GETALD, GETOLD
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 27 Dec 1992   *
!
! Xinghong He 98-08-06
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use parameter_def,    only: nnnw, nnnp
      use coun_c
      use damp_c, only: odamp, cdamp
      use def_c
      use default_c
      use fixd_c, only: lfix
      use iounit_c
      use grid_c, only: h, hp, n, rnt, r
      use invt_c, only: noinvt
      use mcpb_c
      use mpi_c
      use node_c, only: nnodep
      use npar_c
      use npot_c, only: zz
      use npol_c, only: npol_message_len
      use npol_mpi, only: npol_init_mpi, npol_abort_mpi
      use orb_c , only: nak, nh, np, nw, ped
      use orthct_c
      use scf_c, only: scnsty, method
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use setisompi_i
      use setqic_i
      use radgrd_i
      use nucpot_i
      use setrwfmpi_i
      use getaldmpi_i
      use getoldmpi_i
      use convrt_i
      use getrsl_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(out) :: eol
      character  :: isofile*(*)
      character  :: rwffile*(*)
      character(*), intent(in) :: npfile
      character  :: idblk(*)*8
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character*6, parameter :: myname = 'GETSCD'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: np_status
      character(npol_message_len) :: np_message
      integer :: iqadum
      integer , dimension(nnnw) :: indx
      integer :: i, j, iend, ibeg, lenth, nsubs, loc
      real(double) :: odampu, cdampu
      logical :: yes
      character :: record*80, cnum*20
!-----------------------------------------------
!
!   Open, check, load data from, and close the  .iso  file
!
      call setisompi (isofile)  ! NOW ALSO GETS THE RADIAL PARAMETERS, Jon CGrumer (Lund, 2013)
!
!   Set default speed of light and grid parameters
!
      c = cvac
! DEFAULT IS SET ABOVE IN SETISO (from isodata), JG (Lund, 2013)
      if (nparm == 0) then
         rnt = exp((-65.d0/16.d0))/z
         h = 0.5d0**4
         n = min(220,nnnp)
      else
!         default comes here
!CFF     ... RNT should be Z-dependent
         rnt = 2.0d-06/z
         h = 5.d-2
         n = nnnp
      endif
      hp = 0.d0
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      accy = h**6
      if (ndef /= 0) then

                        if (myid .eq. 0) then
                        !--------------------
         write (istde,'(A)')'Change the default speed of',&
            ' light or radial grid parameters?  (y/n) '
         yes = getyn()
         if (yes) then
            write (istde,*) 'Speed of light = ',cvac,';', ' revise ?'
            yes = getyn()
            if (yes) then
               write (istde, *) 'Enter the revised value:'
               read (5, *) c
            endif
!
!   Determine the parameters controlling the radial grid
!
            write (istde, *) 'The default radial grid parameters for ',&
               'this case are:'
            write (istde, *) ' RNT = ', rnt
            write (istde, *) ' H = ', h
            write (istde, *) ' HP = ', hp
            write (istde, *) ' N = ', n
            write (istde, *) ' revise these values?'
            yes = getyn()
            if (yes) then
               write (istde, *) 'Enter RNT:'
               read (5, *) rnt
               write (istde, *) 'Enter H:'
               read (5, *) h
               write (istde, *) 'Enter HP:'
               read (5, *) hp
               write (istde, *) 'Enter N:'
               read (5, *) n
!b Revised grid
               write (istde,*) 'Revised RNT = ', rnt
               write (istde,*) 'Revised H   = ', h
               write (istde,*) 'Revised HP  = ', hp
               write (istde,*) 'Revised N   = ', n
            endif
         endif
!
!b
!b read ACCY on input
!b
         write (istde,*) 'Revise the default ACCY = ', accy
         yes = getyn ()
            if (yes) then
               write (istde,*) 'Enter ACCY:'
               read *, accy
!b Revised ACCY
               write (istde,*) 'Revised ACCY = ', accy
         endif
                        endif ! (myid .EQ. 0)
                        !--------------------

         call mpi_bcast (c,   1, mpi_double_precision, 0,    &
                                 mpi_comm_world, ierr)
         call mpi_bcast (rnt, 1, mpi_double_precision, 0,    &
                                 mpi_comm_world, ierr)
         call mpi_bcast (h,   1, mpi_double_precision, 0,    &
                                 mpi_comm_world, ierr)
         call mpi_bcast (hp,  1, mpi_double_precision, 0,    &
                                 mpi_comm_world, ierr)
         call mpi_bcast (n,   1, mpi_integer, 0, mpi_comm_world, ierr)
!b ACCY
         call mpi_bcast (accy,  1, mpi_double_precision, 0,  &
                                 mpi_comm_world, ierr)
      endif
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
      call npol_init_mpi(npfile, z, n, r, zz, mpi_comm_world, np_status, np_message)
      if (np_status /= 0) call npol_abort_mpi(mpi_comm_world, np_message)
!
!   Load the subshell radial wavefunction estimates
!
      call setrwfmpi (rwffile)
!
!   Set some defaults
!
      thresh = 0.05d0

!      IORDER(I) = I    ! Completely determined in GETOLDMPI
      method(:nw) = 3
      noinvt(:nw) = .true.
!FF      ... Revised Sep. 2014
!CFF      ODAMP(:NW) = 1.D0
      odamp(:nw) = 0.0d0
      scnsty(:nw) = 0.0d0

      where (nak(:nw) < 0)
         nnodep(:nw) = np(:nw) + nak(:nw)
      elsewhere
         nnodep(:nw) = np(:nw) - nak(:nw) - 1
      end where

      if (diag) then
         eol = .false.
         call getaldmpi       ! (E)AL type calculation,
                              ! H(DC) will not be diagonalised
      else if (lfordr) then
         eol = .true.
         call getoldmpi (idblk)   ! EOL type calculation
      else
         if (myid == 0) then
            eol = .true.
         end if
         call mpi_bcast (eol, 1, mpi_logical, 0, mpi_comm_world, ierr)
         if (eol) then
            call getoldmpi (idblk)
         else
            call getaldmpi    ! (E)AL type calculation,
                              ! H(DC) will not be diagonalised
         endif
      endif

      if(myid == 0) then
         write (istde,*) 'Enter the maximum number of SCF cycles:'
         read (*,*) nscf
         if (ndef.eq.0) then
            write(734,*) nscf,'! Number of SCF cycles'
          end if
      end if
      call mpi_bcast (nscf, 1, mpi_integer, 0, mpi_comm_world, ierr)
!
!   Allow the user to modify other defaults
!
      if (ndef /= 0) then
         if(myid == 0) then
            write (istde, '(A)') 'Modify other defaults?  (y/n) '
            yes = getyn()
         end if
         call mpi_bcast (yes, 1, mpi_logical, 0, mpi_comm_world, ierr)
      else
         yes = .false.
      endif

      if (.not.yes) return
!=======================================================================
! From here to end, "other defaults" are handled. For simplicity
! We'll let node-0 do the job and then broadcast results to all
! nodes.
!=======================================================================
                  !-------------------------------------------
      if (myid == 0) then                        ! This is a _big_ IF
                  !-------------------------------------------
!
!   THRESH
!
         write (istde,*)'An oscillation in the large-component of the ',&
            'radial wavefunction is diregarded'
         write (istde, *) 'for the purposes of node counting if its ', &
            'amplitude is less than 1/20 the'
         write (istde, *) 'maximum amplitude.   Revise this?'
         yes = getyn()
         if (yes) then
    3       continue
            write (istde, *) 'Enter the new threshold value:'
            read (*, *) thresh
            if (thresh <= 0.d0) then
               write (istde, *) myname, ': This must exceed 0;'
               go to 3
            endif
         endif
         yes = .false.
!
!   METHOD
!

! Piece only for printing...

         do i = 1, 4

            do j = 1, nw
               if (.not.(method(j)==i .and. .not.lfix(j))) cycle
               write (istde, *) 'Method ', i, ' is used for ', &
                  'integrating the radial differential ', &
                  'equation for subshells'
               go to 9
            end do

            cycle

    9       continue
            iend = 0
            do j = 1, nw
               if (method(j)==i .and. .not.lfix(j)) then
                  ibeg = iend + 1
                  iend = ibeg
                  record(ibeg:iend) = ' '
                  call convrt (np(j), cnum, lenth)
                  ibeg = iend + 1
                  iend = ibeg + lenth - 1
                  record(ibeg:iend) = cnum(1:lenth)
                  ibeg = iend + 1
                  if (nak(j) < 0) then
                     iend = ibeg
                     record(ibeg:iend) = nh(j)(1:1)
                  else
                     iend = ibeg + 1
                     record(ibeg:iend) = nh(j)(1:2)
                  endif
               endif
               if (iend <= 76) cycle
               write (istde, *) record(1:iend)
               iend = 0
            end do
            if (iend<=0 .or. myid/=0) cycle
            write (istde, *) record(1:iend)
         end do

! Reads user inputs and fills array indx(1:nsubs) where nsubs itself
! is an output from getrsl.
! METHOD(1:4) is the only output. indx() and nsubs are discarded

         write (istde, *) 'Select a different integration method for ',&
            'any subshell radial wavefunction?'
         yes = getyn()
         if (yes) then
            do i = 1, 4
               write (istde, *) 'Method ', i, ':'
               call getrsl (indx, nsubs)
               do j = 1, nsubs
                  loc = indx(j)
                  if (lfix(loc)) cycle
                  method(loc) = i
               end do
            end do
         endif
!
!   NOINVT
!
         write (istde, *) 'The first oscillation of the large component'

         do i = 1, nw
            if (.not.(noinvt(i) .and. .not.lfix(i))) cycle
            write (istde, *) 'of the following radial wavefunctions ', &
               'will be required to be positive'
            go to 15
         end do

         write(istde,*)'of all radial wavefunctions will be required ',&
            'to be positive.   Revise this?'
         yes = getyn()
         go to 17
   15    continue
         iend = 0
         do i = 1, nw
            if (noinvt(i) .and. .not.lfix(i)) then
               ibeg = iend + 1
               iend = ibeg
               record(ibeg:iend) = ' '
               call convrt (np(i), cnum, lenth)
               ibeg = iend + 1
               iend = ibeg + lenth - 1
               record(ibeg:iend) = cnum(1:lenth)
               ibeg = iend + 1
               if (nak(i) < 0) then
                  iend = ibeg
                  record(ibeg:iend) = nh(i)(1:1)
               else
                  iend = ibeg + 1
                  record(ibeg:iend) = nh(i)(1:2)
               endif
            endif
            if (iend <= 76) cycle
            write (istde, *) record(1:iend)
            iend = 0
         end do
         if (iend > 0) write (istde, *) record(1:iend)
         write (istde, *) 'Revise this?'
         yes = getyn()
   17    continue
         if (yes) then
            write(istde,*)'Suppressing enforcement of positive first ',&
               'oscillation:'
            call getrsl (indx, nsubs)
            do i = 1, nsubs
               loc = indx(i)
               if (lfix(loc)) cycle
               noinvt(loc) = .true.
            end do
         endif
!
!   ODAMP
!
         do i = 1, nw
            if (.not.(odamp(i)/=0.d0 .and. .not.lfix(i))) cycle
            write (istde, *) 'Subshell accelerating parameters have ', &
               'been set.   Revise these?'
            yes = getyn()
            go to 20
         end do
         write (istde, *) 'Set accelerating parameters for subshell ', &
            'radial wavefunctions?'
         yes = getyn()
   20    continue
         if (yes) then
            write (istde, *) 'Different accelerating parameters for ', &
               'different subshell radial wavefunction?'
            yes = getyn()
            if (yes) then
   21          continue
               write (istde, *) 'Enter an accelerating parameter'
               write (istde, *) ' (0< ODAMP < 1 allows ODAMP to be ', &
                  'reduced as convergence is approached;'
               write (istde, *) ' -1 < ODAMP < 0 implies |ODAMP| is ', &
                  'held constant):'
               read (*, *) odampu
               if (abs(odampu)==0.d0 .or. abs(odampu)>=1.d0) then
                  write (istde, *) myname, ': Value out of range ...'
                  go to 21
               else
                  call getrsl (indx, nsubs)
                  do i = 1, nsubs
                     loc = indx(i)
                     if (lfix(loc)) cycle
                     odamp(loc) = odampu
                  end do
               endif
            else
   23          continue
               write (istde, *) 'Enter the accelerating parameter'
               write (istde, *) ' (0< ODAMP < 1 allows ODAMP to be ', &
                  'reduced as convergence is approached;'
               write (istde, *) ' -1 < ODAMP < 0 implies |ODAMP| is ', &
                  'held constant):'
               read (*, *) odampu
               if (abs(odampu)==0.d0 .or. abs(odampu)>=1.d0) then
                  write (istde, *) myname, ': Value out of range ...'
                  go to 23
               else
                  where (.not.lfix(:nw))
                     odamp(:nw) = odampu
                  end where
               endif
            endif
         endif
!
!   CDAMP
!
         write (istde, *) 'Set accelerating parameters for the ', &
            'eigenvectors?'
         yes = getyn()
         if (yes) then
            write (istde, *) 'Different accelerating parameters ', &
               'for each eigenvector?'
            yes = getyn()
            if (yes) then
               write (istde, *) 'Enter an accelerating parameter for'
               call convrt (ncmin, record, lenth)
               write (istde, *) ' each of the '//record(1:lenth)//' levels :'
               read (*, *) (cdamp(i),i=1,ncmin)
            else
               write (istde, *) 'Enter the accelerating parameter:'
               read (*, *) cdampu
               cdamp(:ncmin) = cdampu
            endif
         endif
!
!   NSIC
!
         write (istde, *) 'Following the improvement of each of the ', &
            'subshell radial wavefunctions in turn, '
         write (istde, *) 'the ', nsic, ' least self-consistent', &
            ' functions will be improved at the'
         write (istde, *) 'end of the first SCF cycle. Revise this ', &
            'setting?'
         yes = getyn()
         if (yes) then
            write (istde, *) 'Enter the number of additional ', 'improvements:'
            read (*, *) nsic
         endif
!
!   NSOLV
!
         write (istde, *) 'The maximum number of cycles in attempting ', &
            'to solve each radial equation is '
         write (istde, *) nsolv, ' times the principal quantum', &
            ' number of the radial'
         write (istde, *) 'wave-function to be estimated.   ', &
            'Revise this setting?'
         yes = getyn()
         if (yes) then
            write (istde, *) 'Enter the factor that multiplies the ', &
               'principal quantum number:'
            read (*, *) nsolv
         endif
!
!   Orthogonalisation
!
         if (orthst) then
            write (istde, *) 'Subshell radial wavefunctions will be ', &
               'Schmidt orthogonalised immediately'
            write (istde, *) 'following their estimation to all ', &
               'functions with poorer self-consistency.'
            write (istde, *) ' Revise this?'
            yes = getyn()
            if (yes) orthst = .false.
         else
            write (istde, *) 'Subshell radial wavefunctions will be ', &
               'Schmidt orthogonalised at the end of'
            write (istde, *) 'each SCF cycle.   Revise this?'
            yes = getyn()
            if (yes) orthst = .true.
         endif
                  !-------------------------------------------
      endif                                      ! end of the _big_ IF
                  !-------------------------------------------

      call mpi_bcast (thresh,  1, mpi_double_precision, 0,            &
                                  mpi_comm_world, ierr)
      call mpi_bcast (method, nw, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast (noinvt, nw, mpi_logical, 0, mpi_comm_world, ierr)
      call mpi_bcast (odamp,  nw, mpi_double_precision, 0,            &
                                  mpi_comm_world, ierr)
      call mpi_bcast (cdamp, ncmin, mpi_double_precision, 0,          &
                                    mpi_comm_world, ierr)
      call mpi_bcast (nsic, 1, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast (nsolv, 1, mpi_integer, 0, mpi_comm_world, ierr)
!cjb  CALL MPI_Bcast (orthst, nw, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast (orthst, 1, mpi_logical, 0, mpi_comm_world, ierr)

      return
      end subroutine getscdmpi
