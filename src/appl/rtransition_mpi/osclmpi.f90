!***********************************************************************
!                                                                      *
      subroutine oscl(name,fullname,tmpdir,startdir,idstring)
      use, intrinsic :: iso_fortran_env, only: real32
!   This routine controls the main sequence of routine calls for the   *
!   calculation  of  data for transitions between multiconfiguration   *
!   Dirac-Fock energy levels.                                          *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use memory_man
      use blk_c
      use biorb_c
      use def_c, ccmps=>speed_of_light_cm_per_s
      use default_c
      use eigv_c
      use orb_c
      use osc_c
      use prnt_c
      use syma_c
      use titl_c
      use wave_c
      use jj2lsjbio_c
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use cpmix_i
      use alcnsa_i
      use alcnta_i
      use connect_i
      use mctout_i
      use readmix_i
      use mctin_i
      use itrig_i
      use bessj_i
      use csfm_i
      use printa_i
      use printals_i

      use suite_system, only: exit => terminate_process
      use ldlbl1_i, only: ldlbl1
      use ldlbl2_i, only: ldlbl2
      use strsum_i, only: strsum
      use sys_chdir_i, only: sys_chdir
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character (len = 24), intent(inout) :: name(2)
      character(len=3)  ::  idstring
      character (len = 128), intent(inout) :: fullname(2)
      character(len=128) :: startdir,tmpdir
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character(len=4), parameter :: iau = 'Hart'
      character(len=4), parameter :: iev = ' eV '
      character(len=4), parameter :: icm = 'Kays'
      character(len=4), parameter :: ihz = ' Hz '
      character(len=4), parameter :: iang = ' A  '
      integer, parameter :: nca = 65536
      integer, parameter :: nfile = 93
      integer, parameter :: nfile1 = 237
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kka, jkp, nfile2, ios, my_nelec, ncftoti, nwi, &
                 nvecsizi, iopar, iblki, ncftotf, nwf,&
                 nvecsizf, iblkf, nvecpr, i, m,       &
                 ielec, levii, levff, itkpo, itest, nlp, lines, &
                 ilbl, inum_ii, inum_ff, icount1, icount2, iii
      real(real64) :: factor, omega, argu, asfa, asfb, asfatmp, asfbtmp
      real(real64), pointer :: et(:), et1(:), ipr(:), ipr1(:), next(:)
      logical :: avail, lsame
      character :: answ*1, iunits*4, g92mix*6
!-----------------------------------------------
!
! NCFI(I): the end position of the Ith block for the initial states in the globle CSF list
! NCFF(I): the end position of the Ith block for the final states in the globle CSF list
!
! for the case that the inital file and final file are same
!
      lsame = trim(fullname(1)) == trim(fullname(2))
      if (lsame.and.mpi_rank.eq.0) call cpmix (fullname, inpci)
      call mpi_barrier(mpi_comm_world,mpi_error_code)
!
! write header for the result file
!
!cjb fort.24
        if (mpi_rank .eq. 0) then
      if (lsame) then
         write(24,*) 'Transition in file:'
         write(24,*) 'f = ',trim(name(1))
      else
         write(24,*) 'Transition between files:'
         write(24,*) 'f1 = ',trim(name(1))
         write(24,*) 'f2 = ',trim(name(2))
      end if
        end if ! myid .EQ. 0
!
      call alcnsa (jja, jjb, hb1, hb2, hc1, hc2, hm1, &
         hm2, lab, nptr, nsdim, 1)
      call alcnta (isldr, isldr1, xsldr, ntdim, 1)
!
!   Make a connection between the orbitals of the
!   merged list and the initial and final state lists
!
      call connect
!
!   Set up units for printing transition energy
!
      if (ltc(1)) then
!
!   Print transition energies in Angstroms
!
         factor = hartree_to_wavenumber_cm
         factor = 1.0e8_real64/factor
         iunits = iang
!
      else if (ltc(2)) then
!
!   Print energies in eV
!
         factor = hartree_to_ev
         iunits = iev
!
      else if (ltc(3)) then
!
!   Print transition energies in Hartree Atomic Units
!
         factor = 1.0_real64
         iunits = iau
!
      else if (ltc(4)) then
!
!   Print transition energies in Hz
!
         factor = hartree_to_wavenumber_cm
         factor = factor*ccmps
         iunits = ihz
!
      else if (ltc(5)) then
!
!   Print transition energies in Kaysers
!
         factor = hartree_to_wavenumber_cm
         iunits = icm
!
      endif
!
!   Select type of transition
!
!   KK  =  0  for electric multipole.
!       =  1  for magnetic multipole.
!
!   CF:  IOPAR  =  (-1)**N      Electric N-pole
!               =  (-1)**(N+1)  Magnetic N-pole.
!                               N > 0
!
      kka = 1
      do jkp = 1, nkp
         nfile2 = nfile1 + jkp
!
! read the head of the file of mixing coef. for initial
!
         if (lsame) fullname(1) = trim(fullname(2))//'_CP'
         if (inpci == 0) then
            open(unit=68,file=trim(fullname(1))//'.cbm',              &
               form='UNFORMATTED',status='OLD')
         else
            open(unit=68, file=trim(fullname(1))//'.bm',              &
               form='UNFORMATTED',status='OLD')
         endif
         if (lsame) fullname(1) = fullname(2)
         read (68, iostat=ios) g92mix
         if (ios/=0 .or. g92mix/='G92MIX') then
            write (*, *) 'Not a GRASP mixing file'
            stop
         endif
         read (68) electron_count, ncftoti, nwi, nvectoti, nvecsizi, nblocki
         write (*, *) '   nelec  = ', my_nelec
         write (*, *) '   ncftoti = ', ncftoti
         write (*, *) '   nwi     = ', nwi
         write (*, *) '   nblocki = ', nblocki
         write (*, *)
!   lbl
         icount1 = 0
         call ldlbl1 (fullname(1))
! If not available generate angular coefficients for all pares of blocks
         call mctout (iopar, jkp, name)
         do iblki = 1, nblocki
            call readmix (name, inpci, 1)
!   lbl
            icount1 = nvecii + icount1
!
! read the head of the file of mixing coef. for final
!
            if (inpci == 0) then
               open(unit=78,file=trim(fullname(2))//'.cbm',          &
                  form='UNFORMATTED',status='OLD')
            else
               open(unit=78, file=trim(fullname(2))//'.bm',          &
                  form='UNFORMATTED',status='OLD')
            endif
            read (78, iostat=ios) g92mix
            if (ios/=0 .or. g92mix/='G92MIX') then
               write (*, *) 'Not a GRASP mixing file'
               stop
            endif
            read (78) electron_count, ncftotf, nwf, nvectotf, nvecsizf, nblockf
            write (*, *) '   nelec  = ', my_nelec
            write (*, *) '   ncftotf = ', ncftotf
            write (*, *) '   nwf     = ', nwf
            write (*, *) '   nblockf = ', nblockf
            write (*, *)
!GG  lbl
            if(iblki .eq. 1) then
               call ldlbl2 (fullname(2))
               if(mpi_rank.eq.0.and.iopen_status1.eq.0.and.iopen_status2.eq.0) then
                  iii = len_trim(startdir)
                  call sys_chdir(trim(startdir),iii,mpi_error_code)
                  if (mpi_error_code.ne.0) then
                     print *, "error changing dir!"
                    call exit(1)
                  endif
                  ilbl = 1
                  call strsum(name,inpci,ilbl)
                  iii = len_trim(tmpdir)
                  call sys_chdir(trim(tmpdir)//'/'//idstring,iii+4,mpi_error_code)
                  if (mpi_error_code.ne.0) then
                     print *, "Error changing dir!"
                     call exit(1)
                  endif
                  write(32,*) 'Transition between files:'
                  write(32,*) trim(name(1))
                  write(32,*) trim(name(2))
               end if
            end if
            icount2 = 0
!GG   lbl end
            do iblkf = 1, nblockf
               call readmix (name, inpci, 2)
!   lbl
               icount2 = nvecff + icount2
!
!   Allocate storage
!
               call alloc (totb, nvecff, 'TOTB', 'OSCL')
               call alloc (totc, nvecff, 'TOTC', 'OSCL')
!
               nvecpr = nvecii*nvecff
               call alloc (et, nvecpr, 'ET', 'OSCL' )
               call alloc (et1, nvecpr, 'ET1','OSCL' )
               call alloc (ipr, nvecpr, 'IPR', 'OSCL' )
               call alloc (ipr1, nvecpr,'IPR1', 'OSCL' )
               call alloc (next, nvecpr, 'NEXT', 'OSCL')
!
!   Initialization for total decay rate
!
               totc(:nvecff) = 0.0_real64
               totb(:nvecff) = 0.0_real64
!
               call mctin (iopar, jkp, name)
!
               if (lk > 0) then
                  ielec = (-1)**lk
                  if (ielec == iopar) then
                     kk = 0
                     kka = 0
                  else
                     kk = 1
                     ielec = -ielec
                  endif
!
!   Set up list of levels for calculation of oscillator strengths
!   sort list into increasing order of energy if option 6 set
!
                  if (mpi_rank.eq.0.and.iblki==1.and.iblkf==1) then
                     if (kk == 0) then
                        write (24, 308) lk
                     else
                        write (24, 309) lk
                     endif
                     write (24, 310)
                     if (ltc(1)) then
                        write (24, 311)
                        if (.not.ltc(7)) then
                           write (24, 312)
                        else
                           write (24, 313)
                        endif
                     else
                        write (24, 314)
                        if (.not.ltc(7)) then
                           write (24, 315) iunits
                        else
                           write (24, 316) iunits
                        endif
                     endif
                  endif !myid=0
!
                  do levii = 1, nvecii
!   lbl
                     inum_ii = icount1 - nvecii + levii
                     do levff = 1, nvecff
!   lbl
                        inum_ff = icount2 - nvecff + levff
!
!   Check for consistent parity and J
!
                        itkpo = lk + lk + 1
                        if (itrig(iatjpoii(levii),iatjpoff(levff),itkpo) == 0) &
                           cycle
                        itest = iasparii(levii)*iasparff(levff)*ielec
                        if (itest < 0) cycle
!
!   Calculate and print transition probability data
!
                        nlp = 70 - 8
                        lines = nlp
!
                        if (lines >= nlp) lines = 0
!
                        m = levff + nvecff*(levii - 1)
!             M = LEVFF+NVECII*(LEVII-1)
                        et(m) = evalff(levff) + eavff - evalii(levii) - eavii
                        if (lsame .and. et(m)<=0.0_real32) cycle
                        omega = -et(m)
                        argu = omega/speed_of_light_au
                        call bessj (argu)
!
!  Calculate oscillator strength between the ASFs
!
                        call csfm (asfatmp,asfbtmp,levii,levff)
                        call mpi_allreduce(asfatmp,asfa,1,mpi_double_precision, &
                             mpi_sum, mpi_comm_world, mpi_error_code)
                        call mpi_allreduce(asfbtmp,asfb,1,mpi_double_precision, &
                             mpi_sum, mpi_comm_world, mpi_error_code)
                        if (mpi_rank .eq. 0) then
                           call printa(asfa,asfb,levii,levff,omega,factor, &
                               lines,lsame)
                           if(iopen_status1.eq.0 .and. iopen_status2 .eq.0) then
                              call printals (inum_ii,inum_ff,asfa,asfb,levii,levff,&
                              omega,factor)
                           end if
                        end if
!              WRITE (24,317)
                     end do
                  end do
               endif
!
!   Deallocate storage; this is local to OSCL
!
               call dalloc (totb, 'TOTB', 'OSCL')
               call dalloc (totc, 'TOTC', 'OSCL')
!
               call dalloc (et, 'ET', 'OSCL')
               call dalloc (et1, 'ET!', 'OSCL')
               call dalloc (ipr, 'IPR', 'OSCL')
               call dalloc (ipr1, 'IPR1', 'OSCL')
               call dalloc (next, 'NEXT', 'OSCL')
!
!   Deallocate storage; this is allocated in READMIX
!
               call dalloc (evalff, 'EVALFF', 'OSCL')
               call dalloc (evecff, 'EVECFF', 'OSCL')
               call dalloc (ivecff, 'IVECFF', 'OSCL')
               call dalloc (iatjpoff, 'IATJPOFF', 'OSCL')
               call dalloc (iasparff, 'IASPARFF', 'OSCL')

            end do
            close(78)
!
!   Deallocate storage; this is allocated in READMIX
!
            call dalloc (evalii, 'EVALII', 'OSCL')
            call dalloc (evecii, 'EVECII', 'OSCL')
            call dalloc (ivecii, 'IVECII', 'OSCL')
            call dalloc (iatjpoii, 'IATJPOII', 'OSCL')
            call dalloc (iasparii, 'IASPARII', 'OSCL')

         end do
         close(68)
         close(nfile2)
      end do
      call dalloc (kp, 'KP', 'OSCL')
!
! close and delete duplicated mixing file
      if (lsame) then
         if (inpci == 0) then
            open(68, file=trim(fullname(2))//'_CP.cbm')
         else
            open(68, file=trim(fullname(2))//'_CP.bm')
         endif
!GG         CLOSE(68, STATUS='delete')
      endif
!
      call alcnsa (jja, jjb, hb1, hb2, hc1, hc2, hm1, hm2, lab, nptr, nsdim, 3)
      call alcnta (isldr, isldr1, xsldr, ntdim, 3)
!
!   This was allocated in LOAD
!
!      CALL DALLOC (PF, 'QF', 'OSLC')
!      CALL DALLOC (QF, 'QF', 'OSCL')
!      CALL DALLOC (IVEC, 'IVEC', 'OSCL')
!      CALL DALLOC (IATJPO, 'IATJPO', 'OSCL')
!      CALL DALLOC (IASPAR, 'IASPAR', 'OSCL')
!      CALL DALLOC (EVAL, 'EVAL', 'OSCL')
!      CALL DALLOC (EVEC, 'EVEC', 'OSCL')
!
!   Close all files
!
      if (mpi_rank .eq. 0) close(24)
!
      return
!
  302 format(/,' ***** Warning *****')
  303 format(/,/,' ***** Error in OSCL *****')
  307 format(/,' Dynamic allocation computed incorrectly: Bug.')
  308 format(/,/,' Electric 2**(',i2,')-pole transitions')
  309 format(/,/,' Magnetic 2**(',i2,')-pole transitions')
  310 format(1x,33('='))
  311 format(/,'   Upper state        Lower state  ',8x,'Gauge',8x,'Wavelength'&
         ,13x,'Einstein coefficients',13x,'Oscillator')
  312 format(81x,'-1',15x,'3 -2 -1',/,' Level  J Parity',4x,'Level  J ',&
         'Parity',21x,'(Angstroms)',10x,'A (s  )',9x,'gB (m s  J  )',7x,&
         'strength gf'/)
  313 format(' Level  J Parity',4x,'Level  J Parity',21x,'(Angstroms)',10x,&
         'A (au)',13x,'gB (au)',10x,'strength gf'/)
  314 format(/,' Upper       Lower ')
  315 format(' Lev  J P',3x,'Lev  J P',7x,'E (',a4,')',9x,'A (s-1)',10x,'gf',&
         12x,'S')
  316 format(' Level  J Parity',4x,'Level  J Parity',23x,'(',a4,')',13x,&
         'A (au)',13x,'gB (au)',10x,'strength gf'/)
  317 format(/,1x,124('+'))
  318 format(/,/,' Radiative lifetimes '/,' ======================='/,/,&
         ' Level      Lifetime  s (-1)')
  319 format(1x,i4,6x,'Coulomb: ',1p,1d20.7)
  320 format(10x,'Babushkin:',1p,1d20.7,/)
  321 format(1x,i4,5x,'Magnetic: ',1p,1d20.7,/)
      return
!
      end subroutine oscl
