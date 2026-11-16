!***********************************************************************
!                                                                      *
      subroutine ris_cal_mpi (name)
!                                                                      *
!   MPI version of RIS_CAL. Each rank accumulates a strided subset of  *
!   configuration-pair contributions and the partial sums are reduced   *
!   before the root rank writes the final output.                       *
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: max_orbitals
      use memory_man
      use debug_c
      use decide_c
      use def_c
      use grid_c
      use npar_c
      use prnt_c
      use syma_c
      use orb_c
      use teilst_c
      use buffer_c
      use ris_c
      use jlabl_c, labj=> jlbr, labp=>jlbp
      use eigv_c
      use mpi_c
      use rintdens_i
      use rinti_i
      use rint_i
      use rintdensvec_i
      use rinti_nms_i
      use cre_i
      use vinti_i
      use rint_sms2_i
      use rint_sms3_i
      use angdata_i
      use getyn_i
      use densread_i
      use densread_seltz_i
      use densnew_seltz_i
      use densread_mpi_i
      use densread_seltz_mpi_i
      use densnew_i
      use densnew_mpi_i
      use densnew_seltz_mpi_i
      use smsread_i
      use smsread_mpi_i
      use smsnew_i
      use smsnew_mpi_i
      use edensityfit_i
      use mpiu_i, only: gdrsummpi_root
      use starttime_i, only: starttime
      use stoptime_i, only: stoptime
      implicit none (type, external)
      character(len=24), intent(in) :: name
      real(real64), pointer :: dint1vec(:, :, :)
      real(real64), pointer   :: dens1vec(:, :)
      real(real64), pointer     :: densfit(:)
      real(real64), pointer   :: fmat(:, :)
      real(real64), pointer     :: rho(:)
      real(real64), pointer     :: res(:)
      real(real64) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
      real(real64) :: dint1(max_orbitals, max_orbitals), dint2(max_orbitals, max_orbitals), dint3(max_orbitals, max_orbitals), &
                                            dint4(max_orbitals, max_orbitals), dint5(max_orbitals, max_orbitals), &
                                               dint6(max_orbitals, max_orbitals), &
                                            dint7(max_orbitals, max_orbitals)
      real(real64) :: au2fm, rcre, hzsmsu, hznmsu, fo90
      logical :: avail_tb, avail_ob, yes2, saveyes, serial_save
      integer :: i, j, ncount1, doit_ob, doit_tb, nrnuc, yes2_mpi
      integer :: avail_ob_mpi, avail_tb_mpi, serial_save_mpi
      if (mpi_rank .eq. 0) then
         write(*,*) '-------------------------------'
         write(*,*) 'RIS_CAL_MPI: Execution Begins ...'
         write(*,*) '-------------------------------'
      endif

      doit_ob = 0
      doit_tb = 0
      yes2_mpi = 0
      avail_ob_mpi = 0
      avail_tb_mpi = 0
      serial_save_mpi = 0

      fo90 = nuclear_model_parameters(1)+2._real64*log(3._real64)*nuclear_model_parameters(2)
      do i=1,max_radial_points
         if(radial_grid(i+1).gt.fo90.and.radial_grid(i).le.fo90) then
            nrnuc = i+1
            exit
         end if
      end do

      allocate( dens1vec(nvec,nrnuc),dint1vec(max_orbitals,max_orbitals,nrnuc) )
      allocate( densfit(nrnuc) )
      dint1vec(:,:,:) = 0._real64
      dens1vec(:,:) = 0._real64
      allocate( fmat(6,nvec),rho(nvec),res(nvec) )

      call alloc (smsc1, nvec,'SMSC1', 'RIS_CAL_MPI')
      call alloc (smsc2, nvec,'SMSC2', 'RIS_CAL_MPI')
      call alloc (dens1, nvec,'DENS1','RIS_CAL_MPI')
      call alloc (dens2, nvec,'DENS2','RIS_CAL_MPI')
      call alloc (dens3, nvec,'DENS3','RIS_CAL_MPI')
      call alloc (dens4, nvec,'DENS4','RIS_CAL_MPI')
      call alloc (dens5, nvec,'DENS5','RIS_CAL_MPI')
      call alloc (dens6, nvec,'DENS6','RIS_CAL_MPI')
      call alloc (dens7, nvec,'DENS7','RIS_CAL_MPI')

      call starttime (ncount1, 'RIS_CAL_MPI')

      do i = 1,nvec
         smsc1(i)  = 0.0_real64
         smsc2(i)  = 0.0_real64
         dens1(i) = 0.0_real64
         dens2(i) = 0.0_real64
         dens3(i) = 0.0_real64
         dens4(i) = 0.0_real64
         dens5(i) = 0.0_real64
         dens6(i) = 0.0_real64
         dens7(i) = 0.0_real64
      enddo

      if (mpi_rank .eq. 0) then
         write(*,*) ' Compute higher order field shift electronic factors?'
         yes2 = getyn ()
         if (yes2) then
            yes2_mpi = 1
         else
            yes2_mpi = 0
         endif
      endif
      call mpi_bcast(yes2_mpi, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
      yes2 = yes2_mpi .eq. 1

      do i = 1,orbital_count
         do j = 1,orbital_count
            if (orbital_kappa(i).eq.orbital_kappa(j)) then
               dint1(i,j) = rintdens(i,j)
               if(yes2) then
                  call rintdensvec(i,j,dint1vec,nrnuc)
               end if
               call rinti_nms(i,j,dint2(i,j),dint7(i,j))
               dint3(i,j) = rint(i,j,1)
               dint4(i,j) = rint(i,j,2)
               dint5(i,j) = rint(i,j,-1)
               dint6(i,j) = rint(i,j,-2)
            else
               dint1(i,j) = 0.0_real64
               dint1vec(i,j,:) = 0.0_real64
               dint2(i,j) = 0.0_real64
               dint3(i,j) = 0.0_real64
               dint4(i,j) = 0.0_real64
               dint5(i,j) = 0.0_real64
               dint6(i,j) = 0.0_real64
               dint7(i,j) = 0.0_real64
            endif
         end do
      end do

      do i = 1,orbital_count
         do j = 1,orbital_count
            if (i.ne.j) then
               rcre = cre(orbital_kappa(i),1,orbital_kappa(j))
               if (dabs(rcre) .gt. cutoff) then
                  vint  (i,j) = vinti(i,j)
                  vint2(i,j) = vint(i,j)                                  &
                     + rint_sms2(i,j)/rcre                                &
                     + rint_sms3(i,j)
               else
                  vint (i,j) = 0.0_real64
                  vint2(i,j) = 0.0_real64
               endif
            else
               vint (i,j) = 0.0_real64
               vint2(i,j) = 0.0_real64
            endif
         end do
      end do

      if (mpi_rank .eq. 0) then
         call angdata(name,avail_ob,1)
         if (avail_ob) close(50)
         call angdata(name,avail_tb,2)
         if (avail_tb) close(51)
         avail_ob_mpi = merge(1, 0, avail_ob)
         avail_tb_mpi = merge(1, 0, avail_tb)
         if ((.not. avail_ob) .and. (.not. avail_tb)) then
            print *,' Save ang. coefficients of one- and two-body op.?'
            saveyes = getyn ()
            print *
            if(saveyes) then
               doit_ob = 1
               doit_tb = 1
            endif
         elseif (.not. avail_ob) then
            print *,' Save ang. coefficients of one-body op. ?'
            saveyes = getyn ()
            print *
            if(saveyes) doit_ob = 1
         elseif (.not. avail_tb) then
            print *,' Save ang. coefficients of two-body op. ?'
            saveyes = getyn ()
            print *
            if(saveyes) doit_tb = 1
         endif
      endif
      call mpi_bcast(avail_ob_mpi, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
      call mpi_bcast(avail_tb_mpi, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
      call mpi_bcast(doit_ob, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
      call mpi_bcast(doit_tb, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
      avail_ob = avail_ob_mpi .eq. 1
      avail_tb = avail_tb_mpi .eq. 1
      if (mpi_rank .eq. 0) then
         ! Saving angular coefficients still has to stay on the root
         ! process because the serial file format is shared.
         serial_save = ((.not. avail_ob .and. doit_ob .eq. 1) .or.      &
                        (.not. avail_tb .and. doit_tb .eq. 1))
         serial_save_mpi = merge(1, 0, serial_save)
      endif
      call mpi_bcast(serial_save_mpi, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
      serial_save = serial_save_mpi .eq. 1

      if (serial_save) then
         call mpi_barrier(mpi_comm_world, mpi_error_code)
         if (mpi_rank .eq. 0) then
            j = index(name,' ')
            if (avail_ob) then
               open(unit=50,file = name(1:j-1)//'.IOB',status='OLD',form='UNFORMATTED')
               if (yes2) then
                  call densread_seltz(dint1,dint2,dint3,dint4,dint5,dint6, &
                       dint7,dint1vec,dens1vec,nrnuc)
               else
                  call densread(dint1,dint2,dint3,dint4,dint5,dint6,dint7)
               endif
               close(50)
            else
               if (doit_ob .eq. 1) then
                  open(unit=50,file = name(1:j-1)//'.IOB',status='UNKNOWN',form='UNFORMATTED')
               endif
               if (yes2) then
                  call densnew_seltz(doit_ob,dint1,dint2,dint3,dint4,dint5,dint6, &
                       dint7,dint1vec,dens1vec,nrnuc)
               else
                  call densnew(doit_ob,dint1,dint2,dint3,dint4,dint5,dint6,dint7)
               endif
            endif

            if (avail_tb) then
               open(unit=51,file = name(1:j-1)//'.ITB',status='OLD',form='UNFORMATTED')
               call smsread(vint,vint2)
               close(51)
            else
               if (doit_tb .eq. 1) then
                  open(unit=51,file = name(1:j-1)//'.ITB',status='UNKNOWN',form='UNFORMATTED')
               endif
               call smsnew(doit_tb,vint,vint2)
            endif
         endif
         call mpi_barrier(mpi_comm_world, mpi_error_code)
         go to 800
      endif

      if (avail_ob) then
         j = index(name,' ')
         open(unit=50,file = name(1:j-1)//'.IOB',status='OLD',form='UNFORMATTED')
         if (yes2) then
            call densread_seltz_mpi(dint1,dint2,dint3,dint4,dint5,dint6, &
                 dint7,dint1vec,dens1vec,nrnuc)
         else
            call densread_mpi(dint1,dint2,dint3,dint4,dint5,dint6,dint7)
         end if
         close(50)
      else
         if(yes2) then
            call densnew_seltz_mpi(doit_ob,dint1,dint2,dint3,dint4,dint5,dint6, &
                 dint7,dint1vec,dens1vec,nrnuc)
         else
            call densnew_mpi(doit_ob,dint1,dint2,dint3,dint4,dint5,dint6,dint7)
         end if
      endif

      if (avail_tb) then
         j = index(name,' ')
         open(unit=51,file = name(1:j-1)//'.ITB',status='OLD',form='UNFORMATTED')
         call smsread_mpi(vint,vint2)
         close(51)
      else
         call smsnew_mpi(doit_tb,vint,vint2)
      endif

      call mpi_barrier(mpi_comm_world, mpi_error_code)

      call gdrsummpi_root(smsc1(:), nvec)
      call gdrsummpi_root(smsc2(:), nvec)
      call gdrsummpi_root(dens1(:), nvec)
      call gdrsummpi_root(dens2(:), nvec)
      call gdrsummpi_root(dens3(:), nvec)
      call gdrsummpi_root(dens4(:), nvec)
      call gdrsummpi_root(dens5(:), nvec)
      call gdrsummpi_root(dens6(:), nvec)
      call gdrsummpi_root(dens7(:), nvec)
      if (yes2) call gdrsummpi_root(dens1vec(:,:), nvec*nrnuc)

  800 continue
      if (mpi_rank .ne. 0) go to 900

      write (24,'(a24,i3)') 'Number of eigenvalues: ',nvec
      write(24,*)
      write(24,*)
      write (24,309)
      do i = 1, nvec
         write (24,323) ivec(i),labj(iatjpo(i)),labp((iaspar(i)+3)/2), &
              hamiltonian_energy_shift+shifted_eigenvalues(i)
      end do
      write (24,308)
      do i = 1,nvec
         write(24,314)
         write (24,313) ivec(i),labj(iatjpo(i)),labp((iaspar(i)+3)/2), &
              dens7(i),(dens2(i)-dens7(i)),dens2(i)
         hznmsu=dens7(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s
         write (24,316) hznmsu*(1.0e-9_real64), &
              (dens2(i)-dens7(i))*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64), &
              dens2(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64)
      end do

      write (24,302)
      do i = 1,nvec
         write(24,314)
         write (24,313) ivec(i),labj(iatjpo(i)),labp((iaspar(i)+3)/2), &
              smsc1(i),smsc2(i)-smsc1(i),smsc2(i)
         hzsmsu=smsc1(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s
         if (dabs(hzsmsu) .le. (10**8)) then
            write (24,315) hzsmsu*(1.0e-6_real64), &
                 (smsc2(i)-smsc1(i))*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-6_real64), &
                 smsc2(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-6_real64)
         else
            write (24,316) hzsmsu*(1.0e-9_real64), &
                 (smsc2(i)-smsc1(i))*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64), &
                 smsc2(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64)
         endif
      end do
      if(yes2) then
         write (24,307)
         do i = 1,nvec
            call edensityfit(radial_grid,dens1vec(i,:),nuclear_charge,nuclear_model_parameters,nrnuc,fmat(:,i),rho(i),res(i))
            write (24,303) ivec(i),labj(iatjpo(i)), &
                 labp((iaspar(i)+3)/2),rho(i)
         end do
         write (24,317)
         do i = 1,nvec
            write (24,324) ivec(i),labj(iatjpo(i)), &
                 labp((iaspar(i)+3)/2),fmat(1,i),fmat(2,i), &
                 fmat(3,i),fmat(4,i),res(i)
         end do
         write (24,327)
         do i = 1,nvec
            write (24,334) ivec(i),labj(iatjpo(i)), &
                 labp((iaspar(i)+3)/2),fmat(5,i),fmat(6,i)
         end do
         write(24,*)
      else
         write (24,306)
         do i = 1,nvec
            write (24,303) ivec(i),labj(iatjpo(i)), &
                 labp((iaspar(i)+3)/2),dens1(i)
         end do
      end if

  900 continue
      call dalloc (smsc1,'SMSC1','RIS_CAL_MPI')
      call dalloc (smsc2,'SMSC2','RIS_CAL_MPI')
      call dalloc (dens1,'DENS1','RIS_CAL_MPI')
      call dalloc (dens2,'DENS2','RIS_CAL_MPI')
      call dalloc (dens3,'DENS3','RIS_CAL_MPI')
      call dalloc (dens4,'DENS4','RIS_CAL_MPI')
      call dalloc (dens5,'DENS5','RIS_CAL_MPI')
      call dalloc (dens6,'DENS6','RIS_CAL_MPI')
      call dalloc (dens7,'DENS7','RIS_CAL_MPI')

      if (mpi_rank .eq. 0) then
         write(*,*) '-------------------------------'
         write(*,*) 'RIS_CAL_MPI: Execution Finished ...'
         write(*,*) '-------------------------------'
      endif
      call stoptime (ncount1, 'RIS_CAL_MPI')
      return
301   format (//' CUTOFF set to ',1pd22.15)
302   format (//' Level  J Parity  Specific mass shift parameter')
303   format (1x,i3,5x,2a4,3x,d20.10)
304   format (1x,i3,5x,2a4,3x,2d20.10)
324   format (1x,i3,5x,2a4,3x,4d20.10,f10.4)
334   format (1x,i3,5x,2a4,3x,2d20.10)
306   format (//' Electron density in atomic units' &
           //' Level  J Parity',8x,'DENS (a.u.)'/)
307   format (//' Level  J Parity  Electron density in atomic units' &
           //24x,'Dens. (a.u.)')
317   format (//' Level  J Parity  Field shift electronic factors and av&
           &erage point discrepancy in fit' &
           //24x, &
           'F0 (GHz/fm^2)',6x,' F2 (GHz/fm^4)', &
           6x,' F4 (GHz/fm^6)',6x,' F6 (GHz/fm^8)', &
           6x,' Disc. (per mille)')
327   format (//' Level  J Parity  Field shift electronic factors (corre&
           &cted for varying density inside nucleus)' &
           //24x, &
           'F0VED0 (GHz/fm^2)   F0VED1 (GHz/fm^4)')
308   format (//' Level  J Parity  Normal mass shift parameter')
309   format (' Level  J Parity  Energy')
313   format (1x,i3,5x,2a4,3x,3d20.10,'  (a.u.)')
333   format (1x,i3,5x,2a4,3x,3d20.10,'  (GHz u)')
323   format (1x,i3,5x,2a4,3x,d20.10,'  (a.u.)')
314   format (/,29x,'<K^1>',13x,'<K^2+K^3>',9x,'<K^1+K^2+K^3>')
315   format (20x,3d20.10 ,2x,'(MHz u)')
316   format (20x,3d20.10,2x,'(GHz u)')
      return
      end subroutine ris_cal_mpi
