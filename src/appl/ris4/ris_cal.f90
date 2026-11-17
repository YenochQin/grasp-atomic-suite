!***********************************************************************
!                                                                      *
      subroutine ris_cal (name)
!                                                                      *
!   This routine controls the main sequence of routine calls for the   *
!   calculation  of the MS parameters, the electron density at the     *
!   origin and radial expectation values.                              *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, GETYN                          *
!               [SMS92]: RINTDENS, VINTI                               *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!                                         Last revision: 10 Nov 1995   *
!                                                                      *
!   Modified by C. Naz\'e  Feb. 2011                                   *
!   Modified by J. Ekman   Jan. 2014                                   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: max_orbitals
      use memory_man
      use debug_c, only: cutoff, double, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use decide_c, only: double, lfordr, lnms, lse, lsms, ltrans, lvp
      use def_c
      use grid_c
      use npar_c, only: double, nuclear_model_parameter_count, nuclear_model_parameters
      use prnt_c, only: byte, double, ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use syma_c, only: double, iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
      use orb_c
      use teilst_c, only: ndtea, ntei, pindte, pvalte, firstt
      use buffer_c
      use ris_c
      use jlabl_c, labj=> jlbr, labp=>jlbp
      use ris_c
      use eigv_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rintdens_i, only: rintdens
      use rinti_i, only: rinti
      use rint_i, only: rint
      use rintdensvec_i, only: rintdensvec
      use rinti_nms_i, only: rinti_nms
      use cre_i, only: cre
      use vinti_i, only: vinti
      use rint_sms2_i, only: rint_sms2
      use rint_sms3_i, only: rint_sms3
      use angdata_i, only: angdata
      use input_io, only: fail_input_io
      use getyn_i, only: getyn
      use densread_i, only: densread
      use densnew_i, only: densnew
      use smsread_i, only: smsread
      use smsnew_i, only: smsnew
      use edensityfit_i, only: edensityfit
      use densnew_seltz_i, only: densnew_seltz
      use densread_seltz_i, only: densread_seltz
      use starttime_i, only: starttime
      use stoptime_i, only: stoptime
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=24), intent(in) :: name
!      REAL(DOUBLE), INTENT(IN) :: DR2
!      INTEGER, INTENT(IN)      :: NOPAR
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64), pointer :: dint1vec(:, :, :)
      real(real64), pointer   :: dens1vec(:, :)
      real(real64), pointer     :: densfit(:)
      real(real64), pointer     :: fmat(:, :)
      real(real64), pointer     :: rho(:)
      real(real64), pointer     :: res(:)


      real(real64)      :: tshell(max_orbitals)
      real(real64) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
      real(real64) :: dint1(max_orbitals, max_orbitals), dint2(max_orbitals, max_orbitals), dint3(max_orbitals, max_orbitals), &
                                            dint4(max_orbitals, max_orbitals), dint5(max_orbitals, max_orbitals), &
                                               dint6(max_orbitals, max_orbitals), &
                                            dint7(max_orbitals, max_orbitals)
      real(real64) :: au2fm, rcre, hzsmsu, hznmsu, fo90
      character(len=11) :: cnum
      character(len=2)  :: ck
      logical :: vsh, nucde, smssh, yes, yes2, avail_tb, avail_ob
      integer :: i, j, l, ncount1, doit_ob, doit_tb, nrnuc, cache_status
!-----------------------------------------------
!

!
!   Matrix elements smaller than CUTOFF are not accumulated
!
!      PARAMETER (CUTOFF = 1.0D-10)

      write(*,*) "-------------------------------"
      write(*,*) "RIS_CAL: Execution Begins ..."
      write(*,*) "-------------------------------"

!     Determine point end point NRNUC for density fit
      fo90 = nuclear_model_parameters(1)+2.0_real64*log(3.0_real64)*nuclear_model_parameters(2)
      !PARM(1) - Fermi parameter c
      !PARM(2) - Fermi parameter a
      do i=1,max_radial_points
         if(radial_grid(i+1)>fo90.and.radial_grid(i)<=fo90) then
            nrnuc = i+1
            exit
         end if
      end do


!   Allocate memory

      allocate( dens1vec(nvec,nrnuc),dint1vec(max_orbitals,max_orbitals,nrnuc) )    ! JE ADD
      allocate( densfit(nrnuc) )                                   ! JE ADD

      dint1vec(:,:,:) = 0.0_real64                                      ! JE ADD
      dens1vec(:,:) = 0.0_real64                                        ! JE ADD

      allocate( fmat(6,nvec),rho(nvec),res(nvec) )                ! JE ADD

!
!   Allocate storage for local arrays
!
      call alloc (smsc1, nvec,"SMSC1", "RIS_CAL")
      call alloc (smsc2, nvec,"SMSC2", "RIS_CAL")
      call alloc (dens1,nvec,"DENS1","RIS_CAL")
      call alloc (dens2,nvec,"DENS2","RIS_CAL")
      call alloc (dens3,nvec,"DENS3","RIS_CAL")
      call alloc (dens4,nvec,"DENS4","RIS_CAL")
      call alloc (dens5,nvec,"DENS5","RIS_CAL")
      call alloc (dens6,nvec,"DENS6","RIS_CAL")
      call alloc (dens7,nvec,"DENS7","RIS_CAL")

      call starttime (ncount1, "RIS_CAL")

!
!   Initialise
!
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
      end do

!      VSH   = .TRUE.
!      SMSSH = .TRUE.
!
!   Calculate all integrals needed for the volume shift calc.
!
      !      IF (VSH) THEN

      write(*,*) " Compute higher order field shift electronic factors?"
      yes2 = getyn ()

      do 5 i = 1,orbital_count
        do 4 j = 1,orbital_count
          if (orbital_kappa(i)==orbital_kappa(j)) then
            dint1(i,j) = rintdens(i,j)
            !GG                 DINT2(I,J) = RINTI(I,J,1)
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
          end if
    4   continue
    5 continue
!      ENDIF
!
!   Calculate and save the Vinti integrals
!
!      IF (SMSSH) THEN
      do 7 i = 1,orbital_count
        do 6 j = 1,orbital_count
          if (i/=j) then
            rcre = cre(orbital_kappa(i),1,orbital_kappa(j))
            if (dabs(rcre) > cutoff) then
              vint  (i,j) = vinti(i,j)
              vint2(i,j) = vint(i,j)                                   &
                         + rint_sms2(i,j)/rcre                         &
                         + rint_sms3(i,j)
            else
             vint (i,j) = 0.0_real64
             vint2(i,j) = 0.0_real64
            end if
          else
            vint (i,j) = 0.0_real64
            vint2(i,j) = 0.0_real64
          end if
    6   continue
    7 continue
!      ENDIF
!
!   See if the appropriate angular data is available. If so,
!   then read the angular files and perform the calculation.
!   If not, the user can choose to save them or not
!
      doit_ob = 0
      doit_tb = 0
      call angdata(name,avail_ob,1,cache_status)
      if (cache_status /= 0) call fail_input_io(cache_status,50,"ris4: one-body angular cache "//trim(name))
      call angdata(name,avail_tb,2,cache_status)
      if (cache_status /= 0) call fail_input_io(cache_status,51,"ris4: two-body angular cache "//trim(name))
      if ((.not. avail_ob) .and. (.not. avail_tb)) then
        print *," Save ang. coefficients of one- and two-body op.?"
        yes = getyn ()
        print *
!C        IF(.NOT. YES) THEN
!C          PRINT *,' Save ang. coefficients of one-body op. ONLY?'
!C          YES = GETYN ()
!C          PRINT *
!C          IF(.NOT. YES) THEN
!C            PRINT *,' Save ang. coefficients of two-body op. ONLY?'
!C            YES = GETYN ()
!C            PRINT *
!C            IF(YES) DOIT_TB = 1
!C          ELSE
!C            DOIT_OB = 1
!C          ENDIF
!C        ELSE
       if(yes) then
          doit_ob = 1
          doit_tb = 1
        end if
      else if (.not. avail_ob) then
        print *," Save ang. coefficients of one-body op. ?"
        yes = getyn ()
        print *
        if(yes) doit_ob = 1
      else if (.not. avail_tb) then
        print *," Save ang. coefficients of two-body op. ?"
        yes = getyn ()
        print *
        if(yes) doit_tb = 1
      end if
      if (avail_ob) then
        j = index(name," ")
        open(unit=50,file = name(1:j-1)//".IOB",status="UNKNOWN"       &
            ,form="UNFORMATTED")
!    Read angular data from file and compute matrix elements
        if(yes2) then
           call densread_seltz(dint1,dint2,dint3, & ! JE ADD
                dint4,dint5,dint6, & ! JE ADD
                dint7,dint1vec,dens1vec,nrnuc) ! JE ADD
        else
           call densread(dint1,dint2,dint3,dint4,dint5,dint6,dint7)
        end if

        !        CALL DENSREAD(DINT1,DINT2,DINT3,             &   ! JE ADD
        !                      DINT4,DINT5,DINT6,             &   ! JE ADD
        !                      DINT7,DINT1VEC,DENS1VEC)           ! JE ADD
     else
!Check if the user wants to save one-body ang. coeff.
        if(doit_ob ==1) then
          j = index(name," ")
          open(unit=50,file = name(1:j-1)//".IOB",status="UNKNOWN"     &
            ,form="UNFORMATTED")
        end if
        if(yes2) then
           call densnew_seltz(doit_ob,dint1,dint2,dint3, & ! JE ADD
                dint4,dint5,dint6,dint7, & ! JE ADD
                dint1vec,dens1vec,nrnuc) ! JE ADD
        else
           call densnew(doit_ob,dint1,dint2,dint3,dint4, &
                dint5,dint6,dint7)
        end if

!        CALL NATORBNEW(NAME)                            ! JE ADD
!        CALL DENSNEW(DOIT_OB,DINT1,DINT2,DINT3,     &   ! JE ADD
!             DINT4,DINT5,DINT6,DINT7,               &   ! JE ADD
!             DINT1VEC,DENS1VEC)                         ! JE ADD
      end if
!
!
!
      if (avail_tb) then
        j = index(name," ")
        open(unit=51,file = name(1:j-1)//".ITB",status="UNKNOWN"       &
          ,form="UNFORMATTED")
        call smsread(vint,vint2)
      else
! Check if the user wants to save two-body ang. coeff.
        if(doit_tb==1) then
          j = index(name," ")
          open(unit=51,file = name(1:j-1)//".ITB",status="UNKNOWN"     &
             ,form="UNFORMATTED")
        end if
         call smsnew(doit_tb,vint,vint2)
      end if
!
!   Printouts
      !

      write (24,"(a24,i3)") "Number of eigenvalues: ",nvec
      write(24,*)
      write(24,*)
      write (24,309)
      do i = 1, nvec
         write (24,323) ivec(i),labj(iatjpo(i)),labp((iaspar(i)+3)/2), &
              hamiltonian_energy_shift+shifted_eigenvalues(i)
      end do
      write (24,308)
      do 18 i = 1,nvec
         write(24,314)
         write (24,313) ivec(i),labj(iatjpo(i)),labp((iaspar(i)+3)/2), &
              dens7(i),(dens2(i)-dens7(i)),dens2(i)
         hznmsu=dens7(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s
!         if (dabs(HzNMSu) .LE. (10**8)) THEN
!             WRITE (24,315) HzNMSu*(1.0D-06),
!     :      (DENS2(I)-DENS7(I))*AUMAMU*AUCM*CCMS*(1.0D-06),
!     :       DENS2(I)*AUMAMU*AUCM*CCMS*(1.0D-06)
!         else
         write (24,316) hznmsu*(1.0e-9_real64), &
              (dens2(i)-dens7(i))*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64), &
              dens2(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64)
!         endif
 18   continue

      write (24,302)
      do 14 i = 1,nvec
         write(24,314)
         write (24,313) ivec(i),labj(iatjpo(i)),labp((iaspar(i)+3)/2), &
              smsc1(i),smsc2(i)-smsc1(i),smsc2(i)
         hzsmsu=smsc1(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s
         if (dabs(hzsmsu) <= (10**8)) then
            write (24,315) hzsmsu*(1.0e-6_real64), &
                 (smsc2(i)-smsc1(i))*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-6_real64), &
                 smsc2(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-6_real64)
         else
            write (24,316) hzsmsu*(1.0e-9_real64), &
                 (smsc2(i)-smsc1(i))*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64), &
                 smsc2(i)*electron_mass_u*hartree_to_wavenumber_cm*speed_of_light_cm_per_s*(1.0e-9_real64)
         end if
 14   continue
      if(yes2) then
         write (24,307)
         do 16 i = 1,nvec
            call edensityfit(radial_grid,dens1vec(i,:),nuclear_charge,nuclear_model_parameters,nrnuc,fmat(:,i),rho(i),res(i))  ! JE ADD
            write (24,303) ivec(i),labj(iatjpo(i)), &
                 labp((iaspar(i)+3)/2),rho(i)
 16      continue
         write (24,317)
         do 17 i = 1,nvec
            !CALL EDENSITYFIT(R,DENS1VEC(I,:),Z,PARM,NRNUC,FMAT,RHO,RES)  ! JE ADD
            write (24,324) ivec(i),labj(iatjpo(i)), &
                 labp((iaspar(i)+3)/2),fmat(1,i),fmat(2,i), &
                 fmat(3,i),fmat(4,i),res(i)
 17      continue
         write (24,327)
         do 19 i = 1,nvec
            !CALL EDENSITYFIT(R,DENS1VEC(I,:),Z,PARM,NRNUC,FMAT,RHO,RES)  ! JE ADD
            write (24,334) ivec(i),labj(iatjpo(i)), &
                 labp((iaspar(i)+3)/2),fmat(5,i),fmat(6,i)
 19      continue
         write(24,*)
      else
         write (24,306)
         do 36 i = 1,nvec
            write (24,303) ivec(i),labj(iatjpo(i)), &
                 labp((iaspar(i)+3)/2),dens1(i)
 36      continue
      end if
!     WRITE (24,309)
!      DO 20 I = 1,NVEC
!     WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
!     :                  DENS3(I)
!     20 CONTINUE
!     WRITE (24,310)
!     DO 22 I = 1,NVEC
!     WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
!     :                  DENS4(I)
!     22 CONTINUE
!     WRITE (24,311)
!     DO 24 I = 1,NVEC
!         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
!     :                  DENS5(I)
!     24 CONTINUE
!     WRITE (24,312)
!     DO 26 I = 1,NVEC
!     WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
!     :                  DENS6(I)
!     26 CONTINUE
!
!   Dealloc
!
      call dalloc (smsc1,"SMSC1","RIS_CAL")
      call dalloc (smsc2,"SMSC2","RIS_CAL")
      call dalloc (dens1,"DENS1","RIS_CAL")
      call dalloc (dens2,"DENS2","RIS_CAL")
      call dalloc (dens3,"DENS3","RIS_CAL")
      call dalloc (dens4,"DENS4","RIS_CAL")
      call dalloc (dens5,"DENS5","RIS_CAL")
      call dalloc (dens6,"DENS6","RIS_CAL")
      call dalloc (dens7,"DENS7","RIS_CAL")

      write(*,*) "-------------------------------"
      write(*,*) "RIS_CAL: Execution Finished ..."
      write(*,*) "-------------------------------"

      call stoptime (ncount1, "RIS_CAL")
      return
!
301   format (//" CUTOFF set to ",1pd22.15)
302   format (//" Level  J Parity  Specific mass shift parameter")

!     Electron density values
303   format (1x,i3,5x,2a4,3x,d20.10)
304   format (1x,i3,5x,2a4,3x,2d20.10)

!     Electronic factor values
324   format (1x,i3,5x,2a4,3x,4d20.10,f10.4)

!     Field shift values
334   format (1x,i3,5x,2a4,3x,2d20.10)

      !     Electron density header
306   format (//" Electron density in atomic units" &
           //" Level  J Parity",8x,"DENS (a.u.)"/)

307   format (//" Level  J Parity  Electron density in atomic units" &
           //24x,"Dens. (a.u.)")

!     Electronic factors header
317   format (//" Level  J Parity  Field shift electronic factors and av&
           &erage point discrepancy in fit" &
           //24x, &
           "F0 (GHz/fm^2)",6x," F2 (GHz/fm^4)", &
           6x," F4 (GHz/fm^6)",6x," F6 (GHz/fm^8)", &
           6x," Disc. (per mille)")

!     Electronic field shift header
327   format (//" Level  J Parity  Field shift electronic factors (corre&
           &cted for varying density inside nucleus)" &
           //24x, &
           "F0VED0 (GHz/fm^2)   F0VED1 (GHz/fm^4)")

308   format (//" Level  J Parity  Normal mass shift parameter")

309   format (" Level  J Parity  Energy")

!     309 FORMAT (//' Radial expectationvalue'
!     :        //' Level  J Parity',8X,'<r> (a.u.)'/)
!     310 FORMAT (//' Radial expectationvalue'
!     :        //' Level  J Parity',8X,'<r2> (a.u.)'/)
!     311 FORMAT (//' Radial expectationvalue'
!     :        //' Level  J Parity',8X,'<r-1> (a.u.)'/)
!     312 FORMAT (//' Radial expectationvalue'
!     :        //' Level  J Parity',8X,'<r-2> (a.u.)'/)
313   format (1x,i3,5x,2a4,3x,3d20.10,"  (a.u.)")
333   format (1x,i3,5x,2a4,3x,3d20.10,"  (GHz u)")
323   format (1x,i3,5x,2a4,3x,d20.10,"  (a.u.)")
314   format (/,29x,"<K^1>",13x,"<K^2+K^3>",9x,"<K^1+K^2+K^3>")
315   format (20x,3d20.10 ,2x,"(MHz u)")
316   format (20x,3d20.10,2x,"(GHz u)")
      !
      return
    end subroutine ris_cal
