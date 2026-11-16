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
      use vast_kind_param,  only: double
      use parameter_def,    only: nnnw
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
      use ris_c
      use eigv_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
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
      use densnew_i
      use smsread_i
      use smsnew_i
      use edensityfit_i
      use densnew_seltz_i, only: densnew_seltz
      use densread_seltz_i, only: densread_seltz
      use starttime_i, only: starttime
      use stoptime_i, only: stoptime
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character*24, intent(in) :: name
!      REAL(DOUBLE), INTENT(IN) :: DR2
!      INTEGER, INTENT(IN)      :: NOPAR
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double), dimension(:,:,:), pointer :: dint1vec
      real(double), dimension(:,:), pointer   :: dens1vec
      real(double), dimension(:), pointer     :: densfit
      real(double), dimension(:,:), pointer     :: fmat
      real(double), dimension(:), pointer     :: rho
      real(double), dimension(:), pointer     :: res


      real(double), dimension(nnnw)      :: tshell
      real(double), dimension(nnnw,nnnw) :: vint, vint2
      real(double), dimension(nnnw,nnnw) :: dint1, dint2, dint3, &
                                            dint4, dint5, dint6, &
                                            dint7
      real(double) :: au2fm, rcre, hzsmsu, hznmsu, fo90
      character*11 :: cnum
      character*2  :: ck
      logical :: vsh, nucde, smssh, yes, yes2, avail_tb, avail_ob
      integer :: i, j, l, ncount1, doit_ob, doit_tb, nrnuc
!-----------------------------------------------
!

!
!   Matrix elements smaller than CUTOFF are not accumulated
!
!      PARAMETER (CUTOFF = 1.0D-10)

      write(*,*) '-------------------------------'
      write(*,*) 'RIS_CAL: Execution Begins ...'
      write(*,*) '-------------------------------'

!     Determine point end point NRNUC for density fit
      fo90 = parm(1)+2.d0*log(3.d0)*parm(2)
      !PARM(1) - Fermi parameter c
      !PARM(2) - Fermi parameter a
      do i=1,nnnp
         if(r(i+1).gt.fo90.and.r(i).le.fo90) then
            nrnuc = i+1
            exit
         end if
      end do


!   Allocate memory

      allocate( dens1vec(nvec,nrnuc),dint1vec(nnnw,nnnw,nrnuc) )    ! JE ADD
      allocate( densfit(nrnuc) )                                   ! JE ADD

      dint1vec(:,:,:) = 0.d0                                      ! JE ADD
      dens1vec(:,:) = 0.d0                                        ! JE ADD

      allocate( fmat(6,nvec),rho(nvec),res(nvec) )                ! JE ADD

!
!   Allocate storage for local arrays
!
      call alloc (smsc1, nvec,'SMSC1', 'RIS_CAL')
      call alloc (smsc2, nvec,'SMSC2', 'RIS_CAL')
      call alloc (dens1,nvec,'DENS1','RIS_CAL')
      call alloc (dens2,nvec,'DENS2','RIS_CAL')
      call alloc (dens3,nvec,'DENS3','RIS_CAL')
      call alloc (dens4,nvec,'DENS4','RIS_CAL')
      call alloc (dens5,nvec,'DENS5','RIS_CAL')
      call alloc (dens6,nvec,'DENS6','RIS_CAL')
      call alloc (dens7,nvec,'DENS7','RIS_CAL')

      call starttime (ncount1, 'RIS_CAL')

!
!   Initialise
!
      do i = 1,nvec
         smsc1(i)  = 0.0d00
         smsc2(i)  = 0.0d00
         dens1(i) = 0.0d00
         dens2(i) = 0.0d00
         dens3(i) = 0.0d00
         dens4(i) = 0.0d00
         dens5(i) = 0.0d00
         dens6(i) = 0.0d00
         dens7(i) = 0.0d00
      enddo

!      VSH   = .TRUE.
!      SMSSH = .TRUE.
!
!   Calculate all integrals needed for the volume shift calc.
!
      !      IF (VSH) THEN

      write(*,*) ' Compute higher order field shift electronic factors?'
      yes2 = getyn ()

      do 5 i = 1,nw
        do 4 j = 1,nw
          if (nak(i).eq.nak(j)) then
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
            dint1(i,j) = 0.0d00
            dint1vec(i,j,:) = 0.0d00
            dint2(i,j) = 0.0d00
            dint3(i,j) = 0.0d00
            dint4(i,j) = 0.0d00
            dint5(i,j) = 0.0d00
            dint6(i,j) = 0.0d00
            dint7(i,j) = 0.0d00
          endif
    4   continue
    5 continue
!      ENDIF
!
!   Calculate and save the Vinti integrals
!
!      IF (SMSSH) THEN
      do 7 i = 1,nw
        do 6 j = 1,nw
          if (i.ne.j) then
            rcre = cre(nak(i),1,nak(j))
            if (dabs(rcre) .gt. cutoff) then
              vint  (i,j) = vinti(i,j)
              vint2(i,j) = vint(i,j)                                   &
                         + rint_sms2(i,j)/rcre                         &
                         + rint_sms3(i,j)
            else
             vint (i,j) = 0.0d00
             vint2(i,j) = 0.0d00
            endif
          else
            vint (i,j) = 0.0d00
            vint2(i,j) = 0.0d00
          endif
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
      call angdata(name,avail_ob,1)
      call angdata(name,avail_tb,2)
      if ((.not. avail_ob) .and. (.not. avail_tb)) then
        print *,' Save ang. coefficients of one- and two-body op.?'
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
        endif
      elseif (.not. avail_ob) then
        print *,' Save ang. coefficients of one-body op. ?'
        yes = getyn ()
        print *
        if(yes) doit_ob = 1
      elseif (.not. avail_tb) then
        print *,' Save ang. coefficients of two-body op. ?'
        yes = getyn ()
        print *
        if(yes) doit_tb = 1
      endif
      if (avail_ob) then
        j = index(name,' ')
        open(unit=50,file = name(1:j-1)//'.IOB',status='UNKNOWN'       &
            ,form='UNFORMATTED')
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
        if(doit_ob .eq.1) then
          j = index(name,' ')
          open(unit=50,file = name(1:j-1)//'.IOB',status='UNKNOWN'     &
            ,form='UNFORMATTED')
        endif
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
      endif
!
!
!
      if (avail_tb) then
        j = index(name,' ')
        open(unit=51,file = name(1:j-1)//'.ITB',status='UNKNOWN'       &
          ,form='UNFORMATTED')
        call smsread(vint,vint2)
      else
! Check if the user wants to save two-body ang. coeff.
        if(doit_tb.eq.1) then
          j = index(name,' ')
          open(unit=51,file = name(1:j-1)//'.ITB',status='UNKNOWN'     &
             ,form='UNFORMATTED')
        endif
         call smsnew(doit_tb,vint,vint2)
      endif
!
!   Printouts
      !

      write (24,'(a24,i3)') 'Number of eigenvalues: ',nvec
      write(24,*)
      write(24,*)
      write (24,309)
      do i = 1, nvec
         write (24,323) ivec(i),labj(iatjpo(i)),labp((iaspar(i)+3)/2), &
              eav+eval(i)
      end do
      write (24,308)
      do 18 i = 1,nvec
         write(24,314)
         write (24,313) ivec(i),labj(iatjpo(i)),labp((iaspar(i)+3)/2), &
              dens7(i),(dens2(i)-dens7(i)),dens2(i)
         hznmsu=dens7(i)*aumamu*aucm*ccms
!         if (dabs(HzNMSu) .LE. (10**8)) THEN
!             WRITE (24,315) HzNMSu*(1.0D-06),
!     :      (DENS2(I)-DENS7(I))*AUMAMU*AUCM*CCMS*(1.0D-06),
!     :       DENS2(I)*AUMAMU*AUCM*CCMS*(1.0D-06)
!         else
         write (24,316) hznmsu*(1.0d-09), &
              (dens2(i)-dens7(i))*aumamu*aucm*ccms*(1.0d-09), &
              dens2(i)*aumamu*aucm*ccms*(1.0d-09)
!         endif
 18   continue

      write (24,302)
      do 14 i = 1,nvec
         write(24,314)
         write (24,313) ivec(i),labj(iatjpo(i)),labp((iaspar(i)+3)/2), &
              smsc1(i),smsc2(i)-smsc1(i),smsc2(i)
         hzsmsu=smsc1(i)*aumamu*aucm*ccms
         if (dabs(hzsmsu) .le. (10**8)) then
            write (24,315) hzsmsu*(1.0d-06), &
                 (smsc2(i)-smsc1(i))*aumamu*aucm*ccms*(1.0d-06), &
                 smsc2(i)*aumamu*aucm*ccms*(1.0d-06)
         else
            write (24,316) hzsmsu*(1.0d-09), &
                 (smsc2(i)-smsc1(i))*aumamu*aucm*ccms*(1.0d-09), &
                 smsc2(i)*aumamu*aucm*ccms*(1.0d-09)
         endif
 14   continue
      if(yes2) then
         write (24,307)
         do 16 i = 1,nvec
            call edensityfit(r,dens1vec(i,:),z,parm,nrnuc,fmat(:,i),rho(i),res(i))  ! JE ADD
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
      call dalloc (smsc1,'SMSC1','RIS_CAL')
      call dalloc (smsc2,'SMSC2','RIS_CAL')
      call dalloc (dens1,'DENS1','RIS_CAL')
      call dalloc (dens2,'DENS2','RIS_CAL')
      call dalloc (dens3,'DENS3','RIS_CAL')
      call dalloc (dens4,'DENS4','RIS_CAL')
      call dalloc (dens5,'DENS5','RIS_CAL')
      call dalloc (dens6,'DENS6','RIS_CAL')
      call dalloc (dens7,'DENS7','RIS_CAL')

      write(*,*) '-------------------------------'
      write(*,*) 'RIS_CAL: Execution Finished ...'
      write(*,*) '-------------------------------'

      call stoptime (ncount1, 'RIS_CAL')
      return
!
301   format (//' CUTOFF set to ',1pd22.15)
302   format (//' Level  J Parity  Specific mass shift parameter')

!     Electron density values
303   format (1x,i3,5x,2a4,3x,d20.10)
304   format (1x,i3,5x,2a4,3x,2d20.10)

!     Electronic factor values
324   format (1x,i3,5x,2a4,3x,4d20.10,f10.4)

!     Field shift values
334   format (1x,i3,5x,2a4,3x,2d20.10)

      !     Electron density header
306   format (//' Electron density in atomic units' &
           //' Level  J Parity',8x,'DENS (a.u.)'/)

307   format (//' Level  J Parity  Electron density in atomic units' &
           //24x,'Dens. (a.u.)')

!     Electronic factors header
317   format (//' Level  J Parity  Field shift electronic factors and av&
           &erage point discrepancy in fit' &
           //24x, &
           'F0 (GHz/fm^2)',6x,' F2 (GHz/fm^4)', &
           6x,' F4 (GHz/fm^6)',6x,' F6 (GHz/fm^8)', &
           6x,' Disc. (per mille)')

!     Electronic field shift header
327   format (//' Level  J Parity  Field shift electronic factors (corre&
           &cted for varying density inside nucleus)' &
           //24x, &
           'F0VED0 (GHz/fm^2)   F0VED1 (GHz/fm^4)')

308   format (//' Level  J Parity  Normal mass shift parameter')

309   format (' Level  J Parity  Energy')

!     309 FORMAT (//' Radial expectationvalue'
!     :        //' Level  J Parity',8X,'<r> (a.u.)'/)
!     310 FORMAT (//' Radial expectationvalue'
!     :        //' Level  J Parity',8X,'<r2> (a.u.)'/)
!     311 FORMAT (//' Radial expectationvalue'
!     :        //' Level  J Parity',8X,'<r-1> (a.u.)'/)
!     312 FORMAT (//' Radial expectationvalue'
!     :        //' Level  J Parity',8X,'<r-2> (a.u.)'/)
313   format (1x,i3,5x,2a4,3x,3d20.10,'  (a.u.)')
333   format (1x,i3,5x,2a4,3x,3d20.10,'  (GHz u)')
323   format (1x,i3,5x,2a4,3x,d20.10,'  (a.u.)')
314   format (/,29x,'<K^1>',13x,'<K^2+K^3>',9x,'<K^1+K^2+K^3>')
315   format (20x,3d20.10 ,2x,'(MHz u)')
316   format (20x,3d20.10,2x,'(GHz u)')
      !
      return
    end subroutine ris_cal
