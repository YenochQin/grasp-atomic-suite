
!======================================================================
! Written by Chongyang Chen, Fudan University, Shanghai, Aug 2024
!
! Report the memory (RAM) used. 
! Intrinsic-array payload sizes use standard Fortran storage inquiries.
! 
!======================================================================
      SUBROUTINE REPORT_MEM(ICALL, STRFRM)
      USE vast_kind_param
      USE suite_storage, ONLY: storage_bytes

      use symexpand_mod
      use symmatrix_mod
      use csfg_scf_C
      use csfg_tv_C

      use orb_C,   ONLY: IQA, NW
      use eigv_C,  ONLY: EVEC
 
! Sizes of variables in csfg_decide_C are all small
      ! use csfg_decide_C

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,       INTENT(IN) :: ICALL
      CHARACTER*(*), INTENT(IN) :: STRFRM
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL :: TOTMEM = 0.0D0
      !DATA TOTMEM /0.0D0/

! Use size information from vast_kind_param:
! integer, parameter :: byte_log     = selected_int_kind(2)
! integer, parameter :: short_log    = selected_int_kind(4)
! integer, parameter :: long_log     = selected_int_kind(18)
! integer, parameter :: byte         = selected_int_kind(2)
! integer, parameter :: short        = selected_int_kind(4)
! integer, parameter :: long         = selected_int_kind(18)
! integer, parameter :: double       = selected_real_kind(13)
! integer, parameter :: extended     = selected_real_kind(30)
! integer, parameter :: double_ext   = selected_real_kind(50)
! integer, parameter :: dble_complex = selected_real_kind(13)
! integer, parameter :: byte_log     = selected_int_kind(2)
! integer, parameter :: short_log    = selected_int_kind(4)

      REAL(DOUBLE), parameter :: F1GB  = 1.0D0/(1024 * 1024 * 1024)
      INTEGER      :: INT4 = 4, INT8 = 8, IBYTE = 1, ISHORT = 2
      INTEGER      :: IFLT = 4, IDBL = 8
!-----------------------------------------------
      REAL         :: TMPMEM, SUMM

      WRITE(*, *)"ReportMem called from " // TRIM(STRFRM)
      IF (ICALL == 1) THEN
        SUMM = 0.D0 
        ! MAP1(:,:)
        TMPMEM = F1GB*storage_bytes(MAP1)
        WRITE(*, *) 'Mem (GB) needed by MAP1     =', TMPMEM
        SUMM   = SUMM   + TMPMEM

        ! NTYPE
        TMPMEM = F1GB*storage_bytes(NTYPE)
        WRITE(*, *) 'Mem (GB) needed by NTPYE    =', TMPMEM
        SUMM   = SUMM   + TMPMEM

        ! ICGPDET, ICGPEND, IRGPDET
        TMPMEM = F1GB*(storage_bytes(ICGPDET)+storage_bytes(ICGPEND)+storage_bytes(IRGPDET))
        WRITE(*, *) 'Mem (GB) for part 1 arrays from mcp files: &
                     &ICGPDET, ICGPEND,  IRGPDET', TMPMEM
        TOTMEM = TOTMEM + TMPMEM

        ! IQA, this array could be deallocated before entering SCF.
        ! Of course, some modifications needed. 2024/08/18
        TMPMEM = F1GB*storage_bytes(IQA)
        WRITE(*, *) 'Mem (GB) needed by IQA_CSFG =', TMPMEM
        SUMM   = SUMM   + TMPMEM

        ! IQA_CSF
        TMPMEM = F1GB*storage_bytes(IQA_CSF)
        WRITE(*, *) 'Mem (GB) needed by IQA_CSF  =', TMPMEM
        SUMM   = SUMM   + TMPMEM

        TOTMEM = TOTMEM + SUMM 
        WRITE(*, *) 'Mem (GB) needed for SCF    (Maintained) :', TOTMEM

      ELSEIF(ICALL == 2) THEN
        SUMM   = 0.D0 

        TMPMEM = F1GB*storage_bytes(LABTH)
        WRITE(*, *) 'Mem (GB) needed by LABTH    =', TMPMEM
        SUMM   = SUMM   + TMPMEM
        
        TMPMEM = F1GB*storage_bytes(TCOEFH)
        WRITE(*, *) 'Mem (GB) needed by TCOEFH   =', TMPMEM
        SUMM   = SUMM   + TMPMEM

        TMPMEM = F1GB*storage_bytes(NVHP  )
        WRITE(*, *) 'Mem (GB) needed by NVHP     =', TMPMEM
        SUMM   = SUMM   + TMPMEM

        TMPMEM = F1GB*storage_bytes(LABVKH)
        WRITE(*, *) 'Mem (GB) needed by LABVKH   =', TMPMEM
        SUMM   = SUMM   + TMPMEM

        TMPMEM = F1GB*storage_bytes(VCOEFH)
        WRITE(*, *) 'Mem (GB) needed by VCOEFH   =', TMPMEM
        SUMM   = SUMM   + TMPMEM

        WRITE(*, *) "LABTH, TCOEFH,  NVHP, LABVKH, VCOEFH    :", SUMM

        TMPMEM = F1GB*storage_bytes(ICOLIDX)
        WRITE(*, *) 'Mem (GB) needed by ICOLIDX              :', TMPMEM
        SUMM   = SUMM   + TMPMEM

        TMPMEM = F1GB*storage_bytes(EVEC)
        WRITE(*, *) 'Mem (GB) needed by EVEC                 :', TMPMEM
        SUMM   = SUMM   + TMPMEM

        TOTMEM = TOTMEM + SUMM
        WRITE(*, *) "Total Mem needed by now    (Maintained) :", TOTMEM

      ELSEIF (ICALL == 3) THEN
        SUMM   = 0.D0
        
        TMPMEM = F1GB*(storage_bytes(NDA) + storage_bytes(DA) + storage_bytes(NXA) +        &
                       storage_bytes(XA)  + storage_bytes(NYA)+ storage_bytes(YA))
        WRITE(*, *) "NDA, DA, NXA, XA, NYA, YA  (Maintained) :", TMPMEM
        SUMM   = SUMM   + TMPMEM

        WRITE(*, *) "NXAopt,XAopt,NYAopt,YAopt (Deallocated) :", &
                     REAL(F1GB*(storage_bytes(NXAOPT) + storage_bytes(XAOPT)+ &
                                storage_bytes(NYAOPT) + storage_bytes(YAOPT)))
        WRITE(*, *) "NXAWRK, XAWRK,NYAWRK,YAWRK (Maintained) :", &
                     REAL(F1GB*(storage_bytes(NXAWRK) + storage_bytes(XAWRK)))

        TMPMEM = F1GB*(storage_bytes(NDAOPT) + storage_bytes(DAOPT) + storage_bytes(NXAWRK)+&
                       storage_bytes(XAWRK)  + storage_bytes(NYAOPT)+ storage_bytes(YAOPT) +&
                       storage_bytes(TMPXYA))
        WRITE(*, *) "NDAOPT,DAOPT, NXAWRK,XAWRK, NYAWRK,YAWRK:", TMPMEM
        SUMM   = SUMM + TMPMEM

        TOTMEM = TOTMEM + SUMM
        WRITE(*, *) "Total Mem needed by now    (Maintained) :", TOTMEM

      ELSEIF (ICALL == 4) THEN
        SUMM   = 0.D0
      ELSE
        STOP "Program stopped at REPORT_MEM ..."
      ENDIF

      RETURN
      END SUBROUTINE REPORT_MEM
