!***********************************************************************
!                                                                      *
      subroutine ldlbl2 (name)
!                                                                      *
!   Open, check and load data from the  .lsj.lbl   file of the         *
!   inital state.                                                      *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC.                                        *
!                                                                      *
!   Written by G. Gaigalas,                                            *
!   NIST                                                  May 2011     *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use array_memory_management, only: alloc
      use transition_lsj_labels, only: rlev_ener_1, rlev_ener_2, string_csf1, string_csf2, nvectoti, nvectotf, &
          iopen_status1, iopen_status2, lev_pos_1, lev_j_1, lev_par_1, lev_pos_2, lev_j_2, lev_par_2
      implicit none (type, external)
      integer :: lsj_label_input_unit
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character, intent(in) :: name*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64) :: weights
      character    :: record*15
      integer      :: icount, io_status, itest, j
!-----------------------------------------------
!
!   Check the first record of the file; if not as expected, try again
!
      j = index(name," ")
      open (newunit=lsj_label_input_unit,file=name(1:j-1)//".lsj.lbl",form="FORMATTED",  &
           status="OLD",iostat = iopen_status2)
      if(iopen_status2 == 0) then
         read (lsj_label_input_unit,"(1A15)",iostat = io_status) record
         if (io_status /= 0) then
            print *, "Not a i *.lsj.lbl  File;"
            close (lsj_label_input_unit)
            stop
         else
            call alloc (lev_pos_2,  nvectotf,"Lev_POS_2",  "LDLBL2")
            call alloc (lev_j_2,    nvectotf,"Lev_J_2",    "LDLBL2")
            call alloc (lev_par_2,  nvectotf,"Lev_Par_2",  "LDLBL2")
            call alloc (rlev_ener_2,nvectotf,"RLev_ENER_2","LDLBL2")
!CPJ        CALL ALLOC (string_CSF2,NVECTOTF,'string_CSF2','LDLBL2')
            allocate(string_csf2(nvectotf))
!
            icount = 1
            read (lsj_label_input_unit,"(1X,I2,1X,A4,5X,A1,8X,F16.9)",iostat = io_status)     &
              lev_pos_2(icount),lev_j_2(icount),lev_par_2(icount),    &
              rlev_ener_2(icount)
            if (io_status /= 0) goto 1
!
            read (lsj_label_input_unit,"(7X,F12.8,17X,A)") weights,string_csf2(icount)
!
    2       read (lsj_label_input_unit,"(1X,I2)",iostat = io_status) itest
            if (io_status /= 0) goto 1
            if (itest == 0) goto 2
            backspace lsj_label_input_unit
            icount = icount + 1
            read (lsj_label_input_unit,"(1X,I2,1X,A4,5X,A1,8X,F16.9)",iostat = io_status)     &
              lev_pos_2(icount),lev_j_2(icount),lev_par_2(icount),    &
              rlev_ener_2(icount)
            if (io_status /= 0) goto 1
            read (lsj_label_input_unit,"(7X,F12.8,17X,A)") weights,string_csf2(icount)
            goto 2
         end if
      end if
    1 continue
      close (lsj_label_input_unit)
      return
      end subroutine ldlbl2
