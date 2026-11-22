module load_initial_lsj_labels_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: load_initial_lsj_labels
contains

    !    Open, check and load data from the  .lsj.lbl   file of the
    !    inital state.
    !
    ! Written by G. Gaigalas,
    !    NIST                                                  May 2011
    !
    ! Translated by Charlotte Froese Fischer
    !                        Gediminas Gaigalas  10/05/17
    subroutine load_initial_lsj_labels (name)
        !  Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use array_memory_management, only: alloc
        use transition_lsj_labels, only: rlev_ener_1, string_csf1, nvectoti, iopen_status1, lev_pos_1, lev_j_1, lev_par_1
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: lsj_label_input_unit
        character(len=*), intent(in) :: name
        real(real64) :: weights
        character    :: record*15
        integer      :: icount, io_status, itest, j
        ! Check the first record of the file; if not as expected, try again
        j = index(name," ")
        open (newunit=lsj_label_input_unit,file=name(1:j-1)//".lsj.lbl",form="FORMATTED",  &
            status="OLD",iostat = iopen_status1)
        if(iopen_status1 == 0) then
            read (lsj_label_input_unit,"(1A15)",iostat = io_status) record
            if (io_status /= 0) then
                print *, "Not a i *.lsj.lbl  File;"
                close (lsj_label_input_unit)
                stop
            else
                call alloc (lev_pos_1,  nvectoti,"Lev_POS_1",  "load_initial_lsj_labels")
                call alloc (lev_j_1,    nvectoti,"Lev_J_1",    "load_initial_lsj_labels")
                call alloc (lev_par_1,  nvectoti,"Lev_Par_1",  "load_initial_lsj_labels")
                call alloc (rlev_ener_1,nvectoti,"RLev_ENER_1","load_initial_lsj_labels")
                ! CPJ        CALL ALLOC (string_CSF1,NVECTOTI,'string_CSF1','LDLBL1')
                allocate(string_csf1(nvectoti))
                icount = 1
                do
                    read (lsj_label_input_unit,"(1X,I2,1X,A4,5X,A1,8X,F16.9)",iostat = io_status)     &
                        lev_pos_1(icount),lev_j_1(icount),lev_par_1(icount),    &
                        rlev_ener_1(icount)
                    if (io_status /= 0) exit
                    read (lsj_label_input_unit,"(7X,F12.8,17X,A)") weights,string_csf1(icount)

                    ! Skip additional composition records before the next level.
                    do
                        read (lsj_label_input_unit,"(1X,I2)",iostat = io_status) itest
                        if (io_status /= 0) exit
                        if (itest /= 0) exit
                    end do
                    if (io_status /= 0) exit
                    backspace lsj_label_input_unit
                    icount = icount + 1
                end do
            end if
        end if
        close (lsj_label_input_unit)
        return
    end subroutine load_initial_lsj_labels
end module load_initial_lsj_labels_interface
