module load_final_lsj_labels_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: load_final_lsj_labels
contains

    ! Read the final state's LSJ labels.
    !
    ! Written by G. Gaigalas,
    !    NIST                                                  May 2011
    !
    ! Translated by Charlotte Froese Fischer
    !                        Gediminas Gaigalas  10/05/17
    subroutine load_final_lsj_labels(name, from_ci)
        !  Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use calculation_state_files, only: lsj_label_input_path
        use array_memory_management, only: alloc
        use transition_lsj_labels, only: rlev_ener_2, string_csf2, nvectotf, iopen_status2, lev_pos_2, lev_j_2, lev_par_2
        implicit none (type, external)

        ! Arguments
        character(*), intent(in) :: name
        logical, intent(in) :: from_ci

        ! Local declarations
        character(:), allocatable :: label_path
        integer :: lsj_label_input_unit, icount, io_status, itest
        real(real64) :: weights
        character(15) :: record

        label_path = lsj_label_input_path(name, from_ci)
        open(newunit=lsj_label_input_unit, file=label_path, form="formatted", status="old", &
            action="read", iostat=iopen_status2)
        if (iopen_status2 /= 0) return
        read (lsj_label_input_unit,"(1A15)",iostat = io_status) record
        if (io_status /= 0) then
            print *, "Cannot read LSJ labels: ", label_path
            close (lsj_label_input_unit)
            stop
        else
            call alloc (lev_pos_2,  nvectotf,"Lev_POS_2",  "load_final_lsj_labels")
            call alloc (lev_j_2,    nvectotf,"Lev_J_2",    "load_final_lsj_labels")
            call alloc (lev_par_2,  nvectotf,"Lev_Par_2",  "load_final_lsj_labels")
            call alloc (rlev_ener_2,nvectotf,"RLev_ENER_2","load_final_lsj_labels")
            ! CPJ        CALL ALLOC (string_CSF2,NVECTOTF,'string_CSF2','LDLBL2')
            allocate(string_csf2(nvectotf))
            icount = 1
            do
                read (lsj_label_input_unit,"(I3,1X,A4,5X,A1,8X,F16.9)",iostat = io_status)     &
                    lev_pos_2(icount),lev_j_2(icount),lev_par_2(icount),    &
                    rlev_ener_2(icount)
                if (io_status /= 0) exit
                read (lsj_label_input_unit,"(7X,F12.8,17X,A)") weights,string_csf2(icount)

                ! Skip additional composition records before the next level.
                do
                    read (lsj_label_input_unit,"(I3)",iostat = io_status) itest
                    if (io_status /= 0) exit
                    if (itest /= 0) exit
                end do
                if (io_status /= 0) exit
                backspace lsj_label_input_unit
                icount = icount + 1
            end do
        end if
        close (lsj_label_input_unit)
        return
    end subroutine load_final_lsj_labels
end module load_final_lsj_labels_interface
