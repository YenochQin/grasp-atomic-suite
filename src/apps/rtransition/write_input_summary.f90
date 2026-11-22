module write_input_summary_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_input_summary
contains

    ! Generates the first part of  oscl92.sum  (on stream 24).
    !
    ! Written by Farid A. Parpia. Last revision: 28 Dec 1992
    subroutine write_input_summary(name, inpci, ilbl)
        ! Module dependencies
        use calculation_file_units, only: lsj_transition_output_unit, summary_output_unit
        use transition_lsj_labels, only: iopen_status1, iopen_status2
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: inpci
        integer :: ilbl
        character, intent(in) :: name(2)*24

        ! Local declarations
        integer :: i, j
        character(len=:), allocatable :: lsj_transition_filename

        i = index(name(1)," ")
        j = index(name(2)," ")
        if(ilbl == 0) then
            if (inpci == 0) then
                open(newunit=summary_output_unit,file=name(1)(1:i-1)//"."//name(2)(1:j-1)//".ct",form=&
                    "FORMATTED", status="UNKNOWN",position="asis")
            else
                open(newunit=summary_output_unit,file=name(1)(1:i-1)//"."//name(2)(1:j-1)//".t",form=&
                    "FORMATTED",status="UNKNOWN",position="asis")
            end if
        else if(ilbl == 1) then
            if(iopen_status1==0 .and. iopen_status2 ==0) then
                if (inpci == 0) then
                    lsj_transition_filename = name(1)(1:i-1)//"."//name(2)(1:j-1)//".ct.lsj"
                else
                    lsj_transition_filename = name(1)(1:i-1)//"."//name(2)(1:j-1)//".t.lsj"
                end if
                if (lsj_transition_output_unit == -1) then
                    open(newunit=lsj_transition_output_unit, file=lsj_transition_filename, &
                        form="FORMATTED", status="UNKNOWN")
                else
                    open(unit=lsj_transition_output_unit, file=lsj_transition_filename, &
                        form="FORMATTED", status="UNKNOWN")
                end if
            end if
        end if
        ! Write a timestamp in the summary header.
        !
        !    Write out the basic dimensions of the initial state electron cloud
        !
        ! //' electrons in the initial state cloud'
        ! //' relativistic CSFs'
        ! //' relativistic subshells.'
        !
        !    Write out the basic dimensions of the final state electron cloud
        !
        ! //' electrons in the final state cloud'
        ! //' relativistic CSFs'
        ! //' relativistic subshells.'
        !
        ! Report the first-order cutoff when the calculation restricts CSF interactions.
        !
        ! //' the zero-order space;'
        !
        ! Report the nuclear model and isotope parameters.
        !
        ! Report the enabled Hamiltonian corrections.
        !
        ! Report the radial grid parameters.
        !
        !    Write out the orbital properties
        !
        ! GAMAII(I),PFII(2,I),QFII(2,I),MFII(I)
        !
        ! GAMAFF(I),PFFF(2,I),QFFF(2,I),MFFF(I)
        !
        !    Write the list of eigenpair indices for the initial state
        !
        !    Write the list of eigenpair indices for the final state
        return

        return
    end subroutine write_input_summary
end module write_input_summary_interface
