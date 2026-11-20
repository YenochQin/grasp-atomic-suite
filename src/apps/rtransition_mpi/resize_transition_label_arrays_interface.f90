module resize_transition_label_arrays_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine resize_transition_label_arrays(ibeg, ilab, ilast, ileft, iptcs, irigh, lblin, &
            lldim, imode)
            implicit none (type, external)

            ! Arguments
            integer, intent(inout) :: lldim
            integer, intent(in) :: imode
            integer, pointer :: ibeg(:), ilab(:), ilast(:), ileft(:), &
                iptcs(:), irigh(:), lblin(:)
            !...This routine performs I/O.
        end subroutine resize_transition_label_arrays
    end interface
end module resize_transition_label_arrays_interface
