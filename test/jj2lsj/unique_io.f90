program unique_io
    use jj2lsj_conversion, only: assign_unique_lsj_labels
    implicit none (type, external)
    integer :: lsj_label_unit, unique_summary_unit, unique_label_unit
    open(newunit=lsj_label_unit, file="input.lbl", status="old", form="formatted")
    open(newunit=unique_summary_unit, file="unique.summary", status="replace")
    open(newunit=unique_label_unit, file="unique.lbl", status="replace")
    call assign_unique_lsj_labels(lsj_label_unit, unique_summary_unit, unique_label_unit)
    close(lsj_label_unit)
    close(unique_summary_unit)
    close(unique_label_unit)
end program unique_io
