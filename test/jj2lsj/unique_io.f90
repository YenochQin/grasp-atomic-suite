program unique_io
    use jj2lsj_code, only: uniquelsj
    implicit none (type, external)
    open(57, file="input.lbl", status="old", form="formatted")
    open(72, file="unique.summary", status="replace")
    open(73, file="unique.lbl", status="replace")
    call uniquelsj
    close(57)
    close(72)
    close(73)
end program unique_io
