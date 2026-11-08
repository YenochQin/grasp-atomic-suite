! Exercise the real parser, LS basis generator and transformation, then expose
! every signed double for bitwise comparison, including the zero truncated tail.
program coefficients
    use jj2lsj_code
    use orb_C, only: NCF
    use prnt_C, only: NVEC
    use eigv_C, only: EVEC
    use blk_C, only: NBLOCK, NCFINBLK
    use factt_I
    use setiso_I
    use vast_kind_param, only: LONG
    implicit none
    integer :: counts(Blocks_number)
    integer, allocatable :: levels(:,:)
    integer :: block, first, last, lev, level, i, outc, outj, unique, outt
    integer, allocatable :: mask(:)
    integer(LONG) :: stride
    real(DOUBLE) :: threshold, weight

    call factt
    call setiso('isodata')
    allocate(levels(Blocks_number,Vectors_number))
    call inscreen(threshold,levels,counts,outc,outj,unique,outt)
    allocate(mask(NCF))
    stride = NCF
    last = 0
    open(80,file='coefficients.raw',status='replace',access='stream',form='unformatted')
    do block = 1,NBLOCK
        first = last+1
        last = first+NCFINBLK(block)-1
        if (counts(block) == 0) cycle
        mask = 0
        do i = first,last
            weight = 0.0_DOUBLE
            do lev = 1,counts(block)
                level = levels(block,lev)
                weight = weight+abs(EVEC(i+(level-1)*stride))
            end do
            if (abs(weight) >= abs(EPSNEW)) mask(i) = 1
        end do
        call setLS(mask,first,last)
        allocate(asf_set_LS%asf(NVEC))
        do i = 1,NVEC
            allocate(asf_set_LS%asf(i)%eigenvector(asf_set_LS%csf_set_LS%nocsf))
            asf_set_LS%asf(i)%eigenvector = 0.0_DOUBLE
        end do
        call asf2ls(1,mask,counts(block),block,levels,first,last,outt)
        write(80) block, counts(block), asf_set_LS%csf_set_LS%nocsf
        do lev = 1,counts(block)
            write(80) asf_set_LS%asf(levels(block,lev))%eigenvector
        end do
        call dallocASFLS(asf_set_LS)
    end do
    close(80)
    print *, 'Coefficient dump complete'
end program coefficients
