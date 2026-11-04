! Stable LS ranking. The input coefficients are never reordered or truncated.
module jj2lsj_postprocess
    use vast_kind_param, only: DOUBLE
    use ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: ls_ranking, rank_ls, post_worker_count
    type :: ls_ranking
        integer, allocatable :: order(:)
        real(DOUBLE) :: norm = 0.0_DOUBLE
        real(DOUBLE) :: seconds = 0.0_DOUBLE
    end type
contains
    integer function post_worker_count(levels) result(workers)
!$      use omp_lib, only: omp_get_max_threads
        integer, intent(in) :: levels
        character(len=16) :: setting
        integer :: status
        workers=1
        ! Keep the conversion team's size between active regions: its SIXJ
        ! threadprivate caches must not depend on a changing team size.
!$      if (levels>1) workers=omp_get_max_threads()
        call get_environment_variable('JJ2LSJ_SERIAL',setting,status=status)
        if (status==0) then
            if (trim(setting)=='1') workers=1
        end if
        call get_environment_variable('JJ2LSJ_POST_SERIAL',setting,status=status)
        if (status==0) then
            if (trim(setting)=='1') workers=1
        end if
    end function

    recursive subroutine rank_ls(values, ranking)
        real(DOUBLE), intent(in) :: values(:)
        type(ls_ranking), intent(inout) :: ranking
        real(DOUBLE), allocatable :: squares(:)
        integer, allocatable :: scratch(:)
        integer :: n, i, width, first, middle, last, left, right, out
        n=size(values)
        if (allocated(ranking%order)) deallocate(ranking%order)
        allocate(ranking%order(n),squares(n),scratch(n))
        ranking%norm=0.0_DOUBLE
        do i=1,n
            if (.not. ieee_is_finite(values(i))) error stop 'jj2lsj: nonfinite LS coefficient'
            squares(i)=values(i)*values(i)
            ! Match the legacy accumulation in original LS order.
            ranking%norm=ranking%norm+squares(i)
            ranking%order(i)=i
        end do
        width=1
        do while (width<n)
            first=1
            do while (first<=n)
                middle=first+min(width,n-first+1)-1
                last=middle+min(width,n-middle)
                left=first; right=middle+1
                do out=first,last
                    if (left>middle) then
                        scratch(out)=ranking%order(right); right=right+1
                    else if (right>last) then
                        scratch(out)=ranking%order(left); left=left+1
                    else if (squares(ranking%order(right))>squares(ranking%order(left))) then
                        scratch(out)=ranking%order(right); right=right+1
                    else
                        ! Strict > reproduces insertion-sort ties, including +/-0.
                        scratch(out)=ranking%order(left); left=left+1
                    end if
                end do
                if (last==n) exit
                first=last+1
            end do
            ranking%order=scratch
            if (width>n/2) exit
            width=width*2
        end do
    end subroutine
end module
