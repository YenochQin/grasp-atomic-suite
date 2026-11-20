module try_allocate_sort_buffers_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: try_allocate_sort_buffers
contains

    !***********************************************************************
    !                                                                      *
    subroutine try_allocate_sort_buffers(n, irez)
        ! Module dependencies
        !                                                                      *
        !   Written by G. Gaigalas            Last revision: 27 October 2017   *
        !                                                                      *
        !                                                                      *
        !***********************************************************************
        !
        use, intrinsic :: iso_fortran_env, only: real64
        use array_memory_management, only: dalloc
        implicit none (type, external)

        ! Arguments
        integer, intent(in)  :: n
        integer, intent(out) :: irez

        ! Local declarations
        real(real64), pointer :: pr(:)
        integer, pointer :: pi(:)
        integer :: error
        !-----------------------------------------------
        irez = 1
        allocate(pr(n),stat=error)
        if (error /= 0) then
            irez = 0
            return
        end if
        allocate(pi(3*n),stat=error)
        if (error /= 0) then
            irez = 0
            call dalloc (pr,"PR","allocCheck")
            return
        end if
        call dalloc (pr,"PR","allocCheck")
        call dalloc (pi,"PI","allocCheck")
        return
    end subroutine try_allocate_sort_buffers
end module try_allocate_sort_buffers_interface
