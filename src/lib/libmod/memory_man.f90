module memory_man
   use vast_kind_param
   interface alloc
      module procedure alloc_0i, alloc_1i, alloc_2i, alloc_3i,     &
                                 alloc_1il,                        &
                       alloc_1r, alloc_2r,                         &
                       alloc_1rl,                                  &
                       alloc_0b, alloc_1b, alloc_2b, alloc_3b,     &
                       alloc_1c
   end interface

   interface ralloc
      module procedure ralloc_0i, ralloc_1i, ralloc_2i, ralloc_3i, &
                       ralloc_1r, ralloc_2r,                       &
                       ralloc_0b, ralloc_1b, ralloc_2b, ralloc_3b, &
                       ralloc_1c
   end interface ralloc

   interface dalloc
      module procedure dalloc_1i, dalloc_2i, dalloc_3i,            &
                       dalloc_1r, dalloc_2r,                       &
                       dalloc_1b, dalloc_2b, dalloc_3b,            &
                       dalloc_1c
   end interface

   contains
   subroutine alloc_0i(p, first, last, var, sub )
     implicit none
     integer, pointer, dimension(:) :: p
     integer, intent(in) :: first, last
     character(len=*), intent(in) :: var, sub
     integer :: error

!    write(777,*) "0i first, last, var, sub",first, last, var,"  ", sub
     nullify(p)
     if (associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        stop
     else
        allocate(p(first:last),stat=error)
        if (error/= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,2I15)') 'Number of elements requested: ',first, last
        stop
        end if
     end if
   end subroutine alloc_0i

   subroutine alloc_1i(p, n, var, sub)
     implicit none
     integer, pointer, dimension(:) :: p
     integer, intent(in) :: n
     character(len=*), intent(in) :: var, sub
     integer :: error

!    write(777,*) "1i n, var, sub",n, var,"  ", sub
     nullify(p)
     if (associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        stop
     else
        allocate(p(n),stat=error)
        if (error/= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15)') 'Number of elements requested: ',n
        stop
        end if
     end if
   end subroutine alloc_1i

   subroutine alloc_2i(p, n1, n2, var, sub)
     implicit none
     integer, pointer, dimension(:,:) :: p
     integer, intent(in) :: n1, n2
     character(len=*), intent(in) :: var, sub
     integer :: error

!    write(777,*) "2i n1, n2, var, sub",n1, n2, var,"  ", sub
     nullify(p)
     if (associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        stop
     else
        allocate(p(n1,n2),stat=error)
        if (error /= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15,I15)') 'Number of elements requested: ',n1,n2
        stop
        end if
     end if
   end subroutine alloc_2i

   subroutine alloc_3i(p, n1, n2, n3, var, sub)
     implicit none
     integer, pointer, dimension(:,:,:) :: p
     integer, intent(in) :: n1, n2, n3
     character(len=*), intent(in) :: var, sub
     integer :: error

!    write(777,*) "3i n1, n2, n3, var, sub",n1, n2, n3, var,"  ", sub
     nullify(p)
     if (associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        stop
     else
        allocate(p(n1,n2,n3),stat=error)
        if (error /= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,3I15)') 'Number of elements requested: ',n1,n2,n3
        stop
        end if
     end if
   end subroutine alloc_3i

   subroutine alloc_1il(p, n, var, sub)
     implicit none
     integer, pointer, dimension(:) :: p
     integer(long), intent(in) :: n
     character(len=*), intent(in) :: var, sub
     integer :: error

!    write(777,*) "1iL n, var, sub",n, var,"  ", sub
     nullify(p)
     if (associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        stop
     else
        allocate(p(n),stat=error)
        if (error/= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15)') 'Number of elements requested: ',n
        stop
        end if
     end if
   end subroutine alloc_1il

   subroutine alloc_1r(p, n, var, sub)
     implicit none
     real(double), pointer, dimension(:) :: p
     integer, intent(in) :: n
     character(len=*), intent(in) :: var, sub
     integer :: error

!    write(777,*) "1r n, var, sub",n, var,"  ", sub
     nullify(p)
     if (associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        stop
     else
        allocate(p(n),stat=error)
        if (error /= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15)') 'Number of elements requested: ',n
        stop
        end if
     end if
   end subroutine alloc_1r

   subroutine alloc_2r(p, n1, n2, var, sub)
     implicit none
     real(double), pointer, dimension(:,:) :: p
     integer, intent(in) :: n1, n2
     character(len=*), intent(in) :: var, sub
     integer :: error

!    write(777,*) "2r n1, n2, var, sub",n1, n2, var,"  ", sub
     nullify(p)
     if (associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        stop
     else
        allocate(p(n1,n2),stat=error)
        if (error /= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15,I15)') 'Number of elements requested: ',n1, n2
        stop
        end if
     end if
   end subroutine alloc_2r

   subroutine alloc_1rl(p, n, var, sub)
     implicit none
     real(double), pointer, dimension(:) :: p
     integer(long), intent(in) :: n
     character(len=*), intent(in) :: var, sub
     integer :: error

!    write(777,*) "1rL n, var, sub",n, var,"  ", sub
     nullify(p)
     if (associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        stop
     else
        allocate(p(n),stat=error)
        if (error /= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15)') 'Number of elements requested: ',n
        stop
        end if
     end if
   end subroutine alloc_1rl

   subroutine alloc_0b(p, first, last, var, sub )
     implicit none
     integer(byte), pointer, dimension(:) :: p
     integer, intent(in) :: first, last
     character(len=*), intent(in) :: var, sub
     integer :: error

!    write(777,*) "0b first, last, var, sub",first, last, var,"  ", sub
     nullify(p)
     if (associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        stop
     else
        allocate(p(first:last),stat=error)
        if (error/= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,2I15)') 'Number of elements requested: ',first, last
        stop
        end if
     end if
   end subroutine alloc_0b

   subroutine alloc_1b(p, n, var, sub)
     implicit none
     integer(byte), pointer, dimension(:) :: p
     integer, intent(in) :: n
     character(len=*), intent(in) :: var, sub
     integer :: error

!    write(777,*) "1b n, var, sub",n, var, "  ", sub
     nullify(p)
     if (associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        stop
     else
        allocate(p(n),stat=error)
        if (error/= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15)') 'Number of elements requested: ',n
        stop
        end if
     end if
   end subroutine alloc_1b

   subroutine alloc_2b(p, n1, n2, var, sub)
     implicit none
     integer(byte), pointer, dimension(:,:) :: p
     integer, intent(in) :: n1, n2
     character(len=*), intent(in) :: var, sub
     integer :: error

!     write(777,*) "2b n1, n2, var, sub",n1, n2, var, "  ", sub
     nullify(p)
     if (associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        stop
     else
        allocate(p(n1,n2),stat=error)
        if (error /= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15,I15)') 'Number of elements requested: ',n1,n2
        stop
        end if
     end if
   end subroutine alloc_2b

   subroutine alloc_3b(p, n1, n2, n3, var, sub)
     implicit none
     integer(byte), pointer, dimension(:,:,:) :: p
     integer, intent(in) :: n1, n2, n3
     character(len=*), intent(in) :: var, sub
     integer :: error

!    write(777,*) "3b n1, n2, n3, var, sub",n1, n2, n3, var, "  ", sub
     nullify(p)
     if (associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        stop
     else
        allocate(p(n1,n2,n3),stat=error)
        if (error /= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,3I15)') 'Number of elements requested: ',n1,n2,n3
        stop
        end if
     end if
   end subroutine alloc_3b

   subroutine alloc_1c(p, n, var, sub)
     implicit none
!CFF     CHARACTER(LEN=256), POINTER, DIMENSION(:) :: p
     character(len=*), pointer, dimension(:) :: p
     integer, intent(in) :: n
     character(len=*), intent(in) :: var, sub
     integer :: error

!    write(777,*) "1c n, var, sub",n, var, "  ", sub
     if (associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        stop
     else
        allocate(p(n),stat=error)
        if (error/= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15)') 'Number of elements requested: ',n
        stop
        end if
     end if
   end subroutine alloc_1c

   subroutine dalloc_1i(p, var, sub)
     implicit none
     integer, pointer, dimension(:) :: p
     character(len=*) :: var, sub
!    write(777,*) "dalloc_1i var, sub",var, "  ", sub
     if (associated(p)) then
        deallocate (p)
     else
        write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     end if
   end subroutine dalloc_1i

   subroutine dalloc_2i(p, var, sub)
     implicit none
     integer, pointer, dimension(:,:) :: p
     character(len=*) :: var, sub
!    write(777,*) "dalloc_2i var, sub",var, "  ", sub
     if (associated(p)) then
        deallocate (p)
     else
        write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     end if
   end subroutine dalloc_2i

   subroutine dalloc_3i(p, var, sub)
     implicit none
     integer, pointer, dimension(:,:,:) :: p
     character(len=*) :: var, sub
!     write(777,*) "dalloc_3i var, sub",var, "  ", sub
     if (associated(p)) then
        deallocate (p)
     else
        write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     end if
   end subroutine dalloc_3i

   subroutine dalloc_1r(p, var, sub)
     implicit none
     real(double), pointer, dimension(:) :: p
     character(len=*) :: var, sub
!    write(777,*) "dalloc_1r var, sub",var, "  ", sub
     if (associated(p)) then
        deallocate (p)
     else
        write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     end if
   end subroutine dalloc_1r

   subroutine dalloc_2r(p, var, sub)
     implicit none
     real(double), pointer, dimension(:,:) :: p
     character(len=*) :: var, sub
!    write(777,*) "dalloc_2r var, sub",var, "  ", sub
     if (associated(p)) then
        deallocate (p)
     else
        write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     end if
   end subroutine dalloc_2r

   subroutine dalloc_1b(p, var, sub)
     implicit none
     integer(byte), pointer, dimension(:) :: p
     character(len=*) :: var, sub
!    write(777,*) "dalloc_1b var, sub",var, "  ", sub
     if (associated(p)) then
        deallocate (p)
     else
        write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     end if
   end subroutine dalloc_1b

   subroutine dalloc_2b(p, var, sub)
     implicit none
     integer(byte), pointer, dimension(:,:) :: p
     character(len=*) :: var, sub
!    write(777,*) "dalloc_2b var, sub",var, "  ", sub
     if (associated(p)) then
        deallocate (p)
     else
        write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     end if
   end subroutine dalloc_2b

   subroutine dalloc_3b(p, var, sub)
     implicit none
     integer(byte), pointer, dimension(:,:,:) :: p
     character(len=*) :: var, sub
!    write(777,*) "dalloc_3b var, sub",var, "  ", sub
     if (associated(p)) then
        deallocate (p)
     else
        write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     end if
   end subroutine dalloc_3b

   subroutine dalloc_1c(p, var, sub)
     implicit none
     character(len=256), pointer, dimension(:) :: p
     character(len=*) :: var, sub
!    write(777,*) "dalloc_1c var, sub",var, "  ", sub
     if (associated(p)) then
        deallocate (p)
     else
        write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     end if
   end subroutine dalloc_1c


   subroutine ralloc_0i(p, first, last, var, sub)
     implicit none
     integer, pointer, dimension(:) :: p, pnew
     integer, intent(in) :: first, last
     character(len=*), intent(in) :: var, sub
     integer :: error, nold
!    write(777,*) "ralloc_0i first, last, var, sub",first, last,var, "  ", sub
     if (.not.  associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        stop
     else
        allocate(pnew(first:last),stat=error)
        if (error/= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to reallocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15)') 'Elements requested for range: ', &
                            first, ' to ', last
        stop
        end if
        nold = min(size(p), last+1) - 1
        pnew(first:nold) = p(first:nold)
        deallocate(p)
        p => pnew
     end if
   end subroutine ralloc_0i


   subroutine ralloc_1i(p, n, var, sub)
     implicit none
     integer, pointer, dimension(:) :: p, pnew
     integer, intent(in) :: n
     character(len=*), intent(in) :: var, sub
     integer :: error, nold
!    write(777,*) "ralloc_1i n, var, sub",n, var, "  ", sub
     if (.not.  associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        stop
     else
        allocate(pnew(n),stat=error)
        if (error/= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to reallocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15)') 'Number of elements requested: ',n
        stop
        end if
        nold = min(size(p),n)
        pnew(1:nold) = p(1:nold)
        deallocate(p)
        p => pnew
     end if
   end subroutine ralloc_1i

   subroutine ralloc_2i(p, n1, n2, var, sub)
     implicit none
     integer, pointer, dimension(:,:) :: p, pnew
     integer, intent(in) :: n1, n2
     character(len=*), intent(in) :: var, sub
     integer :: error, n1_old, n2_old

!    write(777,*) "ralloc_2i n1, n2, var, sub",n1, n2, var, "  ", sub
     if (.not. associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        stop
     else
        allocate(pnew(n1,n2),stat=error)
        if (error /= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15,I15)') 'Number of elements requested: ',n1,n2
        stop
        end if
        n1_old = min(size(p, dim=1),n1)
        n2_old = min(size(p, dim=2),n2)
        pnew(1:n1_old, 1:n2_old) = p(1:n1_old, 1:n2_old)
        deallocate(p)
        p => pnew
     end if
   end subroutine ralloc_2i

   subroutine ralloc_3i(p, n1, n2, n3, var, sub)
     implicit none
     integer, pointer, dimension(:,:,:) :: p, pnew
     integer, intent(in) :: n1, n2, n3
     character(len=*), intent(in) :: var, sub
     integer :: error, n1_old, n2_old, n3_old

!    write(777,*) "ralloc_3i var, sub",var, "  ", sub
     if (.not. associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        stop
     else
        allocate(pnew(n1,n2,n3),stat=error)
        if (error /= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15,I15)') 'Number of elements requested: ',n1,n2,n3
        stop
        end if
        n1_old = min(size(p, dim=1),n1)
        n2_old = min(size(p, dim=2),n2)
        n3_old = min(size(p, dim=3),n3)
        pnew(1:n1_old, 1:n2_old, 1:n3_old) = p(1:n1_old, 1:n2_old, 1:n3_old)
        deallocate(p)
        p => pnew
     end if
   end subroutine ralloc_3i

   subroutine ralloc_1r(p, n, var, sub)
     implicit none
     real(double), pointer, dimension(:) :: p, pnew
     integer, intent(in) :: n
     character(len=*), intent(in) :: var, sub
     integer :: error, nold

!    write(777,*) "ralloc_1r n, var, sub",n, var, "  ", sub
     if (.not. associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        stop
     else
        allocate(pnew(n),stat=error)
        if (error /= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15)') 'Number of elements requested: ',n
        stop
        end if
        nold = min(size(p),n)
        pnew(1:nold) = p(1:nold)
        deallocate(p)
        p => pnew
     end if
   end subroutine ralloc_1r

   subroutine ralloc_2r(p, n1, n2, var, sub)
     implicit none
     real(double), pointer, dimension(:,:) :: p, pnew
     integer, intent(in) :: n1, n2
     character(len=*), intent(in) :: var, sub
     integer :: error, n1_old, n2_old

!    write(777,*) "ralloc_2r n1, n2, var, sub",n1, n2, var, "  ", sub
     if (.not. associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        stop
     else
        allocate(pnew(n1,n2),stat=error)
        if (error /= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15,I15)') 'Number of elements requested: ',n1, n2
        stop
        end if
        n1_old = min(size(p, dim=1),n1)
        n2_old = min(size(p, dim=2),n2)
        pnew(1:n1_old, 1:n2_old) = p(1:n1_old, 1:n2_old)
        deallocate(p)
        p => pnew
     end if
   end subroutine ralloc_2r

   subroutine ralloc_0b(p, first, last, var, sub)
     implicit none
     integer(byte), pointer, dimension(:) :: p, pnew
     integer, intent(in) :: first, last
     character(len=*), intent(in) :: var, sub
     integer :: error, nold
!    write(777,*) "ralloc_0b var, sub",var, "  ", sub
     if (.not.  associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        stop
     else
        allocate(pnew(first:last),stat=error)
        if (error/= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to reallocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15)') 'Elements requested for range: ', &
                            first, ' to ', last
        stop
        end if
        nold = min(size(p), last+1) - 1
        pnew(first:nold) = p(first:nold)
        deallocate(p)
        p => pnew
     end if
   end subroutine ralloc_0b

   subroutine ralloc_1b(p, n, var, sub)
     implicit none
     integer(byte), pointer, dimension(:) :: p, pnew
     integer, intent(in) :: n
     character(len=*), intent(in) :: var, sub
     integer :: error, nold
!    write(777,*) "ralloc_1b var, sub",var, "  ", sub
     if (.not.  associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        stop
     else
        allocate(pnew(n),stat=error)
        if (error/= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to reallocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15)') 'Number of elements requested: ',n
        stop
        end if
        nold = min(size(p),n)
        pnew(1:nold) = p(1:nold)
        deallocate(p)
        p => pnew
     end if
   end subroutine ralloc_1b

   subroutine ralloc_2b(p, n1, n2, var, sub)
     implicit none
     integer(byte), pointer, dimension(:,:) :: p, pnew
     integer, intent(in) :: n1, n2
     character(len=*), intent(in) :: var, sub
     integer :: error, n1_old, n2_old

!    write(777,*) "ralloc_2b var, sub",var, "  ", sub
     if (.not. associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        stop
     else
        allocate(pnew(n1,n2),stat=error)
        if (error /= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15,I15)') 'Number of elements requested: ',n1,n2
        stop
        end if
        n1_old = min(size(p, dim=1),n1)
        n2_old = min(size(p, dim=2),n2)
        pnew(1:n1_old, 1:n2_old) = p(1:n1_old, 1:n2_old)
        deallocate(p)
        p => pnew
     end if
   end subroutine ralloc_2b

   subroutine ralloc_3b(p, n1, n2, n3, var, sub)
     implicit none
     integer(byte), pointer, dimension(:,:,:) :: p, pnew
     integer, intent(in) :: n1, n2, n3
     character(len=*), intent(in) :: var, sub
     integer :: error, n1_old, n2_old, n3_old, i
!     INTEGER :: error, n1_old, n2_old, n3_old

!    write(777,*) "ralloc_3b var, sub",var, "  ", sub
     if (.not. associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        stop
     else
        allocate(pnew(n1,n2,n3),stat=error)
        if (error /= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15,I15)') 'Number of elements requested: ',n1,n2,n3
        stop
        end if
        n1_old = min(size(p, dim=1),n1)
        n2_old = min(size(p, dim=2),n2)
        n3_old = min(size(p, dim=3),n3)
        do i = 1,n2_old
           pnew(1:n1_old, i, 1:n3_old) = p(1:n1_old, i, 1:n3_old)
        end do
!        pnew(1:n1_old, 1:n2_old, 1:n3_old) = p(1:n1_old, 1:n2_old, 1:n3_old)
        deallocate(p)
        p => pnew
     end if
   end subroutine ralloc_3b

   subroutine ralloc_1c(p, n, var, sub)
     implicit none
     character(len=256), pointer, dimension(:) :: p, pnew
     integer, intent(in) :: n
     character(len=*), intent(in) :: var, sub
     integer :: error, nold

!    write(777,*) "ralloc_1c var, sub",var, "  ", sub
     if (.not. associated(p)) then
        write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        stop
     else
        allocate(pnew(n),stat=error)
        if (error /= 0) then
          write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          write(6, '(T10,A,I15)') 'Number of elements requested: ',n
        stop
        end if
        nold = min(size(p),n)
        pnew(1:nold) = p(1:nold)
        deallocate(p)
        p => pnew
     end if
   end subroutine ralloc_1c

end module memory_man
