program csfg_reallocation
    use vast_kind_param, only: double, long
    use csfg_memory_man, only: alloc, ralloc, dalloc
    use symmatrix_restart_C, only: nelmntres
    implicit none
    integer, pointer :: integers(:) => null()
    real(double), pointer :: reals(:) => null()
    integer :: i

    if (storage_size(nelmntres) /= 64) error stop 'restart count width changed'
    call alloc(integers, 5_long, 'integers', 'test')
    integers = [(i, i=1,5)]
    call ralloc(integers, 3_long, 'integers', 'test')
    if (size(integers,kind=long) /= 3_long) error stop 'integer shrink extent'
    if (any(integers /= [1,2,3])) error stop 'integer shrink lost values'
    call ralloc(integers, 7_long, 'integers', 'test')
    if (size(integers,kind=long) /= 7_long) error stop 'integer grow extent'
    if (any(integers(:3) /= [1,2,3])) error stop 'integer grow lost values'
    call dalloc(integers, 'integers', 'test')

    call alloc(reals, 5_long, 'reals', 'test')
    reals = [(real(i,double)/4, i=1,5)]
    call ralloc(reals, 3_long, 'reals', 'test')
    if (size(reals,kind=long) /= 3_long) error stop 'real shrink extent'
    if (any(reals /= [0.25_double,0.5_double,0.75_double])) error stop 'real shrink lost values'
    call ralloc(reals, 7_long, 'reals', 'test')
    if (size(reals,kind=long) /= 7_long) error stop 'real grow extent'
    if (any(reals(:3) /= [0.25_double,0.5_double,0.75_double])) error stop 'real grow lost values'
    call dalloc(reals, 'reals', 'test')
end program csfg_reallocation
