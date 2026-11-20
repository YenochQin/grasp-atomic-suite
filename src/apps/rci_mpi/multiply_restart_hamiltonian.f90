!***********************************************************************
!                                                                      *
subroutine multiply_restart_hamiltonian(n, m, b, c)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    !                                                                      *
    !   Matrix-matrix product: C = AB.  A  sparse  representation of the   *
    !   lower triangle of the  (NxN)  matrix  A  is read from the disk     *
    !                                                                      *
    !   This is an adaptation of  Andreas Stathopulos'  routine  SPSBMV,   *
    !   and is specific to GRASP2 derivatives.                             *
    !                                                                      *
    !   Call(s) to: [AUXBLAS]: DINIT/SINIT;                                *
    !               [SPBLAS]: DAXPYI/SAXPYI, DDOTI/SDOTI.                  *
    !                                                                      *
    !   F A Parpia and A Stathopoulos         Last revision: 13 Oct 1992   *
    !   Block Version by Xinghong He          Last revision: 18 Jun 1998   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
        csf_mixing_coefficients, evecff, evecii
    use ci_restart_units, only: restart_matrix_unit
    use rci_matrix_file_position, only: nposition
    use mpi_state,           only:  mpi_rank, mpi_process_count
    use real_vector_initialization, only: dinit
    use file_positioning, only: posfile
    use sparse_dot_and_scatter, only: dmerge
    use mpi_utilities, only: gdsummpi
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: n
    integer, intent(in) :: m
    real(real64), intent(in) :: b(n, m)
    real(real64), intent(inout) :: c(n, m)

    ! Local declarations
    integer :: irow(n)
    !cjb  INTEGER :: MYID, NPROCS, NCF, NCFDUM, ICCUTDUM, MYIDDUM, NPROCSDUM, ICOL&
    integer ::               csf_count, file_csf_count, file_first_order_cutoff, file_mpi_rank, file_mpi_process_count, icol&
        , nelc, ir, iv
    real(real64) :: emt(n)
    real(real64) :: elsto, diag, dl
    !-----------------------------------------------
    !
    !     !...nposition+1 is the current position of the .res file
    !     !...It is set in matrix and used in maneig, spodmv
    !
    !-----------------------------------------------------------------------
    write (output_unit, *) "Calling spodmv..."
    csf_count = n

    !   Initialise the result matrix; note that this is specific to the
    !   data structure of DVDSON

    call dinit (n*m, 0.0_real64, c, 1)

    !...moved from maneig before "CALL GDVD (SPODMV..."
    call posfile (0, restart_matrix_unit, nposition)

    read (restart_matrix_unit) file_csf_count, file_first_order_cutoff, file_mpi_rank, file_mpi_process_count
    if (csf_count/=file_csf_count .or. mpi_rank/=file_mpi_rank .or. file_mpi_process_count/=mpi_process_count) stop &
        "spodmv: ncf read wrong"


    do icol = mpi_rank + 1, n, mpi_process_count
        read (restart_matrix_unit) nelc, elsto, (emt(ir),ir=1,nelc), (irow(ir),ir=1,nelc)
        do iv = 1, m
            diag = c(icol,iv) + (emt(nelc)-hamiltonian_energy_shift)*b(icol,iv)
            call dmerge (nelc - 1, b(1:n,iv), c(1:n,iv), irow(1:n), &
                emt(1:n), b(icol,iv), dl)
            c(icol,iv) = diag + dl
        end do
    end do

    call gdsummpi (c, n*m)

    return
end subroutine multiply_restart_hamiltonian
