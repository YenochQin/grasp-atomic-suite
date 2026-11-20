!***********************************************************************
!                                                                      *
subroutine apply_occupation_factor(input_mixing_coefficients, csf_count, level_count, i, l, const, &
    & transformed_mixing_coefficients, ntestg)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    !                                                                      *
    !   Calculates the action of the operator                              *
    !   Const ** E(li,li) on a set of vectors                              *
    !                                                                      *
    !   Adapted for GRASP, daughter of TIINI, born February 1996           *
    !                                                                      *
    !   =====                                                              *
    !   Input                                                              *
    !   =====                                                              *
    !   CIIN      : Input CI vectors                                       *
    !   NCSCF     : Length of CI expansion                                 *
    !   NCIV      : Number of CI vectors                                   *
    !   I         : Shell number                                           *
    !   L         : symmetry                                               *
    !   NSHLP(L,K): Gives the shell number as defined in getcsl for        *
    !               the K:th shell with symmetry L                         *
    !   CONST     : The constant                                           *
    !                                                                      *
    !   ======                                                             *
    !   Output                                                             *
    !   ======                                                             *
    !                                                                      *
    !   CIOUT : List of output CI vectors                                  *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: orbital_index_radix
    use mcp_coefficient_storage, only: cnn, intgrl, intptr, jann, jbnn, ncoeff, nintg
    use biorthogonal_subshell_mapping, only: max_orbitals, nlmax, nshlp, nshlpp

    ! Procedure interfaces
    use fill_vector_interface, only: fill_vector
    use write_matrix_interface, only: write_matrix
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: csf_count
    integer, intent(in) :: level_count
    integer, intent(in) :: i
    integer, intent(in) :: l
    integer, intent(in) :: ntestg
    real(real64), intent(in) :: const
    real(real64), intent(in) :: input_mixing_coefficients(csf_count,level_count)
    real(real64)  :: transformed_mixing_coefficients(csf_count,level_count)

    ! Local constants
    integer, parameter :: orbital_label_radix = orbital_index_radix

    ! Local declarations
    integer :: ntestl, ntest, nfound, k, ia, ib, ifirst, ielmnt, ival, ileft&
        , ivec
    real(real64) :: constn
    !-----------------------------------------------
    !
    !
    ntestl = 0
    ntest = max(ntestl,ntestg)


    if (ntest >= 10) write (output_unit, *) " Entering TIINI"

    call fill_vector (transformed_mixing_coefficients, 0.0_real64, csf_count*level_count)
    !
    !.  Obtain address of first coupling coefficient for h(il,il) :  IFIRST
    !.  Obtain number of        coupling coefficient for h(il,il) :  NFOUND
    !.
    !     NFOUND : Number of coefficients obtained
    !     IVAL   : actual RACAH coefficient <CSF(L)!E(il,il)!CSF(R)>
    !     ILEFT  = CSF(L) ?
    !
    nfound = 0
    do k = 1, nintg
        ia = intgrl(k)/orbital_label_radix
        ib = mod(intgrl(k),orbital_label_radix)
        if (nshlp(l,i)/=ia .or. ia/=ib) cycle
        if (k == 1) then
            ifirst = 1
        else
            ifirst = intptr(k - 1) + 1
        end if
        nfound = intptr(k) - ifirst + 1
        exit
    end do
    do ielmnt = 1, nfound
        ! Bug 2011-08-18 Per Jonsson        IVAL = CNN(IFIRST-1+IELMNT)
        ival = idnint(cnn(ifirst-1+ielmnt))
        constn = const**ival
        ileft = jann(ifirst - 1 + ielmnt)
        transformed_mixing_coefficients(ileft,:level_count) = constn*input_mixing_coefficients(ileft,:level_count)
    end do
    !
    !.  The previous provided us with all
    !   terms with nonvanishing occupation.
    !   For terms with vanishing occupation of il,
    !   just copy coefficients, since (x) ** 0 = 1
    !
    where (transformed_mixing_coefficients(:csf_count,:level_count) == 0.0_real64)
        transformed_mixing_coefficients(:csf_count,:level_count) = input_mixing_coefficients(:csf_count,:level_count)
    end where
    !
    if (ntest >= 100) then
        write (output_unit, *)
        write (output_unit, *) " Input and output vectors from TIINI I,L", i, l
        call write_matrix (input_mixing_coefficients, csf_count, level_count, csf_count, level_count)
        write (output_unit, *)
        call write_matrix (transformed_mixing_coefficients, csf_count, level_count, csf_count, level_count)
        write (output_unit, *)
    end if
    !
    if (ntest >= 10) write (output_unit, *) " Leaving  TIINI"

    return
end subroutine apply_occupation_factor
