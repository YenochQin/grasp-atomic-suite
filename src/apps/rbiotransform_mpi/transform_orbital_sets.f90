!***********************************************************************
!                                                                      *
subroutine transform_orbital_sets(pi, qi, nli, ninshli, pf, qf, nlf, ninshlf, &
    ngrid, mxl, scr, lscr, ntestg, cishl, cici, cfshl, cfci)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    !                                                                      *
    ! Generate Matrices for rotating radial functions and                  *
    ! for counter rotating CI coefficients                                 *
    !                                                                      *
    ! This subroutine is a modification of the corresponding MCHF          *
    ! routine                                                              *
    !                                                                      *
    ! Per Jonsson,     Lund June 1996                                      *
    !                                                                      *
    ! This subroutine has been modified to support transformations         *
    ! in cases like 3s2p1d 3s2p                                            *
    !                                                                      *
    ! Per Jonsson,     Lund Feb  1997                                      *
    !                                                                      *
    ! =====                                                                *
    ! Input                                                                *
    ! =====                                                                *
    !                                                                      *
    ! PI,QI   : Input radial shells for initial state                      *
    !         : (large and small component)                                *
    ! NLI     : Number of shells per kappa for initial state               *
    ! NINSHLI : Number of inactive shells  per kappa for initial state     *
    !                                                                      *
    ! PF,QF   : Input radial shells for final   state                      *
    !         : (large and small component)                                *
    ! NLF     : Number of shells per kappa for final   state               *
    ! NINSHLF : Number of inactive shells  per kappa for final state       *
    !                                                                      *
    ! NGRID   : Number of gridpoints                                       *
    ! MXL     : Number of kappa quantum numbers                            *
    ! SCR     :                                                            *
    ! LSCR    : Total length of scratch space.                             *
    ! NTESTG  : Global PRINT flag : = 0 => complete silence                *
    !                               = 1 => test and PRINT overlap matrices,*
    !                                      PRINT header and  t matrix      *
    !                            .gt. 1 => Hope you now what you are doing *
    !                                                                      *
    ! ======                                                               *
    ! Output                                                               *
    ! ======                                                               *
    !                                                                      *
    ! CISHL : Rotation matrix for initial state shells                     *
    ! CICI  : Rotation matrix for initial state CI coefficients            *
    ! CFSHL : Rotation matrix for final state shells                       *
    ! CFCI  : Rotation matrix for final state CI coefficients              *
    ! PI,QI : Initial shells in biorthogonal basis                         *
    !         biorthogonal basis                                           *
    ! PF,QF : Final   shells in biorthogonal basis                         *
    !                                                                      *
    !                                                                      *
    ! Inactive Shells : The functions of the inactive shells must be       *
    !                   be supplied to the program, and NLI,NLF            *
    !                   must refer to the total number of                  *
    !                   occupied shells ( inactive+active)                 *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64

    ! Procedure interfaces
    use find_integer_extremum_interface, only: find_integer_extremum
    use sum_integer_vector_interface, only: sum_integer_vector
    !      USE gets_I
    !      USE wrtmat_I
    !      USE copvec_I
    !      USE invmat_I
    !      USE ulla_I
    !      USE trpmat_I
    !      USE matml4_I
    !      USE scalve_I
    !      USE setvec_I
    !      USE pamtmt_I
    use copy_vector_interface, only: copy_vector
    use calculate_subshell_overlap_interface, only: calculate_subshell_overlap
    use invert_dense_submatrix_interface, only: invert_dense_submatrix
    use multiply_submatrices_interface, only: multiply_submatrices
    use build_subshell_transformation_matrix_interface, only: build_subshell_transformation_matrix
    use scale_vector_interface, only: scale_vector
    use fill_vector_interface, only: fill_vector
    use transpose_matrix_interface, only: transpose_matrix
    use combine_upper_lower_factors_interface, only: combine_upper_lower_factors
    use write_matrix_interface, only: write_matrix
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: ngrid,  mxl
    integer, intent(in) :: lscr,  ntestg
    integer, intent(in) :: nli(mxl), nlf(mxl)
    integer :: ninshli(mxl), ninshlf(mxl)
    !     The following have dimensions (NGRID:*)
    !GG      REAL(DOUBLE), DIMENSION(:,:), pointer :: PI, QI, PF, QF
    real(real64) :: pi(ngrid, *), qi(ngrid, *), pf(ngrid, *), qf(ngrid, *)
    !     The following have fixed dimensions (NGRID, NLMAX)
    real(real64)  :: cishl(*), cfshl(*)
    !     The work array
    !     REAL(DOUBLE), DIMENSION(LWORK1) :: work
    !     Th following have fixed dimensions (20*NLMAX*NLMAX)
    real(real64) :: cici(*), cfci(*)
    real(real64) :: scr(*)

    ! Local declarations
    integer :: ntestl, ntest, ili, ilf, nlimx, nlfmx, nlifmx, nlti, &
        nltf, kfree, klpi, klpf, klpif, klstot, klsif, klsifi, &
        klci, klcf, klscr, l, iioff, ifoff, ni, nf, nifmn, iii,&
        jjj, klpmx, klpmn, nmx, nmn, klcmx, klcmn, ndiff, j, &
        klti, i, kltf
    real(real64) :: tii, tiii
    !-----------------------------------------------
    !
    ntestl = 0
    ntest = max(ntestl,ntestg)
    if (ntest >= 1) then
        write (output_unit, *)
        write (output_unit, *) "                   *************************"
        write (output_unit, *) "                   *   Entering BIOTR1     *"
        write (output_unit, *) "                   *************************"
        write (output_unit, *)
    end if
    !
    !. Scratch should at least be of length
    !
    ili = 1
    ilf = 1
    !
    !. Largest number of shells of a given symmetry
    !
    nlimx = find_integer_extremum(nli,mxl,1)
    nlfmx = find_integer_extremum(nlf,mxl,1)
    nlifmx = max(nlimx,nlfmx)
    if (ntest >= 10) write (output_unit, *) " NLIMX,NLFMX NLIFMX ", nlimx, nlfmx, &
        nlifmx
    !
    !. Total numner of shells
    !
    nlti = sum_integer_vector(nli,mxl)
    nltf = sum_integer_vector(nlf,mxl)
    if (ntest >= 10) write (output_unit, *) " NLTI NLTF", nlti, nltf
    !
    ! Scratch space for orbital rotations
    kfree = 1
    !
    klpi = kfree
    kfree = kfree + nlimx*ngrid
    !      WRITE (6, *) ' In biotrn: KLPI    = ', KLPI
    !      WRITE (6, *) '            KFREE   = ', KFREE
    !
    klpf = kfree
    kfree = kfree + nlfmx*ngrid
    !      WRITE (6, *) ' In biotrn: KLPF    = ', KLPF
    !      WRITE (6, *) '            KFREE   = ', KFREE
    !
    klpif = kfree
    kfree = kfree + nlifmx*ngrid
    !      WRITE (6, *) ' In biotrn: KLPIF   = ', KLPIF
    !      WRITE (6, *) '            KFREE   = ', KFREE
    !
    !. Total overlap matrix
    !
    klstot = kfree
    kfree = kfree + nlti*nltf
    !      WRITE (6, *) ' In biotrn: KLSTOT  = ', KLSTOT
    !      WRITE (6, *) '            KFREE   = ', KFREE
    !
    klsif = kfree
    kfree = kfree + nlifmx**2
    !      WRITE (6, *) ' In biotrn: KLSIF   = ', KLSIF
    !      WRITE (6, *) '            KFREE   = ', KFREE
    !
    klsifi = kfree
    kfree = kfree + nlifmx**2
    !      WRITE (6, *) ' In biotrn: KLSIFI  = ', KLSIFI
    !      WRITE (6, *) '            KFREE   = ', KFREE
    !
    klci = kfree
    kfree = kfree + nlifmx**2
    !      WRITE (6, *) ' In biotrn: KLCI    = ', KLCI
    !      WRITE (6, *) '            KFREE   = ', KFREE
    !
    klcf = kfree
    kfree = kfree + nlifmx**2
    !      WRITE (6, *) ' In biotrn: KLCF    = ', KLCF
    !      WRITE (6, *) '            KFREE   = ', KFREE
    !
    klscr = kfree
    kfree = kfree + nlifmx**2 + nlifmx*(nlifmx + 1)
    !      WRITE (6, *) ' In biotrn: KLSCR   = ', KLSCR
    !      WRITE (6, *) '            KFREE   = ', KFREE
    !      WRITE (6, *) '         =>  FREE =     ', KFREE
    !
    !. Check length of scratch
    !
    if (lscr <= kfree - 1) then
        write (output_unit,*)" BIOTR1 in trouble ! "
        write (output_unit,*)" Increase dimension of scratch before call to BIOTR1"
        write (output_unit,*)" Current and required length (LSCR,KFREE-1)", lscr, &
            kfree - 1
        stop "Increase LWORK before call to BIOTR1"
    end if
    !
    !. Obtain overlap matrix
    !
    call calculate_subshell_overlap (scr(klstot), nlti,nltf)

    do l = 1, mxl
        if (ntest >= 5) then
            write (output_unit, *) "   L = ", l
            write (output_unit, *) "       Orbital rotation..."
        end if
        !
        !. Offsets for given L in shell matrices
        !
        if (l == 1) then
            iioff = 1
            ifoff = 1
        else
            iioff = iioff + nli(l-1)**2
            ifoff = ifoff + nlf(l-1)**2
        end if
        if (ntest >= 1) then
            write (output_unit, *)
            write (output_unit, *) &
                " BIOTRN : Information on transformations of shells with L =",l
            write (output_unit, *)
        end if
        !
        ! =========================================================
        ! 1 : Obtain Biorthogonal forms of initial and final shells
        ! =========================================================
        !
        !. The overlap matrix can be written
        !
        !     * * * * * * *
        !     *       *   *
        !     *   S   * X *
        !     *       *   *
        !     * * * * * * *
        !
        ! where S is the quadratic subblock.
        ! The basis functions corresponding to the quadratic
        ! subblock are made biorthogonal with an UL decomposition.
        ! The remaining basis
        ! functions are made biorthogonal by choosing  the
        ! transformation
        !            * * *
        !            *   *
        !            * Y *
        !            *   *
        !            * * *
        !            * 1 *
        !            * * *
        !
        !
        ! With Y = -S-1*X
        !
        ni = nli(l)
        nf = nlf(l)
        !
        !.1.1 : Obtain shells of given L in proper order for
        !       biorthogonal treatment
        !
        if (l == 1) then
            ili = 1
            ilf = 1
        else
            ili = ili + nli(l-1)
            ilf = ilf + nlf(l-1)
        end if

        !
        ! 1.2 obtain biorthogonal of the first min(ni,nf) shells
        !
        nifmn = min(ni,nf)
        !
        !. Overlap matrix SIF = Integral (PI(I)*PF(J))
        !
        !ww Per change to support cases like 3s2p1d 3s2p. The belonging endif
        !ww Is just at the end

        if (nifmn <= 0) cycle

        do iii = 1, nifmn
            do jjj = 1, nifmn
                scr(klsif+(jjj-1)*nifmn+iii-1) = scr(klstot-1+(jjj+ilf-1-1)*nlti&
                    +iii+ili-1)
            end do
        end do
        if (ntest >= 15) then
            write (output_unit, *) " Overlap matrix "
            call write_matrix (scr(klsif), nifmn, nifmn, nifmn, nifmn)
        end if
        !
        ! Obtain upper triangular CI and CF so CI(T) S CF = 1
        ! or CF CI(T) = S-1, which corresponds to an UL decomposition
        !
        !. Invert S
        !
        call copy_vector (scr(klsif), scr(klsifi), nifmn**2)
        call invert_dense_submatrix (scr(klsifi), scr(klci), nifmn, nifmn)

        !. UL decompose
        call copy_vector (scr(klsifi), scr(klsif), nifmn**2)
        call combine_upper_lower_factors (scr(klsif), scr(klcf), scr(klci), nifmn, scr(klscr))
        call transpose_matrix (scr(klci), nifmn, nifmn, scr(klscr))
        call copy_vector (scr(klscr), scr(klci), nifmn**2)
        !
        !. The transformation matrix between the first NIFMX
        !. shells is now known, biorthogonalize remaining orbitals
        !
        if (ni/=nf .and. ni/=0 .and. nf/=0) then
            if (ni > nf) then
                klpmx = klpi
                klpmn = klpf
                nmx = ni
                nmn = nf
                klcmx = klci
                klcmn = klcf
            else
                klpmx = klpf
                klpmn = klpi
                nmx = nf
                nmn = ni
                klcmx = klcf
                klcmn = klci
            end if
            ndiff = nmx - nmn
            !
            ! Y = -S-1 * X
            !. overlap X between remaining orbitals and the other set
            !
            if (ni > nf) then
                !
                ! I columns F rows
                !
                do iii = nmn + 1, nmx
                    do jjj = 1, nf
                        scr(klsif+(iii-nmn-1)*nf+jjj-1) = scr(klstot-1+(jjj+ilf-1-&
                            1)*nlti+iii+ili-1)
                    end do
                end do
            else if (nf > ni) then
                ! F columns I rows
                do jjj = nmn + 1, nmx
                    do iii = 1, ni
                        scr(klsif+(jjj-nmn-1)*ni+iii-1) = scr(klstot-1+(jjj+ilf-1-&
                            1)*nlti+iii+ili-1)
                    end do
                end do
            end if
            !
            if (ni > nf) then
                call transpose_matrix (scr(klsifi), nmn, nmn, scr(klscr))
                call copy_vector (scr(klscr), scr(klsifi), nmn**2)
            end if
            call multiply_submatrices (scr(klscr), scr(klsifi), scr(klsif), nmn, ndiff, nmn, &
                nmn, nmn, ndiff, 0)
            call scale_vector (scr(klscr), -1.0_real64, nmn*ndiff)
            call copy_vector (scr(klscr), scr(klsif), nmn*ndiff)
            !
            ! Construct complete CMX
            !
            call fill_vector (scr(klscr), 0.0_real64, nmx**2)
            do j = 1, nmx
                if (j <= nifmn) then
                    call copy_vector (scr(klcmx+(j-1)*nifmn), scr(klscr+(j-1)*nmx), &
                        nmn)
                else
                    call copy_vector (scr(klsif+(j-nmn-1)*nmn), scr(klscr+(j-1)*nmx), &
                        nmn)
                    scr(klscr-1+(j-1)*nmx+j) = 1.0_real64
                end if
            end do
            !
            call copy_vector (scr(klscr), scr(klcmx), nmx**2)
        end if
        !ww Pertest
        !         ENDIF
        !
        !. The two upper triangular matrices CI and CF are now known
        !. Transfer to permanent arrays
        !
        call copy_vector (scr(klci), cishl(iioff), ni**2)
        call copy_vector (scr(klcf), cfshl(ifoff), nf**2)
        !
        !. Rotate the large component of the shells
        !
        call copy_vector (pi(1,ili), scr(klpi), ni*ngrid)
        call copy_vector (pf(1,ilf), scr(klpf), nf*ngrid)
        !
        !         WRITE (*, *) 'Transformation matrices initial'
        !         CALL WRTMAT (SCR(KLCI), NI, NI, NI, NI)
        call multiply_submatrices (scr(klpif), scr(klpi), scr(klci), ngrid, ni, ngrid, ni, &
            ni, ni, 0)
        call copy_vector (scr(klpif), pi(1,ili), ni*ngrid)
        !         WRITE (*, *) 'Transformation matrices final'
        !         CALL WRTMAT (SCR(KLCF), NF, NF, NF, NF)
        call multiply_submatrices (scr(klpif), scr(klpf), scr(klcf), ngrid, nf, ngrid, nf, &
            nf, nf, 0)
        call copy_vector (scr(klpif), pf(1,ilf), nf*ngrid)
        !
        !. Rotate the small component of the shells
        !
        call copy_vector (qi(1,ili), scr(klpi), ni*ngrid)
        call copy_vector (qf(1,ilf), scr(klpf), nf*ngrid)
        !
        call multiply_submatrices (scr(klpif), scr(klpi), scr(klci), ngrid, ni, ngrid, ni, &
            ni, ni, 0)
        call copy_vector (scr(klpif), qi(1,ili), ni*ngrid)
        call multiply_submatrices (scr(klpif), scr(klpf), scr(klcf), ngrid, nf, ngrid, nf, &
            nf, nf, 0)
        call copy_vector (scr(klpif), qf(1,ilf), nf*ngrid)
        !
        if (ntest >= 1) then
            write (output_unit, *) " Test of overlap of biorthonormal functions"
            ! F columns I rows
            do jjj = 1, nf
                do iii = 1, ni
                    scr(klsif+(jjj-1)*ni+iii-1) = scr(klstot-1+(jjj+ilf-1-1)*nlti&
                        +iii+ili-1)
                end do
            end do
            call multiply_submatrices (scr(klscr), scr(klci), scr(klsif), ni, nf, ni, ni, ni&
                , nf, 1)
            call multiply_submatrices (scr(klsif), scr(klscr), scr(klcf), ni, nf, ni, nf, nf&
                , nf, 0)
            write (output_unit, *) &
                " new overlap matrix ( should be 1 on diag, 0 elsewhere )"
            call write_matrix (scr(klsif), ni, nf, ni, nf)
        end if

        if (ntest >= 1) then
            write (output_unit, *)
            write (output_unit, *) " Orbital Rotation matrix for I state"
            call write_matrix (cishl(iioff), ni, ni, ni, ni)
            write (output_unit, *) " Orbital Rotation matrix for F state"
            call write_matrix (cfshl(ifoff), nf, nf, nf, nf)
            write (output_unit, *)
        end if
        !
        !. Matrix for counterrotation of CI coefficients, initial state
        !
        klti = klsif
        call build_subshell_transformation_matrix (scr(klci), scr(klti), scr(klscr), ni)
        do i = 1, ni
            tii = scr(klti-1+(i-1)*ni+i)
            tiii = 1.0_real64/tii
            call scale_vector (scr(klti+(i-1)*ni), tiii, i - 1)
        end do
        call copy_vector (scr(klti), cici(iioff), ni*ni)
        !
        !. Matrix for counterrotation of CI coefficients, Final state
        !
        kltf = klsif
        call build_subshell_transformation_matrix (scr(klcf), scr(kltf), scr(klscr), nf)
        do i = 1, nf
            tii = scr(kltf-1+(i-1)*nf+i)
            tiii = 1.0_real64/tii
            call scale_vector (scr(kltf+(i-1)*nf), tiii, i - 1)
        end do
        call copy_vector (scr(kltf), cfci(ifoff), nf*nf)
        if (ntest < 1) cycle
        write (output_unit, *)
        write (output_unit, *) " CI-Rotation matrix for I state"
        call write_matrix (cici(iioff), ni, ni, ni, ni)
        write (output_unit, *) " CI-Rotation matrix for F state"
        call write_matrix (cfci(ifoff), nf, nf, nf, nf)
        write (output_unit, *)

    end do
    !
    return
end subroutine transform_orbital_sets
