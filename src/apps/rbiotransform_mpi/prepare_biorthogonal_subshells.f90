!***********************************************************************
!                                                                      *
subroutine prepare_biorthogonal_subshells(ntestg, ncore1, ncore2)
    ! Module dependencies
    !                                                                      *
    !   This subroutine determines the number of kappa quantum numbers     *
    !   KAMAX together with the number of orbitals of each kappa.          *
    !                                                                      *
    !   This subroutine also checks if the orbital ordering is normal      *
    !   or reversed. If reversed then JA and JB have to be permuted        *
    !   In the subroutine ti1tv                                            *
    !                                                                      *
    !   Written by Per Jonsson                                             *
    !                                                                      *
    !***********************************************************************
    use grasp_parameters,   only: max_orbitals
    use biorthogonal_orbital_sets,         only: nwii,nakii,npii,nwff,nakff,npff
    use biorthogonal_subshell_layout,         only: nlmax,ninii,ninff,nshlii,nshlff,      &
        ikappa,nakinvii,nakinvff,nlmax,       &
        nshlpii,nshlpff,nshlppii,nshlppff,kamax
    use biorthogonal_orbital_order, only: nordff, nordii
    implicit none (type, external)

    ! Arguments
    integer , intent(in) :: ntestg
    integer , intent(in) :: ncore1
    integer , intent(in) :: ncore2

    ! Local declarations
    integer :: j1
    integer :: isort(2*max_orbitals)
    integer :: ntestl, ntest, i, j, k, ital, nref
    logical :: klar
    !-----------------------------------------------
    !
    !
    !
    ntestl = 0
    ntest = max0(ntestl,ntestg)

    isort(:2*max_orbitals) = 100

    ikappa(:nlmax) = 0
    nshlii(:nlmax) = 0
    nshlff(:nlmax) = 0
    ninii(:nlmax) = 0
    ninff(:nlmax) = 0

    nakinvii(:max_orbitals) = 0
    nakinvff(:max_orbitals) = 0

    nshlpii(:nlmax,:nlmax) = 0
    nshlpff(:nlmax,:nlmax) = 0

    nshlppii(:nlmax,:max_orbitals) = 0
    nshlppff(:nlmax,:max_orbitals) = 0
    !
    !   Sort the kappa quantum numbers
    !
    do k = 1, nwii + nwff
        if (k <= nwii) then
            ital = nakii(k)
        else
            ital = nakff(k-nwii)
        end if
        i = k - 1
        klar = .false.
        do while (i>0 .and. .not.klar)
            if (ital <= isort(i)) then
                isort(i+1) = isort(i)
                i = i - 1
            else
                klar = .true.
            end if
        end do
        isort(i+1) = ital
    end do
    !
    !   Determine the unique set of kappa IKAPPA
    !
    kamax = 1
    ikappa(1) = isort(1)
    do k = 1, 2*max_orbitals - 1
        if (isort(k)==isort(k+1) .or. isort(k+1)>=100) cycle
        kamax = kamax + 1
        ikappa(kamax) = isort(k+1)
    end do
    !
    !  Make a connection between each kappa and a number in the
    !  range [1,KAMAX] as to know on which file to dump the data
    !  Determine the number of shells NSHLII for each I in the
    !  range [1,KAMAX]
    !
    if (ntest >= 10) then
        write (*, *) "******************"
        write (*, *) " Entering kapdata"
        write (*, *) "******************"
        write (*, *)
        write (*, *) "There are", nwii, "orbitals in the initial state"
        write (*, *) "with the following n and kappa quantum numbers"
    end if

    do j = 1, nwii
        if (ntest >= 10) write (*, *) "orbital number", j, "n and kappa", npii&
            (j), nakii(j)
        if (j <= ncore1) then
            do i = 1, kamax
                if (ikappa(i) /= nakii(j)) cycle
                nakinvii(j) = i
                ninii(i) = ninii(i) + 1
                nshlii(i) = nshlii(i) + 1
                nshlpii(i,nshlii(i)) = j
                nshlppii(i,j) = nshlii(i)
            end do
        else
            do i = 1, kamax
                if (ikappa(i) /= nakii(j)) cycle
                nakinvii(j) = i
                nshlii(i) = nshlii(i) + 1
                nshlpii(i,nshlii(i)) = j
                nshlppii(i,j) = nshlii(i)
            end do
        end if
    end do

    if (ntest >= 10) then
        write (*, *) "There are", nwff, "orbitals in the final state"
        write (*, *) "with the following n and kappa quantum numbers"
    end if
    do j = 1, nwff
        if (ntest >= 10) write (*, *) "orbital number", j, "n and kappa=", &
            npff(j), nakff(j)
        if (j <= ncore2) then
            do i = 1, kamax
                if (ikappa(i) /= nakff(j)) cycle
                nakinvff(j) = i
                ninff(i) = ninff(i) + 1
                nshlff(i) = nshlff(i) + 1
                nshlpff(i,nshlff(i)) = j
                nshlppff(i,j) = nshlff(i)
            end do
        else
            do i = 1, kamax
                if (ikappa(i) /= nakff(j)) cycle
                nakinvff(j) = i
                nshlff(i) = nshlff(i) + 1
                nshlpff(i,nshlff(i)) = j
                nshlppff(i,j) = nshlff(i)
            end do
        end if
    end do

    if (ntest >= 10) then
        write (*, *) "Total number of different kappa", kamax
        do i = 1, kamax
            write (*, *) "L=", i, "corresponds to kappa=", ikappa(i)
            write (*, *) "nr of init.  orb. with this kappa=", nshlii(i)
            write (*, *) "nr of final. orb. with this kappa=", nshlff(i)
        end do
        do i = 1, kamax
            write (*, *) "Position in initial state list"
            do j = 1, nshlii(i)
                write (*, *) "L=", i, "orb. nr", j, ",position", nshlpii(i,j)
            end do
            write (*, *) "Position in final   state list"
            do j = 1, nshlff(i)
                write (*, *) "L=", i, "orb. nr", j, ",position", nshlpff(i,j)
            end do
            write (*, *) "Relative positions for initial state orbitals"
            do j = 1, nwii
                if (nshlppii(i,j) == 0) cycle
                write (*, *) "Orbital", j, "is nr", nshlppii(i,j), "with kappa"&
                    , ikappa(i)
            end do
            write (*, *) "Relative positions for final state orbitals"
            do j = 1, nwff
                if (nshlppff(i,j) == 0) cycle
                write (*, *) "Orbital", j, "is nr", nshlppff(i,j), "with kappa"&
                    , ikappa(i)
            end do
        end do

    end if
    !
    !  Check if the orbital ordering is normal or reversed.
    !
    !  For normal ordering
    !
    !    NPII(NSHLPII(I,1)) < NPII(NSHLPII(I,2)) < ... < NPII(NSHLPII(I,NSHLII(I))
    !
    !  For reversed ordering
    !
    !    NPII(NSHLPII(I,1)) > NPII(NSHLPII(I,2)) > ... > NPII(NSHLPII(I,NSHLII(I))
    !
    nordii = 0
    do i = 1, kamax
        nref = 0
        do j = 1 + ninii(i), nshlii(i)
            if (npii(nshlpii(i,j)) < nref) nordii = 1
            nref = npii(nshlpii(i,j))
        end do
    end do

    nordff = 0
    do i = 1, kamax
        nref = 0
        do j = 1 + ninff(i), nshlff(i)
            if (npff(nshlpff(i,j)) < nref) nordff = 1
            nref = npff(nshlpff(i,j))
        end do
    end do

    if (nordii /= nordff) then
        write (*, *) " Orbital order of the initial and final states"
        write (*, *) " should be the same. STOP"
        stop
    end if
    !
    !   If not normal order check if reversed order
    !
    if (nordii == 1) then
        do i = 1, kamax
            nref = 0
            do j = nshlii(i), 1 + ninii(i), -1
                if (npii(nshlpii(i,j)) < nref) nordii = 2
                nref = npii(nshlpii(i,j))
            end do
        end do

        do i = 1, kamax
            nref = 0
            do j = nshlff(i), 1 + ninff(i), -1
                if (npff(nshlpff(i,j)) < nref) nordff = 2
                nref = npff(nshlpff(i,j))
            end do
        end do
    end if

    if (nordii==2 .or. nordff==2) then
        write (*, *) " The orbital order is neither normal or reversed"
        write (*, *) " STOP"
        stop
    end if

    !w      write(*,*) 'Give nordii'
    !w      read(*,*) nordii
    !w      nordff = nordii

    if (nordii == 0) then
        write (*, *) " Normal orbital ordering"
    else
        write (*, *) " Reverse orbital ordering"
    end if

    !      WRITE (*, *) '*****************'
    !      WRITE (*, *) ' Leaving kapdata'
    !      WRITE (*, *) '*****************'
    !      WRITE (*, *)

    return
end subroutine prepare_biorthogonal_subshells
