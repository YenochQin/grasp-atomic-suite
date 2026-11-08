! Exhaustive small-integer lookup inputs; compare output with the frozen
! d1260ba implementation, including all three table branches and misses.
program nmtejj_probe
    use nmtejj_I, only: nmtejj
    implicit none
    integer :: j, iq, ij, nk, nd, modes
    do j = 1, 37, 2
        modes = 1
        if (j == 9) modes = 2
        do nk = 1, modes
            do nd = 1, modes
                do iq = 0, (j+1)/2
                    do ij = 0, 99
                        write(*,'(6(I0,1X))') j, iq, ij, nk, nd, nmtejj(iq, ij, j, nk, nd)
                    end do
                end do
            end do
        end do
    end do
    write(*,'(6(I0,1X))') 39, 0, 0, 1, 1, nmtejj(0, 0, 39, 1, 1)
end program nmtejj_probe
