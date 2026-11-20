module satisfies_breit_group1_triangle_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: satisfies_breit_group1_triangle
contains

    !***********************************************************************
    !                                                                      *
    logical function satisfies_breit_group1_triangle (orbital_a,orbital_b,orbital_c,orbital_d,k)
        ! Module dependencies
        !                                                                      *
        !     Evaluate the triangular relation for the Breit integral of       *
        !     See paper by I. Grant and B McKenzie, JPB 13, (1980) 2671-2681   *
        !     formula 5                                                        *
        !                                                                      *
        !   Written by Per Jonsson                                             *
        !                                                                      *
        !***********************************************************************
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        use orbital_csf_basis, only: orbital_l, orbital_two_j
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_a,orbital_b,orbital_c,orbital_d,k

        ! Local declarations
        integer :: la,lb,lc,ld,ja,jb,jc,jd,l
        !-----------------------------------------------
        !
        la = orbital_l(orbital_a)
        lb = orbital_l(orbital_b)
        lc = orbital_l(orbital_c)
        ld = orbital_l(orbital_d)
        ja = orbital_two_j(orbital_a)
        jb = orbital_two_j(orbital_b)
        jc = orbital_two_j(orbital_c)
        jd = orbital_two_j(orbital_d)

        satisfies_breit_group1_triangle = .false.

        do l = k-1,k+1
            if ((mod(k+la+lb,2)==1).and.(abs(ja-jb)<=2*l).and.     &
                (abs(ja+jb)>=2*l)) then
                if ((mod(k+lc+ld,2)==1).and.(abs(jc-jd)<=2*l).and.  &
                    (abs(jc+jd)>=2*l)) then
                    satisfies_breit_group1_triangle = .true.
                    return
                end if
            end if
        end do

        return
    end function satisfies_breit_group1_triangle
end module satisfies_breit_group1_triangle_interface
