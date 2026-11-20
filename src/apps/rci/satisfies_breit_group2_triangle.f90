!***********************************************************************
!                                                                      *
logical function satisfies_breit_group2_triangle (orbital_a,orbital_b,orbital_c,orbital_d,l)
    ! Module dependencies
    !                                                                      *
    !     Evaluate the triangle relations for Breit integrals of type 2    *
    !     See Grant and McKenzier JPB 13 (1980), 2675-2681 formula 5       *
    !                                                                      *
    !   Written by Per Jonsson                                             *
    !                                                                      *
    !***********************************************************************
    !...Translated by Charlotte Froese Fischer
    !                       Gediminas Gaigalas  10/05/17
    use orbital_csf_basis, only: orbital_l, orbital_two_j
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: orbital_a,orbital_b,orbital_c,orbital_d,l

    ! Local declarations
    integer :: la,lb,lc,ld,ja,jb,jc,jd
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

    satisfies_breit_group2_triangle = .false.

    if ((mod(l+la+lb,2)==0).and.(abs(ja-jb)<=2*l).and.           &
        (abs(ja+jb)>=2*l)) then
        if ((mod(l+lc+ld,2)==0).and.(abs(jc-jd)<=2*l).and.        &
            (abs(jc+jd)>=2*l)) then
            satisfies_breit_group2_triangle = .true.
            return
        end if
    end if

    return
end function satisfies_breit_group2_triangle
