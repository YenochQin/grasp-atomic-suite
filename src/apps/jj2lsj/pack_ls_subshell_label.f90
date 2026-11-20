module pack_ls_subshell_label_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: pack_ls_subshell_label
contains

    !
    !***********************************************************************
    !                                                                      *
    subroutine pack_ls_subshell_label(subshell_count, subshell_labels, subshell_occupations, coupling_labels, packed_label)
        ! Procedure interfaces
        !                                                                      *
        !   Subroutine written by Bin LIU                                      *
        !                                                                      *
        !   Rules for encoding                                                 *
        !   1. All blanks deleted                                              *
        !   2. If Qi=1, omit Qi                                                *
        !   3. If Qi=1 or Qi>=4l+1, omit ALFAi                                 *
        !   4. If i=1 or (Qi=4l+2 and i<>m), insert '.'; else _BETAi.          *
        !                                                                      *
        !   Modified by G. Gaigalas,                                May 2011   *
        !                                                                      *
        !***********************************************************************
        !...Switches:
        use orbital_letter_to_l_interface, only: orbital_letter_to_l
        implicit none (type, external)

        ! Arguments
        integer, intent(in)            :: subshell_count
        !GG      CHARACTER(LEN=64), INTENT(OUT) :: STR
        character(len=164), intent(out) :: packed_label
        integer, intent(in)            :: subshell_occupations(*)
        character(len=3), intent(in)   :: subshell_labels(*), coupling_labels(*)

        ! Local declarations
        integer          :: full1, full, const, j, i, k, n, n1, nn, icore
        integer          :: l, previous_q
        character(len=1) :: ch1
        character(len=3) :: ch3
        !-----------------------------------------------
        !  FULL   :  4l+2
        !  CONST  :   constant for converting lowercase to uppercase
        !  CH*    :  temporary variables
        !
        const = ichar("a") - ichar("A")
        packed_label = " "
        j = 1
        icore = 0
        full = 0; n = 0; full1 = 0; n1 = 0
        !
        !   -----  begin to encode  -----
        !
        do i = 1, subshell_count
            !        Fortran does not guarantee short-circuit evaluation of .AND.;
            !        never reference Q(0), even when a condition starts with I /= 1.
            previous_q = 0
            if (i > 1) previous_q = subshell_occupations(i-1)
            if (subshell_occupations(i) == 0) cycle
            if (subshell_labels(i)(1:1) == " ") then
                packed_label(j:j+1) = subshell_labels(i)(2:3)
                k = 2
            else
                packed_label(j:j+2) = subshell_labels(i)
                k = 3
                if (subshell_labels(i)(3:3) == " ") k = 2
            end if
            ch1 = packed_label(j+1:j+1)
            if (ch1>="A" .and. ch1<="Z") packed_label(j+1:j+1) = char(ichar(ch1) + const)
            if(i /= 1) full1 = full
            if(i /= 1) n1 = n
            l = orbital_letter_to_l(ch1)
            full = 4*l + 2
            !
            !  -----  convert Qi into character  -----
            !
            j = j + k
            n = subshell_occupations(i)
            if (n > 14) stop "Too many electrons"
            if(i /= 1 .and. full1 /= n1) icore = icore + full1 - n1
            !
            !  -----  If Qi<>1, add Qi
            !           If Qi<4l+1, add TERMi for the shell -----
            !
            if (n /= 1) then
                packed_label(j:j) = "("
                j = j + 1
                nn = n
                if (n > 9) then
                    packed_label(j:j) = "1"
                    j = j + 1
                    nn = n - 10
                end if
                packed_label(j:j+1) = char(ichar("0") + nn)//")"
                j = j + 2
                if(i/=1 .and. icore<2 .and. previous_q<=1 .and. n==2) then
                    continue
                else if(i/=1 .and. icore==0 .and. previous_q==full1 .and. n<3  &
                    .and. subshell_count==i) then
                    continue
                else if(i/=1 .and. icore==0 .and. previous_q==full1 .and. l<2  &
                    .and. subshell_count==i) then
                    continue
                else if(i/=1 .and. icore==0 .and. previous_q==full1 .and. l==2 &
                    .and. n > 7 .and. subshell_count==i) then
                    continue
                else if(i/=1 .and. icore==0 .and. previous_q==full1 .and. l==3 &
                    .and. n > 11 .and. subshell_count==i) then
                    continue
                else if (n<(full-1) .and. subshell_count/=1) then
                    ch3 = coupling_labels(i)
                    ch1 = ch3(2:2)
                    if (ch1>="a" .and. ch1<="z") ch3(2:2) = char(ichar(ch1) - const)
                    if(len_trim(ch3) == 2) then
                        packed_label(j:j+1) = trim(ch3)
                        j = j + 2
                    else if(l == 2 .and. n < 3) then
                        packed_label(j:j+1) = ch3(1:2)
                        j = j + 2
                    else if(l == 2 .and. n > 7) then
                        packed_label(j:j+1) = ch3(1:2)
                        j = j + 2
                    else if(l == 2 .and. n == 3) then
                        if(ch3(1:1) == "2" .and. ch3(2:2) == "D") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else
                            packed_label(j:j+1) = ch3(1:2)
                            j = j + 2
                        end if
                    else if(l == 2 .and. n == 7) then
                        if(ch3(1:1) == "2" .and. ch3(2:2) == "D") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else
                            packed_label(j:j+1) = ch3(1:2)
                            j = j + 2
                        end if
                    else if(l == 2 .and. n == 4) then
                        if(ch3(1:1) == "3" .and. ch3(2:2) == "F") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "3" .and. ch3(2:2) == "P") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "G") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "D") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "S") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else
                            packed_label(j:j+1) = ch3(1:2)
                            j = j + 2
                        end if
                    else if(l == 2 .and. n == 6) then
                        if(ch3(1:1) == "3" .and. ch3(2:2) == "F") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "3" .and. ch3(2:2) == "P") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "G") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "D") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "S") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else
                            packed_label(j:j+1) = ch3(1:2)
                            j = j + 2
                        end if
                    else if(l == 2 .and. n == 5) then
                        if(ch3(1:1) == "2" .and. ch3(2:2) == "G") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "F") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "D") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else
                            packed_label(j:j+1) = ch3(1:2)
                            j = j + 2
                        end if
                    else if(l == 3 .and. n < 3) then
                        packed_label(j:j+1) = ch3(1:2)
                        j = j + 2
                    else if(l == 3 .and. n > 11) then
                        packed_label(j:j+1) = ch3(1:2)
                        j = j + 2
                    else if(l == 3 .and. n == 3) then
                        if(ch3(1:1) == "2" .and. ch3(2:2) == "H") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "G") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "F") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "D") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else
                            packed_label(j:j+1) = ch3(1:2)
                            j = j + 2
                        end if
                    else if(l == 3 .and. n == 11) then
                        if(ch3(1:1) == "2" .and. ch3(2:2) == "H") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "G") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "F") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "D") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else
                            packed_label(j:j+1) = ch3(1:2)
                            j = j + 2
                        end if
                    else if(l == 3 .and. n == 4) then
                        if(ch3(1:1) == "3" .and. ch3(2:2) == "K") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "3" .and. ch3(2:2) == "I") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "3" .and. ch3(2:2) == "H") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "3" .and. ch3(2:2) == "G") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "3" .and. ch3(2:2) == "F") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "3" .and. ch3(2:2) == "D") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "3" .and. ch3(2:2) == "P") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "L") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "I") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "H") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "G") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "D") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "S") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else
                            packed_label(j:j+1) = ch3(1:2)
                            j = j + 2
                        end if
                    else if(l == 3 .and. n == 10) then
                        if(ch3(1:1) == "3" .and. ch3(2:2) == "K") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "3" .and. ch3(2:2) == "I") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "3" .and. ch3(2:2) == "H") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "3" .and. ch3(2:2) == "G") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "3" .and. ch3(2:2) == "F") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "3" .and. ch3(2:2) == "D") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "3" .and. ch3(2:2) == "P") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "L") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "I") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "H") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "G") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "D") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "1" .and. ch3(2:2) == "S") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else
                            packed_label(j:j+1) = ch3(1:2)
                            j = j + 2
                        end if
                    else if(l == 3 .and. n == 5) then
                        if(ch3(1:1) == "4" .and. ch3(2:2) == "K") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "4" .and. ch3(2:2) == "I") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "4" .and. ch3(2:2) == "H") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "4" .and. ch3(2:2) == "G") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "4" .and. ch3(2:2) == "F") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "4" .and. ch3(2:2) == "D") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "4" .and. ch3(2:2) == "P") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "M") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "L") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "K") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "I") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "H") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "G") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "F") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "D") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "P") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else
                            packed_label(j:j+1) = ch3(1:2)
                            j = j + 2
                        end if
                    else if(l == 3 .and. n == 9) then
                        if(ch3(1:1) == "4" .and. ch3(2:2) == "K") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "4" .and. ch3(2:2) == "I") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "4" .and. ch3(2:2) == "H") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "4" .and. ch3(2:2) == "G") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "4" .and. ch3(2:2) == "F") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "4" .and. ch3(2:2) == "D") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "4" .and. ch3(2:2) == "P") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "M") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "L") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "K") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "I") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "H") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "G") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "F") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "D") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else if(ch3(1:1) == "2" .and. ch3(2:2) == "P") then
                            packed_label(j:j+2) = ch3
                            j = j + 3
                        else
                            packed_label(j:j+1) = ch3(1:2)
                            j = j + 2
                        end if
                    else
                        packed_label(j:j+2) = ch3
                        j = j + 3
                    end if
                end if
            end if
            !
            !  -----  If i=1 or Qi=4l+2 and i<>m,
            !           insert '.'; else _RESULTANTi.  -----
            !
            if(i/=1 .and. previous_q==full1 .and. n==1 .and. subshell_count/=i .and.      &
                icore<2) then
                continue
            else if(i/=1 .and. previous_q==full1 .and. n==full-1 .and.       &
                subshell_count/=i .and. icore<2) then
                continue
            else if(i/=1 .and. previous_q==full1 .and. icore==0 .and.        &
                subshell_count/=i) then
                continue
            else if (i/=1 .and. n/=full .or. subshell_count==i) then
                ch3 = coupling_labels(subshell_count+i-1)
                ch1 = ch3(2:2)
                if (ch1>="a" .and. ch1<="z") ch3(2:2)=char(ichar(ch1)-const)
                k = 2
                if (subshell_count == 1) k = 3
                packed_label(j:j+k) = "_"//ch3(1:k)
                j = j + k + 1
            end if
            if (i == subshell_count) cycle
            packed_label(j:j) = "."
            j = j + 1
        end do
        !
        !>>>>>    Because of a compiler error on the SUN, the following is
        !         needed to have the string printed correctly.
        !         STR = STR(1:J)
        return
    end subroutine pack_ls_subshell_label
end module pack_ls_subshell_label_interface
