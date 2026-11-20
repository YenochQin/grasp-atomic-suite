! Stores intermediate-term index ranges by occupation parity for 2j = 9, using term indices 301 through 306.
module jj_reduced_coefficients_9
    implicit none (type, external)

    ! Module state
    integer :: first_same_parity_terms_two_j_9(6), last_same_parity_terms_two_j_9(6), first_opposite_parity_terms_two_j_9(6), &
        last_opposite_parity_terms_two_j_9(6)

    ! Initial values
    data first_same_parity_terms_two_j_9/301,5*302/
    data last_same_parity_terms_two_j_9/301,5*306/
    data first_opposite_parity_terms_two_j_9/302,5*301/
    data last_opposite_parity_terms_two_j_9/306,5*301/
end module jj_reduced_coefficients_9
