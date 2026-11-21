module rmcdhf_mcp_memory
    ! Module dependencies
    !   Written by  G. Gaigalas               Vilnius, September 2021  *
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64), pointer :: angular_contribution_coefficients_mcp_31(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_32(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_33(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_34(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_35(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_36(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_37(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_38(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_39(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_40(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_41(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_42(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_43(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_44(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_45(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_46(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_47(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_48(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_49(:)
    real(real64), pointer :: angular_contribution_coefficients_mcp_50(:)
    real(real64), pointer :: coeff(:)
    integer, pointer ::imin1_30(:), imax1_30(:)
    integer, pointer ::imin2_30(:), imax2_30(:)
    integer, pointer ::indx(:), iclmn(:)
    integer, pointer ::iendc_30(:), irow_30(:)
    integer, pointer ::lab_31(:),integral_contribution_counts_by_label_mcp_31(:),contribution_csf_column_indices_mcp_31(:), &
        & contribution_sparse_matrix_indices_mcp_31(:)
    integer, pointer ::lab_32(:),integral_contribution_counts_by_label_mcp_32(:),contribution_csf_column_indices_mcp_32(:), &
        & contribution_sparse_matrix_indices_mcp_32(:)
    integer, pointer ::lab_33(:),integral_contribution_counts_by_label_mcp_33(:),contribution_csf_column_indices_mcp_33(:), &
        & contribution_sparse_matrix_indices_mcp_33(:)
    integer, pointer ::lab_34(:),integral_contribution_counts_by_label_mcp_34(:),contribution_csf_column_indices_mcp_34(:), &
        & contribution_sparse_matrix_indices_mcp_34(:)
    integer, pointer ::lab_35(:),integral_contribution_counts_by_label_mcp_35(:),contribution_csf_column_indices_mcp_35(:), &
        & contribution_sparse_matrix_indices_mcp_35(:)
    integer, pointer ::lab_36(:),integral_contribution_counts_by_label_mcp_36(:),contribution_csf_column_indices_mcp_36(:), &
        & contribution_sparse_matrix_indices_mcp_36(:)
    integer, pointer ::lab_37(:),integral_contribution_counts_by_label_mcp_37(:),contribution_csf_column_indices_mcp_37(:), &
        & contribution_sparse_matrix_indices_mcp_37(:)
    integer, pointer ::lab_38(:),integral_contribution_counts_by_label_mcp_38(:),contribution_csf_column_indices_mcp_38(:), &
        & contribution_sparse_matrix_indices_mcp_38(:)
    integer, pointer ::lab_39(:),integral_contribution_counts_by_label_mcp_39(:),contribution_csf_column_indices_mcp_39(:), &
        & contribution_sparse_matrix_indices_mcp_39(:)
    integer, pointer ::lab_40(:),integral_contribution_counts_by_label_mcp_40(:),contribution_csf_column_indices_mcp_40(:), &
        & contribution_sparse_matrix_indices_mcp_40(:)
    integer, pointer ::lab_41(:),integral_contribution_counts_by_label_mcp_41(:),contribution_csf_column_indices_mcp_41(:), &
        & contribution_sparse_matrix_indices_mcp_41(:)
    integer, pointer ::lab_42(:),integral_contribution_counts_by_label_mcp_42(:),contribution_csf_column_indices_mcp_42(:), &
        & contribution_sparse_matrix_indices_mcp_42(:)
    integer, pointer ::lab_43(:),integral_contribution_counts_by_label_mcp_43(:),contribution_csf_column_indices_mcp_43(:), &
        & contribution_sparse_matrix_indices_mcp_43(:)
    integer, pointer ::lab_44(:),integral_contribution_counts_by_label_mcp_44(:),contribution_csf_column_indices_mcp_44(:), &
        & contribution_sparse_matrix_indices_mcp_44(:)
    integer, pointer ::lab_45(:),integral_contribution_counts_by_label_mcp_45(:),contribution_csf_column_indices_mcp_45(:), &
        & contribution_sparse_matrix_indices_mcp_45(:)
    integer, pointer ::lab_46(:),integral_contribution_counts_by_label_mcp_46(:),contribution_csf_column_indices_mcp_46(:), &
        & contribution_sparse_matrix_indices_mcp_46(:)
    integer, pointer ::lab_47(:),integral_contribution_counts_by_label_mcp_47(:),contribution_csf_column_indices_mcp_47(:), &
        & contribution_sparse_matrix_indices_mcp_47(:)
    integer, pointer ::lab_48(:),integral_contribution_counts_by_label_mcp_48(:),contribution_csf_column_indices_mcp_48(:), &
        & contribution_sparse_matrix_indices_mcp_48(:)
    integer, pointer ::lab_49(:),integral_contribution_counts_by_label_mcp_49(:),contribution_csf_column_indices_mcp_49(:), &
        & contribution_sparse_matrix_indices_mcp_49(:)
    integer, pointer ::lab_50(:),integral_contribution_counts_by_label_mcp_50(:),contribution_csf_column_indices_mcp_50(:), &
        & contribution_sparse_matrix_indices_mcp_50(:)
end module rmcdhf_mcp_memory
