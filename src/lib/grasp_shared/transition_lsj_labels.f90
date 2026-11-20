! Stores the LSJ labels and level metadata read for the initial and final transition states.
module transition_lsj_labels
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64), pointer :: rlev_ener_1(:),rlev_ener_2(:)
    character(len=64), pointer :: string_csf1(:)
    character(len=64), pointer :: string_csf2(:)
    integer :: nvectoti,nvectotf,iopen_status1,iopen_status2
    integer, pointer :: lev_pos_1(:),lev_j_1(:),lev_par_1(:)
    integer, pointer :: lev_pos_2(:),lev_j_2(:),lev_par_2(:)
end module transition_lsj_labels
