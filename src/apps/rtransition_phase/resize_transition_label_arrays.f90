!***********************************************************************
!                                                                      *
subroutine resize_transition_label_arrays(ibeg, ilab, ilast, ileft, iptcsf, iright, lblint, &
    lldim, imode)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    !                                                                      *
    !   This  subprogram allocates (IMODE = 1), reallocates (IMODE = 2),   *
    !   and deallocates (IMODE = 3) storage for  certain arrays that are   *
    !   local to SUBROUTINE TRSORT.                                        *
    !                                                                      *
    !   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
    !                                                                      *
    !   Written by Farid A. Parpia            Last revision: 28 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    use array_memory_management, only: alloc, dalloc, ralloc
    implicit none (type, external)

    ! Arguments
    integer  :: lldim
    integer , intent(in) :: imode
    integer, pointer :: ibeg(:), ilab(:), ilast(:), ileft(:), &
        iptcsf(:), iright(:), lblint(:)

    ! Local declarations
    integer :: newsiz
    !-----------------------------------------------
    !
    !
    select case (imode)
    case (1)
        !
        !   Initial array dimension
        !
        lldim = 64
        !
        !   Allocate storage for arrays
        !
        call alloc (ibeg, lldim, "IBEG", "ALCLLA")
        call alloc (ilab, lldim, "ILAB", "ALCLLA")
        call alloc (ilast, lldim, "ILAST", "ALCLLA")
        call alloc (ileft, lldim, "ILEFT", "ALCLLA")
        call alloc (iptcsf, lldim, "IPTCSF", "ALCLLA")
        call alloc (iright, lldim, "IRIGHT", "ALCLLA")
        call alloc (lblint, lldim, "LBLINT", "ALCLLA")
    !
    case (2)
        !
        !   Double the allocation of storage for the arrays
        !
        newsiz = 2*lldim
        !
        call ralloc (ibeg, newsiz, "IBEG", "ALCLLA")
        call ralloc (ilab, newsiz, "ILAB", "ALCLLA")
        call ralloc (ilast, newsiz, "ILAST", "ALCLLA")
        call ralloc (ileft, newsiz, "ILEFT", "ALCLLA")
        call ralloc (iptcsf, newsiz, "IPTCSF", "ALCLLA")
        call ralloc (iright, newsiz, "IRIGHT", "ALCLLA")
        call ralloc (lblint, newsiz, "LBLINT", "ALCLLA")
        !
        lldim = newsiz
    !
    case (3)
        !
        !   Deallocate the storage for the arrays
        !
        call dalloc (ibeg, "IBEG", "ALCLLA")
        call dalloc (ilab, "ILAB", "ALCLLA")
        call dalloc (ilast, "ILAST", "ALCLLA")
        call dalloc (ileft, "ILEFT", "ALCLLA")
        call dalloc (iptcsf, "IPTCSF", "ALCLLA")
        call dalloc (iright, "IRIGHT", "ALCLLA")
        call dalloc (lblint, "LBLINT", "ALCLLA")
    !
    case default
        !
        write (output_unit, *) "ALCLLA: Invalid argument IMODE = ", imode
        stop
        !
    end select
    !
    return
end subroutine resize_transition_label_arrays
