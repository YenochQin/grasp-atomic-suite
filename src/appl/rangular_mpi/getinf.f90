!***********************************************************************
!                                                                      *
      subroutine getinf
!                                                                      *
!   Interactively determines data governing the generation of MCP co-  *
!   efficients.                                                        *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT, GETYN.                                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 19 Dec 1992   *
!   Modified Xinghong He                  Last revision: 30 Jun 1998   *
!                                                                      *
!   File shared (hard link) by mcpmpi, mcpblk                          *
!                                                                      *
!   Updated to treat ICCUT for block                                   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cons_c, only: double, eleven, eps, four, half, one, seven, ten, tenth, three, two, zero
      use orb_c, only: double
      use iccu_c, only: iccut, iccutblk
      use default_c, only: double, inpci, ndef, ndump
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
      use mcp_c, only: diagonal_only, double, max_coulomb_multipole_rank, use_first_order_restriction
      use hblock_c, only: csf_counts_by_block, double, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i, only: getyn
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
!
!   Determine the physical effects specifications
!
      diagonal_only = .false.
      if (diagonal_only) then
         use_first_order_restriction = .false.
         do i = 1,100
            iccut(i) = 0
         end do
      else
         if (ndef /= 0) then
            use_first_order_restriction = .true.
         else
            use_first_order_restriction = .false.
         end if
         if (use_first_order_restriction) then
            write (stderr_unit,*) "The contribution of CSFs 1 -- ICCUT will",&
                            " be treated variationally;"
            write (stderr_unit,*) "the remainder perturbatively; enter ICCUT:"
            do i = 1,symmetry_block_count
              write(stderr_unit,*) "Give ICCUT for block",i
    1         read *, iccut(i)
              if ((iccut(i) <= 1).or.(iccut(i) >= csf_counts_by_block(i))) then
                write (stderr_unit,*) "GETINF: ICCUT must be greater than 1",&
                                " and less than ",csf_counts_by_block(i)
                write (stderr_unit,*) " please reenter ICCUT:"
                goto 1
              end if
            end do
         end if
      end if
      return
      end subroutine getinf
