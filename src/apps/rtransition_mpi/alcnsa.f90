!***********************************************************************
!                                                                      *
      subroutine alcnsa(jja, jjb, hb1, hb2, hc1, hc2, hm1, &
          hm2, lab, nptr, nsdim, imode)
      use, intrinsic :: iso_fortran_env, only: output_unit
!                                                                      *
!   This  subprogram allocates (IMODE = 1), reallocates (IMODE = 2),   *
!   and deallocates (IMODE = 3) storage for  certain arrays that are   *
!   used in 12_oscl.                                                   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
!                                                                      *
!   Farid A. Parpia.                      Last revision: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use array_memory_management, only: alloc, dalloc, ralloc
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: nsdim
      integer, intent(in) :: imode
      integer, pointer :: jja(:), jjb(:), lab(:), nptr(:)
      real(real64), pointer :: hb1(:), hb2(:), hc1(:), hc2(:), hm1(:), hm2(:)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: newsiz
!-----------------------------------------------
!
      select case (imode)
      case (1)
!
!   Initial array dimension
!
         nsdim = 1
!
!   Allocate storage for arrays
!
         call alloc (jja, nsdim, "JJA", "ALCNSA")
         call alloc (jjb, nsdim, "JJB", "ALCNSA" )
         call alloc (hb1, nsdim, "HB1", "ALCNSA")
         call alloc (hb2, nsdim, "HB2", "ALCNSA")
         call alloc (hc1, nsdim, "HC1", "ALCNSA")
         call alloc (hc2, nsdim, "HC2", "ALCNSA")
         call alloc (hm1, nsdim, "HM1", "ALCNSA")
         call alloc (hm2, nsdim, "HM2", "ALCNSA")
         call alloc (lab, nsdim, "LAB", "ALCNSA" )
         call alloc (nptr, nsdim, "NPTR", "ALCNSA" )
!
      case (2)
!
!   Double the allocation of storage for the arrays
!
         newsiz = 2*nsdim
!
         call ralloc (jja, newsiz, "JJA", "ALCNSA")
         call ralloc (jjb, newsiz, "JJB", "ALCNSA" )
         call ralloc (hb1, newsiz, "HB1", "ALCNSA")
         call ralloc (hb2, newsiz, "HB2", "ALCNSA")
         call ralloc (hc1, newsiz, "HC1", "ALCNSA")
         call ralloc (hc2, newsiz, "HC2", "ALCNSA")
         call ralloc (hm1, newsiz, "HM1", "ALCNSA")
         call ralloc (hm2, newsiz, "HM2", "ALCNSA")
         call ralloc (lab, newsiz, "LAB", "ALCNSA" )
         call ralloc (nptr, newsiz, "NPTR", "ALCNSA" )
!
         nsdim = newsiz
!
      case (3)
!
!   Deallocate the storage for the arrays
!
         call dalloc (jja, "JJA", "ALCNSA")
         call dalloc (jjb, "JJB", "ALCNSA" )
         call dalloc (hb1, "HB1", "ALCNSA")
         call dalloc (hb2, "HB2", "ALCNSA")
         call dalloc (hc1, "HC1", "ALCNSA")
         call dalloc (hc2, "HC2", "ALCNSA")
         call dalloc (hm1, "HM1", "ALCNSA")
         call dalloc (hm2, "HM2", "ALCNSA")
         call dalloc (lab, "LAB", "ALCNSA" )
         call dalloc (nptr, "NPTR", "ALCNSA" )
!
      case default
!
         write (output_unit, *) "ALCNSA: Invalid argument IMODE = ", imode
         stop
!
      end select
!
      return
      end subroutine alcnsa
