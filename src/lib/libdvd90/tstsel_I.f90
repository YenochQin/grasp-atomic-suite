      module tstsel_i
      implicit none (type, external)
      interface
      logical function tstsel (kpass, nume, neig, iselec, svec, eigval, icv&
         , crite, critc, rowlast, ind, oldval, nncv, incv)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: kpass
      integer, intent(in) :: nume
      integer, intent(in) :: neig
      integer, intent(in) :: iselec(neig)
      real(real64), intent(in) :: svec(kpass*nume)
      real(real64), intent(in) :: eigval(nume)
      integer, intent(inout) :: icv(nume)
      real(real64), intent(in) :: crite
      real(real64), intent(in) :: critc
      real(real64), intent(inout) :: rowlast(neig)
      integer, intent(inout) :: ind(neig)
      real(real64), intent(in) :: oldval(nume)
      integer, intent(inout) :: nncv
      integer, intent(out) :: incv(neig)
!VAST.../MPI/ MYID(IN)
!VAST...Calls: IDAMAX
!...This routine performs I/O.
      end function
      end interface
      end module
