!
!***********************************************************************
!                                                                      *
      module debug_c
!                                                                      *
!***********************************************************************
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  14:35:02   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: ibug1, ibug2, ibug3, ibug4, ibug5, ibug6
      logical, dimension(5) :: ldbpa
      logical, dimension(5) :: ldbpg
      logical, dimension(30) :: ldbpr
      real(double) :: cutoff  ! used by bioscl
      end module debug_c
