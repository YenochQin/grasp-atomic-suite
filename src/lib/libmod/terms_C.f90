!
!***********************************************************************
!                                                                      *
      module terms_c
!                                                                      *
!***********************************************************************
      implicit none (type, external)
      integer :: itab(39)
      integer :: jtab(40)
      integer :: ntab(357)
      integer :: nrows
      integer, private :: i
      data nrows/ 39/
!
!   A row is defined by a subshell angular momentum and an occupation
!   number
!
!   Each entry ITAB gives the number of terms in a row
!
!   Each entry JTAB gives the starting location -1 of the first triad
!   in a row
!
!   Each triad in NTAB is (v,w,2J+1); here v is the seniority,
!   w resolves any degeneracy in the seniority scheme, and J is the
!   subshell total angular momentum
!
!   Empty subshell or full subshell
!
      data (itab(i),i=1,1)/ 1/
      data (jtab(i),i=1,1)/ 0/
      data (ntab(i),i=1,3)/ 0, 0, 1/
!
!   s, p-   (j = 1/2)
!
      data (itab(i),i=2,2)/ 1/
      data (jtab(i),i=2,2)/ 3/
      data (ntab(i),i=4,6)/ 1, 0, 2/
!
!   p, d-   (j = 3/2)
!
      data (itab(i),i=3,4)/ 1, 2/
      data (jtab(i),i=3,4)/ 6, 9/
      data (ntab(i),i=7,15)/ 1, 0, 4, 0, 0, 1, 2, 0, 5/
!
!  d, f-   (j = 5/2)
!
      data (itab(i),i=5,7)/ 1, 3, 3/
      data (jtab(i),i=5,7)/ 15, 18, 27/
      data (ntab(i),i=16,36)/ 1, 0, 6, 0, 0, 1, 2, 0, 5, 2, 0, 9, 1, 0, 6, 3, 0&
         , 4, 3, 0, 10/
!
!   f, g-   (j = 7/2)
!
      data (itab(i),i=8,11)/ 1, 4, 6, 8/
      data (jtab(i),i=8,11)/ 36, 39, 51, 69/
      data (ntab(i),i=37,93)/ 1, 0, 8, 0, 0, 1, 2, 0, 5, 2, 0, 9, 2, 0, 13, 1, &
         0, 8, 3, 0, 4, 3, 0, 6, 3, 0, 10, 3, 0, 12, 3, 0, 16, 0, 0, 1, 2, 0, 5&
         , 2, 0, 9, 2, 0, 13, 4, 0, 5, 4, 0, 9, 4, 0, 11, 4, 0, 17/
!
!   g, h-   (j = 9/2)
!
      data (itab(i),i=12,16)/ 1, 5, 10, 18, 20/
      data (jtab(i),i=12,16)/ 93, 96, 111, 141, 195/
      data (ntab(i),i=94,255)/ 1, 0, 10, 0, 0, 1, 2, 0, 5, 2, 0, 9, 2, 0, 13, 2&
         , 0, 17, 1, 0, 10, 3, 0, 4, 3, 0, 6, 3, 0, 8, 3, 0, 10, 3, 0, 12, 3, 0&
         , 14, 3, 0, 16, 3, 0, 18, 3, 0, 22, 0, 0, 1, 2, 0, 5, 2, 0, 9, 2, 0, &
         13, 2, 0, 17, 4, 0, 1, 4, 0, 5, 4, 0, 7, 4, 0, 9, 4, 1, 9, 4, 0, 11, 4&
         , 0, 13, 4, 1, 13, 4, 0, 15, 4, 0, 17, 4, 0, 19, 4, 0, 21, 4, 0, 25, 1&
         , 0, 10, 3, 0, 4, 3, 0, 6, 3, 0, 8, 3, 0, 10, 3, 0, 12, 3, 0, 14, 3, 0&
         , 16, 3, 0, 18, 3, 0, 22, 5, 0, 2, 5, 0, 6, 5, 0, 8, 5, 0, 10, 5, 0, &
         12, 5, 0, 14, 5, 0, 16, 5, 0, 18, 5, 0, 20, 5, 0, 26/
!
!   h, i-   (j = 11/2)
!
!   h, i-   (j = 11/2)
!
!   First two rows only
!
      data (itab(i),i=17,18)/ 1, 6/
      data (jtab(i),i=17,19)/ 255, 258, 277/
      data (ntab(i),i=256,276)/ 1, 0, 12, 0, 0, 1, 2, 0, 5, 2, 0, 9, 2, 0, 13, &
         2, 0, 17, 2, 0, 21/
!
!   i, k-   (j = 13/2)
!
!   First two rows only
!
      data (itab(i),i=23,24)/ 1, 7/
      data (jtab(i),i=23,25)/ 276, 279, 301/
      data (ntab(i),i=277,300)/ 1, 0, 14, 0, 0, 1, 2, 0, 5, 2, 0, 9, 2, 0, 13, &
         2, 0, 17, 2, 0, 21, 2, 0, 25/
!
!   k, l-   (j = 15/2)
!
!   First two rows only
!
      data (itab(i),i=30,31)/ 1, 8/
      data (jtab(i),i=30,32)/ 300, 303, 328/
      data (ntab(i),i=301,327)/ 1, 0, 16, 0, 0, 1, 2, 0, 5, 2, 0, 9, 2, 0, 13, &
         2, 0, 17, 2, 0, 21, 2, 0, 25, 2, 0, 29/
!
!   l, m-   (j = 17/2)
!
!   First two rows only
!
      data (itab(i),i=38,39)/ 1, 9/
      data (jtab(i),i=38,40)/ 327, 330, 358/
      data (ntab(i),i=328,357)/ 1, 0, 18, 0, 0, 1, 2, 0, 5, 2, 0, 9, 2, 0, 13, &
         2, 0, 17, 2, 0, 21, 2, 0, 25, 2, 0, 29, 2, 0, 33/
      end module terms_c
