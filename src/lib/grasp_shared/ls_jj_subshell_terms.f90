! Defines LS-coupled subshell terms and the records used to look up LS–jj transformation coefficients.
      module ls_jj_subshell_terms
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none (type, external)
      private :: int64
      ! Subshell orbital l and term labels w, 2Q, 2L, and 2S.
      type :: subshell_term_ls
         integer :: subshell_l
         integer :: w
         integer :: subshell_two_quasispin
         integer :: subshell_two_l
         integer :: subshell_two_s
      end type subshell_term_ls
      ! A lookup record represents factor * sqrt(nom / denom).
      type :: ls_jj_me
         integer       :: w, q, l, s, j, nm, qm, jm, qp, jp
         integer       :: factor
         integer(int64) :: nom, denom
      end type ls_jj_me
      ! Number of tabulated coefficients for each subshell and occupation.
      integer, parameter :: ls_jj_number_p3 =  10, ls_jj_number_p4  =    9, &
                            ls_jj_number_p5 =   2, ls_jj_number_p6  =    1, &
                            ls_jj_number_d3 =  65, ls_jj_number_d4  =  166, &
                            ls_jj_number_d5 = 184, ls_jj_number_d6  =  166, &
                            ls_jj_number_d7 =  65, ls_jj_number_d8  =   19, &
                            ls_jj_number_d9 =   2, ls_jj_number_d10 =    1, &
                            ls_jj_number_f3 = 216, ls_jj_number_f4  = 1210, &
                            ls_jj_number_f5 =3799, ls_jj_number_f6  = 7313, &
                            ls_jj_number_f7 =8003
   ! Antisymmetric LS terms for each nonrelativistic subshell.
   type(subshell_term_ls), parameter ::   &
      term_ls_s(1:2) =[                                       &
      subshell_term_ls(0, 0, 0, 0, 1),                    &
      subshell_term_ls(0, 0, 1, 0, 0)  ]
   type(subshell_term_ls), parameter ::   &
      term_ls_p(1:6) =[                                       &
      subshell_term_ls(1, 0, 0, 0, 3),                    &
      subshell_term_ls(1, 0, 2, 2, 1),                    &
      subshell_term_ls(1, 0, 0, 4, 1),                    &
      subshell_term_ls(1, 0, 3, 0, 0),                    &
      subshell_term_ls(1, 0, 1, 2, 2),                    &
      subshell_term_ls(1, 0, 1, 4, 0)  ]
   type(subshell_term_ls), parameter ::  &
      term_ls_d(1:32) =[                                       &
      subshell_term_ls(2, 0, 0, 0, 5),                    &
      subshell_term_ls(2, 0, 0, 0, 1),                    &
      subshell_term_ls(2, 0, 2, 2, 3),                    &
      subshell_term_ls(2, 0, 2, 2, 1),                    &
      subshell_term_ls(2, 0, 4, 4, 1),                    &
      subshell_term_ls(2, 0, 2, 4, 1),                    &
      subshell_term_ls(2, 0, 0, 4, 3),                    &
      subshell_term_ls(2, 0, 0, 4, 1),                    &
      subshell_term_ls(2, 0, 2, 6, 3),                    &
      subshell_term_ls(2, 0, 2, 6, 1),                    &
      subshell_term_ls(2, 0, 0, 6, 1),                    &
      subshell_term_ls(2, 0, 2, 8, 1),                    &
      subshell_term_ls(2, 0, 0, 8, 3),                    &
      subshell_term_ls(2, 0, 0, 8, 1),                    &
      subshell_term_ls(2, 0, 2,10, 1),                    &
      subshell_term_ls(2, 0, 0,12, 1),                    &
      subshell_term_ls(2, 0, 5, 0, 0),                    &
      subshell_term_ls(2, 0, 1, 0, 0),                    &
      subshell_term_ls(2, 0, 3, 2, 2),                    &
      subshell_term_ls(2, 0, 1, 2, 2),                    &
      subshell_term_ls(2, 0, 1, 4, 4),                    &
      subshell_term_ls(2, 0, 1, 4, 2),                    &
      subshell_term_ls(2, 0, 3, 4, 0),                    &
      subshell_term_ls(2, 0, 1, 4, 0),                    &
      subshell_term_ls(2, 0, 3, 6, 2),                    &
      subshell_term_ls(2, 0, 1, 6, 2),                    &
      subshell_term_ls(2, 0, 1, 6, 0),                    &
      subshell_term_ls(2, 0, 1, 8, 2),                    &
      subshell_term_ls(2, 0, 3, 8, 0),                    &
      subshell_term_ls(2, 0, 1, 8, 0),                    &
      subshell_term_ls(2, 0, 1,10, 2),                    &
      subshell_term_ls(2, 0, 1,12, 0)  ]
   type(subshell_term_ls), parameter :: &
      term_ls_f(1:238) =[                                       &
      subshell_term_ls(3, 0, 0, 0, 7),                    &
      subshell_term_ls(3, 0, 2, 2, 5),                    &
      subshell_term_ls(3, 0, 0, 4, 5),                    &
      subshell_term_ls(3, 0, 2, 6, 5),                    &
      subshell_term_ls(3, 0, 0, 8, 5),                    &
      subshell_term_ls(3, 0, 2,10, 5),                    &
      subshell_term_ls(3, 0, 0,12, 5),                    &
      subshell_term_ls(3, 1, 4, 0, 3),                    &
      subshell_term_ls(3, 2, 0, 0, 3),                    &
      subshell_term_ls(3, 1, 2, 2, 3),                    &
      subshell_term_ls(3, 2, 2, 2, 3),                    &
      subshell_term_ls(3, 1, 4, 4, 3),                    &
      subshell_term_ls(3, 2, 2, 4, 3),                    &
      subshell_term_ls(3, 3, 2, 4, 3),                    &
      subshell_term_ls(3, 4, 0, 4, 3),                    &
      subshell_term_ls(3, 5, 0, 4, 3),                    &
      subshell_term_ls(3, 6, 0, 4, 3),                    &
      subshell_term_ls(3, 1, 4, 6, 3),                    &
      subshell_term_ls(3, 2, 2, 6, 3),                    &
      subshell_term_ls(3, 3, 2, 6, 3),                    &
      subshell_term_ls(3, 4, 2, 6, 3),                    &
      subshell_term_ls(3, 5, 0, 6, 3),                    &
      subshell_term_ls(3, 1, 4, 8, 3),                    &
      subshell_term_ls(3, 2, 2, 8, 3),                    &
      subshell_term_ls(3, 3, 2, 8, 3),                    &
      subshell_term_ls(3, 4, 2, 8, 3),                    &
      subshell_term_ls(3, 5, 0, 8, 3),                    &
      subshell_term_ls(3, 6, 0, 8, 3),                    &
      subshell_term_ls(3, 7, 0, 8, 3),                    &
      subshell_term_ls(3, 1, 2,10, 3),                    &
      subshell_term_ls(3, 2, 2,10, 3),                    &
      subshell_term_ls(3, 3, 2,10, 3),                    &
      subshell_term_ls(3, 4, 0,10, 3),                    &
      subshell_term_ls(3, 5, 0,10, 3),                    &
      subshell_term_ls(3, 1, 4,12, 3),                    &
      subshell_term_ls(3, 2, 2,12, 3),                    &
      subshell_term_ls(3, 3, 2,12, 3),                    &
      subshell_term_ls(3, 4, 0,12, 3),                    &
      subshell_term_ls(3, 5, 0,12, 3),                    &
      subshell_term_ls(3, 1, 2,14, 3),                    &
      subshell_term_ls(3, 2, 2,14, 3),                    &
      subshell_term_ls(3, 3, 0,14, 3),                    &
      subshell_term_ls(3, 1, 2,16, 3),                    &
      subshell_term_ls(3, 2, 0,16, 3),                    &
      subshell_term_ls(3, 3, 0,16, 3),                    &
      subshell_term_ls(3, 0, 2,18, 3),                    &
      subshell_term_ls(3, 0, 0,20, 3),                    &
      subshell_term_ls(3, 1, 0, 0, 1),                    &
      subshell_term_ls(3, 2, 0, 0, 1),                    &
      subshell_term_ls(3, 1, 4, 2, 1),                    &
      subshell_term_ls(3, 2, 2, 2, 1),                    &
      subshell_term_ls(3, 3, 2, 2, 1),                    &
      subshell_term_ls(3, 4, 2, 2, 1),                    &
      subshell_term_ls(3, 5, 0, 2, 1),                    &
      subshell_term_ls(3, 1, 4, 4, 1),                    &
      subshell_term_ls(3, 2, 4, 4, 1),                    &
      subshell_term_ls(3, 3, 2, 4, 1),                    &
      subshell_term_ls(3, 4, 2, 4, 1),                    &
      subshell_term_ls(3, 5, 2, 4, 1),                    &
      subshell_term_ls(3, 6, 0, 4, 1),                    &
      subshell_term_ls(3, 7, 0, 4, 1),                    &
      subshell_term_ls(3, 1, 6, 6, 1),                    &
      subshell_term_ls(3, 2, 4, 6, 1),                    &
      subshell_term_ls(3, 3, 2, 6, 1),                    &
      subshell_term_ls(3, 4, 2, 6, 1),                    &
      subshell_term_ls(3, 5, 2, 6, 1),                    &
      subshell_term_ls(3, 6, 2, 6, 1),                    &
      subshell_term_ls(3, 7, 2, 6, 1),                    &
      subshell_term_ls(3, 8, 0, 6, 1),                    &
      subshell_term_ls(3, 9, 0, 6, 1),                    &
      subshell_term_ls(3,10, 0, 6, 1),                    &
      subshell_term_ls(3, 1, 4, 8, 1),                    &
      subshell_term_ls(3, 2, 4, 8, 1),                    &
      subshell_term_ls(3, 3, 2, 8, 1),                    &
      subshell_term_ls(3, 4, 2, 8, 1),                    &
      subshell_term_ls(3, 5, 2, 8, 1),                    &
      subshell_term_ls(3, 6, 2, 8, 1),                    &
      subshell_term_ls(3, 7, 0, 8, 1),                    &
      subshell_term_ls(3, 8, 0, 8, 1),                    &
      subshell_term_ls(3, 9, 0, 8, 1),                    &
      subshell_term_ls(3,10, 0, 8, 1),                    &
      subshell_term_ls(3, 1, 4,10, 1),                    &
      subshell_term_ls(3, 2, 4,10, 1),                    &
      subshell_term_ls(3, 3, 2,10, 1),                    &
      subshell_term_ls(3, 4, 2,10, 1),                    &
      subshell_term_ls(3, 5, 2,10, 1),                    &
      subshell_term_ls(3, 6, 2,10, 1),                    &
      subshell_term_ls(3, 7, 2,10, 1),                    &
      subshell_term_ls(3, 8, 0,10, 1),                    &
      subshell_term_ls(3, 9, 0,10, 1),                    &
      subshell_term_ls(3, 1, 4,12, 1),                    &
      subshell_term_ls(3, 2, 2,12, 1),                    &
      subshell_term_ls(3, 3, 2,12, 1),                    &
      subshell_term_ls(3, 4, 2,12, 1),                    &
      subshell_term_ls(3, 5, 2,12, 1),                    &
      subshell_term_ls(3, 6, 0,12, 1),                    &
      subshell_term_ls(3, 7, 0,12, 1),                    &
      subshell_term_ls(3, 8, 0,12, 1),                    &
      subshell_term_ls(3, 9, 0,12, 1),                    &
      subshell_term_ls(3, 1, 4,14, 1),                    &
      subshell_term_ls(3, 2, 2,14, 1),                    &
      subshell_term_ls(3, 3, 2,14, 1),                    &
      subshell_term_ls(3, 4, 2,14, 1),                    &
      subshell_term_ls(3, 5, 2,14, 1),                    &
      subshell_term_ls(3, 6, 0,14, 1),                    &
      subshell_term_ls(3, 7, 0,14, 1),                    &
      subshell_term_ls(3, 1, 4,16, 1),                    &
      subshell_term_ls(3, 2, 2,16, 1),                    &
      subshell_term_ls(3, 3, 2,16, 1),                    &
      subshell_term_ls(3, 4, 0,16, 1),                    &
      subshell_term_ls(3, 5, 0,16, 1),                    &
      subshell_term_ls(3, 1, 2,18, 1),                    &
      subshell_term_ls(3, 2, 2,18, 1),                    &
      subshell_term_ls(3, 3, 0,18, 1),                    &
      subshell_term_ls(3, 4, 0,18, 1),                    &
      subshell_term_ls(3, 1, 2,20, 1),                    &
      subshell_term_ls(3, 2, 0,20, 1),                    &
      subshell_term_ls(3, 0, 2,22, 1),                    &
      subshell_term_ls(3, 0, 0,24, 1),                    &
      subshell_term_ls(3, 0, 1, 6, 6),                    &
      subshell_term_ls(3, 1, 3, 4, 4),                    &
      subshell_term_ls(3, 2, 1, 4, 4),                    &
      subshell_term_ls(3, 3, 1, 4, 4),                    &
      subshell_term_ls(3, 1, 3, 6, 4),                    &
      subshell_term_ls(3, 2, 1, 6, 4),                    &
      subshell_term_ls(3, 1, 3, 8, 4),                    &
      subshell_term_ls(3, 2, 1, 8, 4),                    &
      subshell_term_ls(3, 3, 1, 8, 4),                    &
      subshell_term_ls(3, 0, 1, 2, 4),                    &
      subshell_term_ls(3, 1, 1,10, 4),                    &
      subshell_term_ls(3, 2, 1,10, 4),                    &
      subshell_term_ls(3, 0, 3, 0, 4),                    &
      subshell_term_ls(3, 1, 3,12, 4),                    &
      subshell_term_ls(3, 2, 1,12, 4),                    &
      subshell_term_ls(3, 0, 1,14, 4),                    &
      subshell_term_ls(3, 0, 1,16, 4),                    &
      subshell_term_ls(3, 1, 5, 6, 2),                    &
      subshell_term_ls(3, 2, 3, 6, 2),                    &
      subshell_term_ls(3, 6, 1, 6, 2),                    &
      subshell_term_ls(3, 8, 1, 6, 2),                    &
      subshell_term_ls(3, 1, 3, 4, 2),                    &
      subshell_term_ls(3, 2, 3, 4, 2),                    &
      subshell_term_ls(3, 3, 1, 4, 2),                    &
      subshell_term_ls(3, 4, 1, 4, 2),                    &
      subshell_term_ls(3, 3, 3, 6, 2),                    &
      subshell_term_ls(3, 5, 1, 6, 2),                    &
      subshell_term_ls(3, 1, 3, 8, 2),                    &
      subshell_term_ls(3, 2, 3, 8, 2),                    &
      subshell_term_ls(3, 4, 1, 8, 2),                    &
      subshell_term_ls(3, 5, 1, 8, 2),                    &
      subshell_term_ls(3, 5, 1, 4, 2),                    &
      subshell_term_ls(3, 4, 3, 6, 2),                    &
      subshell_term_ls(3, 7, 1, 6, 2),                    &
      subshell_term_ls(3, 9, 1, 6, 2),                    &
      subshell_term_ls(3, 3, 3, 8, 2),                    &
      subshell_term_ls(3, 6, 1, 8, 2),                    &
      subshell_term_ls(3, 7, 1, 8, 2),                    &
      subshell_term_ls(3, 1, 5, 2, 2),                    &
      subshell_term_ls(3, 2, 3, 2, 2),                    &
      subshell_term_ls(3, 3, 3, 2, 2),                    &
      subshell_term_ls(3, 1, 5,10, 2),                    &
      subshell_term_ls(3, 2, 3,10, 2),                    &
      subshell_term_ls(3, 3, 3,10, 2),                    &
      subshell_term_ls(3, 4, 3,10, 2),                    &
      subshell_term_ls(3, 4, 1, 2, 2),                    &
      subshell_term_ls(3, 5, 1,10, 2),                    &
      subshell_term_ls(3, 6, 1,10, 2),                    &
      subshell_term_ls(3, 5, 1, 2, 2),                    &
      subshell_term_ls(3, 6, 1, 2, 2),                    &
      subshell_term_ls(3, 7, 1,10, 2),                    &
      subshell_term_ls(3, 8, 1,10, 2),                    &
      subshell_term_ls(3, 9, 1,10, 2),                    &
      subshell_term_ls(3, 1, 3,12, 2),                    &
      subshell_term_ls(3, 2, 3,12, 2),                    &
      subshell_term_ls(3, 3, 1,12, 2),                    &
      subshell_term_ls(3, 3, 1,12, 2),                    &
      subshell_term_ls(3, 5, 1,12, 2),                    &
      subshell_term_ls(3, 6, 1,12, 2),                    &
      subshell_term_ls(3, 1, 3,14, 2),                    &
      subshell_term_ls(3, 2, 3,14, 2),                    &
      subshell_term_ls(3, 3, 1,14, 2),                    &
      subshell_term_ls(3, 4, 1,14, 2),                    &
      subshell_term_ls(3, 5, 1,14, 2),                    &
      subshell_term_ls(3, 6, 1,14, 2),                    &
      subshell_term_ls(3, 1, 3,16, 2),                    &
      subshell_term_ls(3, 2, 1,16, 2),                    &
      subshell_term_ls(3, 3, 1,16, 2),                    &
      subshell_term_ls(3, 1, 3,18, 2),                    &
      subshell_term_ls(3, 2, 1,18, 2),                    &
      subshell_term_ls(3, 3, 1,18, 2),                    &
      subshell_term_ls(3, 0, 1,20, 2),                    &
      subshell_term_ls(3, 0, 1,22, 2),                    &
      subshell_term_ls(3, 2, 1, 6, 0),                    &
      subshell_term_ls(3, 3, 1, 6, 0),                    &
      subshell_term_ls(3, 4, 1, 6, 0),                    &
      subshell_term_ls(3, 1, 5, 4, 0),                    &
      subshell_term_ls(3, 2, 3, 4, 0),                    &
      subshell_term_ls(3, 3, 3, 4, 0),                    &
      subshell_term_ls(3, 1, 3, 6, 0),                    &
      subshell_term_ls(3, 1, 5, 8, 0),                    &
      subshell_term_ls(3, 2, 3, 8, 0),                    &
      subshell_term_ls(3, 3, 3, 8, 0),                    &
      subshell_term_ls(3, 5, 1, 4, 0),                    &
      subshell_term_ls(3, 5, 1, 8, 0),                    &
      subshell_term_ls(3, 6, 1, 4, 0),                    &
      subshell_term_ls(3, 6, 1, 8, 0),                    &
      subshell_term_ls(3, 7, 1, 8, 0),                    &
      subshell_term_ls(3, 8, 1, 8, 0),                    &
      subshell_term_ls(3, 4, 3, 4, 0),                    &
      subshell_term_ls(3, 4, 3, 8, 0),                    &
      subshell_term_ls(3, 1, 3,10, 0),                    &
      subshell_term_ls(3, 2, 3,10, 0),                    &
      subshell_term_ls(3, 0, 1, 2, 0),                    &
      subshell_term_ls(3, 3, 1,10, 0),                    &
      subshell_term_ls(3, 4, 1,10, 0),                    &
      subshell_term_ls(3, 1, 7, 0, 0),                    &
      subshell_term_ls(3, 1, 5,12, 0),                    &
      subshell_term_ls(3, 2, 3, 0, 0),                    &
      subshell_term_ls(3, 2, 3,12, 0),                    &
      subshell_term_ls(3, 3, 3,12, 0),                    &
      subshell_term_ls(3, 3, 1, 0, 0),                    &
      subshell_term_ls(3, 4, 1,12, 0),                    &
      subshell_term_ls(3, 5, 1,12, 0),                    &
      subshell_term_ls(3, 4, 1, 0, 0),                    &
      subshell_term_ls(3, 6, 1,12, 0),                    &
      subshell_term_ls(3, 7, 1,12, 0),                    &
      subshell_term_ls(3, 1, 3,14, 0),                    &
      subshell_term_ls(3, 2, 1,14, 0),                    &
      subshell_term_ls(3, 3, 1,14, 0),                    &
      subshell_term_ls(3, 1, 3,16, 0),                    &
      subshell_term_ls(3, 2, 3,16, 0),                    &
      subshell_term_ls(3, 3, 1,16, 0),                    &
      subshell_term_ls(3, 4, 1,16, 0),                    &
      subshell_term_ls(3, 1, 1,18, 0),                    &
      subshell_term_ls(3, 2, 1,18, 0),                    &
      subshell_term_ls(3, 1, 3,20, 0),                    &
      subshell_term_ls(3, 2, 1,20, 0),                    &
      subshell_term_ls(3, 0, 1,24, 0)  ]
   type(subshell_term_ls), parameter ::   &
      term_ls_g1(1:1) =[                                      &
      subshell_term_ls(4, 0, 8, 8, 1)  ]
   type(subshell_term_ls), parameter ::   &
      term_ls_g2(1:9) =[                                      &
      subshell_term_ls(4, 0, 9, 0, 0),                    &
      subshell_term_ls(4, 0, 7, 2, 2),                    &
      subshell_term_ls(4, 0, 7, 4, 0),                    &
      subshell_term_ls(4, 0, 7, 6, 2),                    &
      subshell_term_ls(4, 0, 7, 8, 0),                    &
      subshell_term_ls(4, 0, 7,10, 2),                    &
      subshell_term_ls(4, 0, 7,12, 0),                    &
      subshell_term_ls(4, 0, 7,14, 2),                    &
      subshell_term_ls(4, 0, 7,16, 0)  ]
   type(subshell_term_ls), parameter ::   &
      term_ls_h1(1:1) =[                                      &
      subshell_term_ls(5, 0,10,10, 1)  ]
   type(subshell_term_ls), parameter ::  &
      term_ls_h2(1:11) =[                                      &
      subshell_term_ls(5, 0,11, 0, 0),                    &
      subshell_term_ls(5, 0, 9, 2, 2),                    &
      subshell_term_ls(5, 0, 9, 4, 0),                    &
      subshell_term_ls(5, 0, 9, 6, 2),                    &
      subshell_term_ls(5, 0, 9, 8, 0),                    &
      subshell_term_ls(5, 0, 9,10, 2),                    &
      subshell_term_ls(5, 0, 9,12, 0),                    &
      subshell_term_ls(5, 0, 9,14, 2),                    &
      subshell_term_ls(5, 0, 9,16, 0),                    &
      subshell_term_ls(5, 0, 9,18, 2),                    &
      subshell_term_ls(5, 0, 9,20, 0)  ]
   type(subshell_term_ls), parameter ::   &
      term_ls_i1(1:1) =[                                      &
      subshell_term_ls(6, 0,12,12, 1)  ]
   type(subshell_term_ls), parameter ::  &
      term_ls_i2(1:13) =[                                      &
      subshell_term_ls(6, 0,13, 0, 0),                    &
      subshell_term_ls(6, 0,11, 2, 2),                    &
      subshell_term_ls(6, 0,11, 4, 0),                    &
      subshell_term_ls(6, 0,11, 6, 2),                    &
      subshell_term_ls(6, 0,11, 8, 0),                    &
      subshell_term_ls(6, 0,11,10, 2),                    &
      subshell_term_ls(6, 0,11,12, 0),                    &
      subshell_term_ls(6, 0,11,14, 2),                    &
      subshell_term_ls(6, 0,11,16, 0),                    &
      subshell_term_ls(6, 0,11,18, 2),                    &
      subshell_term_ls(6, 0,11,20, 0),                    &
      subshell_term_ls(6, 0,11,22, 2),                    &
      subshell_term_ls(6, 0,11,24, 0)  ]
   type(subshell_term_ls), parameter ::   &
      term_ls_k1(1:1) =[                                      &
      subshell_term_ls(7, 0,14,14, 1)  ]
   type(subshell_term_ls), parameter ::  &
      term_ls_k2(1:15) =[                                      &
      subshell_term_ls(7, 0,15, 0, 0),                    &
      subshell_term_ls(7, 0,13, 2, 2),                    &
      subshell_term_ls(7, 0,13, 4, 0),                    &
      subshell_term_ls(7, 0,13, 6, 2),                    &
      subshell_term_ls(7, 0,13, 8, 0),                    &
      subshell_term_ls(7, 0,13,10, 2),                    &
      subshell_term_ls(7, 0,13,12, 0),                    &
      subshell_term_ls(7, 0,13,14, 2),                    &
      subshell_term_ls(7, 0,13,16, 0),                    &
      subshell_term_ls(7, 0,13,18, 2),                    &
      subshell_term_ls(7, 0,13,20, 0),                    &
      subshell_term_ls(7, 0,13,22, 2),                    &
      subshell_term_ls(7, 0,13,24, 0),                    &
      subshell_term_ls(7, 0,13,26, 2),                    &
      subshell_term_ls(7, 0,13,28, 0) ]
   type(subshell_term_ls), parameter ::   &
      term_ls_l1(1:1) =[                                      &
      subshell_term_ls(8, 0,16,16, 1)  ]
   type(subshell_term_ls), parameter ::  &
      term_ls_l2(1:17) =[                                      &
      subshell_term_ls(8, 0,17, 0, 0),                    &
      subshell_term_ls(8, 0,15, 2, 2),                    &
      subshell_term_ls(8, 0,15, 4, 0),                    &
      subshell_term_ls(8, 0,15, 6, 2),                    &
      subshell_term_ls(8, 0,15, 8, 0),                    &
      subshell_term_ls(8, 0,15,10, 2),                    &
      subshell_term_ls(8, 0,15,12, 0),                    &
      subshell_term_ls(8, 0,15,14, 2),                    &
      subshell_term_ls(8, 0,15,16, 0),                    &
      subshell_term_ls(8, 0,15,18, 2),                    &
      subshell_term_ls(8, 0,15,20, 0),                    &
      subshell_term_ls(8, 0,15,22, 2),                    &
      subshell_term_ls(8, 0,15,24, 0),                    &
      subshell_term_ls(8, 0,15,26, 2),                    &
      subshell_term_ls(8, 0,15,28, 0),                    &
      subshell_term_ls(8, 0,15,30, 2),                    &
      subshell_term_ls(8, 0,15,32, 0) ]
   type(subshell_term_ls), parameter ::   &
      term_ls_m1(1:1) =[                                      &
      subshell_term_ls(9, 0,18,18, 1)  ]
   type(subshell_term_ls), parameter ::  &
      term_ls_m2(1:19) =[                                      &
      subshell_term_ls(9, 0,19, 0, 0),                    &
      subshell_term_ls(9, 0,17, 2, 2),                    &
      subshell_term_ls(9, 0,17, 4, 0),                    &
      subshell_term_ls(9, 0,17, 6, 2),                    &
      subshell_term_ls(9, 0,17, 8, 0),                    &
      subshell_term_ls(9, 0,17,10, 2),                    &
      subshell_term_ls(9, 0,17,12, 0),                    &
      subshell_term_ls(9, 0,17,14, 2),                    &
      subshell_term_ls(9, 0,17,16, 0),                    &
      subshell_term_ls(9, 0,17,18, 2),                    &
      subshell_term_ls(9, 0,17,20, 0),                    &
      subshell_term_ls(9, 0,17,22, 2),                    &
      subshell_term_ls(9, 0,17,24, 0),                    &
      subshell_term_ls(9, 0,17,26, 2),                    &
      subshell_term_ls(9, 0,17,28, 0),                    &
      subshell_term_ls(9, 0,17,30, 2),                    &
      subshell_term_ls(9, 0,17,32, 0),                    &
      subshell_term_ls(9, 0,17,34, 2),                    &
      subshell_term_ls(9, 0,17,36, 0) ]
      end module ls_jj_subshell_terms
