PROGRAM timeing_loop
    IMPLICIT NONE
    INTEGER, PARAMETER :: MK = KIND(1d0)
    INTEGER, PARAMETER :: m = 3
    INTEGER(MK) :: i, j, k, n, n_max
    REAL(MK), DIMENSION(m,m,m) :: f
    REAL ts , te , t_initial, t_flip_i_k, t_flip_i_k_j_k

    t_initial = 0d0
    t_flip_i_k = 0d0
    n_max = 10000000

    CALL INITIAL

    DO n = 1,n_max
        f = 1d0
        CALL CPU_TIME( ts )
        CALL INITIAL
        CALL CPU_TIME( te )
        t_initial = t_initial+te-ts
        !WRITE(*,*) "(sub INITIAL) Wall time = ", t

        f = 1d0
        CALL CPU_TIME( ts )
        CALL FLIP_i_k
        CALL CPU_TIME( te )
        t_flip_i_k = t_flip_i_k+te-ts

        f = 1d0
        CALL CPU_TIME( ts )
        CALL FLIP_i_k_j_k
        CALL CPU_TIME( te )
        t_flip_i_k_j_k = t_flip_i_k_j_k+te-ts

    END DO
    WRITE(*,*) "n_max: ", n_max
    WRITE(*,*) "(sub INITIAL) Avg, Wall time = ", t_initial/n_max, "(tot: ", t_initial,")"
    WRITE(*,*) "(sub FLIP_i_k) Avg, Wall time = ", t_flip_i_k/n_max, "(tot: ", t_flip_i_k, ")"
    WRITE(*,*) "(sub FLIP_i_k_j_k) Avg, Wall time = ", t_flip_i_k_j_k/n_max, "(tot: ", t_flip_i_k_j_k, ")"
    contains
    SUBROUTINE INITIAL
        DO i =1,m
            Do j =1,m
                DO k=1,m
                    f ( i , j , k) = f ( i , j , k) * 2d0
                ENDDO
            ENDDO
        ENDDO
    END SUBROUTINE INITIAL

    SUBROUTINE FLIP_i_k
        DO k =1,m
            Do j =1,m
                DO i=1,m
                    f ( i , j , k) = f ( i , j , k) * 2d0
                ENDDO
            ENDDO
        ENDDO
    END SUBROUTINE FLIP_i_k

    SUBROUTINE FLIP_i_k_j_k
        DO j =1,m
            Do k =1,m
                DO i=1,m
                    f ( i , j , k) = f ( i , j , k) * 2d0
                ENDDO
            ENDDO
        ENDDO
    END SUBROUTINE FLIP_i_k_j_k
END PROGRAM timeing_loop

!! --------------------- Results -------------------------- !!
! Struture: compiler : flags (f95=sun compiler,  n_max:  10000000)
!
! 1.: f95 : -O0
! (sub INITIAL) Avg, Wall time =  1.6093159E-7 (tot:  1.6093159 )
! (sub FLIP_i_k) Avg, Wall time =  1.5586004E-7 (tot:  1.5586004 )
! (sub FLIP_i_k_j_k) Avg, Wall time =  1.555914E-7 (tot:  1.5559139 )
!
! 2.: gfortran : -O0
! (sub INITIAL) Avg, Wall time =    3.48219885E-07 (tot:    3.48219872     )
! (sub FLIP_i_k) Avg, Wall time =    3.39278785E-07 (tot:    3.39278793     )
! (sub FLIP_i_k_j_k) Avg, Wall time =    3.38117019E-07 (tot:    3.38117027     )
!
! 3.: f95 : -O3
! (sub INITIAL) Avg, Wall time =  3.9894818E-8 (tot:  0.3989482 )
! (sub FLIP_i_k) Avg, Wall time =  3.978181E-8 (tot:  0.3978181 )
! (sub FLIP_i_k_j_k) Avg, Wall time =  3.9766167E-8 (tot:  0.3976617 )
!
! 4.: gfortran : -O3
! (sub INITIAL) Avg, Wall time =    2.36039355E-07 (tot:    2.36039352     )
! (sub FLIP_i_k) Avg, Wall time =    2.36034779E-07 (tot:    2.36034775     )
! (sub FLIP_i_k_j_k) Avg, Wall time =    2.35964492E-07 (tot:    2.35964489     )
!
! Conclution: The sun compiler is better than the gfortran compiler with a run time
! that is allmost 1/3 the time for the gfortran compiler


