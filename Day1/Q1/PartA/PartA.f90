PROGRAM timeing_loop
    IMPLICIT NONE
    INTEGER, PARAMETER :: MK = KIND(0d0)
    INTEGER(MK) :: n = 1e8
    REAL(MK), DIMENSION(4) :: a
    REAL ts , te , t

    CALL PartA_int

    CALL CPU_TIME( ts )
    CALL PartA_int
    CALL CPU_TIME( te )
    t = te-ts
    WRITE(*,*) "(Part A - Initial) Wall time = ", t

    CALL CPU_TIME( ts )
    CALL PartA_move_loop
    CALL CPU_TIME( te )
    t = te-ts
    WRITE(*,*) "(Part A - Flip loop) Wall time = ", t


    CALL CPU_TIME( ts )
    CALL PartA_unrole
    CALL CPU_TIME( te )
    t = te-ts
    WRITE(*,*) "(Part A - Unrole) Wall time = ", t

    CONTAINS
        SUBROUTINE PartA_int
        INTEGER(MK) :: i,j

        a = 1d0
        DO i = 1,n
            DO j = 1,4
                a(j) = a(j) + a(j) -2d0
            END DO
        END DO

    END SUBROUTINE PartA_int
    SUBROUTINE PartA_move_loop
        INTEGER(MK) :: i,j

        a = 1d0
        DO j = 1,4
            DO i = 1,n
                a(j) = a(j) + a(j) -2d0
            END DO
        END DO

    END SUBROUTINE PartA_move_loop
    SUBROUTINE PartA_unrole
        INTEGER(MK) :: i

        a = 1d0
        DO i = 1,n
            a(1) = a(1) + a(1) -2d0
            a(2) = a(2) + a(2) -2d0
            a(3) = a(3) + a(3) -2d0
            a(4) = a(4) + a(4) -2d0
        END DO

    END SUBROUTINE PartA_unrole

END PROGRAM timeing_loop