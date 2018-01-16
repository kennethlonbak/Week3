PROGRAM Bandwidth_MPI
    IMPLICIT NONE
    INCLUDE "mpif.h"
    ! Declaring
    INTEGER, PARAMETER ::  MK=KIND(1d0)
    INTEGER(MK), PARAMETER :: msg_len = 200, n=2**18
    INTEGER(MK) :: i
    INTEGER :: rank, Npr, ierror
    INTEGER :: status(MPI_STATUS_SIZE)
    REAL(MK), DIMENSION(n) :: array_send, array_recv

    ! Initilizing

    CALL MPI_Init(ierror)
    CALL MPI_Comm_Rank(MPI_COMM_WORLD,rank,ierror)
    CALL MPI_Comm_Size(MPI_COMM_WORLD,Npr,ierror)
    WRITE(*,*) "Rank:", rank, " (Npr:", Npr, ")"

    DO i=1,n
        IF (rank == 0) THEN
            ! Send
            array_send(i) = MPI_WTIME()
            CALL MPI_SSend(array_send(1:i),i,MPI_DOUBLE_PRECISION,1,i,MPI_COMM_WORLD,ierror)
            ! Recive
            CALL MPI_Recv(array_recv(1:i),i,MPI_DOUBLE_PRECISION,1,i-1,MPI_COMM_WORLD,status,ierror)
            array_send(i) = MPI_WTIME()-array_send(i)
            !WRITE(*,*) "(Rank 0) Last element recv:", array_recv(i)
        ELSE
            ! Recive
            CALL MPI_Recv(array_recv(1:i),i,MPI_DOUBLE_PRECISION,0,i,MPI_COMM_WORLD,status,ierror)
            !WRITE(*,*) "(Rank 1) Last element recv:", array_recv(i)
            ! Send
            CALL MPI_SSend(array_send(1:i),i,MPI_DOUBLE_PRECISION,0,i-1,MPI_COMM_WORLD,ierror)
        END IF
    END DO
    CALL MPI_Barrier(MPI_COMM_WORLD,ierror)

    IF (rank == 0) THEN
        OPEN(61,FILE="TIMING.dat",action="write")
        DO i = 1,n
            WRITE(61,*) array_send(i), i*8
        END DO
    END IF

    CALL MPI_Finalize(ierror)
    CALL SLEEP(1)




END PROGRAM Bandwidth_MPI


