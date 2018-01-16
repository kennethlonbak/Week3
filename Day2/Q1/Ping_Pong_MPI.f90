PROGRAM Ping_pong_MPI
    IMPLICIT NONE
    INCLUDE "mpif.h"
    ! Declaring
    INTEGER, PARAMETER :: msg_len = 200
    INTEGER :: rank, Npr, ierror, i, n=100000
    CHARACTER(LEN=msg_len) :: msg_send, msg_recv
    INTEGER :: status(MPI_STATUS_SIZE)
    REAL :: ts, te

    ! Initilizing

    CALL MPI_Init(ierror)
    CALL MPI_Comm_Rank(MPI_COMM_WORLD,rank,ierror)
    CALL MPI_Comm_Size(MPI_COMM_WORLD,Npr,ierror)
    WRITE(*,*) "Rank:", rank, " (Npr:", Npr, ")"
    ts = MPI_WTIME()
    DO i=1,n
        IF (rank == 0) THEN
            ! Send
            msg_send = "Ping"
            CALL MPI_SSend(msg_send,msg_len,MPI_CHARACTER,1,i,MPI_COMM_WORLD,ierror)
            ! Recive
            CALL MPI_Recv(msg_recv,msg_len,MPI_CHARACTER,1,i+n+1,MPI_COMM_WORLD,status,ierror)
            WRITE(*,*) "i=",i,"Rank 0, got ", TRIM(msg_recv), " from rank 1"
        ELSE
            ! Recive
            CALL MPI_Recv(msg_recv,msg_len,MPI_CHARACTER,0,i,MPI_COMM_WORLD,status,ierror)
            WRITE(*,*) "i=",i,"Rank 1, got ", TRIM(msg_recv), " from rank 0"
            ! Send
            msg_send = "Pong"
            CALL MPI_SSend(msg_send,msg_len,MPI_CHARACTER,0,i+n+1,MPI_COMM_WORLD,ierror)
        END IF
    END DO
    CALL MPI_Barrier(MPI_COMM_WORLD,ierror)
    te = MPI_WTIME()

    CALL MPI_Finalize(ierror)
    CALL SLEEP(1)
    IF (rank == 1) WRITE(*,*) "Avg. transfer time (pr. msg):", &
            (te-ts)/n,"pm", (te-ts)/SQRT(REAL(n))**3 ," (Wall tot:", te-ts, ", n:",n, ")"


END PROGRAM Ping_pong_MPI

! Commenting out print statements to run the test
! Output:
! Rank:           0  (Npr:           2 )
! Rank:           1  (Npr:           2 )
! Avg. transfer time (pr. msg):   1.21093751E-06 pm   3.82932130E-09  (Wall tot:  0.121093750     , n:      100000 )
!
! With print statements:
! i=       99999 Rank 0, got Pong from rank 1
! i=      100000 Rank 0, got Pong from rank 1
! Avg. transfer time (pr. msg):   3.26562513E-05 pm   1.03268142E-07  (Wall tot:   3.26562500     , n:      100000 )


