PROGRAM Bandwidth_MPI
    IMPLICIT NONE
    INCLUDE "mpif.h"
    ! Declaring
    INTEGER, PARAMETER ::  MK=KIND(1d0)
    INTEGER(MK), PARAMETER :: n_array = 2**18
    INTEGER(MK) :: i, k, imax_num = 1e6, imax_array = 1e4
    INTEGER :: rank, Npr, ierror
    INTEGER :: status(MPI_STATUS_SIZE)
    CHARACTER(LEN=100) :: fileout
    REAL(MK), DIMENSION(n_array) :: array_send=1d0, array_recv
    REAL(8) :: ts, te
    INTEGER :: num_send, num_recv

    CALL MPI_Init(ierror)
    CALL MPI_Comm_Rank(MPI_COMM_WORLD,rank,ierror)
    CALL MPI_Comm_Size(MPI_COMM_WORLD,Npr,ierror)

    IF (MOD(Npr,2) /= 0) THEN
        WRITE(*,*) "Npr needs to be an even number! (Given Npr=",Npr,")"
        CALL MPI_Finalize(ierror)
        STOP
    END IF

    ! ------------------------ Latency -------------------------------- !
    ts  = MPI_WTIME()
    DO i=1,imax_num
        IF (MOD(rank,2) == 0) THEN
            ! Send
            !WRITE(*,*) "SEND", rank+1,rank
            CALL MPI_SSend(num_send,1,MPI_INTEGER,rank+1,rank,MPI_COMM_WORLD,ierror)
            ! Recive
            CALL MPI_Recv(num_recv,1,MPI_INTEGER,rank+1,rank+1,MPI_COMM_WORLD,status,ierror)
        ELSE
            ! Recive
            !WRITE(*,*) "RECV", rank-1,rank-1
            CALL MPI_Recv(num_recv,1,MPI_INTEGER,rank-1,rank-1,MPI_COMM_WORLD,status,ierror)
            ! Send
            CALL MPI_SSend(num_send,1,MPI_INTEGER,rank-1,rank,MPI_COMM_WORLD,ierror)
        END IF
    END DO
    te = MPI_WTIME()
    IF (MOD(rank,2) == 0) THEN
        WRITE(*,"(A,E10.4,A,I1,AI1,A,F15.7,A)") "Avg. Latency is: ", (te-ts)/imax_num, &
                " s (Rank ", rank, "<->Rank ",rank+1,", t_tot:", te-ts, ")"
    END IF

    ! ------------------- Bandwidth -------------------------------------- !
    ts  = MPI_WTIME()
    DO i=1,imax_array
        IF (MOD(rank,2) == 0) THEN
            ! Send
            CALL MPI_SSend(array_send,n_array,MPI_DOUBLE_PRECISION,rank+1,rank,MPI_COMM_WORLD,ierror)
            ! Recive
            CALL MPI_Recv(array_recv,n_array,MPI_DOUBLE_PRECISION,rank+1,rank+1,MPI_COMM_WORLD,status,ierror)
        ELSE
            ! Recive
            CALL MPI_Recv(array_recv,n_array,MPI_DOUBLE_PRECISION,rank-1,rank-1,MPI_COMM_WORLD,status,ierror)
            ! Send
            CALL MPI_SSend(array_send,n_array,MPI_DOUBLE_PRECISION,rank-1,rank,MPI_COMM_WORLD,ierror)
        END IF
    END DO
    te = MPI_WTIME()
    IF (MOD(rank,2) == 0) THEN
        WRITE(*,"(A,f8.1,A,I1,AI1,A,F15.7,A)") "Avg. Bandwidth is: ", real(8*n_array*imax_array)/(1d6*(te-ts)), &
                " MB/s (Rank ", rank, "<->Rank ",rank+1,", t_tot:", te-ts, ")"
    END IF


    CALL MPI_Finalize(ierror)


END PROGRAM Bandwidth_MPI

! Running:
!
! mpirun -np 2 Bandwidth_MPI.prog
! Avg. Latency is: 0.1061E-05 s (Rank 0<->Rank 1, t_tot:      1.0606183)
! Avg. Bandwidth is:   7326.6 MB/s (Rank 0<->Rank 1, t_tot:      2.8623627)
!
! mpirun -np 4 Bandwidth_MPI.prog
! Avg. Latency is: 0.1082E-05 s (Rank 2<->Rank 3, t_tot:      1.0818536)
! Avg. Latency is: 0.1092E-05 s (Rank 0<->Rank 1, t_tot:      1.0919312)
! Avg. Bandwidth is:   4733.0 MB/s (Rank 0<->Rank 1, t_tot:      4.4308716)
! Avg. Bandwidth is:   4702.3 MB/s (Rank 2<->Rank 3, t_tot:      4.4598312)
!
! mpirun -np 6 Bandwidth_MPI.prog
! Avg. Latency is: 0.1293E-05 s (Rank 0<->Rank 1, t_tot:      1.2932095)
! Avg. Latency is: 0.1333E-05 s (Rank 2<->Rank 3, t_tot:      1.3325731)
! Avg. Latency is: 0.1333E-05 s (Rank 4<->Rank 5, t_tot:      1.3330252)
! Avg. Bandwidth is:   2980.2 MB/s (Rank 0<->Rank 1, t_tot:      7.0368865)
! Avg. Bandwidth is:   2885.8 MB/s (Rank 2<->Rank 3, t_tot:      7.2672328)
! Avg. Bandwidth is:   2875.8 MB/s (Rank 4<->Rank 5, t_tot:      7.2924287)
!
! mpirun -np 8 Bandwidth_MPI.prog
! Avg. Latency is: 0.1512E-05 s (Rank 2<->Rank 3, t_tot:      1.5117851)
! Avg. Latency is: 0.1552E-05 s (Rank 6<->Rank 7, t_tot:      1.5515365)
! Avg. Latency is: 0.1573E-05 s (Rank 0<->Rank 1, t_tot:      1.5732317)
! Avg. Latency is: 0.1716E-05 s (Rank 4<->Rank 5, t_tot:      1.7155836)
! Avg. Bandwidth is:   2188.3 MB/s (Rank 0<->Rank 1, t_tot:      9.5833186)
! Avg. Bandwidth is:   2169.1 MB/s (Rank 6<->Rank 7, t_tot:      9.6681903)
! Avg. Bandwidth is:   2101.7 MB/s (Rank 4<->Rank 5, t_tot:      9.9785839)
! Avg. Bandwidth is:   2057.2 MB/s (Rank 2<->Rank 3, t_tot:     10.1939995)





