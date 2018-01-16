PROGRAM Hallo_World_MPI
    IMPLICIT NONE
    INCLUDE "mpif.h"
    INTEGER :: rank, Npr, ierror
    CALL MPI_Init(ierror)
    CALL MPI_Comm_Rank(MPI_COMM_WORLD,rank,ierror)
    CALL MPI_Comm_Size(MPI_COMM_WORLD,Npr,ierror)

    IF (rank == 0) THEN
        WRITE(*,*) "Hallo World, from rank:",rank, " (Npr:", Npr, ")"
    ELSE
        WRITE(*,*) "rank:",rank, " (Npr:", Npr, ")"
    END IF

    CALL MPI_Finalize(ierror)
END PROGRAM Hallo_World_MPI

! For: mpirun -np 4 Hallo_World_MPI
! The output is:
! rank:           2  (Npr:           4 )
! rank:           3  (Npr:           4 )
! Hallo World, from rank:           0  (Npr:           4 )
! rank:           1  (Npr:           4 )
!
! Which is the desired output
