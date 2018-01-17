   PROGRAM Non_blocking_communication
   IMPLICIT NONE 
   INCLUDE 'mpif.h'
   INTEGER :: rank, Nproc, ierror, tag, dest, next, prev, iproc
   INTEGER :: count=1, status(MPI_STATUS_SIZE)
   INTEGER :: sendbuf, recvbuf, handle
   DOUBLE PRECISION ::  ts , te

   CALL MPI_INIT(ierror)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, Nproc, ierror)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

!Sequence
   prev=MOD(rank-1,Nproc)
   next=MOD(rank+1,Nproc) 
IF (rank == 0) THEN !This is done so rank 0 does not receive from -1, but from the last rank instead. 
prev=Nproc-1
ENDIF
!PRINT*, rank, 'prev', prev, 'next', next

!Loop
sendbuf=rank !Initialize
DO iproc=0,Nproc-1 
   IF (iproc.LT.Nproc-1) THEN! For all iteratinos except the last, a send-receive is called
   CALL MPI_SENDRECV(sendbuf,1,MPI_INTEGER, next,tag,recvbuf,1,MPI_INTEGER,prev,tag,MPI_COMM_WORLD,status,ierror)
   sendbuf=recvbuf+rank
   ENDIF

   !CALL SLEEP(1)
   !Print*, 'Iteration: ', iproc, 'Processor no:', rank,',Receive buffer' recvbuf, ',Send buffer: ', sendbuf,
ENDDO
CALL SLEEP (1) 
PRINT*, 'Processor no:', rank, 'SUM of ranks is = ', sendbuf
CALL MPI_FINALIZE(ierror)
END PROGRAM


