MODULE m_Par_matmul

    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTERFACE read_vari_arg
        MODULE PROCEDURE read_vari_arg_integer, read_vari_arg_double, read_vari_arg_char, read_vari_arg_logical
    END INTERFACE read_vari_arg

    INTEGER, PARAMETER :: MK = KIND(0d0)

    ! Matrix variabels
    REAL(MK), DIMENSION(:,:), ALLOCATABLE  :: C, A_own, A_recv, B_own, B_recv
    REAL(MK), DIMENSION(:,:), ALLOCATABLE  :: C_tot, A_tot, B_tot

    REAL(MK), DIMENSION(:), ALLOCATABLE  :: sr_vec
    INTEGER :: N, i, j, k, k_bar, is, ie, js, je, isize, jsize, q

    ! MPI varaiabels -------------------------------
    INTEGER, PARAMETER :: ndim = 2
    INTEGER, DIMENSION(4) :: src, dest, handle, tag
    LOGICAL, DIMENSION(ndim) :: wrap = .true.
    INTEGER, DIMENSION(ndim) :: dimsize, coord
    ! Comunicators
    INTEGER :: cart_comm, row_comm
    ! Default arguments
    INTEGER :: ierror, Nproc, rank, rank_row, rank_col, status(MPI_STATUS_SIZE), root_rank

    CONTAINS
    SUBROUTINE initialize
        N = read_vari_arg("N",99)

        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, Nproc, ierror) ! Getting Nproc

        IF (MOD(INT(SQRT(REAL(Nproc))*100),100)==0) THEN
            dimsize = int(SQRT(REAL(Nproc)))
            q = dimsize(1)
            IF (MOD(INT(N/REAL(q)*100),100)/=0) THEN
                PRINT*, "Only work for N/q as a whole number (q=",q," N=",N,INT(REAL(N)/q)*100,")"
                CALL MPI_FINALIZE(ierror)
                STOP
            END IF
        ELSE
            PRINT*, "Only work for sqrt(Nproc) as a whole number"
            CALL MPI_FINALIZE(ierror)
            STOP
        END IF


    END SUBROUTINE initialize

    SUBROUTINE create_cart

        CALL MPI_Cart_create(MPI_COMM_WORLD,ndim,dimsize,wrap,.true.,cart_comm,ierror)
        CALL MPI_COMM_RANK(cart_comm, rank, ierror)
        CALL MPI_Cart_coords(cart_comm,rank,ndim,coord,ierror)

    END SUBROUTINE create_cart

    SUBROUTINE check_cart_grid
        INTEGER, DIMENSION(ndim) :: coord_temp
        INTEGER :: ist,iet,jst,jet,isizet,jsizet
        IF (rank == 0) THEN
            DO i = 0,Nproc-1
                CALL MPI_Cart_coords(cart_comm,i,ndim,coord_temp,ierror)
                CALL get_index_from_coord(coord_temp,ist,iet,jst,jet,isizet,jsizet)
                WRITE(*,"(A,I2,A,I2,A,I2)") "Rank", i, " coord(1)=",coord_temp(1), " coord(2)=", coord_temp(2)
                WRITE(*,"(A,I4,A,I4,A,I4,A,I4)") "is",ist,", ie",iet,", js",jst,", je",jet
                WRITE(*,"(A,I4,A,I4)") "isize=", isizet, ", jsize=", jsizet
                WRITE(*,*)
            END DO
        END IF
    END SUBROUTINE check_cart_grid

    SUBROUTINE get_index_from_coord(coord_in,is_out,ie_out,js_out,je_out,isize,jsize)
        IMPLICIT NONE
        INTEGER, DIMENSION(ndim), INTENT(IN) :: coord_in
        INTEGER, INTENT(OUT) :: is_out,ie_out,js_out,je_out,isize,jsize
        is_out = FLOOR(coord_in(1)*N/REAL(dimsize(1)))+1
        ie_out = FLOOR((coord_in(1)+1)*N/REAL(dimsize(1)))
        js_out = FLOOR(coord_in(2)*N/REAL(dimsize(2)))+1
        je_out = FLOOR((coord_in(2)+1)*N/REAL(dimsize(2)))
        isize = ie_out-is_out+1
        jsize = je_out-js_out+1
    END SUBROUTINE get_index_from_coord

    SUBROUTINE gen_sub_matrix_rank
        C = 0d0
        DO i = 1,isize
            DO j = 1,jsize
                A_own(i,j) = rank
                B_own(i,j) = rank
                A_recv(i,j) = rank
                B_recv(i,j) = rank
            END DO
        END DO
    END SUBROUTINE gen_sub_matrix_rank

    SUBROUTINE gen_sub_matrix_cos
        C = 0d0
        DO i = 1,isize
            DO j = 1,jsize
                A_own(i,j) = cos(REAL((i-1)*N+j-1))
                B_own(i,j) = cos(REAL((i-1)*N+j-1))
                A_recv(i,j) = cos(REAL((i-1)*N+j-1))
                B_recv(i,j) = cos(REAL((i-1)*N+j-1))
            END DO
        END DO
    END SUBROUTINE gen_sub_matrix_cos

    SUBROUTINE broadcast_row_mat
        root_rank = k_bar
        ! Send matrix size
        DO j= 1,jsize
            sr_vec(1:isize) = A_own(1:isize,j)
            CALL MPI_Bcast(sr_vec,isize,MPI_DOUBLE_PRECISION,root_rank,row_comm,ierror)
            A_recv(1:isize,j) = sr_vec(1:isize)
        END DO
        !CALL MPI_Barrier(cart_comm,ierror)
    END SUBROUTINE broadcast_row_mat

    SUBROUTINE test_broadcast(coord_in, pre_string)
        INTEGER, DIMENSION(ndim), INTENT(IN) :: coord_in
        CHARACTER(LEN=*), INTENT(IN) :: pre_string

        IF (coord(1) == coord_in(1).and.coord(2) == coord_in(2)) THEN
            WRITE(*,*) TRIM(pre_string)," root_rank",root_rank," Coord=(",coord,")",", k_bar=",k_bar
            WRITE(*,*) "A_recv(1,1)=",A_recv(1,1),"A_own(1,1)=",A_own(1,1)
            WRITE(*,*)
        END IF
    END SUBROUTINE

    SUBROUTINE send_B
        INTEGER :: temp, rank_above, rank_below
        ! Get rank of above submatrix (to send)
        CALL MPI_Cart_shift(cart_comm,0,-1,temp,rank_above,ierror)
        CALL MPI_Cart_shift(cart_comm,0, 1,temp,rank_below,ierror)

        IF (rank==0) THEN
            !CALL MPI_Cart_coords(cart_comm,rank,ndim,coords,ierror)
        END IF


        ! Reciving mat from below
        IF (rank_col > 0) THEN
            CALL recv_mat(B_recv, cart_comm, rank_below)
        END IF

        ! Sending mat
        CALL send_mat(B_own, isize, jsize, cart_comm, rank_above)

        IF (rank_col == 0) THEN
            CALL recv_mat(B_recv, cart_comm, rank_below)
        END IF

        B_own = B_recv

    END SUBROUTINE send_B

    SUBROUTINE send_mat(mat, isize_in, jsize_in, comm, send_to_rank, is_in, ie_in, js_in, je_in)
        IMPLICIT NONE
        REAL(MK), DIMENSION(:,:), INTENT(OUT) :: mat
        INTEGER, INTENT(IN) :: comm, isize_in, jsize_in, send_to_rank
        INTEGER, INTENT(IN), OPTIONAL :: is_in, ie_in, js_in, je_in
        INTEGER :: ierror, j
        REAL(MK), DIMENSION(:), ALLOCATABLE  :: sr_vec

        CALL MPI_SSend(isize_in,1,MPI_INTEGER,send_to_rank,1,comm,ierror)
        CALL MPI_SSend(jsize_in,1,MPI_INTEGER,send_to_rank,1,comm,ierror)
        ALLOCATE(sr_vec(isize_in))

        IF (PRESENT(je_in)) THEN
            CALL MPI_SSend(is_in,1,MPI_INTEGER,send_to_rank,1,comm,ierror)
            CALL MPI_SSend(ie_in,1,MPI_INTEGER,send_to_rank,1,comm,ierror)
            CALL MPI_SSend(js_in,1,MPI_INTEGER,send_to_rank,1,comm,ierror)
            CALL MPI_SSend(je_in,1,MPI_INTEGER,send_to_rank,1,comm,ierror)
        END IF

        DO j=1,jsize_in
            sr_vec = mat(1:isize_in,j)
            CALL MPI_SSend(sr_vec,isize_in,MPI_DOUBLE_PRECISION,send_to_rank,1,comm,ierror)
        END DO
        DEALLOCATE(sr_vec)
    END SUBROUTINE send_mat

    SUBROUTINE recv_mat(mat, comm, recv_from_rank, is_out, ie_out, js_out, je_out)
        IMPLICIT NONE
        REAL(MK), DIMENSION(:,:), INTENT(OUT) :: mat
        INTEGER, INTENT(IN) :: comm, recv_from_rank
        INTEGER, INTENT(OUT), OPTIONAL :: is_out, ie_out, js_out, je_out
        INTEGER :: rank_int, ierror, isizet,jsizet, status(MPI_STATUS_SIZE), i,j
        REAL(MK), DIMENSION(:), ALLOCATABLE  :: sr_vec

        ! Get isize and jsize from recive process
        CALL MPI_Recv(isizet,1,MPI_INTEGER,recv_from_rank,1,comm,status,ierror)
        CALL MPI_Recv(jsizet,1,MPI_INTEGER,recv_from_rank,1,comm,status,ierror)
        ALLOCATE(sr_vec(isizet))

        IF (PRESENT(je_out)) THEN
            CALL MPI_Recv(is_out,1,MPI_INTEGER,recv_from_rank,1,comm,status,ierror)
            CALL MPI_Recv(ie_out,1,MPI_INTEGER,recv_from_rank,1,comm,status,ierror)
            CALL MPI_Recv(js_out,1,MPI_INTEGER,recv_from_rank,1,comm,status,ierror)
            CALL MPI_Recv(je_out,1,MPI_INTEGER,recv_from_rank,1,comm,status,ierror)
        END IF
        DO j=1,jsizet
            CALL MPI_Recv(sr_vec,isizet,MPI_DOUBLE_PRECISION,recv_from_rank,1,comm,status,ierror)
            mat(1:isizet,j) = sr_vec(1:isizet)
            !print*, j
        END DO
        DEALLOCATE(sr_vec)
    END SUBROUTINE recv_mat

    SUBROUTINE test_b_send(coord_in, pre_string)
        INTEGER, DIMENSION(ndim), INTENT(IN) :: coord_in
        CHARACTER(LEN=*), INTENT(IN) :: pre_string

        IF (coord(1) == coord_in(1).and.coord(2) == coord_in(2)) THEN
            WRITE(*,*) TRIM(pre_string)," Coord=(",coord,")",", k_bar=",k_bar,"B_recv(10,10)=",B_recv(10,10)
        END IF
    END SUBROUTINE test_b_send

    SUBROUTINE collect_C_tot

        ! Allocate C_tot
        ALLOCATE(C_tot(N,N))

        ! Set own part of C_tot
        C_tot(is:ie,js:je) = C

        DO i = 1,Nproc-1
            CALL recv_mat(C, cart_comm, i, is, ie, js, je)
            C_tot(is:ie,js:je) = C
        END DO

    END SUBROUTINE

    SUBROUTINE compare_with_seq
        REAL(MK), DIMENSION(N,N)  :: C_seq
        INTEGER :: is_the_same

        CALL write_matrix(A_tot, "A_tot.dat", .true.)
        CALL write_matrix(B_tot, "B_tot.dat", .true.)
        C_seq = matmul(A_tot,B_tot)
        CALL write_matrix(C_seq, "C_seq.dat", .true.)
        is_the_same = 1
        DO j = 1,N
            DO i = 1,N
                !PRINT*,i,j
                !PRINT*, C_tot(i,j), C_seq(i,j)
                IF (C_tot(i,j) /= C_seq(i,j)) THEN
                    IF (abs(C_tot(i,j)- C_seq(i,j)) < 1d-8) THEN
                        is_the_same = 2
                    ELSE
                        is_the_same = 0
                        !exit
                    END IF
                END IF
            END DO
        END DO
        IF (is_the_same == 1) THEN
            PRINT*, "The parallel and sequential gives exactly the same"
        ELSE IF (is_the_same == 2) THEN
            PRINT*, "The parallel and sequential gives the same to machine presision"
        ELSE
            PRINT*, "The parallel and sequential does NOT give the same"
        END IF

    END SUBROUTINE compare_with_seq

    SUBROUTINE collect_A_tot
        INTEGER :: is_tmp,ie_tmp,js_tmp,je_tmp
        REAL(MK), DIMENSION(isize,jsize) :: A_tmp
        IF (rank==0) THEN
            ! Allocate C_tot
            ALLOCATE(A_tot(N,N))


            ! Set own part of C_tot
            A_tot(is:ie,js:je) = A_own

            DO i = 1,Nproc-1
                CALL recv_mat(A_tmp, cart_comm,i, is_tmp,ie_tmp,js_tmp,je_tmp)
                A_tot(is_tmp:ie_tmp,js_tmp:je_tmp) = A_tmp
            END DO

        ELSE
            CALL send_mat(A_own, isize, jsize, cart_comm, 0, is, ie, js, je)
        END IF
    END SUBROUTINE

    SUBROUTINE collect_B_tot
        INTEGER :: is_tmp,ie_tmp,js_tmp,je_tmp
        REAL(MK), DIMENSION(isize,jsize) :: B_tmp
        IF (rank==0) THEN
            ! Allocate C_tot
            ALLOCATE(B_tot(N,N))


            ! Set own part of C_tot
            B_tot(is:ie,js:je) = B_own

            DO i = 1,Nproc-1
                CALL recv_mat(B_tmp, cart_comm,i, is_tmp,ie_tmp,js_tmp,je_tmp)
                B_tot(is_tmp:ie_tmp,js_tmp:je_tmp) = B_tmp
            END DO

        ELSE
            CALL send_mat(B_own, isize, jsize, cart_comm, 0, is, ie, js, je)
        END IF
    END SUBROUTINE

    ! Matrix to file
    SUBROUTINE write_matrix(mat, filename, write_mat)
        REAL(MK), DIMENSION(:,:), INTENT(IN) :: mat
        CHARACTER(LEN=*), INTENT(IN) :: filename
        LOGICAL :: write_mat
        INTEGER :: i,j

        OPEN(61, FILE=filename, action="write")


        IF (write_mat) THEN
            ! writing out matrix
            DO j = 1,N
                DO i = 1,N
                    WRITE(61,"(E16.8E2,A)",advance="no") mat(i,j), " "
                END DO
                WRITE(61,"(A)")
            END DO
        ELSE
            WRITE(61,*) 0D0, 0D0
        END IF
        CLOSE(61)
    END SUBROUTINE write_matrix

    ! Input reading
    FUNCTION read_vari_arg_integer(var_name,default_val)
        INTEGER :: read_vari_arg_integer
        CHARACTER(LEN=*), INTENT(IN) :: var_name
        INTEGER, INTENT(IN) :: default_val
        INTEGER :: i, i_arg
        CHARACTER(LEN=200) :: arg
        read_vari_arg_integer = default_val
        DO i = 1,iargc()
            CALL GETARG(i,arg)
            i_arg = INDEX(arg,var_name)
            IF (i_arg > 0) THEN
                i_arg = INDEX(arg,"=")
                IF (LEN(var_name) == i_arg-1) THEN
                    arg = arg(i_arg+1:len(arg))
                    READ(arg,*) read_vari_arg_integer
                    WRITE(*,*) "The variable ", var_name, " was set,", var_name,"=",read_vari_arg_integer
                END IF
            end if
        end do
    END FUNCTION read_vari_arg_integer

    FUNCTION read_vari_arg_double(var_name,default_val)
        REAL(KIND(0D0)) :: read_vari_arg_double
        CHARACTER(LEN=*), INTENT(IN) :: var_name
        REAL(KIND(0D0)), INTENT(IN) :: default_val
        INTEGER :: i, i_arg
        CHARACTER(LEN=200) :: arg
        read_vari_arg_double = default_val
        DO i = 1,iargc()
            CALL GETARG(i,arg)
            i_arg = INDEX(arg,var_name)
            IF (i_arg > 0) THEN
                i_arg = INDEX(arg,"=")
                IF (LEN(var_name) == i_arg-1) THEN
                    arg = arg(i_arg+1:len(arg))
                    READ(arg,*) read_vari_arg_double
                    WRITE(*,*) "The variable ", var_name, " was set,", var_name,"=",read_vari_arg_double
                END IF
            end if
        end do
    END FUNCTION read_vari_arg_double

    FUNCTION read_vari_arg_char(var_name,default_val,loc_print)
        CHARACTER(LEN=200) :: read_vari_arg_char
        CHARACTER(LEN=*), INTENT(IN) :: var_name
        CHARACTER(LEN=*), INTENT(IN) :: default_val
        LOGICAL, OPTIONAL :: loc_print
        INTEGER :: i, i_arg
        CHARACTER(LEN=200) :: arg
        IF (.NOT.PRESENT(loc_print)) loc_print = .true.
        read_vari_arg_char = default_val
        DO i = 1,iargc()
            CALL GETARG(i,arg)
            i_arg = INDEX(arg,var_name)
            IF (i_arg > 0) THEN
                i_arg = INDEX(arg,"=")
                IF (LEN(var_name) == i_arg-1) THEN
                    arg = arg(i_arg+1:len(arg))
                    read_vari_arg_char = arg
                    IF (loc_print) THEN
                        WRITE(*,*) "The variable ", var_name, " was set,", var_name,"=",TRIM(read_vari_arg_char)
                    END IF
                END IF
            end if
        end do
    END FUNCTION read_vari_arg_char

    FUNCTION read_vari_arg_logical(var_name,default_val)
        LOGICAL :: read_vari_arg_logical
        CHARACTER(LEN=*), INTENT(IN) :: var_name
        LOGICAL, INTENT(IN) :: default_val
        INTEGER :: i, i_arg
        CHARACTER(LEN=200) :: arg
        read_vari_arg_logical = default_val
        DO i = 1,iargc()
            CALL GETARG(i,arg)
            i_arg = INDEX(arg,var_name)
            IF (i_arg > 0) THEN
                i_arg = INDEX(arg,"=")
                IF (LEN(var_name) == i_arg-1) THEN
                    arg = TRIM(arg(i_arg+1:len(arg)))
                    IF (arg == ".false.") THEN
                        read_vari_arg_logical = .false.
                    ELSE
                        read_vari_arg_logical = .true.
                    END IF
                    WRITE(*,*) "The variable ", var_name, " was set,", var_name,"=",read_vari_arg_logical
                END IF
            end if
        end do
    END FUNCTION read_vari_arg_logical

end module m_Par_matmul