MODULE m_GS_solver

    INCLUDE 'mpif.h'
    INTERFACE read_vari_arg
        MODULE PROCEDURE read_vari_arg_integer, read_vari_arg_double, read_vari_arg_char, read_vari_arg_logical
    END INTERFACE read_vari_arg

    INTEGER, PARAMETER :: MK = KIND(0d0)
    ! Variabels that can be changed by user
    INTEGER :: N, k_max, k, N_th
    REAL(MK) :: d_min, d_loc, d
    CHARACTER(LEN=200) :: filename
    ! Variabels to be used by the solver
    REAL(MK), DIMENSION(:,:), ALLOCATABLE :: uk, fdx2
    REAL(MK), DIMENSION(:), ALLOCATABLE :: bus,bds,bls,brs
    REAL(MK), DIMENSION(:), ALLOCATABLE :: bur,bdr,blr,brr
    REAL(MK) :: dx, dx2
    REAL(8) :: wall_time
    INTEGER :: solver_type ! Solver type: 1= Jacobi, 2=Gauss-Seidel
    INTEGER :: mod_state, problem
    LOGICAL :: show_state, write_mat

    ! MPI variabels
    INTEGER, PARAMETER :: ndim = 2
    INTEGER :: rank, Nproc, ierror, CART_COMM, i_dir,  i_dest
    INTEGER, DIMENSION(4) :: src, dest, handle, tag
    LOGICAL, DIMENSION(ndim) :: wrap = .false.
    INTEGER, DIMENSION(ndim) :: dimsize, coords, dims, coor_temp
    INTEGER :: status(MPI_STATUS_SIZE), j_u, j_d, i_l, i_r
    INTEGER :: i, j, rb_i, i_min, i_max, j_min, j_max, i_size, j_size, is,ie,js,je

    CONTAINS
    SUBROUTINE RB_GS_PAL
        IMPLICIT NONE
        ! Initialize MPI
        CALL MPI_INIT(ierror)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, Nproc, ierror)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

        ! Making grid layout
        SELECT CASE (Nproc)
            CASE (1)
                PRINT *, "It is not possible to use Nproc=1. Use solver_type=2 or 3 for that"
                CALL MPI_FINALIZE(ierror)
                STOP
            CASE (2)
                dimsize(1) = 2
                dimsize(2) = 1
            CASE (4)
                dimsize(1) = 2
                dimsize(2) = 2
            CASE (6)
                dimsize(1) = 3
                dimsize(2) = 2
            CASE (8)
                dimsize(1) = 4
                dimsize(2) = 2
            CASE (9)
                dimsize(1) = 3
                dimsize(2) = 3
            CASE DEFAULT
                PRINT *, "Only Nproc=2,4 is supported"
                CALL MPI_FINALIZE(ierror)
                STOP
        END SELECT

        ! Creating cart communicator
        CALL MPI_Cart_create(MPI_COMM_WORLD,ndim,dimsize,wrap,.true.,CART_COMM,ierror)

        ! Getting the coordinates for each process
        CALL MPI_Cart_get(cart_comm,ndim,dimsize,wrap,coords,ierror)
        CALL MPI_Cart_rank(cart_comm,coords,rank,ierror)
        !PRINT*,"Coords:",coords, " dimsize:", dimsize

        ! Make sub grid limits from coords
        i_min = FLOOR(coords(2)*N/REAL(dimsize(2)))+1
        i_max = FLOOR((coords(2)+1)*N/REAL(dimsize(2)))
        j_min = FLOOR(coords(1)*N/REAL(dimsize(1)))+1
        j_max = FLOOR((coords(1)+1)*N/REAL(dimsize(1)))
        i_size = i_max-i_min+1
        j_size = j_max-j_min+1

        ! Getting naibour location
        CALL MPI_Cart_shift(cart_comm, 0, 1, src(1), dest(1), ierror) ! Up
        CALL MPI_Cart_shift(cart_comm, 0,-1, src(2), dest(2), ierror) ! Down
        CALL MPI_Cart_shift(cart_comm, 1,-1, src(3), dest(3), ierror) ! Left
        CALL MPI_Cart_shift(cart_comm, 1, 1, src(4), dest(4), ierror) ! Right

        ! Checking if process is at global bounday
        j_u = logical2integer((dest(1) == MPI_PROC_NULL))
        j_d = logical2integer((dest(2) == MPI_PROC_NULL))
        i_l = logical2integer((dest(3) == MPI_PROC_NULL))
        i_r = logical2integer((dest(4) == MPI_PROC_NULL))

        ! Making run index
        is = i_min+i_l
        ie = i_max-i_r
        js = j_min+j_d
        je = j_max-j_u

        ! Allocating bond array
        ALLOCATE(bus(i_size),bds(i_size),bls(j_size),brs(j_size))
        ALLOCATE(bur(i_size),bdr(i_size),blr(j_size),brr(j_size))

        CALL PRINT_RANK_AND_NAI(1)
        DO k = 1,k_max

            DO j= js,je
                DO i = is,ie
                    IF (MOD(i+j,2)==0) THEN
                        uk(i,j) = (uk(i,j-1)+uk(i,j+1)+uk(i-1,j)+uk(i+1,j)+fdx2(i,j))*25d-2
                    END IF
                END DO
            END DO

            ! Update Ghost points
            CALL UPDATE_GHOST

            ! BLACK grid ---------------------------------------------------- !
            ! Send and recive bond on BLACK grid
            DO j= js,je
                DO i = is,ie
                    IF (MOD(i+j,2)/=0) THEN
                        uk(i,j) = (uk(i,j-1)+uk(i,j+1)+uk(i-1,j)+uk(i+1,j)+fdx2(i,j))*25d-2
                    END IF
                END DO
            END DO

            ! Update Ghost points
            CALL UPDATE_GHOST

            ! Calculating residual
            d_loc = get_residual()
            CALL MPI_Allreduce(d_loc, d, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cart_comm, ierror)

            ! Build convergence cretia
            IF (d < d_min) THEN
                !print*,d_loc
                exit
            end if

        end do

        IF (rank==0) THEN
            IF (k > k_max) THEN
                !WRITE(*,*)
                WRITE(*,*) "!! WARNING: The solver did NOT converge within k_max:", k_max, " !!"
                WRITE(*,"(A,ES8.2E2,A,ES8.2E2,A)") " d: ", d, " (d_min:",d_min,")"
                !WRITE(*,*)
            ELSE
                WRITE(*,*) "The solver converged after: ", k, "Iterations (k_max: ", k_max, ")"
                WRITE(*,"(A,ES8.2E2,A,ES8.2E2,A)") " d: ", d, " (d_min:",d_min,")"
            end if
        END IF

        ! Collecting matrix
        IF (rank == 0) THEN
            DO i = 1,Nproc-1
                ! Get i,j range (i_min, i_max, j_min, j_max)
                CALL MPI_Recv(i_min,1,MPI_INTEGER,i,1,cart_comm,status,ierror)
                CALL MPI_Recv(i_max,1,MPI_INTEGER,i,1,cart_comm,status,ierror)
                CALL MPI_Recv(j_min,1,MPI_INTEGER,i,1,cart_comm,status,ierror)
                CALL MPI_Recv(j_max,1,MPI_INTEGER,i,1,cart_comm,status,ierror)
                i_size = i_max-i_min+1

                ! Send each vector in
                DO j = j_min,j_max
                    CALL MPI_Recv(bls, i_size,MPI_DOUBLE_PRECISION,i,1,cart_comm,status,ierror)
                    uk(i_min:i_max,j)=bls
                END DO
            END DO
        ELSE
            ! Send i,j range (i_min, i_max, j_min, j_max)
            CALL MPI_Ssend(i_min, 1,MPI_INTEGER, 0, 1 ,cart_comm, ierror)
            CALL MPI_Ssend(i_max, 1,MPI_INTEGER, 0, 1 ,cart_comm, ierror)
            CALL MPI_Ssend(j_min, 1,MPI_INTEGER, 0, 1 ,cart_comm, ierror)
            CALL MPI_Ssend(j_max, 1,MPI_INTEGER, 0, 1 ,cart_comm, ierror)

            ! Send rows in matrix
            DO j = j_min,j_max
                bls = uk(i_min:i_max,j)
                CALL MPI_Ssend(bls, i_size,MPI_DOUBLE_PRECISION, 0, 1 ,cart_comm, ierror)
            END DO
        END IF
        CALL MPI_FINALIZE(ierror)
    END SUBROUTINE RB_GS_PAL

    SUBROUTINE UPDATE_GHOST
            ! SEND -------------------------------------------- !
            ! Up
            bus(1:i_size) = uk(i_min:i_max,j_max)
            CALL MPI_ISsend(bus(1:i_size), i_size,MPI_DOUBLE_PRECISION, dest(1),1 ,cart_comm, handle(1), ierror)
            ! Down
            bds(1:i_size) = uk(i_min:i_max,j_min)
            CALL MPI_ISsend(bds(1:i_size), i_size,MPI_DOUBLE_PRECISION, dest(2),1 ,cart_comm, handle(2), ierror)
            ! Left
            bls(1:j_size) = uk(i_min,j_min:j_max)
            CALL MPI_ISsend(bls(1:j_size), j_size,MPI_DOUBLE_PRECISION, dest(3),1 ,cart_comm, handle(3), ierror)
            ! Right
            brs(1:j_size) = uk(i_max,j_min:j_max)
            CALL MPI_ISsend(brs(1:j_size), j_size,MPI_DOUBLE_PRECISION, dest(4),1 ,cart_comm, handle(4), ierror)

            ! Recive ------------------------------------------ !
            ! Up
            CALL MPI_IRECV(bur, j_size,MPI_DOUBLE_PRECISION, dest(1), 1, cart_comm, status, ierror)
            ! Down
            CALL MPI_IRECV(bdr, j_size,MPI_DOUBLE_PRECISION, dest(2), 1, cart_comm, status, ierror)
            ! Left
            CALL MPI_IRECV(blr, i_size,MPI_DOUBLE_PRECISION, dest(3), 1, cart_comm, status, ierror)
            ! Down
            CALL MPI_IRECV(brr, i_size,MPI_DOUBLE_PRECISION, dest(4), 1, cart_comm, status, ierror)

            ! Wait ------------------------------------------- !
            ! Up
            CALL MPI_Wait(handle(1),status,ierror)
            ! Down
            CALL MPI_Wait(handle(2),status,ierror)
            ! Left
            CALL MPI_Wait(handle(3),status,ierror)
            ! Right
            CALL MPI_Wait(handle(4),status,ierror)

            ! Setting Ghost points --------------------------- !
            ! Up
            IF (dest(1) /= MPI_PROC_NULL) THEN
                uk(i_min:i_max,j_max+1) = bur(1:i_size)
            END IF
            ! Down
            IF (dest(2) /= MPI_PROC_NULL) THEN
                uk(i_min:i_max,j_min-1) = bdr(1:i_size)
            END IF
            ! Left
            IF (dest(3) /= MPI_PROC_NULL) THEN
                uk(i_min-1,j_min:j_max) = blr(1:i_size)
            END IF
            ! Right
            IF (dest(4) /= MPI_PROC_NULL) THEN
                uk(i_max+1,j_min:j_max) = brr(1:i_size)
            END IF
    END SUBROUTINE UPDATE_GHOST

    SUBROUTINE RB_GS_SEQ
        i_min=1; i_max = N
        i_bu=1; i_bd=1
        j_min=1; j_max = N
        j_bl=1; j_br=1
        WRITE(*,*)
        WRITE(*,*) "Starting RB GS SEQ iterations. (k_max=",k_max," N=",N," d_min=",d_min,")"

        DO k = 1,k_max
            d = 0d0

            DO j= j_min+j_l,j_max-j_r
                DO i = i_min+i_d,i_max-i_u
                    IF (MOD(i+j,2)==0) THEN
                        uk(i,j) = (uk(i,j-1)+uk(i,j+1)+uk(i-1,j)+uk(i+1,j)+fdx2(i,j))*25d-2
                    END IF
                END DO
            END DO

            ! BLACK grid ---------------------------------------------------- !
            ! Send and recive bond on BLACK grid
            DO j= j_min+j_l,j_max-j_r
                DO i = i_min+i_d,i_max-i_u
                    IF (MOD(i+j,2)/=0) THEN
                        uk(i,j) = (uk(i,j-1)+uk(i,j+1)+uk(i-1,j)+uk(i+1,j)+fdx2(i,j))*25d-2
                    END IF
                END DO
            END DO

            ! Calculating residual
            d = get_residual()

            ! Build convergence cretia
            IF (d < d_min.and.(k > 10)) THEN
                WRITE(*,*) "The solver converged after: ", k, "Iterations (k_max: ", k_max, ")"
                WRITE(*,"(A,ES8.2E2,A,ES8.2E2,A)") " d: ", d, " (d_min:",d_min,")"
                exit
            end if

            IF (MOD(k,mod_state)==0.and.show_state) THEN
                WRITE(*,"(A,I5,A,ES8.2E2,A)") " Solution is not converged yet. (k= ",k,", d= ",d,")"
            end if
        end do



        IF (k > k_max) THEN
            WRITE(*,*)
            WRITE(*,*) "!! WARNING: The solver did NOT converge within k_max:", k_max, " !!"
            WRITE(*,"(A,ES8.2E2,A,ES8.2E2,A)") " d: ", d, " (d_min:",d_min,")"
            WRITE(*,*)
        end if

    END SUBROUTINE RB_GS_SEQ

    SUBROUTINE GS_SEQ
        i_min=1; i_max = N
        i_bu=1; i_bd=1
        j_min=1; j_max = N
        j_bl=1; j_br=1

        WRITE(*,*)
        WRITE(*,*) "Starting GS SEQ iterations. (k_max=",k_max," N=",N," d_min=",d_min,")"

        DO k = 1,k_max
            DO j=2,N-1
                DO i=2,N-1
                    uk(i,j) = (uk(i,j-1)+uk(i,j+1)+uk(i-1,j)+uk(i+1,j)+fdx2(i,j))*25d-2
                END DO
            END DO

            ! Calculating residual
            d = get_residual()

            ! Build convergence cretia
            IF (d < d_min.and.(k > 10)) THEN
                WRITE(*,*) "The solver converged after: ", k, "Iterations (k_max: ", k_max, ")"
                WRITE(*,"(A,ES8.2E2,A,ES8.2E2,A)") " d: ", d, " (d_min:",d_min,")"
                exit
            end if

            IF (MOD(k,mod_state)==0.and.show_state) THEN
                WRITE(*,"(A,I4,A,ES8.2E2,A)") " Solution is not converged yet. (k= ",k,", d= ",d,")"
            end if
        end do

        IF (k > k_max) THEN
            WRITE(*,*)
            WRITE(*,*) "!! WARNING: The solver did NOT converge within k_max:", k_max, " !!"
            WRITE(*,"(A,ES8.2E2,A,ES8.2E2,A)") " d: ", d, " (d_min:",d_min,")"
            WRITE(*,*)
        end if

    END SUBROUTINE GS_SEQ

    SUBROUTINE initilize
        CHARACTER(LEN=200) :: temp_name

        ! Default MPI values (Only used in seq mode)
        rank = 0

        ! Set N, k_max, d_min (Reads from argument or sets default value)
        N = read_vari_arg("N",100)
        k_max = read_vari_arg("k_max",100)
        d_min = read_vari_arg("d_min",1d-2)

        ! Setting problem
        problem = read_vari_arg("problem",1)

        ! Setting Solver Type
        solver_type = read_vari_arg("solver_type",1)

        ! Setting output filename
        filename = "A"
        IF (filename /= read_vari_arg_char("filesub",filename,.false.)) THEN
            temp_name = read_vari_arg_char("filesub",filename,.true.)
            WRITE(filename,"(A,A,A,I4.4,A)") "DATA/",TRIM(temp_name),"_N",N,".dat"

        ELSE
            WRITE(filename,"(A,I4.4,A)") "DATA/DATA_N",N,".dat"
        end if
        filename = read_vari_arg("filename",filename,.true.)

        ! Show state
        show_state = read_vari_arg("show_state",.true.)

        ! Setting mod_state
        mod_state = read_vari_arg("mod_state",20)

        ! Writing the full matrix in output file
        write_mat = read_vari_arg("write_mat",.true.)

        ! Grid spacing
        dx = 2d0/(N-1)
        dx2 = (dx)**2

        ! Allocate fields (The fiels are initilized when setting bounday conditions)
        ALLOCATE(uk(N,N))
        ALLOCATE(fdx2(N,N))

    END SUBROUTINE initilize

    SUBROUTINE set_bond
        uk = -1d0
        uk(1:N,1) = 0d0
        uk(1:N,N) = 0d0
        uk(1,1:N) = 0d0
        uk(N,1:N) = 0d0
    END SUBROUTINE set_bond

    SUBROUTINE set_f
        fdx2 = -dx2
    END SUBROUTINE set_f

    FUNCTION get_residual()
        REAL(MK) :: get_residual
        d = 0d0
        DO j= js,je
            DO i = is,ie
                d = d + abs((uk(i,j-1)+uk(i,j+1)+uk(i-1,j)+uk(i+1,j)-4d0*uk(i,j))/dx2-1d0)
            END DO
        END DO
        get_residual = d
    END FUNCTION get_residual

    SUBROUTINE write_matrix(mat, filename, write_mat)
        REAL(MK), DIMENSION(:,:), INTENT(IN) :: mat
        CHARACTER(LEN=*), INTENT(IN) :: filename
        LOGICAL :: write_mat
        INTEGER :: i,j

        OPEN(61, FILE=filename, action="write")

        ! Writing header
        WRITE(61,*) "N=", N
        WRITE(61,*) "k_max=", k_max
        WRITE(61,*) "k=", k
        WRITE(61,*) "d_min=", d_min
        WRITE(61,*) "d=", d
        WRITE(61,*) "Wall_time=", wall_time

        IF (k>k_max) THEN
            WRITE(61,*) "State= NOT converged"
        ELSE
            WRITE(61,*) "State= Converged"
        end if

        IF (solver_type == 2) THEN
            WRITE(61,*) "solver_type= Red-Black GS SEQ"
        ELSE IF (solver_type == 3) THEN
            WRITE(61,*) "solver_type= GS SEQ"
        ELSE IF (solver_type == 1) THEN
            WRITE(61,*) "solver_type= Red-Black GS Pal"
        end if

        IF (problem == 2) THEN
            WRITE(61,*) "problem= Harmonic"
        ELSE
            WRITE(61,*) "problem= Problem"
        end if

        IF (write_mat) THEN
            WRITE(61,*) "write_mat= True"
        ELSE
            WRITE(61,*) "write_mat= False"
        end if

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

    SUBROUTINE PRINT_RANK_AND_NAI(rank_in)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: rank_in
        IF (rank == rank_in) THEN
            PRINT *, "Rank",rank, " Coords", coords, ": i_min", i_min,&
                    " i_max", i_max,"j_min",j_min,"j_max",j_max
            PRINT*,"j_u",j_u, "j_d", j_d, "i_l",i_l,"i_r",i_r
            PRINT*,"is",is, " ie",ie, " js", js, " je", je
            ! Up
            IF (dest(1) /= MPI_PROC_NULL) THEN
                CALL MPI_Cart_coords(cart_comm,dest(1),ndim,coor_temp,ierror)
                PRINT*,"Up, rank:", dest(1), " coo:", coor_temp
            END IF
            ! Down
            IF (dest(2) /= MPI_PROC_NULL) THEN
                CALL MPI_Cart_coords(cart_comm,dest(2),ndim,coor_temp,ierror)
                PRINT*,"Down, rank:", dest(2), " coo:", coor_temp
            END IF
            ! Left
            IF (dest(3) /= MPI_PROC_NULL) THEN
                CALL MPI_Cart_coords(cart_comm,dest(3),ndim,coor_temp,ierror)
                PRINT*,"Left, rank:", dest(3), " coo:", coor_temp
            END IF
            ! Right
            IF (dest(4) /= MPI_PROC_NULL) THEN
                CALL MPI_Cart_coords(cart_comm,dest(4),ndim,coor_temp,ierror)
                PRINT*,"Right, rank:", dest(4), " coo:", coor_temp
            END IF

        END IF
    END SUBROUTINE PRINT_RANK_AND_NAI

    FUNCTION logical2integer(loc_in)
        INTEGER :: logical2integer
        LOGICAL, INTENT(IN) :: loc_in
        IF (loc_in) THEN
            logical2integer = 1
        ELSE
            logical2integer = 0
        END IF
    END FUNCTION logical2integer
end module m_GS_solver