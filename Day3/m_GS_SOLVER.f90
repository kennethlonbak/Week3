MODULE m_GS_solver
    INTERFACE read_vari_arg
        MODULE PROCEDURE read_vari_arg_integer, read_vari_arg_double, read_vari_arg_char, read_vari_arg_logical
    END INTERFACE read_vari_arg

    INTEGER, PARAMETER :: MK = KIND(0d0)
    ! Variabels that can be changed by user
    INTEGER :: N, k_max, k, N_th
    REAL(MK) :: d_min, d
    CHARACTER(LEN=200) :: filename
    ! Variabels to be used by the solver
    REAL(MK), DIMENSION(:,:), ALLOCATABLE :: uk, fdx2
    REAL(MK) :: dx, dx2, wall_time
    INTEGER :: solver_type ! Solver type: 1= Jacobi, 2=Gauss-Seidel
    INTEGER :: mod_state, problem
    LOGICAL :: show_state, write_mat

    CONTAINS
    SUBROUTINE RB_GS_PAL(N,k_max,d_min,uk,fdx2,wall_time,k,d,show_state,mod_state)

        IMPLICIT NONE
        INCLUDE 'mpif.h'
        INTEGER, INTENT(IN) :: N, k_max
        REAL(MK), DIMENSION(N,N), INTENT(INOUT) :: uk
        REAL(MK), DIMENSION(N,N), INTENT(IN) :: fdx2
        REAL(MK), INTENT(IN) :: d_min
        INTEGER, INTENT(OUT), OPTIONAL :: k
        REAL(MK), INTENT(OUT), OPTIONAL :: d, wall_time
        LOGICAL :: show_state
        INTEGER :: mod_state
        INTEGER :: i, j, rb_i, rank, Nproc, ierror

        ! Initialize MPI
        CALL MPI_INIT(ierror)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, Nproc, ierror)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

        ! Determine grid layout
        SELECT CASE (Nproc)
            CASE (1)
                WRITE(*,*) "Can not run use Nproc=1, program is terminated"
                CALL MPI_FINALIZE(ierror)
                STOP
            CASE DEFAULT
                WRITE(*,*) "NOT implemented yet, Rank ", rank," Nproc ", Nproc
                CALL MPI_FINALIZE(ierror)
                STOP
        END SELECT


        WRITE(*,*)
        WRITE(*,*) "Starting RB GS PAL iterations. (Rank",rank, "k_max=",k_max," N=",N," d_min=",d_min,")"

        DO k = 1,k_max
            d = 0d0

            DO i=2,N-1
                IF (MOD(i,2) == 0) THEN
                    rb_i = 0
                ELSE
                    rb_i = 1
                END IF

                DO j=2+rb_i,N-1+rb_i,2
                    uk(i,j) = (uk(i,j-1)+uk(i,j+1)+uk(i-1,j)+uk(i+1,j)+fdx2(i,j))*25d-2
                    d = d + (uk(i,j)-uk(i,j))**2
                END DO
            END DO
            DO i=2,N-1
                IF (MOD(i,2) == 0) THEN
                    rb_i = 1
                ELSE
                    rb_i = 0
                END IF
                DO j=2+rb_i,N-1+rb_i,2
                    uk(i,j) = (uk(i,j-1)+uk(i,j+1)+uk(i-1,j)+uk(i+1,j)+fdx2(i,j))*25d-2
                    d = d + (uk(i,j)-uk(i,j))**2
                END DO
            END DO
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
        WRITE(*,*) "Wall time=", wall_time, " secounds"
        CALL MPI_FINALIZE(ierror)
    END SUBROUTINE RB_GS_PAL

    SUBROUTINE RB_GS_SEQ
        INTEGER :: i, j, rb_i


        WRITE(*,*)
        WRITE(*,*) "Starting RB GS SEQ iterations. (k_max=",k_max," N=",N," d_min=",d_min,")"

        DO k = 1,k_max
            d = 0d0

            DO i=2,N-1
                IF (MOD(i,2) == 0) THEN
                    rb_i = 0
                ELSE
                    rb_i = 1
                END IF

                DO j=2+rb_i,N-1+rb_i,2
                    uk(i,j) = (uk(i,j-1)+uk(i,j+1)+uk(i-1,j)+uk(i+1,j)+fdx2(i,j))*25d-2
                    d = d + (uk(i,j)-uk(i,j))**2
                END DO
            END DO
            DO i=2,N-1
                IF (MOD(i,2) == 0) THEN
                    rb_i = 1
                ELSE
                    rb_i = 0
                END IF
                DO j=2+rb_i,N-1+rb_i,2
                    uk(i,j) = (uk(i,j-1)+uk(i,j+1)+uk(i-1,j)+uk(i+1,j)+fdx2(i,j))*25d-2
                    d = d + (uk(i,j)-uk(i,j))**2
                END DO
            END DO
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
        WRITE(*,*) "Wall time=", wall_time, " secounds"

    END SUBROUTINE RB_GS_SEQ

    SUBROUTINE GS_SEQ
        INTEGER :: i, j


        WRITE(*,*)
        WRITE(*,*) "Starting GS SEQ iterations. (k_max=",k_max," N=",N," d_min=",d_min,")"

        DO k = 1,k_max
            DO j=2,N-1
                DO i=2,N-1
                    uk(i,j) = (uk(i,j-1)+uk(i,j+1)+uk(i-1,j)+uk(i+1,j)+fdx2(i,j))*25d-2
                END DO
            END DO
            d = 0d0
            DO i=2,N-1
                DO j=2,N-1
                    d = d + abs((uk(i,j-1)+uk(i,j+1)+uk(i-1,j)+uk(i+1,j)-4d0*uk(i,j))/dx2-1)
                END DO
            END DO
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
        WRITE(*,*) "Wall time=", wall_time, " secounds"

    END SUBROUTINE GS_SEQ

    SUBROUTINE initilize
        CHARACTER(LEN=200) :: temp_name

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
            DO i = 1,N
                DO j = 1,N
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
end module m_GS_solver