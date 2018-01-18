program GS_SOLVER
    ! Declearing variabels
    USE m_gs_solver

    ! Initlize variabels
    CALL initilize

    ! Set boundary condition
    IF (problem == 1) THEN
        CALL set_bond
    end if

    ! Set source field
    IF (problem == 1) THEN
        CALL set_f
    end if

    ! Solve via Jacobi solver
    wall_time = MPI_WTIME()
    IF (solver_type == 2) THEN
        CALL RB_GS_SEQ
    ELSE IF (solver_type == 3) THEN
        CALL GS_SEQ
    ELSE IF (solver_type == 1) THEN
        CALL RB_GS_PAL
    end if
    WRITE(*,*) "Wall time=", MPI_WTIME()-wall_time, " secounds"
    ! Write out soulution to text file
    IF (rank == 0) THEN
       CALL write_matrix(uk,filename,write_mat)
    END IF


end program GS_SOLVER