program GS_SOLVER
    ! Declearing variabels
    USE m_gs_solver


    WRITE(*,*) "--------------------------- STARTING SOLVER ------------------------------------"
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
    IF (solver_type == 2) THEN
        CALL RB_GS_SEQ
    ELSE IF (solver_type == 3) THEN
        CALL GS_SEQ
    ELSE IF (solver_type == 1) THEN
        CALL RB_GS_PAL(N,k_max,d_min,uk,fdx2,wall_time,k,d,show_state,mod_state)
    end if

    ! Write out soulution to text file
    CALL write_matrix(uk,filename,write_mat)

    WRITE(*,*) "------------------------------ ENDING SOLVER -----------------------------------"
    WRITE(*,*)
end program GS_SOLVER