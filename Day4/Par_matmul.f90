program Par_matmul
    ! Declaring variabels
    USE m_Par_matmul

    ! Start MPI
    CALL MPI_INIT(ierror)

    ! Get input parameters (N,)
    CALL initialize

    ! Generate Caratisian grid ( Checking if sqrt(N) is a whole number)
    CALL create_cart

    ! Create process index
    CALL get_index_from_coord(coord,is,ie,js,je,isize,jsize)

    ! Allocate sub matrices
    ALLOCATE(C(isize,jsize),A_own(isize,jsize),B_own(isize,jsize))
    ALLOCATE(A_recv(isize,jsize),B_recv(isize,jsize))
    ALLOCATE(sr_vec(isize))

    ! Checking the cartisian grid (output to terminal)
    CALL check_cart_grid

    ! Make row comunicator
    CALL MPI_Cart_sub(cart_comm,(/.false.,.true./),row_comm,ierror)
    CALL MPI_COMM_RANK(row_comm, rank_row, ierror)

    ! Make colum comunicator
    CALL MPI_Cart_sub(cart_comm,(/.true.,.false./),col_comm,ierror)
    CALL MPI_COMM_RANK(col_comm, rank_col, ierror)

    ! Generate/load sub matrix
    CALL gen_sub_matrix_cos
    !CALL gen_sub_matrix_rank
    CALL collect_A_tot
    CALL collect_B_tot

    ! Make k do loop
    DO k = 0, q-1
        ! Calculate k_bar
        k_bar = MOD(coord(1)+k,q)

        ! Brodcast A(i,k_bar) to each row
        !CALL test_broadcast((/0,0/), "Before")
        CALL broadcast_row_mat
        !CALL test_broadcast((/0,0/), "After")

        ! Compute matmul(A,B) and add to C
        C = C + matmul(A_recv,B_recv)

        ! Send B to matrix above
        !CALL test_b_send((/0,0/), "Before")
        CALL send_B
        !CALL test_b_send((/0,0/), "After")
    END DO

    ! Collect C at rank 0 and print it out.
    IF (rank==0) THEN
        CALL collect_C_tot
        CALL compare_with_seq
        CALL write_matrix(C_tot, "C_tot.dat", .true.)
    ELSE
        CALL send_mat(C, isize, jsize, cart_comm, 0, is, ie, js, je)
    END IF

    CALL MPI_FINALIZE(ierror)

end program Par_matmul