program SPMDriver

    use final
    use license
    use mod_idle
    use file_data
    use mpidata
    use error_data
    use precision_data
    use mod_simulation

    implicit none

    integer :: m

    real( tkind ) :: time_st, time_ed

    ! Check licensing if any
    call init_license

    ! Start up MPI
    call mpi_init( ierr )
    call mpi_comm_size( mpi_comm_world, nproc, ierr )
    call mpi_comm_rank( mpi_comm_world, iam, ierr )

    ! Redefine MPI types for the application
    SPModel_comm_world = mpi_comm_world
    int_type = mpi_integer
    dbl_type = mpi_double_precision
    lgc_type = mpi_logical
    max_type = mpi_max
    min_type = mpi_min
    sum_type = mpi_sum

    ! Get experiment file name
    if( iam == 0 ) then
        !write(6,*)'Author Tristan Salles'
        m = iargc()
        fin_noflat = ' '
        if(m < 1) then
            write(6,*)'use: mpirun -np X ./lecode <input-file-name> '
            write(6,*)'where X refers to the number of processors for the run. '
            attempt = ARG_FILE
        elseif(m == 1)then
            call getarg(1,fin_noflat)
            call read_input_to_flatten
        endif
    endif
    call completion

    time_st = mpi_wtime( )

    ! Initialisation driver
    call spm_initial
    if( iam == 0 )then
        write(6,*)'Initialisation phase done ...'
        write(6,*)
    endif

    ! Run driver
    call spm_run
    if( iam == 0 )then
        write(6,*)
        write(6,*)'Processes phase done ...'
        write(6,*)
    endif

    ! Finalisation driver
    call spm_final
    if( iam == 0 )then
        write(6,*)'Finalisation phase done ...'
        write(6,*)
    endif

    ! Output simulation duration
    time_ed = mpi_wtime( )
    if( iam == 0 ) print*,'Time elapse:',time_ed - time_st

    ! Finalise mpi
    call mpi_finalize( ierr )

end program SPMDriver
