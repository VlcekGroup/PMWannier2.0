program wannier
      use OMP_LIB
      use commvar
      implicit none
      
      call getarg(1,fname1)
      call print_start_time
      call write_header
      start_all = omp_get_wtime()
      call read_cnt
      call read_input
      call read_bin
      call write_read_info
      call mat_alloc1
      call OMP_set_num_threads(nthreads)
      call wannier_core
      finish_all = omp_get_wtime()
      call print_end_time
      write(6,*) '--------------------------------------------------------------'; call flush(6)
      write(6,*) "Job done and total wall time:", finish_all-start_all; call flush(6)
 end program wannier
