subroutine wannier_core
  use OMP_LIB
  use commvar
  implicit none
  real*8 :: finish_prep

  call awf_prep
  call mat_dealloc1
  finish_prep = omp_get_wtime()
  write(6,*) 'max_iter and delta_t from input:', max_iter, delta_t; call flush(6)
  write(6,*) 'Time spent preparing for steepest descent:', finish_prep-start_all; call flush(6)

  if(seq_wan.or.sto_wan) then
    if(seq_wan.and.sto_wan) stop "error: seq_wan can only do either deterministic or stochastic search!"
    if(sto_wan) then
      call sto_search
    else
      call det_search
    end if
  else !not seq_wan or sto_wan
    call simple_search
  end if
end subroutine wannier_core

subroutine det_search
  use OMP_LIB
  use commvar
  implicit none
  real*8 :: start_m, finish_m
  real*8 :: dOBJ
  integer :: i_frag, i_step_main, i_hit

  i_hit = 0
    if(frag_wan) then
      call frag_wan_prep
      if(n_frag.gt.1) stop 'error: seq_wan applies to one fragment only!'
    else
      if(.not.loc_wan) stop 'error: seq_wan applies to either loc_wan or frag_wan!'
    end if

      i_step_main = istep_0
      do
        start_m = omp_get_wtime()
        call make_CO_det(i_step_main)
        call mat_alloc2(1)
        write(6,*) '--------------------------------------------------------------'; call flush(6)
        write(6,*) 'Localizing states on selected atoms at i_step:', i_step_main; call flush(6)
        write(6,*) 'Start steepest descent process with the number of states:', ns; call flush(6)
        if(frag_wan) then
          call steepest_descent_fw_sw(i_step_main)
        elseif(loc_wan) then
          call steepest_descent_lw_sw(i_step_main)
        end if
        call update_core_rest(i_step_main)
        call mat_dealloc2
        if(mod(i_step_main,step_prnt_tmp)==0) call write_wf_tmp(i_step_main)
        finish_m = omp_get_wtime()
        write(6,*) 'wtime at i_step:', finish_m-start_m, i_step_main; call flush(6)

        if(i_step_main.gt.istep_0) then
          write(6,*) 'OBJ value at i_step:', OBJ_m1, i_step_main; call flush(6)
          write(6,*) '--------------------------------------------------------------'; call flush(6)
          dOBJ = OBJ_m1 - OBJ_m0
          OBJ_m0 = OBJ_m1
          write(6,*) 'Objective function value change and i_step:', dOBJ, i_step_main; call flush(6)
          write(6,*) '--------------------------------------------------------------'; call flush(6)

          if(dOBJ.gt.d_OBJ_max) then
            d_OBJ_max = dOBJ
          end if

          if(i_step_main.ge.nblock) then
            if(mod(i_step_main,nblock)==0) then
              if(d_OBJ_max.lt.con_thres_outer) then
                write(6,*) 'Objective function value is converged at i_step:', i_step_main; call flush(6)
                write(6,*) 'The maximum Objective function change at convergence:', d_OBJ_max; call flush(6)
                exit
              end if
            end if
          end if

          if(i_step_main.gt.nblock) then
            if(mod(i_step_main,nblock)==0) d_OBJ_max = -1d0
          end if

          if(i_step_main.ge.total_step) then
             write(6,*) 'Iteration reaches maximum steps'; call flush(6)
             exit
          end if

          i_step_main = i_step_main + 1

        else !i_step_main
          write(6,*) 'OBJ value at i_step:', OBJ_m0, i_step_main; call flush(6)
          write(6,*) '--------------------------------------------------------------'; call flush(6)
          if(i_step_main.ge.total_step) then
              write(6,*) 'Iteration reaches maximum steps'; call flush(6)
              exit
          end if
          i_step_main = i_step_main + 1
        end if !i_step_main

      end do
      write(6,*) 'Total iterations and total effective iterations', iter_acc, iter_eff_acc; call flush(6)
      call write_wf_full_seq(1)
      call write_wf_trun_seq(1)
      call mat_dealloc3(1)
end subroutine det_search

subroutine sto_search
    use OMP_LIB
    use commvar
    implicit none
    real*8 :: start_mc, finish_mc
    real*8 :: dOBJ
    integer :: i_frag, i_mc_main, i_hit

    i_hit = 0
    if(frag_wan) then
      call frag_wan_prep
      if(n_frag.gt.1) stop 'error: sto_wan applies to one fragment only!'
    else
      if(.not.loc_wan) stop 'error: sto_wan applies to either loc_wan or frag_wan!'
    end if

      i_mc_main = i_mc_step0
      do
        start_mc = omp_get_wtime()
        call make_CO_sto(i_mc_main)
        call mat_alloc2(1)
        write(6,*) '--------------------------------------------------------------'; call flush(6)
        write(6,*) 'Localizing states on selected atoms at MC_step:', i_mc_main; call flush(6)
        write(6,*) 'Start steepest descent process with the number of states:', ns; call flush(6)
        if(frag_wan) then
          call steepest_descent_fw_sw(i_mc_main)
        elseif(loc_wan) then
          call steepest_descent_lw_sw(i_mc_main)
        end if
        call mat_dealloc2
        if(mod(i_mc_main,step_prnt_tmp)==0) call write_wf_tmp_sto(i_mc_main)
        finish_mc = omp_get_wtime()
        write(6,*) 'wtime at MC_step:', finish_mc-start_mc, i_mc_main; call flush(6)

        if(i_mc_main.gt.i_mc_step0) then
          write(6,*) 'OBJ value at MC_step:', OBJ_m1, i_mc_main; call flush(6)
          write(6,*) '--------------------------------------------------------------'; call flush(6)
          dOBJ = OBJ_m1 - OBJ_m0
          write(6,*) 'Objective function value change and MC_step:', dOBJ, i_mc_main; call flush(6)
          write(6,*) '--------------------------------------------------------------'; call flush(6)
          if(dOBJ.gt.0d0) then
            OBJ_m0 = OBJ_m1
            if(dOBJ.le.con_thres_outer) then
              i_hit = i_hit + 1
              if(i_hit==nhit_cong_out) then
                write(6,*) 'Objective function value is converged at MC_step:', i_mc_main; call flush(6)
                exit
              else
                write(6,*) 'Convergence hitting times and MC_step:', i_hit, i_mc_main; call flush(6)
              end if
            elseif(i_mc_main.ge.mc_step) then
              write(6,*) 'Iteration reaches maximum mc_steps'; call flush(6)
              exit
            else
              i_hit = 0
            end if
            i_mc_main = i_mc_main + 1
            LO_sv = LO
          else !dOBJ<0
            i_mc_main = i_mc_main + 1
            LO = LO_sv
          end if !dOBJ
        else !i_mc==i_mc_step
          write(6,*) 'OBJ value at MC_step:', OBJ_m0, i_mc_main; call flush(6)
          write(6,*) '--------------------------------------------------------------'; call flush(6)
          if(i_mc_main.ge.mc_step) then
              write(6,*) 'Iteration reaches maximum mc_steps'; call flush(6)
              exit
          end if
          LO_sv = LO
          i_mc_main = i_mc_main + 1
        end if !i_mc

      end do
        write(6,*) 'Total iterations and total effective iterations', iter_acc, iter_eff_acc; call flush(6)
        call write_wf_full(1)
        call write_wf_trun(1)
        call mat_dealloc3(1)
end subroutine sto_search

subroutine simple_search
    use OMP_LIB
    use commvar
    implicit none
    integer :: i_frag
    if(frag_wan) then
      call frag_wan_prep
      do i_frag = 1, n_frag
        call make_CO(i_frag)
        call mat_alloc2(i_frag)
        write(6,*) '--------------------------------------------------------------'; call flush(6)
        write(6,*) 'Solving the wannier function localized on fragment', i_frag; call flush(6)
        write(6,*) 'Start steepest descent process with the number of states:', ns; call flush(6)

        call steepest_descent_fw(i_frag)
        call mat_dealloc2
        call write_wf_full(i_frag)
        call write_wf_trun(i_frag)
        call mat_dealloc3(i_frag)
      end do
    else !not frag_wan
      call make_CO(1)
      call mat_alloc2(1)
      write(6,*) '--------------------------------------------------------------'; call flush(6)
      write(6,*) 'Start steepest descent process with the number of states:', ns; call flush(6)

      if(loc_wan) then
        if(ns_lw.lt.0) stop 'not providing the number of states to search for in loc_wan!'
        call steepest_descent_lw
      else
        call steepest_descent_full
      end if
      call mat_dealloc2
      call write_wf_full(1)
      if(loc_wan) call write_wf_trun(1)
      call mat_dealloc3(1)
    end if
end subroutine simple_search 
