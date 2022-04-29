!                     
!                     
!    The routine(s) in this file are a part of the  
!                     PMWannier2.0                 
!    suite, developed 2020-2022, and copyrighted    
!    to the authors: Guorong Weng and Vojtech Vlcek,
!    at the University of California, Santa Barbara.
!                                                   
!                                                   
!   If you use or modify any part of this routine   
!   the header should be kept and unmodified.          
!                                                   
!                                                   
! 
subroutine steepest_descent_fw_sw(i_step)
      use OMP_LIB
      use commvar
      implicit none
      integer::i,j,k,l,m,count_iter, n_sub
      real*8::sum_itmd, OBJ_0, OBJ_1, d_OBJ
      real*8::start,finish, start1
      real*8:: start_head, finish_tail
      real*8 :: delta_t_sw
      integer :: iter_eff, i_step, ns_fw2, st, i_hit
      integer, allocatable :: orb_indx0(:)
      real*8, allocatable :: Q_A0_fw(:,:), Q_A1_fw(:,:), U_dag(:,:), U_itmd(:,:)
     
      sum_itmd = 0d0
      OBJ_0 = 0d0
      OBJ_1 = 0d0
      d_OBJ = 0d0
      count_iter = 0
      delta_t_sw = delta_t

      iter_eff = 0
      ns_fw2 = ns_fw(1)
      i_hit = 0

      if(sto_wan) then
        n_sub = n_det
      else
        n_sub = n_c
      end if

      allocate(Q_A0_fw(ns,ns),stat=st)
      if(st/=0) stop "error: allocation of Q_A0_fw!"
      allocate(Q_A1_fw(ns,ns),stat=st)
      if(st/=0) stop "error: allocation of Q_A1_fw!"
      allocate(U_dag(ns,ns),stat=st)
      if(st/=0) stop "error: allocation of U_dag!"
      allocate(U_itmd(ns,ns),stat=st)
      if(st/=0) stop "error: allocation of U_itmd matrix!"

      start_head = omp_get_wtime()

      if(ns_fw2.gt.ns) stop "asking for too many states of interest."

      do
        if(count_iter.eq.0) then
          start1 = omp_get_wtime()
          !$OMP parallel
          !$OMP do schedule(dynamic)
            do j=1,ns
              do k=1,ns
                 Q_A0_fw(j,k) = sum(CO(:,j)*FWF(:,1)*CO(:,k))*dV
              end do
            end do
          !$OMP end do
          !$OMP end parallel
          finish = omp_get_wtime()
          !write(6,*) "time spent on Q_A0:", finish-start1; call flush(6)
   
          !$OMP parallel
          !$OMP do
            do i=1,ns
             orb_con(i) = Q_A0_fw(i,i)
            end do
          !$OMP end do
          !$OMP end parallel

          if(.not.allocated(orb_indx0)) then
            allocate(orb_indx0(n_sub),stat=st)
            if(st/=0) stop "error: allocation of orb_indx0"
          end if

          do i=1,n_sub
            forb_con(i) = maxval(orb_con)
            orb_indx0(i) = maxloc(orb_con,1)
            orb_con(orb_indx0(i)) = 0d0
          end do
   
          !$OMP parallel shared(OBJ_0)
          !$OMP do reduction(+:OBJ_0) schedule(dynamic)
            do j=1,ns_fw2
               OBJ_0 = OBJ_0 + Q_A0_fw(orb_indx0(j),orb_indx0(j))**2d0
            end do
          !$OMP end do
          !$OMP end parallel
          
          finish = omp_get_wtime()
          write(6,*) 'step, OBJ, dOBJ, wall_time, and delta_t:'; call flush(6)
          write(6,*) count_iter,OBJ_0,0,finish-start1,delta_t_sw; call flush(6)
        end if

        start1 = omp_get_wtime()
        !$OMP parallel
        !$OMP do schedule(dynamic)
         do j=1,ns
           do i=1,ns
              A(i,j) = delta_t_sw * (Q_A0_fw(j,i)*(Q_A0_fw(j,j)-Q_A0_fw(i,i))-Q_A0_fw(i,j)*(Q_A0_fw(i,i)-Q_A0_fw(j,j)))
           end do
         end do
        !$OMP end do
        !$OMP end parallel
        finish = omp_get_wtime()
        !write(6,*) "time spent on A:", finish-start1; call flush(6)

        start = omp_get_wtime()
        !get the unitary matrix U and U dagger
        call r8mat_expm1 (ns, A, U)

        U_dag = transpose(U)
        finish = omp_get_wtime()
        !write(6,*) "time spent on U and U_dag:", finish-start; call flush(6)

        start = omp_get_wtime()
        !get the Q_A1 matrix
        Q_A1_fw = matmul(Q_A0_fw,U)
        Q_A1_fw = matmul(U_dag,Q_A1_fw)
        finish = omp_get_wtime()
        !write(6,*) "time spent on Q_A1:", finish-start; call flush(6)

        !$OMP parallel
        !$OMP do
          do i=1,ns
             orb_con(i) = Q_A1_fw(i,i)
          end do
        !$OMP end do
        !$OMP end parallel

        do i=1,n_sub
          forb_con(i) = maxval(orb_con)
          orb_indx(i) = maxloc(orb_con,1)
          orb_con(orb_indx(i)) = -1d0
        end do

        !calculate the objective function OBJ_1
        !$OMP parallel shared(OBJ_1)
        !$OMP do reduction(+:OBJ_1) schedule(dynamic)
          do j=1,ns_fw2
             OBJ_1 = OBJ_1 + Q_A1_fw(orb_indx(j),orb_indx(j))**2d0
          end do
        !$OMP end do
        !$OMP end parallel

        finish = omp_get_wtime()
        d_OBJ = OBJ_1-OBJ_0
        count_iter = count_iter + 1

      if(mod(count_iter,step_prnt_info)==0) then
        write(6,*) count_iter,OBJ_1,d_OBJ,finish-start1,delta_t_sw; call flush(6)
      end if

      if(mod(count_iter,sfw)==0) then
        write(6,*) '--------------------------------------------------------------'
        write(6,*) 'orbital_index, locality'; call flush(6)
        do i=1,ns_fw2
          write(6,*) orb_indx(i),forb_con(i); call flush(6)
        end do
        write(6,*) '--------------------------------------------------------------'
      end if

      if(count_iter.eq.max_iter) then
        if(d_OBJ.gt.0) then
          iter_eff = iter_eff + 1
          if(iter_eff==1) then
            U_itmd = U_dag
          else
            U_itmd = matmul(U_dag,U_itmd)
          end if

          if(i_step==istep_0) then
            OBJ_m0 = OBJ_1
          else
            OBJ_m1 = OBJ_1
          end if
        else
          if(i_step==istep_0) then
            OBJ_m0 = OBJ_0
          else
            OBJ_m1 = OBJ_0
          end if
          orb_indx(:) = orb_indx0(:)
        end if
        exit
      end if

      if(d_OBJ.gt.0) then
        iter_eff = iter_eff + 1
        if(iter_eff==1) then
          U_itmd = U_dag
        else
          U_itmd = matmul(U_dag,U_itmd)
        end if

        if(d_OBJ.lt.con_thres_inner) then
          if(i_step==istep_0) then
            OBJ_m0 = OBJ_1
          else
            OBJ_m1 = OBJ_1
          end if
          i_hit = i_hit + 1
          if(i_hit==nhit_cong_in) exit
        else
          i_hit = 0
        end if

        Q_A0_fw = Q_A1_fw
        orb_indx0(:) = orb_indx(:)
        OBJ_0 = OBJ_1
        OBJ_1 = 0d0
      else !d_OBJ<0
        delta_t_sw = delta_t_sw/corr_fac
        OBJ_1 = 0d0
      end if

 end do
 
 if(allocated(orb_indx0)) deallocate(orb_indx0)
 if(allocated(U_dag)) deallocate(U_dag)
 if(allocated(Q_A0_fw)) deallocate(Q_A0_fw)
 if(allocated(Q_A1_fw)) deallocate(Q_A1_fw)

 iter_acc = iter_acc + count_iter
 iter_eff_acc = iter_eff_acc + iter_eff
 finish_tail = omp_get_wtime()
 
write(6,*) '--------------------------------------------------------------'; call flush(6)
write(6,*) "End of iterations!"; call flush(6)
write(6,*) "Wannier function index and its locality"; call flush(6)

do i=1,ns_fw2
  write(6,*) orb_indx(i),forb_con(i); call flush(6)
end do

write(6,*) '--------------------------------------------------------------'; call flush(6)
write(6,*) "i_step, iter_tot, iter_eff, and wtime:",i_step, count_iter, iter_eff, finish_tail-start_head; call flush(6)

write(6,*) '--------------------------------------------------------------'; call flush(6)
write(6,*) "Transforming into Localized States!"; call flush(6)
if(iter_eff.gt.0) then
  !$OMP parallel
  !$OMP do
  do j=1,nx*ny*nz
    do i=1,ns
       LO(j,i) = sum(U_itmd(i,:)*CO(j,:))
    end do
  end do
  !$OMP end do
  !$OMP end parallel
else
  LO = CO
end if

if(allocated(U_itmd)) deallocate(U_itmd)
end subroutine steepest_descent_fw_sw
