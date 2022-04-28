subroutine steepest_descent_full
      use OMP_LIB
      use commvar
      implicit none
      integer::i,j,k,l,m,count_iter
      integer::st, indx, iter_eff, i_hit
      real*8::sum_itmd, OBJ_0, OBJ_1, d_OBJ
      real*8::start,finish
      real*8 :: start_head, finish_tail
      real*8, allocatable :: Q1_itmd(:,:), U_dag(:,:), Q0_itmd(:,:), U_itmd(:,:)
      
      sum_itmd = 0d0
      OBJ_0 = 0d0
      OBJ_1 = 0d0
      d_OBJ = 0d0
      count_iter = 0
      iter_eff = 0
      i_hit = 0
      allocate(Q1_itmd(ns,ns),stat=st)
      if(st/=0) stop "error: allocation of Q1_itmd!"
      allocate(Q0_itmd(ns,ns),stat=st)
      if(st/=0) stop "error: allocation of Q0_itmd!"
      allocate(U_dag(ns,ns),stat=st)
      if(st/=0) stop "error: allocation of U_dag!"
      allocate(U_itmd(ns,ns),stat=st)
      if(st/=0) stop "error: allocation of U_itmd!"

      start_head = omp_get_wtime()

      do
        if(count_iter.eq.0) then
          start = omp_get_wtime()
          !$OMP parallel
          !$OMP do schedule(dynamic)
           do j=1,ns
             do k=1,ns
               do i=1,noa
                 Q_A0(i,j,k) = sum(CO(:,j)*AWF(:,i)*CO(:,k))*dV
               end do
             end do
           end do
          !$OMP end do
          !$OMP end parallel
     
          !$OMP parallel shared(OBJ_0)
          !$OMP do reduction(+:OBJ_0) schedule(dynamic)
              do i=1,ns
                 OBJ_0 = OBJ_0 + sum(Q_A0(:,i,i)**2d0)
              end do
          !$OMP end do
          !$OMP end parallel
        
          !IF(allocated(AWF)) Deallocate(AWF)
          finish = omp_get_wtime()

          write(6,*) 'step, OBJ, dOBJ, wall_time, and delta_t:'; call flush(6)
          write(6,*) count_iter,OBJ_0,0,finish-start,delta_t; call flush(6)
        end if

        start = omp_get_wtime()
        !$OMP parallel
        !$OMP do schedule(dynamic)
         do j=1,ns
           do i=1,ns
             A(i,j) = delta_t*sum(Q_A0(:,j,i)*(Q_A0(:,j,j)-Q_A0(:,i,i))-Q_A0(:,i,j)*(Q_A0(:,i,i)-Q_A0(:,j,j)))
           end do
         end do
        !$OMP end do
        !$OMP end parallel

        call r8mat_expm1 (ns, A, U)
        U_dag = transpose(U)

        !get Q_A1
        do i=1, noa
          Q0_itmd(:,:) = Q_A0(i,:,:)
          Q1_itmd = matmul(Q0_itmd, U)
          Q1_itmd = matmul(U_dag, Q1_itmd)
          Q_A1(i,:,:) = Q1_itmd(:,:)
        end do


        !$OMP parallel shared(OBJ_1)
        !$OMP do reduction(+:OBJ_1) schedule(dynamic)
          do i=1,ns
             OBJ_1 = OBJ_1 + sum(Q_A1(:,i,i)**2d0)
          end do
        !$OMP end do
        !$OMP end parallel

        d_OBJ = OBJ_1-OBJ_0
        count_iter = count_iter + 1
        
        finish = omp_get_wtime() 
        write(6,*) count_iter,OBJ_1,d_OBJ,finish-start,delta_t; call flush(6)

        if(count_iter.eq.max_iter) then
          if(d_OBJ.gt.0) then
            iter_eff = iter_eff + 1
            if(iter_eff==1) then
              U_itmd = U_dag
            else
              U_itmd = matmul(U_dag,U_itmd)
            end if
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
            i_hit = i_hit + 1
            if(i_hit==nhit_cong_in) exit
          else
            i_hit = 0
          end if

          Q_A0 = Q_A1
          OBJ_0 = OBJ_1
          OBJ_1 = 0d0

          write(fname2,'(a,i3.3,a)') 'U_tmp.txt'
          open(222,file=trim(fname2),form='formatted',status='replace',action='write')
            do i=1,ns
              write(222,*) U_itmd(i,:)
            end do
          close(222)

        else
         delta_t = delta_t/corr_fac
         OBJ_1 = 0d0
        end if
 end do

if(d_OBJ.gt.0) then
  !$OMP parallel
  !$OMP do
  do i=1, ns
    orb_con(i) = sum(Q_A1(:,i,i))
  end do
  !$OMP end do
  !$OMP end parallel
else
  !$OMP parallel
  !$OMP do
  do i=1, ns
    orb_con(i) = sum(Q_A0(:,i,i))
  end do
  !$OMP end do
  !$OMP end parallel
end if

finish_tail = omp_get_wtime()
if(allocated(Q1_itmd)) deallocate(Q1_itmd)
if(allocated(Q0_itmd)) deallocate(Q0_itmd)
if(allocated(U_dag)) deallocate(U_dag)
 
write(6,*) '--------------------------------------------------------------'; call flush(6)
write(6,*) "End of iterations!"; call flush(6)
write(6,*) "Wannier function index and its locality"; call flush(6)
do i=1, ns
  indx = maxloc(orb_con,1)
  write(6,*) indx, orb_con(indx); call flush(6)
  orb_con(indx) = 0d0
end do

write(6,*) '--------------------------------------------------------------'; call flush(6)
write(6,*) "iter_tot, iter_eff, and wtime:",count_iter, iter_eff, finish_tail-start_head; call flush(6)

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

end subroutine steepest_descent_full
