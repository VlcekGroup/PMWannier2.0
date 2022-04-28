subroutine make_CO_sto(i_mc)
      use OMP_LIB
      use commvar
      implicit none
      integer :: i, st, i_sto, i_mc, j
      integer :: n_seed, ns_gs
      integer, allocatable :: seed_no(:)
      real*8, allocatable :: rand_r8(:,:)
      real*8 :: norm_fac, ovlp
      real*8, allocatable :: ovlp_check(:,:), orb_sto(:), orb_sub(:)
      real*8 :: start_sto, end_sto
      
      if(i_mc==i_mc_step0) then
        call read_AO2
        if(ntot_sto.lt.0) stop "number of states in total for sto_wan must be provided for this version!"
        if(ntot_sto.gt.nstates2) stop "asking for too many states for stochastic wannierization!"
        if(n_det.lt.0) then
          if(loc_wan) then
            n_det = ns_lw
          elseif(frag_wan) then
            n_det = ns_fw(1)
          end if
        end if
        n_tot2 = n_det + n_sto
        ns = n_tot2
      end if

      allocate(orb_sto(nx*ny*nz), stat=st)
      if(st/=0) stop "error: allocation of orb_sto array"
      allocate(orb_sub(nx*ny*nz), stat=st)
      if(st/=0) stop "error: allocation of orb_sub1 array"     


      allocate(CO(nx*ny*nz,ns),stat=st)
      if(st/=0) stop "error: allocation of CO matrix!"
      !allocate(rand_st(nx*ny*nz,n_sto),stat=st)
      !if(st/=0) stop "error: allocation of rand_st matrix!"

      if(i_mc==i_mc_step0) then
        allocate(det_st(nx*ny*nz,n_det),stat=st)
        if(st/=0) stop "error: allocation of det_st matrix!"
        !write(6,*) "calling get_det!"; call flush(6)
        call get_det
        do i=1, n_det
          CO(:,i) = det_st(:,i)
        end do
        deallocate(det_st)
      else
        !call orb_indx_LO
        do i=1, n_det
          CO(:,i) = LO(:,orb_indx(i))
        end do
      end if

      write(6,*) "done preparing the deterministic states in CO"; call flush(6)
      
      start_sto = omp_get_wtime()
      call random_seed(size=n_seed)
      allocate(seed_no(n_seed),stat=st)
      if(st/=0) stop "error: allocation of seed_no array!"
      write(6,*) "Done allocation seed_no!"; call flush(6)

      allocate(rand_r8(ntot_sto,n_sto),stat=st)
      if(st/=0) stop "error: allocation of rand_r8 matrix!"
      seed_no = n_det * 933273 + i_mc
      call random_seed(put=seed_no)
      call random_number(rand_r8)

      ns_gs = n_det
      do i_sto = 1, n_sto
        orb_sto(:) = 0d0
        !$OMP parallel shared(orb_sto)
        !$OMP do reduction(+:orb_sto) schedule(dynamic)
        do j = 1, ntot_sto
          orb_sto(:) = orb_sto(:) + rand_r8(j,i_sto) * AO2(:,j)
        end do
        !$OMP end do
        !$OMP end parallel
        orb_sub(:) = 0d0
        !$OMP parallel private(ovlp) shared(orb_sub)
        !$OMP do reduction(+:orb_sub) schedule(dynamic)
        do j = 1, ns_gs
          ovlp = sum(CO(:,j)*orb_sto(:)) / sum(abs(CO(:,j))**2d0)
          orb_sub(:) = orb_sub(:) + ovlp*CO(:,j)
        end do
        !$OMP end do
        !$OMP end parallel
        orb_sto(:) = orb_sto(:) - orb_sub(:)
        norm_fac = sum(abs(orb_sto(:))**2d0)*dV
        orb_sto(:) = (1d0/dsqrt(norm_fac)) * orb_sto(:)
        if(any(orb_sto/=orb_sto)) stop "NaN problem of orb_sto!"
        CO(:,n_det+i_sto) = orb_sto(:)
        ns_gs = ns_gs + 1
      end do
      end_sto = omp_get_wtime()

      write(6,*) "done preparing and orthgonalizing the stochastic states in CO"; call flush(6)
      write(6,*) "Time spent in preparing stochastic states:", end_sto-start_sto; call flush(6)
      
      if(allocated(orb_sub)) deallocate(orb_sub)
      if(allocated(seed_no)) deallocate(seed_no)
      if(allocated(rand_r8)) deallocate(rand_r8)
      if(allocated(orb_sto)) deallocate(orb_sto)
      
      if(check_CO) then
        if(i_mc==i_mc_step0) then
          !check orthonormality of CO
          allocate(ovlp_check(ns,ns),stat=st)
          if(st/=0) stop "allocation of ovlp_check_matrix!"
          do i=1, ns
            do j=1, ns
              ovlp_check(i,j) = sum(CO(:,i)*CO(:,j))*dV
            end do
          end do
          write(123,*) "MC_step:", i_mc; call flush(123)
          do i=1, ns
            write(123,*) ovlp_check(i,:); call flush(123)
          end do

          if(allocated(ovlp_check)) deallocate(ovlp_check)
        end if
      end if
end subroutine make_CO_sto

subroutine get_det
     use OMP_LIB
     use commvar
     implicit none
     integer :: i,j,k,st
     integer :: ns_sub, na
     real*8,allocatable :: locality(:), Q_diag(:,:)
     integer,allocatable :: label(:)
     
     if(restart) then
       ns_sub = nstates
     else
       ns_sub = ntot_sto
     end if

     if(ns_sub.lt.n_det) stop "not enough states to construct det_st matrix!"

     if(loc_wan) then
       na = noa_lw
     elseif(frag_wan) then
       na = 1
     end if

     allocate(label(na),stat=st)
     if(st/=0) stop "error: allocation of the label array"
     allocate(locality(ns_sub),stat=st)
     if(st/=0) stop "error: allcatiion of locality matrix"
     allocate(Q_diag(na,ns_sub),stat=st)
     if(st/=0) stop "error: allocation of Q_diag matrix"

     if(loc_wan) then
       label = atom_label
     else
       do i=1, na
         label(i) = i
       end do
     end if
     write(6,*) "done with label"; call flush(6)
     
     if(frag_wan) then
       !$OMP parallel
       !$OMP do schedule(dynamic)
       do j=1, ns_sub
          Q_diag(1,j) = sum(AO(:,j)*FWF(:,1)*AO(:,j))*dV
       end do
       !$OMP end do
       !$OMP end parallel
     elseif(loc_wan) then
       !$OMP parallel
       !$OMP do schedule(dynamic)
       do j=1, ns_sub
          do i=1, na
            Q_diag(i,j) = sum(AO(:,j)*AWF(:,label(i))*AO(:,j))*dV
          end do
       end do
       !$OMP end do
       !$OMP end parallel
     end if

     write(6,*) "done with Q_diag"; call flush(6)

     !$OMP parallel
     !$OMP do
       do i=1, ns_sub
         locality(i) = sum(Q_diag(:,i))
       end do
     !$OMP end do
     !$OMP end parallel
     write(6,*) "done with locality"; call flush(6)

     allocate(state_indx(n_det),stat=st)
     if(st/=0) stop "error: allocation of state_indx array"

     do i=1, n_det
       state_indx(i) = maxloc(locality,1)
       locality(state_indx(i)) = -1d0
     end do

     k=1
     do i=1, ns_sub
       if(any(state_indx==i)) then
         det_st(:,k) = AO(:,i)
         k = k + 1
       end if
     end do
     if((k-1)/=n_det) stop "error: generating the det_st matrix"

     write(6,*) "done with det_st"; call flush(6)

     if(allocated(locality)) deallocate(locality)
     if(allocated(Q_diag)) deallocate(Q_diag)
     if(allocated(label)) deallocate(label)
     if(allocated(AO)) deallocate(AO)
     if(allocated(state_indx)) deallocate(state_indx)
end subroutine get_det
