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
subroutine make_CO_det(i_step)
      use OMP_LIB
      use commvar
      implicit none
      integer :: i, st, i_step, j
      integer :: indx_start1, indx_end1
      integer :: i_blk, i_blk_enter
      real*8, allocatable :: ovlp_check(:,:)
      real*8 :: start_co, end_co
      
      start_co = omp_get_wtime()

      if(i_step==istep_0) then
        if(ntot_seq.lt.0) stop "number of states in total for seq_wan must be provided!"
        if(ntot_seq.gt.nstates) stop "asking for too many states for sequential wannierization!"
        if(n_c.lt.0) then
          if(loc_wan) then
            n_c = ns_lw
          elseif(frag_wan) then
            n_c = ns_fw(1)
          end if
        end if
      end if


      if(i_step==istep_0) then
        allocate(core_st(nx*ny*nz,n_c),stat=st)
        if(st/=0) stop "error: allocation of core_st matrix!"
        n_rest = ntot_seq - n_c
        allocate(rest_st(nx*ny*nz,n_rest),stat=st)
        if(st/=0) stop "error: allocation of rest_st matrix!"
        call get_core_rest(i_step)

        nblock = int(n_rest/n_r)
        if(mod(n_rest,n_r)/=0) nblock = nblock + 1
        write(6,*) "Number of blocks in total:", nblock; call flush(6)
        
        allocate(nst_block(nblock),stat=st)
        if(st/=0) stop "error: allocation of nst_block!"
        
        nst_block(:) = n_r
        if(mod(n_rest,n_r)/=0) nst_block(nblock) = mod(n_rest,n_r)
        
        write(6,*) "Number of states in each block:"; call flush(6)
        write(6,*) nst_block(:); call flush(6)
      end if

      if(i_step.gt.nblock) then
        i_blk = mod(i_step,nblock)
      else
        i_blk = i_step
      end if

      n_tot2 = n_c + nst_block(i_blk)
      ns = n_tot2

      allocate(CO(nx*ny*nz,ns),stat=st)
      if(st/=0) stop "error: allocation of CO matrix!"

      indx_start1 = (i_blk-1)*n_r + 1
      indx_end1 = indx_start1 + nst_block(i_blk) - 1
      
      do i=1, n_c
        CO(:,i) = core_st(:,i)
      end do
      
      j = indx_start1
      do i = n_c+1, ns
        CO(:,i) = rest_st(:,j)
        j = j + 1
      end do
      if((j-1)/=indx_end1) stop "error: assigning rest_st to CO!"

      write(6,*) "done preparing CO"; call flush(6)

      end_co = omp_get_wtime()

      write(6,*) "Time spent in preparing CO:", end_co-start_co; call flush(6)
      
      if(check_CO) then
        if(i_step==istep_0) then
          !check orthonormality of CO
          allocate(ovlp_check(ns,ns),stat=st)
          if(st/=0) stop "allocation of ovlp_check_matrix!"
          do i=1, ns
            do j=1, ns
              ovlp_check(i,j) = sum(CO(:,i)*CO(:,j))*dV
            end do
          end do
          write(123,*) "i_step:", i_step; call flush(123)
          do i=1, ns
            write(123,*) ovlp_check(i,:); call flush(123)
          end do

          if(allocated(ovlp_check)) deallocate(ovlp_check)
        end if
      end if
end subroutine make_CO_det

subroutine get_core_rest(i_step_gdr)
     use OMP_LIB
     use commvar
     implicit none
     integer :: i,j,k,st
     integer :: ns_sub, na, i_step_gdr
     real*8,allocatable :: locality(:), Q_diag(:,:), AO_rod(:,:)
     integer,allocatable :: label(:)
     integer :: orb_indx_dr
     integer, allocatable :: state_indx2(:)
     
     if(restart) then
       ns_sub = nstates
     else
       ns_sub = ntot_seq
     end if

     if(ns_sub/=ntot_seq) stop "for restarting seq_wan nstates has to be ntot_seq!"

     if(ns_sub.lt.n_c) stop "not enough states to construct core_st matrix!"

     if(restart) then
       do i=1, n_c
         core_st(:,i) = AO(:,i)
       end do

       do i=1, n_rest
         rest_st(:,i) = AO(:,i+n_c)
       end do

     else !not restart
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
       allocate(AO_rod(nx*ny*nz,ns_sub),stat=st)
       if(st/=0) stop "error: allocation of AO_rod matrix"

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

       do i=1, ns_sub
         orb_indx_dr = maxloc(locality,1)
         locality(orb_indx_dr) = -1d0
         AO_rod(:,i) = AO(:,orb_indx_dr)
       end do

       do i=1, n_c
         core_st(:,i) = AO_rod(:,i)
       end do

       do j=1, n_rest
         rest_st(:,j) = AO_rod(:,j+n_c)
       end do

       write(6,*) "done with core_st and rest_st"; call flush(6)

       if(allocated(locality)) deallocate(locality)
       if(allocated(Q_diag)) deallocate(Q_diag)
       if(allocated(label)) deallocate(label)
       if(allocated(AO_rod)) deallocate(AO_rod)
       if(allocated(AO)) deallocate(AO)
     end if
end subroutine get_core_rest

subroutine update_core_rest(i_step_up)
     use commvar
     implicit none
     integer :: i,j,k,st
     integer :: ns_back, i_blk, i_step_up, na
     integer :: indx_start2, indx_end2, orb_indx_dr
     real*8 :: tmp(nx*ny*nz)
     real*8, allocatable :: Q_diag(:,:), rest_rod(:,:), locality(:)

     if(i_step_up.gt.nblock) then
       i_blk = mod(i_step_up,nblock)
     else
       i_blk = i_step_up
     end if

     !update core_st and rest_st based on previous run
     do i = 1, n_c
       core_st(:,i) = LO(:,orb_indx(i))
     end do

     ns_back = n_c + nst_block(i_blk)
     indx_start2 = (i_blk-1)*n_r + 1
     indx_end2 = indx_start2 + nst_block(i_blk) - 1

     j = indx_start2
     do i = 1, ns_back
       if(.not.(any(orb_indx==i))) then
         rest_st(:,j) = LO(:,i)
         j = j + 1
       end if
     end do
     if((j-1)/=indx_end2) stop "error: updating the rest_st matrix!"
     if(allocated(LO)) deallocate(LO)

     write(6,*) "Finish updating core_st and rest_st!"; call flush(6)

     if(mod(i_step_up,nblock)==0) then
       if(reorder_rest) then
         write(6,*) "Reordering rest_st for the next round!"; call flush(6)
         
         if(frag_wan) then
           na = 1 
         elseif(loc_wan) then
           na = noa_lw
         end if

         allocate(Q_diag(na,n_rest),stat=st)
         if(st/=0) stop "error: allocation of Q_diag matrix"
         allocate(rest_rod(nx*ny*nz,n_rest),stat=st)
         if(st/=0) stop "error: allocation of rest_rod matrix"
         allocate(locality(n_rest),stat=st)
         if(st/=0) stop "error: allocation of locality matrix"

         if(frag_wan) then
           !$OMP parallel
           !$OMP do schedule(dynamic)
           do j=1, n_rest
              Q_diag(1,j) = sum(rest_st(:,j)*FWF(:,1)*rest_st(:,j))*dV
           end do
           !$OMP end do
           !$OMP end parallel
         elseif(loc_wan) then
           !$OMP parallel
           !$OMP do schedule(dynamic)
           do i=1, na
             do j=1, n_rest
                Q_diag(i,j) = sum(rest_st(:,j)*AWF(:,atom_label(i))*rest_st(:,j))*dV
             end do
           end do
           !$OMP end do
           !$OMP end parallel
         end if

         write(6,*) "done with Q_diag"; call flush(6)

         if(loc_wan) then
           !$OMP parallel
           !$OMP do schedule(dynamic)
             do i=1, n_rest
               locality(i) = sum(Q_diag(:,i))
             end do
           !$OMP end do
           !$OMP end parallel
         elseif(frag_wan) then
           locality(:) = Q_diag(1,:)
         end if

         do i=1, n_rest
           orb_indx_dr = maxloc(locality,1)
           locality(orb_indx_dr) = -1d0
           rest_rod(:,i) = rest_st(:,orb_indx_dr)
         end do

         rest_st = rest_rod

         if(allocated(Q_diag)) deallocate(Q_diag)
         if(allocated(rest_rod)) deallocate(rest_rod)
         if(allocated(locality)) deallocate(locality)
       end if
     end if
             
end subroutine update_core_rest
