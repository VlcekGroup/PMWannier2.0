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
subroutine make_CO(i_frag)
      use commvar
      implicit none
      integer :: i, st, n_rd, i_frag
      real*8 :: emax

      if(occ) then
        fname2 = 'wannier_occ_tmp.bin'
        fname3 = 'wannier_occ.bin'
        if(full) then
          if(restart) then
            ns = nstates
          else
            ns = n_occ
          endif
          allocate(CO(nx*ny*nz,ns),stat=st)
          if(st/=0) stop "error: allocation of CO matrix"
          do i=1, ns
            CO(:,i) = AO(:,i)
          end do
        elseif(sub) then
          call get_indx_occ(i_frag)
          write(6,*) 'done getting state index'; call flush(6)
          allocate(CO(nx*ny*nz,ns),stat=st)
          if(st/=0) stop "error: allocation of CO matrix"
          do i=1, ns
            CO(:,i) = AO(:,state_indx(i))
          end do
        end if
      elseif(unocc) then
        fname2 = 'wannier_unocc_tmp.bin'
        fname3 = 'wanner_unocc.bin'
        if(full) then
          ns = n_unocc
          allocate(CO(nx*ny*nz,ns),stat=st)
          if(st/=0) stop "error: allocation of CO matrix"
          if(restart) then
            CO = AO
          else
            do i=1, ns
              CO(:,i) = AO(:,i+n_occ)
            end do
          end if
        elseif(sub) then
          if(energy) then
            if(restart) stop "energy window subspace unavailable for restart"
            emax = evls(n_occ) + evwd
            n_rd = 0
            do i=1,n_unocc
              if(evls(i+n_occ).lt.emax) then
                n_rd = n_rd + 1
              else 
                exit
              end if
            end do
            ns = n_rd
            write(6,*) "number of states used upon the eigen value:", ns; call flush(6)
            if(ns==0) stop "no unoccupied state is within the eigen value window"
            allocate(CO(nx*ny*nz,ns),stat=st)
            if(st/=0) stop "error: allocation of CO matrix"
            do i=1, ns
              CO(:,i) = AO(:,i+n_occ)
            end do
          else
            call get_indx_unocc(i_frag)
            allocate(CO(nx*ny*nz,ns),stat=st)
            if(st/=0) stop "error: allocation of CO matrix"
            do i=1, ns
              CO(:,i) = AO(:,state_indx(i))
            end do
          end if
        end if
      elseif(random) then
        fname3 = 'wannier.bin'
        fname2 = 'wannier_tmp.bin'
        allocate(CO(nx*ny*nz,ns),stat=st)
        if(st/=0) stop "error: allocation of CO matrix"
        do i=1, ns
          CO(:,i) = AO(:,state_indx(i))
        end do
      end if
      
      if((.not.frag_wan).or.(i_frag==n_frag)) then
        if(allocated(AO)) deallocate(AO)
      end if
      if(allocated(state_indx)) deallocate(state_indx)
end subroutine make_CO

subroutine get_indx_occ(i_frag)
     use OMP_LIB
     use commvar
     implicit none
     integer :: i,j,st,indx,i_frag
     integer :: ns_sub, n_rd, na
     real*8,allocatable :: locality(:), Q_diag(:,:)
     integer,allocatable :: label(:)

     if(restart) then
       ns_sub = nstates
     else
       ns_sub = n_occ
     end if

     if(loc_wan) then
       na = noa_lw
     elseif(frag_wan) then
       na = 1
     else
       na = noa
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
            do i=1, na
              Q_diag(i,j) = sum(AO(:,j)*FWF(:,i_frag)*AO(:,j))*dV
            end do
          end do
       !$OMP end do
       !$OMP end parallel
     else
       !$OMP parallel
       !$OMP do schedule(dynamic)
          do j=1, ns_sub
            do i=1, na
              Q_diag(i,j) = sum(AO(:,j)*AWF(:,label(i))*AO(:,j))*dV
            end do
          end do
       !$OMP end do
       !$OMP end parallel
     endif
     write(6,*) "done with Q_diag"; call flush(6)

     !$OMP parallel
     !$OMP do
       do i=1, ns_sub
         locality(i) = sum(Q_diag(:,i))
       end do
     !$OMP end do
     !$OMP end parallel
     write(6,*) "done with locality"; call flush(6)

     if(n_state) then
       if(ns.gt.ns_sub) stop "error: asking for too many states for occupied subspace calculation"
       write(6,*) "number of states used upon your request:", ns; call flush(6)
       allocate(state_indx(ns),stat=st)
       if(st/=0) stop "error: allocation of state_indx array"
       do i=1, ns
         state_indx(i) = maxloc(locality,1)
         locality(state_indx(i)) = -1d0
       end do
     elseif(loc_cut) then
       n_rd = 0
       do i=1, ns_sub
         if(locality(i).ge.loc_thres) n_rd = n_rd + 1
       end do
       ns = n_rd
       write(6,*) "number of states used upon the locality:", ns; call flush(6)
       allocate(state_indx(ns),stat=st)
       if(st/=0) stop "error: allocation of state_indx array"
       j=1
       do i=1, ns_sub
         if(locality(i).ge.loc_thres) then
           state_indx(j) = i
           j = j+1
         endif
       end do
       if(j/=(ns+1)) stop "error: generating state_indx array"
     elseif(energy) then
       stop "energy window subspace unavailable for occupied space"
     end if

     if(allocated(locality)) deallocate(locality)
     if(allocated(Q_diag)) deallocate(Q_diag)
     if(allocated(label)) deallocate(label)
end subroutine get_indx_occ

subroutine get_indx_unocc(i_frag)
     use OMP_LIB
     use commvar
     implicit none
     integer :: i,j,st,indx, i_frag
     integer :: ns_sub, n_rd, na
     real*8,allocatable :: locality(:), Q_diag(:,:)
     integer,allocatable :: label(:)

     ns_sub = n_unocc

     if(loc_wan) then
       na = noa_lw
     elseif(frag_wan) then
       na = 1
     else
       na = noa
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
     
     if(restart) then
       if(frag_wan) then
         !$OMP parallel
         !$OMP do schedule(dynamic)
           do j=1, ns_sub
             do i=1, na
               Q_diag(i,j) = sum(AO(:,j)*FWF(:,i_frag)*AO(:,j))*dV
             end do
           end do
         !$OMP end do
         !$OMP end parallel
       else
         !$OMP parallel
         !$OMP do schedule(dynamic)
           do j=1, ns_sub
             do i=1, na
               Q_diag(i,j) = sum(AO(:,j)*AWF(:,label(i))*AO(:,j))*dV
             end do
           end do
         !$OMP end do
         !$OMP end parallel
       endif
     else
       if(frag_wan) then
         !$OMP parallel
         !$OMP do schedule(dynamic)
            do j=1, ns_sub
              do i=1, na
                Q_diag(i,j) = sum(AO(:,j+n_occ)*FWF(:,i_frag)*AO(:,j+n_occ))*dV
              end do
            end do
         !$OMP end do
         !$OMP end parallel
       else
         !$OMP parallel
         !$OMP do schedule(dynamic)
            do j=1, ns_sub
              do i=1, na
                Q_diag(i,j) = sum(AO(:,j+n_occ)*AWF(:,label(i))*AO(:,j+n_occ))*dV
              end do
            end do
         !$OMP end do
         !$OMP end parallel
       end if
     end if

     !$OMP parallel
     !$OMP do
       do i=1, ns_sub
         locality(i) = sum(Q_diag(:,i))
       end do
     !$OMP end do
     !$OMP end parallel

     if(n_state) then
       if(ns.gt.ns_sub) stop "error: asking for too many states for unoccupied subspace calculation"
       write(6,*) "number of states used upon your request:", ns; call flush(6)
       allocate(state_indx(ns),stat=st)
       if(st/=0) stop "error: allocation of state_indx array"
       do i=1, ns
         state_indx(i) = maxloc(locality,1)
         locality(state_indx(i)) = -1d0
       end do

       if(.not.restart) then
         state_indx(:) = state_indx(:) + n_occ
       endif

     elseif(loc_cut) then
       n_rd = 0
       do i=1, ns_sub
         if(locality(i).ge.loc_thres) n_rd = n_rd + 1
       end do
       ns = n_rd
       write(6,*) "number of states used upon the locality:", ns; call flush(6)
       allocate(state_indx(ns),stat=st)
       if(st/=0) stop "error: allocation of state_indx array"
       j=1
       do i=1, ns_sub
         if(locality(i).ge.loc_thres) then
           state_indx(j) = i
           j = j+1
         endif
       end do
       if(j/=(ns+1)) stop "error: generating state_indx array"

       if(.not.restart) then
         state_indx(:) = state_indx(:) + n_occ
       end if
     end if

     if(allocated(locality)) deallocate(locality)
     if(allocated(Q_diag)) deallocate(Q_diag)
     if(allocated(label)) deallocate(label)
end subroutine get_indx_unocc
