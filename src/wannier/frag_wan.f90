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
subroutine frag_wan_prep
      use commvar
      implicit none
      integer :: st, i, j, inpstat
      
      write(6,*) '--------------------------------------------------------------'; call flush(6)
      write(6,*) 'Preparing the fragment weight function with the number of fragment', n_frag; call flush(6)
      if(n_frag.lt.0) stop "not providing the number of fragment for localization!"

        allocate(FWF(nx*ny*nz,n_frag),stat=st)
        if(st/=0) stop "error: allocation of the fragment weight function matrix"
        allocate(ns_fw(n_frag),stat=st)
        if(st/=0) stop "error: allocation of the fragment weight function matrix"
         
        FWF = 0d0
        open(222,file='group_atom',form='formatted',status='old',action='read')
            read(222,*)
          do j=1, n_frag
            read(222,*,iostat=inpstat) ns_fw(j)
            if(inpstat/=0) stop "ns_fw value missing in the group_atom file"
            read(222,*,iostat=inpstat) noa_fw
            if(inpstat/=0) stop "noa_fw value missing in the group_atom file"
            allocate(frag_label(noa_fw),stat=st)
            if(st/=0) stop "error: allocation of the host label matrix"
            read(222,*,iostat=inpstat) frag_label(:)
            if(inpstat/=0) stop "frag_label value missing in the group_atom file" 
            do i=1, noa_fw
              FWF(:,j) = FWF(:,j) + AWF(:,frag_label(i))
            end do
            if(allocated(frag_label)) deallocate(frag_label)
          end do
        close(222)

end subroutine frag_wan_prep
