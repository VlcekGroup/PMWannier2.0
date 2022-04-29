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
subroutine read_input
      use commvar
      implicit none
      character(100) :: rch
      integer :: inpstat, st, i
      integer :: first_indx, last_indx


      write(6,*) '--------------------------------------------------------------'; call flush(6)
      write(6,*) "Start reading input file."; call flush(6)
      
      !general
      restart = .false.
      nthreads = 1
      max_iter = 500
      delta_t = 5.0
      corr_fac = 1.1
      con_thres_inner = 1E-7
      nhit_cong_in = 3
      n_unocc = -1
      rcut = 7.18096
      space = 'occ'
      occ = .true.
      unocc = .false.
      random = .false.
      continuous = .false.
      space_size = 'full'
      full = .true.
      sub = .false.
      def_sub = 'loc_cut'
      loc_cut = .true.
      loc_thres = 1E-3
      energy = .false.
      evwd = -1d0
      n_state = .false.
      atom_space = 'full'
      full_atom = .true.
      !loc_wan
      loc_wan = .false.
      noa_lw = -1
      ns_lw = -1
      sfw = 50
      !frag_wan
      frag_wan = .false.
      n_frag = -1
      noa_fw = -1
      exhaust_space = 'simple'
      sto_wan = .false.
      seq_wan = .false.
      simp_wan = .true.
      !sequential wannierzation
      check_CO = .false.
      ntot_seq = -1
      ntot_sto = -1
      n_tot2 = -1
      nstates2 = -1
      step_prnt_info = 50
      step_prnt_tmp = 100
      iter_acc = 0
      iter_eff_acc = 0
      con_thres_outer = 5E-7
      OBJ_m0 = 0d0
      OBJ_m1 = 0d0
      !stochastic search
      n_det = -1
      n_sto = 1
      mc_step = 1000
      i_mc_step0 = 1
      nhit_cong_out = 5
      !deterministic search
      total_step = 1000
      n_r = 1
      n_c = -1 
      istep_0 = 1
      reorder_rest = .true.
      d_OBJ_max = -1d0

      open(111,file='input',form='formatted',status='old',action='read')
      do
        read(111,*,iostat=inpstat) rch
        if(inpstat/=0) exit
        select case(rch)
        case('restart')
          read(111,*,iostat=inpstat) restart
          if(inpstat/=0) stop "restart value missing in input"
        case('max_iter')
          read(111,*,iostat=inpstat) max_iter
          if(inpstat/=0) stop "max_iter value missing in input"
        case('delta_t')
          read(111,*,iostat=inpstat) delta_t
          if(inpstat/=0) stop "delta_t value missing in input"
        case('delta_t_correction')
          read(111,*,iostat=inpstat) corr_fac
          if(inpstat/=0) stop "corr_fac value missing in input"
        case('convergence_threshold_inner')
          read(111,*,iostat=inpstat) con_thres_inner
          if(inpstat/=0) stop "con_thres_inner value missing in input"
        case('nhit_convergence_inner')
          read(111,*,iostat=inpstat) nhit_cong_in
          if(inpstat/=0) stop "nhit_cong_in value missing in input" 
        case('num_of_threads')
          read(111,*,iostat=inpstat) nthreads
          if(inpstat/=0) stop "nthreads value missing in input"
        case('atom_space')
          read(111,*,iostat=inpstat) atom_space
          if(inpstat/=0) stop "atom_space value missing input"
          if(atom_space=='local') then
            full_atom = .false.
            loc_wan = .true.
          elseif(atom_space=='fragment') then
            full_atom = .false.
            frag_wan = .true.
          elseif(atom_space=='full') then
            full_atom = .true.
          else
            stop "error: unrecognized atom space type!"
          end if
        case('space')
          read(111,*,iostat=inpstat) space
          if(inpstat/=0) stop "space value missing in input"
          if(space=='occ') then
            occ = .true.
          elseif(space=='unocc') then
            occ = .false.
            unocc = .true.
          elseif(space=='random') then
            occ = .false.
            random = .true.
          else
            stop "error: unrecognized space type in input"
          end if
          if(random) then
            read(111,*,iostat=inpstat) ns
            if(inpstat/=0) stop "ns value missing in input"
            allocate(state_indx(ns),stat=st)
            if(st/=0) stop "error: allocation of the state_indx array"
            read(111,*,iostat=inpstat) continuous
            if(inpstat/=0) stop "continuous value missing in input"
            if(continuous) then
              read(111,*,iostat=inpstat) first_indx, last_indx
              if(inpstat/=0) stop "first and last indx value missing in input"
              do i=1, ns
                state_indx(i) = first_indx
                first_indx = first_indx + 1
              end do
              if((first_indx-1)/=last_indx) stop "error in reading continuous state indx"
            else
              read(111,*,iostat=inpstat) state_indx(:)
              if(inpstat/=0) stop "state_indx value missing in input"
            end if
          end if
        case('space_size')
          read(111,*,iostat=inpstat) space_size
          if(inpstat/=0) stop "space_size value missing in input"
          if(space_size=='full') then
            full = .true.
          elseif(space_size=='sub') then
            full = .false.
            sub = .true.
          else
            stop "error: unrecognized space size!"
          end if
        case('define_subspace')
          read(111,*,iostat=inpstat) def_sub
          if(inpstat/=0) stop "def_sub value missing in input"
          if(def_sub=='loc_cut') then
            loc_cut = .true.
            read(111,*,iostat=inpstat) loc_thres
            if(inpstat/=0) stop "loc_thres value missing in input"
          elseif(def_sub=='nstates') then
            loc_cut = .false.
            n_state = .true.
            read(111,*,iostat=inpstat) ns
            if(inpstat/=0) stop "ns value missing in input"
          elseif(def_sub=='energy') then
             loc_cut = .false.
             energy = .false.
             read(111,*,iostat=inpstat) evwd
             if(inpstat/=0) stop "evwd value missing in input"
           else
              stop "error: unrecognized def_sub variable"
           end if
        case('local_wannierization')
          read(111,*,iostat=inpstat) ns_lw
          if(inpstat/=0) stop "ns_lw value missing in input"
          read(111,*,iostat=inpstat) noa_lw
          if(inpstat/=0) stop "noa_lw value missing in input"
          allocate(atom_label(noa_lw),stat=st)
          if(st/=0) stop "error: allocation of atom_label matrix!"
          read(111,*,iostat=inpstat) atom_label(:)
          if(inpstat/=0) stop "atom_label value missing in input"
        case('step_findwf')
          read(111,*,iostat=inpstat) sfw
          if(inpstat/=0) stop "sfw value missing in input"
        case('fragment_wannierization')
          read(111,*,iostat=inpstat) n_frag
          if(inpstat/=0) stop "n_frag value missing in input"
        case('exhaust_space')
          read(111,*,iostat=inpstat) exhaust_space
          if(inpstat/=0) stop "exhaust space value missing in input"
          if(exhaust_space=='sequential') then
            simp_wan = .false.
            seq_wan = .true.
          elseif(exhaust_space=='stochastic') then
            simp_wan = .false.
            sto_wan = .true.
          elseif(exhaust_space=='simple') then
            simp_wan = .true.
          else
            stop "error: unrecognized orbital space search"
          end if
        case('convergence_threshold_outer')
          read(111,*,iostat=inpstat) con_thres_outer
          if(inpstat/=0) stop "con_thres_outer value missing in input"
        case('nhit_cong_outer')
          read(111,*,iostat=inpstat) nhit_cong_out
          if(inpstat/=0) stop "nhit_cong_outer value missing in input"
        case('step_prnt_info')
          read(111,*,iostat=inpstat) step_prnt_info
          if(inpstat/=0) stop "step_prnt_info value missing in input"
        case('step_prnt_tmp')
          read(111,*,iostat=inpstat) step_prnt_tmp
          if(inpstat/=0) stop "step_prnt_tmp value missing in input"
        case('check_CO')
          read(111,*,iostat=inpstat) check_CO
          if(inpstat/=0) stop "check_CO value missing in input"
        case('nstates_seqwan')
          read(111,*,iostat=inpstat) ntot_seq
          if(inpstat/=0) stop "nstates_seqwan value missing in inpt"
        case('nstates_core')
          read(111,*,iostat=inpstat) n_c
          if(inpstat/=0) stop "n_c value missing in input"
        case('nstates_r')
          read(111,*,iostat=inpstat) n_r
          if(inpstat/=0) stop "n_r value missing in input"
        case('total_step')
          read(111,*,iostat=inpstat) total_step
          if(inpstat/=0) stop "total_step value missing in input"
        case('istep_0')
          read(111,*,iostat=inpstat) istep_0
          if(inpstat/=0) stop "i_mc_step value missing in input"
        case('reorder_rest')
          read(111,*,iostat=inpstat) reorder_rest
          if(inpstat/=0) stop "reorder_rest value missing in input"
        case('nstates_stowan')
          read(111,*,iostat=inpstat) ntot_sto
          if(inpstat/=0) stop "ntot_sto value missing in input"
        case('nstates_det')
          read(111,*,iostat=inpstat) n_det
          if(inpstat/=0) stop "n_det value missing in input"
        case('nstates_sto')
          read(111,*,iostat=inpstat) n_sto
          if(inpstat/=0) stop "n_sto value missing in input"
        case('mc_step')
          read(111,*,iostat=inpstat) mc_step
          if(inpstat/=0) stop "mc_step value missing in input"
        case('i_mc_step0')
          read(111,*,iostat=inpstat) i_mc_step0
          if(inpstat/=0) stop "i_mc_step0 missing in input"
        case default
          write(6,*) "unrecognized input variable", rch
          stop
        end select
      end do

      close(111)

 end subroutine read_input

 subroutine read_cnt
     use commvar
     implicit none
     character(50)::rch,ccha
     integer::i,st,inpstat
     integer::cha, zcha, nvele

     write(6,*) '--------------------------------------------------------------'; call flush(6)
     write(6,*) 'Start reading cnt.ini file.'; call flush(6)

     i=0
     open(111,file='cnt.ini',form='formatted',status='old',action='read')
       do
         read(111,*,iostat=inpstat) rch
         if(inpstat/=0) exit
         i = i+1
       enddo

       noa = i

       rewind(111)

       allocate(COA(noa,3),stat=st)
       if(st/=0) stop 'error: allcation of atomic coordinate matrix'
       allocate(Nel(noa),stat=st)
       if(st/=0) stop 'error: allocation of atomic effective electron array'
       allocate(Nvel(noa),stat=st)

       do i=1,noa
         read(111,*) ccha,COA(i,:)
         if(ccha=='Si') then
          cha = 14
          zcha = 4
         else if(ccha=='P') then
            cha = 15
            zcha = 5
         else if(ccha=='H ') then
            cha = 1
            zcha = 1
         else if(ccha=='C ') then
            cha = 6
            zcha = 4
         else if(ccha=='B ') then
            cha = 5
            zcha = 3
         elseif(ccha=='Li') then
            cha = 3
            zcha = 1
         else if(ccha=='N ') then
            cha = 7
            zcha = 5
         else if(ccha=='O ') then
            cha = 8
            zcha = 6
         else if(ccha=='F ') then
            cha = 9
            zcha = 7
         else if(ccha=='Zn') then
            cha = 30
            zcha = 12
         else if(ccha=='Ni') then
            cha = 28
            zcha = 10
         else if(ccha=='W ') then
            cha = 74
            zcha = 6
         else if(ccha=='Au') then
            cha = 79
            zcha = 11
         else if(ccha=='S ') then
            cha = 16
            zcha = 6
         else if(ccha=='He') then
            cha = 2
            zcha = 2
         else if(ccha=='Se') then
            cha = 34
            zcha = 6
         else if(ccha=='Ar') then
            cha = 18
            zcha = 8
         else if(ccha=='As') then
            cha = 33
            zcha = 5
         else if(ccha=='Cu') then
            cha = 29
            zcha = 11
         else if(ccha=='Cl') then
            cha = 17
            zcha = 7
         else if(ccha=='Br') then
            cha = 35
            zcha = 7
         else if(ccha=='I ') then
            cha = 53
            zcha = 7
         else if(ccha=='Be') then
            cha = 4
            zcha = 2
         else if(ccha=='Al') then
            cha = 13
            zcha = 3
         else if(ccha=='Rb') then
            cha = 37
            zcha = 1
         else if(ccha=='Na') then
            cha = 11
            zcha = 1
         else if(ccha=='K ') then
            cha = 19
            zcha = 1
         else if(ccha=='Ga') then
            cha = 31
            zcha = 3
         else if(ccha=='Ge') then
            cha = 32
            zcha = 4
         else if(ccha=='Kr') then
            cha = 36
            zcha = 8
         else if(ccha=='Ne') then
            cha = 10
            zcha = 8
         else if(ccha=='Ag') then
            cha = 47
            zcha = 11
         else if(ccha=='Ti') then
            cha = 22
            zcha = 4
         else if(ccha=='Xe') then
            cha = 54
            zcha = 8d0
         else if(ccha=='Mg') then
            cha = 12
            zcha = 2
         else
            write(6,*)' stopping, since we need to include also the case of atom ',ccha
            stop
         end if
         Nel(i) = cha
         Nvel(i) = zcha
     end do
     close(111)

     nvele = 0

     do i=1,noa
       nvele = nvele + Nvel(i)
     end do

     n_occ = nvele/2

     write(6,*) 'Number of atoms in total:', noa; call flush(6)
     write(6,*) 'Number of occupied states in total:', n_occ; call flush(6)
       
end subroutine read_cnt

subroutine read_bin
        use commvar
        implicit none
        integer::i,st
        character(9) :: ch

          write(6,*) '--------------------------------------------------------------'; call flush(6)
          write(6,*) 'Start reading bin file.'; call flush(6)

          open(111,file=trim(fname1),form='unformatted',status='old',action='read')
          read(111) ch,nx
          read(111) ch,ny
          read(111) ch,nz
          read(111) ch,dx
          read(111) ch,dy
          read(111) ch,dz
          read(111) ch,nsp
          read(111) ch,nstates
          write(6,*) 'Number of states in total:', nstates; call flush(6)
          allocate(evls(nstates),stat=st)
          if(st/=0) stop 'error:allocation of eigenvalue array'
          read(111) !skip a line
          read(111) evls(:)
          read(111) !skip a line

          if(unocc) then
            if(restart) then
              n_unocc = nstates
            else
              n_unocc = nstates - n_occ
            endif
          write(6,*) 'Number of unoccupied states in total:', n_unocc; call flush(6)
          endif
          
          write(6,*) 'Number of grid points:', nx*ny*nz; call flush(6)
          allocate(AO(nx*ny*nz,nstates),stat=st)
          if(st/=0) stop 'error:allocation of all orbital matrix'
            
          do i=1,nstates
            read(111) !skip a line
            read(111) AO(:,i)
          end do 

          close(111)

        dV = dx*dy*dz
end subroutine read_bin

subroutine read_AO2
        use commvar
        implicit none
        integer::i,st
        character(9) :: ch

          write(6,*) '--------------------------------------------------------------'; call flush(6)
          write(6,*) 'Start reading wf.bin file.'; call flush(6)

          open(111,file='wf.bin',form='unformatted',status='old',action='read')
          read(111) ch,nx
          read(111) ch,ny
          read(111) ch,nz
          read(111) ch,dx
          read(111) ch,dy
          read(111) ch,dz
          read(111) ch,nsp
          read(111) ch,nstates2
          if(allocated(evls)) deallocate(evls)
          allocate(evls(nstates2),stat=st)
          if(st/=0) stop 'error:allocation of eigenvalue array'
          read(111) !skip a line
          read(111) evls(:)
          read(111) !skip a line

          allocate(AO2(nx*ny*nz,nstates2),stat=st)
          if(st/=0) stop 'error:allocation of all orbital matrix'

          do i=1,nstates2
            read(111) !skip a line
            read(111) AO2(:,i)
          end do

          close(111)

end subroutine read_AO2
