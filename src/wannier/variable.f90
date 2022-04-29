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
module commvar
      implicit none
      save
      real*8,allocatable::CO(:,:),LO(:,:),evls(:),AO(:,:)
      real*8,allocatable::Q_A0(:,:,:),Q_A1(:,:,:)
      real*8,allocatable::AWF(:,:),GR(:,:),Nel(:)
      real*8,allocatable::COA(:,:),DOA(:,:)
      real*8,allocatable::A(:,:),U(:,:)
      integer,allocatable::Nvel(:)
      integer,allocatable::state_indx(:)
      integer:: nstates, n_occ, n_unocc
      integer:: noa
      integer:: nx,ny,nz,nsp
      integer:: nthreads
      integer:: max_iter
      integer:: ns
      real*8,parameter::gamma_a=0.944863,pi=4*atan(1.0_8)
      real*8::sqrt_2pi=dsqrt(2d0*pi),delta_t,corr_fac, rcut
      real*8::dV,dx,dy,dz
      real*8::start_all,finish_all
      real*8::con_thres_inner
      real*8::evwd, loc_thres
      logical::restart
      character(16) :: space, space_size, def_sub
      character(50) :: fname1, fname2, fname3
      logical::occ, unocc, random
      logical::full, sub
      logical::n_state, loc_cut, energy
      logical::continuous
      character(16)::atom_space, exhaust_space
      logical::frag_wan, full_atom, loc_wan
      logical::sto_wan, seq_wan, simp_wan
      !loc_wan
      integer :: noa_lw, ns_lw, sfw
      integer,allocatable::atom_label(:),orb_indx(:)
      real*8,allocatable::orb_con(:),forb_con(:)
      !frag_wan
      integer :: noa_fw, n_frag
      real*8,allocatable::FWF(:,:)
      integer,allocatable::frag_label(:)
      integer,allocatable::ns_fw(:)    
      !sequential wannierization
      logical :: check_CO
      integer :: ntot_seq, ntot_sto, n_tot2, nstates2
      integer :: step_prnt_info, step_prnt_tmp, nhit_cong_in, iter_acc, iter_eff_acc
      real*8 :: con_thres_outer, OBJ_m0, OBJ_m1
      !stochastic wannierization
      integer :: n_det, n_sto
      integer :: mc_step, i_mc_step0, nhit_cong_out
      real*8, allocatable :: AO2(:,:), det_st(:,:), LO_sv(:,:)
      !deterministic sequential
      integer :: n_c, n_r
      integer :: nblock, n_rest, total_step, istep_0
      integer, allocatable :: nst_block(:)
      real*8, allocatable :: rest_st(:,:), core_st(:,:)
      logical :: reorder_rest
      real*8 :: d_OBJ_max
end module commvar
