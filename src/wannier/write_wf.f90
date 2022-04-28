subroutine write_wf_full(i_frag)
  use commvar
  implicit none
  integer::i, i_frag
   
   if(frag_wan) then
     write(fname3,'(a,i3.3,a)') 'wannier',i_frag,'.bin'
   else
     write(fname3,'(a,i3.3,a)') 'wannier.bin'
   end if

   open(222,file=trim(fname3),form='unformatted',status='replace',action='write')
   write(222) 'nx       ', nx
   write(222) 'ny       ', ny
   write(222) 'nz       ', nz
   write(222) 'dx       ', dx
   write(222) 'dy       ', dy
   write(222) 'dz       ', dz
   write(222) 'nsp      ', nsp
   write(222) 'nstates  ', ns
   write(222) 'evls     '
   write(222) evls(:)
   write(222) 'orbitals '
   do i=1,ns
    write(222) i, 1
    write(222) LO(:,i)
   end do
  close(222)
end subroutine write_wf_full

subroutine write_wf_trun(i_frag)
  use commvar
  implicit none
  character(50) :: fname4
  integer::i, ns_trun, i_frag

  if(loc_wan) then
    ns_trun = ns_lw
    write(fname4,'(a,i3.3,a)') 'wannier_trun.bin'
  elseif(frag_wan) then
    ns_trun = ns_fw(i_frag)
    write(fname4,'(a,i3.3,a)') 'wannier',i_frag,'_trun.bin'
  end if
   
   open(222,file=trim(fname4),form='unformatted',status='replace',action='write')
   write(222) 'nx       ', nx
   write(222) 'ny       ', ny
   write(222) 'nz       ', nz
   write(222) 'dx       ', dx
   write(222) 'dy       ', dy
   write(222) 'dz       ', dz
   write(222) 'nsp      ', nsp
   write(222) 'nstates  ', ns_trun
   write(222) 'evls     '
   write(222) evls(:)
   write(222) 'orbitals '
   do i=1,ns_trun
    write(222) i, 1
    write(222) LO(:,orb_indx(i))
   end do
  close(222)
end subroutine write_wf_trun

subroutine write_wf_tmp(i_step_wf)
  use commvar
  implicit none
  character(50) :: fname4
  integer::i, ns_tmp, i_step_wf, i_step_wf2

   i_step_wf2 = i_step_wf
   i_step_wf2 = mod(i_step_wf,1000)
   write(fname4,'(a,i3.3,a)') 'step',i_step_wf2,'.bin'

   open(222,file=trim(fname4),form='unformatted',status='replace',action='write')
   write(222) 'nx       ', nx
   write(222) 'ny       ', ny
   write(222) 'nz       ', nz
   write(222) 'dx       ', dx
   write(222) 'dy       ', dy
   write(222) 'dz       ', dz
   write(222) 'nsp      ', nsp
   write(222) 'nstates  ', ntot_seq
   write(222) 'evls     '
   write(222) evls(:)
   write(222) 'orbitals '
   do i=1, n_c
     write(222) i, 1
     write(222) core_st(:,i)
   end do
   do i=1, n_rest
    write(222) i+n_c, 1
    write(222) rest_st(:,i)
   end do
 close(222)
end subroutine write_wf_tmp

subroutine write_wf_full_seq
  use commvar
  implicit none
  integer :: i
   
   write(fname3,'(a,i3.3,a)') 'wannier.bin'

   open(222,file=trim(fname3),form='unformatted',status='replace',action='write')
   write(222) 'nx       ', nx
   write(222) 'ny       ', ny
   write(222) 'nz       ', nz
   write(222) 'dx       ', dx
   write(222) 'dy       ', dy
   write(222) 'dz       ', dz
   write(222) 'nsp      ', nsp
   write(222) 'nstates  ', ntot_seq
   write(222) 'evls     '
   write(222) evls(:)
   write(222) 'orbitals '
   do i=1, n_c
     write(222) i, 1
     write(222) core_st(:,i)
   end do
   do i=1, n_rest
    write(222) i+n_c, 1
    write(222) rest_st(:,i)
   end do
 close(222)
end subroutine write_wf_full_seq

subroutine write_wf_trun_seq(i_frag)
  use commvar
  implicit none
  character(50) :: fname4
  integer::i, ns_trun, i_frag

  if(loc_wan) then
    write(fname4,'(a,i3.3,a)') 'wannier_trun.bin'
  elseif(frag_wan) then
    write(fname4,'(a,i3.3,a)') 'wannier',i_frag,'_trun.bin'
  end if

  ns_trun = n_c
   
   open(222,file=trim(fname4),form='unformatted',status='replace',action='write')
   write(222) 'nx       ', nx
   write(222) 'ny       ', ny
   write(222) 'nz       ', nz
   write(222) 'dx       ', dx
   write(222) 'dy       ', dy
   write(222) 'dz       ', dz
   write(222) 'nsp      ', nsp
   write(222) 'nstates  ', ns_trun
   write(222) 'evls     '
   write(222) evls(:)
   write(222) 'orbitals '
   do i=1,ns_trun
    write(222) i, 1
    write(222) core_st(:,i)
   end do
  close(222)
end subroutine write_wf_trun_seq

subroutine write_wf_tmp_sto(i_mc)
  use commvar
  implicit none
  character(50) :: fname4
  integer::i, ns_tmp, i_mc, i_mc2

   ns_tmp = ns
   i_mc2 = i_mc
   i_mc2 = mod(i_mc2,1000)
   write(fname4,'(a,i3.3,a)') 'mcstep',i_mc2,'.bin'

   open(222,file=trim(fname4),form='unformatted',status='replace',action='write')
   write(222) 'nx       ', nx
   write(222) 'ny       ', ny
   write(222) 'nz       ', nz
   write(222) 'dx       ', dx
   write(222) 'dy       ', dy
   write(222) 'dz       ', dz
   write(222) 'nsp      ', nsp
   write(222) 'nstates  ', ns_tmp
   write(222) 'evls     '
   write(222) evls(:)
   write(222) 'orbitals '
   do i=1,ns_tmp
    write(222) i, 1
    write(222) LO(:,i)
   end do
  close(222)
end subroutine write_wf_tmp_sto
