subroutine mat_alloc1
        use commvar
        implicit none
        integer::st

        allocate(DOA(nx*ny*nz,noa),stat=st)
        if(st/=0) stop 'error: allocation of density of atom matrix'
        allocate(GR(nx*ny*nz,3),stat=st)
        if(st/=0) stop 'error: allocation of grid matrix'
        allocate(AWF(nx*ny*nz,noa),stat=st)
        if(st/=0) stop 'error: allocation of atomic weight function matrix'
end subroutine mat_alloc1

subroutine mat_alloc2(i_frag)
        use commvar
        implicit none
        integer :: st, i_frag, n_sub

        if(sto_wan) then
          n_sub = n_det
        else
          n_sub = n_c
        end if

        !write(6,*) "Squre of ns:", ns**2; call flush(6)
        allocate(A(ns,ns),stat=st)
        if(st/=0) stop 'error:allocation of A matrix'
        allocate(U(ns,ns),stat=st)
        if(st/=0) stop 'error: allocation of U matrix'
        allocate(orb_con(ns),stat=st)
        if(st/=0) stop 'error:allocation in orb_con matrix'
        if(.not.allocated(LO)) then
          allocate(LO(nx*ny*nz,ns),stat=st)
          if(st/=0) stop 'error: allocation of localized orbital matrix'
        end if


        if(seq_wan.or.sto_wan) then
          if(.not.allocated(orb_indx)) then
            allocate(orb_indx(n_sub),stat=st)
            if(st/=0) stop 'error: allocation in orb_indx array!'
          endif

          if(.not.allocated(forb_con)) then
            allocate(forb_con(n_sub),stat=st)
            if(st/=0) stop 'error: allocation in forb_con array1'
          end if

          if(sto_wan) then
            if(.not.allocated(LO_sv)) then
              allocate(LO_sv(nx*ny*nz,ns),stat=st)
              if(st/=0) stop 'error: allocation of localized orbital matrix!'
            end if
          end if

          if(loc_wan) then
            allocate(Q_A0(noa_lw,ns,ns),stat=st)
            if(st/=0) stop 'error:allocation in Q_A0 matrix'
            allocate(Q_A1(noa_lw,ns,ns),stat=st)
            if(st/=0) stop 'error:allocation in Q_A1 matrix'
          end if
        else !not seq_wan or sto_wan
          if(loc_wan) then
            allocate(Q_A0(noa_lw,ns,ns),stat=st)
            if(st/=0) stop 'error:allocation in Q_A0 matrix'
            allocate(Q_A1(noa_lw,ns,ns),stat=st)
            if(st/=0) stop 'error:allocation in Q_A1 matrix'
            allocate(orb_indx(ns_lw),stat=st)
            if(st/=0) stop 'error:allocation in orb_indx matrix'
            allocate(forb_con(ns_lw),stat=st)
            if(st/=0) stop 'error:allocation in found orb_con matrix'
          elseif(frag_wan) then
            allocate(orb_indx(ns_fw(i_frag)),stat=st)
            if(st/=0) stop 'error:allocation in orb_indx matrix'
            allocate(forb_con(ns_fw(i_frag)),stat=st)
            if(st/=0) stop 'error:allocation in found orb_con matrix'
          else
            allocate(Q_A0(noa,ns,ns),stat=st)
            if(st/=0) stop 'error:allocation in Q_A0 matrix'
            allocate(Q_A1(noa,ns,ns),stat=st)
            if(st/=0) stop 'error:allocation in Q_A1 matrix'
          end if
        end if
end subroutine mat_alloc2

subroutine mat_dealloc1
      use commvar, only: COA,GR,Nel,DOA
      implicit none

      If(allocated(COA)) Deallocate(COA)
      If(allocated(GR)) Deallocate(GR)
      If(allocated(Nel)) Deallocate(Nel)
      If(allocated(DOA)) Deallocate(DOA)
end subroutine mat_dealloc1

subroutine mat_dealloc2
      use commvar, only: Q_A0,Q_A1,A,U,orb_con,forb_con,CO

      If(allocated(Q_A0)) Deallocate(Q_A0)
      If(allocated(Q_A1)) Deallocate(Q_A1)
      If(allocated(A)) Deallocate(A)
      If(allocated(U)) Deallocate(U)
      If(allocated(CO)) Deallocate(CO)
      If(allocated(orb_con)) Deallocate(orb_con)
      If(allocated(forb_con)) Deallocate(forb_con)
end subroutine mat_dealloc2

subroutine mat_dealloc3(i_frag)
      use commvar
      implicit none
      integer :: i_frag

      If((.not.frag_wan).or.(i_frag==n_frag)) then
        If(allocated(evls)) Deallocate(evls)
      end if
      If(allocated(LO)) Deallocate(LO)
      If(allocated(LO_sv)) deallocate(LO_sv)
      If(allocated(orb_indx)) Deallocate(orb_indx)
end subroutine mat_dealloc3
