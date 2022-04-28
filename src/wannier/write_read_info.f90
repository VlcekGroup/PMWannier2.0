subroutine write_read_info
      use commvar
      implicit none

      write(6,*) '--------------------------------------------------------------'; call flush(6)
      if(restart) then
        write(6,*) 'Wannierization restarts from previous wannier functions.'; call flush(6)
      else
        write(6,*) 'Wannierization starts from Kohn-Sham eigen states.'; call flush(6)
      endif

      write(6,*) '--------------------------------------------------------------'; call flush(6)
      if(loc_wan) then
        write(6,*) 'Performing Wannierization on selected atoms.'; call flush(6)
      elseif(frag_wan) then
        write(6,*) 'Performing Wannierization on selected fragments.'; call flush(6)
      else
        write(6,*) 'Performing Wannierization on all atoms.'; call flush(6)
      end if

      write(6,*) '--------------------------------------------------------------'; call flush(6)
      if(occ) then
        write(6,*) 'Performing Wannierization in the occupied space.'; call flush(6)
      elseif(unocc) then
        write(6,*) 'Performing Wannierization in the unoccupied states.'; call flush(6)
      elseif(random) then
        write(6,*) 'Performing Wannierization of random specified state.'; call flush(6)
      end if

      if(seq_wan.or.sto_wan) then
        write(6,*) 'Performing sequential Wannierization of specified state.'; call flush(6)
      end if

end subroutine write_read_info
