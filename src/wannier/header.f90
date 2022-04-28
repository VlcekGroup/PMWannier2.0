subroutine write_header
  implicit none
     write(6,*)' ************************************************************** '
     write(6,*)
     write(6,*)' Pipek-Mezey Wannier Code, v3.0 (Feb/2022) '
     write(6,*)
     write(6,*)' ************************************************************** '
     write(6,*)
     call flush(6)
end subroutine write_header
