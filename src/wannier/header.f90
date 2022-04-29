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
subroutine write_header
  implicit none
     write(6,*)' ************************************************************** '
     write(6,*)
     write(6,*)' Pipek-Mezey Wannier Code, v2.0 (Feb/2022) '
     write(6,*)
     write(6,*)' ************************************************************** '
     write(6,*)
     call flush(6)
end subroutine write_header
