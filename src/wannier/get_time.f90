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
subroutine print_start_time

  implicit none
  integer :: values(8)
  integer :: date

  call date_and_time(VALUES=values)

  date = values(2)*10**8 + values(3)*10**6 + values(5)*10**4 + values(6)*10**2 + values(7)

  write(6,*) "Job starts at the following time:", date; call flush(6)

end subroutine print_start_time

subroutine print_end_time

  implicit none
  integer :: values(8)
  integer :: date

  call date_and_time(VALUES=values)
  write(6,*) "Finish call date and time!"; flush(6)

  date = values(2)*10**8 + values(3)*10**6 + values(5)*10**4 + values(6)*10**2 + values(7)

  write(6,*) "Job ends at the following time:", date; call flush(6)

      
end subroutine print_end_time

