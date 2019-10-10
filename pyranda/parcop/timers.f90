module LES_timers
  use iso_c_binding

  real(c_double) :: compute_time = 0.0D0
  real(c_double) :: comm_time = 0.0D0
  real(c_double) :: total_time = 0.0D0

  real(c_double) :: compute_time_start
  real(c_double) :: comm_time_start
  real(c_double) :: total_time_start
  

contains
  
  subroutine startCPU()
    call cpu_time(compute_time_start)    
  end subroutine startCPU


  subroutine endCPU()
    real(c_double) :: finish
    call cpu_time(finish)
    compute_time = finish - compute_time_start    
  end subroutine endCPU



  subroutine startCOMM()
    call cpu_time(comm_time_start) 
  end subroutine startCOMM

  subroutine endCOMM()
    real(c_double) :: finish
    call cpu_time(finish)
    comm_time = finish - comm_time_start  + comm_time
  end subroutine endCOMM


end module LES_timers
