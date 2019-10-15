module LES_timers
  use iso_c_binding

  real(c_double) :: compute_time = 0.0D0
  real(c_double) :: comm_time = 0.0D0
  real(c_double) :: total_time = 0.0D0

  real(c_double) :: compute_time_start
  real(c_double) :: comm_time_start
  real(c_double) :: total_time_start
  
  real(c_double),dimension(10) :: custom_time_start
  real(c_double),dimension(10) :: custom_time = 0.0

  !custom_time = 0.0
  
contains
  
  subroutine startCPU()
    call cpu_time(compute_time_start)    
  end subroutine startCPU


  subroutine endCPU()
    real(c_double) :: finish
    call cpu_time(finish)
    compute_time = finish - compute_time_start + compute_time
  end subroutine endCPU



  subroutine startCOMM()
    call cpu_time(comm_time_start) 
  end subroutine startCOMM

  subroutine endCOMM()
    real(c_double) :: finish
    call cpu_time(finish)
    comm_time = finish - comm_time_start  + comm_time
  end subroutine endCOMM


  subroutine startCUSTOM(index)
    integer, intent(in) :: index
    call cpu_time(custom_time_start(index)) 
  end subroutine startCUSTOM

  subroutine endCUSTOM(index)
    integer, intent(in) :: index
    real(c_double) :: finish
    call cpu_time(finish)
    custom_time(index) = finish - custom_time_start(index)  + custom_time(index)
  end subroutine endCUSTOM
  

end module LES_timers
