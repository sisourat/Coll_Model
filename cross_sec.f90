program bProb
  implicit none

  real(kind=8), allocatable :: para_b(:) ! parameter
  real(kind=8), allocatable :: Prob_state(:,:) ! probability
  real(kind=8), allocatable :: eigval(:) ! energy
  real(kind=8), allocatable :: sigma(:) ! cross section of each state( upper limit)
  integer                   :: num_bound,    num_single,   num_double! number of Exc, S, D
  real(kind=8)              :: sum_bound,    sum_single,   sum_double! tot cross sections of Exc, S, D( upper limit)
  real(kind=8)              :: P_of_b_ex,    P_of_b_SI,    P_of_b_DI ! (num_b, tot p)  probability of each b
  real(kind=8)              :: P_of_b_ex_dy, P_of_b_SI_dy, P_of_b_DI_dy  ! (num_b, tot p)  probability of each b
  real(kind=8), allocatable :: sum_norm_Dyson(:) ! SND 
  real(kind=8), allocatable :: norm_Dyson(:,:) ! each P_kj in Dyson_norms.txt
  real(kind=8)              :: db, normtot
  real(kind=8)              :: sum_bound_new, sum_single_new, sum_double_new !!cross section (Dyson)
  ! tot cross sections of each, using the Dyson model

  real(kind=8)              :: thresh_bou_sin, thresh_sin_dou !to sort the energy

  integer                   :: i, j, k
  integer                   :: num_b, num_state, cation_state, dyson_num! the number of b, states(he), states(he+) 

  !real(kind=8), parameter   :: db=0.2D0 ! depend on the coll_input.
  !!! you can also use 
  real(kind=8), parameter   :: pi=DACOS(-1.0D0)

  thresh_bou_sin= -2.0D0
  thresh_sin_dou=  0.0D0

  open(unit=110,file="Prop_collision.out",status="old")
  read(110,*) num_b, num_state

  allocate(para_b(1:num_b))
  allocate(Prob_state(1:num_b,1:num_state))
  allocate(eigval(1:num_state))
  allocate(sigma(1:num_state))

  do i=1,num_state
    read(110,*) eigval(i)
  end do
  do i=1,num_b
    read(110,*) para_b(i), Prob_state(i,1:num_state)
    write(*,*) para_b(i), Prob_state(i,1:num_state)
  end do
  close(110)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(unit=130,file="cross_section_up.dat",status="replace")
  ! cross section integral.
  do i=1,num_state
    sigma(i)=0.0D0
    do j=1,num_b-1
      db=para_b(j+1)-para_b(j)
      sigma(i)=sigma(i) + (para_b(j)*Prob_state(j,i)+para_b(j+1)*Prob_state(j+1,i))*db/2.0D0
    end do
    sigma(i)=2.0D0*pi*sigma(i)
    write(130,"(I5,F18.12,2F18.10)") i,eigval(i),sigma(i),sigma(i)/3.57D0
    write(*,"(I5,F18.12,2F18.10)") i,eigval(i),sigma(i),sigma(i)/3.57D0
  end do
  close(130)

end program
