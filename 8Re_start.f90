subroutine binresult
use constant
use variable
use ANCF
Use Local
implicit none
character(12):: fn
integer:: i,j,k,num
real(8):: f0
integer,allocatable:: i4a(:)
allocate(i4a(ne_orig))
i4a=0
write(fn,'(i8.8,a)') nstep, '.b'
open(unit=100,file=fn,form='unformatted')
write(100) t,dt,x,y,u0,v0,pn,nstep,x0_ale,y0_ale,um_ale,vm_ale,DELTX,DELTY,acc0_ang,vel0_ang,dis0_ang,acc_ang,vel_ang,dis_ang,&
           & ACC0_Y,VEL0_Y,DIS0_Y,ACC_Y,VEL_Y,DIS_Y,ACC0_X,VEL0_X,DIS0_X,ACC_X,VEL_X,DIS_X,tor0,tor,CL0,CL,CD0,CD,CL1,CD1,CL2, & 
           & CD2,ACC0_Y_UP,VEL0_Y_UP,DIS0_Y_UP,ACC_Y_UP,VEL_Y_UP,DIS_Y_UP,AF_e0, AF_e,AF_ed0, AF_ed,AF_edd0, AF_edd!,Delta_x510,delta_y510,Delta_x0510,delta_y0510

end subroutine binresult
!
!
!
!
subroutine restart
use constant
use variable
use ANCF
Use Local
implicit none
integer:: i,j,k,nn_read,ne_read,nff
logical:: l
real(8):: f0,temp
integer,allocatable:: ele_new_old_read(:)
character(12):: fn,fn1,add
!open(unit=1000,file='nstep.txt')
!read(1000,*)fn
!close(1000)
!write(*,*)fn

open(unit=101,file=file_name,form='unformatted')
read(101) t,dt,x,y,u0,v0,pn,nstep,x0_ale,y0_ale,um_ale,vm_ale,DELTX,DELTY,acc0_ang,vel0_ang,dis0_ang,acc_ang,vel_ang,dis_ang,&
           & ACC0_Y,VEL0_Y,DIS0_Y,ACC_Y,VEL_Y,DIS_Y,ACC0_X,VEL0_X,DIS0_X,ACC_X,VEL_X,DIS_X,tor0,tor,CL0,CL,CD0,CD,CL1,CD1,CL2, & 
           & CD2,ACC0_Y_UP,VEL0_Y_UP,DIS0_Y_UP,ACC_Y_UP,VEL_Y_UP,DIS_Y_UP,AF_e0, AF_e,AF_ed0, AF_ed,AF_edd0, AF_edd
write(*,*)
write(*,*) '********** congratulations ! **********'
write(*,*) 
write(*,*) '********** restaring ! **********'
write(*,*) 


end subroutine restart