 
SUBROUTINE GAMBIT_READ    !KS! 改写
USE CONSTANT
USE VARIABLE
IMPLICIT NONE
INTEGER I,J,TEMP
REAL(8)::XL,YU,YD,XR

OPEN(UNIT=15,FILE=MeshFileName)

READ(15,*) NN_ORIG,NE_ORIG
WRITE(*,*) NN_ORIG,NE_ORIG

ALLOCATE(X_ORIG(NN_ORIG))
ALLOCATE(Y_ORIG(NN_ORIG))
ALLOCATE(NODEELE_ORIG(NE_ORIG,4))

DO I=1,NN_ORIG
    READ(15,*) TEMP,X_ORIG(I),Y_ORIG(I)
END DO

DO I=1,NE_ORIG
    READ(15,*) TEMP,(NODEELE_ORIG(I,J),J=1,4)
END DO

CLOSE(15)

OPEN(UNIT=10,FILE='MESH_ORIG.DAT')
WRITE(10,*) 'TITLE="A"'
WRITE(10,*) 'VARIABLES="X","Y" '

WRITE(10,'(1X,A7,I6,A4,I6,A31)')  'ZONE N=' , NN_ORIG, ', E=',NE_ORIG, ',F=FEPOINT, ET=QUADRILATERAL'

WRITE(10,*) 'DT=(SINGLE SINGLE )'

DO I=1,NN_ORIG
    WRITE(10,*) X_ORIG(I),Y_ORIG(I)
END DO

DO I=1,NE_ORIG
    WRITE(10,*) (NODEELE_ORIG(I,J),J=1,4)
END DO
CLOSE(10)

!========================================================================! ALE会用到
NN_ALE=NN_ORIG                                                                  
NE_ALE=NE_ORIG                                                                  
ALLOCATE(X_ALE(NN_ALE),Y_ALE(NN_ALE),NODEELE_ALE(NE_ALE,4),NODETYPE(NN_ALE))    
X_ALE=X_ORIG                                                                    
Y_ALE=Y_ORIG                                                                    
NODEELE_ALE=NODEELE_ORIG                                                        
NODETYPE=0  !在下面会更新出在 A 以及 B区的NODETYPE为 3和2，剩下的不在 A 或者B 的就是0
!========================================================================
DO I=1,NN_ALE

    IF((X_ALE(I)*X_ALE(I)+Y_ALE(I)*Y_ALE(I))>=(spring_in**2+0.02).AND.(X_ALE(I)*X_ALE(I)+Y_ALE(I)*Y_ALE(I))<(spring_out**2+0.02))THEN      ! 论文里面的B区
        NODETYPE(I)=2
    ELSE IF((X_ALE(I)**2+Y_ALE(I)**2)<(spring_in**2+0.02))THEN                                                             ! 论文里面的A区
        NODETYPE(I)=3
    END IF

END DO
    END SUBROUTINE GAMBIT_READ
    
 
    SUBROUTINE Flexible_MeshREAD    !For Ks Case
    USE CONSTANT
    USE VARIABLE
    IMPLICIT NONE

    INTEGER I,J,TEMP
    REAL(8)::XL,YU,YD,XR

    INTEGER::ICINDEX,ISINDEX,JINDEXK
    REAL(8)::XSVOTEXK,YSVOTEXK,XVOTEXK,YVOTEXK,DIS_ALEK,DIS_ALEMINK


    OPEN(UNIT=15,FILE=MeshFileName)

    READ(15,*) NN_ORIG,NE_ORIG
    WRITE(*,*) NN_ORIG,NE_ORIG

    ALLOCATE(X_ORIG(NN_ORIG))
    ALLOCATE(Y_ORIG(NN_ORIG))
    ALLOCATE(NODEELE_ORIG(NE_ORIG,4))

    DO I=1,NN_ORIG
        READ(15,*) TEMP,X_ORIG(I),Y_ORIG(I)
    END DO

    DO I=1,NE_ORIG
        READ(15,*) TEMP,(NODEELE_ORIG(I,J),J=1,4)
    END DO

    CLOSE(15)

    OPEN(UNIT=10,FILE='MESH_ORIG.DAT')
    WRITE(10,*) 'TITLE="A"'
    WRITE(10,*) 'VARIABLES="X","Y" '

    WRITE(10,'(1X,A7,I6,A4,I6,A31)')  'ZONE N=' , NN_ORIG, ', E=',NE_ORIG, ',F=FEPOINT, ET=QUADRILATERAL'

    WRITE(10,*) 'DT=(SINGLE SINGLE )'

    DO I=1,NN_ORIG
        WRITE(10,*) X_ORIG(I),Y_ORIG(I)
    END DO

    DO I=1,NE_ORIG
        WRITE(10,*) (NODEELE_ORIG(I,J),J=1,4)
    END DO
    CLOSE(10)

    NN_ALE=NN_ORIG
    NE_ALE=NE_ORIG
    ALLOCATE(X_ALE(NN_ALE),Y_ALE(NN_ALE),NODEELE_ALE(NE_ALE,4),NODETYPE(NN_ALE))
    X_ALE=X_ORIG
    Y_ALE=Y_ORIG
    NODEELE_ALE=NODEELE_ORIG
    NODETYPE=0  ! Boundary    : 1
                ! First Layer : 78
                ! Do nothing  : 79
                ! Spring Method Mesh Transform domain : 0

    DO I=1,NN_ALE

        IF((X_ALE(I)*X_ALE(I)+Y_ALE(I)*Y_ALE(I))>=(SPRING_OUT**2))THEN     !.or.X_ALE(I)<0.0D0  !.or.(X_ALE(I)<0.3.or.(DABS(Y_ALE(I))>2))
            NODETYPE(I)= 79 ! Do Nothing 
        END IF
    END DO
    
    
    END SUBROUTINE Flexible_MeshREAD
    
    
    SubRoutine FirstMeshLayer_Adj_Structure
    USE CONSTANT
    USE VARIABLE
    IMPLICIT NONE

    INTEGER::ICINDEX,ISINDEX,JINDEXK,I
    REAL(8)::XSVOTEXK,YSVOTEXK,XVOTEXK,YVOTEXK,DIS_ALEK,DIS_ALEMINK,FDIS1,FDIS2,FDIS3,FDIS4,FDIS5

    DO I=1,NN_ALE
        DIS_ALEK=0.0D0
        DIS_ALEMINK=1E6
        
        XSVOTEXK=X_ALE(I)
        YSVOTEXK=Y_ALE(I)

        DO ICINDEX=1,NPBN_ALE
            JINDEXK=PBN_ALE(ICINDEX)
            XVOTEXK=X_ALE(JINDEXK)
            YVOTEXK=Y_ALE(JINDEXK)
            DIS_ALEK=SQRT((XSVOTEXK-XVOTEXK)**2+(YSVOTEXK-YVOTEXK)**2)

            IF(DIS_ALEMINK>=DIS_ALEK)THEN
                DIS_ALEMINK=DIS_ALEK
            END IF
        END DO

        DO ISINDEX=1,NSBN_ALE
            JINDEXK=SBN_ALE(ISINDEX)
            XVOTEXK=X_ALE(JINDEXK)
            YVOTEXK=Y_ALE(JINDEXK)
            DIS_ALEK=SQRT((XSVOTEXK-XVOTEXK)**2+(YSVOTEXK-YVOTEXK)**2)

            IF(DIS_ALEMINK>=DIS_ALEK)THEN
                DIS_ALEMINK=DIS_ALEK
            END IF
        END DO

        FDIS1 = 0.02
        FDIS2 = 0.04
        FDIS3 = 0.08
        FDIS4 = 0.16
        FDIS5 = 0.32
        
        IF((DIS_ALEMINK<=FDIS5).and.(DIS_ALEMINK>FDIS4).and.(DIS_ALEMINK /= 0).and.(NodeType(i)== 0))THEN
            NODETYPE(I)=75 
        ElseIf((DIS_ALEMINK<=FDIS4).and.(DIS_ALEMINK>FDIS3).and.(DIS_ALEMINK /= 0).and.(NodeType(i)== 0)) Then
            NODETYPE(I)=74 
        ElseIf((DIS_ALEMINK<=FDIS3).and.(DIS_ALEMINK>FDIS2).and.(DIS_ALEMINK /= 0).and.(NodeType(i)== 0)) Then
            NODETYPE(I)=73 
        ElseIf((DIS_ALEMINK<=FDIS2).and.(DIS_ALEMINK>FDIS1).and.(DIS_ALEMINK /= 0).and.(NodeType(i)== 0)) Then
            NODETYPE(I)=72 
        ElseIf((DIS_ALEMINK<=FDIS1).and.(DIS_ALEMINK /= 0).and.(NodeType(i)== 0)) Then
            NODETYPE(I)=71 
        END IF
    END DO
    
    
    
    End SubRoutine FirstMeshLayer_Adj_Structure
    
subroutine MESH_ALE
use constant
use variable
implicit none
integer i,j


   
    I=1
	IF (I==1) THEN
	OPEN(UNIT=10,FILE='MESH_ALE.DAT')
	WRITE(10,*) ' TITLE=" MESHEXTWOWN " '
	WRITE(10,*) ' VARIABLES="x","y"'
	WRITE(10,*) ' ZONE N=', NN_ORIG,', E=',NE_ORIG,', F=FEPOINT, ET=QUADRILATERAL'
	WRITE(10,*) ' DT=(SINGLE SINGLE )'
	
    DO I=1,NN_ORIG
	WRITE(10,*) X(I),Y(I)
    END DO
	
    DO I=1,NE_ORIG
	WRITE(10,*) (NodeEle(I,J),J=1,4)
	END DO
	CLOSE(10)
	END IF
  
end subroutine MESH_ALE




SUBROUTINE ANS_MESH_READ              !没看
	USE VARIABLE
	USE CONSTANT
    IMPLICIT NONE
    
    INTEGER:: I,J,K

	INTEGER,ALLOCATABLE :: NODEELE_LIU   (:,:)
	REAL(8),ALLOCATABLE :: COORNODE (:,:)
	REAL(8)             :: TEMP
      
                  
	!***********************************************************************************
	ALLOCATE(COORNODE(NN_orig,6))

11	FORMAT(BN,1/)
12	FORMAT(BN,2/)
14	FORMAT(BN,4/)
19	FORMAT(BN,9/)  ! BLANK FOR NONE LINES.
21	FORMAT(BN,11/)
23	FORMAT(BN,13/)	

	OPEN(UNIT=10,FILE='NLIST.lis',STATUS='OLD')
	READ(10,14)
 
	K=0
		
	DO WHILE (.TRUE.)

		DO I=1,20
		K=K+1
		IF (K>NN_orig) EXIT

		READ(10,*,END=100) TEMP,COORNODE (K,1),COORNODE (K,2),TEMP

		END DO

		IF (K>=NN_orig) EXIT

		READ(10,11,END=100) 

	END DO
100   CLOSE(10)


	!---------------------------------------------

	

	ALLOCATE(NODEELE_LIU(NE_ORIG,4))

	OPEN(UNIT=10,FILE='ELIST.lis',STATUS='OLD')
	READ(10,14)
	K=0
	DO WHILE(.TRUE.)
		DO I=1,20	
		K=K+1
		IF (K>NE_ORIG) EXIT	
		READ(10,*,END=101) TEMP,TEMP,TEMP,TEMP,TEMP,TEMP,NODEELE_LIU(k,1),NODEELE_LIU(k,2),NODEELE_LIU(k,3),NODEELE_LIU(k,4)
		
		END DO
		IF (K>=NE_ORIG) EXIT

		READ(10,12,END=101)

	END DO
101   CLOSE(10)
	


	
	OPEN(UNIT=10,FILE='COORNODE.TXT')


	
	DO I=1,nn_orig
	WRITE(10,*) (COORNODE(I,J),J=1,2)
	END DO
	CLOSE(10)

	OPEN(UNIT=10,FILE='NODEELE.TXT')

	DO I=1,NE_ORIG
	WRITE(10,*) (NODEELE_LIU(I,J),J=1,4)
	END DO

	CLOSE(10)

	open(unit=10,file='mesh_tec.dat')
	WRITE(10,*) 'TITLE="a"'
	WRITE(10,*) 'VARIABLES="X","Y" '

	WRITE(10,'(1X,A7,I6,A4,I6,A31)')  'ZONE N=' , NN_ORIG, ', E=',NE_ORIG, ',F=FEPOINT, ET=QUADRILATERAL'
	
	WRITE(10,*) 'DT=(SINGLE SINGLE )' 

	DO I=1,NN_ORIG
	WRITE(10,*) (COORNODE(I,J),J=1,2)
	END DO
	DO I=1,NE_ORIG
	WRITE(10,*) (NODEELE_LIU(I,J),J=1,4)
	END DO
	close(10)

END  SUBROUTINE ANS_MESH_READ
      
      
      


SUBROUTINE Read_Grid
    use CONSTANT
    use VARIABLE
    IMPLICIT NONE
    integer i,j
     
      
    ALLOCATE(NODEELE_orig(Ne_orig,4))
    ALLOCATE(x_orig(nn_orig))
    ALLOCATE(y_orig(nn_orig))
      
     
    OPEN(UNIT=10,FILE='COORNODE.TXT')
	
	DO I=1,nn_orig
	READ(10,*) x_orig(I),y_orig(i)
	END DO
	CLOSE(10)

	OPEN(UNIT=10,FILE='NODEELE.TXT')

	DO I=1,NE_ORIG
	READ(10,*) (NodeEle_ORIG(I,J),J=1,4)
	END DO

	CLOSE(10)
      
        
    open(unit=10,file='mesh_ORIG.dat')
	WRITE(10,*) 'TITLE="a"'
	WRITE(10,*) 'VARIABLES="X","Y" '

	WRITE(10,'(1X,A7,I6,A4,I6,A31)')  'ZONE N=' , NN_ORIG, ', E=',NE_ORIG, ',F=FEPOINT, ET=QUADRILATERAL'
	
	WRITE(10,*) 'DT=(SINGLE SINGLE )' 

	DO I=1,NN_ORIG
	WRITE(10,*) X_ORIG(I),Y_ORIG(I)
	END DO
	DO I=1,NE_ORIG
	WRITE(10,*) (NodeEle_ORIG(I,J),J=1,4)
	END DO
	close(10)
      
!========================================================================           ! ale
	nn_ale=nn_orig                                                      !
	ne_ale=ne_orig                                                      !  
	allocate(x_ale(nn_ale),y_ale(nn_ale),nodeele_ale(ne_ale,4))         !
	x_ale=x_orig                                                        !
	y_ale=y_orig                                                        !
	nodeele_ale=nodeele_orig                                            !
!=======================================================================!   
      

END  SUBROUTINE Read_Grid



SUBROUTINE READ_LAPLACE
USE CONSTANT
USE VARIABLE
IMPLICIT NONE
integer i,j
     
    ALLOCATE(NODEELE_orig(Ne_orig,4))
    ALLOCATE(x_orig(nn_orig))
    ALLOCATE(y_orig(nn_orig))
      

      
    open(unit=10,file='LAPLACE_NODE.TXT',STATUS='OLD')

    do i=1,nn_orig
    read(10,*) x_orig(I),y_orig(i)
    end do

    CLOSE(10)
    
            
    open(unit=10,file='LAPLACE_ELE.TXT',status='old')
    do i=1,ne_orig
    read(10,*) nodeele_orig(i,1),nodeele_orig(i,2),nodeele_orig(i,3),nodeele_orig(i,4)
    end do
    CLOSE(10)

  
    END SUBROUTINE READ_LAPLACE
    
subroutine gmsh_m_read
use constant
use variable
implicit none
integer i,j !JY!i,j循环变量 temp临时存储(不用管)
real(8) temp,temp2
integer nc !JY!num of curve
!character*2 temp_c

open(unit=15,file='roundedSquareCylinderWithoutplate.m')  

read(15,*)
    read(15,*) nc,nn_orig,ne_orig
    write(*,*)  'zone nn_orig=' , nn_orig, ', ne_orig=',ne_orig
    allocate(nodeele_orig(ne_orig,4))
    allocate(x_orig(nn_orig))
    allocate(y_orig(nn_orig))
    
    DO I=1,6
        read(15,*)     
    END DO
    
    !读取各节点坐标
	do i=1,nn_orig    
	    read(15,*) x_orig(i),y_orig(i)
    end do
    
    do i=1,4+nc+2
    read(15,*)
    end do

    !JY!读取各节点编号
	do i=1,ne_orig 
	read(15,*) (nodeele_orig(i,j),j=4,1,-1)    
    end do
    
close(15)

!JY!输出初始网格，可以用tecplot画
    open(unit=10,file='MESH_ORIG.dat')
	write(10,*) 'title="a"'
	write(10,*) 'variables="x","y" '

    !JY!输出外网格
	write(10,'(1x,a7,i6,a4,i6,a31)')  'zone n=' , nn_orig, ', e=',ne_orig, ',f=fepoint, et=quadrilateral'
	
	write(10,*) 'dt=(single single )' 

	do i=1,nn_orig
	write(10,*) x_orig(i),y_orig(i)
    end do
    
	do i=1,ne_orig
	write(10,*) (nodeele_orig(i,j),j=1,4)
    end do
    CLOSE(10)
    
    !========================================================================           ! ale会用到
	nn_ale=nn_orig                                                                  !
	ne_ale=ne_orig                                                                  !  
	allocate(x_ale(nn_ale),y_ale(nn_ale),nodeele_ale(ne_ale,4),nodetype(nn_ale))    !
	x_ale=x_orig                                                                    !
	y_ale=y_orig                                                                    !
	nodeele_ale=nodeele_orig                                                        !
	nodetype=0                                                                      !
!================================================================== ====   XHL  更新后的网格       !
    do i=1,nn_ale                  
       
        if((x_ale(i)*x_ale(i)+y_ale(i)*y_ale(i))>=spring_in**2+0.02.and.(x_ale(i)*x_ale(i)+y_ale(i)*y_ale(i))<spring_out**2+0.02)then
            nodetype(i)=2 !JY!nodetype()在spring_method里会出现，对应弹簧网格不同的处理方法!JY!注意以后该成滑移网格时会用到
        else if((x_ale(i)**2+y_ale(i)**2)<spring_in**2+0.02)then
            nodetype(i)=3
        end if
        
    end do
!===========================================================================!XHL 更新后的网格
      
end subroutine gmsh_m_read