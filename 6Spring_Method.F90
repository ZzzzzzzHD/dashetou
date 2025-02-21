
    SUBROUTINE NODE_ADJ_NODE_LMM
    USE CONSTANT
    USE VARIABLE
    IMPLICIT NONE

    INTEGER:: I,J,K,NUM,NODE1,NODE2,ii,nk,m
    ALLOCATE(NODE_NODE_LMM(NN_ALE,20))

    ii=0
    if(ii==1) then

        DO I=1,NN_ALE
            NUM=0
            DO J=1,NE_ALE
                DO K=1,4
                    IF(K<=3)THEN
                        NODE1=NODEELE_ORIG(J,K)
                        NODE2=NODEELE_ORIG(J,K+1)
                    ELSE
                        NODE1=NODEELE_ORIG(J,4)
                        NODE2=NODEELE_ORIG(J,1)                     !与节点i相连的节点
                    END IF

                    IF(NODE1==I)THEN
                        NUM=NUM+1
                        NODE_NODE_LMM(I,NUM)=NODE2
                    ELSEIF(NODE2==I)THEN
                        NUM=NUM+1
                        NODE_NODE_LMM(I,NUM)=NODE1
                    END IF
                END DO
            END DO
        END DO

    end if

    do i=1,nn_ale
        nk=0
        do j=1,30
            k=node_ele(i,j)                              !节点i所属的单元k
            if (k/=0) then
                do m=1,4

                    if(m<4) then
                        node1=nodeele_orig(k,m)
                        node2=nodeele_orig(k,m+1)
                    else
                        node1=nodeele_orig(k,4)
                        node2=nodeele_orig(k,1)
                    end if

                    if(node1==i) then
                        nk=nk+1
                        if (nk>20) THEN
                            WRITE(*,*)'ghv669'
                            PAUSE
                            stop
                        END IF
                        NODE_NODE_LMM(i,nk)= node2
                    else if (node2==i) then
                        nk=nk+1
                        if (nk>20) THEN
                            WRITE(*,*)'ghv668'
                            PAUSE
                            stop
                        END IF
                        NODE_NODE_LMM(i,nk)= node1
                    end if
                end do
            end if
        end do
    end do


    DO I=1,NN_ALE
        DO J=1,19
            DO K=J+1,20
                IF( NODE_NODE_LMM(I,J)==NODE_NODE_LMM(I,K))THEN
                    NODE_NODE_LMM(I,K)=0
                END IF
            END DO
        END DO
    END DO

    END SUBROUTINE NODE_ADJ_NODE_LMM


    SUBROUTINE For_FindMinELEArea
    USE CONSTANT
    USE VARIABLE
    IMPLICIT NONE

    INTEGER:: I,J,K,NUM,NODE1,NODE2,nk,m,Net(4),Index,ii,jj,qq
    Integer,Allocatable::StoreELES(:,:)
    Integer,Allocatable::RELES0(:),RELES1(:)

    Allocate(NODE_NODE_LMM_ELE(NN_ALE,MAXNUM_NODE_NODE,2))
    ALLOCATE(StoreELES(NN_ALE,20))
    NODE_NODE_LMM_ELE = 0
    StoreELES = 0

    Allocate(RELES0(20),RELES1(20))

    Do i=1,nn_ale
        Index = 1
        Do j = 1,Ne_ale
            Net(:) = NodeELE_ALE(j,1:4)
            Do k = 1,4
                Node1 = Net(k)
                If (Node1 == i) then
                    StoreELES(I,Index) = j
                    Index = Index +1
                    If (Index > 20) Then
                        Write(*,*) 'How will There be a Line Connects More Than 20 Elements !'
                    End If
                End IF
            End Do
        End Do
    end do

    Do I = 1, NN_ale
        RELES0(:) = StoreELES(i,:)
        Do j = 1, MAXNUM_NODE_NODE
            qq=1
            K = NODE_NODE_LMM(I,J)
            IF( K > 0)THEN
                RELES1(:) = StoreELES(k,:)
            End if
            Do ii = 1 ,size(RELES0)
                Do jj = 1,size(RELES1)
                    IF (RELES0(ii)/=0.and.(RELES1(jj)/=0.and.RELES0(ii)==RELES1(jj))) then
                        NODE_NODE_LMM_ELE(i,j,qq) = RELES0(ii)
                        qq = qq + 1
                        If (qq > 3) Then
                            Write(*,*) 'Something Wrong !'
                        End If
                    End if
                End Do
            End do
        End do
    End do

    !DO I=1,NN_ALE
    !    DO J=1,19
    !        DO K=J+1,20
    !            IF( NODE_NODE_LMM_ELE(I,J,:)==NODE_NODE_LMM_ELE(I,K,:))THEN
    !            NODE_NODE_LMM_ELE(I,K,:)=0
    !            END IF
    !        END DO
    !    END DO
    !END DO

    DEAllocate(RELES0,RELES1)
    DEALLOCATE(StoreELES)

    END SUBROUTINE For_FindMinELEArea


    SUBROUTINE Iter_Area
    USE VARIABLE
    USE CONSTANT

    INTEGER:: I,J,K,NODE1,NODE2,IP,K0
    REAL(8):: S,XX2,YY2,XX1,YY1
    REAL(8):: NODEELE_INVERSE(4)

    IF (ALLOCATED(SELE_Iter)) DEALLOCATE(SELE_Iter)

    ALLOCATE(SELE_Iter(NE_ALE))
    SELE_Iter(NE_ALE) = 0.0d0

    DO I=1,NE_ALE
        S=0.0D0
        IP=IPE(I)
        K0=IP-1
        DO J=1,IP
            IF (J<=K0) THEN
                NODE1=NODEELE_ALE(I,J)
                NODE2=NODEELE_ALE(I,J+1)
            ELSE
                NODE1=NODEELE_ALE(I,4)
                NODE2=NODEELE_ALE(I,1)
            END IF
            XX1=X0_ALE(NODE1)
            YY1=Y0_ALE(NODE1)
            XX2=X0_ALE(NODE2)
            YY2=Y0_ALE(NODE2)
            S=S+(XX1*YY2-XX2*YY1)
        END DO

        IF (S>0.0D0) THEN
            SELE_Iter(I)=S*0.5D0
        ELSE
            WRITE(*,*) ' Why Still Negative Area ?'
            WRITE(*,*) ' THE ELEMENT IS A :', I
            WRITE(*,*) S
            PAUSE
        END IF
    END DO
    END SUBROUTINE


    SUBROUTINE Iter_Angle
    USE VARIABLE
    USE CONSTANT

    INTEGER::i,j,KNode,Index,qq,NodeTemp,Node(2),Net(4),EleIndex,KEle
    REAL(8)::Angle(4),aaLength,bbLength,ccLength,ATempData
    REAL(8):: pi,Error

    pi=dacos(-1.0d0)
    Do I = 1,NN_ALE
        Do j = 1, MAXNUM_NODE_NODE
            KNode = NODE_NODE_LMM(I,J)
            IF (Knode > 0) Then
                Do EleIndex = 1,2
                    KEle  = NODE_NODE_LMM_ELE(I,J,EleIndex)
                    IF( KEle /= 0)THEN
                        Net(:) = NodeELE_ALE(KEle,1:4)
                        qq=1
                        Do ii = 1,4
                            NodeTemp = Net(ii)
                            If (NodeTemp/=I.and.NodeTemp/=Knode) Then
                                Node(qq) = NodeTemp
                                qq = qq+1
                                If(qq > 3)Then
                                    Write(*,*) 'Something Bad !', I,KNode,KEle
                                End If
                            End If
                        End Do

                        aaLength = (X0_Ale(Node(1))-X0_Ale(KNode))**2+(Y0_Ale(Node(1))-Y0_Ale(KNode))**2
                        bbLength = (X0_Ale(I)-X0_Ale(Node(1)))**2+(Y0_Ale(I)-Y0_Ale(Node(1)))**2
                        ccLength = (X0_Ale(I)-X0_Ale(KNode))**2+(Y0_Ale(I)-Y0_Ale(KNode))**2
                        ATempData = (bbLength+ccLength-aaLength)/(2*dsqrt(bbLength)*dsqrt(ccLength))
                        Angle(1) =Dacos(ATempData) ! I 点和Node(1)连接

                        aaLength = (X0_Ale(I)-X0_Ale(Node(2)))**2+(Y0_Ale(I)-Y0_Ale(Node(2)))**2
                        bbLength = (X0_Ale(KNode)-X0_Ale(Node(2)))**2+(Y0_Ale(KNode)-Y0_Ale(Node(2)))**2
                        ccLength = (X0_Ale(I)-X0_Ale(KNode))**2+(Y0_Ale(I)-Y0_Ale(KNode))**2
                        ATempData = (bbLength+ccLength-aaLength)/(2.0d0*dsqrt(bbLength)*dsqrt(ccLength))
                        Angle(2) =Dacos(ATempData) ! knode 和 Node(2) 相连

                        IF ( Angle(1)==0.or. Angle(2)==0) then
                            Write(*,*) I,J,Knode
                            Pause
                        End if

                        aaLength = (X0_Ale(I)-X0_Ale(Node(2)))**2+(Y0_Ale(I)-Y0_Ale(Node(2)))**2
                        bbLength = (X0_Ale(I)-X0_Ale(Node(1)))**2+(Y0_Ale(I)-Y0_Ale(Node(1)))**2
                        ccLength = (X0_Ale(Node(1))-X0_Ale(Node(2)))**2+(Y0_Ale(Node(1))-Y0_Ale(Node(2)))**2
                        ATempData = (bbLength+ccLength-aaLength)/(2.0d0*dsqrt(bbLength)*dsqrt(ccLength))
                        Angle(3) =Dacos(ATempData)

                        aaLength = (X0_Ale(Node(1))-X0_Ale(KNode))**2+(Y0_Ale(Node(1))-Y0_Ale(KNode))**2
                        bbLength = (X0_Ale(KNode)-X0_Ale(Node(2)))**2+(Y0_Ale(KNode)-Y0_Ale(Node(2)))**2
                        ccLength = (X0_Ale(Node(1))-X0_Ale(Node(2)))**2+(Y0_Ale(Node(1))-Y0_Ale(Node(2)))**2
                        ATempData = (bbLength+ccLength-aaLength)/(2.0d0*dsqrt(bbLength)*dsqrt(ccLength))
                        Angle(4) =Dacos(ATempData)

                        Error = Dabs(Sum(Angle(1:4))-2*pi)
                        !Write(*,*) I,KNODE,KELE,ERROR,Node(1:2)
                        !Pause

                        If ( Error < 1e-6) Then
                            If (EleIndex == 1)then
                                NODE_NODE_LMM_Angle(I,J,1) = Angle(1)
                                NODE_NODE_LMM_Angle(I,J,2) = Angle(2)
                            Else
                                NODE_NODE_LMM_Angle(I,J,3) = Angle(1)
                                NODE_NODE_LMM_Angle(I,J,4) = Angle(2)
                            End If
                        Else
                            aaLength = (X0_Ale(Node(2))-X0_Ale(KNode))**2+(Y0_Ale(Node(2))-Y0_Ale(KNode))**2
                            bbLength = (X0_Ale(I)-X0_Ale(Node(2)))**2+(Y0_Ale(I)-Y0_Ale(Node(2)))**2
                            ccLength = (X0_Ale(I)-X0_Ale(KNode))**2+(Y0_Ale(I)-Y0_Ale(KNode))**2
                            ATempData = (bbLength+ccLength-aaLength)/(2.0d0*dsqrt(bbLength)*dsqrt(ccLength))
                            Angle(1) =Dacos(ATempData) ! I 点和Node(2)连接

                            aaLength = (X0_Ale(I)-X0_Ale(Node(1)))**2+(Y0_Ale(I)-Y0_Ale(Node(1)))**2
                            bbLength = (X0_Ale(KNode)-X0_Ale(Node(1)))**2+(Y0_Ale(KNode)-Y0_Ale(Node(1)))**2
                            ccLength = (X0_Ale(I)-X0_Ale(KNode))**2+(Y0_Ale(I)-Y0_Ale(KNode))**2
                            ATempData = (bbLength+ccLength-aaLength)/(2.0d0*dsqrt(bbLength)*dsqrt(ccLength))
                            Angle(2) =Dacos(ATempData) ! knode 和 Node(1) 相连

                            aaLength = (X0_Ale(I)-X0_Ale(Node(1)))**2+(Y0_Ale(I)-Y0_Ale(Node(1)))**2
                            bbLength = (X0_Ale(I)-X0_Ale(Node(2)))**2+(Y0_Ale(I)-Y0_Ale(Node(2)))**2
                            ccLength = (X0_Ale(Node(1))-X0_Ale(Node(2)))**2+(Y0_Ale(Node(1))-Y0_Ale(Node(2)))**2
                            ATempData = (bbLength+ccLength-aaLength)/(2.0d0*dsqrt(bbLength)*dsqrt(ccLength))
                            Angle(3) =Dacos(ATempData)

                            aaLength = (X0_Ale(Node(2))-X0_Ale(KNode))**2+(Y0_Ale(Node(2))-Y0_Ale(KNode))**2
                            bbLength = (X0_Ale(Node(1))-X0_Ale(Node(2)))**2+(Y0_Ale(Node(1))-Y0_Ale(Node(2)))**2
                            ccLength = (X0_Ale(Node(1))-X0_Ale(KNode))**2+(Y0_Ale(Node(1))-Y0_Ale(KNode))**2
                            ATempData = (bbLength+ccLength-aaLength)/(2.0d0*dsqrt(bbLength)*dsqrt(ccLength))
                            Angle(4) =Dacos(ATempData)

                            Error = Dabs(Sum(Angle(1:4))-2*pi)
                            !WRITE(*,*) 'ERRORRRRRRRRRRRR',I,KNODE,KELE,ERROR,Node(1:2)
                            !Pause
                            IF(Error < 1e-6 ) then
                                If (EleIndex == 1)then
                                    NODE_NODE_LMM_Angle(I,J,1) = Angle(1)
                                    NODE_NODE_LMM_Angle(I,J,2) = Angle(2)
                                Else
                                    NODE_NODE_LMM_Angle(I,J,3) = Angle(1)
                                    NODE_NODE_LMM_Angle(I,J,4) = Angle(2)
                                End If
                            Else
                                Write(*,*) 'How ? Precision ?' ,Error
                                Pause
                                Stop
                            End If

                        End If
                    End if
                End Do
            End If

            If(NodeType(i)==1.and.NodeType(j)==1) Then
                NODE_NODE_LMM_Angle(I,J,3) = pi/2
                NODE_NODE_LMM_Angle(I,J,4) = pi/2
            End IF
        End do
    End do

    END SUBROUTINE



    SUBROUTINE SPRING_METHOD
    Use FlexibleOnly
    USE CONSTANT
    USE VARIABLE
    IMPLICIT NONE

    INTEGER::NUM_LMM,NUM_LIU,ElementIndex,FNUM,Types_NodeTypes(6,2),FIIJK,FIJK
    INTEGER,ALLOCATABLE::NODE_INDEX_TYPES(:)
    REAL(8)::ttt0,ttt1,Flo_xi,Flo_S(4)
    REAL(8)::KIJ(maxnum_node_node),KPLUS,PI,X0,Y0
    INTEGER::I,J,N,k,num_it,num
    REAL(8)::SMALL=1.0d-8,theta1,theta2,theta3,R00,dd
    real(8)::spring_ring_width,InternalAngel,FA_SPR71,FA_SPR72,FA_SPR73,FA_SPR74,FA_SPR75,FA_SPR1
    REAL(8),ALLOCATABLE:: DELTX00(:),DELTY00(:)

    INTEGER::JJ,NKindSN,Index_X
    REAL(8),ALLOCATABLE::SORTED_X_INDEX_ARR(:,:),AddiDis(:)
    REAL(8)::xxFtemp,xFtemp,FDELTA_X,FDELTA_Y
    REAL(8)::X_INDEX_ARR(NSBN_ALE,2)
    REAL(8)::FOmega,SCL
    
    ALLOCATE(deltx00(nn_ale),delty00(nn_ale))
    

    xxFtemp = 0
    xFtemp  = 0
    
    spring_ring_width=spring_out-spring_in

    pi=dacos(-1.0d0)
    num=0


    IF (IROTA==1) THEN               !圆柱与板绕圆心以角度ACOS(WT+PI/2)旋转
        !====================================================================更新后的网格 XHL
        DO I=1,NN_ALE
            IF(NODETYPE(I)==3)THEN	      !网格跟随系统旋转

                X0=X_ALE(I)
                Y0=Y_ALE(I)
                DD=DSQRT(X0*X0+Y0*Y0)
                THETA2=DASIN(Y0/DD)
                IF (X0<0.0D0)   THEN
                    THETA2=PI-THETA2
                END IF

                IF(T>T_START)THEN
                    DELTX(I)=DD*DCOS(DIS_ANG+THETA2)-DD*DCOS(THETA2+DIS0_ANG)
                    DELTY(I)=DD*DSIN(DIS_ANG+THETA2)-DD*DSIN(THETA2+DIS0_ANG)                         !+DIS_Y-DIS0_Y
                ELSE
                    DELTX(I)=0.0D0        !JY!ORIG:DELTX(J)=0.0D0
                    DELTY(I)=0.0D0        !JY!ORIG:DELTY(J)=0.0D0
                END IF

            ELSE IF(NODETYPE(I)==2)THEN     !处理网格旋转变形

                X0=X_ALE(I)
                Y0=Y_ALE(I)

                DD=DSQRT(X0*X0+Y0*Y0)
                THETA2=DASIN(Y0/DD)
                IF (X0<0.0D0)   THEN
                    THETA2=PI-THETA2
                END IF
                IF( T>T_START ) THEN
                    DELTX(I)=DD*DCOS(DIS_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH)+THETA2)-DD*DCOS(THETA2+DIS0_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH))
                    DELTY(I)=DD*DSIN(DIS_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH)+THETA2)-DD*DSIN(THETA2+DIS0_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH))           !+DIS_Y-DIS0_Y
                ELSE
                    DELTX(I)=0.0D0        !JY!ORIG:DELTX(J)=0.0D0
                    DELTY(I)=0.0D0        !JY!ORIG:DELTY(J)=0.0D0
                END IF

            ENDIF
        END DO
        !===================================================================================================更新后的网格 XHL

        IF(TT0>=T_START)THEN
            !THETA1=AMC_ALE*DCOS(OMG_ALE*(T-T_START)+PI/2.0D0)
            !THETA3=AMC_ALE*DCOS(OMG_ALE*(TT0-T_START)+PI/2.0D0)     !THETA3是THETA1的前一步
            DIS_ANG=AMC_ALE*DCOS(OMG_ALE*(T-T_START)+PI/2.0D0)        !PDD
            DIS0_ANG=AMC_ALE*DCOS(OMG_ALE*(TT0-T_START)+PI/2.0D0)     !PDD
        ELSE
            ! THETA1=0.0D0
            ! THETA3=0.0D0
            DIS_ANG=0.0D0      !PDD
            DIS0_ANG=0.0D0     !PDD
            GOTO 100
        END IF

        DO I=1,NSBN_ALE
            NUM=NUM+1
            J=SBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)

            !DELTX(J)=X0*(DCOS(THETA1)-DCOS(THETA3))-Y0*(DSIN(THETA1)-DSIN(THETA3))
            !DELTY(J)=X0*(DSIN(THETA1)-DSIN(THETA3))+Y0*(DCOS(THETA1)-DCOS(THETA3))
            DELTX(J)=X0*(DCOS(DIS_ANG)-DCOS(DIS0_ANG))-Y0*(DSIN(DIS_ANG)-DSIN(DIS0_ANG))    !PDD
            DELTY(J)=X0*(DSIN(DIS_ANG)-DSIN(DIS0_ANG))+Y0*(DCOS(DIS_ANG)-DCOS(DIS0_ANG))    !PDD
        END DO

        DO I=1,NPBN_ALE
            NUM=NUM+1
            J=PBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)

            DD=DSQRT(X0*X0+Y0*Y0)
            THETA2=DASIN(Y0/DD)
            IF (X0<0.0D0)   THEN
                THETA2=PI-THETA2
            END IF

            !DELTX(J)=DD*DCOS(THETA1+THETA2)-DD*DCOS(THETA2+THETA3)
            !DELTY(J)=DD*DSIN(THETA1+THETA2)-DD*DSIN(THETA2+THETA3)
            DELTX(J)=DD*DCOS(DIS_ANG+THETA2)-DD*DCOS(THETA2+DIS0_ANG)  !PDD
            DELTY(J)=DD*DSIN(DIS_ANG+THETA2)-DD*DSIN(THETA2+DIS0_ANG)  !PDD
        END DO


        !ELSE IF (IROTA==10) THEN          !圆柱固定，板绕圆心以板端Y=ASIN(WT)旋转，与参考文献73验证
        !
        !      IF(TT0>=T_START)THEN
        !         THETA1=DASIN(AMC_ALE*DSIN(OMG_ALE*(T-T_START)))
        !         THETA3=DASIN(AMC_ALE*DSIN(OMG_ALE*(TT0-T_START)))
        !      ELSE
        !         THETA1=0.0D0
        !         THETA3=0.0D0
        !         GOTO 100
        !      END IF
        !
        !
        !      DO I=1,NSBN_ALE
        !      NUM=NUM+1
        !   J=SBN_ALE(I)
        !   X0=X_ALE(J)
        !   Y0=Y_ALE(J)
        !
        !   R00=X0-X00
        !   DELTX(J)=R00*(DCOS(THETA1)-DCOS(THETA3))
        !   DELTY(J)=R00*(DSIN(THETA1)-DSIN(THETA3))
        !      END DO
        !
        !      DO I=1,NPBN_ALE
        !	   NUM=NUM+1
        !   J=PBN_ALE(I)
        !   DELTX(J)=0.0D0
        !   DELTY(J)=0.0D0
        !      END DO
        !
        !
        !ELSE IF(IROTA==11)THEN            !圆柱固定，板绕圆心以角度ACOS(WT+PI/2)旋转，与参考文献72验证
        !
        !      IF(TT0>=T_START)THEN
        !         THETA1=AMC_ALE*DCOS(OMG_ALE*(T-T_START)+PI/2.0D0)
        !         THETA3=AMC_ALE*DCOS(OMG_ALE*(TT0-T_START)+PI/2.0D0)
        !      ELSE
        !         THETA1=0.0D0
        !         THETA3=0.0D0
        !         GOTO 100
        !      END IF
        !
        !
        !      DO I=1,NSBN_ALE
        !      NUM=NUM+1
        !   J=SBN_ALE(I)
        !   X0=X_ALE(J)
        !   Y0=Y_ALE(J)
        !
        !   R00=X0-X00
        !   DELTX(J)=R00*(DCOS(THETA1)-DCOS(THETA3))
        !   DELTY(J)=R00*(DSIN(THETA1)-DSIN(THETA3))
        !      END DO
        !
        !      DO I=1,NPBN_ALE
        ! 	   NUM=NUM+1
        !   J=PBN_ALE(I)
        !   DELTX(J)=0.0D0
        !   DELTY(J)=0.0D0
        !      END DO

    ELSE IF(IROTA==0)THEN            !圆柱与板固定

        DO I=1,NSBN_ALE
            NUM=NUM+1
            J=SBN_ALE(I)
            DELTX(J)=0.0D0
            DELTY(J)=0.0D0
        END DO

        DO I=1,NPBN_ALE
            NUM=NUM+1
            J=PBN_ALE(I)
            DELTX(J)=0.0D0
            DELTY(J)=0.0D0
        END DO

    ELSE IF(IROTA==2)THEN            !圆柱与板绕圆心自由旋转
        !====================================================================更新后的网格 XHL
        DO I=1,NN_ALE
            IF(NODETYPE(I)==3)THEN	      !网格跟随系统旋转

                X0=X_ALE(I)
                Y0=Y_ALE(I)
                DD=DSQRT(X0*X0+Y0*Y0)
                THETA2=DASIN(Y0/DD)
                IF (X0<0.0D0)   THEN
                    THETA2=PI-THETA2
                END IF

                IF(T>T_START)THEN
                    DELTX(I)=DD*DCOS(DIS_ANG+THETA2)-DD*DCOS(THETA2+DIS0_ANG)
                    DELTY(I)=DD*DSIN(DIS_ANG+THETA2)-DD*DSIN(THETA2+DIS0_ANG)                         !+DIS_Y-DIS0_Y

                ELSE
                    DELTX(I)=0.0D0        !JY!ORIG:DELTX(J)=0.0D0
                    DELTY(I)=0.0D0        !JY!ORIG:DELTY(J)=0.0D0
                END IF

            ELSE IF(NODETYPE(I)==2)THEN     !处理网格旋转变形

                X0=X_ALE(I)
                Y0=Y_ALE(I)

                DD=DSQRT(X0*X0+Y0*Y0)
                THETA2=DASIN(Y0/DD)
                IF (X0<0.0D0)   THEN
                    THETA2=PI-THETA2
                END IF

                !JY!ORIG
                !IF( T>T_START ) THEN
                ! DELTX(I)=DD*DCOS(DIS_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH)+THETA2)-DD*DCOS(THETA2+DIS0_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH))
                ! DELTY(I)=DD*DSIN(DIS_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH)+THETA2)-DD*DSIN(THETA2+DIS0_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH))
                !    DELTX(J)=0.0D0
                !    DELTY(J)=0.0D0
                !END IF

                !JY!NEW
                IF( T>T_START ) THEN
                    DELTX(I)=DD*DCOS(DIS_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH)+THETA2)-DD*DCOS(THETA2+DIS0_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH))
                    DELTY(I)=DD*DSIN(DIS_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH)+THETA2)-DD*DSIN(THETA2+DIS0_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH))
                ELSE
                    DELTX(I)=0.0D0        !JY!ORIG:DELTX(J)=0.0D0
                    DELTY(I)=0.0D0        !JY!ORIG:DELTY(J)=0.0D0
                END IF

            ENDIF
        END DO
        !===================================================================================================更新后的网格 XHL
        DO I=1,NSBN_ALE
            NUM=NUM+1
            J=SBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)

            IF( T>T_START ) THEN
                DELTX(J)=X0*(DCOS(DIS_ANG)-DCOS(DIS0_ANG))-Y0*(DSIN(DIS_ANG)-DSIN(DIS0_ANG))
                DELTY(J)=X0*(DSIN(DIS_ANG)-DSIN(DIS0_ANG))+Y0*(DCOS(DIS_ANG)-DCOS(DIS0_ANG))
            ELSE
                DELTX(J)=0.0D0
                DELTY(J)=0.0D0
            END IF

        END DO

        DO I=1,NPBN_ALE
            NUM=NUM+1
            J=PBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)

            DD=DSQRT(X0*X0+Y0*Y0)
            THETA2=DASIN(Y0/DD)
            IF (X0<0.0D0)   THEN
                THETA2=PI-THETA2
            END IF

            IF( T>T_START ) THEN
                DELTX(J)=DD*DCOS(DIS_ANG+THETA2)-DD*DCOS(THETA2+DIS0_ANG)
                DELTY(J)=DD*DSIN(DIS_ANG+THETA2)-DD*DSIN(THETA2+DIS0_ANG)
            ELSE
                DELTX(J)=0.0D0
                DELTY(J)=0.0D0
            END IF

        END DO

    ELSE IF(IROTA==200.OR.IROTA==201.OR.IROTA==202.OR.IROTA==2012)THEN            !JY!2DOF
        !====================================================================更新后的网格 XHL
        DO I=1,NN_ALE
            IF(NODETYPE(I)==3)THEN	      !网格跟随系统旋转

                X0=X_ALE(I)
                Y0=Y_ALE(I)
                DD=DSQRT(X0**2+Y0**2)
                THETA2=DASIN(Y0/DD)
                IF (X0<0.0D0)   THEN
                    THETA2=PI-THETA2
                END IF

                IF(T>T_START)THEN
                    DELTX(I)=DD*DCOS(DIS_ANG+THETA2)-DD*DCOS(THETA2+DIS0_ANG)
                    DELTY(I)=DD*DSIN(DIS_ANG+THETA2)-DD*DSIN(THETA2+DIS0_ANG)+DIS_Y-DIS0_Y                         !+DIS_Y-DIS0_Y

                ELSE
                    DELTX(I)=0.0D0        !JY!ORIG:DELTX(J)=0.0D0
                    DELTY(I)=0.0D0        !JY!ORIG:DELTY(J)=0.0D0
                END IF

            ELSE IF(NODETYPE(I)==2)THEN     !处理网格旋转变形

                X0=X_ALE(I)
                Y0=Y_ALE(I)

                DD=DSQRT(X0**2+Y0**2)
                THETA2=DASIN(Y0/DD)
                IF (X0<0.0D0)   THEN
                    THETA2=PI-THETA2
                END IF
                IF( T>T_START ) THEN
                    DELTX(I)=DD*DCOS(DIS_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH)+THETA2)-DD*DCOS(THETA2+DIS0_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH))
                    DELTY(I)=DD*DSIN(DIS_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH)+THETA2)-DD*DSIN(THETA2+DIS0_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH))+DIS_Y-DIS0_Y
                ELSE
                    DELTX(I)=0.0D0        !JY!ORIG:DELTX(J)=0.0D0
                    DELTY(I)=0.0D0        !JY!ORIG:DELTY(J)=0.0D0
                END IF

            ENDIF
        END DO
        !===================================================================================================更新后的网格 XHL
        DO I=1,NSBN_ALE
            NUM=NUM+1
            J=SBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)

            IF( T>T_START ) THEN
                DELTX(J)=X0*(DCOS(DIS_ANG)-DCOS(DIS0_ANG))-Y0*(DSIN(DIS_ANG)-DSIN(DIS0_ANG))
                DELTY(J)=X0*(DSIN(DIS_ANG)-DSIN(DIS0_ANG))+Y0*(DCOS(DIS_ANG)-DCOS(DIS0_ANG))+DIS_Y-DIS0_Y
            ELSE
                DELTX(J)=0.0D0
                DELTY(J)=0.0D0
            END IF

        END DO

        DO I=1,NPBN_ALE
            NUM=NUM+1
            J=PBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)

            DD=DSQRT(X0**2+Y0**2)
            THETA2=DASIN(Y0/DD)
            IF (X0<0.0D0)   THEN
                THETA2=PI-THETA2
            END IF

            IF( T>T_START ) THEN
                DELTX(J)=DD*DCOS(DIS_ANG+THETA2)-DD*DCOS(THETA2+DIS0_ANG)
                DELTY(J)=DD*DSIN(DIS_ANG+THETA2)-DD*DSIN(THETA2+DIS0_ANG)+DIS_Y-DIS0_Y
            ELSE
                DELTX(J)=0.0D0
                DELTY(J)=0.0D0
            END IF

        END DO



    ELSE IF(IROTA==300)THEN            !JY!3DOF
        !====================================================================更新后的网格 XHL
        DO I=1,NN_ALE
            IF(NODETYPE(I)==3)THEN	      !网格跟随系统旋转

                X0=X_ALE(I)
                Y0=Y_ALE(I)
                DD=DSQRT(X0**2+Y0**2)
                THETA2=DASIN(Y0/DD)
                IF (X0<0.0D0)   THEN
                    THETA2=PI-THETA2
                END IF

                IF(T>T_START)THEN
                    DELTX(I)=DD*DCOS(DIS_ANG+THETA2)-DD*DCOS(THETA2+DIS0_ANG)+DIS_X-DIS0_X
                    DELTY(I)=DD*DSIN(DIS_ANG+THETA2)-DD*DSIN(THETA2+DIS0_ANG)+DIS_Y-DIS0_Y                         !+DIS_Y-DIS0_Y

                ELSE
                    DELTX(I)=0.0D0        !JY!ORIG:DELTX(J)=0.0D0
                    DELTY(I)=0.0D0        !JY!ORIG:DELTY(J)=0.0D0
                END IF

            ELSE IF(NODETYPE(I)==2)THEN     !处理网格旋转变形

                X0=X_ALE(I)
                Y0=Y_ALE(I)

                DD=DSQRT(X0**2+Y0**2)
                THETA2=DASIN(Y0/DD)
                IF (X0<0.0D0)   THEN
                    THETA2=PI-THETA2
                END IF
                IF( T>T_START ) THEN
                    DELTX(I)=DD*DCOS(DIS_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH)+THETA2)-DD*DCOS(THETA2+DIS0_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH))+DIS_X-DIS0_X
                    DELTY(I)=DD*DSIN(DIS_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH)+THETA2)-DD*DSIN(THETA2+DIS0_ANG*((-1.0/SPRING_RING_WIDTH)*DD+SPRING_OUT/SPRING_RING_WIDTH))+DIS_Y-DIS0_Y
                ELSE
                    DELTX(I)=0.0D0        !JY!ORIG:DELTX(J)=0.0D0
                    DELTY(I)=0.0D0        !JY!ORIG:DELTY(J)=0.0D0
                END IF

            ENDIF
        END DO
        !===================================================================================================更新后的网格 XHL
        DO I=1,NSBN_ALE
            NUM=NUM+1
            J=SBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)

            IF( T>T_START ) THEN
                DELTX(J)=X0*(DCOS(DIS_ANG)-DCOS(DIS0_ANG))-Y0*(DSIN(DIS_ANG)-DSIN(DIS0_ANG))+DIS_X-DIS0_X
                DELTY(J)=X0*(DSIN(DIS_ANG)-DSIN(DIS0_ANG))+Y0*(DCOS(DIS_ANG)-DCOS(DIS0_ANG))+DIS_Y-DIS0_Y
            ELSE
                DELTX(J)=0.0D0
                DELTY(J)=0.0D0
            END IF

        END DO

        DO I=1,NPBN_ALE
            NUM=NUM+1
            J=PBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)

            DD=DSQRT(X0**2+Y0**2)
            THETA2=DASIN(Y0/DD)
            IF (X0<0.0D0)   THEN
                THETA2=PI-THETA2
            END IF

            IF( T>T_START ) THEN
                DELTX(J)=DD*DCOS(DIS_ANG+THETA2)-DD*DCOS(THETA2+DIS0_ANG)+DIS_X-DIS0_X
                DELTY(J)=DD*DSIN(DIS_ANG+THETA2)-DD*DSIN(THETA2+DIS0_ANG)+DIS_Y-DIS0_Y
            ELSE
                DELTX(J)=0.0D0
                DELTY(J)=0.0D0
            END IF

        END DO

    ELSE IF(IROTA==5107)THEN ! FORCE OSCILLATION
        
        
        SCL = 0.1
        
        DO I=1,NPBN_ALE
            NUM=NUM+1
            J=PBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)

            DD=DSQRT(X0**2+Y0**2)
            THETA2=DASIN(Y0/DD)
            IF (X0<0.0D0)   THEN
                THETA2=PI-THETA2
            END IF

            IF( T>T_START ) THEN
                DELTX(J)=0.0D0
                DELTY(J)=0.6*SIN(T*SCL*PI)-0.6*SIN((T-DT)*SCL*PI)
            ELSE
                DELTX(J)=0.0D0
                DELTY(J)=0.0D0
            END IF
        END DO

        DO I=1,NSBN_ALE
            NUM=NUM+1
            J=SBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)

            X0=X0-0.5D0
            
            FOmega = 1.5
            
            
            IF( T>T_START ) THEN
                DELTX(J)=0.0D0
                !DELTY(J)=TAN(0.6*SIN(T*0.5*PI)-0.6*SIN((T-DT)*0.5*PI))*X0+0.6*SIN(T*0.5*PI)-0.6*SIN((T-DT)*0.5*PI)
                DELTY(J) = 0.6*SIN(T*SCL*PI)-0.6*SIN((T-DT)*SCL*PI)+(Dsin(T*SCL)-Dsin(SCL*(T-Dt)))*(Dcosh(FOmega*X0)-Dcos(FOmega*X0)+(Dcos(FOmega)+Dcosh(FOmega))/(Dsin(FOmega)+Dsinh(FOmega))*(Dsin(FOmega*X0)-Dsinh(FOmega*X0)))
                
                X_INDEX_ARR(I,2) = 0.6*SIN(T*SCL*PI)+Dsin(T*SCL)*(Dcosh(FOmega*X0)-Dcos(FOmega*X0)+(Dcos(FOmega)+Dcosh(FOmega))/(Dsin(FOmega)+Dsinh(FOmega))*(Dsin(FOmega*X0)-Dsinh(FOmega*X0)))
            ELSE
                DELTX(J)=0.0D0
                DELTY(J)=0.0D0
            END IF
            
            X_INDEX_ARR(I,1) = X_ALE(J)
        END DO

        DO I =1,NSBN_ALE
            xFtemp=X_INDEX_ARR(I,1)
            Do J = I+1,NSBN_ALE
                xxFtemp =X_INDEX_ARR(J,1)
                IF (xFtemp==xxFtemp) THEN
                    X_INDEX_ARR(J,1)=0
                END IF
            End DO
        End do

        NKindSN=0
        DO I =1,NSBN_ALE
            xFtemp=X_INDEX_ARR(I,1)
            IF(xFtemp /= 0) then
                NKindSN =  NKindSN + 1
            End if
        End do

        !PRINT*, NKindSN
        !PAUSE
        
        Allocate(SORTED_X_INDEX_ARR(NKindSN,3),AddiDis(NKindSN))
        SORTED_X_INDEX_ARR = 0.0D0
        AddiDis = 0.0
        
        J = 1
        DO I = 1, NSBN_ALE
            IF (X_INDEX_ARR(I,1)/=0)THEN
                SORTED_X_INDEX_ARR(J,1)=X_INDEX_ARR(I,1)
                SORTED_X_INDEX_ARR(J,2)=X_INDEX_ARR(I,2)
                J=J+1
            END IF
        END DO

        Call Sort_Ascend_Multi(SORTED_X_INDEX_ARR,NKindSN,3) ! 按照第一列的大小进行排列

        PLATE_DS = 0.0
        DO I = NKINDSN,2,-1
            FDELTA_Y = SORTED_X_INDEX_ARR(I,2) - SORTED_X_INDEX_ARR(I-1,2)
            PLATE_DS = X_Plate_Orig(I,1) - X_Plate_Orig(I-1,1)
            FA_THETA(I) = ASIN(FDELTA_Y / PLATE_DS)
            
            FDELTA_X =  (PLATE_DS * DCOS(FA_THETA(I)) - (X_CA_ORR(I,1) - X_CA_ORR(I-1,1)))
            SORTED_X_INDEX_ARR(I,3) = FDELTA_X
            
            !PRINT*, FDELTA_X
        END DO
        SORTED_X_INDEX_ARR(1,3) = 0.0
        
        
        !DO I = 1, NSBN_ALE
        !    J = SBN_ALE(I)
        !    X0 = X_ALE(J)
        !    Y0 = Y_ALE(J)
        !    Index_X = findloc(X_Plate_Orig(:,1),X0,1)
        !    IF(Y0 >= 0.0) THEN
        !        DELTX(J)=DELTX(J) + SUM(SORTED_X_INDEX_ARR(1:Index_X,3)) - Y0 * (DSIN(FA_THETA(Index_X))-DSIN(FA_THETA0(Index_X)))
        !        DELTY(J)=DELTY(J) + Y0 * ((DCOS(FA_THETA(Index_X)) - 1.0)-(DCOS(FA_THETA0(Index_X)) - 1.0))
        !    ELSE 
        !        DELTX(J)=DELTX(J) + SUM(SORTED_X_INDEX_ARR(1:Index_X,3)) - Y0 * (DSIN(FA_THETA(Index_X))-DSIN(FA_THETA0(Index_X)))
        !        DELTY(J)=DELTY(J) - Y0 * ((DCOS(FA_THETA(Index_X)) - 1.0)-(DCOS(FA_THETA0(Index_X)) - 1.0))
        !    END IF
        !    AddiDis(Index_X) = DELTX(J)
        !    
        !END DO
        
        X_CA_ORR(:,1) = X_CA_ORR(:,1) + AddiDis(:)
        
        FA_THETA0 = FA_THETA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
    ELSE IF(IROTA==510)THEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !CALL ITER_AREA
        !CALL ITER_ANGLE

        DO I=1,NPBN_ALE
            NUM=NUM+1
            J=PBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)

            DD=DSQRT(X0**2+Y0**2)
            THETA2=DASIN(Y0/DD)
            IF (X0<0.0D0)   THEN
                THETA2=PI-THETA2
            END IF

            IF( T>T_START ) THEN
                DELTX(J)=DD*DCOS(DIS_ANG+THETA2)-DD*DCOS(THETA2+DIS0_ANG)
                DELTY(J)=DD*DSIN(DIS_ANG+THETA2)-DD*DSIN(THETA2+DIS0_ANG)+DIS_Y-DIS0_Y
            ELSE
                DELTX(J)=0.0D0
                DELTY(J)=0.0D0
            END IF
        END DO

        DO I=1,NSBN_ALE
            !首先要把计算出来的位移进行插值到网格节点上
            NUM=NUM+1
            J=SBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)

            X0=X0-0.5D0
            ELEMENTINDEX=FLOOR(X0/FL_PARAMS_L)

            IF (ELEMENTINDEX<FL_PARAMS_NE)THEN
                FLO_XI=(X0-ELEMENTINDEX*FL_PARAMS_L)/FL_PARAMS_L  !局部坐标
                FLO_S=0.0D0
                CALL SHAPE_FUN(FLO_S,4,FLO_XI,0) !S现在是4*1，转置后，1*4，乘以4*1的位移阵，就得到了节点的挠度

                IF( T>T_START ) THEN
                    DELTX(J)=1.0*(DCOS(DIS_ANG)-DCOS(DIS0_ANG))!两部分，一部分是由于圆柱运动造成的整体向前的位移，另一部分是由于纵向变形造成的X轴的收缩
                    DELTY(J)=FLO_S(1)*FDIS_Y(ELEMENTINDEX+1)+FLO_S(2)*FDIS_ANG(ELEMENTINDEX+1)+ &
                        & FLO_S(3)*FDIS_Y(ELEMENTINDEX+2)+FLO_S(4)*FDIS_ANG(ELEMENTINDEX+2) - &
                        & (FLO_S(1)*FDIS0_Y(ELEMENTINDEX+1)+FLO_S(2)*FDIS0_ANG(ELEMENTINDEX+1) + &
                        &  FLO_S(3)*FDIS0_Y(ELEMENTINDEX+2)+FLO_S(4)*FDIS0_ANG(ELEMENTINDEX+2)) + DIS_Y-DIS0_Y
                    X_INDEX_ARR(I,2) = FLO_S(1)*FDIS_Y(ELEMENTINDEX+1)+FLO_S(2)*FDIS_ANG(ELEMENTINDEX+1)+ &
                        & FLO_S(3)*FDIS_Y(ELEMENTINDEX+2)+FLO_S(4)*FDIS_ANG(ELEMENTINDEX+2)
                ELSE
                    DELTX(J)=0.0D0
                    DELTY(J)=0.0D0
                END IF
            ELSE
                FLO_XI=(X0-(ELEMENTINDEX-1)*FL_PARAMS_L)/FL_PARAMS_L
                FLO_S=0.0D0
                CALL SHAPE_FUN(FLO_S,4,FLO_XI,0)

                IF( T>T_START ) THEN
                    DELTX(J)=1.0*(DCOS(DIS_ANG)-DCOS(DIS0_ANG))!两部分，一部分是由于圆柱运动造成的整体向前的位移，另一部分是由于纵向变形造成的X轴的收缩
                    DELTY(J)=FLO_S(1)*FDIS_Y(ELEMENTINDEX)+FLO_S(2)*FDIS_ANG(ELEMENTINDEX)+ &
                        & FLO_S(3)*FDIS_Y(ELEMENTINDEX+1)+FLO_S(4)*FDIS_ANG(ELEMENTINDEX+1) - &
                        & (FLO_S(1)*FDIS0_Y(ELEMENTINDEX)+FLO_S(2)*FDIS0_ANG(ELEMENTINDEX) + &
                        & FLO_S(3)*FDIS0_Y(ELEMENTINDEX+1)+FLO_S(4)*FDIS0_ANG(ELEMENTINDEX+1)) + DIS_Y - DIS0_Y
                    
                    X_INDEX_ARR(I,2) = FLO_S(1)*FDIS_Y(ELEMENTINDEX)+FLO_S(2)*FDIS_ANG(ELEMENTINDEX)+ &
                        & FLO_S(3)*FDIS_Y(ELEMENTINDEX+1)+FLO_S(4)*FDIS_ANG(ELEMENTINDEX+1)
                ELSE
                    DELTX(J)=0.0D0
                    DELTY(J)=0.0D0
                END IF
            END IF
            !先插值出来Y方向的位移再说
            !还需要角度THETA确定单元的旋转
            !!!!!!!!!!!!!!!!在这里，记录网格点前后的数据是有必要的，考虑创建一个单独的数组进行储存!!!!!!!
            !!!!!!暂时在这里没有考虑THETA角造成的旋转，以及纵向变形对横向变形的影响!!!!!!!!!!!

            X_INDEX_ARR(I,1) = X_ALE(J)    ! 存储的是节点的x坐标
            
        END DO

        DO I =1,NSBN_ALE
            xFtemp=X_INDEX_ARR(I,1)
            Do J = I+1,NSBN_ALE
                xxFtemp =X_INDEX_ARR(J,1)
                IF (xFtemp==xxFtemp) THEN
                    X_INDEX_ARR(J,1)=0
                END IF
            End DO
        End do

        NKindSN=0
        DO I =1,NSBN_ALE
            xFtemp=X_INDEX_ARR(I,1)
            IF(xFtemp /= 0) then
                NKindSN =  NKindSN + 1
            End if
        End do

        !PRINT*, NKindSN
        !PAUSE
        
        Allocate(SORTED_X_INDEX_ARR(NKindSN,3),AddiDis(NKindSN))
        SORTED_X_INDEX_ARR = 0.0D0 ! 1：原始的轴线x坐标，2：原始的轴线x坐标对应的
        AddiDis = 0.0 
        
        J = 1
        DO I = 1, NSBN_ALE
            IF (X_INDEX_ARR(I,1)/=0)THEN
                SORTED_X_INDEX_ARR(J,1)=X_INDEX_ARR(I,1)
                SORTED_X_INDEX_ARR(J,2)=X_INDEX_ARR(I,2)  
                J=J+1
            END IF
        END DO

        Call Sort_Ascend_Multi(SORTED_X_INDEX_ARR,NKindSN,3) ! 按照第一列的大小进行排列

        PLATE_DS = 0.0
        DO I = NKINDSN,2,-1
            FDELTA_Y = SORTED_X_INDEX_ARR(I,2) - SORTED_X_INDEX_ARR(I-1,2)
            PLATE_DS = X_Plate_Orig(I,1) - X_Plate_Orig(I-1,1)
            FA_THETA(I) = ASIN(FDELTA_Y / PLATE_DS)
            FDELTA_X =  (PLATE_DS * DCOS(FA_THETA(I)) - (X_CA_ORR(I,1) - X_CA_ORR(I-1,1)))
            SORTED_X_INDEX_ARR(I,3) = FDELTA_X
            
            !SORTED_X_INDEX_ARR(I,4) = FA_THETA(I) - FA_THETA0(I)
            
            !PRINT*, FDELTA_X
        END DO
        SORTED_X_INDEX_ARR(1,3) = 0.0
        SORTED_X_INDEX_ARR(1,4) = 0.0
    
        !DO I = 1, NSBN_ALE
        !    J = SBN_ALE(I)
        !    X0 = X_ALE(J)
        !    Y0 = Y_ALE(J)
        !    Index_X = findloc(SORTED_X_INDEX_ARR(:,1),X0,1)
        !    IF(Y0 >= 0.0) THEN
        !        DELTX(J)=DELTX(J) + SUM(SORTED_X_INDEX_ARR(1:Index_X,3)) - Y0 * (DSIN(FA_THETA(Index_X))-DSIN(FA_THETA0(Index_X)))
        !        DELTY(J)=DELTY(J) + Y0 * ((DCOS(FA_THETA(Index_X)) - 1.0)-(DCOS(FA_THETA0(Index_X)) - 1.0))
        !        
        !    ELSE 
        !        DELTX(J)=DELTX(J) + SUM(SORTED_X_INDEX_ARR(1:Index_X,3)) - Y0 * (DSIN(FA_THETA(Index_X))-DSIN(FA_THETA0(Index_X)))
        !        DELTY(J)=DELTY(J) - Y0 * ((DCOS(FA_THETA(Index_X)) - 1.0)-(DCOS(FA_THETA0(Index_X)) - 1.0))
        !        
        !    END IF
        !
        !    AddiDis(Index_X) = DELTX(J)
        !
        !END DO
        
        X_CA_ORR(:,1) = X_CA_ORR(:,1) + AddiDis(:)
        
        FA_THETA0 = FA_THETA

    ELSE
        WRITE(*,*)"ERROR IN THE IROTA"
        PAUSE
        STOP
    END IF

    DO I=1,NLBN_ale
        num=num+1
        j=lbn_ale(I)
        DELTX(j)=0.0D0
        DELTY(j)=0.0D0
    END DO


    DO I=1,Ndbn_ale
        num=num+1
        j=dbn_ale(I)
        DELTX(j)=0.0D0
        DELTY(j)=0.0D0
    END DO


    DO I=1,Nubn_ale
        num=num+1
        j=ubn_ale(I)
        DELTX(j)=0.0D0
        DELTY(j)=0.0D0
    END DO


    do  i=1,nrbn_ale
        num=num+1
        j=rbn_ale(I)
        DELTX(j)=0.0D0
        DELTY(j)=0.0D0
    end do


    if(num/=nbn_ale) then
        print*, ' wrong',num,nbn_ale
        PAUSE
        stop
    end if

    FA_SPR1  = 320000
    FA_SPR71 = FA_SPR1 / 10
    FA_SPR72 = FA_SPR71 / 10
    FA_SPR73 = FA_SPR72 / 5
    FA_SPR74 = FA_SPR73 / 3
    FA_SPR75 = FA_SPR74 / 2

    Types_NodeTypes(1,1) = 71
    Types_NodeTypes(2,1) = 72
    Types_NodeTypes(3,1) = 73
    Types_NodeTypes(4,1) = 74
    Types_NodeTypes(5,1) = 75
    Types_NodeTypes(6,1) = 0
    
    Types_NodeTypes(:,2) = 0
    
    DO N = 1,NN_ALE
        Select case(Nodetype(N))
        CASE(71)
            Types_NodeTypes(1,2) = Types_NodeTypes(1,2)+1
        CASE(72)
            Types_NodeTypes(2,2) = Types_NodeTypes(2,2)+1
        CASE(73)
            Types_NodeTypes(3,2) = Types_NodeTypes(3,2)+1
        CASE(74)
            Types_NodeTypes(4,2) = Types_NodeTypes(4,2)+1
        CASE(75)
            Types_NodeTypes(5,2) = Types_NodeTypes(5,2)+1
        CASE(0)
            Types_NodeTypes(6,2) = Types_NodeTypes(6,2)+1
        END SELECT
    END DO
    
    ALLOCATE(NODE_INDEX_TYPES(SUM(Types_NodeTypes(:,2))))
    
    K = 1
    DO FNUM = 1,size(Types_NodeTypes(:,1))
        DO N = 1,NN_ALE
            IF (Nodetype(N) == Types_NodeTypes(FNUM,1)) THEN
                NODE_INDEX_TYPES(K) = N
                K = K +1
            END IF
        END DO
    END DO

    IF(K/=1+SUM(Types_NodeTypes(:,2)))THEN
        PRINT*,'NUMBER DO NOT EQUAL'
    END IF
    

    IF( T>T_START ) THEN
        DO FNUM = 1,size(Types_NodeTypes(:,1))
            DO N=1,ITER_NUM                        !以固定迭代次数计算
                FIJK = Types_NodeTypes(FNUM,2)
                DO  FIIJK = SUM(Types_NodeTypes(1:FNUM,2))- Types_NodeTypes(FNUM,2)+1,SUM(Types_NodeTypes(1:FNUM,2))                  !对ALE区域的所有坐标进行计算
                    I = NODE_INDEX_TYPES(FIIJK)
                    KPLUS=0.0                      !刚度赋初值
                    NUM_LMM=0
                    NUM_LIU=0
                    DO  J=1,MAXNUM_NODE_NODE                  !与本节点相连的最多节点个数
                        K=NODE_NODE_LMM(I,J)       !节点周围的节点，有哪几个
                        IF(K>0)THEN               !除去重合的点
                            NUM_LIU=NUM_LIU+1
                        END IF
                    END DO

                    DELTX00(I)=0.0D0
                    DELTY00(I)=0.0D0
                    DO J=1,MAXNUM_NODE_NODE
                        K=NODE_NODE_LMM(I,J)
                        IF(K>0)THEN                                                                     !除去重合的点。
                            Select case(Nodetype(i))
                            case(71)
                                Select case(Nodetype(K))
                                Case (1)
                                    FA_SPR = FA_SPR1
                                Case (71)
                                    FA_SPR = FA_SPR71
                                Case (72)
                                    FA_SPR = FA_SPR72
                                End Select
                            case(72)
                                Select case(Nodetype(K))
                                Case (1)
                                    FA_SPR = FA_SPR71
                                Case (71)
                                    FA_SPR = FA_SPR72
                                Case (72)
                                    FA_SPR = FA_SPR72
                                Case (73)
                                    FA_SPR = FA_SPR73
                                End Select
                            case(73)
                                Select case(Nodetype(K))
                                Case (71)
                                    FA_SPR = FA_SPR72
                                Case (72)
                                    FA_SPR = FA_SPR73
                                Case (73)
                                    FA_SPR = FA_SPR73
                                Case (74)
                                    FA_SPR = FA_SPR73
                                End Select
                            case(74)
                                Select case(Nodetype(K))
                                Case (72)
                                    FA_SPR = FA_SPR73
                                Case (73)
                                    FA_SPR = FA_SPR73
                                Case (74)
                                    FA_SPR = FA_SPR74
                                Case (75)
                                    FA_SPR = FA_SPR75
                                End Select
                                FA_SPR = 500D0
                            case(75)
                                Select case(Nodetype(K))
                                Case (73)
                                    FA_SPR = FA_SPR74
                                Case (74)
                                    FA_SPR = FA_SPR75
                                Case (75)
                                    FA_SPR = FA_SPR75
                                Case (0)
                                    FA_SPR = FA_SPR75
                                End Select
                            case(0)
                                Select case(Nodetype(K))
                                Case (74)
                                    FA_SPR = FA_SPR74
                                Case (75)
                                    FA_SPR = FA_SPR75
                                Case (0)
                                    FA_SPR = 1.0d0
                                End Select
                            End Select

                            NUM_LMM=NUM_LMM+1
                            KIJ(NUM_LMM)=FA_SPR*((X0_ALE(I)-X0_ALE(K))**2+(Y0_ALE(I)-Y0_ALE(K))**2)**FUSAI_SPR

                            KPLUS=KPLUS+KIJ(NUM_LMM)
                            DELTX00(I)=DELTX00(I)+KIJ(NUM_LMM)*DELTX(K)
                            DELTY00(I)=DELTY00(I)+KIJ(NUM_LMM)*DELTY(K)

                            !IF ((NodeType(I)==78).OR.(NodeType(K)==78)) Then
                            !    !FA_SPR=20.0d0
                            !    NUM_LMM=NUM_LMM+1
                            !
                            !    !InternalAngel = Max(Dabs(NODE_NODE_LMM_Angle(I,J,1)-pi/2),&
                            !    !    &Dabs(NODE_NODE_LMM_Angle(I,J,2)-pi/2),&
                            !    !    &Dabs(NODE_NODE_LMM_Angle(I,J,3)-pi/2),&
                            !    !    &Dabs(NODE_NODE_LMM_Angle(I,J,4)-pi/2))
                            !    !
                            !    !!Write(*,*) InternalAngel,NODE_NODE_LMM_Angle(I,J,:)
                            !    !!Pause
                            !    !
                            !    !KIJ(NUM_LMM)=25*(1/min(SELE_Iter(NODE_NODE_LMM_ELE(I,J,1)),SELE_Iter(NODE_NODE_LMM_ELE(I,J,2)))&
                            !    !    & + 1 / Dsin(pi/2-InternalAngel)**2 &
                            !    !    &+((X0_ALE(I)-X0_ALE(K))**2+(Y0_ALE(I)-Y0_ALE(K))**2)**FUSAI_SPR)
                            !    !
                            !    !write(*,*) KIJ(NUM_LMM),i,k
                            !
                            !    KIJ(NUM_LMM)=FA_SPR*((X0_ALE(I)-X0_ALE(K))**2+(Y0_ALE(I)-Y0_ALE(K))**2)**FUSAI_SPR
                            !
                            !    KPLUS=KPLUS+KIJ(NUM_LMM)                                                 !总的刚度
                            !    DELTX00(I)=DELTX00(I)+KIJ(NUM_LMM)*DELTX(K)                              !本节点的位移等于与本节点相连的节点的刚度乘以其上一步的位移
                            !    DELTY00(I)=DELTY00(I)+KIJ(NUM_LMM)*DELTY(K)
                            !Else
                            !    NUM_LMM=NUM_LMM+1
                            !    !InternalAngel = Max(Dabs(NODE_NODE_LMM_Angle(I,J,1)-pi/2),&
                            !    !    &Dabs(NODE_NODE_LMM_Angle(I,J,2)-pi/2),&
                            !    !    &Dabs(NODE_NODE_LMM_Angle(I,J,3)-pi/2),&
                            !    !    &Dabs(NODE_NODE_LMM_Angle(I,J,4)-pi/2))
                            !    !KIJ(NUM_LMM)=5.0*(1/min(SELE_Iter(NODE_NODE_LMM_ELE(I,J,1)),SELE_Iter(NODE_NODE_LMM_ELE(I,J,2)))&
                            !    !    & + 1 / Dsin(pi/2-InternalAngel)**2 &
                            !    !    &+((X0_ALE(I)-X0_ALE(K))**2+(Y0_ALE(I)-Y0_ALE(K))**2)**FUSAI_SPR)
                            !    !
                            !    KIJ(NUM_LMM)=FA_SPR*((X0_ALE(I)-X0_ALE(K))**2+(Y0_ALE(I)-Y0_ALE(K))**2)**FUSAI_SPR
                            !    !write(*,*) KIJ(NUM_LMM),i,k
                            !
                            !    KPLUS=KPLUS+KIJ(NUM_LMM)                                                 !总的刚度
                            !    DELTX00(I)=DELTX00(I)+KIJ(NUM_LMM)*DELTX(K)                              !本节点的位移等于与本节点相连的节点的刚度乘以其上一步的位移
                            !    DELTY00(I)=DELTY00(I)+KIJ(NUM_LMM)*DELTY(K)
                            !End if
                        END IF
                    END DO
                    DELTX(I)=DELTX00(I)
                    DELTY(I)=DELTY00(I)
                    DELTX(I)=DELTX(I)/KPLUS                                                                !节点的位移等于总位移除以总刚度
                    DELTY(I)=DELTY(I)/KPLUS
                ENDDO
            ENDDO
        End Do
    ELSE
        DELTX=0.0D0
        DELTY=0.0D0
    END IF

100 DO  I=1,NN_ALE
        X1_ALE(I)=X0_ALE(i)+DELTX(i)
        Y1_ALE(I)=Y0_ALE(i)+DELTY(i)
    ENDDO

    write(*,*)

    DEALLOCATE(DELTX00,DELTY00)

    END SUBROUTINE SPRING_METHOD

    SUBROUTINE SPRING_METHOD_Test
    Use FlexibleOnly
    USE CONSTANT
    USE VARIABLE
    IMPLICIT NONE

    INTEGER::NUM_LMM,NUM_LIU,ElementIndex
    REAL(8)::ttt0,ttt1,Flo_xi,Flo_S(4)
    REAL(8)::KIJ(maxnum_node_node),KPLUS,PI,X0,Y0
    INTEGER::I,J,N,k,num_it,num
    REAL(8)::SMALL=1.0d-8,theta1,theta2,theta3,R00,dd
    real(8)::spring_ring_width


    REAL(8),ALLOCATABLE:: DELTX00(:),DELTY00(:)
    ALLOCATE(deltx00(nn_ale),delty00(nn_ale))

    spring_ring_width=spring_out-spring_in

    pi=dacos(-1.0d0)
    num=0


    if(IROTA==2012)then
        DO I=1,NPBN_ALE
            NUM=NUM+1
            J=PBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)

            DD=DSQRT(X0**2+Y0**2)
            THETA2=DASIN(Y0/DD)
            IF (X0<0.0D0)   THEN
                THETA2=PI-THETA2
            END IF

            IF( T>T_START ) THEN
                DELTX(J)=DD*DCOS(DIS_ANG+THETA2)-DD*DCOS(THETA2+DIS0_ANG)
                DELTY(J)=DD*DSIN(DIS_ANG+THETA2)-DD*DSIN(THETA2+DIS0_ANG)+DIS_Y-DIS0_Y
            ELSE
                DELTX(J)=0.0D0
                DELTY(J)=0.0D0
            END IF

            write(*,*) DELTY(J)
        END DO


        DO I=1,NSBN_ALE

            !首先要把计算出来的位移进行插值到网格节点上
            NUM=NUM+1
            J=SBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)

            X0=X0-0.5D0
            ElementIndex=FLOOR(X0/FL_Params_L)

            IF (ElementIndex<FL_Params_ne)THEN
                Flo_xi=(X0-ElementIndex*FL_Params_L)/FL_Params_L  !局部坐标
                Flo_S=0.0D0
                call Shape_Fun(Flo_S,4,Flo_xi,0) !S现在是4*1，转置后，1*4，乘以4*1的位移阵，就得到了节点的挠度

                If (ElementIndex==0) Then
                    ElementIndex = 1
                    IF( T>T_START ) THEN
                        DELTX(J)=1.0*(dcos(Dis_ang)-dcos(Dis0_ang))!两部分，一部分是由于圆柱运动造成的整体向前的位移，另一部分是由于纵向变形造成的X轴的收缩
                        DELTY(J)=Flo_S(1)*Fdis_Y(ElementIndex)+Flo_S(2)*Fdis_ang(ElementIndex)+ &
                            & Flo_S(3)*Fdis_Y(ElementIndex+1)+Flo_S(4)*Fdis_ang(ElementIndex+1) - &
                            & (Flo_S(1)*Fdis0_Y(ElementIndex)+Flo_S(2)*Fdis0_ang(ElementIndex) + &
                            & Flo_S(3)*Fdis0_Y(ElementIndex+1)+Flo_S(4)*Fdis0_ang(ElementIndex+1)) + DIS_y-Dis0_y
                    Else
                        DELTX(J)=0.0d0
                        DELTY(J)=0.0d0
                    End if
                End if

            ELSE
                Flo_xi=(X0-(ElementIndex-1)*FL_Params_L)/FL_Params_L
                Flo_S=0.0D0
                call Shape_Fun(Flo_S,4,Flo_xi,0)

                IF( T>T_START ) THEN
                    DELTX(J)=1.0*(dcos(Dis_ang)-dcos(Dis0_ang))!两部分，一部分是由于圆柱运动造成的整体向前的位移，另一部分是由于纵向变形造成的X轴的收缩
                    DELTY(J)=Flo_S(1)*Fdis_Y(ElementIndex-1)+Flo_S(2)*Fdis_ang(ElementIndex-1)+ &
                        & Flo_S(3)*Fdis_Y(ElementIndex)+Flo_S(4)*Fdis_ang(ElementIndex) - &
                        & (Flo_S(1)*Fdis0_Y(ElementIndex-1)+Flo_S(2)*Fdis0_ang(ElementIndex-1) + &
                        & Flo_S(3)*Fdis0_Y(ElementIndex)+Flo_S(4)*Fdis0_ang(ElementIndex)) + DIS_y-Dis0_y

                Else
                    DELTX(J)=0.0d0
                    DELTY(J)=0.0d0
                End if
            End if

            !先插值出来y方向的位移再说
            !还需要角度theta确定单元的旋转

            !!!!!!!!!!!!!!!!在这里，记录网格点前后的数据是有必要的，考虑创建一个单独的数组进行储存!!!!!!!
            !!!!!!暂时在这里没有考虑theta角造成的旋转，以及纵向变形对横向变形的影响!!!!!!!!!!!

            write(*,*) I,DELTY(J)

        END DO

    else

        write(*,*)"error in the irota"
        PAUSE
        stop

    END IF


    DO I=1,NLBN_ale
        num=num+1
        j=lbn_ale(I)
        DELTX(j)=0.0D0
        DELTY(j)=0.0D0
    END DO


    DO I=1,Ndbn_ale
        num=num+1
        j=dbn_ale(I)
        DELTX(j)=0.0D0
        DELTY(j)=0.0D0
    END DO


    DO I=1,Nubn_ale
        num=num+1
        j=ubn_ale(I)
        DELTX(j)=0.0D0
        DELTY(j)=0.0D0
    END DO


    do  i=1,nrbn_ale
        num=num+1
        j=rbn_ale(I)
        DELTX(j)=0.0D0
        DELTY(j)=0.0D0
    end do


    if(num/=nbn_ale) then
        print*, ' wrong',num,nbn_ale
        PAUSE
        stop
    end if





    DEALLOCATE(DELTX00,DELTY00)

    END SUBROUTINE
