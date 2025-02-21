    SubRoutine RBFMeshUpdate
    Use FlexibleOnly
    Use ANCF
    Use Variable
    Use CONSTANT
    use lapack95
    use blas95
    Use Local
    Use Omp_lib

    Implicit NONE
    Integer::I,J,K,II,JJ,KK,fksSize,Parallel_Switch
    Real(8)::X0,Y0,SCL,FOmega,PI,eee,sss

    Integer,Allocatable::ipiv(:)

    Real(8),ALLOCATABLE::RBF_Delta(:,:),RBF_POS_Delta(:,:)
    REAL(8)::Lo_xi,Lo_S(4),theta1,theta2,DD
    REAL(8)::xxFtemp,xFtemp
    Integer::ElementIndex,Is,Ie,NUM
    REAL(8)::Lo_SI(AFP%Dim,AFP%NumE),Lo_SI_Tran(AFP%NumE,AFP%Dim),Coord(AFP%Dim),eele(AFP%NumE)
    REAL(8)::Delta_x,Delta_Y,FNX,FNY,RBFDelta_X,RBFDelta_Y,SR
    REAL(8)::CPU_START,CPU_FINISH
    Character(8)::ZoneName

    !Main Function: Calculate Update X0_ALE to X1_ALE
    fksSize = NBN_ALE_KS
    Allocate(RBF_Delta(Mesh_Update_Num,2),RBF_POS_Delta(fksSize,2))

    
    PI=DACOS(-1.0D0)

    RBF_Delta       = 0.0
    RBF_POS_Delta   = 0.0

    !Step1 define Position of Boundary Points
    If (IROTA == 5107) Then
        FOmega = 1.1
        SCL = 0.5
        DO I = 1,NPBN_ALE
            IF( T>T_START ) THEN
                J = Pbn_ALE(I)
                X1_ALE(J) = X0_ALE(J)  !X0_ALE(J)
                Y1_ALE(j) = Y0_ALE(j)  !+(SIN(T) - SIN(T-DT))
            Else
                X1_ALE(J) = X0_ALE(J)
                Y1_ALE(j) = Y0_ALE(J)
            End If
        End Do
        Do I = 1,NSBN_ALE
            IF( T>T_START ) THEN
                J = Sbn_ALE(I)
                X1_ALE(J)=X0_ALE(J)
                X0 = X_ALE(J) -0.5
                !y0 = (Dsin(T*SCL)-Dsin((T-DT)*SCL))*((Dcosh(FOmega*X0))-DCOS(FOmega*X0)+(Dcos(FOmega)+Dcosh(FOmega))/(Dsin(FOmega)+Dsinh(FOmega))*(Dsin(FOmega*X0)-Dsinh(FOmega*X0)))
                Y1_ALE(J)=Y_ALE(J) + (Dsin(T*3.1415926))*SCL*((Dcosh(FOmega*X0))-DCOS(FOmega*X0)+(Dcos(FOmega * 3.5D0)+Dcosh(FOmega * 3.5D0)) &
                    & /(Dsin(FOmega * 3.5D0)+Dsinh(FOmega * 3.5D0))*(Dsin(FOmega*X0)-Dsinh(FOmega*X0)))
                !write(*,*) (Dsin(T/2))*SCL*((Dcosh(FOmega*X0))-DCOS(FOmega*X0)+(Dcos(FOmega * 3.5D0)+Dcosh(FOmega * 3.5D0)) &
                !& /(Dsin(FOmega * 3.5D0)+Dsinh(FOmega * 3.5D0))*(Dsin(FOmega*X0)-Dsinh(FOmega*X0)))
                if ((abs(X0-3.5d0)<0.0001).and.(abs(Y_ALE(J)-0.1d0)<0.0001))Then
                    write(*,*) Y1_ALE(J)
                End if
            ELSE
                X1_ALE(J) = X0_ALE(J)
                Y1_ALE(j) = Y0_ALE(J)
            End if
        End Do

    ELSE IF((IROTA==510).OR.(IROTA==5102).OR.(Irota==5101).OR.(IROTA==520))THEN
        IF((IROTA==510).OR.(IROTA==5101))Then
            DO I=1,NPBN_ALE !Cylinder Remain unmoved
                J=PBN_ALE(I)
                IF( T>T_START ) THEN
                    X1_ALE(J)= X0_ALE(J)
                    Y1_ALE(J)= Y0_ALE(J)
                ELSE
                    X1_ALE(J)= X0_ALE(J)
                    Y1_ALE(J)= Y0_ALE(J)
                END IF
            END DO
        Else IF ((Irota==5102).OR.(IROTA==520)) Then
            DO I=1,NPBN_ALE !Cylinder Oscillate
                J=PBN_ALE(I)
                IF( T>T_START ) THEN
                    X1_ALE(J)= X_ALE(J)
                    Y1_ALE(J)= Y_ALE(J) + AF_e(2)
                ELSE
                    X1_ALE(J)= X_ALE(J)
                    Y1_ALE(J)= Y_ALE(J)
                END IF
            END DO
        End If

        Sec_Rota0 = Sec_Rota
        Delta_x0510  = Delta_x510
        Delta_y0510  = Delta_y510

        DO I=1,NSBN_ALE
            NUM=NUM + 1
            J=SBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)
            X0=X0-0.5D0
            ELEMENTINDEX=FLOOR(X0/AFP%LE)
            eele = 0.0
            IF (ELEMENTINDEX<AFP%NE)THEN
                Is = AFP%NumD*ELEMENTINDEX  + 1
                Ie = AFP%NumD*(ELEMENTINDEX+2)
                eele = AF_E(Is:Ie)! 12 Coordinate of the Element
                Lo_xi=(X0-ELEMENTINDEX*AFP%LE)/AFP%LE
                Lo_S=0.0D0
                Call ANCF_Shape_Fun(Lo_S,Lo_SI,Lo_xi,0)
                IF( T>T_START ) THEN
                    Coord = 0.0
                    Call S_e(Lo_S,eele,Coord)
                    X1_ALE(J)= 0.5D0 + Coord(1)
                    Y1_ALE(J)= Y0 + Coord(2)
                ELSE
                    X1_ALE(J)=X0_ALE(J)
                    Y1_ALE(J)=Y0_ALE(J)
                END IF

                Lo_S = 0.0D0
                Call ANCF_Shape_Fun(Lo_S,Lo_SI,Lo_xi,1)
                IF( T>T_START ) THEN
                    Call S_e(Lo_S,eele,Sec_Rota(I,:))
                    Sec_Rota(I,:) = Sec_Rota(I,:) / AFP%LE
                Else
                    Sec_Rota(I,:) = 0.0D0
                End If

            ELSE
                Is = AFP%NumD*ELEMENTINDEX - AFP%NumD + 1
                Ie = AFP%NumD*(ELEMENTINDEX + 1 )
                eele = AF_E(Is:Ie)
                Lo_xi=(X0-(ELEMENTINDEX-1)*AFP%LE)/AFP%LE
                Lo_S = 0.0D0
                Call ANCF_Shape_Fun(Lo_S,Lo_SI,Lo_xi,0)
                IF( T>T_START ) THEN
                    Coord = 0.0
                    Call S_e(Lo_S,eele,Coord)
                    X1_ALE(J)= 0.5D0 + Coord(1)
                    Y1_ALE(J)= Y0 + Coord(2)
                ELSE
                    X1_ALE(J)=X0_ALE(J)
                    Y1_ALE(J)=Y0_ALE(J)
                END IF

                Lo_S = 0.0D0
                Call ANCF_Shape_Fun(Lo_S,Lo_SI,Lo_xi,1)
                IF( T>T_START ) THEN
                    Call S_e(Lo_S,eele,Sec_Rota(I,:))
                    Sec_Rota(I,:) = Sec_Rota(I,:) / AFP%LE
                Else
                    Sec_Rota(I,:) = 0.0D0
                End If
            END IF
        END DO

        ! 截面旋转补正
        KK = MaxVal(PlatePointIndex(:,1))
        DO I=1,NSBN_ALE
            !J = SBN_ALE(I)
            K = PlatePointIndex(I,1)
            IF (K/=KK) Then
                X_Neutal(K,1) = 0.5D0 * Sum(X1_Ale(PlatePointIndex_UDF(K,:)))
                Y_Neutal(K,1) = 0.5D0 * Sum(Y1_Ale(PlatePointIndex_UDF(K,:)))
            Else
                X_Neutal(K,1) = 1.0D0 / Size(PlatePointIndex_Tail(:,1)) * Sum(X1_Ale(PlatePointIndex_Tail(:,1)))
                Y_Neutal(K,1) = 1.0D0 / Size(PlatePointIndex_Tail(:,1)) * Sum(Y1_Ale(PlatePointIndex_Tail(:,1)))
            End IF
        END DO

        DO I=1,NSBN_ALE
            J = SBN_ALE(I)
            K = PlatePointIndex(I,1)
            Y0 = Y_ALE(J)
            If (K > 1) Then
                Call KNXNY(X_Neutal(K,1),Y_Neutal(K,1),X_Neutal(K-1,1),Y_Neutal(K-1,1),FNX,FNY)!永远后面点减前面点
                Delta_x510(K) = - Y0 * FNX
                Delta_y510(K) = - Y0 * ( 1.0D0 - FNY)
            Else
                Delta_x510(K)= 0.0D0
                Delta_y510(K)= 0.0D0
            End If
            X1_ALE(J) = X1_ALE(J) + Delta_x510(K)
            Y1_ALE(J) = Y1_ALE(J) + Delta_y510(K)
        END DO


    ELSE IF(IROTA==200.OR.IROTA==201)THEN            !JY!2DOF
        DO I=1,NSBN_ALE
            NUM=NUM+1
            J=SBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)
            IF( T>T_START ) THEN
                X1_ALE(J)=X0_ALE(J) + X0*(DCOS(DIS_ANG)-DCOS(DIS0_ANG))-Y0*(DSIN(DIS_ANG)-DSIN(DIS0_ANG))
                Y1_ALE(J)=Y0_ALE(J) + X0*(DSIN(DIS_ANG)-DSIN(DIS0_ANG))+Y0*(DCOS(DIS_ANG)-DCOS(DIS0_ANG))+DIS_Y-DIS0_Y
            ELSE
                X1_ALE(J)=X0_ALE(J)
                Y1_ALE(J)=Y0_ALE(J)
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
                X1_ALE(J)=X0_ALE(J) + DD*DCOS(DIS_ANG+THETA2)-DD*DCOS(THETA2+DIS0_ANG)
                Y1_ALE(J)=Y0_ALE(J) + DD*DSIN(DIS_ANG+THETA2)-DD*DSIN(THETA2+DIS0_ANG)+DIS_Y-DIS0_Y
            ELSE
                X1_ALE(J)=X0_ALE(J)
                Y1_ALE(J)=Y0_ALE(J)
            END IF
        END DO
    ELSE IF(IROTA==2201)THEN            
        DO I=1,NSBN_ALE       !DO
            NUM=NUM+1
            J=SBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)
            IF( T>T_START ) THEN
                X1_ALE(J)=X_ALE(J) 
                Y1_ALE(J)=Y_ALE(J) + DIS_Y
            ELSE
                X1_ALE(J)=X0_ALE(J)
                Y1_ALE(J)=Y0_ALE(J)
            END IF
        END DO
        DO I=1,NPBN_ALE      !UP
            NUM=NUM+1
            J=PBN_ALE(I)
            IF( T>T_START ) THEN
                X1_ALE(J)=X_ALE(J) 
                Y1_ALE(J)=Y_ALE(J) + DIS_Y_UP
            ELSE
                X1_ALE(J)=X0_ALE(J)
                Y1_ALE(J)=Y0_ALE(J)
            END IF
        END DO  
    ELSE IF(IROTA==2101)THEN
        DO I=1,NSBN_ALE       !DO
            NUM=NUM+1
            J=SBN_ALE(I)
            X0=X_ALE(J)
            Y0=Y_ALE(J)
            IF( T>T_START ) THEN
                X1_ALE(J)=X_ALE(J)
                Y1_ALE(J)=Y_ALE(J) + DIS_Y
            ELSE
                X1_ALE(J)=X0_ALE(J)
                Y1_ALE(J)=Y0_ALE(J)
            END IF
        END DO
        DO I=1,NPBN_ALE      !UP
            NUM=NUM+1
            J=PBN_ALE(I)
            IF( T>T_START ) THEN
                X1_ALE(J)=X_ALE(J) 
                Y1_ALE(J)=Y_ALE(J)
            ELSE
                X1_ALE(J)=X0_ALE(J)
                Y1_ALE(J)=Y0_ALE(J)
            END IF
        END DO  
    END IF

    DO I=1,NLBN_ALE
        J=LBN_ALE(I)
        X1_ALE(J)= X0_ALE(J)
        Y1_ALE(J)= Y0_ALE(J)
    END DO

    DO I=1,NDBN_ALE
        J=DBN_ALE(I)
        X1_ALE(J)= X0_ALE(J)
        Y1_ALE(J)= Y0_ALE(J)
    END DO

    DO I=1,NUBN_ALE
        J=UBN_ALE(I)
        X1_ALE(J)= X0_ALE(J)
        Y1_ALE(J)= Y0_ALE(J)
    END DO

    DO  I=1,NRBN_ALE
        J=RBN_ALE(I)
        X1_ALE(J)= X0_ALE(J)
        Y1_ALE(J)= Y0_ALE(J)
    END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Step 2 Calculte w!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Droped Now No Need!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Step 3 Calculate X1_Ale Of Inner Filed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    RBF_POS_Delta(:,1)  = X1_ALE(BN_ALE_KS(:))-X_ALE(BN_ALE_KS(:))
    RBF_POS_Delta(:,2)  = Y1_ALE(BN_ALE_KS(:))-Y_ALE(BN_ALE_KS(:))

    Call gemm(Left_Constant,RBF_POS_Delta,RBF_Delta)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MKL3

    If((RE_START==1).or.(Nstep==1))Then
        X1_Ale = X0_Ale
        Y1_Ale = Y0_Ale
    End IF

    X1_Ale(Mesh_Update_Index(1:Mesh_Update_Num)) = X_Ale(Mesh_Update_Index(1:Mesh_Update_Num)) + RBF_Delta(:,1)
    Y1_Ale(Mesh_Update_Index(1:Mesh_Update_Num)) = Y_Ale(Mesh_Update_Index(1:Mesh_Update_Num)) + RBF_Delta(:,2)

    DO I=1,NDBN_ALE
        J=DBN_ALE(I)
        X1_ALE(J)= X0_ALE(J)
        Y1_ALE(J)= Y0_ALE(J)
    END DO

    DO I=1,NUBN_ALE
        J=UBN_ALE(I)
        X1_ALE(J)= X0_ALE(J)
        Y1_ALE(J)= Y0_ALE(J)
    END DO

    DO  I=1,NRBN_ALE
        J=RBN_ALE(I)
        X1_ALE(J)= X0_ALE(J)
        Y1_ALE(J)= Y0_ALE(J)
    END DO


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Output!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    End SubRoutine RBFMeshUpdate


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    SubRoutine GatherSp
    !20211211 把ALE区域边界的位移也考虑进来，对于较小计算域的算例重要
    Use VARIABLE
    Use FlexibleOnly
    Use Local

    Implicit None
    Integer:: I,II,J,JJ,K,Num,NBN_ALE_KS0
    Real(8)::X0,Y0,XX0,YY0
    Real(8),Allocatable::temp(:,:),TempTail(:,:)
    Integer,Allocatable::SPBN_ALE_KS0(:),BN_ALE_KS0(:),BoundType0(:)

    NSPBN_ALE_KS = NSBN_ALE+NPBN_ALE
    Allocate(SPBN_ALE_KS(NSPBN_ALE_KS),temp(NSPBN_ALE_KS,2),SPBN_ALE_KS0(NSPBN_ALE_KS))

    Num = 1
    Do I = 1,NPBN_ALE
        SPBN_ALE_KS(Num) = PBN_ALE(I)
        Num = Num +1
    End Do

    Do I = 1,NSBN_ALE
        SPBN_ALE_KS(Num) = SBN_ALE(I)
        Num = Num +1
    End Do

    If((Num - 1)/=(NSBN_ALE+NPBN_ALE))Then
        Print*, 'Num Do Not Equal To N S P '
        Pause
        Stop
    End If

    Temp = 9999.00
    K = 0
    Do I = 1 ,NSPBN_ALE_KS
        J = SPBN_ALE_KS(I)
        X0 = X_ALE(J)
        Y0 = Y_ALE(J)
        If (Y0>=0.0)then
            K = K + 1
            temp(I,1) = X0
            temp(I,2) = J
        End If
    End do

    Call Sort_Ascend_Multi(temp,NSPBN_ALE_KS,2)

    X0 = temp(K,1)
    K = 0
    JJ = 0
    Do I = NSPBN_ALE_KS,1,-1
        If (Temp(I,1)==X0)Then
            K = K + 1
            If ( I > JJ ) Then
                JJ = I
            End IF
        End IF
    End Do

    Allocate(TempTail(K,2))
    J = 0
    Do I = JJ,JJ- K + 1,-1
        J = J + 1
        TempTail(J,1) = Y_ALE(Temp(I,2))
        TempTail(J,2) = temp(I,2)!Index
    End Do

    Call Sort_Descend_Multi(TempTail,K,2)

    J = 0
    Do I = JJ- K + 1,JJ
        J = J + 1
        temp(I,1) = X_ALe(TempTail(J,2))
        temp(I,2) = TempTail(J,2)
    End Do

    SPBN_ALE_KS0 = 0
    k = 1
    Do I = 1,NSPBN_ALE_KS
        If (Temp(I,1)/=9999.00)Then
            SPBN_ALE_KS0(I) = Temp(I,2)
            K = K + 1
        End If
    End Do

    Temp = 9999.00
    Do I = 1 ,NSPBN_ALE_KS
        J = SPBN_ALE_KS(I)
        X0 = X_ALE(J)
        Y0 = Y_ALE(J)
        If (Y0<0.0)then
            temp(I,1) = X0
            temp(I,2) = J
        End If
    End do

    call Sort_Ascend_Multi(temp,NSPBN_ALE_KS,2)

    Do I = NSPBN_ALE_KS,1,-1
        If (Temp(I,1)/=9999.00)Then
            SPBN_ALE_KS0(K) = Temp(I,2)
            K = K +1
        End If

    End Do

    If (k-1/=NSPBN_ALE_KS) then
        print*,'aya'
    end if

    SPBN_ALE_KS = SPBN_ALE_KS0

    !!!!!!!!!!!!!!!!!!!!!!!!Add Boundary Dis!!!!!!!!!!!!!!!!!!!!!!!

    NBN_ALE_KS0 = NSBN_ALE+NPBN_ALE+NRBN_ALE!+NDBN_ALE+NLBN_ALE+NUBN_ALE
    Allocate(BN_ALE_KS0(NBN_ALE_KS0),BoundType0(NBN_ALE_KS0))

    Num = 1
    Do I = 1,NSPBN_ALE_KS
        BN_ALE_KS0(I) = SPBN_ALE_KS(I)
        BoundType0(I) = 1
        Num = Num +1
    End Do
    !
    !Do I = 1,NLBN_ALE
    !    BN_ALE_KS0(Num) = LBN_ALE(I)
    !    BoundType0(Num) = 2
    !    Num = Num +1
    !End Do
    !
    !Do I = 1,NDBN_ALE
    !    BN_ALE_KS0(Num) = DBN_ALE(I)
    !    BoundType0(Num) = 2
    !    Num = Num +1
    !End Do
    !
    Do I = 1,NRBN_ALE
        BN_ALE_KS0(Num) = RBN_ALE(I)
        BoundType0(Num) = 2
        Num = Num +1
    End Do
    !
    !Do I = 1,NUBN_ALE
    !    BN_ALE_KS0(Num) = UBN_ALE(I)
    !    BoundType0(Num) = 2
    !    Num = Num +1
    !End Do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Remove Repeated Values!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Num = 0
    DO I =1,NBN_ALE_KS0
        J = BN_ALE_KS0(I)
        If((J/=0))Then
            If(I<=NBN_ALE_KS0)Then
                X0 = X_ALE(J)
                Y0 = Y_ALE(J)
                Do II = I+1,NBN_ALE_KS0
                    JJ = BN_ALE_KS0(II)
                    If(JJ/=0)Then
                        XX0 = X_ALE(JJ)
                        YY0 = Y_ALE(JJ)
                        IF ((X0==XX0).And.(Y0==YY0)) Then
                            BN_ALE_KS0(II) = 0
                        End If
                    End IF
                End Do
            End IF
            Num = Num + 1
        End IF
    End Do

    NBN_ALE_KS = Num
    Allocate(BN_ALE_KS(NBN_ALE_KS),BoundType(NBN_ALE_KS))

    Num = 0
    DO I =1,NBN_ALE_KS0
        J = BN_ALE_KS0(I)
        K = BoundType0(I)
        If(J/=0)Then
            Num = Num + 1
            BN_ALE_KS(Num) = J
            BoundType(num) = K
        End IF
    End Do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Remove Repeated Values!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DeAllocate(temp,SPBN_ALE_KS0,BN_ALE_KS0)
    End SubRoutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SubRoutine CalPlatePointIndex

    Use VARIABLE,Only:X_Ale,Y_Ale,NSBN_ALE,SBN_ALE,X0_Ale,Y0_Ale
    Use ANCF
    Use Local

    Implicit None
    !!!!Local!!!!
    Integer::I,J,JJ,K,KK,Num,IndexTemp(NSBN_ALE,1)
    !!!!Local!!!!

    Allocate(PlatePointIndex(NSBN_ALE,2))
    PlatePointIndex = 0.0D0
    IndexTemp = 0

    Do I = 1,NSBN_ALE
        J = SBN_AlE(I)
        PlatePointIndex(I,1) = X_ALE(J)
        PlatePointIndex(I,2) = I!每个点对应的index
    End Do

    Call Sort_Ascend_Multi_V2(PlatePointIndex,NSBN_ALE,2,1)

    Num = 1
    Do I = 2,NSBN_ALE
        If (PlatePointIndex(I-1,1) < PlatePointIndex(I,1)) Then
            Num = Num + 1
            IndexTemp(I,1) = Num!“单元”号
        ElseIf(PlatePointIndex(I-1,1) == PlatePointIndex(I,1)) Then
            IndexTemp(I,1) = Num
        ElseIf(PlatePointIndex(I-1,1) > PlatePointIndex(I,1)) Then
            Print*,'Wrong! Code#666102'
            Pause
            Stop
        End If
    End Do
    IndexTemp(1,1) = 1
    PlatePointIndex(:,1) = IndexTemp(:,1)
    Allocate(PlatePointIndex_UDF(Num-1,2))
    K  = MaxVal(PlatePointIndex(:,1))

    Call Sort_Ascend_Multi_V2(PlatePointIndex,NSBN_ALE,2,2)

    Num = 1
    JJ = 1
    KK = 0
    PlatePointIndex_UDF(:,1) = 123456789
    Do I = 1,NSBN_ALE
        J = SBN_ALE(I)
        JJ = PlatePointIndex(I,1)!Ele Index
        If ( JJ /= K )Then
            If (PlatePointIndex_UDF(JJ,1)==123456789) Then
                PlatePointIndex_UDF(JJ,1) = J
            Else
                PlatePointIndex_UDF(JJ,2) = J
            End If
        Else
            KK = KK + 1
        End If
    End Do

    Allocate(PlatePointIndex_Tail(KK,1))

    Num = 0
    Do I = 1,NSBN_ALE
        IF(PlatePointIndex(I,1) == K )Then
            Num = Num + 1
            PlatePointIndex_Tail(Num,1) = SBN_ALE(I)
        End If
    End Do

    I = MaxVal(PlatePointIndex(:,1))
    Allocate(X_Neutal(I,1))
    Allocate(Y_Neutal(I,1))
    X_Neutal=0.0D0
    Y_Neutal=0.0D0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!在此定义

    KK = MaxVal(PlatePointIndex(:,1))
    DO I=1,NSBN_ALE
        !J = SBN_ALE(I)
        K = PlatePointIndex(I,1)
        IF (K/=KK) Then
            X_Neutal(K,1) = 0.5D0 * Sum(X0_Ale(PlatePointIndex_UDF(K,:)))
            Y_Neutal(K,1) = 0.5D0 * Sum(Y0_Ale(PlatePointIndex_UDF(K,:)))
        Else
            X_Neutal(K,1) = 1.0D0 / Size(PlatePointIndex_Tail(:,1)) * Sum(X0_Ale(PlatePointIndex_Tail(:,1)))
            Y_Neutal(K,1) = 1.0D0 / Size(PlatePointIndex_Tail(:,1)) * Sum(Y0_Ale(PlatePointIndex_Tail(:,1)))
        End IF
    END DO

    I = MaxVal(PlatePointIndex(:,1))
    Allocate(Delta_x510(I),delta_y510(I),Delta_x0510(I),delta_y0510(I))
    Delta_x510   = 0.0D0
    delta_y510   = 0.0D0
    Delta_x0510  = 0.0D0
    delta_y0510  = 0.0D0

    End SubRoutine


    SubRoutine RBF_Greedy
    Use FlexibleOnly,Only:NBN_ALE_KS,BN_ALE_KS
    Use Variable
    Use ANCF
    use Lapack95
    Use Blas95
    Use Local

    Integer::I,Num,II,JJ,MorePoint
    Real(8)::IFStop_Sel,sss,eee,aaad
    Real(8),Allocatable::Red_Index(:),Err_Red(:,:)
    Integer,Allocatable::TempIndex(:),MaxIndex(:),AddedIndex(:)

    Reduction_Size = 25
    MorePoint = 25

    IF(Allocated(Err_X_Red)) DEAllocate(Err_X_Red)
    IF(Allocated(Err_Y_Red)) DEAllocate(Err_Y_Red)
    IF(Allocated(Err_Red)) DEAllocate(Err_Red)
    IF(Allocated(BN_ALE_KS_Red)) DEAllocate(BN_ALE_KS_Red)

    Allocate(Err_X_Red(NBN_ALE_KS,2))
    Allocate(Err_Y_Red(NBN_ALE_KS,2))
    Allocate(Err_Red(NBN_ALE_KS,2))
    Allocate(BN_ALE_KS_Red(Reduction_Size))

    Allocate(MaxIndex(2))
    Allocate(Red_Index(Reduction_Size),TempIndex(Reduction_Size))

    AddedIndex = 0

    Do I = 1,Reduction_Size
        Call RANDOM_NUMBER(Red_Index(I))
        Red_Index(I) = Floor(Red_Index(I)*(NBN_ALE_KS-1)) + 1!避免0
        !BN_ALE_KS_Red(I) = BN_ALE_KS(Red_Index(I))
        BN_ALE_KS_Red(I) = Red_Index(I)
    End Do

    DO I =1,Reduction_Size-1
        J = BN_ALE_KS_Red(I)
        If((J/=0))Then
            Do II = I+1,Reduction_Size
                JJ = BN_ALE_KS_Red(II)
                If(JJ/=J)Then
                    Continue
                Else
                    BN_ALE_KS_Red(II) = 0
                End IF
            End Do
        End IF
    End Do

    Num = 0
    DO I =1,Reduction_Size
        J = BN_ALE_KS_Red(I)
        If(J/=0)Then
            Num = Num + 1
        End IF
    End Do

    DeAllocate(TempIndex)
    Allocate(TempIndex(Num))
    Num = 0
    DO I =1,Reduction_Size
        J = BN_ALE_KS_Red(I)
        If(J/=0)Then
            Num = Num + 1
            TempIndex(Num) = J
        End IF
    End Do
    DeAllocate(BN_ALE_KS_Red)
    Allocate(BN_ALE_KS_Red(Num))
    BN_ALE_KS_Red = TempIndex

    Call SelectMorePoint

    Err_Red(:,2) =(/(i,i=1,NBN_ALE_KS)/)
    Err_Red(:,1) =DSQRT( Err_X_Red(:,1)*Err_X_Red(:,1) + Err_Y_Red(:,1) * Err_Y_Red(:,1) )
    IFStop_Sel = Sum(Err_Red(:,1))/NBN_ALE_KS


    Allocate(AddedIndex(MorePoint))
    AddedIndex = 0
    Do I =1,MorePoint
        MaxIndex = MaxLoc(Err_Red,1)
        Err_Red(MaxIndex(1),1) = 0.0D0
        AddedIndex(I) = Err_Red(MaxIndex(1),2)
    End Do
    Num =0
    aaad = 0.0D0

    Do While(IFStop_Sel>1.0e-7)
        Num = Num + 1

        Reduction_Size = Size(BN_ALE_KS_Red) + MorePoint
        DeAllocate(BN_ALE_KS_Red)
        Allocate(BN_ALE_KS_Red(Reduction_Size))
        BN_ALE_KS_Red(1:Reduction_Size-MorePoint) = TempIndex(:)
        BN_ALE_KS_Red(Reduction_Size-MorePoint+1:Reduction_Size) = AddedIndex(1:MorePoint)
        !
        Num = 0
        DO I =1,Reduction_Size-1
            J = BN_ALE_KS_Red(I)
            If((J/=0))Then
                Do II = I+1,Reduction_Size
                    JJ = BN_ALE_KS_Red(II)
                    If(JJ/=J)Then
                        Continue
                    Else
                        BN_ALE_KS_Red(II) = 0
                    End IF
                End Do
            End IF
        End Do

        Num = 0
        DO I =1,Reduction_Size
            J = BN_ALE_KS_Red(I)
            If(J/=0)Then
                Num = Num + 1
            End IF
        End Do

        DeAllocate(TempIndex)
        Allocate(TempIndex(Num))
        Num = 0
        DO I =1,Reduction_Size
            J = BN_ALE_KS_Red(I)
            If(J/=0)Then
                Num = Num + 1
                TempIndex(Num) = J
            End IF
        End Do

        DeAllocate(BN_ALE_KS_Red)
        Allocate(BN_ALE_KS_Red(Num))
        BN_ALE_KS_Red = TempIndex

        Call CPU_TIME(sss)
        Call SelectMorePoint
        Call CPU_TIME(eee)
        aaad = aaad + eee-sss

        Err_Red(:,1) =DSQRT( Err_X_Red(:,1)*Err_X_Red(:,1) + Err_Y_Red(:,1) * Err_Y_Red(:,1) )
        IFStop_Sel = Sum(Err_Red(:,1))/NBN_ALE_KS

        Do I =1,MorePoint
            MaxIndex = MaxLoc(Err_Red,1)
            Err_Red(MaxIndex(1),1) = 0.0D0
            AddedIndex(I) = Err_Red(MaxIndex(1),2)
        End Do

        !Write(*,*) 'IFStop_Sel = ',IFStop_Sel,Num,Reduction_Size
    End Do

    Write(*,*) 'MaxErr = ',Err_Red(MaxIndex(1),1),Num,Reduction_Size,aaad

    Call CPU_TIME(sss)
    Call UpdataUseLastSelect
    Call CPU_TIME(eee)
    PRINT*,'Update', eee-sss

    End SubRoutine RBF_Greedy



    SubRoutine SelectMorePoint
    use Lapack95
    use blas95
    Use FlexibleOnly,Only:NBN_ALE_KS,BN_ALE_KS
    Use Variable
    Use ANCF
    Use Local

    Integer::I,J,K,II,JJ,Red_Size
    Real(8)::X0,Y0,RBFDelta_X_Red,RBFDelta_Y_Red

    Real(8),Allocatable::RBF_POS_X_Red(:),RBF_POS_Y_Red(:),RBF_PhiMA_Red(:,:),RBF_PhiMA_Red_Work(:,:)
    Real(8),Allocatable::OLRBF_PhiMA(:,:)
    Real(8),Allocatable::RBF_POS_X_Delta_Red(:),RBF_POS_Y_Delta_Red(:)
    Real(8),Allocatable::RBF_W_Y_Red(:,:)
    Real(8),Allocatable::RBF_W_X_Red(:,:)
    Real(8),Allocatable::RBF_PhiVe(:)
    Real(8),Allocatable::RBF_PhiMatrix_Red(:,:),Temp(:,:)

    Red_Size = Size(BN_ALE_KS_Red)
    Allocate(RBF_POS_X_Red(Red_Size))
    Allocate(RBF_POS_Y_Red(Red_Size))
    Allocate(RBF_W_X_Red(Red_Size,1),RBF_W_Y_Red(Red_Size,1))
    Allocate(RBF_PhiMA_Red(Red_Size,Red_Size),RBF_PhiMA_Red_Work(Red_Size,Red_Size))
    Allocate(OLRBF_PhiMA(Red_Size,Red_Size))
    Allocate(RBF_POS_X_Delta_Red(Red_Size))
    Allocate(RBF_POS_Y_Delta_Red(Red_Size))
    Allocate(RBF_PhiVe(Red_Size))
    Allocate(RBF_PhiMatrix_Red(NBN_ALE_KS,Red_Size),Temp(NBN_ALE_KS,1))

    RBF_POS_X_Red = 0.0D0
    RBF_POS_Y_Red = 0.0D0
    RBF_PhiMA_Red = 0.0D0
    RBF_PhiMA_Red_Work = 0.0D0
    OLRBF_PhiMA   = 0.0D0
    RBF_POS_X_Delta_Red = 0.0D0
    RBF_POS_Y_Delta_Red = 0.0D0
    RBF_W_Y_Red = 0.0D0
    RBF_W_X_Red = 0.0D0
    RBF_PhiMatrix_Red = 0.0D0
    Temp = 0.0D0

    Do I = 1,Red_Size
        J = BN_ALE_KS_Red(I)
        DO II = 1,Red_Size
            JJ =  BN_ALE_KS_Red(II)
            RBF_PhiMA_Red(I,II) = RBF_PhiMA(J,JJ)
        End Do
        RBF_POS_X_Delta_Red(I) = RBF_POS_X_Delta_Into_Red(J)
        RBF_POS_Y_Delta_Red(I) = RBF_POS_Y_Delta_Into_Red(J)
    End Do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MKL
    OLRBF_PhiMA = RBF_PhiMA_Red
    call SYSV(OLRBF_PhiMA,RBF_POS_Y_Delta_Red,info = ii)
    IF (II/=0)THEN
        PRINT*,'OPS SYSV CAN NOT SLOVE THE PROBLEM1 Red'
        PAUSE
        STOP
    END IF

    OLRBF_PhiMA = RBF_PhiMA_Red
    call SYSV(OLRBF_PhiMA,RBF_POS_X_Delta_Red,info = ii)
    IF (II/=0)THEN
        PRINT*,'OPS SYSV CAN NOT SLOVE THE PROBLEM2 Red'
        PAUSE
        STOP
    END IF

    RBF_W_Y_Red(:,1) = RBF_POS_Y_Delta_Red(:)
    RBF_W_X_Red(:,1) = RBF_POS_X_Delta_Red(:)


    Do I = 1, NBN_ALE_KS
        J= BN_ALE_KS(I)

        RBF_PhiVe = 0.0
        RBFDelta_X_Red = 0.0
        RBFDelta_Y_Red = 0.0
        Do II = 1,Red_Size
            JJ = BN_ALE_KS_Red(II)
            RBF_PhiMatrix_Red(I,II) = RBF_PhiMA(I,JJ)
        End Do
    End Do

    RBF_PhiMA_Red_Work = RBF_PhiMatrix_Red
    Call gemm(RBF_PhiMA_Red_Work,RBF_W_X_Red,Temp)
    Err_X_Red(:,1)= Temp(:,1) - RBF_POS_X_Delta_Into_Red(:)

    Temp = 0.0D0
    RBF_PhiMA_Red_Work = RBF_PhiMatrix_Red
    Call gemm(RBF_PhiMA_Red_Work,RBF_W_Y_Red,Temp)
    Err_Y_Red(:,1) = Temp(:,1) - RBF_POS_Y_Delta_Into_Red(:)


    Err_X_Red(:,2) = (/(i,i=1,NBN_ALE_KS)/)
    Err_Y_Red(:,2) = (/(i,i=1,NBN_ALE_KS)/)


    End SubRoutine


    SubRoutine UpdataUseLastSelect
    Use FlexibleOnly,Only:NBN_ALE_KS,BN_ALE_KS
    Use Variable
    Use ANCF
    use Lapack95
    Use Blas95
    Use Local

    Integer::I,J,K,II,JJ,Red_Size
    Real(8)::X0,Y0,RBFDelta_X_Red,RBFDelta_Y_Red

    Real(8),Allocatable::RBF_POS_X_Red(:),RBF_POS_Y_Red(:),RBF_PhiMA_Red(:,:)
    Real(8),Allocatable::OLRBF_PhiMA(:,:)
    Real(8),Allocatable::RBF_POS_X_Delta_Red(:),RBF_POS_Y_Delta_Red(:)
    Real(8),Allocatable::RBF_W_Y_Red(:)
    Real(8),Allocatable::RBF_W_X_Red(:)
    Real(8),Allocatable::RBF_PhiVe(:)

    Red_Size = Size(BN_ALE_KS_Red)
    Allocate(RBF_POS_X_Red(Red_Size))
    Allocate(RBF_POS_Y_Red(Red_Size))
    Allocate(RBF_PhiMA_Red(Red_Size,Red_Size))
    Allocate(OLRBF_PhiMA(Red_Size,Red_Size))
    Allocate(RBF_POS_X_Delta_Red(Red_Size))
    Allocate(RBF_POS_Y_Delta_Red(Red_Size))
    Allocate(RBF_PhiVe(Red_Size))

    RBF_POS_X_Red = 0.0D0
    RBF_POS_Y_Red = 0.0D0
    RBF_PhiMA_Red = 0.0D0
    OLRBF_PhiMA   = 0.0D0
    RBF_POS_X_Delta_Red = 0.0D0
    RBF_POS_Y_Delta_Red = 0.0D0
    RBF_W_Y_Red = 0.0D0
    RBF_W_X_Red = 0.0D0


    Do I = 1,Red_Size
        J = BN_ALE_KS(BN_ALE_KS_Red(I))
        RBF_POS_X_Red(I) = X_ALE(J)
        RBF_POS_Y_Red(I) = Y_ALE(J)

        DO II = 1,Red_Size
            JJ =  BN_ALE_KS(BN_ALE_KS_Red(II))
            X0 = X_ALE(JJ)
            y0 = Y_ALE(JJ)
            CALL phifun(RBF_POS_X_Red(I),RBF_POS_Y_Red(I),X0,y0,RBF_PhiMA_Red(I,II))
        End Do

        RBF_POS_X_Delta_Red(I) = X1_ALE(J)-X_ALE(J)
        RBF_POS_Y_Delta_Red(I) = Y1_ALE(J)-Y_ALE(J)
    End Do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MKL
    OLRBF_PhiMA = RBF_PhiMA_Red
    call SYSV(OLRBF_PhiMA,RBF_POS_Y_Delta_Red,info = ii)
    IF (II/=0)THEN
        PRINT*,'OPS SYSV CAN NOT SLOVE THE PROBLEM1 Red'
        PAUSE
        STOP
    END IF

    OLRBF_PhiMA = RBF_PhiMA_Red
    call SYSV(OLRBF_PhiMA,RBF_POS_X_Delta_Red,info = ii)
    IF (II/=0)THEN
        PRINT*,'OPS SYSV CAN NOT SLOVE THE PROBLEM2 Red'
        PAUSE
        STOP
    END IF

    RBF_W_Y_Red = RBF_POS_Y_Delta_Red
    RBF_W_X_Red = RBF_POS_X_Delta_Red


    Do I  = 1, NN_ALe
        X0 = X_ALE(I)
        Y0 = Y_ALE(I)
        If((X0>=-1.0D0).And.(Dabs(Y0)<=3.5D0))Then
            RBF_PhiVe = 0.0
            RBFDelta_X_Red = 0.0
            RBFDelta_Y_Red = 0.0
            Do II = 1,Red_Size
                CALL phifun(X0,Y0,RBF_POS_X_Red(II),RBF_POS_Y_Red(II),RBF_PhiVe(II))
                RBFDelta_X_Red = RBFDelta_X_Red + RBF_PhiVe(II) * RBF_W_X_Red(II)
                RBFDelta_Y_Red = RBFDelta_Y_Red + RBF_PhiVe(II) * RBF_W_Y_Red(II)
            End Do
            X1_ALE(I) = X_ALE(I) + RBFDelta_X_Red
            Y1_ALE(I) = Y_ALE(I) + RBFDelta_Y_Red
        Else
            X1_ALE(I) = X_ALE(I)
            Y1_ALE(I) = Y_ALE(I)
        End If
    End Do

    End SubRoutine UpdataUseLastSelect

    SubRoutine RBF_INIT
    Use FlexibleOnly,Only:NBN_ALE_KS,BN_ALE_KS
    Use Variable
    Use ANCF
    use Lapack95
    Use Blas95
    Use Local

    Real(8)::X0,Y0,XX0,YY0
    Integer::I,J,II,Num
    Integer,ALLOCATABLE::IPIV(:)

    Num = 0
    Allocate(Mesh_Update_Index(NN_ALE))
    Mesh_Update_Index = 0

    Do I = 1, NN_ALe
        X0 = X_ALE(I)
        Y0 = Y_ALE(I)
        If(X0>-5.0D0.And.DABS(Y0)<6.0D0)Then
            Num = Num + 1
            Mesh_Update_Index(Num) = I
        End If
    End Do

    Mesh_Update_Num = Num

    Allocate(RBF_PhiMA_Constant(Mesh_Update_Num,NBN_ALE_KS),Left_Constant(Mesh_Update_Num,NBN_ALE_KS))
    RBF_PhiMA_Constant = 0.0D0
    Left_Constant = 0.0D0

    Do I = 1, Mesh_Update_Num
        X0 = X_ALE(Mesh_Update_Index(I))
        Y0 = Y_ALE(Mesh_Update_Index(I))
        Do J = 1,NBN_ALE_KS
            II = BN_ALE_KS(J)
            XX0 = X_ALE(II)
            YY0 = Y_ALE(II)
            CALL phifun(X0,Y0,XX0,YY0,RBF_PhiMA_Constant(I,J))
        End Do
    End Do

    Allocate(RBF_PhiMA(NBN_ALE_KS,NBN_ALE_KS),IPIV(NBN_ALE_KS))
    RBF_PhiMA = 0.0D0

    Do I = 1,NBN_ALE_KS
        J = BN_ALE_KS(I)
        X0 = X_ALE(J)
        Y0 = Y_ALE(J)
        DO II = 1,NBN_ALE_KS
            JJ =  BN_ALE_KS(II)
            XX0 = X_ALE(JJ)
            YY0 = Y_ALE(JJ)
            CALL phifun(X0,Y0,XX0,YY0,RBF_PhiMA(I,II))
        End Do
    End Do

    call getrf(RBF_PhiMA,ipiv=ipiv)
    Call getri(RBF_PhiMA,ipiv)
    call gemm(RBF_PhiMA_Constant,RBF_PhiMA,Left_Constant)

    End SubRoutine RBF_INIT



