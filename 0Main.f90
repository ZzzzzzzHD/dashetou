

    PROGRAM CIRCULARCYLINDER_WITH_RIGIDPLATE
    USE CONSTANT
    USE VARIABLE
    USE FlexibleOnly
    USE Local

    IMPLICIT NONE
    INTEGER:: I,J,K,ISURFACEOUT,COUNT,int_CPU_START,int_CPU_FINISH
    CHARACTER(12):: FN,FORCE_FILENAME,ANGLE_FILENAME,Y_FILENAME

    REAL(8):: PI
    REAL(8)::CPU_START,CPU_FINISH,CPU_S,CPU_MIN,CPU_H!JY!记录计算时间用
    ALLOCATE(THETA(2))
    THETA(1)=0.8D0;       THETA(2)=1.0D0

    Call RANDOM_SEED()
    CALL CPU_TIME(CPU_START)!JY!记录计算时间用
    CALL README

    Call Flexible_MeshREAD

    CALL SURROUND_ELE_ORIG	                        !找单元周围的单元
    CALL POLYBN_0226(1)                             !找到边界单元节点
    CALL NODE_ADJ_NODE_LMM                          !找到节点周围的节点

    CALL EXTEND_MESH_RIGHT_0226                     !延伸右网格
    CALL EXTEND_MESH_LEFT_0226                      !延伸左网格
    CALL SURROUND_ELE_ORIG
    CALL POLYBN_0226(2)
    CALL EXTEND_MESH_UP_0226
    CALL EXTEND_MESH_DOWN_0226

    CALL SURROUND_ELE_ORIG



    CALL ELEMENTTYPE                              ! 单元类型（四边形）
    CALL CENTRE_ELE                               ! 求单元中心坐标

    CALL Clock2AntiClock_Total
    CALL POLYBN_0226(3)

    CALL SURROUND_ELE                             ! 找单元周围的单元
    CALL NODE_ADJ_NODE                            ! 找到节点周围的节点
    !Call For_FindMinELEArea                      ! 仅柔性用到

    CALL C_ISOPE                                  ! 形函数及求导

    IF (ALLOCATED(X))		DEALLOCATE(X)
    IF (ALLOCATED(Y))		DEALLOCATE(Y)
    IF (ALLOCATED(NODEELE))         DEALLOCATE(NODEELE)

    NN=NN_ORIG
    NE=NE_ORIG

    ALLOCATE(X(NN),Y(NN),NODEELE(NE,4))
    X=X_ORIG
    Y=Y_ORIG
    NODEELE=NODEELE_ORIG

    DEALLOCATE(X_ORIG,Y_ORIG,NODEELE_ORIG)

    IF (ALLOCATED(LBN)) DEALLOCATE(LBN)
    IF (ALLOCATED(RBN)) DEALLOCATE(RBN)
    IF (ALLOCATED(UBN)) DEALLOCATE(UBN)
    IF (ALLOCATED(DBN)) DEALLOCATE(DBN)
    IF (ALLOCATED(PBN)) DEALLOCATE(PBN)
    IF (ALLOCATED(SBN)) DEALLOCATE(SBN)

    NLBN=NLBN_ORIG
    NRBN=NRBN_ORIG
    NUBN=NUBN_ORIG
    NDBN=NDBN_ORIG
    NPBN=NPBN_ORIG
    NSBN=NSBN_ORIG
    ALLOCATE(LBN(NLBN),RBN(NRBN),UBN(NUBN),DBN(NDBN),PBN(NPBN),SBN(NSBN))

    LBN=LBN_ORIG
    RBN=RBN_ORIG
    UBN=UBN_ORIG
    DBN=DBN_ORIG
    PBN=PBN_ORIG
    SBN=SBN_ORIG
    DEALLOCATE(LBN_ORIG,RBN_ORIG,UBN_ORIG,DBN_ORIG,PBN_ORIG,SBN_ORIG)

    IF(ALLOCATED(RBE)) DEALLOCATE(RBE)
    IF(ALLOCATED(PBE)) DEALLOCATE(PBE)
    IF(ALLOCATED(SBE)) DEALLOCATE(SBE)

    NRBE=NRBE_ORIG
    NPBE=NPBE_ORIG
    NSBE=NSBE_ORIG
    ALLOCATE(RBE(NRBE,2),PBE(NPBE,2),SBE(NSBE,2))

    RBE=RBE_ORIG
    PBE=PBE_ORIG
    SBE=SBE_ORIG
    DEALLOCATE(RBE_ORIG,PBE_ORIG,SBE_ORIG)


    ALLOCATE(XBC_ALE(NBN_ALE),YBC_ALE(NBN_ALE))
    ALLOCATE(X0_ALE(NN_ALE),Y0_ALE(NN_ALE),X1_ALE(NN_ALE),Y1_ALE(NN_ALE))
    X0_ALE=X_ALE
    Y0_ALE=Y_ALE
    !*********************************************************************!
    !*********************************************************************!
    T=0.0D0
    DT=DT000

    ALLOCATE(U0(NN),V0(NN),PN(NN),PN1(NN),VOR(NN))
    U0=0.0D0
    V0=0.0D0
    PN=0.0D0

    ALLOCATE(UM_ALE(NN),VM_ALE(NN))
    UM_ALE=0.0D0
    VM_ALE=0.0D0
    ALLOCATE(DELTX(NN_ALE),DELTY(NN_ALE))
    DELTX=0.0
    DELTY=0.0
    !*********************************************************************!
    !*********************************************************************!

    CALL INDEX00
    CALL INDEXALE

    ALLOCATE(U  (NN));	ALLOCATE(V  (NN));	ALLOCATE(UB (NN));ALLOCATE(UB2 (NN)); ALLOCATE(VB2 (NN))
    ALLOCATE(VB (NN));	ALLOCATE(CLE(NE));	ALLOCATE(AM (NN))
    ALLOCATE(AP (N1DA), BP(NN))

    CALL MIN_X0_X00

    CAll Clock2AntiClock_ALE

    IF (MESH_TRANSFORM==1) THEN
        ALLOCATE(AP_ALE(N1DA_ALE),BP_ALE(NN_ALE),APX_ALE(N1DA_ALE),APY_ALE(N1DA_ALE))
    END IF

    CALL FirstMeshLayer_Adj_Structure

    IF (MESH_TRANSFORM==1) THEN
        AP_ALE=0.0D0
        CALL COEFFICIENTMATRIXALE	! ASSEMBLED AT X_ALE AND Y_ALE
    ElseIf(MESH_TRANSFORM==3) Then
        Call CalPlatePointIndex
        Call GatherSp
        Call RBF_INIT
    END IF

    If((IROTA==510).OR.(IROTA==5102).OR.(IROTA==5101).OR.(IROTA==520))Then
        ALLOCATE(Sec_Rota(NSBN_ALE,3),Sec_Rota0(NSBN_ALE,3))
        Sec_Rota = 0.0D0
        Sec_Rota0 = Sec_Rota
    End if

    !=========================================================================!
    DIS0_ANG=0.0D0; VEL0_ANG=0.0D0;  ACC0_ANG=0.0D0
    DIS_ANG =0.0D0;  VEL_ANG=0.0D0;  ACC_ANG =0.0D0
    tor0=0.0D0; tor=0.0d0

    !========JY========
    DIS0_X=0.0D0; VEL0_X=0.0D0;  ACC0_X=0.0D0
    DIS_X =0.0D0;  VEL_X=0.0D0;  ACC_X =0.0D0

    DIS0_Y=0.0D0; VEL0_Y=0.0D0;  ACC0_Y=0.0D0
    DIS_Y =0.0D0;  VEL_Y=0.0D0;  ACC_Y =0.0D0
    
    DIS0_Y_UP=0.0D0; VEL0_Y_UP=0.0D0;  ACC0_Y_UP=0.0D0
    DIS_Y_UP =0.0D0;  VEL_Y_UP=0.0D0;  ACC_Y_UP =0.0D0

    CL0=0.0D0;CL=0.0D0;CL1=0.0D0;CL2=0.0D0
    CD0=0.0D0;CD=0.0D0;CD1=0.0D0;CD2=0.0D0

    WORK=0.0D0;WORK_TOR=0.0D0;WORK_F_Y=0.0D0;WORK_F_X=0.0D0

    FORCE_X_CENTER_0=0.0D0

    IF(CSD_METHOD == 1)THEN
        !====================================Flexible=====================================!
        Allocate(Fdis_ang(FL_Params_ne+1),Fdis_Y(FL_Params_ne+1),Fdis0_ang(FL_Params_ne+1),Fdis0_Y(FL_Params_ne+1))
        Allocate(FVel_ang(FL_Params_ne+1),FVel_Y(FL_Params_ne+1),FVel0_ang(FL_Params_ne+1),FVel0_Y(FL_Params_ne+1))
        Allocate(FAcc_ang(FL_Params_ne+1),FAcc_Y(FL_Params_ne+1),FAcc0_ang(FL_Params_ne+1),FAcc0_Y(FL_Params_ne+1))

        Fdis_ang    =0.0d0
        Fdis_Y      =0.0d0
        Fdis0_ang   =0.0d0
        Fdis0_Y     =0.0d0
        FVel_ang    =0.0d0
        FVel_Y      =0.0d0
        FVel0_ang   =0.0d0
        FVel0_Y     =0.0d0
        FAcc_ang    =0.0d0
        FAcc_Y      =0.0d0
        FAcc0_ang   =0.0d0
        FAcc0_Y     =0.0d0

        FL_L   =  eta_l
        FL_Params_rho =   MASSR * DENSITY
        FL_Params_n  =    2*(FL_Params_ne+1)
        FL_Params_L  =    FL_L / FL_Params_ne
        Fixed = [3,0]

        Allocate(M_MASSM(FL_Params_n,FL_Params_n))
        M_MASSM = 0.0d0
        Call Make_Mass_Matrix

        Allocate(K_StiffM(FL_Params_n,FL_Params_n))
        K_StiffM = 0.0d0

        Call Make_Stiff_Matrix

        Allocate(FEULERBER(FL_Params_n))
        FEULERBER    = 0.0d0

        Allocate(Keff(FL_Params_n,FL_Params_n))
        Keff= 0.0d0

        Allocate(Feff(FL_Params_n))
        Feff= 0.0d0

        Allocate(Fq(FL_Params_n))
        Allocate(Fqd(FL_Params_n))
        Allocate(Fqdd(FL_Params_n))
        Fq   = 0.0d0
        Fqd  = 0.0d0
        Fqdd = 0.0d0

        Allocate(Foq(FL_Params_n))
        Allocate(Foqd(FL_Params_n))
        Allocate(Foqdd(FL_Params_n))
        Foq     = 0.0d0
        Foqd    = 0.0d0
        Foqdd   = 0.0d0

    ELSEIF(CSD_METHOD == 2)THEN
        CALL ANCF_INIT
    ELSE
        PRINT*,'INVALID CSD METHOD--ERROR CODE#002C'
    END IF

    !====================================Flexible=====================================!

    !IF(IROTA==201)THEN
    !    WRITE(FORCE_FILENAME,'(a,i2.2,a)') 'FOECE_201_',FRE_Y, '.DAT'
    !    WRITE(ANGLE_FILENAME,'(a,i2.2,a)') 'ANGLE_201_'FRE_Y, '.DAT'
    !    WRITE(Y_FILENAME,'(a,i2.2,a)') 'Y_201_',FRE_Y, '.DAT'
    !ELSE IF(IROTA==202)THEN
    !    WRITE(FORCE_FILENAME,'(a,i2.2,ai2.2,a,)') 'FOECE_202_U_ry',FRE_Y,'U_rr',FRE_R, '.DAT'
    !    WRITE(ANGLE_FILENAME,'(a,i2.2,ai2.2,a,)') 'ANGLE_202_U_ry',FRE_Y,'U_rr',FRE_R, '.DAT'
    !    WRITE(Y_FILENAME,'(a,i2.2,a,i2.2,a,)') 'Y_202_U_ry',FRE_Y,'U_rr',FRE_R, '.DAT'
    !ELSE IF (IROTA==2012)THEN
    !    WRITE(FORCE_FILENAME,'(a,i2.2,a)') 'FOECE_2012_',FRE_R, '.DAT'
    !    WRITE(ANGLE_FILENAME,'(a,i2.2,a)') 'ANGLE_2012_'FRE_R, '.DAT'
    !    WRITE(Y_FILENAME,'(a,i2.2,a)') 'Y_2012_',FRE_R, '.DAT'
    !END IF

    IF(RE_START==0)THEN
        OPEN(UNIT=999,FILE='FORCE.DAT')
        WRITE(999,'(A99)') ' VARIABLES = "T" "CD" "CL" "CD1" "CL1" "CD2" "CL2"  '
        OPEN(UNIT=1025,FILE='WORK.DAT')
        OPEN(UNIT=1026,FILE='WORK_TOR.DAT')
        OPEN(UNIT=1027,FILE='WORK_F_Y.DAT')
        OPEN(UNIT=1028,FILE='WORK_F_X.DAT')
        WRITE(1025,'(A99)') ' VARIABLES = "T+0.5DT" "TOTAL_WORK" "DW"'
        WRITE(1026,'(A99)') ' VARIABLES = "T+0.5DT" "TOR_TOTAL_WORK" "TOR_DW"'
        WRITE(1027,'(A99)') ' VARIABLES = "T+0.5DT" "CL_TOTAL_WORK" "CL_DW"'
        WRITE(1028,'(A99)') ' VARIABLES = "T+0.5DT" "CD_TOTAL_WORK" "CD_DW"'
        !IF(IROTA==1)  THEN
        !OPEN(UNIT=777,FILE='TORQUE.DAT')  !zhangcheng
        !OPEN(UNIT=888,FILE='ANGLE.DAT')
        !WRITE(777,'(A99)') ' VARIABLES = "T" "TOR" "TOR1" "TORP" "TORS" "TOR_PIPE_TAU"  '
        !WRITE(888,'(A99)') ' VARIABLES = "T" "DIS_ANG" "VEL_ANG" "ACC_ANG"  '
        !END IF
        IF(IROTA==300)THEN
            OPEN(UNIT=666,FILE='Y.DAT')
            OPEN(UNIT=444,FILE='X.DAT')
            OPEN(UNIT=233,FILE='X_Y.DAT')
            WRITE(666,'(A55)') ' VARIABLES = "T" "DIS_Y" "VEL_Y" "ACC_Y" '
            WRITE(444,'(A55)') ' VARIABLES = "T" "DIS_X" "VEL_X" "ACC_X" '
            WRITE(233,'(A55)') ' VARIABLES = "T" "X" "Y" '
        ELSE IF(IROTA==200.or.IROTA==201.or.IROTA==202.or.IROTA==2012) THEN
            OPEN(UNIT=666,FILE='Y.DAT')
            WRITE(666,'(A55)') ' VARIABLES = "T" "DIS_Y" "VEL_Y" "ACC_Y" '
        ELSE IF(IROTA==2201) THEN
            OPEN(UNIT=6661,FILE='Y_UP.DAT')
            OPEN(UNIT=666,FILE='Y_DO.DAT')
            WRITE(6661,'(A55)') ' VARIABLES = "T" "DIS_Y_UP" "VEL_Y_UP" "ACC_Y_UP" '
            WRITE(666,'(A55)') ' VARIABLES = "T" "DIS_Y" "VEL_Y" "ACC_Y" '         
        ELSE IF(IROTA==2101) THEN
            OPEN(UNIT=666,FILE='Y_DO.DAT')
            WRITE(666,'(A55)') ' VARIABLES = "T" "DIS_Y" "VEL_Y" "ACC_Y" '         
        END IF
        OPEN(UNIT=2331,FILE='FORCE_CENTER.DAT')
        WRITE(2331,'(A55)') ' VARIABLES = "T" "X" "Y" '
        OPEN(UNIT=2332,FILE='FORCE_CENTER_RELA.DAT')
        WRITE(2332,'(A55)') ' VARIABLES = "T" "X" "Y" '

    ELSE IF(RE_START==1) THEN
        CALL RESTART
        OPEN(UNIT=999,FILE='FORCE.DAT',POSITION='APPEND')
        !OPEN(UNIT=777,FILE='TORQUE.DAT',POSITION='APPEND')
        !OPEN(UNIT=888,FILE='ANGLE.DAT',POSITION='APPEND')
        OPEN(UNIT=1025,FILE='WORK.DAT',POSITION='APPEND')
        OPEN(UNIT=1026,FILE='WORK_TOR.DAT',POSITION='APPEND')
        OPEN(UNIT=1027,FILE='WORK_F_Y.DAT',POSITION='APPEND')
        OPEN(UNIT=1028,FILE='WORK_F_X.DAT',POSITION='APPEND')
        IF(IROTA==300)THEN
            OPEN(UNIT=666,FILE='Y.DAT',POSITION='APPEND')
            OPEN(UNIT=444,FILE='X.DAT',POSITION='APPEND')
            OPEN(UNIT=233,FILE='X_Y.DAT',POSITION='APPEND')
        ELSE IF(IROTA==200.or.IROTA==201.or.IROTA==202.or.IROTA==2012) THEN
            OPEN(UNIT=666,FILE='Y.DAT',POSITION='APPEND')
        ELSE IF(IROTA==2201) THEN
            OPEN(UNIT=6661,FILE='Y_UP.DAT',POSITION='APPEND')
            OPEN(UNIT=666,FILE='Y_DO.DAT',POSITION='APPEND')   
        ELSE IF(IROTA==2101) THEN
            OPEN(UNIT=666,FILE='Y_DO.DAT',POSITION='APPEND')
        END IF
        OPEN(UNIT=2331,FILE='FORCE_CENTER.DAT',POSITION='APPEND')
        OPEN(UNIT=2332,FILE='FORCE_CENTER_RELA.DAT',POSITION='APPEND')
        T_00=T
    END IF
    !==========================================================================


    NSTEP=NSTEP+1
    PI=DACOS(-1.0D0)
    COUNT=1

    DO WHILE (NSTEP <= MAXSTEP)
        TT0=T
        T=T+DT

        !jy!test==============
        IF(DT<0.00015)THEN
            IF (COUNT.LE.5)THEN
                CALL  OUTPUT
                COUNT=COUNT+1
            ELSE
                PAUSE
                STOP
            END IF
        END IF

        IF(DT<0.000001.and.nstep<100)THEN
            dt=0.001
        ELSE IF(DT<0.000001.and.nstep.GT.100) THEN
            Call OutPut
            WRITE(*,*) 'DT SMALLER THEN 1E-6 !'
            PAUSE
            STOP
        END IF
        !====================

        IF(T>time_stop+0.0009)THEN
            !CALL  OUTPUT_ale_right
            CALL CPU_TIME(CPU_FINISH)
            CPU_S=CPU_FINISH-CPU_START
            CPU_MIN=CPU_S/60.0
            CPU_H=CPU_MIN/60.0

            !JY!输出计算时间
            OPEN(UNIT=1357,FILE='TMIE.TXT')
            WRITE(1357,*)'USING TIME:  ',CPU_S,'S    ',CPU_MIN,'MIN    ',CPU_H,'H    '
            CLOSE(1357)
            WRITE(*,*)'USING TIME:  ',CPU_S,'S    ',CPU_MIN,'MIN    ',CPU_H,'H    '

            IF(DABS(T*DT000-NSTEP)>1E-6) THEN!JY!如果动态时间步长起了作用，那么在结束运行时进行输出，若需要续算，直接调用最后的b文件就行，不需要再删除多余数据
                CALL OUTPUT
                CALL BINRESULT
            END IF

            STOP
        END IF

        WRITE(*,'(A10,I12,10X,A10,F12.6)')'NUM STEP =',NSTEP,'TIME     =',T
        WRITE(*,'(A10,I12,10X,A10,I12  )')'NUM NODE =',NN   ,'NUM ELEM =',NE

        CALL CHARACTERISTIC_LENGTH

        AM =0.0D0
        AP =0.0D0
        CALL COEFFICIENTMATRIX                            !形成单元有限元方程和压力泊松方程的左端项,组装压力泊松方程左端项，N-S方程左端项

        U =   U0;		V =   V0
        UB=0.0D0;		VB=0.0D0
        CALL STEP1
        !CALL STEP1_SOLUTION

        !IF (NSTEP>10)THEN
        !    CALL OutFlow_Boundary!Newman Type Velovity Boundary Condition->Open Boundary Condition(OBC)
        !END IF

        Do i=1,nn
            Ub(I)=Ub(I)*(1.0D0/AM(I))
            Vb(I)=Vb(I)*(1.0D0/AM(I))
        End Do

        Do i=1,nn
            U(I)=U0(I)+Ub(I)
            V(I)=V0(I)+Vb(I)
        End Do

        CALL IMPOSE_VELOCITY_BC

        BP=0.0D0
        CALL STEP2                                !组装压力泊松方程右端项
        CALL BOUNDARY_INTEGRAL_PREEQU             !压力泊松方程边界积分项
        CALL PRESSURE_1BC                         !压力右边界条件
        CALL E_ET_Q                               !与 共轭梯度解方程相关
        PN1=PN

        CALL SOLUTION_BICGS(AP,BP,PN1)	          !解泊松方程

        UB2=0.0D0
        VB2=0.0D0
        CALL STEP3
        CALL STEP3_SOLUTION

        PN=PN1
        DO I=1,NN
            U(I)=U(I)+UB2(I)
            V(I)=V(I)+VB2(I)
        END DO

        CALL IMPOSE_VELOCITY_BC

        U0=U                                            !将每一步计算出来的U V作为下一步计算的初值，也就是u^n (^表示上标)
        V0=V

        IF (MOD(NSTEP,NW)==0.OR.NSTEP<5.OR.(NSTEP>=NW_START.and.NSTEP<=NW_END.and.MOD(NSTEP,nw)==0)) THEN
            CALL  OUTPUT
        END IF

        IF(T>TIME_CALL_FORCE) 	CALL FORCE


        IF(IROTA==200)THEN
            IF (T>T_START)  THEN
                !==================================dis_ang=======================================================================
                DXRK=(2.0D0*PI*FRE_R)**2+8.0D0*PI*FRE_R*DAMPING_R/DT+4.0D0/DT/DT    !等效刚度

                DXRF=50.0D0*TOR/(25.0D0*PI+2.0D0)/MASSR+&
                    &-(48.0D0/(75.0D0*PI+52.0D0))*(-ACC0_X*DSIN(DIS0_ANG)+ACC0_Y*DCOS(DIS0_ANG))&
                    &+4.0D0*PI*FRE_R*DAMPING_R*(2.0D0*DIS0_ANG/DT+VEL0_ANG)+4.0D0*DIS0_ANG/DT/DT+4.0D0*VEL0_ANG/DT+ACC0_ANG    !等效荷载

                DIS_ANG=DXRF/DXRK
                VEL_ANG=2.0D0/DT*(DIS_ANG-DIS0_ANG)-VEL0_ANG
                ACC_ANG=4.0D0/DT/DT*(DIS_ANG-DIS0_ANG)-4.0D0/DT*VEL0_ANG-ACC0_ANG
                !==================================dis_ang=======================================================================

                !==================================dis_y==========================================================================
                DXYK=(2.0D0*PI*FRE_Y)**2+8.0D0*PI*FRE_Y*DAMPING_Y/DT+4.0D0/DT/DT    !等效刚度
                DXYF=50.0D0*CL/(25.0D0*PI+2.0D0)/MASSR-(2.0D0/(25.0D0*PI+2))*(ACC0_ANG*DCOS(DIS0_ANG)-VEL0_ANG**2*DSIN(DIS0_ANG))+4.0D0*PI*FRE_Y*DAMPING_Y*(2.0D0*DIS0_Y/DT+VEL0_Y)+4.0D0*DIS0_Y/DT/DT+4.0D0*VEL0_Y/DT+ACC0_Y    !等效荷载
                DIS_Y=DXYF/DXYK
                VEL_Y=2.0D0/DT*(DIS_Y-DIS0_Y)-VEL0_Y
                ACC_Y=4.0D0/DT/DT*(DIS_Y-DIS0_Y)-4.0D0/DT*VEL0_Y-ACC0_Y
                !==================================dis_y==========================================================================

                IF(MOD(NSTEP,5)==0) THEN
                    !WRITE(888,'(4E18.10)')T,DIS_ANG,VEL_ANG,ACC_ANG
                    WRITE(666,'(4E18.10)')T,DIS_Y,VEL_Y,ACC_Y
                END IF

            END IF
        END IF

        IF(IROTA==201)THEN
            IF (T>T_START)  THEN
                !==================================dis_y==========================================================================
                DXYK=(2.0D0*PI*FRE_Y)**2+8.0D0*PI*FRE_Y*DAMPING_Y/DT+4.0D0/DT/DT    !等效刚度
                DXYF=CL/2/(1.0D0-(4.0D0-PI)*ETA_R**2.0D0+1.0D0*ETA_T*ETA_L)/MASSR+4.0D0*PI*FRE_Y*DAMPING_Y*(2.0D0*DIS0_Y/DT+VEL0_Y)+4.0D0*DIS0_Y/DT/DT+4.0D0*VEL0_Y/DT+ACC0_Y    !等效荷载
                DIS_Y=DXYF/DXYK
                VEL_Y=2.0D0/DT*(DIS_Y-DIS0_Y)-VEL0_Y
                ACC_Y=4.0D0/DT/DT*(DIS_Y-DIS0_Y)-4.0D0/DT*VEL0_Y-ACC0_Y
                !==================================dis_y==========================================================================
                IF(MOD(NSTEP,5)==0) THEN
                    !WRITE(888,'(4E18.10)')T,DIS_ANG,VEL_ANG,ACC_ANG
                    WRITE(666,'(4E18.10)')T,DIS_Y,VEL_Y,ACC_Y
                END IF
            END IF
        END IF

        IF(IROTA==202)THEN
            IF (T>T_START)  THEN
                !==================================dis_ang=======================================================================
                DXRK=(2.0D0*PI*FRE_R)**2.0D0+8.0D0*PI*FRE_R*DAMPING_R/DT+4.0D0/DT/DT    !等效刚度
                DXRF=TOR/2.0D0/(1.0D0-(4.0D0-PI)*ETA_R**2.0D0+1.0D0*ETA_T*ETA_L)/MASSR &
                    &-ETA_T*ETA_L*(1.0D0+ETA_L)*(ACC0_Y*DCOS(DIS0_ANG)-ACC0_X*DSIN(DIS0_ANG)) &
                    &/(2.0D0*(1.0D0/6.0D0-2.0D0*eta_R**4.0D0/3.0D0-2.0D0*ETA_R**2.0D0*(1.0D0-ETA_R)**2.0D0+PI*ETA_R**4.0D0/2.0D0-32.0D0*ETA_R**4.0D0/9.0D0/PI &
                    &+2.0D0*PI*ETA_R**2.0D0*(1.0D0/2.0D0-ETA_R+4.0D0*ETA_R/3.0D0/PI)**2.0D0+ETA_T*ETA_L**3.0D0/12.0D0 +ETA_T*ETA_L*(1+ETA_L)**2/4))&
                    &+4.0D0*PI*FRE_R*DAMPING_R*(2.0D0*DIS0_ANG/DT+VEL0_ANG)+4.0D0*DIS0_ANG/DT/DT+4.0D0*VEL0_ANG/DT+ACC0_ANG    !等效荷载

                DIS_ANG=DXRF/DXRK
                VEL_ANG=2.0D0/DT*(DIS_ANG-DIS0_ANG)-VEL0_ANG
                ACC_ANG=4.0D0/DT/DT*(DIS_ANG-DIS0_ANG)-4.0D0/DT*VEL0_ANG-ACC0_ANG
                !==================================dis_ang=======================================================================

                !==================================dis_y==========================================================================
                DXYK=(2.0D0*PI*FRE_Y)**2.0D0+8.0D0*PI*FRE_Y*DAMPING_Y/DT+4.0D0/DT/DT    !等效刚度
                DXYF=CL/2/(1.0D0-(4.0D0-PI)*ETA_R**2.0D0+1.0D0*ETA_T*ETA_L)/MASSR &
                    &-(1+ETA_L)*ETA_T*ETA_L/2.0D0/(1.0D0-(4.0D0-PI)*ETA_R**2.0D0+1.0D0*ETA_T*ETA_L)*(ACC0_ANG*DCOS(DIS0_ANG)-VEL0_ANG**2*DSIN(DIS0_ANG))&
                    &+4.0D0*PI*FRE_Y*DAMPING_Y*(2.0D0*DIS0_Y/DT+VEL0_Y)+4.0D0*DIS0_Y/DT/DT+4.0D0*VEL0_Y/DT+ACC0_Y    !等效荷载
                DIS_Y=DXYF/DXYK
                VEL_Y=2.0D0/DT*(DIS_Y-DIS0_Y)-VEL0_Y
                ACC_Y=4.0D0/DT/DT*(DIS_Y-DIS0_Y)-4.0D0/DT*VEL0_Y-ACC0_Y
                !==================================dis_y==========================================================================

                IF(MOD(NSTEP,5)==0) THEN
                    !WRITE(888,'(4E18.10)')T,DIS_ANG,VEL_ANG,ACC_ANG
                    WRITE(666,'(4E18.10)')T,DIS_Y,VEL_Y,ACC_Y
                END IF

            END IF
        END IF

        IF(IROTA==2012)THEN
            IF (T>T_START)  THEN
                !==================================dis_ang=======================================================================
                DXRK=(2.0D0*PI*FRE_R)**2+8.0D0*PI*FRE_R*DAMPING_R/DT+4.0D0/DT/DT    !等效刚度
                DXRF=TOR/2.0D0/(1.0D0-(4.0D0-PI)*ETA_R**2.0D0+1.0D0*ETA_T*ETA_L)/MASSR &
                    &+4.0D0*PI*FRE_R*DAMPING_R*(2.0D0*DIS0_ANG/DT+VEL0_ANG)+4.0D0*DIS0_ANG/DT/DT+4.0D0*VEL0_ANG/DT+ACC0_ANG

                DIS_ANG=DXRF/DXRK
                VEL_ANG=2.0D0/DT*(DIS_ANG-DIS0_ANG)-VEL0_ANG
                ACC_ANG=4.0D0/DT/DT*(DIS_ANG-DIS0_ANG)-4.0D0/DT*VEL0_ANG-ACC0_ANG
                !==================================dis_ang=======================================================================

                IF(MOD(NSTEP,5)==0) THEN
                    !WRITE(888,'(4E18.10)')T,DIS_ANG,VEL_ANG,ACC_ANG
                    WRITE(666,'(4E18.10)')T,DIS_Y,VEL_Y,ACC_Y
                END IF

            END IF
        END IF
        
        IF(IROTA==2201)THEN
            IF (T>T_START)  THEN
                !下游圆柱
                !==================================dis_y==========================================================================
                DXYK=(2.0D0*PI*FRE_Y)**2+8.0D0*PI*FRE_Y*DAMPING_Y/DT+4.0D0/DT/DT    !等效刚度
                DXYF=CL2*2/PI/MASSR+4.0D0*PI*FRE_Y*DAMPING_Y &
                *(2.0D0*DIS0_Y/DT+VEL0_Y)+4.0D0*DIS0_Y/DT/DT+4.0D0*VEL0_Y/DT+ACC0_Y    !等效荷载
                DIS_Y=DXYF/DXYK
                VEL_Y=2.0D0/DT*(DIS_Y-DIS0_Y)-VEL0_Y
                ACC_Y=4.0D0/DT/DT*(DIS_Y-DIS0_Y)-4.0D0/DT*VEL0_Y-ACC0_Y             
                !==================================dis_y==========================================================================
 
                
                !上游圆柱
                !==================================dis_y==========================================================================
                DXYK=(2.0D0*PI*FRE_Y)**2+8.0D0*PI*FRE_Y*DAMPING_Y/DT+4.0D0/DT/DT    !等效刚度
                DXYF=CL1*2/PI/MASSR+4.0D0*PI*FRE_Y*DAMPING_Y &
                *(2.0D0*DIS0_Y_UP/DT+VEL0_Y_UP)+4.0D0*DIS0_Y_UP/DT/DT+4.0D0*VEL0_Y_UP/DT+ACC0_Y_UP    !等效荷载
                DIS_Y_UP=DXYF/DXYK
                VEL_Y_UP=2.0D0/DT*(DIS_Y_UP-DIS0_Y_UP)-VEL0_Y_UP
                ACC_Y_UP=4.0D0/DT/DT*(DIS_Y_UP-DIS0_Y_UP)-4.0D0/DT*VEL0_Y_UP-ACC0_Y_UP             
                !==================================dis_y==========================================================================
                
 
                IF(MOD(NSTEP,5)==0) THEN
                WRITE(666,'(4E18.10)')T,DIS_Y,VEL_Y,ACC_Y
                WRITE(6661,'(4E18.10)')T,DIS_Y_UP,VEL_Y_UP,ACC_Y_UP
                END IF

		    END IF
        END IF

        IF(IROTA==2101)THEN
            IF (T>T_START)  THEN
                !下游圆柱
                !==================================dis_y==========================================================================
                DXYK=(2.0D0*PI*FRE_Y)**2+8.0D0*PI*FRE_Y*DAMPING_Y/DT+4.0D0/DT/DT    !等效刚度
                DXYF=CL2*2/PI/MASSR+4.0D0*PI*FRE_Y*DAMPING_Y &
                    *(2.0D0*DIS0_Y/DT+VEL0_Y)+4.0D0*DIS0_Y/DT/DT+4.0D0*VEL0_Y/DT+ACC0_Y    !等效荷载
                DIS_Y=DXYF/DXYK
                VEL_Y=2.0D0/DT*(DIS_Y-DIS0_Y)-VEL0_Y
                ACC_Y=4.0D0/DT/DT*(DIS_Y-DIS0_Y)-4.0D0/DT*VEL0_Y-ACC0_Y
                !==================================dis_y==========================================================================


                !上游圆柱固定
                !==================================dis_y==========================================================================
                
                !==================================dis_y==========================================================================


                IF(MOD(NSTEP,5)==0) THEN
                    WRITE(666,'(4E18.10)')T,DIS_Y,VEL_Y,ACC_Y
                END IF

            END IF
        END IF
        

        IF((IROTA == 510 ).Or.(IROTA == 5102 ).OR.(IROTA==5101).OR.(IROTA==520))THEN        !5107 For Flexible Forced Oscillation Test
            !510 different inlet; 5101 fixed cylinder oscillating plate;
            IF ((T>T_START).AND.(CSD_METHOD == 1))  THEN
                CALL FLEXIBLEMAIN
            ELSEIF((T>T_START).AND.(CSD_METHOD == 2)) THEN
                CALL ANCFMAIN
            END IF
        END IF

        If(T>T_Start)Then
            IF (MESH_TRANSFORM==1) THEN !MeshUpdate
                CALL DISPLACEMENT_ALE_BC
                BP_ALE =0.0D0
                APX_ALE=AP_ALE
                X1_ALE =X0_ALE
                CALL APPLYING_BC1( APX_ALE,BP_ALE,XBC_ALE )
                CALL SOLUTION_BICGS_ALE(APX_ALE,BP_ALE,X1_ALE)
                BP_ALE =0.0D0
                APY_ALE=AP_ALE
                Y1_ALE =Y0_ALE
                CALL APPLYING_BC1( APY_ALE,BP_ALE,YBC_ALE )
                CALL SOLUTION_BICGS_ALE(APY_ALE,BP_ALE,Y1_ALE)
            ELSEIF(MESH_TRANSFORM==2)THEN
                CALL SPRING_METHOD
            ELSE
                CALL system_clock(count = int_CPU_START)
                CALL RBFMeshUpdate
                CALL system_clock(count = int_CPU_FINISH)
                Write(*,'(A10,f12.6)') 'RBF Time =',(int_CPU_FINISH-int_CPU_START)/10000.0D0!不清楚为什么会差10倍
            END IF
        End If

        IF ( T>T_START ) THEN
            DIS0_ANG=DIS_ANG
            VEL0_ANG=VEL_ANG
            ACC0_ANG=ACC_ANG

            DIS0_X=DIS_X
            VEL0_X=VEL_X
            ACC0_X=ACC_X

            DIS0_Y=DIS_Y
            VEL0_Y=VEL_Y
            ACC0_Y=ACC_Y
            
            DIS0_Y_UP=DIS_Y_UP
            VEL0_Y_UP=VEL_Y_UP
            ACC0_Y_UP=ACC_Y_UP 

            tor0=tor
            CL0=CL
            CD0=CD
        ELSE
            X1_ALE = X0_ALE
            Y1_ALE = Y0_ALE
        END IF

        UM_ALE=0.0D0
        VM_ALE=0.0D0

        UM_ALE(1:NN_ALE)=( X1_ALE(1:NN_ALE)-X0_ALE(1:NN_ALE) )/DT
        VM_ALE(1:NN_ALE)=( Y1_ALE(1:NN_ALE)-Y0_ALE(1:NN_ALE) )/DT

        X0_ALE=X1_ALE
        Y0_ALE=Y1_ALE

        X(1:NN_ALE)=X1_ALE(1:NN_ALE)
        Y(1:NN_ALE)=Y1_ALE(1:NN_ALE)
        
        IF (MOD(NSTEP,NW)==0.OR.(NSTEP>=NW_START.and.NSTEP<=NW_END.and.MOD(NSTEP,nw)==0)) THEN
            CALL BINRESULT
        END IF


        CALL DYNAMIC_DT

        WRITE(*,'(A54)') '******************************************************'
        WRITE(*,*)
        NSTEP=NSTEP+1

    END DO
    STOP
    END PROGRAM CIRCULARCYLINDER_WITH_RIGIDPLATE



    SUBROUTINE README
    USE CONSTANT
    USE VARIABLE
    Use FlexibleOnly
    Use ANCF
    IMPLICIT NONE

    OPEN(UNIT=10,FILE='README.TXT')
    READ(10,*) TIME_STOP              !JY!结束时间
    READ(10,*) DT000
    READ(10,*) TIME_CALL_FORCE        !JY!计算受力的开始时间
    READ(10,*) SPRING_OUT, SPRING_IN  !JY!弹簧法外侧和内侧的半径
    READ(10,*) VISCOSITY
    READ(10,*) CS
    READ(10,*) RIGHTLAYER, LEFTLAYER, UPLAYER, DOWNLAYER
    READ(10,*) HH
    READ(10,*) MASSR
    READ(10,*) FRE_X,DAMPING_X
    READ(10,*) FRE_Y,DAMPING_Y
    READ(10,*) FRE_R,DAMPING_R
    READ(10,*) IROTA
    READ(10,*) AMC_ALE
    READ(10,*) OMG_ALE
    READ(10,*) NN_ORIG
    READ(10,*) NE_ORIG
    READ(10,*) MESH_TRANSFORM,CSD_METHOD
    READ(10,*) ITER_NUM
    READ(10,*) T_START
    READ(10,*) NW,NW_START,NW_END
    READ(10,*) CFL
    READ(10,*) FA_SPR
    READ(10,*) FUSAI_SPR
    READ(10,*) RE_START
    READ(10,*) FILE_NAME
    READ(10,*) LENGTH_LMM
    READ(10,*) MASS_ADD
    READ(10,*) F_0
    READ(10,*) F_1
    READ(10,*) T_LMM
    READ(10,*) LL0
    READ(10,*) EPSILON_ANG
    READ(10,*) CL_CUT
    READ(10,*) MESHFILENAME
    READ(10,*) ETA_R,ETA_T,ETA_L
    READ(10,*) AFP%E,AFP%A,AFP%I,AFP%RHO
    READ(10,*) AFP%NE
    READ(10,*) DampingCoeiff
    READ(10,*) DS_R,DS_L,DS_U,DS_D
    READ(10,*) R_pg,L_pg,U_pg,D_pg

    CLOSE(10)

    IF(IROTA==11)THEN
        AMC_ALE=(AMC_ALE/180.0D0)*(DACOS(-1.0D0))
        OMG_ALE=OMG_ALE/(2.0D0*DSIN(AMC_ALE))
    END IF

    IF(IROTA==1)THEN
        AMC_ALE=(AMC_ALE/180.0D0)*(DACOS(-1.0D0))       !转化成弧度
    END IF

    OMG_ALE=2.0D0*(DACOS(-1.0D0))*OMG_ALE               !w=2*pi*fn

    FRE_X=1.0D0/FRE_X
    FRE_Y=1.0D0/FRE_Y
    FRE_R=1.0D0/FRE_R
    VISCOSITY=1.0D0/VISCOSITY
    
    END SUBROUTINE README