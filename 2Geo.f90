
    subroutine grid_generation_orig
    use constant
    use VARIABLE
    implicit none

    integer:: i,j,k,ii,jj,kk,num
    real(8):: pi,da,ang0,cdr,dr,dr0,lensq,lensp
    real(8),allocatable:: angle(:),distance(:)

    real(8):: small,eps,cmax,cmin,ss
    real(8),allocatable::rr(:)


    real(8):: scale_lensp_cr0=4.0D0, cr0=0.5D0,rr00=0.5d0
    real(8):: xl,xr,yu,yd

    pi=acos(-1.0d0)

    lensp=scale_lensp_cr0*CR0

    lensq=2.0D0* (cr0+lensp)
    if(lensq<lensp) then
        write(*,*) 'hy99'
        PAUSE
        stop
    end if

    xl=-(cr0+lensp)
    xr= (cr0+lensp)
    yu= (cr0+lensp)
    yd=-(cr0+lensp)

    i=mod(ncc,8)
    if(i/=0) then
        write(*,*) ' sq22'
        PAUSE
        stop
    end if


    allocate(angle(ncc))
    angle=0.0d0

    da=(0.25d0*pi)/float(ncc/8)

    ang0=0.0d0
    num=1

    do i=1,ncc/8
        num=num+1
        angle(num)=ang0+float(i)*da
    end do
    angle(num)=pi/4.0d0
    ang0=pi/4.0d0

    da=(0.5d0*pi)/float(ncc/8*2)
    do i=1,(ncc/8)*2
        num=num+1
        angle(num)=ang0+float(i)*da
    end do
    angle(num)=(3.0d0/4.0d0)*pi
    ang0=(3.0D0/4.0d0)*pi

    do i=1,(ncc/8)*2
        num=num+1
        angle(num)=ang0+float(i)*da
    end do
    angle(num)=(5.0d0/4.0d0)*pi
    ang0=(5.0d0/4.0d0)*pi

    do i=1,(ncc/8)*2
        num=num+1
        angle(num)=ang0+float(i)*da
    end do
    angle(num)=(7.0d0/4.0d0)*pi
    ang0=(7.0d0/4.0d0)*pi

    da=((2.0d0*pi)-(7.0d0*pi/4.0d0))/float(ncc/8)
    do i=1,(ncc/8)-1
        num=num+1
        angle(num)=ang0+float(i)*da
    end do
    allocate(distance(ncc))

    do i=1,ncc
        ang0=angle(i)

        if (ang0<=0.25d0*pi ) then
            distance(i)=0.5d0*lensq/cos(ang0)

        else if	(ang0>0.25d0*pi.and.ang0<=0.75d0*pi ) then
            distance(i)=0.5d0*lensq/sin(ang0)

        else if (ang0>0.75d0*pi.and.ang0<=1.25d0*pi) then
            distance(i)=-0.5d0*lensq/cos(ang0)

        else if (ang0>1.25d0*pi.and.ang0<=1.75d0*pi) then

            distance(i)=-0.5d0*lensq/sin(ang0)
        else

            distance(i)=0.5d0*lensq/cos(ang0)

        end if
    end do


    do i=1,ncc
        distance(i)=distance(i)-cr0
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dr0= cr0 * ( 2.0d0*pi/float(ncc) )
    dr0=dr0*0.25d0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if(dr0*float(ndr)>(lensq/2.0d0-cr0)) then
        dr0=((lensq/2.0d0)-cr0)/float(ndr)
        dr0=dr0/2.0d0
    end if
    allocate(rr(ndr+1))
    rr=0.0d0

    nn_orig=(ncc)*(ndr+1)

    allocate( x_orig( nn_orig ), y_orig( nn_orig ) )
    x_orig=0.0d0
    y_orig=0.0d0

    num=0
    do i=1,ncc

        cmin=0.0d0;	cdr=cmin;	dr=dr0;		kk=0
        do
            kk=kk+1
            if(kk>1000) then
                write(*,*) 'de55';		PAUSE;stop
            end if

            ss=0.0d0
            dr=dr0
            do k=1,ndr
                ss=ss+dr;	dr=dr*cdr
            end do

            if(ss<=distance(i)) then
                cmax=cmax+1.0d0;	cdr=cmax
            else
                exit
            end if
        end do

        kk=0;	cdr=cmax
        do
            kk=kk+1
            if(kk>10000) then
                write(*,*) ' jy00'
                PAUSE
                stop
            end if

            ss=0.0d0;		dr=dr0
            do k=1,ndr
                ss=ss+dr;		dr=cdr*dr
            end do
            small=(ss-distance(i))

            if(dabs(small)<=1.0e-8) then
                exit
            else
                if(small>0.0d0) then
                    cmax=cdr;		cdr=0.5d0*(cmax+cmin)
                else
                    cmin=cdr;		cdr=0.5d0*(cmax+cmin)
                end if
            end if
        end do

        if(i==1.or.i==ncc+1) then
            rr(1)=cr0*cos(angle(i))
        else
            rr(1)=cr0
        end if

        dr=dr0
        do k=2,ndr
            rr(k)=rr(k-1)+dr
            dr=dr*cdr
        end do

        if(i==1.or.i==ncc+1) then
            rr(ndr+1)=xr
        else
            rr(ndr+1)=distance(i)+cr0
        end if

        do j=1,ndr+1
            num=num+1
            x_orig(num)=rr(j)*cos(angle(i))
            y_orig(num)=rr(j)*sin(angle(i))
        end do
    end do


    ne_orig=ncc*ndr

    allocate(nodeele_orig(ne_orig,4))
    nodeele_orig=0
    num=0
    do i=1,ncc-1
        do j=1,ndr
            num=num+1
            k=(i-1)*(ndr+1)+j
            nodeele_orig(num,1)=k
            nodeele_orig(num,2)=k+1

            nodeele_orig(num,3)=k+1+(ndr+1)
            nodeele_orig(num,4)=k+(ndr+1)
        end do
    end do

    k=(ndr+1)*(ncc-1)
    do i=1,ndr

        num=num+1
        nodeele_orig(num,1)=k+i
        nodeele_orig(num,2)=k+1+i

        nodeele_orig(num,3)=i+1
        nodeele_orig(num,4)=i
    end do



    write(*,*) ' nodes   :',nn_orig
    write(*,*) ' elements:',ne_orig
    write(*,*) ' *************** ============== ***************'

    !========================================================================           ! ALE
    NN_ALE=NN_ORIG                                                      !
    NE_ALE=NE_ORIG                                                      !
    ALLOCATE(X_ALE(NN_ALE),Y_ALE(NN_ALE),NODEELE_ALE(NE_ALE,4))        !
    X_ALE=X_ORIG                                                        !
    Y_ALE=Y_ORIG                                                        !
    NODEELE_ALE=NODEELE_ORIG                                            !
    !=======================================================================!
    i=0
    if(i==1) then

        open(unit=10,file='mesh_orig.dat')
        WRITE(10,'(A15)') ' mesh_orig'
        WRITE(10,*) 'TITLE=" LL " '
        WRITE(10,*) 'VARIABLES="X","Y" '
        WRITE(10,'(1X,A7,I6,A4,I6,A31)') 'ZONE N=',NN_orig, ', E=',NE_orig,',F=FEPOINT, ET=QUADRILATERAL'
        WRITE(10,*) 'DT=(SINGLE SINGLE )'
        DO I=1,NN_orig
            WRITE(10,*)	X_orig(I),Y_orig(I)
        END DO
        DO I=1,NE_orig
            WRITE(10,*) (NODEELE_orig(I,J),J=1,4)
        END DO
        CLOSE(10)
    end if

    End subroutine grid_generation_orig
    !
    !
    !
    !
    !
    SUBROUTINE SURROUND_ELE_ORIG
    USE CONSTANT
    USE VARIABLE
    IMPLICIT NONE

    INTEGER:: I,J,K,NODE,NUM
    INTEGER:: ELE_ELE(30),M,N

    Maxnum_Ele_Adj=30

    ALLOCATE(ELE_ADJ(NE_orig,Maxnum_Ele_Adj))
    ALLOCATE(ARRAY_LARGE(NN_orig))
    ALLOCATE(NODE_ELE(NN_orig,30))

    ARRAY_LARGE=0
    NODE_ELE=0
    ELE_ADJ=0

    DO I=1,NE_orig     !k:节点出现的次数    NODE_ELE所属单元
        DO J=1,4
            NODE=NODEELE_orig(I,J)
            ARRAY_LARGE(NODE)=ARRAY_LARGE(NODE)+1
            K=ARRAY_LARGE(NODE)
            IF (K > 30) THEN
                WRITE(*,*) 'MORE ADJ ELES, ERROR IN SURROUND ELE'
                WRITE(*,*) 'ADJ ELES AROUND THE 4 NODES OF ELE:',K
                PAUSE
                STOP
            END IF
            NODE_ELE(NODE,K)=I
        END DO
    END DO

    DO I=1,NE_orig
        NUM=0
        ELE_ELE=0
        DO J=1,4
            NODE=NODEELE_orig(I,J)
            DO K=1,30
                M=NODE_ELE(NODE,K)
                IF (M.NE.0) THEN
                    NUM=NUM+1
                    ELE_ELE(NUM)=M
                ELSE
                    EXIT
                END IF
            END DO
        END DO


        DO J=1,30
            IF (ELE_ELE(J).NE.0) THEN
                IF (ELE_ELE(J)==I) THEN
                    ELE_ELE(J)=0
                END IF
            ELSE
                EXIT
            END IF
        END DO

        DO J=1,29
            M=ELE_ELE(J)
            IF (M.NE.0) THEN
                DO K=J+1,30
                    N=ELE_ELE(K)
                    IF (M==N) THEN
                        ELE_ELE(K)=0
                    END IF
                END DO
            END IF
        END DO

        N=0
        DO J=1,30
            M=ELE_ELE(J)
            IF (M.NE.0) THEN
                N=N+1
                IF (N > Maxnum_Ele_Adj) THEN
                    WRITE(*,*) 'MORE IN SUBROUTNE SURROUND ELE'
                    PAUSE
                    STOP
                END IF
                ELE_ADJ(I,N)=M
            END IF
        END DO

    END DO
    DEALLOCATE(ARRAY_LARGE)

    END SUBROUTINE SURROUND_ELE_ORIG
    !
    !
    !
    !
    !
    subroutine polybn_0226(indexa)
    use constant
    use VARIABLE
    implicit none

    integer:: i,j,k,ii,jj,kk,num,node1,node2,node3,node4
    integer:: e1(100),marker,nmarker,nsl,ele,m,indexa
    real(8):: x1,y1,x2,y2,small,xl,xr,yu,yd

    small=1.0d-6

    xl=minval(x_orig)
    xr=maxval(x_orig)
    yu=maxval(y_orig)
    yd=minval(y_orig)

    WRITE(*,*)
    WRITE(*,*) ' INDEXA:',INDEXA

    ALLOCATE(BE_ORIG(NE_ORIG,3))
    be_orig=0
    nbe_orig=0

    do i=1,ne_orig
        num=0
        do j=1,Maxnum_Ele_Adj
            k=ele_adj(i,j)          !i单元周围的单元k

            if(k/=0) then
                num=num+1
                e1(num)=k               !num周围单元数
            end if
            if(j<Maxnum_Ele_Adj.and.ele_adj(i,j+1)==0) exit
        end do

        if(num==0) then
            write(*,*) ' cna'
            PAUSE
            stop
        end if

        do j=1,4
            marker=0
            nmarker=0
            if(j<=3) then
                node1=nodeele_orig(i,j)
                node2=nodeele_orig(i,j+1)
            else
                node1=nodeele_orig(i,4)
                node2=nodeele_orig(i,1)
            end if

            do k=1,num
                kk=e1(k)
                do ii=1,4
                    if(ii<=3) then
                        node3=nodeele_orig(kk,ii)
                        node4=nodeele_orig(kk,ii+1)
                    else
                        node3=nodeele_orig(kk,4)
                        node4=nodeele_orig(kk,1)
                    end if

                    if( ((node1==node3).and.(node2==node4)).or.((node1==node4).and.(node2==node3)) ) then
                        !    &	                有重合   则不是边界
                        marker=1
                    end if
                end do
                if(marker==1) exit
            end do

            if(marker==0) then
                nbe_orig=nbe_orig+1                   !------XHL----------边界的边，则边界单元+1
                be_orig(nbe_orig,1)=i                 !单元
                be_orig(nbe_orig,2)=j                 !1-4     第j个单元节点
            end if
        end do
    end do

    WRITE(*,*)
    WRITE(*,*) 'NUM OF BE=',NBE_ORIG
    !*	----------------------------------------------------------------
    NLBE_ORIG=0
    NUBE_ORIG=0
    NRBE_ORIG=0
    NDBE_ORIG=0
    NPBE_ORIG=0
    NSBE_ORIG=0

    DO I=1,NBE_ORIG
        ELE=BE_ORIG(I,1)
        nsl=be_orig(i,2)

        if(nsl<=3) then
            node1=nodeele_orig(ele,nsl  )
            node2=nodeele_orig(ele,nsl+1)
        else
            node1=nodeele_orig(ele,4)
            node2=nodeele_orig(ele,1)
        end if

        X1=X_ORIG(NODE1)
        Y1=Y_ORIG(NODE1)
        X2=X_ORIG(NODE2)
        Y2=Y_ORIG(NODE2)

        IF		((ABS(X1-XL)<SMALL).AND.(ABS(X2-XL)<SMALL).AND.(ABS(X1-X2)<SMALL) )THEN
            !     &
            NLBE_ORIG=NLBE_ORIG+1
            BE_ORIG(I,3)=1 ! 1:	L
            if(indexa==1) then
                nodetype(node1)=1             !找到ALE边界，在边弹簧法中指定边界条件
                nodetype(node2)=1
            end if
            !*
        ELSE IF ((ABS(Y1-YU)<SMALL).AND.(ABS(Y2-YU)<SMALL).AND. (ABS(Y1-Y2)<SMALL)) THEN

            NUBE_ORIG=NUBE_ORIG+1
            BE_ORIG(I,3)=2 ! 2: U
            if(indexa==1) then
                nodetype(node1)=1
                nodetype(node2)=1
            end if
            !*
        ELSE IF ((ABS(X1-XR)<SMALL).AND.(ABS(X2-XR)<SMALL).AND.(ABS(X1-X2)<SMALL)) THEN

            NRBE_ORIG=NRBE_ORIG+1
            BE_ORIG(I,3)=3 ! 3: R
            if(indexa==1) then
                nodetype(node1)=1
                nodetype(node2)=1
            end if

        ELSE IF ((ABS(Y1-YD)<SMALL).AND.(ABS(Y2-YD)<SMALL).AND.(ABS(Y1-Y2)<SMALL)) THEN

            NDBE_ORIG=NDBE_ORIG+1
            BE_ORIG(I,3)=4 ! 4: D
            if(indexa==1) then
                nodetype(node1)=1
                nodetype(node2)=1
            end if
            !*
        ELSE IF (IROTA==2201 .OR. IROTA==2101 .OR. IROTA==0)THEN
            IF ((dabs(x1*x1+y1*y1-0.25d0)<small).and.(dabs(x2*x2+y2*y2-0.25d0)<small)) then
        
            NpBE_ORIG=NpBE_ORIG+1   !上游圆柱
            BE_ORIG(I,3)=5 ! 5:	P
            if(indexa==1) then
                nodetype(node1)=1
                nodetype(node2)=1
            end if
        
            ELSE 
        
            nsbe_orig=nsbe_orig+1    !下游圆柱
            be_orig(i,3)=6 ! 6: S
            if(indexa==1) then
                nodetype(node1)=1
                nodetype(node2)=1
            end if
        
            END IF
            
        ELSE IF (IROTA==200 .OR. IROTA==201 .OR. IROTA==510 .OR. IROTA==5102 .OR. Irota==5101 .OR. IROTA==520)THEN    
            IF((X1-0.50D0)>SMALL.AND.(X2-0.50D0)>SMALL.AND.(ABS(X1)<5.0D0).AND.(ABS(X2)<5.0D0).AND.(ABS(Y1)<0.50D0).AND.(ABS(Y2)<0.50D0))THEN

            NSBE_ORIG=NSBE_ORIG+1
            BE_ORIG(I,3)=6 ! 6: S
            !write(*,*) x1,y1
            !write(*,*) x2,y2
            IF(INDEXA==1) THEN
                NODETYPE(NODE1)=1
                NODETYPE(NODE2)=1
            END IF

            ELSE
            NPBE_ORIG=NPBE_ORIG+1
            BE_ORIG(I,3)=5 ! 5:	P
            !write(*,*) x1,y1
            !write(*,*) x2,y2
            IF(INDEXA==1) THEN
                NODETYPE(NODE1)=1
                NODETYPE(NODE2)=1
            END IF
            END IF
            
         END IF
        !ELSE IF ((dabs(x1*x1+y1*y1-0.25d0)<small).and.(dabs(x2*x2+y2*y2-0.25d0)<small)) then
        !
        !    NpBE_ORIG=NpBE_ORIG+1
        !    BE_ORIG(I,3)=5 ! 5:	P
        !    if(indexa==1) then
        !        nodetype(node1)=1
        !        nodetype(node2)=1
        !    end if
        !
        !ELSE
        !
        !    nsbe_orig=nsbe_orig+1
        !    be_orig(i,3)=6 ! 6: S
        !    if(indexa==1) then
        !        nodetype(node1)=1
        !        nodetype(node2)=1
        !    end if
        !
        !END IF
    END DO

    WRITE(*,*)
    WRITE(*,*) ' NLBE_ORIG=',NLBE_ORIG
    WRITE(*,*) ' NUBE_ORIG=',NUBE_ORIG
    WRITE(*,*) ' NRBE_ORIG=',NRBE_ORIG
    WRITE(*,*) ' NDBE_ORIG=',NDBE_ORIG
    WRITE(*,*) ' NPBE_ORIG=',NPBE_ORIG
    WRITE(*,*) ' NSBE_ORIG=',NSBE_ORIG

    WRITE(*,*) ' SUM   BES=',NLBE_ORIG+NUBE_ORIG+NRBE_ORIG+NDBE_ORIG+NPBE_ORIG+NSBE_ORIG

    WRITE(*,*)

    IF	(ALLOCATED(LBE_ORIG))	DEALLOCATE(LBE_ORIG)
    IF	(ALLOCATED(RBE_ORIG))	DEALLOCATE(RBE_ORIG)
    IF	(ALLOCATED(UBE_ORIG))	DEALLOCATE(UBE_ORIG)
    IF	(ALLOCATED(PBE_ORIG))	DEALLOCATE(PBE_ORIG)
    IF	(ALLOCATED(DBE_ORIG))	DEALLOCATE(DBE_ORIG)
    IF	(ALLOCATED(SBE_ORIG))	DEALLOCATE(SBE_ORIG)


    ALLOCATE(LBE_ORIG(NLBE_ORIG,2))
    ALLOCATE(UBE_ORIG(NUBE_ORIG,2))
    ALLOCATE(RBE_ORIG(NRBE_ORIG,2))
    ALLOCATE(DBE_ORIG(NDBE_ORIG,2))
    ALLOCATE(PBE_ORIG(NPBE_ORIG,2))
    ALLOCATE(SBE_ORIG(NSBE_ORIG,2))

    nlbe_orig=0
    nube_orig=0
    nrbe_orig=0
    ndbe_orig=0
    npbe_orig=0
    nsbe_orig=0

    do i=1,nbe_orig
        j=be_orig(i,3)
        if(j==1) then
            nlbe_orig=nlbe_orig+1
            lbe_orig(nlbe_orig,:)=be_orig(i,1:2)
        else if (j==2) then
            nube_orig=nube_orig+1
            ube_orig(nube_orig,:)=be_orig(i,1:2)
        else if (j==3) then
            nrbe_orig=nrbe_orig+1
            rbe_orig(nrbe_orig,:)=be_orig(i,1:2)
        else if (j==4) then
            ndbe_orig=ndbe_orig+1
            dbe_orig(ndbe_orig,:)=be_orig(i,1:2)
        else if (j==5) then
            npbe_orig=npbe_orig+1
            pbe_orig(npbe_orig,:)=be_orig(i,1:2)
        else if (j==6) then
            nsbe_orig=nsbe_orig+1
            sbe_orig(nsbe_orig,:)=be_orig(i,1:2)

        end if
    end do

    !C	----------------------------------------------------------------
    !C     CLASSFYING BOUNDARY NODES
    !C	----------------------------------------------------------------
    IF(ALLOCATED(I4ARRAY)) DEALLOCATE(I4ARRAY)

    ALLOCATE(I4ARRAY(NLBE_ORIG*2))
    i4array=0
    DO I=1,NLBE_ORIG
        J=LBE_ORIG(I,1)                              !单元
        k=LBE_ORIG(I,2)                              !节点
        IF(k<=3) THEN
            NODE1=k;		NODE2=k+1
        ELSE
            NODE1=4;		NODE2=1
        END IF

        I4ARRAY(2*I-1)=nodeele_ORIG(j,NODE1)         !节点
        I4ARRAY(2*I  )=nodeele_ORIG(j,NODE2)
    END DO

    nlbn_orig=0
    DO I=1,NLBE_ORIG*2-1
        K=I4ARRAY(I)         !节点               去重复
        if(k/=0) then
            nlbn_orig=nlbn_orig+1
            DO J=I+1,NLBE_ORIG*2
                M=I4ARRAY(J)
                IF (K==M) THEN
                    I4ARRAY(J)=0
                END IF
            END DO
        end if
    END DO

    if(i4array(nlbe_orig*2)/=0) then
        nlbn_orig=nlbn_orig+1
    end if

    if(allocated(lbn_orig)) deallocate(lbn_orig)
    allocate(lbn_orig(nlbn_orig))

    lbn_orig=0

    NLBN_ORIG=0
    DO I=1,NLBE_ORIG*2
        J=I4ARRAY(I)
        IF (J/=0) THEN
            NLBN_ORIG=NLBN_ORIG+1
            lbn_orig(nlbn_orig)=i4array(i)
        END IF
    END DO
    deallocate(i4array)
    !**********************************************************************************
    !**********************************************************************************
    allocate(i4array(nube_orig*2))
    i4array=0
    DO I=1,NuBE_ORIG
        J=uBE_ORIG(I,1)
        k=uBE_ORIG(I,2)
        IF(k<=3) THEN
            NODE1=k;		NODE2=k+1
        ELSE
            NODE1=4;		NODE2=1
        END IF

        i4array(2*I-1)=nodeele_ORIG(j,NODE1)
        i4array(2*I  )=nodeele_ORIG(j,NODE2)
    END DO

    nubn_orig=0
    DO I=1,NuBE_ORIG*2-1
        K=I4ARRAY(I)
        if(k/=0) then
            nubn_orig=nubn_orig+1
            DO J=I+1,NuBE_ORIG*2
                M=I4ARRAY(J)
                IF (K==M) THEN
                    i4array(J)=0
                END IF
            END DO
        end if
    END DO

    if(i4array(nube_orig*2)/=0) then
        nubn_orig=nubn_orig+1
    end if

    if (allocated(ubn_orig)) deallocate(ubn_orig)
    allocate(ubn_orig(nubn_orig))

    ubn_orig=0

    NuBN_ORIG=0
    DO I=1,NuBE_ORIG*2
        J=I4ARRAY(I)
        IF (J/=0) THEN
            NuBN_ORIG=NuBN_ORIG+1
            ubn_orig(nubn_orig)=i4array(i)
        END IF
    END DO
    deallocate(i4array)
    !*
    !*
    allocate(i4array(nrbe_orig*2))
    i4array=0
    DO I=1,NrBE_ORIG
        J=rBE_ORIG(I,1)
        k=rBE_ORIG(I,2)
        IF(k<=3) THEN
            NODE1=k;		NODE2=k+1
        ELSE
            NODE1=4;		NODE2=1
        END IF

        i4array(2*I-1)=nodeele_ORIG(j,node1)
        i4array(2*I  )=nodeele_ORIG(j,node2)
    END DO

    nrbn_orig=0
    DO I=1,nrbe_orig*2-1
        K=i4array(I)
        if(k/=0) then
            nrbn_orig=nrbn_orig+1
            DO J=I+1,nrbe_orig*2
                M=i4array(J)
                IF (K==M) THEN
                    i4array(J)=0
                END IF
            END DO
        end if
    END DO

    if(i4array(nrbe_orig*2)/=0) then
        nrbn_orig=nrbn_orig+1
    end if

    if(allocated(rbn_orig)) deallocate(rbn_orig)
    allocate(rbn_orig(nrbn_orig))

    rbn_orig=0

    nrbn_orig=0
    DO I=1,nrbe_orig*2
        J=i4array(I)
        IF (J/=0) THEN
            nrbn_orig=nrbn_orig+1
            rbn_orig(nrbn_orig)=i4array(i)
        END IF
    END DO
    deallocate(i4array)
    !*
    !*
    allocate(i4array(ndbe_orig*2))
    i4array=0
    DO I=1,ndbe_orig
        J=dbe_orig(I,1)
        k=dbe_orig(I,2)
        IF(k<=3) THEN
            node1=k;		node2=k+1
        ELSE
            node1=4;		node2=1
        END IF

        i4array(2*I-1)=nodeele_orig(j,node1)
        i4array(2*I  )=nodeele_orig(j,node2)
    END DO

    ndbn_orig=0
    DO I=1,ndbe_orig*2-1
        K=i4array(I)
        if(k/=0) then
            ndbn_orig=ndbn_orig+1
            DO J=I+1,ndbe_orig*2
                M=i4array(J)
                IF (K==M) THEN
                    i4array(J)=0
                END IF
            END DO
        end if
    END DO

    if(i4array(ndbe_orig*2)/=0) then
        ndbn_orig=ndbn_orig+1
    end if


    if(allocated(dbn_orig)) deallocate(dbn_orig)
    allocate(dbn_orig(ndbn_orig))
    dbn_orig=0

    ndbn_orig=0
    DO I=1,ndbe_orig*2
        J=i4array(I)
        IF (J/=0) THEN
            ndbn_orig=ndbn_orig+1
            dbn_orig(ndbn_orig)=i4array(i)
        END IF
    END DO
    deallocate(i4array)
    !*
    !*
    ALLOCATE(I4ARRAY(NpBE_ORIG*2))
    i4array=0
    DO I=1,NpBE_ORIG
        J=pBE_ORIG(I,1)      !单元
        k=pBE_ORIG(I,2)      !节点
        IF(k<=3) THEN
            NODE1=k;		NODE2=k+1
        ELSE
            NODE1=4;		NODE2=1
        END IF

        I4ARRAY(2*I-1)=nodeele_ORIG(j,NODE1)
        I4ARRAY(2*I  )=nodeele_ORIG(j,NODE2)
    END DO

    npbn_orig=0
    DO I=1,NpBE_ORIG*2-1
        K=I4ARRAY(I)
        if(k/=0) then
            npbn_orig=npbn_orig+1
            DO J=I+1,NpBE_ORIG*2
                M=I4ARRAY(J)
                IF (K==M) THEN
                    I4ARRAY(J)=0
                END IF
            END DO
        end if
    END DO

    if(i4array(npbe_orig*2)/=0) then
        npbn_orig=npbn_orig+1
    end if

    if (allocated(pbn_orig)) deallocate(pbn_orig)
    allocate(pbn_orig(npbn_orig))

    pbn_orig=0

    NpBN_ORIG=0
    DO I=1,NpBE_ORIG*2
        J=I4ARRAY(I)
        IF (J/=0) THEN
            NpBN_ORIG=NpBN_ORIG+1
            pbn_orig(npbn_orig)=i4array(i)
        END IF
    END DO
    deallocate(i4array)
    !*
    !*
    ALLOCATE(I4ARRAY(NsBE_ORIG*2))
    i4array=0
    DO I=1,NsBE_ORIG
        J=sBE_ORIG(I,1)
        k=sBE_ORIG(I,2)
        IF(k<=3) THEN
            NODE1=k;		NODE2=k+1
        ELSE
            NODE1=4;		NODE2=1
        END IF

        I4ARRAY(2*I-1)=nodeele_ORIG(j,NODE1)
        I4ARRAY(2*I  )=nodeele_ORIG(j,NODE2)
    END DO

    nsbn_orig=0
    DO I=1,NsBE_ORIG*2-1
        K=I4ARRAY(I)
        if(k/=0) then
            nsbn_orig=nsbn_orig+1
            DO J=I+1,NsBE_ORIG*2
                M=I4ARRAY(J)
                IF (K==M) THEN
                    I4ARRAY(J)=0
                END IF
            END DO
        end if
    END DO

    if(i4array(nsbe_orig*2)/=0) then
        nsbn_orig=nsbn_orig+1
    end if

    if (allocated(sbn_orig)) deallocate(sbn_orig)
    allocate(sbn_orig(nsbn_orig))

    sbn_orig=0

    NsBN_ORIG=0
    DO I=1,NsBE_ORIG*2
        J=I4ARRAY(I)
        IF (J/=0) THEN
            NsBN_ORIG=NsBN_ORIG+1
            sbn_orig(nsbn_orig)=i4array(i)
        END IF
    END DO
    deallocate(i4array)

    WRITE(*,*) ' NLBN_ORIG =', NLBN_ORIG
    WRITE(*,*) ' NUBN_ORIG =', NUBN_ORIG
    WRITE(*,*) ' NRBN_ORIG =', NRBN_ORIG
    WRITE(*,*) ' NDBN_ORIG =', NDBN_ORIG
    WRITE(*,*) ' NPBN_ORIG =', NPBN_ORIG
    WRITE(*,*) ' NSBN_ORIG =', NSBN_ORIG
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (indexa==1) then                               ! ALE
        nlbn_ALE=nlbn_orig
        nrbn_ALE=nrbn_orig
        nubn_ALE=nubn_orig
        ndbn_ALE=ndbn_orig
        npbn_ALE=npbn_orig
        nsbn_ALE=nsbn_orig

        allocate(lbn_ALE(nlbn_ALE),rbn_ALE(nrbn_ALE),dbn_ALE(ndbn_ALE),ubn_ALE(nubn_ALE),pbn_ALE(npbn_ALE),sbn_ALE(nsbn_ALE))
        lbn_ALE=lbn_orig
        rbn_ALE=rbn_orig
        ubn_ALE=ubn_orig
        dbn_ALE=dbn_orig
        pbn_ALE=pbn_orig
        sbn_ALE=sbn_orig

        nbn_ale=	nlbn_ALE+	nrbn_ALE+	nubn_ALE+	ndbn_ALE+	npbn_ALE+	nsbn_ALE
    end if

    End Subroutine polybn_0226
    !
    !
    !
    !
    SUBROUTINE EXTEND_MESH_RIGHT_0226
    USE CONSTANT
    USE VARIABLE
    !	USE IMSL, dbn=>kmcdf09
    IMPLICIT NONE

    INTEGER:: I,J,K,NET(4),rmre,m
    REAL(8):: DS0,ds,x1,y1,yu,yd,xr

    REAL(8),ALLOCATABLE:: RA(:),RB(:)
    INTEGER,ALLOCATABLE:: IPERM(:)

    xr=maxval(x_orig)
    yu=maxval(y_orig)
    yd=minval(y_orig)

    !DS0=(YU-YD)/(FLOAT(NRBE_ORIG+1))
    !ds=ds0
    ds0 = DS_R
    ds=DS_R

    ALLOCATE(NODEELE((NRBE_ORIG)*RIGHTlayer,4))
    ALLOCATE(RA(NRBN_ORIG),RB(NRBN_ORIG))
    ALLOCATE(IPERM(NRBN_ORIG),I4ARRAY(NRBN_ORIG))

    DO I=1,NRBN_ORIG
        J=RBN_ORIG(I)
        I4ARRAY(I)=J
        IPERM(I)=I
        RA(I)=Y_ORIG(J)
    END DO
    CALL DSVRGP (NRBN_ORIG, RA, RB, IPERM)

    !	do i=1,nrbn_orig
    !	write(*,*) ra(i),y_orig(i4array(iperm(i))),rb(iperm(i))
    !	end do
    !	pause
    !   stop

    DO I=1,NRBN_ORIG
        J=IPERM(I)
        RBN_ORIG(I)=I4ARRAY(J)
    END DO

    k=0
    DO I=1,RIGHTlayer
        DO J=1,NRBE_ORIG !+1
            k=k+1
            IF(I==1) THEN
                NODEELE(K,1)=RBN_ORIG(J)

                NODEELE(K,2)=RBN_ORIG(J+1)
            ELSE
                NODEELE(K,1)=NN_ORIG+(I-2)*NRBN_ORIG+J
                NODEELE(K,2)=NN_ORIG+(I-2)*NRBN_ORIG+J+1
            END IF
            NODEELE(K,3)=NN_ORIG+(I-1)*NRBN_ORIG+J+1
            NODEELE(K,4)=NN_ORIG+(I-1)*NRBN_ORIG+J
        END DO
    END DO

    DO I=1,NRBE_ORIG*RIGHTlayer
        DO J=1,4
            NET(J)=NODEELE(I,J)
        END DO
        NODEELE(I,1)=NET(4)
        NODEELE(I,2)=NET(3)
        NODEELE(I,3)=NET(2)
        NODEELE(I,4)=NET(1)
    END DO
    !*	----------------------------------------------------------
    DO I=1,NRBN_ORIG
        RBN_ORIG(I)=NN_ORIG+NRBN_ORIG*(RIGHTlayer-1)+I
    END DO

    ALLOCATE (x_temp(NN_ORIG+NRBN_ORIG*RIGHTlayer),y_temp(NN_ORIG+NRBN_ORIG*RIGHTlayer), nodeele_temp(NE_ORIG+NRBE_ORIG*RIGHTlayer,4))


    x_temp(1:NN_ORIG)=X_ORIG(1:NN_ORIG)
    y_temp(1:NN_ORIG)=Y_ORIG(1:NN_ORIG)
    nodeele_temp(1:NE_ORIG,:)=NODEELE_ORIG(1:NE_ORIG,:)

    DO I=1,RIGHTlayer
        xr=xr+ds
        DO  J=1,NRBN_ORIG
            K=NN_ORIG+(I-1)*NRBN_ORIG+J
            y_temp(K)=RB(J)
            x_temp(K)=XR
        END DO
        ds=ds*R_pg
        if(ds>5.0d0*ds0) ds=5.0d0*ds0
    END DO

    DO I=1,RIGHTlayer*NRBE_ORIG
        K=NE_ORIG+I
        nodeele_temp(K,:)=NODEELE(I,:)
    END DO


    NE_ORIG=NE_ORIG+NRBE_ORIG*RIGHTlayer
    NN_ORIG=NN_ORIG+NRBN_ORIG*RIGHTlayer

    DEALLOCATE(RA,RB,IPERM,NODEELE_ORIG,X_ORIG,Y_ORIG,I4ARRAY)
    ALLOCATE(NODEELE_ORIG(NE_ORIG,4),X_ORIG(NN_ORIG),Y_ORIG(NN_ORIG))

    NODEELE_ORIG=nodeele_temp
    X_ORIG=x_temp
    Y_ORIG=y_temp
    DEALLOCATE(nodeele_temp,x_temp,y_temp,NODEELE)


    I=0
    IF(I==1) THEN
        OPEN(UNIT=10,FILE='MESHEXTRIGTH.DAT')
        WRITE(10,*) ' TITLE=" MESHEXTRIGTH " '
        WRITE(10,*) ' VARIABLES="x","y"'
        WRITE(10,*) ' ZONE N=', NN_ORIG,', E=',NE_ORIG,',F=FEPOINT, ET=QUADRILATERAL'

        WRITE(10,*) ' DT=(SINGLE SINGLE )'
        DO I=1,NN_ORIG
            WRITE(10,*) X_ORIG(I),Y_ORIG(I)
        END DO
        DO I=1,NE_ORIG
            WRITE(10,*) (NodeEle_ORIG(I,J),J=1,4)
        END DO
        CLOSE(10)
        pause
        write(*,*) ' fgb5v'
        PAUSE
        stop
    END IF

    END SUBROUTINE EXTEND_MESH_RIGHT_0226



    SUBROUTINE EXTEND_MESH_LEFT_0226
    USE CONSTANT
    USE VARIABLE
    !	USE IMSL
    IMPLICIT NONE

    INTEGER:: I,J,K,NET(4),rmre,m
    REAL(8):: DS0,ds,x1,y1

    REAL(8),ALLOCATABLE:: RA(:),RB(:)
    INTEGER,ALLOCATABLE:: IPERM(:)

    real(8):: yu,yd,xl


    yu=maxval(y_orig)
    yd=minval(y_orig)
    XL=MINVAL(X_ORIG)


    !DS0=(YU-YD)/(FLOAT(NLBE_ORIG))
    !DS=-DS0
    DS0 = - DS_L
    DS = - DS_L

    ALLOCATE(NODEELE(NlBE_ORIG*leftlayer,4))
    ALLOCATE(RA(NlBN_ORIG),RB(NlBN_ORIG))
    ALLOCATE(IPERM(NlBN_ORIG),I4ARRAY(NlBN_ORIG))

    DO I=1,NLbN_ORIG
        J=lBN_ORIG(I)
        I4ARRAY(I)=J
        IPERM(I)=I
        RA(I)=Y_ORIG(J)
    END DO
    CALL DSVRGP (NlBN_ORIG, RA, RB, IPERM)

    DO I=1,NLBN_ORIG
        J=IPERM(I)
        lBN_ORIG(I)=I4ARRAY(J)
    END DO

    K=0

    DO I=1,LeftLayer
        DO J=1,NLBE_ORIG
            K=K+1
            IF(I==1) THEN
                NODEELE(K,1)=lBN_ORIG(J)

                NODEELE(K,2)=lBN_ORIG(J+1)
            ELSE
                NODEELE(K,1)=NN_ORIG+(I-2)*NlBN_ORIG+J
                NODEELE(K,2)=NN_ORIG+(I-2)*NlBN_ORIG+J+1
            END IF
            NODEELE(K,3)=NN_ORIG+(I-1)*NlBN_ORIG+J+1
            NODEELE(K,4)=NN_ORIG+(I-1)*NlBN_ORIG+J
        END DO
    END DO

    !*	----------------------------------------------------------

    DO I=1,NlBN_ORIG
        lBN_ORIG(I)=NN_ORIG+NlBN_ORIG*(leftlayer-1)+I
    END DO

    ALLOCATE (x_temp(NN_ORIG+NlBN_ORIG*leftlayer),y_temp(NN_ORIG+NlBN_ORIG*leftlayer),nodeele_temp(NE_ORIG+NlBE_ORIG*leftlayer,4))

    x_temp(1:NN_ORIG)=X_ORIG(1:NN_ORIG)
    y_temp(1:NN_ORIG)=Y_ORIG(1:NN_ORIG)
    nodeele_temp(1:NE_ORIG,:)=NODEELE_ORIG(1:NE_ORIG,:)

    DO I=1,leftlayer
        xL=xL+ds
        DO J=1,NlBN_ORIG
            K=NN_ORIG+(I-1)*NLBN_ORIG+J
            y_temp(K)=RB(J)
            x_temp(K)=XL
        END DO

        ds=ds*L_pg
        if(DABS(ds)>DABS(5.0d0*ds0)) ds=5.0d0*ds0

    END DO

    DO I=1,LeftLayer*NLBE_ORIG
        K=NE_ORIG+I
        nodeele_temp(K,:)=NODEELE(I,:)
    END DO


    NE_ORIG=NE_ORIG+NLBE_ORIG*LEFTlayer
    NN_ORIG=NN_ORIG+NLBN_ORIG*LEFTlayer

    DEALLOCATE(RA,RB,IPERM,NODEELE_ORIG,X_ORIG,Y_ORIG,I4ARRAY)
    ALLOCATE(NODEELE_ORIG(NE_ORIG,4),X_ORIG(NN_ORIG),Y_ORIG(NN_ORIG))

    NODEELE_ORIG=nodeele_temp
    X_ORIG=x_temp
    Y_ORIG=y_temp
    DEALLOCATE(nodeele_temp,x_temp,y_temp,NODEELE)

    DEALLOCATE(LBE_ORIG,RBE_ORIG,DBE_ORIG,UBE_ORIG,pbe_orig,sbe_orig)
    DEALLOCATE(LBN_ORIG,RBN_ORIG,DBN_ORIG,UBN_ORIG,pbn_orig,sbn_orig)
    DEALLOCATE(be_orig,ELE_ADJ,NODE_ELE)

    I=0
    IF(I==1) THEN

        OPEN(UNIT=10,FILE='MESHEXTLEFT.DAT')
        WRITE(10,*) ' TITLE=" MESHEXLEFT " '
        WRITE(10,*) ' VARIABLES="x","y"'
        WRITE(10,*) ' ZONE N=', NN_ORIG,', E=',NE_ORIG,', F=FEPOINT, ET=QUADRILATERAL'
        WRITE(10,*) ' DT=(SINGLE SINGLE )'
        DO I=1,NN_ORIG
            WRITE(10,*) X_ORIG(I),Y_ORIG(I)
        END DO
        DO I=1,NE_ORIG
            WRITE(10,*) (NodeEle_ORIG(I,J),J=1,4)
        END DO
        CLOSE(10)
        PAUSE
        stop
    END IF

    END SUBROUTINE EXTEND_MESH_LEFT_0226



    SUBROUTINE EXTEND_MESH_UP_0226
    USE CONSTANT
    USE VARIABLE
    !	USE IMSL
    IMPLICIT NONE

    INTEGER:: I,J,K
    REAL(8):: DS
    REAL(8),ALLOCATABLE:: RA(:),RB(:)
    INTEGER,ALLOCATABLE:: IPERM(:)


    real(8):: xl,xr,yu

    xl=minval(x_orig)
    xr=maxval(x_orig)
    yu=maxval(y_orig)


    ALLOCATE(NODEELE(NUBE_ORIG*uplayer,4))
    !DS=(XR-XL)/(FLOAT(NUBE_ORIG))
    DS = DS_U


    ALLOCATE(RA(NUBN_ORIG),RB(NUBN_ORIG))
    ALLOCATE(IPERM(NUBN_ORIG),I4ARRAY(NUBN_ORIG))

    DO I=1,NUBN_ORIG
        J=UBN_ORIG(I)
        I4ARRAY(I)=J
        IPERM(I)=I
        RA(I)=X_ORIG(J)
    END DO
    CALL DSVRGP (NUBN_ORIG, RA, RB, IPERM)

    DO I=1,NUBN_ORIG
        J=IPERM(I)
        UBN_ORIG(I)=I4ARRAY(J)
    END DO

    DO I=1,uplayer
        DO J=1,NUBE_ORIG
            K=(I-1)*NUBE_ORIG+J
            IF(I==1) THEN
                NODEELE(K,1)=UBN_ORIG(J)
                NODEELE(K,2)=UBN_ORIG(J+1)
            ELSE
                NODEELE(K,1)=NN_ORIG+(I-2)*NUBN_ORIG+J
                NODEELE(K,2)=NN_ORIG+(I-2)*NUBN_ORIG+J+1
            END IF
            NODEELE(K,3)=NN_ORIG+(I-1)*NUBN_ORIG+J+1
            NODEELE(K,4)=NN_ORIG+(I-1)*NUBN_ORIG+J
        END DO
    END DO

    ALLOCATE (x_temp(NN_ORIG+NUBN_ORIG*uplayer),y_temp(NN_ORIG+NUBN_ORIG*uplayer),nodeele_temp(NE_ORIG+NUBE_ORIG*uplayer,4))

    x_temp(1:NN_ORIG)=X_ORIG(1:NN_ORIG)
    y_temp(1:NN_ORIG)=Y_ORIG(1:NN_ORIG)
    nodeele_temp(1:NE_ORIG,:)=NODEELE_ORIG(1:NE_ORIG,:)

    DO I=1,uplayer
        DS=DS*U_pg
        YU=yu+DS
        DO J=1,NUBN_ORIG
            K=NN_ORIG+(I-1)*NUBN_ORIG+J
            x_temp(K)=RB(J)
            y_temp(K)=YU
        END DO
    END DO

    DO I=1,uplayer*NUBE_ORIG
        K=NE_ORIG+I
        nodeele_temp(K,:)=NODEELE(I,:)
    END DO


    NE_ORIG=NE_ORIG+NUBE_ORIG*(float(uplayer))
    NN_ORIG=NN_ORIG+NUBN_ORIG*(float(uplayer))

    DEALLOCATE(RA,RB,IPERM,NODEELE_ORIG,X_ORIG,Y_ORIG,I4ARRAY)

    ALLOCATE(NODEELE_ORIG(NE_ORIG,4),X_ORIG(NN_ORIG),Y_ORIG(NN_ORIG))

    NODEELE_ORIG=nodeele_temp
    X_ORIG=x_temp
    Y_ORIG=y_temp
    DEALLOCATE(nodeele_temp,x_temp,y_temp,NODEELE)


    I=0
    IF (I==1) THEN
        OPEN(UNIT=10,FILE='MESHEXTUP.DAT')
        WRITE(10,*) ' TITLE=" MESHEXTUP " '
        WRITE(10,*) ' VARIABLES="x","y"'
        WRITE(10,*) ' ZONE N=', NN_ORIG,', E=',NE_ORIG,', F=FEPOINT, ET=QUADRILATERAL'
        WRITE(10,*) ' DT=(SINGLE SINGLE )'
        DO I=1,NN_ORIG
            WRITE(10,*) X_ORIG(I),Y_ORIG(I)
        END DO
        DO I=1,NE_ORIG
            WRITE(10,*) (NodeEle_ORIG(I,J),J=1,4)
        END DO
        CLOSE(10)
    END IF


    END SUBROUTINE EXTEND_MESH_UP_0226




    SUBROUTINE EXTEND_MESH_DOWN_0226
    USE CONSTANT
    USE VARIABLE
    !	USE IMSL
    IMPLICIT NONE

    INTEGER:: I,J,K,NET(4)
    REAL(8):: DS
    REAL(8),ALLOCATABLE:: RA(:),RB(:)
    INTEGER,ALLOCATABLE:: IPERM(:)

    real(8):: xl,xr,yd

    xl=minval(x_orig)
    xr=maxval(x_orig)
    yd=minval(y_orig)


    ALLOCATE(NODEELE(NDBE_ORIG*downlayer,4))

    !DS=(XR-XL)/(FLOAT(NdBE_ORIG))
    !DS=-DS
    DS = - DS_D
    
    ALLOCATE(RA(NdBN_ORIG),RB(NdBN_ORIG))
    ALLOCATE(IPERM(NdBN_ORIG),I4ARRAY(NdBN_ORIG))

    DO I=1,NdBN_ORIG
        J=dBN_ORIG(I)
        I4ARRAY(I)=J
        IPERM(I)=I
        RA(I)=X_ORIG(J)
    END DO
    CALL DSVRGP (NdBN_ORIG, RA, RB, IPERM)


    DO I=1,NdBN_ORIG
        J=IPERM(I)
        dBN_ORIG(I)=I4ARRAY(J)
    END DO

    DO I=1,downlayer
        DO J=1,NdBE_ORIG
            K=(I-1)*NdBE_ORIG+J
            IF(I==1) THEN
                NODEELE(K,1)=dBN_ORIG(J)
                NODEELE(K,2)=dBN_ORIG(J+1)
            ELSE
                NODEELE(K,1)=NN_ORIG+(I-2)*ndBN_ORIG+J
                NODEELE(K,2)=NN_ORIG+(I-2)*ndBN_ORIG+J+1
            END IF
            NODEELE(K,3)=NN_ORIG+(I-1)*NdBN_ORIG+J+1
            NODEELE(K,4)=NN_ORIG+(I-1)*NdBN_ORIG+J
        END DO
    END DO

    DO I=1,NDBE_ORIG*downlayer
        NET(:)=NODEELE(I,:)
        NODEELE(I,1)=NET(4)
        NODEELE(I,2)=NET(3)
        NODEELE(I,3)=NET(2)
        NODEELE(I,4)=NET(1)
    END DO


    ALLOCATE (x_temp(NN_ORIG+NdBN_ORIG*downlayer),y_temp(NN_ORIG+NdBN_ORIG*downlayer),nodeele_temp(NE_ORIG+NdBE_ORIG*downlayer,4))

    x_temp(1:NN_ORIG)=X_ORIG(1:NN_ORIG)
    y_temp(1:NN_ORIG)=Y_ORIG(1:NN_ORIG)
    nodeele_temp(1:NE_ORIG,:)=NODEELE_ORIG(1:NE_ORIG,:)

    DO I=1,downlayer
        DS=DS*D_pg
        Yd=yd+DS
        DO J=1,NdBN_ORIG
            K=NN_ORIG+(I-1)*NdBN_ORIG+J
            x_temp(K)=RB(J)
            y_temp(K)=YD
        END DO
    END DO

    DO I=1,downlayer*NdBE_ORIG
        K=NE_ORIG+I
        nodeele_temp(K,:)=NODEELE(I,:)
    END DO


    NE_ORIG=NE_ORIG+NdBE_ORIG *(float(downlayer))
    NN_ORIG=NN_ORIG+NdBN_ORIG*(float(downlayer))

    DEALLOCATE(RA,RB,IPERM,NODEELE_ORIG,X_ORIG,Y_ORIG,I4ARRAY)

    ALLOCATE(NODEELE_ORIG(NE_ORIG,4),X_ORIG(NN_ORIG),Y_ORIG(NN_ORIG))

    NODEELE_ORIG=nodeele_temp
    X_ORIG=x_temp
    Y_ORIG=y_temp
    DEALLOCATE(nodeele_temp,x_temp,y_temp,NODEELE)


    I=1
    IF (I==1) THEN
        OPEN(UNIT=10,FILE='MESH.DAT')
        WRITE(10,*) ' TITLE=" MESHEXTWOWN " '
        WRITE(10,*) ' VARIABLES="x","y"'
        WRITE(10,*) ' ZONE N=', NN_ORIG,', E=',NE_ORIG,', F=FEPOINT, ET=QUADRILATERAL'
        WRITE(10,*) ' DT=(SINGLE SINGLE )'
        DO I=1,NN_ORIG
            WRITE(10,*) X_ORIG(I),Y_ORIG(I)
        END DO
        DO I=1,NE_ORIG
            WRITE(10,*) (NodeEle_ORIG(I,J),J=1,4)
        END DO
        CLOSE(10)
    END IF

    DEALLOCATE(LBE_ORIG,RBE_ORIG,DBE_ORIG, UBE_ORIG,PBE_ORIG,SBE_ORIG)
    DEALLOCATE(LBN_ORIG,RBN_ORIG,DBN_ORIG, UBN_ORIG,PBN_ORIG,SBN_ORIG)
    DEALLOCATE(be_orig,ELE_ADJ,NODE_ELE)

    END SUBROUTINE EXTEND_MESH_DOWN_0226



    subroutine ElementType          !目的：单元是否是四边形
    use constant
    use variable

    integer:: i,j,k,net(4)
    allocate(ipe(ne_orig))

    do i=1,ne_orig
        net(:)=nodeele_orig(i,:)
        j=net(3)
        k=net(4)
        if(j/=k) then
            ipe(i)=4
        else
            ipe(i)=3
        end if
    end do
    end subroutine ElementType





    subroutine centre_ele           !目的：求该单元中心坐标
    use constant
    use variable
    implicit none

    integer:: i,j,k,net(4),num
    real(8):: x0,y0
    allocate(cxye(ne_orig,2))

    do i=1,ne_orig
        net(:)=nodeele_orig(i,:)
        x0=0.0d0
        y0=0.0d0

        num=ipe(i)
        do k=1,num
            j=net(k)
            x0=x0+x_orig(j)
            y0=y0+y_orig(j)
        end do

        cxye(i,1)=x0/float(num)
        cxye(i,2)=y0/float(num)
    end do

    end subroutine centre_ele





    SUBROUTINE NODE_ADJ_NODE        !结点周围的节点
    USE CONSTANT
    USE VARIABLE
    IMPLICIT NONE
    INTEGER:: I,J,K,II,JJ,KK,NUM,MAXNUM


    ALLOCATE(I4MATRIX (maxnum_node_ele*4, NN_ORIG) )
    I4MATRIX=0

    MAXNUM=0
    DO I=1,NN_ORIG
        NUM=0
        DO J=1,MAXNUM_NODE_ELE
            JJ=NODE_ELE(I,J)   !单元
            IF(JJ/=0 .AND. JJ<=NE_ORIG) THEN
                DO K=1,ipe(jj) !4
                    NUM=NUM+1
                    I4MATRIX(NUM,I)=NODEELE_ORIG(JJ,K)!结点i周围的结点
                END DO
            END IF
        END DO
        IF(NUM>MAXNUM) THEN
            MAXNUM=NUM
        END IF
    END DO

    DO I=1,NN_ORIG
        do j=1,MAXNUM
            if(I4MATRIX(J,I)==I) then ! Excluding itself
                I4MATRIX(J,I)=0
            end if
        end do

        DO J=1,MAXNUM-1
            JJ=I4MATRIX(J,I)
            IF(JJ/=0) THEN
                DO K=J+1,MAXNUM
                    KK=I4MATRIX(K,I)
                    IF(JJ==KK) THEN ! Excluding repeated nodes
                        I4MATRIX(K,I)=0
                    END IF
                END DO
            END IF
        END DO

    END DO

    NUM=0
    DO I=1,NN_ORIG
        K=0
        DO J=1,MAXNUM
            JJ=I4MATRIX(J,I)
            IF(JJ/=0) THEN
                K=K+1
                I4MATRIX(K,I)=JJ
                IF(J>K) THEN
                    I4MATRIX(J,I)=0
                END IF
            END IF
        END DO

        IF(K>NUM) THEN
            NUM=K
        END IF

    END DO

    MAXNUM_NODE_NODE =NUM

    ALLOCATE(NODE_NODE(NN_ORIG,MAXNUM_NODE_NODE))
    DO I=1,NN_ORIG
        NODE_NODE(I,:)=I4MATRIX(1:MAXNUM_NODE_NODE,I)
    END DO
    DEALLOCATE(I4MATRIX)

    WRITE(*,*) ' MAXNUM_NODE_NODE=',MAXNUM_NODE_NODE
    WRITE(*,*)


    END SUBROUTINE NODE_ADJ_NODE



    SUBROUTINE SURROUND_ELE         !目的：找单元周围的单元
    USE CONSTANT
    USE VARIABLE
    IMPLICIT NONE

    INTEGER:: I,J,K,NODE,NUM,nmax,ip
    INTEGER:: ELE_ELE_TEMP(30),M,N

    if (allocated(array_large)) deallocate(array_large)
    ALLOCATE(ARRAY_LARGE(NN_orig))

    ARRAY_LARGE=0

    DO I=1,NE_orig
        ip=ipe(i)
        DO J=1,ip
            NODE=NODEELE_orig(I,J)
            ARRAY_LARGE(NODE)=ARRAY_LARGE(NODE)+1
        END DO

    END DO

    maxnum_node_ele=maxval(array_large)      !即maxnum_node_ele=4     但是运行出来=5      有重合点????但是前边已删除重复的结点
    !???????????????????????????????????????????????????????????????

    if (allocated(NODE_ELE)) deallocate(NODE_ELE)
    ALLOCATE(NODE_ELE(NN_orig,maxnum_node_ele))

    NODE_ELE=0

    ARRAY_LARGE=0

    DO I=1,NE_orig
        ip=ipe(i)
        DO J=1,ip
            NODE=NODEELE_orig(I,J)
            ARRAY_LARGE(NODE)=ARRAY_LARGE(NODE)+1
            K=ARRAY_LARGE(NODE)
            NODE_ELE(NODE,K)=I      !---------------XHL----------------第k个结点node属于单元i
        END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    nmax=0

    if (allocated(i4matrix)) deallocate(i4matrix)
    ALLOCATE(i4matrix(NE_orig,30))

    i4matrix=0

    DO I=1,NE_orig
        NUM=0
        ELE_ELE_TEMP=0
        ip=ipe(i)
        DO J=1,ip
            NODE=NODEELE_orig(I,J)
            DO K=1,maxnum_node_ele
                M=NODE_ELE(NODE,K)
                IF (M.NE.0) THEN
                    NUM=NUM+1
                    ELE_ELE_TEMP(NUM)=M
                ELSE
                    EXIT
                END IF
            END DO
        END DO

        DO J=1,30
            IF (ELE_ELE_TEMP(J).NE.0) THEN
                IF (ELE_ELE_TEMP(J)==I) THEN
                    ELE_ELE_TEMP(J)=0
                END IF
            ELSE
                EXIT
            END IF
        END DO

        DO J=1,29
            M=ELE_ELE_TEMP(J)
            IF (M.NE.0) THEN
                DO K=J+1,30
                    N=ELE_ELE_TEMP(K)
                    IF (M==N) THEN
                        ELE_ELE_TEMP(K)=0
                    END IF
                END DO
            END IF
        END DO

        N=0
        DO J=1,30
            M=ELE_ELE_TEMP(J)
            IF (M.NE.0) THEN
                N=N+1
                IF (N > 30) THEN
                    WRITE(*,*) 'MORE IN SUBROUTNE SURROUND ELE'
                    PAUSE
                    STOP
                END IF
                i4matrix(I,N)=M              !放置单元周围的单元
            END IF
        END DO

        if(n>nmax) then
            nmax=n                            !----------------XHL------------------单元i周围有n个单元
        end if

    END DO


    MAXNUM_ELE_ADJ=NMAX                    !------------XHL--------------MAXNUM_ELE_ADJ=N

    if (allocated(ELE_ADJ)) deallocate(ELE_ADJ)
    ALLOCATE(ELE_ADJ(NE_orig,MAXNUM_ELE_ADJ))

    ELE_ADJ(:,1:MAXNUM_ELE_ADJ)=I4MATRIX(:,1:MAXNUM_ELE_ADJ)     !----------------XHL------------------单元i周围的单元

    DEALLOCATE(ARRAY_LARGE,I4MATRIX)

    END SUBROUTINE SURROUND_ELE
