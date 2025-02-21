    Module Local
    Implicit None

    Integer::NSPBN_ALE_KS
    Real(8)::U_inf = 1.0D0
    Real(8),Allocatable::SPBN_ALE_KS(:),Sec_Rota(:,:),Sec_Rota0(:,:),Delta_x510(:),delta_y510(:),Delta_x0510(:),delta_y0510(:)
    Integer,Allocatable::PlatePointIndex_UDF(:,:),PlatePointIndex_Tail(:,:),BoundType(:)
    Real(8)::SupportR = 6.0D0
    Real(8),ALLOCATABLE::RBF_POS_X_Delta_Into_Red(:)
    Real(8),ALLOCATABLE::RBF_POS_Y_Delta_Into_Red(:)
    Integer::Reduction_Size,Mesh_Update_Num
    Real(8),Allocatable::Err_X_Red(:,:),Err_Y_Red(:,:)
    Integer,Allocatable::BN_ALE_KS_Red(:),Mesh_Update_Index(:)
    Real(8),ALLOCATABLE::RBF_PhiMA(:,:)!For Red 
    
    Real(8),ALLOCATABLE::RBF_PhiMA_Constant(:,:),Left_Constant(:,:)
    
    End Module


    MODULE ANCF
    Implicit None
    TYPE ANCF_PARAMS
        Integer::NumD       = 6
        Integer::DIM        = 3
        Integer::NumE       = 12!elements of a ele
        REAL(8)::L
        INTEGER::NE             !Num of cal ele
        INTEGER::NNODE
        REAL(8)::LE
        REAL(8)::E
        REAL(8)::A
        REAL(8)::I
        REAL(8)::RHO
        Integer::fixed(2)   = [1,0]
    END TYPE ANCF_PARAMS

    TYPE(ANCF_PARAMS)::AFP

    REAL(8),ALLOCATABLE::AF_MASS(:,:),AF_Fex(:),AF_FeL(:)   !Fex : Force External,Fel: Force Elastic

    REAL(8),ALLOCATABLE::AF_e(:),AF_ed(:),AF_edd(:)
    REAL(8),ALLOCATABLE::AF_e0(:),AF_ed0(:),AF_edd0(:)

    Real(8)            ::ForceXFlexi,ForceYFlexi,MomentFlexi !Cylinder Force
    Real(8)            ::AFFex2,AFFel2
    Real(8),Allocatable::AFMass2(:,:)
    Real(8),Allocatable::PlatePointIndex(:,:)
    REAL(8),ALLocatable::X_Neutal(:,:),Y_Neutal(:,:)
    Integer,Allocatable::SBE2SBN_Index(:)
    Real(8),Allocatable::NumEMeshPoint(:,:)

    Integer::NGQ1,NGQ2
    Real(8),Allocatable::GQ3xi(:),GQ3wei(:),GQ5xi(:),GQ5wei(:)

    Contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SubRoutine phifun(x1,y1,x0,y0,fks_phi)
    Use Local
    Real(8)::x1,y1,x0,y0,fks_phi
    Real(8)::EulerDis,fks_xi

    fks_xi = Dsqrt((x1-x0)**2+(y1-y0)**2) / SupportR

    If (fks_xi<=1.0) Then
        fks_phi = (1-fks_xi)**4.0*(4*fks_xi+1)
    Else
        fks_phi = 0.0
    End if

    End SubRoutine phifun
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SubRoutine phifun_V2(x1,y1,x0,y0,fks_phi,SR)
    Real(8)::x1,y1,x0,y0,fks_phi,SR
    Real(8)::EulerDis,SupportR,fks_xi

    SupportR = SR
    EulerDis = Dsqrt((x1-x0)**2+(y1-y0)**2)
    fks_xi = EulerDis / SupportR

    If (fks_xi<=1.0) Then
        fks_phi = (1-fks_xi)**4.0*(4*fks_xi+1)
    Else
        fks_phi = 0.0
    End if

    End SubRoutine phifun_V2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    END MODULE

