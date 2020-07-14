float calculateLocalD(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, YE, 2, 1, m));
    row1.push_back(calcularTenedor(e, ZETA, 2, 1, m));

    row2.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, ZETA, 3, 1, m));

    row3.push_back(calcularTenedor(e, EQUIS, 4, 1, m));
    row3.push_back(calcularTenedor(e, YE, 4, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

void calculateAlpha(int i,Matrix &A,mesh m){
    zeroes(A,3);
    element e = m.getElement(i);
    
    A.at(0).at(0) = OperarRestaTenedor(e, YE, ZETA, 3, 4, m);
    A.at(0).at(1) = OperarRestaTenedor(e, YE, ZETA, 4, 2, m);
    A.at(0).at(2) = OperarRestaTenedor(e, YE, ZETA, 2, 3, m);

    A.at(1).at(0) = OperarRestaTenedor(e, EQUIS, ZETA, 4, 3, m);
    A.at(1).at(1) = OperarRestaTenedor(e, EQUIS, ZETA, 2, 4, m);
    A.at(1).at(2) = OperarRestaTenedor(e, EQUIS, ZETA, 3, 2, m);

    A.at(2).at(0) = OperarRestaTenedor(e, EQUIS, YE, 3, 4, m);
    A.at(2).at(1) = OperarRestaTenedor(e, EQUIS, YE, 4, 2, m);
    A.at(2).at(2) = OperarRestaTenedor(e, EQUIS, YE, 2, 3, m);

}

void calculateBeta(Matrix &B){
    zeroes(B,3,12);

    B.at(0).at(0) = -1; 
    B.at(0).at(1) =  1; 
    B.at(0).at(2) =  0; 
    B.at(0).at(3) =  0; 
    B.at(0).at(4) = -1; 
    B.at(0).at(5) =  1; 
    B.at(0).at(6) =  0; 
    B.at(0).at(7) =  0; 
    B.at(0).at(8) = -1; 
    B.at(0).at(9) =  1; 
    B.at(0).at(10) =  0; 
    B.at(0).at(11) =  0; 
    
    B.at(1).at(0) = -1; 
    B.at(1).at(1) = 0; 
    B.at(1).at(2) = 1;
    B.at(1).at(3) = 0;
    B.at(1).at(4) = -1; 
    B.at(1).at(5) = 0; 
    B.at(1).at(6) = 1;
    B.at(1).at(7) = 0;
    B.at(1).at(8) = -1; 
    B.at(1).at(9) = 0; 
    B.at(1).at(10) = 1;
    B.at(1).at(11) = 0;

    B.at(2).at(0) = -1;
    B.at(2).at(1) = 0;
    B.at(2).at(2) = 0;
    B.at(2).at(3) = 1;
    B.at(2).at(4) = -1;
    B.at(2).at(5) = 0;
    B.at(2).at(6) = 0;
    B.at(2).at(7) = 1;
    B.at(2).at(8) = -1;
    B.at(2).at(9) = 0;
    B.at(2).at(10) = 0;
    B.at(2).at(11) = 1;
}

void calculateOmega(Matrix &C){
    zeroes(C,3,4);
    C.at(0).at(0) = -1; C.at(0).at(1) = 1; C.at(0).at(2) = 0; C.at(0).at(3) = 0;
    C.at(1).at(0) = -1; C.at(1).at(1) = 0; C.at(1).at(2) = 1; C.at(1).at(3) = 0;
    C.at(2).at(0) = -1; C.at(2).at(1) = 0; C.at(2).at(2) = 0; C.at(2).at(3) = 1;
}


void calculateGamma(Matrix &m){
	zeroes(m,12,3);

	m.at(0).at(0) = 1;   m.at(0).at(1) = 0;   m.at(0).at(2) = 0;
	m.at(1).at(0) = 1;   m.at(1).at(1) = 0;   m.at(1).at(2) = 0; 
    m.at(2).at(0) = 1;   m.at(2).at(1) = 0;   m.at(2).at(2) = 0;
	m.at(3).at(0) = 1;   m.at(3).at(1) = 0;   m.at(3).at(2) = 0; 
    m.at(4).at(0) = 0;   m.at(4).at(1) = 1;   m.at(4).at(2) = 0;
	m.at(5).at(0) = 0;   m.at(5).at(1) = 1;   m.at(5).at(2) = 0; 
    m.at(6).at(0) = 0;   m.at(6).at(1) = 1;   m.at(6).at(2) = 0;
	m.at(7).at(0) = 0;   m.at(7).at(1) = 1;   m.at(7).at(2) = 0; 
    m.at(8).at(0) = 0;   m.at(8).at(1) = 0;   m.at(8).at(2) = 1;
	m.at(9).at(0) = 0;   m.at(9).at(1) = 0;   m.at(9).at(2) = 1; 
    m.at(10).at(0) = 0;  m.at(10).at(1) = 0;  m.at(10).at(2) = 1;
	m.at(11).at(0) = 0;  m.at(11).at(1) = 0;  m.at(11).at(2) = 1; 
	
}

float calculateLocalJ(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 4, 1, m));

    row2.push_back(calcularTenedor(e, YE, 2, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 4, 1, m));

    row3.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

Matrix createLocalM(int e,mesh &m){
     Matrix matrixG,matrixC,matrixF;
    float J,Determinant;
    
    /* [  0  -G ]
       [  C  -F ]
    */

    //Matrix G

    Matrix Alpha, Beta;
    Matrix Alpha_t,Beta_t;

    calculateAlpha(e,Alpha,m);
    calculateBeta(Beta);

    transpose(Alpha,Alpha_t);
    transpose(Beta,Beta_t);

    Determinant = calculateLocalD(e,m);
    J = calculateLocalJ(e,m);

    if(Determinant == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }

    element elem = m.getElement(e);

    node n1 = m.getNode(elem.getNode1()-1);
    node n2 = m.getNode(elem.getNode2()-1);
    node n3 = m.getNode(elem.getNode3()-1);
    node n4 = m.getNode(elem.getNode4()-1);
    
    float real_a = (float) (-3*(n1.getX()+n2.getX()+n3.getX()+n4.getX())*J)/(24*pow(Determinant,2));
    
    productRealMatrix(real_a,productMatrixMatrix(Beta_t,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,4),3,3,4),4,3,4),matrixG);


    //Matrix C

    
    Matrix matrixH;
    zeroes(matrixH,4,3);

    matrixH.at(0).at(0) = 2*n1.getX()+n2.getX()+n3.getX()+n4.getX();
    matrixH.at(0).at(1) = n1.getX()+2*n2.getX()+n3.getX()+n4.getX();
    matrixH.at(0).at(2) = n1.getX()+n2.getX()+2*n3.getX()+n4.getX();
    matrixH.at(0).at(0) = n1.getX()+n2.getX()+n3.getX()+2*n4.getX();

    matrixH.at(1).at(0) = 2*n1.getX()+n2.getX()+n3.getX()+n4.getX();
    matrixH.at(1).at(1) = n1.getX()+2*n2.getX()+n3.getX()+n4.getX();
    matrixH.at(1).at(2) = n1.getX()+n2.getX()+2*n3.getX()+n4.getX();
    matrixH.at(1).at(0) = n1.getX()+n2.getX()+n3.getX()+2*n4.getX();

    matrixH.at(2).at(0) = 2*n1.getX()+n2.getX()+n3.getX()+n4.getX();
    matrixH.at(2).at(1) = n1.getX()+2*n2.getX()+n3.getX()+n4.getX();
    matrixH.at(2).at(2) = n1.getX()+n2.getX()+2*n3.getX()+n4.getX();
    matrixH.at(2).at(0) = n1.getX()+n2.getX()+n3.getX()+2*n4.getX();

    float real_k = (float) J/(120*Determinant);

    productRealMatrix(real_k,productMatrixMatrix(matrixH,productMatrixMatrix(Alpha,Beta,3,3,4),4,3,4), matrixC);

    
    //Matrix F
    Matrix matrixI;
    zeroes(matrixI,4,3);

    matrixI.at(0).at(0) = -(2*n1.getY()+n3.getY()+n4.getY());
    matrixI.at(0).at(1) = 2*n2.getY()+n3.getY()+n4.getY();
    matrixI.at(0).at(2) = -(n1.getY()-n2.getY());
    matrixI.at(0).at(0) = -(n1.getY()-n2.getY());

    matrixI.at(1).at(0) = -(2*n1.getY()+n2.getY()+n4.getY());
    matrixI.at(1).at(1) = -(n1.getY()-n3.getY());
    matrixI.at(1).at(2) = n2.getY()+2*n3.getY()+n4.getY();
    matrixI.at(1).at(0) = -(n1.getY()-n3.getY());

    matrixI.at(2).at(0) = -(2*n1.getY()+n2.getY()+n3.getY());
    matrixI.at(2).at(1) = -(n1.getY()-n4.getY());
    matrixI.at(2).at(2) = -(n1.getY()-n4.getY());
    matrixI.at(2).at(0) = n2.getY()+n3.getY()+2*n4.getY();

    float real_g = (float) (J/(24*Determinant*Determinant));

    productRealMatrix(real_g,productMatrixMatrix(matrixI,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,4),3,3,4),4,3,4),matrixF);

    //Matrix M
     Matrix M, fill;
    zeroes(M,8);
    zeroes(fill,4);
    for (int i = 0; i < 4; i++){
        for (int e = 0; e < 4; e++){ 
            fill.at(i).at(e) += static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        }
        
    }
    ubicarSubMatriz(M,0,3,0,3,fill);
    ubicarSubMatriz(M,0,3,4,7,matrixG);
    ubicarSubMatriz(M,4,7,0,3,matrixC);
    ubicarSubMatriz(M,4,7,4,7,matrixF);
    
    return M;

    return M;
}

Vector createLocalb(int e,mesh &m){
    float J;
    float tempx = EXTERNAL_TEMP_X, tempy = EXTERNAL_TEMP_Y,tempz = EXTERNAL_TEMP_Z;
    Vector b,vectorB, vectorE;

    J = calculateLocalJ(e,m);

    if(J == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }

    float real_h = ((tempx+tempy+tempz)*J)/24;

    Vector unitario;
    zeroes(unitario,4);

    unitario.at(0) += 1;
    unitario.at(1) += 1;
    unitario.at(2) += 1;
    unitario.at(3) += 1;

    productRealVector(real_h,unitario,vectorB);

    float real_i = (-273.15*J)/24;

    productRealVector(real_i,unitario,vectorE);

    zeroes(b,8);
    b.at(0) = vectorB.at(0);
    b.at(1) = vectorB.at(1);
    b.at(2) = vectorB.at(2);
    b.at(3) = vectorB.at(3);
    b.at(4) = vectorE.at(0);
    b.at(5) = vectorE.at(1);
    b.at(6) = vectorE.at(2);
    b.at(7) = vectorE.at(3);
    
    
    return b;
}
