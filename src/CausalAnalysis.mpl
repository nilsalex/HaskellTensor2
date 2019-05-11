#function to compute the symbol of a rank deficient matrix 
CausalAnalysis := module()

export evalRand, randSubMatrix, linPoly, randSubMatrixQuad, evalRandQuad, prodTrace, quadPoly, solveMatrixEqns;

option package;

#evaluate all constants randomly
evalRand := proc(M::Matrix)
    uses LinearAlgebra;
    vars := convert(indets(M) minus {k[0],k[1],k[2],k[3]},list);
    fRand := rand(-1000..1000);
    evalL := zip((a,b) -> a = b, vars, [seq(fRand(), i = 1..nops(vars))]);
    (simplify(subs(evalL,M)), evalL);
    end proc;

#construct a random SubMatrix 
randSubMatrix := proc(M::Matrix)
    uses LinearAlgebra, combinat;
    n := RowDimension(M);
    rowList := [seq(1..n)];
    M2 := SubMatrix(M,randcomb(rowList,n-4),randcomb(rowList,n-4));
    if Rank(M2) = n-4 then M2; else randSubMatrix(M); end if;
    end proc;

#compute the principal polynomial for one matrix (linear order) and one rand combination 
linPoly := proc(M::Matrix)
    uses LinearAlgebra;
    n := RowDimension(M);
    (MRand, evalL) := evalRand(M);
    SubM := randSubMatrix(MRand);
    Pol := factor(Determinant(SubM, method=multivar));
    end proc;

#construct the subMatrices for the linear Matrix and the list of quadratic matrices
#linear subMatrix must have full rank
randSubMatrixQuad := proc(Lin::Matrix, Quad::list)
    uses LinearAlgebra, combinat;
    n := RowDimension(Lin);
    rowList := [seq(1..n)];
    rowComb := randcomb(rowList,n-4);
    colComb := randcomb(rowList,n-4);
    QuadL := map(x -> SubMatrix(x,rowComb,colComb), Quad);
    LinM := SubMatrix(Lin,rowComb,colComb);
    if  Rank(LinM) < n-4 then 
        randSubMatrixQuad(Lin,Quad); 
        else (LinM, QuadL); 
    end if;
    end proc;

#evaluate randomly in the linear constants as this is probably 
evalRandQuad := proc(Lin::Matrix, Quad::list)
    uses LinearAlgebra;
    varsLin := convert(indets(Lin) minus {k[0],k[1],k[2],k[3]},list);
    fRand := rand(-1000..1000);
    evalL := zip((a,b) -> a = b, varsLin, [seq(fRand(), i = 1..nops(varsLin))]);
    QuadL := map(x -> simplify(subs(evalL,x)), Quad);
    LinM := simplify(subs(evalL, Lin)); 
    (LinM, QuadL);
    end proc;

#compute the trace of a mutrix product 
prodTrace := proc(M::Matrix, Q::Matrix)
    uses LinearAlgebra;
    size := min(RowDimension(M),ColumnDimension(Q));
    rowsM := [Row(M,[seq(1..size)])];
    colsQ := [Column(Q,[seq(1..size)])];
    l := zip((x,y) -> Multiply(x,y), rowsM, colsQ);
    factor(add(l));
    end proc;

quadPoly := proc(M::Matrix, Q::list)
    uses LinearAlgebra;
    (randM, randQ) := evalRandQuad(M,Q);
    (randSubM, randSubQ) := randSubMatrixQuad(randM, randQ);
    subMInv := MatrixInverse(randSubM, method = polynom);
    polyL := zip((x,i) -> H__i * prodTrace(subMInv,x), randSubQ, [seq(1..21)]);
    add(polyL);
    end proc;

solveMatrixEqns := proc(M::Matrix)
    uses LinearAlgebra;
    colsM := ColumnDimension(M);
    rowsM := RowDimension(M);
    zeroVec := ZeroVector(rowsM);
    sol := convert(LinearSolve(M,zeroVec), list);
    vars := [seq(x[i],i=1..colsM)];
    evalL := zip((x,y) -> x = y, vars, sol); 
    end proc;

end module;
