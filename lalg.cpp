#include "lalg.h"

void Init(const mat_ZZ& B, const long m, vec_long& P, mat_ZZ& Prod, mat_RR & MU) // init all parameters
{
    RR tr1, tr2;

    for(long i = 0; i < m; i++){
        P[i] = i;
        MU[i][i] = 1;
    }

    // Prod[i][i] = B[i]^2
    for(long i = 0; i < m; i++){
        InnerProduct(Prod[i][i], B[i], B[i]);
    }

    for(long i = 1; i < m; i++){
        for(long j = 0; j < i; j++){
            InnerProduct(Prod[i][j], B[i], B[j]);
            Prod[j][i] = Prod[i][j];

            // MU[i][j] = Prod[i][j]/Prod[j][j]
            conv(tr1, Prod[i][j]);
            conv(tr2, Prod[j][j]);
            div(tr1, tr1, tr2);
            MU[i][j] = tr1;
            abs(tr1, tr1);

            // MU[j][i] = tr1 = Prod[i][j]/Prod[i][i]
            conv(tr1, Prod[j][i]);
            conv(tr2, Prod[i][i]);
            div(tr1, tr1, tr2);
            MU[j][i] = tr1;
            abs(tr1, tr1);
        }
    }
}

void Rearrange(const long& m, vec_long& P, const mat_ZZ& Prod)
{
    ZZ tz;
    long tl, i, j;
    for(i = 1; i < m; i++){
        tz = Prod[P[i]][P[i]];
        tl = P[i];
        for(j = i - 1; j >= 0; j--){
            if(tz < Prod[P[j]][P[j]]){
                P[j + 1] = P[j];
            }else{ break; }
        }
        P[j + 1] = tl;
    }
}

void Alter(mat_ZZ& B, const long& m, mat_ZZ& Prod, mat_RR& MU, const long& row, const long& j, const long& q)
{
    if(q){
        ZZ tz;
        RR tr, tr1;
        long index;

        // B[row] = B[row] + q * B[j]
        for(index = 0; index < m; index++){
            mul(tz, q, B[j][index]);
            add(B[row][index], B[row][index], tz);
        }

        // Prod[row][row] = Prod[row][row]^2 + q^2 Prod[j][j] + 2q Prod[row][j]
        mul(tz, q * q, Prod[j][j]);
        add(Prod[row][row], Prod[row][row], tz);
        mul(tz, 2 * q, Prod[row][j]);
        add(Prod[row][row], Prod[row][row], tz);

        // Prod[r][row] = Prod[row][r] = Prod[row][r] + q * Prod[j][r]
        for(index = 0; index < row; index++){
            mul(tz, q, Prod[j][index]);
            add(Prod[row][index], Prod[row][index], tz);
            Prod[index][row] = Prod[row][index];
        }
        for(index = row + 1; index < m; index++){
            mul(tz, q, Prod[j][index]);
            add(Prod[row][index], Prod[row][index], tz);
            Prod[index][row] = Prod[row][index];
        }

        for(index = 0; index < m; index++){
            // MU[row][r] = MU[row][r] + q MU[j][r]
            mul(tr, q, MU[j][index]);
            add(MU[row][index], MU[row][index], tr);
            abs(tr, MU[row][index]);

            // MU[r][row] = Prod[r][row]/Prod[row][row]
            conv(tr, Prod[index][row]);
            conv(tr1, Prod[row][row]);
            div(MU[index][row], tr, tr1);
        }    
    }
}

void MaxRow(const long m, const vec_RR& vecrow, const long& row, RR& maxrow)
{
    RR tr;
    maxrow = 0;

    for(long col = 0; col < row; col++){
        abs(tr, vecrow[col]);
        if(maxrow < tr){ maxrow = tr; }
    }
    for(long col = row + 1; col < m; col++){
        abs(tr, vecrow[col]);
        if(maxrow < tr){ maxrow = tr; }
    }
}

void FactorRange(const long& m, const mat_RR& MU, const long& row, const long& currentrow,
        const RR& maxrow, long& maxfactor, long& minfactor)
{
    RR tr, a, b;
    long col;
    // calculate the range of factor
    for(col = 0; col < row; col++){
        a = MU[row][col];
        b = MU[currentrow][col];
        if(b == 0){
            // do nothing
        }else if(b > 0){
            sub(tr, maxrow, a);
            div(tr, tr, b);
            round(tr, tr);
            if(maxfactor > tr){ conv(maxfactor, tr); }
            negate(tr, maxrow);
            sub(tr, tr, a);
            div(tr, tr, b);
            round(tr, tr);
            if(minfactor < tr){ conv(minfactor, tr); }
        }else if(b < 0){
            sub(tr, maxrow, a);
            div(tr, tr, b);
            round(tr, tr);
            if(minfactor < tr){ conv(minfactor, tr); }

            negate(tr, maxrow);
            sub(tr, tr, a);
            div(tr, tr, b);
            round(tr, tr);
            if(maxfactor > tr){ conv(maxfactor, tr); }
        }
    }
    for(col = row + 1; col < m; col++){
        a = MU[row][col];
        b = MU[currentrow][col];
        if(b == 0){
            // do nothing
        }else if(b > 0){
            sub(tr, maxrow, a);
            div(tr, tr, b);
            round(tr, tr);
            if(maxfactor > tr){ conv(maxfactor, tr); }

            negate(tr, maxrow);
            sub(tr, tr, a);
            div(tr, tr, b);
            round(tr, tr);
            if(minfactor < tr){ conv(minfactor, tr); }
        }else if(b < 0){
            sub(tr, maxrow, a);
            div(tr, tr, b);
            round(tr, tr);
            if(minfactor < tr){ conv(minfactor, tr); }

            negate(tr, maxrow);
            sub(tr, tr, a);
            div(tr, tr, b);
            round(tr, tr);
            if(maxfactor > tr){ conv(maxfactor, tr); }
        }
    }
}

void BestFactor(const long& m, const mat_RR& MU, const long& row, const long& currentrow, RR& oldmax,
        const long& minfactor, const long& maxfactor, long& j, long& q)
{
    RR newmax, tr;
    long col;
    for(long factor = minfactor; factor <= maxfactor; factor++){
        newmax = 0;

        for(col = 0; col < row; col++){
            mul(tr, factor, MU[currentrow][col]);
            add(tr, tr, MU[row][col]);
            abs(tr, tr);
            if(newmax < tr){ newmax = tr; }
        }
        for(col = row + 1; col < m; col++){
            mul(tr, factor, MU[currentrow][col]);
            add(tr, tr, MU[row][col]);
            abs(tr, tr);
            if(newmax < tr){ newmax = tr; }
        }

        if(newmax < oldmax){
            oldmax = newmax;
            j = currentrow;
            q = factor;
        }
    }
}

void RowReduce(const long m, const mat_RR& MU, const long& row, long& j, long& q)
{
    RR maxrow, oldmax;
    long currentrow, minfactor, maxfactor;
    minfactor = -999999;
    maxfactor = 999999;

    // find max absolute value of current row
    MaxRow(m, MU[row], row, maxrow);
    oldmax = maxrow;

    for(currentrow = 0; currentrow < row; currentrow++){
        FactorRange(m, MU, row, currentrow, maxrow, maxfactor, minfactor);

        if(minfactor < maxfactor){
            BestFactor(m, MU, row, currentrow, oldmax, minfactor, maxfactor, j, q);
        }
    }
}

void MUMax(const long& m, const mat_RR& MU, RR& mu)
{
    mu = 0;
    RR tr;
    for(long i = 0; i < m; i++){
        for (long j = 0; j < i; j++){
            abs(tr, MU[i][j]);
            if(mu < tr){ mu = tr; }
            abs(tr, MU[j][i]);
            if(mu < tr){ mu = tr; }
        }
    }
}

static float LReduction(mat_ZZ& B)
{
    assert( B.NumRows() == B.NumCols());
    const long m = B.NumRows();

    long row, j, q, loop;
    float result;
    RR mu, mu_old, detmu_start, detmu_end, optrate;

    // permutation of B
    vec_long P;
    P.SetLength(m);
    // innerproduct of vectors
    mat_ZZ Prod;
    Prod.SetDims(m, m);
    // \mu matrix
    mat_RR MU;
    MU.SetDims(m, m);

    Init(B, m, P, Prod, MU);
    MUMax(m, MU, mu);
    determinant(detmu_start, MU);

    std::cout
        << "basis:" << std::endl << B << std::endl
        << "innerproduct matrix:" << std::endl << Prod << std::endl
        << "mu max: " << mu << std::endl
        << "det(mu): " << detmu_start << std::endl << std::endl;

    loop = 0;
    while(loop < 10){
        for(row = 0; row < m; row++){
            Rearrange(m, P, Prod);
            RowReduce(m, MU, P[row], j, q);
            Alter(B, m, Prod, MU, row, j, q);
        }
        loop++;
    }

    MUMax(m, MU, mu);
    determinant(detmu_end, MU);
    div(optrate, detmu_end, detmu_start);
    conv(result, optrate);

    std::cout
        << "basis:" << std::endl << B << std::endl
        << "innerproduct matrix:" << std::endl << Prod << std::endl
        << "mu max: " << mu << std::endl
        << "det(mu): " << detmu_end << std::endl
        << "optimize rate: " << optrate << std::endl;

    return result;
}
