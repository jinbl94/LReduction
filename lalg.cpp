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
        for(index = 0; index < m; index++){
            if(index == row){ continue; }
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
    for(long col = 0; col < m; col++){
        if(col == row){ continue; }
        abs(tr, vecrow[col]);
        if(maxrow < tr){ maxrow = tr; }
    }
}

void FactorRange(const long& m, const mat_RR& MU, const long& row, const long& currentrow,
        const RR& mu, long& maxfactor, long& minfactor)
{
    RR tr, tr1;
    maxfactor = 9999;
    minfactor = -9999;
    long col;
    for(col = 0; col < m; col++){
        if(col == row){ continue; }
        if(MU[currentrow][col] == 0){ }else{
            sub(tr, mu, MU[row][col]);
            div(tr, tr, MU[currentrow][col]);
            round(tr, tr);

            add(tr1, mu, MU[row][col]);
            negate(tr1, tr1);
            div(tr1, tr1, MU[currentrow][col]);
            round(tr1, tr1);
            if(tr > tr1){
                if(maxfactor > tr){ conv(maxfactor, tr); }
                if(minfactor < tr1){ conv(minfactor, tr1); }
                // std::cout << "range(" << tr1 << ", " << tr << ")\n";
            }else{
                if(maxfactor > tr1){ conv(maxfactor, tr1); }
                if(minfactor < tr){ conv(minfactor, tr); }
                // std::cout << "range(" << tr << ", " << tr1 << ")\n";
            }
        }
    }
    // std::cout << row << ", " << currentrow << ": min: " << minfactor << ", max: " << maxfactor << std::endl;
}

void BestFactor(const long& m, const mat_RR& MU, const long& row, const long& currentrow, RR& oldmax,
        const long& minfactor, const long& maxfactor, long& bestindex, long& bestfactor)
{
    RR tr, tr1;
    for(long i = minfactor; i <= maxfactor; i++){
        if(i == 0){ continue; }
        tr = 0;
        for(long col = 0; col < m; col++){
            if(col == row){ continue; }
            mul(tr1, i, MU[currentrow][col]);
            add(tr1, tr1, MU[row][col]);
            abs(tr1, tr1);
            if(tr < tr1){ tr = tr1; }
        }
        if(tr < oldmax){
            oldmax = tr;
            bestindex = currentrow;
            bestfactor = i;
        }
    }
}

void RowReduce(const long m, const mat_RR& MU, const RR& mu, const long& row, long& bestindex, long& bestfactor)
{
    RR maxrow;
    long minfactor, maxfactor;
    // find max absolute value of current row
    bestindex = 0;
    bestfactor = 0;
    maxrow = mu;
    for(long currentrow = 0; currentrow < m; currentrow++){
        if(currentrow == row){ continue; }
        FactorRange(m, MU, row, currentrow, mu, maxfactor, minfactor);
        if(minfactor > maxfactor){ continue; }
        BestFactor(m, MU, row, currentrow, maxrow, minfactor, maxfactor, bestindex, bestfactor);
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

    long bestindex, bestfactor, loop;
    RR mu, muold, detmustart, detmuend, optrate;

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
    determinant(detmustart, MU);

    // Reduce while mu < muold
    loop = 0;
    bool s = true;
    while(s){
        s = false;
        Rearrange(m, P, Prod);
        MUMax(m, MU, mu);
        for(long i = 0; i < m; i++){
            RowReduce(m, MU, mu, P[i], bestindex, bestfactor);
            Alter(B, m, Prod, MU, P[i], bestindex, bestfactor);
            if(bestfactor != 0){ s = true; }
        }
        loop++;
    }

    determinant(detmuend, MU);
    div(optrate, detmuend, detmustart);
    float result;
    conv(result, optrate);

    return result;
}
