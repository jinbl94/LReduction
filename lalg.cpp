#include "lalg.h"

void Enumerate(const mat_RR& MU, const long& m, const vec_long& P, const long& beta, vec_long& indexes, vec_long& coeff)
{
    vec_long ind, indbound, coe, lowbound, upbound;
    vec_RR tvecrr, tvecrr1;
    RR minmax, tminmax;
    conv(minmax, 1);

    coe.SetLength(beta);
    ind.SetLength(beta);
    indbound.SetLength(beta);
    lowbound.SetLength(beta);
    upbound.SetLength(beta);
    tvecrr.SetLength(m - 1);
    tvecrr1.SetLength(m - 1);

    for(long i = 0; i < beta; i++){
        ind[i] = i + 1;
        indbound[i] = m - beta + i;
    }

    // Run over all possible indexes
    ind[beta - 1] -= 1;
    while(ind != indbound){
        long i;
        for(i = beta - 1; i >= 0; i--){
            if(ind[i] != indbound [i]) break;
        }
        ind[i] += 1;
        ++i;
        for(; i < beta; i++){
            ind[i] = ind[i - 1] + 1;
        }
        // todo: Calculate enumeration bound MU[P[ind[i]]]
        /*****
         * todo
         */

        // Enumerate MU[P[ind]]
        coe = lowbound;
        coe[0] -= 1;
        // Run over all possible combinations
        while(coe != upbound){
            for(i = 0; i < beta; i ++){
                if(coe[i] != upbound[i]) break;
            }
            coe[i] += 1;
            --i;
            for(; i >= 0; i--){
                coe[i] = lowbound[i];
            }
            // Calculate the length of current combination
            tvecrr = MU[P[0]];
            for(i = 0; i < beta; i++){
                mul(tvecrr1, coe[i], MU[P[ind[i]]]);
                add(tvecrr, tvecrr, tvecrr1);
            }
            MaxRow(m - 1, tvecrr, tminmax);
            if(tminmax < minmax){
                minmax = tminmax;
                for(i = 0; i < beta; i++){
                    indexes[i] = P[ind[i]];
                }
                coeff = coe;
            }
        }
    }
}

void Init(const mat_ZZ& B, const long& m, vec_long& P, mat_ZZ& Prod, mat_RR & MU) // init all parameters
{
    RR tr1, tr2;

    for(long i = 0; i < m; i++){
        P[i] = i;
        MU[i][i] = 1;
        InnerProduct(Prod[i][i], B[i], B[i]); // Prod[i][i] = B[i]^2
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

void Alter(const long& m, mat_ZZ& B, mat_ZZ& Prod, mat_RR& MU, const long& row, const long& j, const long& q, mat_ZZ& U)
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

        // U[row] = U[row] + qU[j]
        for(index = 0; index < m; index++){
            mul(tz, q, U[j][index]);
            add(U[row][index], U[row][index], tz);
        }
    }
}
void MaxRow(const long m, const vec_RR& vecrow, RR& maxrow)
{
    RR tr;
    maxrow = 0;
    for(long col = 0; col < m; col++){
        abs(tr, vecrow[col]);
        if(maxrow < tr){ maxrow = tr; }
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
        const RR& maxrow, long& maxfactor, long& minfactor)
{
    RR tr, tr1;
    maxfactor = 9999;
    minfactor = -9999;
    long col;
    for(col = 0; col < m; col++){
        if(col == row){ continue; }
        if(MU[currentrow][col] == 0){ }else{
            sub(tr, maxrow, MU[row][col]);
            div(tr, tr, MU[currentrow][col]);
            round(tr, tr);

            add(tr1, maxrow, MU[row][col]);
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

void BestFactor(const long& m, const mat_RR& MU, const long& row, const long& currentrow, RR& maxrow,
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
        if(tr < maxrow){
            maxrow = tr;
            bestindex = currentrow;
            bestfactor = i;
        }
    }
}

void RowReduce(const long& m, const mat_RR& MU, const long& row, long& bestindex, long& bestfactor)
{
    RR maxrow;
    long minfactor, maxfactor;
    bestindex = 0;
    bestfactor = 0;
    // find max absolute value of current row
    // MaxRow(m, MU[row], row, maxrow);
    MUMax(m, MU, maxrow);
    for(long currentrow = 0; currentrow < m; currentrow++){
        if(currentrow == row){ continue; }
        FactorRange(m, MU, row, currentrow, maxrow, maxfactor, minfactor);
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

static float LReduction(mat_ZZ& B, mat_ZZ& U) // B' = U * B
{
    assert( B.NumRows() == B.NumCols()); // make sure B is a square matrix
    const long m = B.NumRows();

    RR detmustart, detmuend, optrate;
    vec_long P; // permutation of B
    mat_ZZ Prod; // innerproduct of vectors
    mat_RR MU; // mu matrix

    // Initialize all parameters
    P.SetLength(m);
    Prod.SetDims(m, m);
    MU.SetDims(m, m);
    Init(B, m, P, Prod, MU);
    determinant(detmustart, MU);

    bool s = true;
    long factor;
    RR tr;
    while(s){
        s = false;
        Rearrange(m, P, Prod);
        for(long i = 1; i < m; i++){
            for(long j = 0; j < i; j++){
                round(tr, MU[P[i]][P[j]]);
                if(tr == 0){ continue; }
                negate(tr, tr);
                conv(factor, tr);
                Alter(m, B, Prod, MU, P[i], P[j], factor, U);
                s = true;
            }
        }
    }

    long loop = 0;
    s = true;
    long bestindex, bestfactor;
    while(s){ // Reduce until no changes happen
        s = false;
        Rearrange(m, P, Prod);
        for(long i = 0; i < m; i++){
            RowReduce(m, MU, P[i], bestindex, bestfactor);
            Alter(m, B, Prod, MU, P[i], bestindex, bestfactor, U);
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
