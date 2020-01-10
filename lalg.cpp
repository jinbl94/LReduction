#include "lalg.h"

static void MaxRow(const long& m, const vec_RR& vecrow, RR& maxrow)
{
    RR tr;
    maxrow = 0;
    for(long col = 0; col < m; col++){
        abs(tr, vecrow[col]);
        if(maxrow < tr){ maxrow = tr; }
    }
}

static void MaxRowR(const long& m, const long& row, const vec_RR& vecrow, RR& maxrow)
{
    RR tr;
    maxrow = 0;
    for(long col = 0; col < m; col++){
        if(col == row){ continue; }
        abs(tr, vecrow[col]);
        if(maxrow < tr){ maxrow = tr; }
    }
}

static void InnerProductR(const long& m, const long& r, RR& result, const vec_RR& vec1, const vec_RR& vec2)
{
    RR t;
    conv(result, 0);
    for(long i = 0; i < m; i++){
        if(i == r) continue;
        mul(t, vec1[i], vec2[i]);
        add(result, result, t);
    }
}

static void Enumerate(const long& m, const mat_RR& MU, const long& r, const long& beta, vec_long& indices, vec_long& coeff)
{
    assert(indices.length() == beta);
    assert(coeff.length() == beta);

    // indices, bound vector for indices, coefficient vector, lower and upper bound for coefficient vector
    // deltax, delta2x, vecl, vecc are temp variables use in enum
    vec_long ind, indbound, coe, P, deltax, delta2x; // lowbound, upbound;
    mat_RR GSO, GSOmu, matv; // GSO basis realted
    vec_RR GSOsquare, tvecrr, tvecrr1, vecl, vecc, vecy;
    RR radius, minmax, tr, tr1; // Temp values
    long i, j; // temp variable used as index

    P.SetLength(m);
    for(i = 1; i < m; i++){ P[i] = i; }
    P[0] = r; P[r] = 0;

    tvecrr.SetLength(m);
    tvecrr1.SetLength(m);

    coe.SetLength(beta);
    ind.SetLength(beta);
    indbound.SetLength(beta);
    GSOsquare.SetLength(beta + 1);

    vecy.SetLength(beta);
    vecl.SetLength(beta + 1);
    vecc.SetLength(beta);
    deltax.SetLength(beta);
    delta2x.SetLength(beta);

    GSO.SetDims(beta + 1, m);
    GSOmu.SetDims(beta + 1, beta + 1);
    matv.SetDims(beta + 1, m);

    // Target at index beta
    GSO[beta] = MU[P[0]];
    MaxRowR(m, r, GSO[beta], minmax);
    sqr(radius, minmax);
    mul(radius, m - 1, radius);
    //InnerProductR(m, r, radius, GSO[0], GSO[0]);

    // ind = [1, 2, ..., beta], indbound = [m - beta, m - beta + 1, ..., m - 1]
    for(i = 0; i < beta; i++){
        ind[i] = i + 1;
        indbound[i] = m - beta + i;
        coeff[i] = 0;
    }

    // Run over all possible indices
    ind[beta - 1]--;
    while(ind != indbound){
        for(i = beta - 1; i >= 0; i--){ if(ind[i] != indbound [i]) break; }
        ind[i]++; ++i;
        for(; i < beta; i++){ ind[i] = ind[i - 1] + 1; }

        /* Schnorr's enumeration algorithm
         * calculate gram-schmit, and related parameters
         * enumerate them and store the shortest one
         */
        InnerProductR(m, r, GSOsquare[0], GSO[0], GSO[0]);
        for(i = 0; i < beta; i++){
            GSO[i] = MU[P[ind[i]]];
            for(j = 0; j < i; j++){
                InnerProductR(m, r, GSOmu[i][j], GSO[i], GSO[j]);
                div(GSOmu[i][j], GSOmu[i][j], GSOsquare[j]);
                mul(tvecrr, GSOmu[i][j], GSO[j]);
                sub(GSO[i], GSO[i], tvecrr);
            }
            InnerProductR(m, r, GSOsquare[i], GSO[i], GSO[i]);
        }

        // GSO[beta]
        InnerProductR(m, r, vecl[0], GSO[beta], GSO[beta]);
        for(j = 0; j < beta; j++){
            InnerProductR(m, r, GSOmu[beta][j], GSO[beta], GSO[j]);
            div(GSOmu[beta][j], GSOmu[beta][j], GSOsquare[j]);
            mul(tvecrr, GSOmu[beta][j], GSO[j]);
            sub(GSO[beta], GSO[beta], tvecrr);
            InnerProductR(m, r, vecl[j + 1], GSO[beta], GSO[beta]);
        }
        GSOsquare[beta] = vecl[beta];

        // Initialize all variables
        for(i = 0; i < beta; i++){
            coe[i] = 0;
            deltax[i] = 0;
            delta2x[i] = -1;
            negate(vecc[i], GSOmu[beta][i]);
            vecy[i] = 0;
            matv[i] = MU[P[0]];
        }
        matv[beta] = MU[P[0]];
        deltax[0] = 1;
        delta2x[0] = 1;

        i = 0;
        while(true){
            // Length of current combination
            sub(vecy[i], coe[i], vecc[i]);
            abs(vecy[i], vecy[i]);
            mul(tr, vecy[i], vecy[i]);
            mul(tr, tr, GSOsquare[i]);
            add(vecl[i], vecl[i + 1], tr);
            mul(tvecrr, coe[i], MU[P[ind[i]]]);
            add(matv[i], matv[i + 1], tvecrr);
            MaxRowR(m, r, matv[i], tr);

            // If current combination is smller in l_\infnity norm
            if(tr < minmax && i == 1){
                for(j = 0; j < beta; j++) indices[j] = P[ind[j]];
                coeff = coe;
                minmax = tr;
                sqr(radius, minmax);
                mul(radius, m - 1, radius);
            }

            if(vecl[i] <= radius && i > 1){
                --i;
                negate(vecc[i], GSOmu[beta][i]);
                for(j = i + 1; j < beta; j++){
                    mul(tr, coe[j], GSOmu[j][i]);
                    sub(vecc[i], vecc[i], tr);
                }
                round(tr, vecc[i]);
                conv(coe[i], tr);
                deltax[i] = 0;
                delta2x[i] = (vecc[i] < coe[i]) ? 1 : -1;
            }else if(vecl[i] > radius && i == (beta - 1)){ break; }
            else{
                ++i;
                delta2x[i] = - delta2x[i];
                deltax[i] = - deltax[i] + delta2x[i];
                coe[i] = coe[i] + deltax[i];
            }
        }
    }

    tvecrr.kill();
    tvecrr1.kill();
    coe.kill();
    ind.kill();
    indbound.kill();
    GSOsquare.kill();
    vecy.kill();
    vecl.kill();
    vecc.kill();
    deltax.kill();
    delta2x.kill();
    GSO.kill();
    GSOmu.kill();
    matv.kill();
}

static void Init(const mat_ZZ& B, const long& m, vec_long& P, mat_ZZ& Prod, mat_RR & MU)
{
    RR tr, tr1;

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
            conv(tr, Prod[i][j]);
            conv(tr1, Prod[j][j]);
            div(tr, tr, tr1);
            MU[i][j] = tr;
            abs(tr, tr);
            // MU[j][i] = tr = Prod[i][j]/Prod[i][i]
            conv(tr, Prod[j][i]);
            conv(tr1, Prod[i][i]);
            div(tr, tr, tr1);
            MU[j][i] = tr;
            abs(tr, tr);
        }
    }
}

static void Rearrange(const long& m, vec_long& P, const mat_ZZ& Prod)
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

static bool Alter(const long& m, mat_ZZ& B, mat_ZZ& Prod, mat_RR& MU, const long& r, const long& beta, const vec_long& indices, const vec_long& coeff, mat_ZZ* U)
{
    long j;
    for(j = 0; j < beta; j++){ if(coeff[j]) break; }
    if(j != beta){
        vec_ZZ tveczz;
        tveczz.SetLength(m);
        RR tr, tr1;
        long index;
        // B[r] = B[r] + B[indices] * coeff
        for(index = 0; index < beta; index++){
            mul(tveczz, coeff[index], B[indices[index]]);
            add(B[r], B[r], tveczz);
        }
        // Prod[r][index] = Prod[index][r] = B[r] * B[index]
        for(index = 0; index < m; index++){
            InnerProduct(Prod[r][index], B[r], B[index]);
            Prod[index][r] = Prod[r][index];
        }
        for(index = 0; index < m; index++){
            // MU[index][r] = Prod[index][r] / Prod[r][r]
            conv(tr, Prod[r][index]);
            conv(tr1, Prod[r][r]);
            div(MU[index][r], tr, tr1);
            // MU[r][index] = Prod[r][index] / Prod[index][index]
            conv(tr1, Prod[index][index]);
            div(MU[r][index], tr, tr1);
        }    
        if(U){
            // U[r] = U[r] + U[indices] * coeff
            for(index = 0; index < beta; index++){
                mul(tveczz, coeff[index], (*U)[indices[index]]);
                add((*U)[r], (*U)[r], tveczz);
            }
        }
        return true;
    }
    return false;
}

static bool Alter(const long& m, mat_ZZ& B, mat_ZZ& Prod, mat_RR& MU, const long& row, const long& j, const long& q, mat_ZZ* U)
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
            // MU[r][row] = Prod[r][row]/Prod[row][row]
            conv(tr, Prod[index][row]);
            conv(tr1, Prod[row][row]);
            div(MU[index][row], tr, tr1);
        }    

        if(U){
            // U[row] = U[row] + qU[j]
            for(index = 0; index < m; index++){
                mul(tz, q, (*U)[j][index]);
                add((*U)[row][index], (*U)[row][index], tz);
            }
        }
        return true;
    }
    return false;
}

static void FactorRange(const long& m, const mat_RR& MU, const long& row, const long& currentrow,
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
            }else{
                if(maxfactor > tr1){ conv(maxfactor, tr1); }
                if(minfactor < tr){ conv(minfactor, tr); }
            }
        }
    }
}

static void BestFactor(const long& m, const mat_RR& MU, const long& row, const long& currentrow, RR& maxrow,
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

static void RowReduce(const long& m, const mat_RR& MU, const long& row, long& bestindex, long& bestfactor)
{
    RR maxrow;
    long minfactor, maxfactor;
    bestindex = 0;
    bestfactor = 0;
    // find max absolute value of current row
    // MaxRowR(m, row, MU[row], maxrow);
    MUMax(m, MU, maxrow);
    for(long currentrow = 0; currentrow < m; currentrow++){
        if(currentrow == row){ continue; }
        FactorRange(m, MU, row, currentrow, maxrow, maxfactor, minfactor);
        if(minfactor > maxfactor){ continue; }
        BestFactor(m, MU, row, currentrow, maxrow, minfactor, maxfactor, bestindex, bestfactor);
    }
}

static void MUMax(const long& m, const mat_RR& MU, RR& mu)
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

void LReduction(mat_ZZ& B, const long& beta, mat_ZZ* U) // B' = U * B
{
    assert( B.NumRows() == B.NumCols()); // make sure B is a square matrix
    const long m = B.NumRows();

    vec_long P; // permutation of B
    mat_ZZ Prod; // innerproduct of vectors
    mat_RR MU; // mu matrix

    // Initialize all parameters
    P.SetLength(m);
    Prod.SetDims(m, m);
    MU.SetDims(m, m);
    Init(B, m, P, Prod, MU);

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
                s = Alter(m, B, Prod, MU, P[i], P[j], factor, U);
            }
        }
    }

    s = true;
    long loop = 0;
    const long LoopBound = 50;
    vec_long indices, coeff;
    indices.SetLength(beta);
    coeff.SetLength(beta);
    // Reduce until no change happened
    while(s && loop < LoopBound){
        s = false;
        // Rearrange(m, P, Prod);
        for(long i = 0; i < m; i++){
            Enumerate(m, MU, i, beta, indices, coeff);
            if(Alter(m, B, Prod, MU, i, beta, indices, coeff, U)) s = true;
        }
        loop++;
    }
}
