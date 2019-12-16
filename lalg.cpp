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
            }else{
                break;
            }
        }
        P[j + 1] = tl;
    }
}

void Alter(mat_ZZ& B, const long& m, const vec_long& P, mat_ZZ& Prod, mat_RR& MU, const long& k, const long& j, const long& q)
{
    if(q){
        ZZ tz;
        RR tr, tr1;
        long r;
        // B[P[k]] = B[P[k]] + q * B[P[j]]
        for(r = 0; r < m; r++){
            mul(tz, q, B[P[j]][P[r]]);
            sub(B[P[k]][P[r]], B[P[k]][P[r]], tz);
        }

        // Prod[P[k]][P[k]] = Prod[P[k]][P[k]]^2 + q^2 Prod[P[j]][P[j]] + 2q Prod[P[k]][P[j]]
        mul(tz, q * q, Prod[P[j]][P[j]]);
        add(Prod[P[k]][P[k]], Prod[P[k]][P[k]], tz);
        mul(tz, 2 * q, Prod[P[k]][P[j]]);
        sub(Prod[P[k]][P[k]], Prod[P[k]][P[k]], tz);

        // Prod[P[r]][P[k]] = Prod[P[k]][P[r]] = Prod[P[k]][P[r]] + q * Prod[P[j]][P[r]]
        for(r = 0; r < k; r++){
            mul(tz, q, Prod[P[j]][P[r]]);
            sub(Prod[P[k]][P[r]], Prod[P[k]][P[r]], tz);
            Prod[P[r]][P[k]] = Prod[P[k]][P[r]];
        }
        for(r = k + 1; r < m; r++){
            mul(tz, q, Prod[P[j]][P[r]]);
            sub(Prod[P[k]][P[r]], Prod[P[k]][P[r]], tz);
            Prod[P[r]][P[k]] = Prod[P[k]][P[r]];
        }

        // MU[P[k]][P[r]] = MU[P[k]][P[r]] + q MU[P[j]][P[r]]
        // MU[P[r]][P[k]] = Prod[P[r]][P[k]]/Prod[P[k]][P[k]]
        for(r = 0; r < m; r++){
            mul(tr, q, MU[P[j]][P[r]]);
            sub(MU[P[k]][P[r]], MU[P[k]][P[r]], tr);
            abs(tr, MU[P[k]][P[r]]);

            conv(tr, Prod[P[r]][P[k]]);
            conv(tr1, Prod[P[k]][P[k]]);
            div(MU[P[r]][P[k]], tr, tr1);
            abs(tr, MU[P[r]][P[k]]);
        }    
    }else{
        // do nothing
    }
}

void reduce(const long m, vec_long& P, const mat_RR& MU, const long& k, const long& r, long& j, long& q)
{
    j = 1;
    q = 0;
    long tq, l;
    RR alpha = abs(MU[P[k]][P[r]]), talpha, tr;

    div(tr, MU[P[k]][P[r]], MU[P[r]][P[r]]);
    round(tr, tr);
    conv(tq, tr);
    mul(tr, tr, MU[P[r]][P[r]]);
    sub(tr, MU[P[k]][P[r]], tr);
    abs(talpha, tr);
    if(alpha > talpha){
        j = r;
        q = tq;
        alpha = talpha;
    }

    for(l = k + 1; l < m; l++){
        if(MU[P[l]][P[r]] != 0){
            div(tr, MU[P[k]][P[r]], MU[P[l]][P[r]]);
            round(tr, tr);
            conv(tq, tr);
            mul(tr, tr, MU[P[l]][P[r]]);
            sub(tr, MU[P[k]][P[r]], tr);
            abs(talpha, tr);
            if(alpha > talpha){
                j = l;
                q = tq;
                alpha = talpha;
            }
        }else{
            // do nothing
        }
    }
}

void MaxRow(const long m, const vec_RR& row, const long& r, RR& maxrow)
{
    RR tr;

    maxrow = 0;
    for(long col = 0; col < r; col++){
        abs(tr, row[col]);
        if(maxrow < tr){
            maxrow = tr;
        }
    }
    for(long col = r + 1; col < m; col++){
        abs(tr, row[col]);
        if(maxrow < tr){
            maxrow = tr;
        }
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
            if(maxfactor > tr){
                conv(maxfactor, tr);
            }
            negate(tr, maxrow);
            sub(tr, tr, a);
            div(tr, tr, b);
            round(tr, tr);
            if(minfactor < tr){
                conv(minfactor, tr);
            }
        }else if(b < 0){
            sub(tr, maxrow, a);
            div(tr, tr, b);
            round(tr, tr);
            if(minfactor < tr){
                conv(minfactor, tr);
            }
            negate(tr, maxrow);
            sub(tr, tr, a);
            div(tr, tr, b);
            round(tr, tr);
            if(maxfactor > tr){
                conv(maxfactor, tr);
            }
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
            if(maxfactor > tr){
                conv(maxfactor, tr);
            }
            negate(tr, maxrow);
            sub(tr, tr, a);
            div(tr, tr, b);
            round(tr, tr);
            if(minfactor < tr){
                conv(minfactor, tr);
            }
        }else if(b < 0){
            sub(tr, maxrow, a);
            div(tr, tr, b);
            round(tr, tr);
            if(minfactor < tr){
                conv(minfactor, tr);
            }
            negate(tr, maxrow);
            sub(tr, tr, a);
            div(tr, tr, b);
            round(tr, tr);
            if(maxfactor > tr){
                conv(maxfactor, tr);
            }
        }
    }
}

void BestFactor(const long& m, const mat_RR& MU, const long& row, const long& currentrow, RR& oldmax,
        const long& minfactor, const long& maxfactor, long& j, long& q)
{
    RR newmax, a, b, tr;
    long col;
    for(long factor = minfactor; factor <= maxfactor; factor++){
        newmax = 0;

        for(col = 0; col < row; col++){
            a = MU[row][col];
            b = MU[currentrow][col];
            mul(tr, factor, b);
            add(tr, tr, a);
            abs(tr, tr);
            if(newmax < tr){
                newmax = tr;
            }
        }
        for(col = row + 1; col < m; col++){
            a = MU[row][col];
            b = MU[currentrow][col];
            mul(tr, factor, b);
            add(tr, tr, a);
            abs(tr, tr);
            if(newmax < tr){
                newmax = tr;
            }
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
            if(mu < tr){
                mu = tr;
            }
            abs(tr, MU[j][i]);
            if(mu < tr){
                mu = tr;
            }
        }
    }
}

static float LReduction(mat_ZZ& B)
{
    const long m = B.NumRows();
    const long n = B.NumCols();
    long r, k, j, q, loop = 0;
    float result = 0;

    assert( m == n);

    // permutation of B
    vec_long P;
    P.SetLength(m);
    // innerproduct of vectors
    mat_ZZ Prod;
    Prod.SetDims(m, m);
    // \mu matrix
    mat_RR MU;
    MU.SetDims(m, m);
    // \mu_{max};
    RR mu, mu_old, detmu_start, detmu_end, optrate;

    Init(B, m, P, Prod, MU);
    MUMax(m, MU, mu);
    determinant(detmu_start, MU);

    std::cout
        << "basis:" << std::endl << B << std::endl
        << "innerproduct matrix:" << std::endl << Prod << std::endl
        << "mu max: " << mu << std::endl
        << "det(mu): " << detmu_start << std::endl << std::endl;

    while(loop < 10){
        for(r = 0; r < m; r++){
            Rearrange(m, P, Prod);
            RowReduce(m, MU, P[r], j, q);
            Alter(B, m, P, Prod, MU, k, j, q);
            //for(k = r + 1; k < m; k++){
            //    reduce(m, P, MU, k, r, j, q);
            //    Alter(B, m, P, Prod, MU, k, j, q);
            //}
        }
        loop++;
    }

    determinant(detmu_end, MU);
    MUMax(m, MU, mu);
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
