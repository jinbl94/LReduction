#include <iostream>
#include <fstream>
#include <string>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>

using namespace NTL;

void init(const mat_ZZ& B, const long m, vec_long& P, mat_ZZ& Prod, mat_RR & MU, RR& mu);

void reduce(const long m, vec_long& P, const mat_RR& MU, const long& k, const long& r, long& j, long& q);

void alter(mat_ZZ& B, const long& m, const vec_long& P, mat_ZZ& Prod, mat_RR& MU, const long& k, const long& j, const long& q);

void rearrange(const long& m, vec_long& P, const mat_RR& MU, const long& r);

void mumax(const long& m, const mat_RR& MU, RR& mu);

static float lreduction(mat_ZZ& B);

int main()//int argc, char* argv[])
{
    long m = 5;
    ZZ bound;
    bound = 10;
    mat_ZZ B;
    B.SetDims(m, m);

    for(long i = 0; i < m; i++){
        for(long j = 0; j < m; j++){
            RandomBnd(B[i][j], bound);
        }
    }

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

    init(B, m, P, Prod, MU, mu);
    determinant(detmu_start, MU);

    std::cout
        << "basis:\n" << B << std::endl
        << "innerproduct:\n" << Prod << std::endl
        << "mu matrix:\n" << MU << std::endl
        << "\n\nmu max:\n" << mu << std::endl
        << std::endl;

    long r, k, j, q;
    for(long loop = 0; loop < 3; loop++){
        for(r = 0; r < m; r++){
            for(k = r + 1; k < m; k++){
                reduce(m, P, MU, k, r, j, q);
                alter(B, m, P, Prod, MU, k, j, q);
            }
        }
    }

    determinant(detmu_end, MU);
    div(optrate, detmu_end, detmu_start);
    mumax(m, MU, mu);

    std::cout << "------------------------" << std::endl;

    std::cout
        //<< "basis:\n" << B
        << "\n\ninnerproduct:\n" << Prod
        << "\n\nmu matrix:\n" << MU
        << "\n\nmu max:\n" << mu
        << std::endl;

    std::cout << "optimize rate: " << optrate << std::endl;

    //std::ofstream outfile;
    //outfile.open("output");
    //std::ifstream infile;
    //infile.open("input");

    //long m, n, q;
    //float p;

    /* m: challenge lattice dimension
n: reference dimension
q: modulus
B: challenge lattice basis */
    //infile >> m >> n >> q;
    //mat_ZZ B;
    //B.SetDims(m, m);
    //infile >> B;
    //infile.close();

    //p = lreduction(B);

    //outfile.close();

    return 0;
}

void init(const mat_ZZ& B, const long m, vec_long& P, mat_ZZ& Prod, mat_RR & MU, RR& mu) // init all parameters
{
    mu = 0;
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
            div(tr1, tr2, tr1);
            MU[i][j] = tr1;
            abs(tr1, tr1);
            if(mu < tr1){
                mu = tr1;
            }

            // MU[i][j] = tr1 = Prod[i][j]/Prod[j][j]
            conv(tr1, Prod[j][i]);
            conv(tr2, Prod[i][i]);
            div(tr1, tr2, tr1);
            MU[j][i] = tr1;
            abs(tr1, tr1);
            if(mu < tr1){
                mu = tr1;
            }
        }
    }
}

void rearrange(const long& m, vec_long& P, const mat_RR& MU, const long& r)
{
    RR tr;
    long tl, i, j;
    for(i = r + 2; i < m; i++){
        tr = MU[P[i]][P[r]];
        tl = P[i];
        for(j = i - 1; j > r; j--){
            if(tr < MU[P[j]][P[r]]){
                P[j + 1] = P[j];
            }else{
                break;
            }
        }
        P[j + 1] = tl;
    }
}

void alter(mat_ZZ& B, const long& m, const vec_long& P, mat_ZZ& Prod, mat_RR& MU, const long& k, const long& j, const long& q)
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
        //mul(tz, q * q, Prod[P[j]][P[j]]);
        //add(Prod[P[k]][P[k]], Prod[P[k]][P[k]], tz);
        //mul(tz, 2 * q, Prod[P[k]][P[j]]);
        //add(Prod[P[k]][P[k]], Prod[P[k]][P[k]], tz);
        InnerProduct(Prod[P[k]][P[k]], B[P[k]], B[P[k]]);

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
    long tq;
    RR alpha = abs(MU[P[k]][P[r]]), talpha, tr;

    for(long l = 0; l < k; l++){
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
        }else{}
    }
    for(long l = k + 1; l < m; l++){
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

void mumax(const long& m, const mat_RR& MU, RR& mu)
{
    mu = 0;
    RR tr;
    for(long i = 0; i < m; i++){
        for (long j = 0; j < i; j++){
            abs(tr, MU[i][j]);
            if(mu < tr){
                mu = tr;
            }
        }
    }
}

static float lreduction(mat_ZZ& B)
{
    const long m = B.NumRows();
    const long n = B.NumCols();

    float result = 0;

    if(m == n){
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

        init(B, m, P, Prod, MU, mu);
        determinant(detmu_start, MU);

        long loop = 0;
        long r, k, j, q;
        while(mu < mu_old){
            mu_old = mu;
            for(r = 0; r < m; r++){
                for(k = r + 1; k < m; k++){
                    reduce(m, P, MU, k, r, j, q);
                    alter(B, m, P, Prod, MU, k, j, q);
                }
                rearrange(m, P, MU, r);
            }
            mumax(m, MU, mu);
            loop++;
        }

        determinant(detmu_end, MU);
        div(optrate, detmu_end, detmu_start);
        conv(result, optrate);
    }else{
        std::cerr << "Input matirx is not square" << std::endl;
    }

    return result;
}
