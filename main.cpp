#include "lalg.h"

void test();

int main()
{
    //std::ifstream infile;
    //infile.open("input");
    //long m, n, q;
    //float oprate;

    //m: challenge lattice dimension
    //n: reference dimension
    //q: modulus
    //B: challenge lattice basis

    //infile >> m >> n >> q;
    //mat_ZZ B;
    //infile >> B;
    //infile.close();

    //oprate = LReduction(B);

    //outfile.close();

    test();

    return 0;
}

void test()
{
    long m = 100;
    mat_ZZ B, U;
    B.SetDims(m, m);
    ident(U, m);
    // Generate random latice
    ZZ bound;
    bound = 100;
    long sign;
    for(long i = 0; i < m; i++){
        for(long j = 0; j < m; j++){
            RandomBnd(B[i][j], bound);
            RandomBnd(sign, 2);
            if(sign){
                negate(B[i][j], B[i][j]);
            }
        }
    }

    RR mu, detmustart, detmuend, optrate;
    vec_long P; // permutation of B
    mat_ZZ Prod; // innerproduct of vectors
    mat_RR MU, subMU; // mu matrix

    // Initialize all parameters
    P.SetLength(m);
    Prod.SetDims(m, m);
    MU.SetDims(m, m);
    subMU.SetDims(m, m);
    Init(B, m, P, Prod, MU);
    determinant(detmustart, MU);
    MUMax(m, MU, mu);

    std::cout
        << "---- Before Reduction ----\n"
        << "\nmu max: " << mu << "\tdet(mu): " << detmustart << std::endl
        << "---------------\n";

    bool s = true;
    long change = 0;
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
                change++;
            }
        }
    }

    determinant(detmuend, MU);
    div(optrate, detmuend, detmustart);
    MUMax(m, MU, mu);
    std::cout
        << "\n---- After OReduction ----\n"
        << "change: " << change 
        << "\nmu max: " << mu << "\tdet(mu): " << detmuend << "\toptimize rate: " << optrate << std::endl
        << "---------------\n";

    change = 0;
    long j;
    long beta = 4;
    vec_long indices, coeff;
    indices.SetLength(beta);
    coeff.SetLength(beta);
    s = true;
    // Reduce until no changes happen
    while(s){
        s = false;
        Rearrange(m, P, Prod);
        for(long i = 0; i < m; i++){
            Enumerate(m, MU, P[i], beta, indices, coeff);
            Alter(m, B, Prod, MU, P[i], beta, indices, coeff, U);
            for(j = 0; j < beta; j++){ if(coeff[j]) break; }
            if(j != beta){ s = true; change++; }
        }
    }

    determinant(detmuend, MU);
    div(optrate, detmuend, detmustart);
    MUMax(m, MU, mu);
    std::cout
        << "\n---- After LReduction ----\n"
        << "change: " << change 
        << "\nmu max: " << mu << "\tdet(mu): " << detmuend << "\toptimize rate: " << optrate << std::endl
        << "---------------\n";
}
