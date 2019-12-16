#include "lalg.h"

int main()//int argc, char* argv[])
{
    // Test part
    long m = 10;
    ZZ bound;
    bound = 10;

    int data[10][10] = 
    {
        {-8,-5,-4,-6,-4,-5,-1,-8,-9,5},
        {-7,6,7,-5,-1,-4,1,-1,0,-9},
        {-2,2,5,3,7,-5,-3,7,-9,-4},
        {-8,0,7,4,8,1,0,-9,3,-5},
        {9,4,-7,4,7,5,-5,-5,2,6},
        {9,-2,2,7,-7,9,9,-6,7,1},
        {-2,6,-5,-6,2,7,-8,0,9,-6},
        {3,-7,7,8,-6,0,-7,0,-8,6},
        {7,7,-9,5,-2,-8,0,-8,6,-5},
        {-3,-6,5,-8,8,0,7,4,3,-6}
    };

    long sign;
    mat_ZZ B;
    B.SetDims(m, m);
    // Generate random latice
    for(long i = 0; i < m; i++){
        for(long j = 0; j < m; j++){
            conv(B[i][j], data[i][j]);
            //RandomBnd(B[i][j], bound);
            //RandomBnd(sign, 2);
            //if(sign){
            //    negate(B[i][j], B[i][j]);
            //}
        }
    }

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
    MUMax(m, MU, mu);

    std::cout
        << "\n---------------\n"
        << "basis:" << std::endl << B << std::endl
        << "innerproduct matrix:" << std::endl << Prod << std::endl
        << "MU matrix:" << std::endl << MU << std::endl
        << "mu max: " << mu << std::endl
        << "det(mu): " << detmustart << std::endl 
        << "---------------\n";

    // Reduce while mu < muold
    long loop = 0;
    long bestindex, bestfactor;
    bool s = true;
    while(s && loop < 1000){
        s = false;
        Rearrange(m, P, Prod);
        MUMax(m, MU, mu);
        for(long i = 0; i < m; i++){
            RowReduce(m, MU, mu, P[i], bestindex, bestfactor);
            Alter(B, m, Prod, MU, P[i], bestindex, bestfactor);
            if(bestfactor != 0){ 
                s = true;
                std::cout << "B[" << P[i] << "] += " << bestfactor << "B[" << bestindex << "]\n";
            }
        }
        std::cout << "------------------\n\n";
        loop++;
    }

    MUMax(m, MU, mu);
    determinant(detmuend, MU);
    div(optrate, detmuend, detmustart);
    std::cout
        << "\n---------------\n"
        << "loop count: " << loop << std::endl
        << "basis:" << std::endl << B << std::endl
        << "innerproduct matrix:" << std::endl << Prod << std::endl
        << "MU matrix:" << std::endl << MU << std::endl
        << "mu max: " << mu << std::endl
        << "det(mu): " << detmuend << std::endl 
        << "optimize rate: " << optrate << std::endl
        << "---------------\n";
    // Test part end

    //std::ofstream outfile;
    //outfile.open("output");
    //std::ifstream infile;
    //infile.open("input");

    //long m, n, q;
    //float oprate;

    /* m: challenge lattice dimension
n: reference dimension
q: modulus
B: challenge lattice basis
*/
    //infile >> m >> n >> q;
    //mat_ZZ B;
    //infile >> B;
    //infile.close();

    //oprate = LReduction(B);

    //outfile.close();

    return 0;
}
