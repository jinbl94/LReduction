#include "main.h"

int main()
{
    //std::ofstream datafile;
    //datafile.open("input");
    //if(datafile){ GenerateTestData(datafile, 5); }
    //else{ exit(1); }
    //datafile.close();

    std::cout << "***** Test start *****" << std::endl;
    Test();
    std::cout << "***** Test end   *****" << std::endl;
    return 0;
}

void Test()
{
    std::ifstream data;
    std::ofstream outlrd;
    data.open("input");
    outlrd.open("output.lrd");

    long m;
    mat_ZZ B, B1, Prod;
    mat_RR MU;
    ZZ det;
    RR small, detmu;

    outlrd << "Dim\tOrigi\t\tLLL\t\tLRD" << std::endl;
    data >> m;
    B.SetDims(m, m);
    B1.SetDims(m, m);
    Prod.SetDims(m, m);
    MU.SetDims(m, m);
    data >> B;
    B1 = B;

    while(!data.eof()){
        outlrd << m << "\t";

        // Before reduction
        Update(B, m, Prod, MU, small);
        SqrRoot(small, small);
        determinant(detmu, MU);
        outlrd << small << "/" << detmu << "\t";

        /* Reduce basis B 
         * record the length of shortest vector 
         * and the product of squares of all vectors
         */

        // LLL algorithm part
        LLL(det, B1, 0);
        Update(B1, m, Prod, MU, small);
        SqrRoot(small, small);
        determinant(detmu, MU);
        outlrd << small << "/" << detmu << "\t";

        // LRD
        LReduction(B, 6, 0);
        Update(B, m, Prod, MU, small);
        SqrRoot(small, small);
        determinant(detmu, MU);
        outlrd << small << "/" << detmu << std::endl;

        outlrd.flush();
        std::cout << m << " ";
        std::cout.flush();

        // Release resources
        B.kill();
        B1.kill();
        Prod.kill();
        MU.kill();

        // Read basis dimension m, and the basis B
        data >> m;
        B.SetDims(m, m);
        B1.SetDims(m, m);
        Prod.SetDims(m, m);
        MU.SetDims(m, m);
        data >> B;
        B1 = B;
    }

    std::cout << std::endl;
    data.close();
    outlrd.close();
}

// Generate random lattice bases
void GenerateTestData(std::ofstream& data, const long number)
{
    long vecm[1] = {20};
    long i, j, k, sign;
    ZZ det, bound;
    bound = 50;
    mat_ZZ B;

    for(auto m : vecm){
        for(i = 0; i < number; i++){
            B.SetDims(m, m);
            for(j = 0; j < m; j++){
                for(k = 0; k < m; k++){
                    RandomBnd(B[j][k], bound);
                    RandomBnd(sign, 2);
                    if(sign){
                        negate(B[j][k], B[j][k]);
                    }
                }
            }
            determinant(det, B);
            // If B is not m 
            if(det == 0) { i--; continue; }
            data << m << std::endl << B << std::endl;
        }
    }
}

void Update(const mat_ZZ& B, const long& m, mat_ZZ& Prod, mat_RR& MU, RR& small)
{
    RR tr, tr1, t;

    MU[0][0] = 1;
    InnerProduct(Prod[0][0], B[0], B[0]); // Prod[i][i] = B[i]^2
    conv(small, Prod[0][0]);
    for(long i = 1; i < m; i++){
        MU[i][i] = 1;
        InnerProduct(Prod[i][i], B[i], B[i]); // Prod[i][i] = B[i]^2
        conv(t, Prod[i][i]);
        if(t < small) small = t;
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
