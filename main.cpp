#include <NTL/LLL.h>
#include "lalg.h"

void Test();

void GenerateTestData(std::ofstream& data, const long number = 50);

int main()
{
    Test();
    return 0;
}

void Test()
{
    std::ofstream datafile;
    datafile.open("input");
    if(datafile){
        GenerateTestData(datafile);
    }else{
        exit(1);
    }
    datafile.close();

    std::ifstream data;
    std::ofstream outlrd;
    std::ofstream outlll;
    data.open("input");
    outlrd.open("output.lrd");
    outlll.open("output.lll");

    long m, i;
    mat_ZZ B, B1;
    ZZ det, prod, small, t;
    RR smallr;

    while(!data.eof()){
        // Read basis dimension m, and the basis B
        data >> m;
        B.SetDims(m, m);
        B1.SetDims(m, m);
        data >> B;
        B1 = B;
        
        /* Reduce basis B 
         * record the length of shortest vector 
         * and the product of squares of all vectors
        */

        // LReduction part
        LReduction(B, 0);

        InnerProduct(small, B[0], B[0]);
        prod = small;
        for(i = 1; i < m; i++){
            InnerProduct(t, B[i], B[i]);
            mul(prod, prod, t);
            if(t < small) small = t;
        }
        conv(smallr, small);
        SqrRoot(smallr, smallr);
        outlrd << smallr << "\t" << prod << std::endl;

        // LLL algorithm part
        LLL(det, B1, 0);

        InnerProduct(small, B[0], B[0]);
        prod = small;
        for(i = 1; i < m; i++){
            InnerProduct(t, B[i], B[i]);
            mul(prod, prod, t);
            if(t < small) small = t;
        }
        conv(smallr, small);
        SqrRoot(smallr, smallr);
        outlll << smallr << "\t" << prod << std::endl;
    }

    data.close();
    outlrd.close();
    outlll.close();
}

// Generate random lattice bases
void GenerateTestData(std::ofstream& data, const long number)
{
    long vecm[5] = {25, 50, 75, 85, 100};
    long i, j, k, sign;
    ZZ det, bound;
    bound = 100;
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
