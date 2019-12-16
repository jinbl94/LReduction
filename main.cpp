#include "lalg.h"

int main()//int argc, char* argv[])
{
    long m = 5;
    long r = 3;
    ZZ bound;
    bound = 10;

    //std::ofstream outfile;
    //outfile.open("output");
    //std::ifstream infile;
    //infile.open("input");

    //long m, n, q;
    float oprate;

    /* m: challenge lattice dimension
       n: reference dimension
       q: modulus
       B: challenge lattice basis
       */
    //infile >> m >> n >> q;
    long sign;
    mat_ZZ B;
    B.SetDims(m, m);
    for(long i = 0; i < m; i++){
        for(long j = 0; j < m; j++){
            //RandomBnd(B[i][j], bound);
            RandomBnd(B[i][j], bound);
            RandomBnd(sign, 2);
            if(sign){
                negate(B[i][j], B[i][j]);
            }
        }
    }

    //infile >> B;
    //infile.close();

    //oprate = LReduction(B);

    //outfile.close();

    return 0;
}
