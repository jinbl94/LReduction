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
    mat_RR B;
    B.SetDims(m, m);
    for(long i = 0; i < m; i++){
        for(long j = 0; j < m; j++){
            //RandomBnd(B[i][j], bound);
            random(B[i][j]);
        }
    }
    B[r][r] = 1;
    RR maxrow;
    MaxRow(m, B[r], r, maxrow);
    //infile >> B;
    //infile.close();

    //oprate = LReduction(B);

    //outfile.close();

    return 0;
}
