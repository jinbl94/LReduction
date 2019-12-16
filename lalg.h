#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>

using namespace NTL;

void init(const mat_ZZ& B, const long m, vec_long& P, mat_ZZ& Prod, mat_RR & MU);

void reduce(const long m, vec_long& P, const mat_RR& MU, const long& k, const long& r, long& j, long& q);

void alter(mat_ZZ& B, const long& m, const vec_long& P, mat_ZZ& Prod, mat_RR& MU, const long& k, const long& j, const long& q);

void rearrange(const long& m, vec_long& P, const mat_ZZ& Prod);

void mumax(const long& m, const mat_RR& MU, RR& mu);

void MaxRow(const long m, const vec_RR& row, const long& r, RR& maxrow);

void FactorRange(const long& m, const mat_RR& MU, const long& row, const long& currentrow,
        const RR& maxrow, long& maxfactor, long& minfactor);

void BestFactor(const long& m, const mat_RR& MU, const long& row, const long& currentrow, RR& oldmax,
        const long& minfactor, const long& maxfactor, long& j, long& q);

static float lreduction(mat_ZZ& B);
