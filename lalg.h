#ifndef _HEAD_LALG
#define _HEAD_LALG
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/mat_RR.h>

using namespace NTL;

// Find the greatest absolute value of elements of a single row (besides 1, the rth element)
void MaxRow(const long& m, const vec_RR& vecrow, RR& maxrow);
void MaxRowR(const long& m, const long& row, const vec_RR& vecrow, RR& maxrow);

// Innerproduct without the rth element
void InnerProductR(const long& m, const long& r, RR& result, const vec_RR& vec1, const vec_RR& vec2);

// Enumerate all possible combinations of some basis vectors, find the shortest one, and store this combination
// MU : the basis vectors
// s : dimension of the vectors
// l : block size of enumeration
void EnumerateCPMU(const long& m, const mat_RR& MU, const long& beta, const vec_long& indices, vec_long& coeff);
void Enumerate(const long& m, const mat_RR& MU, const long& r, const long& beta, vec_long& indices, vec_long& coeff);

// Initialize all parameters
void Init(const mat_ZZ& B, const long& m, vec_long& P, mat_ZZ& Prod, mat_RR & MU);

// Apply the combinatino to the basis and update related parameters
bool Alter(const long& m, mat_ZZ& B, mat_ZZ& Prod, mat_RR& MU, const long& r, const long& beta, const vec_long& indices, const vec_long& coeff, mat_ZZ* U);
bool Alter(const long& m, mat_ZZ& B, mat_ZZ& Prod, mat_RR& MU, const long& row, const long& j, const long& q, mat_ZZ* U);

// Rearrange the basis according to some standard
void Rearrange(const long& m, vec_long& P, const mat_ZZ& Prod);

// Find the greatest absolute value of elements of MU matrix (besides 1)
void MUMax(const long& m, const mat_RR& MU, RR& mu);

// According to the greatest absolute value, find acceptable range for linear combinations
void FactorRange(const long& m, const mat_RR& MU, const long& row, const long& currentrow,
        const RR& maxrow, long& maxfactor, long& minfactor);

// Enumerate all acceptable combinations, and find the one minimize the greatest absolute value
void BestFactor(const long& m, const mat_RR& MU, const long& row, const long& currentrow, RR& oldmax,
        const long& minfactor, const long& maxfactor, long& bestindex, long& bestfactor);

// Find the beset index and factor
void RowReduce(const long& m, const mat_RR& MU, const long& row, long& bestindex, long& bestfactor);

// L-Reduction main procedure
static float LReduction(mat_ZZ& B, mat_ZZ* U, const long beta = 6);

#else
#endif
