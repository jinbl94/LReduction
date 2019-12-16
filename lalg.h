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
#include <NTL/mat_RR.h>

using namespace NTL;

// L-Reduction main procedure
static float LReduction(mat_ZZ& B);

// Initialize all parameters
void Init(const mat_ZZ& B, const long m, vec_long& P, mat_ZZ& Prod, mat_RR & MU);

// L-Reduction core procedure, wich finds the best linear combination according to some parameters
// Deprecated
void Reduce(const long m, vec_long& P, const mat_RR& MU, const long& k, const long& r, long& j, long& q);

// Apply the combinatino to the basis and update related parameters
void Alter(mat_ZZ& B, const long& m, mat_ZZ& Prod, mat_RR& MU, const long& row, const long& j, const long& q);

// Rearrange the basis according to some standard
void Rearrange(const long& m, vec_long& P, const mat_ZZ& Prod);

// Find the greatest absolute value of elements of MU matrix (besides 1)
void MUMax(const long& m, const mat_RR& MU, RR& mu);

// Find the greatest absolute value of elements of a single row (besides 1, the rth element)
void MaxRow(const long m, const vec_RR& row, const long& r, RR& maxrow);

// According to the greatest absolute value, find acceptable range for linear combinations
void FactorRange(const long& m, const mat_RR& MU, const long& row, const long& currentrow,
        const RR& maxrow, long& maxfactor, long& minfactor);

// Enumerate all acceptable combinations, and find the one minimize the greatest absolute value
void BestFactor(const long& m, const mat_RR& MU, const long& row, const long& currentrow, RR& oldmax,
        const long& minfactor, const long& maxfactor, long& j, long& q);
#else
#endif
