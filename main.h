#include <iostream>
#include <fstream>
#include <NTL/LLL.h>
#include "lalg.h"

void Test();
void GenerateTestData(std::ofstream& data, const long number = 5);
void Update(const mat_ZZ& B, const long& m, mat_ZZ& Prod, mat_RR& MU, RR& small);
