#ifndef DistAlnF
#define DistAlnF

#include <stdlib.h>
#include <stdio.h>
#include "Alignment.h"
#include "Distance.h"

typedef TypeFloat TypeDistFunctionAln(TypeAlignment aln);

TypeFloat computeMatchesAln(TypeAlignment aln);
TypeFloat computeNorm1Aln(TypeAlignment aln);
TypeFloat computeNorm2Aln(TypeAlignment aln);
void printMatchesAln(TypeAlignment aln);

TypeDistance computeWholeDistancePairAln(TypeAlignment aln, TypeDistFunctionAln *distfunc);

#endif
