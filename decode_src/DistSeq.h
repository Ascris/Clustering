#ifndef DistSeqF
#define DistSeqF

#include <stdlib.h>
#include <stdio.h>
#include "Sequences.h"
#include "Distance.h"

typedef TypeFloat TypeDistFunction(TypeSetOfSequences*);

TypeFloat computeAlex(TypeSetOfSequences *set);
TypeDistance computeWholeDistancePair(TypeSetOfSequences set, TypeDistFunction *distfunc);
TypeDistance computeWholeDistanceGlobal(TypeSetOfSequences *set);
TypeDistance computeWholeDistanceGlobalBis(TypeSetOfSequences set, TypeDistFunction *distfunc);
TypeDistance computeMSMDistance(TypeSetOfSequences *set);
TypeDistance computeWholeDistanceDec(TypeSetOfSequences *dec);

#endif
