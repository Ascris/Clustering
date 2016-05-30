#ifndef SuffixTreeF
#define SuffixTreeF

#include <stdlib.h>
#include <stdio.h>
#include "Sequences.h"
#include "MarkovEmbed.h"

#define INC_SIZE_CODE 100

typedef struct SUFFIXNODE {
	TypePosition start, end, trans, next, suffix;
	TypeNumber n;
	char okk;
} TypeSuffixNode;

typedef struct INFO {
	TypePosition *occ, max, sum;
} TypeInfo;

typedef struct SUFFIXTREE {
	TypeSuffixNode *nodes;
	TypeSetOfSequences *set;
	TypePosition root, bottom, size;/*, sizeInfo, *info;*/
} TypeSuffixTree;

TypeSuffixTree *computeSuffixTree(TypeSetOfSequences *set);
void freeSuffixTree(TypeSuffixTree *suffixTree);
void printSuffixTree(FILE *f, TypeSuffixTree *suffixTree);
TypeSuffixTree *getSuffixTree(TypeSetOfSequences *set);
TypePosition findTransition(TypePosition s, TypeSymbol symbol, TypeSuffixTree *suffixTree);




#endif
