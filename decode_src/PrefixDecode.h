#ifndef PrefixDecodeF
#define PrefixDecodeF

#include <stdlib.h>
#include <stdio.h>
#include "Sequences.h"
#include "SuffixTree.h"
#include "PrefixCode.h"

typedef struct ATOMNODE {
	TypePosition parent, child, sibling, f, rchild, rnext, rprec;
	int state;
} TypeAtomNode;

typedef struct ATOMTREE {
	TypeAtomNode *nodes;
	TypePosition size;
} TypeAtomTree;

TypeSetOfSequences *computeAtomTree(TypePosition *sortCode, TypePosition *sizeCode, TypePosition prefsize, TypeSetOfSequences *set);
TypePosition *getTableCode(TypePosition *length, TypePosition size);
TypeSetOfSequences *getDecodedFromScheme(TypeCodeScheme *scheme);
TypeSetOfSequences *getDecodedFromThreshold(TypeSetOfSequences *set, double t, TypeMarkovModel *model);
TypeSetOfSequences *getDecodedFromLength(TypeSetOfSequences *set, TypePosition l);


#endif