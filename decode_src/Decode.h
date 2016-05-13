#ifndef DecodeF
#define DecodeF

#include <stdlib.h>
#include <stdio.h>
#include "Sequences.h"

typedef struct SUFFIXNODE {
	struct SUFFIXNODE *trans, *next, *suffix;
	TypePosition start, end;
} TypeSuffixNode;

typedef struct SUFFIXTREE {
	TypeSuffixNode *root, *bottom;
	TypeSymbol *sequence;
} TypeSuffixTree;

typedef struct LINK {
	struct LINK *next;
	TypePosition position;
} TypeLink;

typedef struct LOCALNODE {
	struct LOCALNODE *trans, *next, *ancestor, *nlink, *ptrans, *pnext;
	int state;
	TypePosition position;
} TypeLocalNode;

typedef struct LOCALTREE {
	TypeLocalNode **leaf;
	TypePosition size;
} TypeLocalTree;

TypeSetOfSequences localDecode(TypeSetOfSequences set,  TypePosition order);
TypeOneSequence localDecodeOneSequence(TypeOneSequence seq, TypePosition order);
TypeSetOfSequences slideDecode(TypeSetOfSequences set,  TypePosition order);
TypeOneSequence slideDecodeOneSequence(TypeOneSequence seq, TypePosition order);

void printDecodedSequencesFasta(FILE *f, TypeSetOfSequences set, TypeSetOfSequences dec, int sizeLine);

TypeOneSequence computeSequence(TypeLocalTree *local);
TypeSuffixTree *computeSuffixTree(TypeOneSequence seq);
TypeLocalTree *initLocalTree(TypeSuffixTree *suffixTree, TypeOneSequence seq, TypePosition order);
TypeLocalTree *computeLocalTree(TypeSuffixTree *suffixTree, TypeOneSequence seq, TypePosition order);
void freeSuffixTree(TypeSuffixTree *suffixTree);
void freeLocalTree(TypeLocalTree *localTree);
void printSuffixTree(FILE *f, TypeSuffixTree *suffixTree);
void printLocalTree(FILE *f, TypeLocalTree *localTree);
void printSequenceDebug(FILE *f, TypeOneSequence s, int sizeLine);

#endif
