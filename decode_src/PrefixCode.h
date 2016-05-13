#ifndef PrefixCodeF
#define PrefixCodeF

#include <stdlib.h>
#include <stdio.h>
#include "Sequences.h"
#include "SuffixTree.h"
#include "MarkovEmbed.h"

typedef struct CODESCHEME {
	TypeSuffixTree *suffixTree;
	TypePosition cardCode, *lengthCode, *code, buffSize;
} TypeCodeScheme;


TypeSetOfSequences *computeCoded(TypeSetOfSequences *set);
void encode(TypeSetOfSequences *out, TypeCodeScheme *scheme);
void printDecodedSequencesFasta(FILE *f, TypeSetOfSequences *set, TypeSetOfSequences *dec, int sizeLine);
void printDecodedSequencesFastaAlt(FILE *f, TypeSetOfSequences *set, TypeSetOfSequences *dec, int sizeLine);
TypeSetOfSequences *computeAntecedent(TypeSetOfSequences *dec, TypeSuffixTree *suffixTree);
void setThreshold(double thr);
TypeCodeScheme *getCodeThreshold(TypeSuffixTree *suffixTree, double threshold, TypeMarkovModel *model);
TypeCodeScheme *getCodeLength(TypeSuffixTree *suffixTree, TypePosition length);
void freeScheme(TypeCodeScheme *scheme);
void buildCodeLength(TypePosition l, TypePosition s, TypePosition depth, TypeCodeScheme *scheme);
void buildCodeThreshold(double t, TypePosition s, TypePosition depth, double init, TypeMarkovModel *model, TypeCodeScheme *scheme);
void printLengthDistribution(FILE *f, TypePosition *length, TypePosition size);

#endif
