#ifndef MarkovEmbedF
#define MarkovEmbedF

#include <stdlib.h>
#include <stdio.h>
#include "Sequences.h"

typedef struct MARKOV_MODEL {
	TypeSymbol cardinal;
	double *init, **trans;
	TypeAmbiguity ambiguity;
} TypeMarkovModel;

typedef struct MARKOV_EMBED {
	TypePosition size, no, **bord;
	TypeMarkovModel model;
	TypeSymbol *word;
	double **state, sink;
} TypeMarkovEmbed;

TypeMarkovModel *estimateMarkovModel(TypeSetOfSequences *set);
TypePosition **getNext(TypeSymbol *w, TypePosition size, TypeSymbol cardinal);
double getProb(TypeSymbol *w, TypePosition size, TypePosition nocc, TypeSetOfSequences *set, TypeMarkovModel *m);
double getProbMat(TypeSymbol *w, TypePosition size, TypePosition nocc, TypeSetOfSequences *set, TypeMarkovModel *m);
void getBounds(TypePosition *start, TypePosition *end, TypePosition nocc, TypeSetOfSequences *set, TypePosition *sorted, double threshold, TypeMarkovModel *m);
double getProbApp(TypeSymbol *w, TypePosition size, TypePosition nocc, TypeSetOfSequences *set, TypeMarkovModel *m);
void freeModel(TypeMarkovModel *model);

#endif
