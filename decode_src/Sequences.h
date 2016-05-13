
#ifndef SequencesF
#define SequencesF

#include <stdlib.h>
#include <stdio.h>

typedef int TypePosition;
typedef int TypeSymbol;
typedef short TypeNumber;
/*Structure storing a set of sequence*/
typedef struct AMBIGUITY {
	char *car;
	TypeSymbol **set, *size, number;
} TypeAmbiguity;

/*Structure storing one sequence*/
typedef struct ONESEQUENCE {
/*cardinal of the alphabet*/
	TypePosition size;
	TypeSymbol cardinal;
	TypeSymbol *sequence;
	TypeAmbiguity ambiguity;
} TypeOneSequence;

/*Structure storing a set of sequence*/
typedef struct SEQUENCES {
	char **name, *table, *title;
/*number of TypeSetOfSequences include in the set*/
	TypeNumber number;
/*cardinal of the alphabet*/
	TypeSymbol cardinal;
/*size[i] size of the sequence i*/
	TypePosition *size, totSize;
/*table of all the sequences*/
	TypeSymbol **sequence;
	TypeAmbiguity ambiguity;
} TypeSetOfSequences;

typedef struct REFERENCE {
	TypePosition pos;
	TypeNumber n;
} TypeRef;

double meanLength(TypeSetOfSequences set);
TypePosition maxLength(TypeSetOfSequences set);
TypePosition minLength(TypeSetOfSequences set);
TypePosition totalLength(TypeSetOfSequences set);

void fixSequencesAmbiguity(TypeSetOfSequences *set);

TypeOneSequence toOne(TypeSetOfSequences set);
void printSequencesFasta(FILE *f, TypeSetOfSequences s, int sizeLine);
TypeSetOfSequences readSequencesFasta(FILE *f, char *table, int canInc);
TypeSetOfSequences readSequencesRaw(FILE *f);
void writeSequences(FILE *f, TypeSetOfSequences seq);
void freeSequences(TypeSetOfSequences seq);
TypeAmbiguity getXNAAmbiguity();
TypeAmbiguity getProteinAmbiguity();
TypeSetOfSequences clone(TypeSetOfSequences s);
int comparePosition(void const *a, void const *b);

#define SEP -1
#define FIN -1
#define DNA "ACGTYRMKWSBDHVN" /*"ACGT" + "YRMKWSBDHVN"*/
#define RNA "ACGUYRMKWSBDHVN" /*"ACGU" + "YRMKWSBDHVN"*/
#define PRO "DEGNQCSTYAVLIPFMWKRHBZX" /*"DEGNQCSTYAVLIPFMWKRH" + "X"*/

#endif
