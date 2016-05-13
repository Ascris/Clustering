#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <gsl/gsl_cdf.h>
#include "Sequences.h"
#include "Utils.h"
#include "PrefixCode.h"

#define CHAR_END -10
#define THRESHOLD 0.00000000001

static TypePosition progress;
typedef struct STACK_ELT {
	TypePosition symb, pos;
} TypeStackElt;
static double threshold = THRESHOLD;


static void setCode(TypePosition s, TypePosition code, TypeCodeScheme *scheme);
static double updateProb(TypeMarkovModel *model, double init, TypeSymbol *w, TypePosition size);



void setThreshold(double thr) {
	threshold = thr;
}

void printLengthDistribution(FILE *f, TypePosition *length, TypePosition size) {
	TypePosition min, max, pos, *occ, tot = 0;

	min = length[0];
	max = length[0];
	for(pos=0; pos<size; pos++) {
		if(length[pos]<min)
			min = length[pos];
		if(length[pos]>max)
			max = length[pos];
	}
	occ = (TypePosition*) monmalloc((max-min+1)*sizeof(TypePosition));
	for(pos=0; pos<=max-min; pos++)
		occ[pos] = 0;
	for(pos=0; pos<size; pos++) {
		occ[length[pos]-min]++;
		tot++;
	}
	for(pos=0; pos<=max-min; pos++)
		fprintf(f, "%ld\t%ld\t%lE\n", pos+min, occ[pos], ((double)occ[pos])/((double)tot));
	monfree((void*)occ);
}
	

/*fail if scheme->code[c]>=0*/
void encode(TypeSetOfSequences *out, TypeCodeScheme *scheme) {
	TypeNumber n;
	TypePosition i, j, tot, code;
	for(n=0; n<scheme->suffixTree->set->number; n++) {
		TypePosition s, c;
		j = 0;
		s = scheme->suffixTree->root;
		c = findTransition(s, scheme->suffixTree->set->sequence[n][j], scheme->suffixTree);
//printf("pos %d-%d code %d (depth %d)\n", n, i, suffixTree->nodes[c].code, j);
		for(i=0; i<scheme->suffixTree->set->size[n]; i++) {
			while(c >= 0 && scheme->code[c]<0) {
//printf(" c code %d_%d (depth %d)\n", c, suffixTree->nodes[c].code, j);
				s = c;
				j += scheme->suffixTree->nodes[s].end-scheme->suffixTree->nodes[s].start+1;
				c = findTransition(s, scheme->suffixTree->set->sequence[n][j], scheme->suffixTree);
			}
if(c<0)
	exitProg(ExitOk,"oooops");
//printf("pos %d-%d code %d_%d (depth %d)\n", n, i, c, suffixTree->nodes[c].code, j);
			out->sequence[n][i] = scheme->code[c];
			s = scheme->suffixTree->nodes[s].suffix;
			c = findTransition(s, scheme->suffixTree->set->sequence[n][j], scheme->suffixTree);
		}
	}
}


/*
TypeSetOfSequences *computeCoded(TypeSetOfSequences *set) {
	TypeSetOfSequences *res, *dec;
	TypeNumber n;
	TypeCodeScheme scheme;

	res = (TypeSetOfSequences*) monmalloc(sizeof(TypeSetOfSequences));
	res->name = set->name;
	res->table = set->table;
	res->number = set->number;
	res->size = set->size;
	res->sequence = (TypeSymbol**) monmalloc(res->number*sizeof(TypeSymbol*));
	for(n=0; n<res->number; n++)
		res->sequence[n] = (TypeSymbol*) monmalloc((res->size[n]+1)*sizeof(TypeSymbol));
	scheme.suffixTree = computeSuffixTree(set);
	res->cardinal = suffixTree->cardCode;
printf("Encode\n");
	encode(res, suffixTree);
printf("Decode\n");
	dec = computeAntecedent(res, suffixTree);
	freeSuffixTree(suffixTree);
	for(n=0; n<res->number; n++)
		monfree((void*) res->sequence[n]);
	monfree((void*) res->sequence);
	monfree((void*) res);
	return dec;
}
*/

TypeCodeScheme *getCodeThreshold(TypeSuffixTree *suffixTree, double threshold, TypeMarkovModel *model) {
	setThreshold(threshold);
	double *prob;
	TypeCodeScheme *scheme;

//suffixTree->model = estimateMarkovModel(set);

	scheme = (TypeCodeScheme*) monmalloc(sizeof(TypeCodeScheme));
	scheme->suffixTree = suffixTree;
	scheme->code = (TypePosition*) monmalloc(suffixTree->size*sizeof(TypePosition));
	scheme->code[scheme->suffixTree->bottom] = -1;
	scheme->code[scheme->suffixTree->root] = -1;
	scheme->buffSize = INC_SIZE_CODE;
	scheme->lengthCode = (TypePosition*) monmalloc(scheme->buffSize*sizeof(TypePosition));
	scheme->cardCode = 0;
progress = 0;
initProgress("code", scheme->suffixTree->size);
	buildCodeThreshold(threshold, scheme->suffixTree->root, 0, 1., model, scheme);
flushProgress();
	scheme->lengthCode = (TypePosition*) monrealloc(scheme->lengthCode, scheme->cardCode*sizeof(TypePosition));
}

TypeCodeScheme *getCodeLength(TypeSuffixTree *suffixTree, TypePosition length) {
	TypeCodeScheme *scheme;

	scheme = (TypeCodeScheme*) monmalloc(sizeof(TypeCodeScheme));
	scheme->suffixTree = suffixTree;
	scheme->code = (TypePosition*) monmalloc(suffixTree->size*sizeof(TypePosition));
	scheme->code[scheme->suffixTree->bottom] = -1;
	scheme->code[scheme->suffixTree->root] = -1;
	scheme->buffSize = INC_SIZE_CODE;
	scheme->lengthCode = (TypePosition*) monmalloc(scheme->buffSize*sizeof(TypePosition));
	scheme->cardCode = 0;
progress = 0;
initProgress("code", suffixTree->size);
	buildCodeLength(length, suffixTree->root, 0, scheme);
flushProgress();
	scheme->lengthCode = (TypePosition*) monrealloc(scheme->lengthCode, scheme->cardCode*sizeof(TypePosition));
	return scheme;
}


void buildCodeLength(TypePosition l, TypePosition s, TypePosition depth, TypeCodeScheme *scheme) {
	TypePosition c;
	if(s<0)
		return;
	updateProgress(++progress);
	if(scheme->suffixTree->nodes[s].trans<0 || (depth>=l)) {
		scheme->code[s] = scheme->cardCode++;
		if(scheme->code[s]>=scheme->buffSize) {
			scheme->buffSize += INC_SIZE_CODE;
			scheme->lengthCode = (TypePosition*) monrealloc(scheme->lengthCode, scheme->buffSize*sizeof(TypePosition));
		}
		scheme->lengthCode[scheme->code[s]] =  l;
		for(c=scheme->suffixTree->nodes[s].trans; c>=0; c=scheme->suffixTree->nodes[c].next)
			setCode(c, scheme->cardCode, scheme);
	} else {
		scheme->code[s] = -1;
		for(c=scheme->suffixTree->nodes[s].trans; c>=0; c=scheme->suffixTree->nodes[c].next)
			buildCodeLength(l, c, depth+scheme->suffixTree->nodes[c].end-scheme->suffixTree->nodes[c].start+1, scheme);
	}
}

double updateProb(TypeMarkovModel *model, double init, TypeSymbol *w, TypePosition size) {
	TypePosition position;
	double res = init;
	for(position=0; position<size; position++)
		res *= model->trans[w[position-1]][w[position]];
	return res;
}

void buildCodeThreshold(double t, TypePosition s, TypePosition depth, double init, TypeMarkovModel *model, TypeCodeScheme *scheme) {
	TypePosition c;
	double prob;
	if(s<0)
		return;
	updateProgress(++progress);
//printf("init %lE, %ld %ld\n", init, scheme->suffixTree->nodes[s].end-scheme->suffixTree->nodes[s].start+1, depth);
	if(scheme->suffixTree->nodes[s].trans>=0) {
		if(depth>1)
			prob = updateProb(model, init, &(scheme->suffixTree->set->sequence[scheme->suffixTree->nodes[s].n][scheme->suffixTree->nodes[s].start]), scheme->suffixTree->nodes[s].end-scheme->suffixTree->nodes[s].start+1);
		else {
			if(depth<1)
				prob = 1.;
			else
				prob = model->init[scheme->suffixTree->set->sequence[scheme->suffixTree->nodes[s].n][scheme->suffixTree->nodes[s].start]];
		}
	}
//getProbApp(&(scheme->suffixTree->set->sequence[scheme->suffixTree->nodes[s].n][scheme->suffixTree->nodes[s].end-depth+1]), depth, 2, scheme->suffixTree->set, model);
//printf("probit\t%lE\t%ld\n\n", prob, scheme->suffixTree->set->totSize);
//if(prob<=1)
//printf("depth %ld\tprob %lE vs %lE\n", depth, getProbApp(&(scheme->suffixTree->set->sequence[scheme->suffixTree->nodes[s].n][scheme->suffixTree->nodes[s].end-depth+1]), depth, 2, scheme->suffixTree->set, model), gsl_cdf_binomial_Q(1, prob, scheme->suffixTree->set->totSize-scheme->suffixTree->set->number*(depth-1)));
//	if(suffixTree->nodes[s].trans<0 || ((getProbApp(&(suffixTree->set->sequence[suffixTree->nodes[s].n][suffixTree->nodes[s].end-depth+1]), depth, 2, suffixTree->set, suffixTree->model)<t))) {
	if(scheme->suffixTree->nodes[s].trans<0 || ((gsl_cdf_binomial_Q(1, prob, scheme->suffixTree->set->totSize-scheme->suffixTree->set->number*(depth-1))<t))) {
		scheme->code[s] = scheme->cardCode++;
		if(scheme->code[s]>=scheme->buffSize) {
			scheme->buffSize += INC_SIZE_CODE;
			scheme->lengthCode = (TypePosition*) monrealloc(scheme->lengthCode, scheme->buffSize*sizeof(TypePosition));
		}
//		scheme->lengthCode[scheme->code[s]] =  depth-scheme->suffixTree->nodes[s].end+scheme->suffixTree->nodes[s].start-1;
		scheme->lengthCode[scheme->code[s]] =  depth-1;
		if(scheme->lengthCode[scheme->code[s]]<1)
			scheme->lengthCode[scheme->code[s]] = 1;
		for(c=scheme->suffixTree->nodes[s].trans; c>=0; c=scheme->suffixTree->nodes[c].next)
			setCode(c, scheme->cardCode, scheme);
	} else {
		scheme->code[s] = -1;
		for(c=scheme->suffixTree->nodes[s].trans; c>=0; c=scheme->suffixTree->nodes[c].next)
			buildCodeThreshold(t, c, depth+scheme->suffixTree->nodes[c].end-scheme->suffixTree->nodes[c].start+1, prob, model, scheme);
	}
}

void freeScheme(TypeCodeScheme *scheme) {
	monfree((void*)scheme->code);
	monfree((void*)scheme->lengthCode);
	freeSuffixTree(scheme->suffixTree);
	monfree((void*)scheme);
}

/*
void buildCode(TypePosition s, TypePosition *code, TypePosition depth, TypeSuffixTree *suffixTree) {
	TypePosition c;
	int m = 10, M = 15, d = 12;
	TypePosition sumc = -1;
	if(s<0)
		return;
	updateProgress(++progress);
//	if(s->max == 1) {
//	if(s->max == 1 && s->sum<=m && (depth-s->end+s->start) >= d) {
//	if(((depth-s->end+s->start) >= d || s->trans == NULL)) {
//	if((suffixTree->info[s].max == 1 && suffixTree->info[s].sum<=m && (depth) >= d) || suffixTree->nodes[s].trans<0) {
//	if((suffixTree->info[s].max == 1 && (depth-suffixTree->nodes[s].end+suffixTree->nodes[s].start) >= d) || suffixTree->nodes[s].trans<0) {
//	if(suffixTree->nodes[s].trans<0 || (suffixTree->info[s].max == 1 && (depth > M || (depth>m && getProbMat(&(suffixTree->set->sequence[suffixTree->nodes[s].n][suffixTree->nodes[s].end-depth+1]), depth, 2, suffixTree->set, suffixTree->model)<threshold)))) {
//	if(suffixTree->info[s].max == 1 || suffixTree->nodes[s].trans<0) {
//	if(suffixTree->nodes[s].trans<0 || (suffixTree->info[s].max == 1 && suffixTree->info[s].sum<=m && (depth) >= d)) {
//	if(suffixTree->nodes[s].trans<0 || (suffixTree->info[s].max == 1 &&  (depth) >= d)) {
//if(suffixTree->nodes[s].trans<0 || (suffixTree->info[s].max == 1 && (depth > M || (depth>m && getProbMat(&(suffixTree->set->sequence[suffixTree->nodes[s].n][suffixTree->nodes[s].end-depth+1]), depth, 2, suffixTree->set, suffixTree->model)<threshold)))) {
//if(suffixTree->nodes[s].trans<0 || (suffixTree->nodes[s].okk >= 0 && (depth > M || (depth>m && getProbApp(&(suffixTree->set->sequence[suffixTree->nodes[s].n][suffixTree->nodes[s].end-depth+1]), depth, 2, suffixTree->set, suffixTree->model)<threshold)))) {
//if(suffixTree->nodes[s].trans<0 || (suffixTree->nodes[s].okk >= 0 && (getProbApp(&(suffixTree->set->sequence[suffixTree->nodes[s].n][suffixTree->nodes[s].end-depth+1]), depth, 2, suffixTree->set, suffixTree->model)<threshold))) {
if(suffixTree->nodes[s].trans<0 || ((getProbApp(&(suffixTree->set->sequence[suffixTree->nodes[s].n][suffixTree->nodes[s].end-depth+1]), depth, 2, suffixTree->set, suffixTree->model)<threshold))) {
/*
if(suffixTree->nodes[s].trans>=0) {
printf("depth %d %.2E/%2E\n", depth, getProbApp(&(suffixTree->set->sequence[suffixTree->nodes[s].n][suffixTree->nodes[s].end-depth+1]), depth, 2, suffixTree->set, suffixTree->model), threshold);
if(depth<9){
TypePosition x;
for(x=0; x<depth; x++)
printf("%d ", suffixTree->set->sequence[suffixTree->nodes[s].n][suffixTree->nodes[s].end-depth+1+x]);
printf("\n");
}
}

		suffixTree->nodes[s].code = (*code)++;
		if(suffixTree->nodes[s].code>=suffixTree->cardCode) {
			suffixTree->cardCode += INC_SIZE_CODE;
			suffixTree->lengthCode = (TypePosition*) monrealloc(suffixTree->lengthCode, suffixTree->cardCode*sizeof(TypePosition));
		}
//		suffixTree->lengthCode[suffixTree->nodes[s].code] = depth;
//		suffixTree->lengthCode[suffixTree->nodes[s].code] = MIN(depth, MAX(d, depth-suffixTree->nodes[s].end+suffixTree->nodes[s].start));
		suffixTree->lengthCode[suffixTree->nodes[s].code] =  depth-suffixTree->nodes[s].end+suffixTree->nodes[s].start-1;
		if(suffixTree->lengthCode[suffixTree->nodes[s].code]<1)
			suffixTree->lengthCode[suffixTree->nodes[s].code] = 1;
//printf("l%d = %d\n", suffixTree->nodes[s].code, suffixTree->lengthCode[suffixTree->nodes[s].code]);
		for(c=suffixTree->nodes[s].trans; c>=0; c=suffixTree->nodes[c].next)
			setCode(c, *code, suffixTree);
	} else
		for(c=suffixTree->nodes[s].trans; c>=0; c=suffixTree->nodes[c].next)
			buildCode(c, code, depth+suffixTree->nodes[c].end-suffixTree->nodes[c].start+1, suffixTree);
}
*/

void setCode(TypePosition s, TypePosition code, TypeCodeScheme *scheme) {
	TypePosition c;
	if(s<0)
		return;
	scheme->code[s] = code; updateProgress(++progress);
	for(c=scheme->suffixTree->nodes[s].trans; c>=0; c=scheme->suffixTree->nodes[c].next)
		setCode(c, code, scheme);
}


/*print a Fasta format file including all the sequences of the set - each decoded symbol is designed by 'X<number>'
where X is the symbol in the same position in the non-decoded sequence and <number> an identifiant number*/
void printDecodedSequencesFasta(FILE *f, TypeSetOfSequences *set, TypeSetOfSequences *dec, int sizeLine) {
	TypeNumber n;
	TypeSymbol *nIdent, *sIdent, s;
	nIdent = (TypeSymbol*) monmalloc(dec->cardinal*sizeof(TypeSymbol));
	for(s=0; s<dec->cardinal; s++)
		nIdent[s] = SPECIAL;
	sIdent = (TypeSymbol*) monmalloc((set->cardinal+set->ambiguity.number)*sizeof(TypeSymbol));
	for(s=0; s<set->cardinal; s++)
		sIdent[s] = 0;
	for(n=0; n<set->number; n++) {
		TypePosition position;
		if(set->name != NULL)
			fprintf(f, ">%s\n", set->name[n]);
		else
			fprintf(f, ">Sequence %d\n", n+1);
		for(position=0; position<dec->size[n]; position += sizeLine) {
			TypePosition j, maxLine;
			maxLine = position+sizeLine;
			maxLine = maxLine>dec->size[n]?dec->size[n]:maxLine;
			if(set->table != NULL) {
				for(j=position; j<maxLine; j++) {
//if(dec->sequence[n][j]>=dec->cardinal || set->sequence[n][j] >= set->cardinal)
//printf("%d:%d/%d d %d/%d i %d/%d nident = %d\n", n, j, dec->size[n], dec->sequence[n][j], dec->cardinal, set->sequence[n][j], set->cardinal, nIdent[dec->sequence[n][j]]);
					if(nIdent[dec->sequence[n][j]] == SPECIAL)
						nIdent[dec->sequence[n][j]] = sIdent[set->sequence[n][j]]++;
					fprintf(f, "%c%ld ", set->table[set->sequence[n][j]], nIdent[dec->sequence[n][j]]);
				}
			} else {
				for(j=position; j<maxLine; j++) {
					if(nIdent[dec->sequence[n][j]] == SPECIAL)
						nIdent[dec->sequence[n][j]] = sIdent[set->sequence[n][j]]++;
					fprintf(f, "%ld-%ld ", set->sequence[n][j], nIdent[set->sequence[n][j]]);
				}
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\n\n");
	}
	monfree((void*) nIdent);
	monfree((void*) sIdent);
}
/*
/*Renvoie l'antecedent maximal a l'ordre N de "seq"*/
/*TypeSetOfSequences *computeAntecedent(TypeSetOfSequences *dec, TypeSuffixTree *suffixTree) {
	TypeRef **next, *first;
	TypeNumber n;
	TypePosition symb, i, tot, card, last, **proj, **far;
	TypeStackElt *stack;
	TypeSetOfSequences *res;

	first = (TypeRef*) monmalloc(suffixTree->cardCode*sizeof(TypeRef));
	far = (TypePosition**) monmalloc(dec->number*sizeof(TypePosition*));
	for(n=0; n<dec->number; n++)
		far[n] = (TypePosition*) monmalloc(dec->size[n]*sizeof(TypePosition));
	for(n=0; n < dec->number; n++) {
		TypePosition best, last;
		best = 0; last = 0;
		while(last<dec->size[n]) {
			for(;best<dec->size[n] && best+suffixTree->lengthCode[dec->sequence[n][best]]-1<last; best++)
			;
			for(; last<dec->size[n] && last<(best+suffixTree->lengthCode[dec->sequence[n][best]]); last++)
				far[n][last] = best;
		}
	}
	next = (TypeRef**) monmalloc(dec->number*sizeof(TypeRef*));
	for(n=0; n<dec->number; n++)
		next[n] = (TypeRef*) monmalloc(dec->size[n]*sizeof(TypeRef));
	for(symb=0; symb<suffixTree->cardCode; symb++)
		first[symb].n = CHAR_END;
	for(n=dec->number-1; n >= 0; n--) {
		for(i=dec->size[n]-1; i>=0; i--) {
			next[n][i] = first[dec->sequence[n][i]];
			first[dec->sequence[n][i]].n = n;
			first[dec->sequence[n][i]].pos = i;
		}
	}
	tot = 0;
	for(symb=0; symb<suffixTree->cardCode; symb++) {
		tot += suffixTree->lengthCode[symb];
	}
	stack = (TypeStackElt*) monmalloc(tot*sizeof(TypeStackElt));
	proj = (TypePosition**) monmalloc(suffixTree->cardCode*sizeof(TypePosition*));
	for(symb=0; symb<suffixTree->cardCode; symb++) {
		proj[symb] = (TypePosition*) monmalloc(suffixTree->lengthCode[symb]*sizeof(TypePosition));
		for(i=0; i<suffixTree->lengthCode[symb]; i++)
			proj[symb][i] = -1;
	}
	card = 0;
	last=-1;
	for(symb=0; symb<suffixTree->cardCode; symb++) {
		for(i=0; i<suffixTree->lengthCode[symb]; i++) {
//printf("\ns %d i %d/%d p %d card %d f %d\n", symb, i, suffixTree->lengthCode[symb], proj[symb][i], card, first[symb].n);
			if(proj[symb][i] < 0) {
				stack[++last].symb = symb;
				stack[last].pos = i;
				proj[symb][i] = card;
				while(last >= 0) {
					TypeStackElt cur;
					TypeRef ref;
					cur = stack[last--];
					ref = first[cur.symb];
					while(ref.n != CHAR_END) {
						TypePosition j, k;
						for(j=far[ref.n][ref.pos+cur.pos]; j<dec->size[ref.n] && j<=ref.pos+cur.pos; j++) {
							if((j+suffixTree->lengthCode[dec->sequence[ref.n][j]]) > ref.pos+cur.pos) {
								if(proj[dec->sequence[ref.n][j]][ref.pos+cur.pos-j] == -1) {
//printf("cur %d:%d last %d j %d (%d - %d) stack (%d,%d)\n", ref.n, ref.pos, last, j, ref.pos+cur.pos, far[ref.n][ref.pos+cur.pos], dec->sequence[ref.n][j], ref.pos+cur.pos-j);
									stack[++last].symb = dec->sequence[ref.n][j];
									stack[last].pos = ref.pos+cur.pos-j;
									proj[dec->sequence[ref.n][j]][ref.pos+cur.pos-j] = card;
								} else {
									if(proj[dec->sequence[ref.n][j]][ref.pos+cur.pos-j] != card) {
										printf("Erreur A\nref %d:%d symb %d seq %d card %d proj %d l %d/%d\n", ref.n, ref.pos, cur.symb, dec->sequence[ref.n][j], card, proj[dec->sequence[ref.n][j]][ref.pos+cur.pos-j], ref.pos+cur.pos-j, suffixTree->lengthCode[dec->sequence[ref.n][j]]);
printf("ref %d:%d cur %d j %d (%d - %d) stack (%d,%d)\n", ref.n, ref.pos, cur.pos, j, ref.pos+cur.pos, far[ref.n][ref.pos+cur.pos], dec->sequence[ref.n][j], ref.pos+cur.pos-j);
										exitProg(ExitOk,"Bye");
									}
								}
							}
						}
						ref = next[ref.n][ref.pos];
					}
				}
				card++;
			}
		}
	}
	monfree((void*) stack);
	monfree((void*)first);
	for(n=0; n<dec->number; n++)
		monfree((void*)next[n]);
	monfree((void*)next);
	for(n=0; n<dec->number; n++)
		monfree((void*)far[n]);
	monfree((void*)far);
	res = (TypeSetOfSequences*) monmalloc(sizeof(TypeSetOfSequences));
	res->name = dec->name;
	res->table = dec->table;
	res->number = dec->number;
	res->size = (TypePosition*) monmalloc(res->number*sizeof(TypePosition));
	res->sequence = (TypeSymbol**) monmalloc(res->number*sizeof(TypeSymbol*));
	for(n=0; n<res->number; n++) {
		res->size[n] = dec->size[n];
		res->sequence[n] = (TypeSymbol*) monmalloc((res->size[n])*sizeof(TypeSymbol));
	}
	res->cardinal = card;
	for(n=0; n<res->number; n++)
		for(i=0; i<dec->size[n]; i++)
			res->sequence[n][i] = proj[dec->sequence[n][i]][0];
	for(symb=0; symb<suffixTree->cardCode; symb++)
		monfree((void*)proj[symb]);
	monfree((void*)proj);
	return res;
}
*/

/*print a Fasta format file including all the sequences of the set - each decoded symbol is designed by 'X<number>'
where X is the symbol in the same position in the non-decoded sequence and <number> an identifiant number*/
void printDecodedSequencesFastaAlt(FILE *f, TypeSetOfSequences *set, TypeSetOfSequences *dec, int sizeLine) {
	TypeNumber n;
	TypeSymbol *nIdent, *sIdent, s;
	TypePosition *occ;

	occ = (TypePosition*) monmalloc(dec->cardinal*sizeof(TypePosition));
	for(s=0; s<dec->cardinal; s++)
		occ[s] = 0;
	for(n=0; n<set->number; n++) {
		TypePosition position;
		for(position=0; position<dec->size[n]; position ++)
			occ[dec->sequence[n][position]]++;
	}
	nIdent = (TypeSymbol*) monmalloc(dec->cardinal*sizeof(TypeSymbol));
	for(s=0; s<dec->cardinal; s++)
		nIdent[s] = SPECIAL;
	sIdent = (TypeSymbol*) monmalloc((set->cardinal+set->ambiguity.number)*sizeof(TypeSymbol));
	for(s=0; s<set->cardinal; s++)
		sIdent[s] = 0;
	for(n=0; n<set->number; n++) {
		TypePosition position;
		if(set->name != NULL)
			fprintf(f, ">%s\n", set->name[n]);
		else
			fprintf(f, ">Sequence %d\n", n+1);
		for(position=0; position<dec->size[n]; position += sizeLine) {
			TypePosition j, maxLine;
			maxLine = position+sizeLine;
			maxLine = maxLine>dec->size[n]?dec->size[n]:maxLine;
			if(set->table != NULL) {
				for(j=position; j<maxLine; j++) {
//if(dec->sequence[n][j]>=dec->cardinal || set->sequence[n][j] >= set->cardinal)
//printf("%d:%d/%d d %d/%d i %d/%d nident = %d\n", n, j, dec->size[n], dec->sequence[n][j], dec->cardinal, set->sequence[n][j], set->cardinal, nIdent[dec->sequence[n][j]]);
					if(occ[dec->sequence[n][j]]>1) {
						if(nIdent[dec->sequence[n][j]] == SPECIAL)
							nIdent[dec->sequence[n][j]] = sIdent[set->sequence[n][j]]++;
						fprintf(f, "%c%ld ", set->table[set->sequence[n][j]], nIdent[dec->sequence[n][j]]);
					} else
						fprintf(f, "%c ", tolower(set->table[set->sequence[n][j]]));

				}
			} else {
				for(j=position; j<maxLine; j++) {
					if(nIdent[dec->sequence[n][j]] == SPECIAL)
						nIdent[dec->sequence[n][j]] = sIdent[set->sequence[n][j]]++;
					fprintf(f, "%ld-%ld ", set->sequence[n][j], nIdent[set->sequence[n][j]]);
				}
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\n\n");
	}
	monfree((void*) nIdent);
	monfree((void*) sIdent);
}


