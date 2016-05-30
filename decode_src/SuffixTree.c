#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "Sequences.h"
#include "Utils.h"
#include "SuffixTree.h"

//#define INFTY LONG_MAX
#define INFTY INT_MAX

#define LEFT 1
#define RIGHT 1 << 1
#define DONE 1 << 2
#define NEW 1 << 3


static TypePosition ident, progress;
static TypeNumber N;


static void addTransition(TypePosition s, TypePosition toadd, TypeSuffixTree *suffixTree);

static void canonize(TypePosition s, TypeNumber n, TypePosition k, TypePosition p, TypePosition *sres, TypePosition *kres, TypeSuffixTree *suffixTree);
static void update(TypePosition s, TypeNumber n, TypePosition k, TypePosition i, TypePosition *sres, TypePosition *kres, TypeSuffixTree *suffixTree);
static int test_and_split(TypePosition s, TypeNumber n, TypePosition k, TypePosition p, TypeSuffixTree *suffixTree, TypePosition *sres);
static TypePosition newSuffixNode(TypeNumber n, TypePosition start, TypePosition end, TypeSuffixTree *suffixTree);

static void printSuffixNode(FILE *f, TypePosition s, TypePosition depth, TypePosition de, TypeSuffixTree *suffixTree);
static void fixEnd(TypeSuffixTree *suffixTree);
static void canonizeEnd(TypePosition s, TypeNumber n, TypePosition *sres, TypePosition *kres, TypeSuffixTree *suffixTree);
static void fillInfoNode(TypePosition s, TypeSuffixTree *suffixTree);
static void fillNodeInfo(TypePosition s, TypeNumber n, TypeSuffixTree *suffixTree);

static TypeSuffixTree *computeSuffixTreeMUM(TypeSetOfSequences *set);

/*Compute local graph of the given order of positions of sequence seq
by computing the suffix tree of seq with Ukkonen's algorithm*/
TypeSuffixTree *computeSuffixTreeMUM(TypeSetOfSequences *set) {
	TypePosition i, k, tot, size, s;
	TypeSuffixTree *suffixTree;
	TypeNumber n;

	ident = 0;
	size = 0;
	for(n=0; n<set->number; n++)
		size += set->size[n]+1;
	suffixTree = (TypeSuffixTree*) monmalloc(sizeof(TypeSuffixTree));

	suffixTree->set = set;
	suffixTree->size = 0;
	suffixTree->nodes = (TypeSuffixNode*) monmalloc(2*size*sizeof(TypeSuffixNode));
	suffixTree->bottom = newSuffixNode(0, 0, -1, suffixTree);
	suffixTree->root = newSuffixNode(0, -1, -1, suffixTree);
	suffixTree->nodes[suffixTree->root].suffix = suffixTree->bottom;
	for(n=0; n<set->number; n++) {
		s = suffixTree->root;
		k = 0;
		for(i=0; i<=set->size[n]; i++) {
			update(s, n, k, i, &s, &k, suffixTree);
			canonize(s, n, k, i, &s, &k, suffixTree);
		}
	}
	fixEnd(suffixTree);
	suffixTree->nodes = (TypeSuffixNode*) monrealloc((void*)suffixTree->nodes, suffixTree->size*sizeof(TypeSuffixNode));
	return suffixTree;
}




/*Compute local graph of the given order of positions of sequence seq
by computing the suffix tree of seq with Ukkonen's algorithm
TypeSuffixTree *computeSuffixTree(TypeSetOfSequences *set) {
	TypePosition i, k, tot, size, code, s;
	TypeSuffixTree *suffixTree;
	TypeNumber n;
	
	ident = 0;
	size = 0;
	for(n=0; n<set->number; n++)
		size += set->size[n]+1;
	suffixTree = (TypeSuffixTree*) monmalloc(sizeof(TypeSuffixTree));
	suffixTree->set = set;
	suffixTree->model = estimateMarkovModel(set);
	suffixTree->size = 0;
	suffixTree->nodes = (TypeSuffixNode*) monmalloc(2*size*sizeof(TypeSuffixNode));
	suffixTree->bottom = newSuffixNode(0, 0, -1, suffixTree);
	suffixTree->root = newSuffixNode(0, -1, -1, suffixTree);
	suffixTree->nodes[suffixTree->root].suffix = suffixTree->bottom;
	suffixTree->cardCode = INC_SIZE_CODE;
	suffixTree->lengthCode = (TypePosition*) monmalloc(suffixTree->cardCode*sizeof(TypePosition));
	N = set->number;
progress = 0;
initProgress("suffix tree", size);
	for(n=0; n<set->number; n++) {
		s = suffixTree->root;
		k = 0;
		for(i=0; i<=set->size[n]; i++) {
			update(s, n, k, i, &s, &k, suffixTree);
			canonize(s, n, k, i, &s, &k, suffixTree);
updateProgress(progress++);
		}
		fillNodeInfo(suffixTree->root, n, suffixTree);
	}
	fixEnd(suffixTree);
flushProgress();
	suffixTree->nodes = (TypeSuffixNode*) monrealloc((void*)suffixTree->nodes, suffixTree->size*sizeof(TypeSuffixNode));

	suffixTree->cardCode = 0;
	code = 0;
progress = 0;
initProgress("code", suffixTree->size);
	buildCode(suffixTree->root, &code, 0, suffixTree);
flushProgress();
	suffixTree->cardCode = code;
	suffixTree->lengthCode = (TypePosition*) monrealloc(suffixTree->lengthCode, suffixTree->cardCode*sizeof(TypePosition));
	return suffixTree;
}
*/

TypeSuffixTree *getSuffixTree(TypeSetOfSequences *set) {
	TypePosition i, k, tot, size, s;
	TypeSuffixTree *suffixTree;
	TypeNumber n;

	size = 0;
	for(n=0; n<set->number; n++)
		size += set->size[n]+1;
	suffixTree = (TypeSuffixTree*) monmalloc(sizeof(TypeSuffixTree));
	suffixTree->set = set;
	suffixTree->size = 0;
	suffixTree->nodes = (TypeSuffixNode*) monmalloc(2*size*sizeof(TypeSuffixNode));
	suffixTree->bottom = newSuffixNode(0, 0, -1, suffixTree);
	suffixTree->root = newSuffixNode(0, -1, -1, suffixTree);
	suffixTree->nodes[suffixTree->root].suffix = suffixTree->bottom;
progress = 0;
initProgress("suffix tree", size);
	for(n=0; n<set->number; n++) {
		s = suffixTree->root;
		k = 0;
		for(i=0; i<=set->size[n]; i++) {
			update(s, n, k, i, &s, &k, suffixTree);
			canonize(s, n, k, i, &s, &k, suffixTree);
updateProgress(progress++);
		}
//		fillNodeInfo(suffixTree->root, n, suffixTree);
	}
	fixEnd(suffixTree);
	suffixTree->nodes = (TypeSuffixNode*) monrealloc((void*)suffixTree->nodes, suffixTree->size*sizeof(TypeSuffixNode));
flushProgress();
	return suffixTree;
}

void canonizeEnd(TypePosition s, TypeNumber n, TypePosition *sres, TypePosition *kres, TypeSuffixTree *suffixTree) {
	TypePosition s2;
	*sres = s;
	s2 = findTransition(*sres, suffixTree->set->sequence[n][*kres], suffixTree);
	while(*kres<=suffixTree->set->size[n] && suffixTree->nodes[s2].suffix >= 0) {
		*sres = s2;
		*kres += suffixTree->nodes[*sres].end-suffixTree->nodes[*sres].start+1;
		s2 = findTransition(*sres, suffixTree->set->sequence[n][*kres], suffixTree);
	}
}

void fixEnd(TypeSuffixTree *suffixTree) {
	TypeNumber n;
	for(n=0; n<suffixTree->set->number; n++) {
		TypePosition s, r;
		TypePosition k = 0, i = 0;
		canonizeEnd(suffixTree->root, n, &s, &k, suffixTree);
		r = findTransition(s, suffixTree->set->sequence[n][k], suffixTree);
		while(k<=suffixTree->set->size[n] && suffixTree->nodes[r].suffix < 0) {
			if(suffixTree->nodes[r].end == INFTY)
				suffixTree->nodes[r].end = suffixTree->set->size[n]-1;
			suffixTree->nodes[r].suffix = findTransition(suffixTree->nodes[s].suffix, suffixTree->set->sequence[n][k], suffixTree);
			canonizeEnd(suffixTree->nodes[s].suffix, n, &s, &k, suffixTree);
			r = findTransition(s, suffixTree->set->sequence[n][k], suffixTree);
		}
	}
}

void fillNodeInfo(TypePosition s, TypeNumber n, TypeSuffixTree *suffixTree) {
	TypePosition c;

	if(s<0)
		return;
	if(suffixTree->nodes[s].trans<0) {
		if(suffixTree->nodes[s].n == n)
			suffixTree->nodes[s].okk = 1;
		else
			suffixTree->nodes[s].okk = 0;
	} else {
		for(c=suffixTree->nodes[s].trans; c>=0; c=suffixTree->nodes[c].next)
			fillNodeInfo(c, n, suffixTree);
		if(suffixTree->nodes[s].okk >= 0)
			suffixTree->nodes[s].okk = 0;
		for(c=suffixTree->nodes[s].trans; c>=0 && suffixTree->nodes[s].okk >= 0 ; c=suffixTree->nodes[c].next) {
			if(suffixTree->nodes[c].okk<0)
				suffixTree->nodes[s].okk = -1;
			else {
				suffixTree->nodes[s].okk += suffixTree->nodes[c].okk;
				if(suffixTree->nodes[s].okk > 1)
					suffixTree->nodes[s].okk = -1;
			}
		}
	}
}

/*
void fillInfoNode(TypePosition s, TypeSuffixTree *suffixTree) {
	TypePosition c;
	TypeNumber n;
	if(s<0)
		return;
	if(suffixTree->nodes[s].trans<0) {
		suffixTree->info[s*suffixTree->sizeInfo+suffixTree->nodes[s].n] = 1;
		suffixTree->info[s*suffixTree->sizeInfo+suffixTree->sizeInfo-1] = 1;
	} else {
		for(c=suffixTree->nodes[s].trans; c>=0; c=suffixTree->nodes[c].next) {
			fillInfoNode(c, suffixTree);
			for(n=0; n<suffixTree->set->number; n++)
				suffixTree->info[s*suffixTree->sizeInfo+n] += suffixTree->info[c*suffixTree->sizeInfo+n];
		}
		suffixTree->info[s*suffixTree->sizeInfo+suffixTree->sizeInfo-1] = suffixTree->info[s*suffixTree->sizeInfo];
		for(n=1; n<suffixTree->set->number; n++) {
			if(suffixTree->info[s*suffixTree->sizeInfo+n]>suffixTree->info[s*suffixTree->sizeInfo+suffixTree->sizeInfo-1])
				suffixTree->info[s*suffixTree->sizeInfo+suffixTree->sizeInfo-1] = suffixTree->info[s*suffixTree->sizeInfo+n];
		}
	}
	if((suffixTree->info[s*suffixTree->sizeInfo+suffixTree->sizeInfo-1] <= 1 && suffixTree->nodes[s].okk == -1) || (suffixTree->info[s*suffixTree->sizeInfo+suffixTree->sizeInfo-1] > 1 && suffixTree->nodes[s].okk >= 0))
		printf("Problem node %d max = %d vs last = %d\n", s, suffixTree->info[s*suffixTree->sizeInfo+suffixTree->sizeInfo-1], suffixTree->nodes[s].okk);
}
*/


/*test whether or not a state with canonical reference pair (s(k,p)) is the end-point (has a t transition i positions later)
unchanged from Ukkonen's algorithm*/
int test_and_split(TypePosition s, TypeNumber n, TypePosition k, TypePosition p, TypeSuffixTree *suffixTree, TypePosition *sres) {
	if(k<=p) {
		TypePosition *prev, s2;
		for(prev = &(suffixTree->nodes[s].trans); *prev >= 0 && suffixTree->set->sequence[suffixTree->nodes[*prev].n][suffixTree->nodes[*prev].start] < suffixTree->set->sequence[n][k]; prev = &(suffixTree->nodes[*prev].next))
		;
		s2 = *prev;
		if(suffixTree->set->sequence[suffixTree->nodes[s2].n][suffixTree->nodes[s2].start+p-k+1] == suffixTree->set->sequence[n][p+1]) {
			*sres = s;
			return 1;
		} else {
			*prev = newSuffixNode(suffixTree->nodes[s2].n, suffixTree->nodes[s2].start, suffixTree->nodes[s2].start+p-k, suffixTree);
			suffixTree->nodes[*prev].next = suffixTree->nodes[s2].next;
			suffixTree->nodes[s2].next = -1;
			suffixTree->nodes[s2].start += p-k+1;
			suffixTree->nodes[*prev].trans = s2;
			*sres = *prev;
			return 0;
		}
	} else {
		*sres = s;
		return (findTransition(s, suffixTree->set->sequence[n][p+1], suffixTree) >= 0);
	}		
}

/*return the canonical reference (sres,(kres,p)) of the state represented by the reference pair (s,(k,p))
unchanged from Ukkonen's algorithm*/
void canonize(TypePosition s, TypeNumber n, TypePosition k, TypePosition p, TypePosition *sres, TypePosition *kres, TypeSuffixTree *suffixTree) {
	*sres = s;
	*kres = k;
//printf("canonize (%d, %d, %d)\n", s, k, p);
	if(p>=k) {
		TypePosition s2;
		s2 = findTransition(*sres, suffixTree->set->sequence[n][*kres], suffixTree);
		while(suffixTree->nodes[s2].end-suffixTree->nodes[s2].start <= p-*kres) {
//printf("(%d, %d)\n", *sres, *kres);
			*kres += suffixTree->nodes[s2].end-suffixTree->nodes[s2].start+1;
			*sres = s2;
			if(*kres<=p)
				s2 = findTransition(*sres, suffixTree->set->sequence[n][*kres], suffixTree);
		}
	}
//printf("fin\n");
}

/*update the suffix suffixTree by adding the ith caracter of the sequence.
(s,(k,i-1)) is the canonical reference pair of the active point,
(sres, (kres,i)) return a reference pair of the next active point
suffixTree stores basic information relative to the sequence and the suffix suffixTree
depth is the depth of the active point (the length of the string spelled out from the root to this one)
order is the given order of decoding
link stores the local graph of this order
only changes from Ukkonen's algorithm are lines starting with /'star' 'star'/ 
concerning update of depth and local graph by adding links*/	
void update(TypePosition s, TypeNumber n, TypePosition k, TypePosition i, TypePosition *sres, TypePosition *kres, TypeSuffixTree *suffixTree) {
	TypePosition oldr,  r;
	int end_point;
	oldr = -1;
	end_point = test_and_split(s, n, k, i-1, suffixTree, &r);
	while(!end_point) {
		addTransition(r, newSuffixNode(n, i, INFTY, suffixTree), suffixTree);
		if(oldr != -1)
			suffixTree->nodes[oldr].suffix = r;
		oldr = r;
		canonize(suffixTree->nodes[s].suffix, n, k, i-1, &s, &k, suffixTree);
		end_point = test_and_split(s, n, k, i-1, suffixTree, &r);
	}
	if(oldr >= 0)
		suffixTree->nodes[oldr].suffix = s;
	*sres = s; *kres = k;
}

/*Basic utilities*/

/*return the son of node with label starting by symbol if it exists, NULL otherwise*/
TypePosition findTransition(TypePosition s, TypeSymbol symbol, TypeSuffixTree *suffixTree) {
	TypePosition tmp;
	if(s == suffixTree->bottom)
		return suffixTree->root;
	for(tmp=suffixTree->nodes[s].trans; tmp >= 0 && suffixTree->set->sequence[suffixTree->nodes[tmp].n][suffixTree->nodes[tmp].start] < symbol; tmp=suffixTree->nodes[tmp].next)
	;
	if(tmp >= 0 && suffixTree->set->sequence[suffixTree->nodes[tmp].n][suffixTree->nodes[tmp].start] == symbol)
		return tmp;
	return -1;		
}

/*insert toadd as son of node in lexicographic order - assume node and toadd aren't NULL*/
void addTransition(TypePosition s, TypePosition toadd, TypeSuffixTree *suffixTree) {
	TypePosition cur;
	if(suffixTree->nodes[s].trans < 0 || suffixTree->set->sequence[suffixTree->nodes[suffixTree->nodes[s].trans].n][suffixTree->nodes[suffixTree->nodes[s].trans].start]>suffixTree->set->sequence[suffixTree->nodes[toadd].n][suffixTree->nodes[toadd].start]) {
		suffixTree->nodes[toadd].next = suffixTree->nodes[s].trans;
		suffixTree->nodes[s].trans = toadd;
		return;
	}
	for(cur = suffixTree->nodes[s].trans; suffixTree->nodes[cur].next >= 0 && suffixTree->set->sequence[suffixTree->nodes[suffixTree->nodes[cur].next].n][suffixTree->nodes[suffixTree->nodes[cur].next].start]<suffixTree->set->sequence[suffixTree->nodes[toadd].n][suffixTree->nodes[toadd].start]; cur=suffixTree->nodes[cur].next)
	;
	suffixTree->nodes[toadd].next = suffixTree->nodes[cur].next;
	suffixTree->nodes[cur].next = toadd;
}	


/*allocated and fill a new node*/
TypePosition newSuffixNode(TypeNumber n, TypePosition start, TypePosition end, TypeSuffixTree *suffixTree) {
	TypeNumber k;
	suffixTree->nodes[suffixTree->size].n = n;
	suffixTree->nodes[suffixTree->size].start = start;
	suffixTree->nodes[suffixTree->size].end = end;
	suffixTree->nodes[suffixTree->size].trans = -1;
	suffixTree->nodes[suffixTree->size].next = -1;
	suffixTree->nodes[suffixTree->size].suffix = -1;
	suffixTree->nodes[suffixTree->size].okk = 0;
/*for(k=0; k<suffixTree->sizeInfo; k++)
	suffixTree->info[suffixTree->size*suffixTree->sizeInfo+k] = 0;
*/	suffixTree->size++;
	return suffixTree->size-1;
}



void freeSuffixTree(TypeSuffixTree *suffixTree) {
	TypePosition i;
	monfree((void*)suffixTree->nodes);
	monfree(suffixTree);
}


/*print the suffix tree*/
void printSuffixTree(FILE *f, TypeSuffixTree *suffixTree) {
	printSuffixNode(f, suffixTree->root, 0, 0, suffixTree);
}


/*print recursively the suffix node*/
void printSuffixNode(FILE *f, TypePosition s, TypePosition depth, TypePosition de, TypeSuffixTree *suffixTree) {
	TypePosition d = depth;
	TypeSymbol t;
	TypePosition tmp;
	
	if(s < 0)
		return;
	if(depth>=0) {
		int i;
		TypePosition d = depth;
		/*Print the branches coming from higher nodes.*/
		for(d=0; d<depth; d++)
			fprintf(f, "|");
		fprintf(f, "+");
		fprintf(f, "(%d,%d,", suffixTree->nodes[s].n, suffixTree->nodes[s].start);
		if(suffixTree->nodes[s].end != INFTY)
			fprintf(f, "%d)", suffixTree->nodes[s].end);
		else
			fprintf(f, "infty)");
		if(suffixTree->nodes[s].suffix >= 0)
			fprintf(f, " %d>%d", s, suffixTree->nodes[s].suffix);
		else
			fprintf(f, " %d>?",s);
		fprintf(f, "(%d,%d)\n", suffixTree->nodes[s].okk,de);
	}
	for(tmp=suffixTree->nodes[s].trans; tmp >= 0; tmp = suffixTree->nodes[tmp].next)
		printSuffixNode(f, tmp, depth+1, de+suffixTree->nodes[tmp].end-suffixTree->nodes[tmp].start+1, suffixTree);
}

/*print all the sequences in a fasta file f*/
void printSequenceDebug(FILE *f, TypeOneSequence s, int sizeLine) {
		TypePosition position;
		for(position=0; position<s.size; position += sizeLine) {
			TypePosition j, maxLine;
			maxLine = position+sizeLine;
			maxLine = maxLine>s.size?s.size:maxLine;
			for(j=position; j<maxLine; j++)
				fprintf(f, " %ld", j);
			fprintf(f, "\n");
			for(j=position; j<maxLine; j++) {
				if(j>9)
					fprintf(f, " ");
				fprintf(f, " %ld", s.sequence[j]);
			}
			fprintf(f, "\n\n");
		}
		fprintf(f, "\n\n");
}


