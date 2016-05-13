#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "Sequences.h"
#include "Utils.h"
#include "SuffixTree.h"
#include "MSM.h"

#define SIZE_BUFFER_MUM 1000

typedef struct MUM {
	TypePosition start0, start1, size;
} TypeMUM;


static void computeMUM(TypePosition s, TypePosition min, TypePosition depth, TypePosition **pos, TypePosition *size, TypeMUM **res, TypePosition *sm, TypePosition *bm, TypeSuffixTree *suffixTree);


double getMSMDist(TypeSetOfSequences *set, TypePosition min) {
	TypeSuffixTree *suffixTree;
	TypeMUM *res;
	TypePosition sm, bm, *(pos[2]), size[2], i, sum;
	set->number = 2;
	pos[0] = (TypePosition*) monmalloc(set->size[0]*sizeof(TypePosition));
	pos[1] = (TypePosition*) monmalloc(set->size[1]*sizeof(TypePosition));
	bm = SIZE_BUFFER_MUM; sm = 0;
	res = (TypeMUM *) monmalloc(bm*sizeof(TypeMUM));
	suffixTree = getSuffixTree(set);
	computeMUM(suffixTree->root, min, 0, pos, size, &res, &sm, &bm, suffixTree);
	freeSuffixTree(suffixTree);
	monfree((void*)pos[0]);
	monfree((void*)pos[1]);
	sum = 0;
	for(i=0; i<sm; i++)
		sum += res[i].size;
	monfree((void*)res);
	return 1.-((double) sum)/((double) MIN(set->size[0], set->size[1]));
}

void getMUM(TypeSetOfSequences *set, TypePosition min) {
	TypeSuffixTree *suffixTree;
	TypeMUM *res;
	TypePosition sm, bm, *(pos[2]), size[2], i;
	set->number = 2;
	pos[0] = (TypePosition*) monmalloc(set->size[0]*sizeof(TypePosition));
	pos[1] = (TypePosition*) monmalloc(set->size[1]*sizeof(TypePosition));
	bm = SIZE_BUFFER_MUM; sm = 0;
	res = (TypeMUM *) monmalloc(bm*sizeof(TypeMUM));
	suffixTree = getSuffixTree(set);
	computeMUM(suffixTree->root, min, 0, pos, size, &res, &sm, &bm, suffixTree);
	freeSuffixTree(suffixTree);
	monfree((void*)pos[0]);
	monfree((void*)pos[1]);
	res = (TypeMUM *) monrealloc((void*) res, bm*sizeof(TypeMUM));
	for(i=0; i<sm; i++)
		printf("%d\t%d\t%d\n", res[i].start0, res[i].start1, res[i].size);
}

void computeMUM(TypePosition s, TypePosition min, TypePosition depth, TypePosition **pos, TypePosition *size, TypeMUM **res, TypePosition *sm, TypePosition *bm, TypeSuffixTree *suffixTree) {
	TypePosition c;
	if(s<0)
		return;
	depth += suffixTree->nodes[s].end-suffixTree->nodes[s].start+1;
	if(depth >= min) {
		if(suffixTree->nodes[s].trans<0) {
			pos[suffixTree->nodes[s].n][size[suffixTree->nodes[s].n]++] = suffixTree->set->size[suffixTree->nodes[s].n]-depth;
		} else {
			TypePosition si[2], sc[2];
 			si[0] = size[0];
			si[1] = size[1];
			c=suffixTree->nodes[s].trans;
			computeMUM(c, min, depth, pos, size, res, sm, bm, suffixTree);
			for(c=suffixTree->nodes[c].next; c>=0; c=suffixTree->nodes[c].next) {
				TypePosition i, j;
				sc[0] = size[0];
				sc[1] = size[1];
				computeMUM(c, min, depth, pos, size, res, sm, bm, suffixTree);
				for(i=si[0]; i<sc[0]; i++)
					for(j=sc[1]; j<size[1]; j++) {
						if((pos[0][i]>0 && pos[1][j]>0 && suffixTree->set->sequence[0][pos[0][i]-1] != suffixTree->set->sequence[1][pos[1][j]-1]) ) {
							if(*sm >= *bm) {
								*bm += SIZE_BUFFER_MUM;
								*res = (TypeMUM *) monrealloc((void*) *res, *bm*sizeof(TypeMUM));
							}
							(*res)[*sm].start0 = pos[0][i];
							(*res)[*sm].start1 = pos[1][j];
							(*res)[*sm].size = depth;
							(*sm)++;
						}
					}
				for(j=si[1]; j<sc[1]; j++)
					for(i=sc[0]; i<size[0]; i++) {
						if((pos[0][i]>0 && pos[1][j]>0 && suffixTree->set->sequence[0][pos[0][i]-1] != suffixTree->set->sequence[1][pos[1][j]-1]) ) {
							if(*sm >= *bm) {
								*bm += SIZE_BUFFER_MUM;
								*res = (TypeMUM *) monrealloc((void*) *res, *bm*sizeof(TypeMUM));
							}
							(*res)[*sm].start0 = pos[0][i];
							(*res)[*sm].start1 = pos[1][j];
							(*res)[*sm].size = depth;
							(*sm)++;
						}
					}
			}
		}
	} else {
		for(c=suffixTree->nodes[s].trans; c>=0; c=suffixTree->nodes[c].next) {
			size[0] = 0;
			size[1] = 0;
			computeMUM(c, min, depth, pos, size, res, sm, bm, suffixTree);
		}
	}
}

