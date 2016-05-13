#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "Utils.h"
#include "SuffixTree.h"
#include "PrefixDecode.h"
#include "DistSeq.h"
#include "MSM.h"


#define SPECIAL -1
static TypePosition fillOcc(TypeSymbol *seq, TypePosition size, TypeSymbol card, TypePosition *occ);
static double getDist(TypePosition *occ0, TypePosition l0, TypePosition *occ1, TypePosition l1, TypeSymbol card);



TypePosition fillOcc(TypeSymbol *seq, TypePosition size, TypeSymbol card, TypePosition *occ) {
	TypePosition p, l=0;
	TypeSymbol c;
	for(c=0; c<card; c++)
		occ[c] = 0;
	for(p=0; p<size; p++)
		if(seq[p]<card) {
			l++;
			occ[seq[p]]++;
		}
	return l;
}

double getDist(TypePosition *occ0, TypePosition l0, TypePosition *occ1, TypePosition l1, TypeSymbol card) {
	double dist = 0.;
	TypeSymbol c;
	for(c=0; c<card; c++)
		dist += MIN(occ0[c], occ1[c]);
	return 1.-dist/((double) MIN(l0, l1));
}

TypeFloat computeAlex(TypeSetOfSequences *set) {
	TypePosition *occ0, *occ1, p, l0=0, l1=0, dist = 0, count = 0,
	*nocc, minocc, maxocc, tmp;
	TypeSymbol c;
	
	occ0 =(TypePosition*) monmalloc(set->cardinal*sizeof(TypePosition));
	occ1 =(TypePosition*) monmalloc(set->cardinal*sizeof(TypePosition));
	for(c=0; c<set->cardinal; c++) {
		occ0[c] = 0;
		occ1[c] = 0;
	}
	for(p=0; p<set->size[0]; p++)
		if(set->sequence[0][p]<set->cardinal) {
			l0++;
			occ0[set->sequence[0][p]]++;
		}
	for(p=0; p<set->size[1]; p++)
		if(set->sequence[1][p]<set->cardinal) {
			l1++;
			occ1[set->sequence[1][p]]++;
		}
	for(c=0; c<set->cardinal; c++) {

		dist += MIN(occ0[c], occ1[c]);
	}
	monfree((void*)occ0);
	monfree((void*)occ1);
	return 1.-dist/((double) MIN(l0, l1));
}




TypeDistance computeMSMDistance(TypeSetOfSequences *set) {
	TypeDistance dist;
	TypeNumber i, j;
	TypeSymbol c;
	TypePosition p, *occi, *occj;
	TypeSetOfSequences stmp;

	dist.number = set->number;
	dist.name = set->name;	
	dist.table = (TypeFloat*) monmalloc(((dist.number*(dist.number-1))/2)*sizeof(TypeFloat));

	stmp.number = 2;
	stmp.size = (TypePosition*) monmalloc(2*sizeof(TypePosition));;
	stmp.sequence = (TypeSymbol **) monmalloc(2*sizeof(TypeSymbol*));
	stmp.cardinal = set->cardinal;
	stmp.ambiguity = set->ambiguity;
	stmp.table = set->table;
	occi =(TypePosition*) monmalloc(stmp.cardinal*sizeof(TypePosition));
	occj =(TypePosition*) monmalloc(stmp.cardinal*sizeof(TypePosition));
	for(i=1; i<dist.number; i++) {
		TypePosition p, li;
		stmp.sequence[0] = set->sequence[i]; stmp.size[0] = set->size[i];
		for(c=0; c<stmp.cardinal; c++)
			occi[c] = 0;
		li = 0;
		for(p=0; p<stmp.size[0]; p++)
			if(stmp.sequence[0][p]<stmp.cardinal) {
				li++;
				occi[stmp.sequence[0][p]]++;
			}
		for(j=0; j<i; j++) {
			TypePosition lj, lmin;
			long ind = (i*(i-1))/2+j;
			double pmatch;
			stmp.sequence[1] = set->sequence[j]; stmp.size[1] = set->size[j];
			for(c=0; c<stmp.cardinal; c++)
				occj[c] = 0;
			lj = 0;
			for(p=0; p<stmp.size[1]; p++)
				if(stmp.sequence[1][p]<stmp.cardinal) {
					lj++;
					occj[stmp.sequence[1][p]]++;
				}
			pmatch = 0.;
			for(c=0; c<stmp.cardinal; c++)
				pmatch += ((double) occi[c]*occj[c])/((double)li*lj);
			lmin = 1 + ((TypePosition) (-log((double) ((1.-pmatch)*li*lj))/log(pmatch)+.5));
			dist.table[ind] = getMSMDist(&stmp, lmin);
		}
	}
	monfree((void*)occi);
	monfree((void*)occj);
	monfree((void*)stmp.size);
	monfree((void*)stmp.sequence);
	return dist;
}

/*
TypeDistance computeWholeDistancePair(TypeSetOfSequences set, TypeDistFunction *distfunc) {
	TypeDistance dist;
	TypeNumber i, j;
	TypeSetOfSequences stmp, *dec;

//printf("****************set.cardinal = %ld\n", set.cardinal);
//exitProg(0, NULL);

	dist.number = set.number;
	dist.name = set.name;	
	dist.table = (TypeFloat*) monmalloc(((dist.number*(dist.number-1))/2)*sizeof(TypeFloat));

	stmp.number = 2;
	stmp.size = (TypePosition*) monmalloc(2*sizeof(TypePosition));;
	stmp.sequence = (TypeSymbol **) monmalloc(2*sizeof(TypeSymbol*));
	stmp.cardinal = set.cardinal;
	stmp.ambiguity = dec->ambiguity;
	stmp.table = dec->table;

//printf("order %ld - %ld\n", orderStart, orderEnd);

	for(i=1; i<dist.number; i++) {
printf("\n\n%s %d/%d\n", dist.name[i], i+1, dist.number);
		for(j=0; j<i; j++) {
			TypeNumber ind = (i*(i-1))/2+j;
//printf("%d/%d\n", j, i);
			stmp.sequence[0] = set.sequence[j]; stmp.size[0] = set.size[j];
			stmp.sequence[1] = set.sequence[i]; stmp.size[1] = set.size[i];
			dec = computeCoded(&stmp);
			dist.table[ind] = distfunc(dec);
			monfree((void*)dec->sequence[0]);
			monfree((void*)dec->sequence[1]);
			monfree((void*)dec->sequence);
			monfree((void*)dec->size);
			monfree((void*)dec);
		}
	}
	monfree((void*)stmp.size);
	monfree((void*)stmp.sequence);
	return dist;
}

TypeDistance computeWholeDistanceGlobal(TypeSetOfSequences *set) {
	TypeDistance dist;
	TypeNumber i, j;
	long ind ;
	TypeSetOfSequences stmp, *dec;
	TypeSuffixTree *suffixTree;
	TypeOneSequence one;
	TypePosition o, p, *occi, *occj, li, lj;
	TypeNumber n;
	TypeSymbol c, *listi, *listj, ni, nj;


	dist.number = set->number;
	dist.name = set->name;
	dist.table = (TypeFloat*) monmalloc(((dist.number*(dist.number-1))/2)*sizeof(TypeFloat));
	for(ind=((dist.number*(dist.number-1))/2)-1; ind>=0; ind--)
		dist.table[ind] = POS_INFTY;

	dec = set;
	listi =(TypePosition*) monmalloc(dec->cardinal*sizeof(TypeSymbol));
	listj =(TypePosition*) monmalloc(dec->cardinal*sizeof(TypeSymbol));
	occi =(TypePosition*) monmalloc(dec->cardinal*sizeof(TypePosition));
	occj =(TypePosition*) monmalloc(dec->cardinal*sizeof(TypePosition));
	for(c=0; c<dec->cardinal; c++) {
		occi[c] = 0;
		occj[c] = 0;
	}
	initProgress("distance", (dist.number*(dist.number-1))/2);
	for(i=1; i<dist.number; i++) {
		ni=0;
		li=0;
		for(p=0; p<dec->size[i]; p++)
			if(dec->sequence[i][p]<dec->cardinal) {
				if(occi[dec->sequence[i][p]] == 0)
					listi[ni++] = dec->sequence[i][p];
				li++;
				occi[dec->sequence[i][p]]++;
			}
		for(j=0; j<i; j++) {
			long ind = (i*(i-1))/2+j;
			TypePosition sum = 0;
			nj=0;
			lj=0;
			for(p=0; p<dec->size[j]; p++)
				if(dec->sequence[j][p]<dec->cardinal) {
					if(occj[dec->sequence[j][p]] == 0)
						listj[nj++] = dec->sequence[j][p];
					lj++;
					occj[dec->sequence[j][p]]++;
				}
			for(c=0; c<ni; c++) {
				sum +=  MIN(occi[listi[c]], occj[listi[c]]);
				occj[listi[c]] = 0;
			}
			for(c=0; c<nj; c++) {
				if(occj[listj[c]]) {
					sum +=  MIN(occi[listj[c]], occj[listj[c]]);
					occj[listj[c]] = 0;
				}
			}
			dist.table[ind] = 1.-((double)sum)/((double) MIN(li, lj));
			updateProgress(ind);
		}
		for(c=0; c<ni; c++)
			occi[listi[c]] = 0;
	}
	printf("\n");
	monfree((void*) occi);
	monfree((void*) occj);
	monfree((void*) listi);
	monfree((void*) listj);
	for(n=0; n<dec->number; n++)
		monfree((void*) dec->sequence[n]);
	monfree((void*) dec->sequence);
	monfree((void*) dec->size);
	monfree((void*) dec);
	return dist;
}

*/

TypeDistance computeWholeDistanceDec(TypeSetOfSequences *dec) {
	TypeDistance dist;
	TypeNumber i, j;
	long ind ;
	TypeSetOfSequences stmp;
	TypeSuffixTree *suffixTree;
	TypeOneSequence one;
	TypePosition o, p, *occi, *occj, li, lj;
	TypeNumber n;
	TypeSymbol c, *listi, *listj, ni, nj;


	dist.number = dec->number;
	dist.name = dec->name;
	dist.table = (TypeFloat*) monmalloc(((dist.number*(dist.number-1))/2)*sizeof(TypeFloat));
	for(ind=((dist.number*(dist.number-1))/2)-1; ind>=0; ind--)
		dist.table[ind] = POS_INFTY;

	listi =(TypePosition*) monmalloc(dec->cardinal*sizeof(TypeSymbol));
	listj =(TypePosition*) monmalloc(dec->cardinal*sizeof(TypeSymbol));
	occi =(TypePosition*) monmalloc(dec->cardinal*sizeof(TypePosition));
	occj =(TypePosition*) monmalloc(dec->cardinal*sizeof(TypePosition));
	for(c=0; c<dec->cardinal; c++) {
		occi[c] = 0;
		occj[c] = 0;
	}
	initProgress("distance", (dist.number*(dist.number-1))/2);
	for(i=1; i<dist.number; i++) {
		ni=0;
		li=0;
		for(p=0; p<dec->size[i]; p++)
			if(dec->sequence[i][p]<dec->cardinal) {
				if(occi[dec->sequence[i][p]] == 0)
					listi[ni++] = dec->sequence[i][p];
				li++;
				occi[dec->sequence[i][p]]++;
			}
		for(j=0; j<i; j++) {
			long ind = (i*(i-1))/2+j;
			TypePosition sum = 0;
			nj=0;
			lj=0;
			for(p=0; p<dec->size[j]; p++)
				if(dec->sequence[j][p]<dec->cardinal) {
					if(occj[dec->sequence[j][p]] == 0)
						listj[nj++] = dec->sequence[j][p];
					lj++;
					occj[dec->sequence[j][p]]++;
				}
			for(c=0; c<ni; c++) {
				sum +=  MIN(occi[listi[c]], occj[listi[c]]);
				occj[listi[c]] = 0;
			}
			for(c=0; c<nj; c++) {
				if(occj[listj[c]]) {
					sum +=  MIN(occi[listj[c]], occj[listj[c]]);
					occj[listj[c]] = 0;
				}
			}
			dist.table[ind] = 1.-((double)sum)/((double) MIN(li, lj));
			updateProgress(ind);
		}
		for(c=0; c<ni; c++)
			occi[listi[c]] = 0;
	}
	monfree((void*) occi);
	monfree((void*) occj);
	monfree((void*) listi);
	monfree((void*) listj);
	return dist;
}











/*
TypeDistance computeWholeDistanceGlobal(TypeSetOfSequences *set) {
	TypeDistance dist;
	TypeNumber i, j;
	long ind ;
	TypeSetOfSequences stmp, *dec;
	TypeSuffixTree *suffixTree;
	TypeOneSequence one;
	TypePosition o, **occ, *l;
	TypeNumber n;
	TypeSymbol c;


	dist.number = set->number;
	dist.name = set->name;
	dist.table = (TypeFloat*) monmalloc(((dist.number*(dist.number-1))/2)*sizeof(TypeFloat));
	for(ind=((dist.number*(dist.number-1))/2)-1; ind>=0; ind--)
		dist.table[ind] = POS_INFTY;

	dec = computeCoded(set);
	l = (TypePosition*) monmalloc(dec->number*sizeof(TypePosition));
	occ = (TypePosition**) monmalloc(dec->number*sizeof(TypePosition*));
	initProgress("distance", (dist.number));
	for(i=0; i<dist.number; i++) {
		occ[i] =(TypePosition*) monmalloc(dec->cardinal*sizeof(TypePosition));
		l[i] = fillOcc(dec->sequence[i], dec->size[i], dec->cardinal, occ[i]);
		updateProgress(i);
	}
	printf("\n");
	initProgress("distance", (dist.number*(dist.number-1))/2);
	for(i=1; i<dist.number; i++) {
		for(j=0; j<i; j++) {
			long ind = (i*(i-1))/2+j;
			dist.table[ind] = getDist(occ[i], l[i], occ[j], l[j], dec->cardinal);
			updateProgress(ind);
		}
	}
	printf("\n");
	for(n=0; n<dec->number; n++) {
		monfree((void*) occ[n]);
		monfree((void*) dec->sequence[n]);
	}
	monfree((void*) occ);
	monfree((void*) l);
	monfree((void*) dec->sequence);
	monfree((void*) dec->size);
	monfree((void*) dec);
	return dist;
}

TypeDistance computeWholeDistanceGlobalBis(TypeSetOfSequences set, TypeDistFunction *distfunc) {
	TypeDistance dist;
	TypeNumber i, j;
	long ind ;
	TypeSetOfSequences stmp, *dec;
	TypeSuffixTree *suffixTree;
	TypeOneSequence one;
	TypePosition o;
	TypeNumber n;

	dist.number = set.number;
	dist.name = set.name;	
	dist.table = (TypeFloat*) monmalloc(((dist.number*(dist.number-1))/2)*sizeof(TypeFloat));
	for(ind=((dist.number*(dist.number-1))/2)-1; ind>=0; ind--)
		dist.table[ind] = POS_INFTY;
	
	dec = computeCoded(&set);

	stmp.number = 2;
	stmp.size = (TypePosition*) monmalloc(2*sizeof(TypePosition));;
	stmp.sequence = (TypeSymbol **) monmalloc(2*sizeof(TypeSymbol*));
	stmp.cardinal = dec->cardinal;
//	stmp.ambiguity = dec->ambiguity;
//	stmp.table = dec->table;

	initProgress("distance", (dist.number*(dist.number-1))/2);
	for(i=1; i<dist.number; i++) {
//printf("\n\n%s %d/%d\n", dist.name[i], i+1, dist.number);
		for(j=0; j<i; j++) {
			long ind = (i*(i-1))/2+j;
//printf("%d/%d\n", j, i);
			updateProgress(ind);
			stmp.sequence[0] = dec->sequence[j]; stmp.size[0] = dec->size[j];
			stmp.sequence[1] = dec->sequence[i]; stmp.size[1] = dec->size[i];
			dist.table[ind] = distfunc(&stmp);
		}
	}
	printf("\n");
	monfree((void*)stmp.size);
	monfree((void*)stmp.sequence);
	for(n=0; n<dec->number; n++)
		monfree((void*) dec->sequence[n]);
	monfree((void*) dec->sequence);
	monfree((void*) dec->size);
	monfree((void*) dec);
	return dist;
}

*/
