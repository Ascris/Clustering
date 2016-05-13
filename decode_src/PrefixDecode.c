#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "Utils.h"
#include "PrefixDecode.h"
#include "PrefixCode.h"

#define bitNew 0x01
#define bitSpe 0x02
#define bitSta 0x04
#define bitIni 0x08
#define bitDon 0x10
#define bitLas 0x20

static void setCode(TypePosition code, TypePosition n, TypeAtomTree *atomTree, TypePosition *tmp);
static void initAtomNode(TypeAtomNode *n);

TypeSetOfSequences *getDecodedFromThreshold(TypeSetOfSequences *set, double t, TypeMarkovModel *model) {
	TypeSetOfSequences *dec;
	TypeSuffixTree *suffixTree;
	TypeCodeScheme *scheme;

	suffixTree = getSuffixTree(set);
	scheme = getCodeThreshold(suffixTree, t, model);
	dec = getDecodedFromScheme(scheme);
	freeScheme(scheme);
	return dec;
}

TypeSetOfSequences *getDecodedFromLength(TypeSetOfSequences *set, TypePosition l) {
	TypeSetOfSequences *dec;
	TypeSuffixTree *suffixTree;
	TypeCodeScheme *scheme;

	suffixTree = getSuffixTree(set);
	scheme = getCodeLength(suffixTree, l);
	dec = getDecodedFromScheme(scheme);
	freeScheme(scheme);
	return dec;
}

TypeSetOfSequences *getDecodedFromScheme(TypeCodeScheme *scheme) {
	TypeSetOfSequences *res, *dec;
	TypeNumber n;
	TypeSetOfSequences *set;
	TypePosition *tab;
	set = scheme->suffixTree->set;
	res = (TypeSetOfSequences*) monmalloc(sizeof(TypeSetOfSequences));
	res->name = set->name;
	res->table = set->table;
	res->number = set->number;
	res->size = set->size;
	res->sequence = (TypeSymbol**) monmalloc(res->number*sizeof(TypeSymbol*));
	for(n=0; n<res->number; n++)
		res->sequence[n] = (TypeSymbol*) monmalloc((res->size[n]+1)*sizeof(TypeSymbol));
	res->cardinal = scheme->cardCode;
	encode(res, scheme);
	tab = getTableCode(scheme->lengthCode, scheme->cardCode);
//printf("min %ld max %ld\n", scheme->lengthCode[tab[scheme->cardCode-1]], scheme->lengthCode[tab[0]]);
	dec = computeAtomTree(tab, scheme->lengthCode, scheme->cardCode, res);
	monfree((void*) tab);
	for(n=0; n<res->number; n++)
		monfree((void*) res->sequence[n]);
	monfree((void*) res->sequence);
	monfree((void*) res);
	return dec;
}

TypeSetOfSequences *computeAtomTree(TypePosition *sortCode, TypePosition *sizeCode, TypePosition prefsize, TypeSetOfSequences *set) {
	TypeAtomTree *atomTree;
	TypePosition p, k, i, size, sizeBis, *offset, *prec, *last, *stack, sizeSta, c, *special, sizeSpe, *new, *tmp, sizeNew, *not, sizeNot;
	TypeNumber n;
	TypeSetOfSequences *res;

	size = 0;
	for(n=0; n<set->number; n++)
		size += set->size[n];
	atomTree = (TypeAtomTree*) monmalloc(sizeof(TypeAtomTree));
	atomTree->nodes = (TypeAtomNode*) monmalloc(2*size*sizeof(TypeAtomNode));
	atomTree->size = 0;
	last = (TypePosition*) monmalloc(prefsize*sizeof(TypePosition));
	prec = (TypePosition*) monmalloc(size*sizeof(TypePosition));
	for(c=0; c<prefsize; c++)
		last[c] = -1;
	for(n=0; n<set->number; n++) {
		initAtomNode(&(atomTree->nodes[atomTree->size]));
		prec[atomTree->size] = last[set->sequence[n][0]];
		last[set->sequence[n][0]] = atomTree->size;
		atomTree->size++;
		for(p=1; p<set->size[n]; p++) {
			initAtomNode(&(atomTree->nodes[atomTree->size]));
			prec[atomTree->size] = last[set->sequence[n][p]];
			last[set->sequence[n][p]] = atomTree->size;
			atomTree->nodes[atomTree->size].rchild = atomTree->size-1;
			atomTree->nodes[atomTree->size-1].f = atomTree->size;
			atomTree->size++;
		}
	}
	sizeBis = size/2;
	stack = (TypePosition*) monmalloc(sizeBis*sizeof(TypePosition));
	sizeSta = 0;
	special = (TypePosition*) monmalloc(sizeBis*sizeof(TypePosition));
	sizeSpe = 0;
	not = (TypePosition*) monmalloc(sizeBis*sizeof(TypePosition));
	sizeNot = 0;
	new = (TypePosition*) monmalloc(sizeBis*sizeof(TypePosition));
	sizeNew = 0;
	c = 0;
	k = sizeCode[sortCode[c]];
initProgress("decode", sizeCode[sortCode[0]]);
	while(k>=1) {
		TypePosition stop, x;
		k = sizeCode[sortCode[c]];
		do {
			TypePosition p, d;
			initAtomNode(&(atomTree->nodes[atomTree->size]));
			for(p=last[sortCode[c]]; p>=0; p=prec[p]) {
				for(d=p; !(atomTree->nodes[d].state&bitIni) && (atomTree->nodes[d].parent>=0); d=atomTree->nodes[d].parent)
					atomTree->nodes[d].state |= bitIni;
				if(!(atomTree->nodes[d].state&bitIni)) {
					atomTree->nodes[d].state |= bitIni;
					atomTree->nodes[d].parent = atomTree->size;
					atomTree->nodes[d].sibling = atomTree->nodes[atomTree->size].child;
					atomTree->nodes[atomTree->size].child = d;
				}
			}
			if(atomTree->nodes[atomTree->size].child>=0) {
				if(atomTree->nodes[atomTree->nodes[atomTree->size].child].sibling>=0) {
					for(d=atomTree->nodes[atomTree->size].child; d>=0; d=atomTree->nodes[d].sibling) {
						TypePosition z;
						for(z=atomTree->nodes[d].rchild; z>=0;) {
							TypePosition y = atomTree->nodes[z].rnext;
							if(atomTree->nodes[z].parent<0 || (atomTree->nodes[atomTree->nodes[z].parent].state&bitNew)) {
								atomTree->nodes[z].f = atomTree->size;
								atomTree->nodes[z].rnext = atomTree->nodes[atomTree->size].rchild;
								atomTree->nodes[z].rprec = -1;
								if(atomTree->nodes[atomTree->size].rchild >= 0)
									atomTree->nodes[atomTree->nodes[atomTree->size].rchild].rprec = z;
								atomTree->nodes[atomTree->size].rchild = z;
							} else {
								if(atomTree->nodes[atomTree->nodes[z].parent].f<0) {
									atomTree->nodes[atomTree->nodes[z].parent].f = atomTree->size;
									atomTree->nodes[atomTree->nodes[z].parent].rnext = atomTree->nodes[atomTree->size].rchild;
									atomTree->nodes[atomTree->nodes[z].parent].rprec = -1;
									if(atomTree->nodes[atomTree->size].rchild >= 0)
										atomTree->nodes[atomTree->nodes[atomTree->size].rchild].rprec = atomTree->nodes[z].parent;
									atomTree->nodes[atomTree->size].rchild = atomTree->nodes[z].parent;
								}
							}
							z = y;
						}
					}
					atomTree->nodes[atomTree->size].state |= bitNew;
					new[sizeNew++] = atomTree->size++;
				} else {
					atomTree->nodes[atomTree->nodes[atomTree->size].child].parent = -1;
				}
			}
			c++;
		} while(c<prefsize && sizeCode[sortCode[c]] == sizeCode[sortCode[c-1]]);
		if(c<prefsize)
			stop = sizeCode[sortCode[c]];
		else
			stop = 0;
		do {
			TypePosition d, x;
updateProgress(sizeCode[sortCode[0]]-k);
//printf("%ld/%ld\n", sizeCode[sortCode[0]]-k, sizeCode[sortCode[0]]);
			while(sizeSpe>0) {
				d = special[--sizeSpe];
				atomTree->nodes[d].state &= ~bitSpe;
				if(atomTree->nodes[d].f < 0) {
					initAtomNode(&(atomTree->nodes[atomTree->size]));
					atomTree->nodes[atomTree->size].state |= bitNew;
					new[sizeNew++] = atomTree->size;
					atomTree->nodes[d].f = atomTree->size;
					atomTree->nodes[d].rnext = atomTree->nodes[atomTree->size].rchild;
					if(atomTree->nodes[atomTree->size].rchild >= 0)
						atomTree->nodes[atomTree->nodes[atomTree->size].rchild].rprec = d;
					atomTree->nodes[atomTree->size].rchild = d;
					stack[sizeSta++] = d;
					while(sizeSta>0) {
						TypePosition a,b;
						a = stack[--sizeSta];
						for(b=atomTree->nodes[a].child; b>=0; b=atomTree->nodes[b].sibling) {
							if(atomTree->nodes[b].f>=0 && !(atomTree->nodes[atomTree->nodes[b].f].state & bitDon)) {
								TypePosition z;
								atomTree->nodes[atomTree->nodes[b].f].state |= bitDon;
								atomTree->nodes[atomTree->nodes[b].f].parent = atomTree->size;
								atomTree->nodes[atomTree->nodes[b].f].sibling = atomTree->nodes[atomTree->size].child;
								atomTree->nodes[atomTree->size].child = atomTree->nodes[b].f;
								for(z=atomTree->nodes[atomTree->nodes[b].f].rchild; z>=0;) {
									TypePosition y = atomTree->nodes[z].rnext;
									if(atomTree->nodes[z].parent<0 || (atomTree->nodes[atomTree->nodes[z].parent].state&bitNew)) {
										atomTree->nodes[z].f = atomTree->size;
										atomTree->nodes[z].rnext = atomTree->nodes[atomTree->size].rchild;
										if(atomTree->nodes[atomTree->size].rchild >= 0)
											atomTree->nodes[atomTree->nodes[atomTree->size].rchild].rprec = z;
										atomTree->nodes[atomTree->size].rchild = z;
									} else {
										if(atomTree->nodes[atomTree->nodes[z].parent].f<0) {
											atomTree->nodes[atomTree->nodes[z].parent].f = atomTree->size;
											atomTree->nodes[atomTree->nodes[z].parent].rnext = atomTree->nodes[atomTree->size].rchild;
											if(atomTree->nodes[atomTree->size].rchild >= 0)
												atomTree->nodes[atomTree->nodes[atomTree->size].rchild].rprec = atomTree->nodes[z].parent;
											atomTree->nodes[atomTree->size].rchild = atomTree->nodes[z].parent;
											if(atomTree->nodes[atomTree->nodes[z].parent].state&bitSpe)
												stack[sizeSta++] = atomTree->nodes[z].parent;
										}
									}
									z = y;
								}
							}
						}
					}
					atomTree->size++;
				}
			}
			for(x=0; x<sizeNot; x++) {
				if(atomTree->nodes[not[x]].f<0) {
					TypePosition c;
					atomTree->nodes[not[x]].f = atomTree->nodes[atomTree->nodes[not[x]].child].f;
					atomTree->nodes[not[x]].rnext = atomTree->nodes[atomTree->nodes[not[x]].f].rchild;
					if(atomTree->nodes[atomTree->nodes[not[x]].f].rchild >= 0)
						atomTree->nodes[atomTree->nodes[atomTree->nodes[not[x]].f].rchild].rprec = not[x];
					atomTree->nodes[atomTree->nodes[not[x]].f].rchild = not[x];
					for(c=atomTree->nodes[not[x]].child; c>=0; c=atomTree->nodes[c].sibling)
						atomTree->nodes[atomTree->nodes[c].rprec].rnext = atomTree->nodes[c].rnext;
				}
			}
			sizeNot = 0;	
			for(x=0; x<sizeNew; x++) {
				TypePosition ref, a;
				ref = atomTree->nodes[atomTree->nodes[new[x]].child].f;
				for(a=atomTree->nodes[atomTree->nodes[new[x]].child].sibling; a>=0 && atomTree->nodes[a].f==ref; a=atomTree->nodes[a].sibling);
				if(a>=0) {
					special[sizeSpe++] = new[x];
					atomTree->nodes[new[x]].state |= bitSpe;
				} else {
					not[sizeNot++] = new[x];
				}
			}
			for(x=0;x<sizeNew; x++)
				atomTree->nodes[new[x]].state &= ~bitNew;
			sizeNew = 0;
			k--;
		} while(k>stop && sizeSpe>0);
	}
flushProgress();
	monfree((void*)last);
	monfree((void*)prec);
	monfree((void*)stack);
	monfree((void*)special);
	monfree((void*)not);
	monfree((void*)new);
	res = (TypeSetOfSequences*) monmalloc(sizeof(TypeSetOfSequences));
	res->name = set->name;
	res->table = set->table;
	res->number = set->number;
	res->cardinal = 0;
	tmp = (TypePosition*) monmalloc(size*sizeof(TypePosition));
	for(i=0; i<size; i++)
		tmp[i] = -1;
	for(i=0; i<size; i++)
		if(tmp[i] < 0) {
			TypePosition a;
			for(a=i; atomTree->nodes[a].parent>=0; a=atomTree->nodes[a].parent);
			setCode(res->cardinal++, a, atomTree, tmp);
		}
	monfree((void*)atomTree->nodes);
	monfree((void*)atomTree);
	res->size = (TypePosition*) monmalloc(res->number*sizeof(TypePosition));
	res->sequence = (TypeSymbol**) monmalloc(res->number*sizeof(TypeSymbol*));
	i = 0;
	for(n=0; n<res->number; n++) {
		TypePosition p;
		res->size[n] = set->size[n];
		res->sequence[n] = (TypeSymbol*) monmalloc((res->size[n])*sizeof(TypeSymbol));
		for(p=0; p<set->size[n]; p++)
			res->sequence[n][p] = tmp[i++];
	}
	monfree((void*)tmp);
	return res;
}



TypePosition *getTableCode(TypePosition *length, TypePosition size) {
	TypePosition i, j, *tab, *tmp, *off, min, max;
	min = length[0];
	max = length[0];
	for(i=1; i<size; i++) {
		if(length[i] < min)
			min = length[i];
		if(length[i] > max)
			max = length[i];
	}
	tmp = (TypePosition*) monmalloc((max-min+1)*sizeof(TypePosition));
	off = (TypePosition*) monmalloc((max-min+1)*sizeof(TypePosition));
	for(j=0; j<=max-min; j++)
		tmp[j] = 0;
	for(i=0; i<size; i++)
		tmp[length[i]-min]++;
	off[0] = 0;
	for(j=1; j<=max-min; j++)
		off[j] = off[j-1]+tmp[j-1];
	for(j=0; j<=max-min; j++)
		tmp[j] = 0;
	tab = (TypePosition*) monmalloc(size*sizeof(TypePosition));
	for(i=0; i<size; i++) {
		if((length[i]-min) == 0) {
			tab[size-1-tmp[length[i]-min]++] = i;
		} else {
			tab[size-1-(off[length[i]-min]+tmp[length[i]-min]++)] = i;
		}
	}
	monfree((void*)tmp);
	monfree((void*)off);
	return tab;
}

void setCode(TypePosition code, TypePosition n, TypeAtomTree *atomTree, TypePosition *tmp) {
	if(atomTree->nodes[n].child < 0)
		tmp[n] = code;
	else {
		TypePosition c;
		for(c=atomTree->nodes[n].child; c>=0; c=atomTree->nodes[c].sibling)
			setCode(code, c, atomTree, tmp);
	}
}

void initAtomNode(TypeAtomNode *n) {
	n->parent = -1;
	n->child = -1;
	n->sibling = -1;
	n->f = -1;
	n->rchild = -1;
	n->rnext = -1;
	n->rprec = -1;
	n->state = 0;
}





