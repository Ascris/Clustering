#include <math.h>
#include <gsl/gsl_cdf.h>
#include "MarkovEmbed.h"
#include "Utils.h"

static void productMat(double **a, double **b, double **c, int size);
static void setMat(double **a, double **b, int size);
static void product(double **a, double *x, double *y, int size);
static void printMat(double **a,int size);
static void printVect(double *a,int size);

TypeMarkovModel *estimateMarkovModel(TypeSetOfSequences *set) {
	TypePosition pos, ntot = 0;
	double *occ;
	TypeSymbol symb;
	TypeMarkovModel *m;
	TypeNumber n;

	m = (TypeMarkovModel*) monmalloc(sizeof(TypeMarkovModel));
	m->cardinal = set->cardinal;
	m->init = (double*) monmalloc(m->cardinal*sizeof(double));
	m->trans = (double**) monmalloc(m->cardinal*sizeof(double*));
	occ = monmalloc(m->cardinal*sizeof(double));
	for(symb=0; symb<m->cardinal; symb++) {
		TypeSymbol dest;
		occ[symb] = 0.;
		m->init[symb] = 0.;
		m->trans[symb] = (double*) monmalloc(m->cardinal*sizeof(double));
		for(dest=0; dest<m->cardinal; dest++)
			m->trans[symb][dest] = 0.;
	}
	for(n=0; n<set->number; n++) {
		TypeSymbol last;
		TypePosition start;
		for(start=0; start<set->size[n] && set->sequence[n][start]>=m->cardinal; start++);
		if(start<set->size[n]) {
			last = set->sequence[n][start];
			start++;
		}
		for(pos=start; pos<set->size[n]; pos++) {
			if(set->sequence[n][pos]<m->cardinal) {
				occ[last]++;
				ntot++;
				(m->trans[last][set->sequence[n][pos]])++;
				last = set->sequence[n][pos];
			} else {
				for(; pos<set->size[n] && set->sequence[n][pos]>=m->cardinal; pos++);
				if(pos<set->size[n])
					last = set->sequence[n][pos];
			}
		}
	}
	for(symb=0; symb<m->cardinal; symb++) {
		TypeSymbol dest;
		m->init[symb] = ((double)occ[symb])/((double)ntot);
		for(dest=0; dest<m->cardinal; dest++)
			m->trans[symb][dest] /= occ[symb];
	}
	monfree((void*)occ);
	return m;
}

void freeModel(TypeMarkovModel *model) {
	TypeSymbol symb;
	for(symb=0; symb<model->cardinal; symb++)
		monfree((void*) model->trans[symb]);
	monfree((void*)model->trans);
	monfree((void*)model->init);
	monfree((void*)model);
}


TypePosition **getNext(TypeSymbol *w, TypePosition size, TypeSymbol cardinal) {
	TypePosition **next, *bord, p, q;
	TypeSymbol a;

	next = (TypePosition**) monmalloc(cardinal*sizeof(TypePosition*));
	bord = (TypePosition*) monmalloc((size+1)*sizeof(TypePosition));
	for(a=0; a<cardinal; a++)
		next[a] = (TypePosition*) monmalloc((size+1)*sizeof(TypePosition));
	for(a=0; a<cardinal; a++)
		next[a][0] = -1;
	bord[0] = -1;
	for(p=0; p<size; p++) {
		for(q=bord[p]; q >= 0 && w[q] != w[p]; q=bord[q]);
		bord[p+1] = q+1;
	}
	for(p=1; p<size; p++) {
		if(w[p]<cardinal) {
			for(a=0; a<w[p]; a++)
				if(bord[p]>0)
					next[a][p-1] = next[a][bord[p]-1];
				else
					next[a][p-1] = (a == w[0])?0:-1;
			next[w[p]][p-1] = p;
			for(a=w[p]+1; a<cardinal; a++)
				if(bord[p]>0)
					next[a][p-1] = next[a][bord[p]-1];
				else
					next[a][p-1] = (a == w[0])?0:-1;
		} else
			for(a=0; a<cardinal; a++)
				next[a][p-1] = -1;
	}
	for(a=0; a<cardinal; a++)
		if(bord[size]>0)
			next[a][size-1] = next[a][bord[p]-1];
		else
			next[a][size-1] = (a == w[0])?0:-1;
	monfree((void*) bord);
	return next;
}	

double getProb(TypeSymbol *w, TypePosition size, TypePosition nocc, TypeSetOfSequences *set, TypeMarkovModel *m) {
	TypePosition **next, o, p;
	TypeNumber n;
	TypeSymbol a,b;
	double **state0, **state1, **tmp, **zero0, **zero1, *start, sink, sum;
	next = getNext(w, size, set->cardinal);
	start = (double*) monmalloc(nocc*sizeof(double));
	state0 = (double**) monmalloc(nocc*sizeof(double*));
	state1 = (double**) monmalloc(nocc*sizeof(double*));
	zero0 = (double**) monmalloc(nocc*sizeof(double*));
	zero1 = (double**) monmalloc(nocc*sizeof(double*));
	for(o=0; o<nocc; o++) {
		state0[o] = (double*) monmalloc(size*sizeof(double));
		state1[o] = (double*) monmalloc(size*sizeof(double));
		zero0[o] = (double*) monmalloc(set->cardinal*sizeof(double));
		zero1[o] = (double*) monmalloc(set->cardinal*sizeof(double));
	}
	start[0] = 1.;
	for(o=1; o<nocc; o++)
		start[o] = 0.;
	sink = 0.;
	for(n=0; n<set->number; n++) {
		TypePosition o, i;
		for(o=0; o<nocc; o++) {
			for(p=1; p<size; p++)
				state0[o][p] = 0.;
			for(a=0; a<w[0]; a++)
				zero0[o][a] = m->init[a]*start[o];
			state0[o][0] = m->init[w[0]]*start[o];
			zero0[o][w[0]] = 0.;
			for(a=w[0]+1; a<set->cardinal; a++)
				zero0[o][a] = m->init[a]*start[o];
		}
		for(i=1; i<set->size[n]; i++) {
			for(o=0; o<nocc; o++) {
				for(p=0; p<size; p++)
					state1[o][p] = 0.;
				for(a=0; a<set->cardinal; a++)
					zero1[o][a] = 0.;
			}
			for(o=0; o<nocc; o++) {
				for(a=0; a<set->cardinal; a++) {
					for(b=0; b<w[0]; b++)
						zero1[o][b] += zero0[o][a]*m->trans[a][b];
					state1[o][0] += zero0[o][a]*m->trans[a][w[0]];
					for(b=w[0]+1; b<set->cardinal; b++)
						zero1[o][b] += zero0[o][a]*m->trans[a][b];
				}
				for(p=0; p<(size-1); p++)
					for(a=0; a<set->cardinal; a++)
						if(next[a][p] >= 0)
							state1[o][next[a][p]] += state0[o][p]*m->trans[w[p]][a];
						else
							zero1[o][a] += state0[o][p]*m->trans[w[p]][a];
			}
			for(o=0; o<nocc-1; o++) {
				for(a=0; a<set->cardinal; a++)
					if(next[a][size-1] >= 0)
						state1[o+1][next[a][size-1]] += state0[o][size-1]*m->trans[w[size-1]][a];
					else
						zero1[o+1][a] += state0[o][size-1]*m->trans[w[size-1]][a];
			}
			sink += state0[nocc-1][size-1];
			tmp = state0; state0 = state1; state1 = tmp;
			tmp = zero0; zero0 = zero1; zero1 = tmp;
		}
		for(o=0; o<nocc; o++) {
			start[o] = 0;
			for(a=0; a<set->cardinal; a++)
				start[o] += zero0[o][a];
			for(p=0; p<(size-1); p++)
				start[o] += state0[o][p];
		}
		for(o=0; o<nocc-1; o++)
			start[o+1] += state0[o][size-1];
		sink += state0[nocc-1][size-1];
	}
	for(o=0; o<nocc; o++) {
		monfree((void*)state0[o]);
		monfree((void*)state1[o]);
		monfree((void*)zero0[o]);
		monfree((void*)zero1[o]);
	}
	monfree((void*)state0);
	monfree((void*)state1);
	monfree((void*)zero0);
	monfree((void*)zero1);
	monfree((void*)start);
	for(a=0; a<set->cardinal; a++)
		monfree((void*)next[a]);
	monfree((void*)next);
	return sink;
}

void productMat(double **a, double **b, double **c, int size) {
	int i, j, k;
	for(i=0; i<size; i++)
		for(j=0; j<size; j++) {
			c[i][j] = 0.;
			for(k=0; k<size; k++)
				c[i][j] += a[i][k]*b[k][j];
		}
}
void setMat(double **a, double **b, int size) {
	int i, j;
	for(i=0; i<size; i++)
		for(j=0; j<size; j++)
			b[i][j] = a[i][j];
}

void product(double **a, double *x, double *y, int size) {
	int i, j;
	for(i=0; i<size; i++) {
		y[i] = 0.;
		for(j=0; j<size; j++)
			y[i] += a[i][j]*x[j];
	}
}


void printMat(double **a,int size) {
	int i, j;
	for(i=0; i<size; i++) {
		for(j=0; j<size; j++)
			printf("%.2lE\t", a[i][j]);
		printf("\n");
	}
}

void printVect(double *a,int size) {
	int i, j;
	for(i=0; i<size; i++)
			printf("%.2lE\n", a[i]);
}

double getProbMat(TypeSymbol *w, TypePosition size, TypePosition nocc, TypeSetOfSequences *set, TypeMarkovModel *m) {
	TypePosition **next, o, p, max, eMax, e;
	TypeNumber n;
	TypeSymbol a,b;
	double *start, sink, sum, ***mat, **ms0, **ms1, **mstmp, *x, *y;
	int **stateInd, **zeroInd, maxInd, i, j;
	next = getNext(w, size, set->cardinal);
	start = (double*) monmalloc(nocc*sizeof(double));
	stateInd = (int**) monmalloc(nocc*sizeof(int*));
	zeroInd = (int**) monmalloc(nocc*sizeof(int*));
	maxInd = 0;
	for(o=0; o<nocc; o++) {
		zeroInd[o] = (int*) monmalloc(set->cardinal*sizeof(int*));
		stateInd[o] = (int*) monmalloc(size*sizeof(int));
		for(a=0; a<set->cardinal; a++)
			zeroInd[o][a] = maxInd++;
		for(p=0; p<size; p++)
			stateInd[o][p] = maxInd++;
	}
	max = set->size[0];
	for(n=1; n<set->number; n++)
		if(set->size[n]>max)
			max = set->size[n];
	eMax = ((TypePosition) floor(log((double)max)/log(2.)))+1;
	mat = (double***) monmalloc(eMax*sizeof(double**));
	for(e=0; e<eMax; e++) {
		mat[e] = (double**) monmalloc((maxInd+1)*sizeof(double*));
		for(i=0; i<=maxInd; i++)
			mat[e][i] = (double*) monmalloc((maxInd+1)*sizeof(double));
	}
	ms0 = (double**) monmalloc((maxInd+1)*sizeof(double*));
	ms1 = (double**) monmalloc((maxInd+1)*sizeof(double*));
	x = (double*) monmalloc((maxInd+1)*sizeof(double));
	y = (double*) monmalloc((maxInd+1)*sizeof(double));
	for(i=0; i<=maxInd; i++) {
		ms0[i] = (double*) monmalloc((maxInd+1)*sizeof(double));
		ms1[i] = (double*) monmalloc((maxInd+1)*sizeof(double));
	}
	for(i=0; i<=maxInd; i++)
		for(j=0; j<=maxInd; j++)
			mat[0][i][j] = 0.;
	for(o=0; o<nocc; o++) {
		for(a=0; a<set->cardinal; a++) {
			if(w[0]<set->cardinal) {
				for(b=0; b<w[0]; b++) {
					mat[0][zeroInd[o][b]][zeroInd[o][a]] = m->trans[a][b];
//				printf("o %d set zero %d %d [%d,%d] = %.2lE\n", o, b, a, zeroInd[o][b], zeroInd[o][a], m->trans[a][b]);
				}
				mat[0][stateInd[o][0]][zeroInd[o][a]] = m->trans[a][w[0]];
				for(b=w[0]+1; b<set->cardinal; b++) {
					mat[0][zeroInd[o][b]][zeroInd[o][a]] = m->trans[a][b];
//					printf("o %d set zero %d %d [%d,%d] = %.2lE\n", o, b, a, zeroInd[o][b], zeroInd[o][a], m->trans[a][b]);
				}
			} else {
				for(b=0; b<set->cardinal; b++)
					mat[0][zeroInd[o][b]][zeroInd[o][a]] = m->trans[a][b];
			}
		}
		for(p=0; p<(size-1); p++)
			for(a=0; a<set->cardinal; a++)
				if(w[p]<set->cardinal) {
					if(next[a][p] >= 0)
						mat[0][stateInd[o][next[a][p]]][stateInd[o][p]] += m->trans[w[p]][a];
					else
						mat[0][zeroInd[o][a]][stateInd[o][p]] += m->trans[w[p]][a];
				} else
					mat[0][zeroInd[o][a]][stateInd[o][p]] += m->init[a];
	}
	for(o=0; o<nocc-1; o++) {
		for(a=0; a<set->cardinal; a++)
			if(w[size-1]<set->cardinal) {
				if(next[a][size-1] >= 0)
					mat[0][stateInd[o+1][next[a][size-1]]][stateInd[o][size-1]] = m->trans[w[size-1]][a];
				else
					mat[0][zeroInd[o+1][a]][stateInd[o][size-1]] = m->trans[w[size-1]][a];
			} else
				mat[0][zeroInd[o][a]][stateInd[o][size-1]] = m->init[a];
	}
	mat[0][maxInd][stateInd[nocc-1][size-1]] = 1.;
	mat[0][maxInd][maxInd] = 1.;
	for(e=1; e<eMax; e++)
		productMat(mat[e-1], mat[e-1], mat[e], maxInd+1);
/*	for(e=0; e<eMax; e++) {
		printf("\n\n %d\n", e+1);
		printMat(mat[e], maxInd+1);
	}
*/	start[0] = 1.;
	for(o=1; o<nocc; o++)
		start[o] = 0.;
	sink = 0.;
	for(n=0; n<set->number; n++) {
		TypePosition o, length;
		int k;
		for(o=0; o<nocc; o++) {
			for(p=1; p<size; p++)
				x[stateInd[o][p]] = 0.;
			if(w[0]<set->cardinal) {
				for(a=0; a<w[0]; a++)
					x[zeroInd[o][a]] = m->init[a]*start[o];
				x[stateInd[o][0]] = m->init[w[0]]*start[o];
				x[zeroInd[o][w[0]]] = 0.;
				for(a=w[0]+1; a<set->cardinal; a++)
					x[zeroInd[o][a]] = m->init[a]*start[o];
			} else
				for(a=0; a<set->cardinal; a++)
					x[zeroInd[o][a]] = m->init[a]*start[o];
		}
		x[maxInd] = sink;
		length = set->size[n];
		for(e=0; e<eMax && length%2 == 0; e++)
			length /= 2;
		setMat(mat[e], ms0, maxInd+1);
		length /= 2; e++;
		for(; e<eMax && length > 0; e++) {
			if(length%2 != 0) {
				productMat(ms0, mat[e], ms1, maxInd+1);
				mstmp = ms1; ms1 = ms0; ms0 = mstmp;
			}
			length /= 2;
		}
/*		printf("length %d\n", length);
		printf("\n\n ms0 %d\n", n);
		printMat(ms0, maxInd+1);
*/		product(ms0, x, y, maxInd+1);
/*		printf("\n\n x %d\n", n);
		printVect(x, maxInd+1);
		printf("\n\n y %d\n", n);
		printVect(y, maxInd+1);
		printf("\n\n");
*/		for(o=0; o<nocc; o++) {
			start[o] = 0;
			for(a=0; a<set->cardinal; a++)
				start[o] += y[zeroInd[o][a]];
			for(p=0; p<(size-1); p++)
				start[o] += y[stateInd[o][p]];
		}
		for(o=0; o<nocc-1; o++)
			start[o+1] += y[stateInd[o][size-1]];
		sink = y[maxInd] + y[stateInd[nocc-1][size-1]];
	}
	for(o=0; o<nocc; o++) {
		monfree((void*)stateInd[o]);
		monfree((void*)zeroInd[o]);
	}
	monfree((void*)stateInd);
	monfree((void*)zeroInd);
	for(i=0; i<=maxInd; i++) {
		monfree((void*)ms0[i]);
		monfree((void*)ms1[i]);
	}
	monfree((void*)ms0);
	monfree((void*)ms1);
	monfree((void*)x);
	monfree((void*)y);
	for(e=0; e<eMax; e++) {
		for(i=0; i<=maxInd; i++)
			monfree((void*)mat[e][i]);
		monfree((void*)mat[e]);
	}
	monfree((void*)mat);
	monfree((void*)start);
	for(a=0; a<set->cardinal; a++)
		monfree((void*)next[a]);
	monfree((void*)next);
	return sink;
}

void getBounds(TypePosition *start, TypePosition *end, TypePosition nocc, TypeSetOfSequences *set, TypePosition *sorted, double threshold, TypeMarkovModel *m) {
	TypePosition **next, o, p, eMax, e, prec, u, v, cardinal = 2;
	TypeNumber n;
	TypeSymbol a,b;
	double *cumul, sink, sum, ***mat, **ms0, **ms1, **mstmp, *x, *y, min, max;
	int **stateInd, **zeroInd, maxInd, i, j;
	cumul = (double*) monmalloc(nocc*sizeof(double));
	stateInd = (int**) monmalloc(nocc*sizeof(int*));
	zeroInd = (int**) monmalloc(nocc*sizeof(int*));
	maxInd = 0;
	max = 0.; min = 1.;
	u = *start; v = *end;
	for(a=0; a<set->cardinal; a++)
		for(b=0; b<set->cardinal; b++) {
			if(m->trans[a][b]>max)
				max = m->trans[a][b];
			if(m->trans[a][b]<min)
				min = m->trans[a][b];
		}
	next = (TypePosition**) monmalloc(cardinal*sizeof(TypePosition*));
	for(v=*start; v<*end; v++) {
	for(a=0; a<cardinal; a++)
		next[a] = (TypePosition*) monmalloc((v+1)*sizeof(TypePosition));
	for(a=0; a<cardinal; a++)
		next[a][0] = -1;
	for(p=1; p<v; p++) {
		next[0][p-1] = p;
		next[1][p-1] = -1;
	}
	next[0][v-1] = v-2;
	next[1][v-1] = -1;

	for(o=0; o<nocc; o++) {
		zeroInd[o] = (int*) monmalloc(set->cardinal*sizeof(int*));
		stateInd[o] = (int*) monmalloc(v*sizeof(int));
		for(a=0; a<set->cardinal; a++)
			zeroInd[o][a] = maxInd++;
		for(p=0; p<v; p++)
			stateInd[o][p] = maxInd++;
	}
	eMax = ((TypePosition) floor(log((double)sorted[set->number-1])/log(2.)))+1;
	mat = (double***) monmalloc(eMax*sizeof(double**));
	for(e=0; e<eMax; e++) {
		mat[e] = (double**) monmalloc((maxInd+1)*sizeof(double*));
		for(i=0; i<=maxInd; i++)
			mat[e][i] = (double*) monmalloc((maxInd+1)*sizeof(double));
	}
	ms0 = (double**) monmalloc((maxInd+1)*sizeof(double*));
	ms1 = (double**) monmalloc((maxInd+1)*sizeof(double*));
	x = (double*) monmalloc((maxInd+1)*sizeof(double));
	y = (double*) monmalloc((maxInd+1)*sizeof(double));
	for(i=0; i<=maxInd; i++) {
		ms0[i] = (double*) monmalloc((maxInd+1)*sizeof(double));
		ms1[i] = (double*) monmalloc((maxInd+1)*sizeof(double));
	}
	for(i=0; i<=maxInd; i++)
		for(j=0; j<=maxInd; j++)
			mat[0][i][j] = 0.;
	for(o=0; o<nocc; o++) {
		mat[0][zeroInd[o][0]][zeroInd[o][0]] = max;
		mat[0][zeroInd[o][0]][zeroInd[o][1]] = max;
		mat[0][zeroInd[o][1]][zeroInd[o][0]] = 1.-max;
		mat[0][zeroInd[o][1]][zeroInd[o][1]] = 1.-max;
		for(p=0; p<(v-1); p++) {
			mat[0][stateInd[o][p+1]][stateInd[o][p]] = max;
			mat[0][zeroInd[o][1]][stateInd[o][p]] = 1.-max;
		}
	}
	for(o=0; o<nocc-1; o++) {
		mat[0][stateInd[o+1][v-2]][stateInd[o][v-1]] = max;
		mat[0][zeroInd[o+1][1]][stateInd[o][v-1]] = 1.-max;
	}
	mat[0][maxInd][stateInd[nocc-1][v-1]] = 1.;
	mat[0][maxInd][maxInd] = 1.;
	for(e=1; e<eMax; e++)
		productMat(mat[e-1], mat[e-1], mat[e], maxInd+1);
	cumul[0] = 1.;
	for(o=1; o<nocc; o++)
		cumul[o] = 0.;
	sink = 0.;
	prec = 0;
	for(n=0; n<set->number; n++) {
		TypePosition o, length;
		int k;
		for(o=0; o<nocc; o++) {
			for(p=1; p<v; p++)
				x[stateInd[o][p]] = 0.;
			x[zeroInd[o][0]] = 0.;
			x[zeroInd[o][1]] = (1.-max)*cumul[o];
			x[stateInd[o][0]] = max*cumul[o];
		}
		x[maxInd] = sink;
		if((length = sorted[n]-prec)>0) {
			if(prec == 0) {
				for(e=0; e<eMax && length%2 == 0; e++)
					length /= 2;
				setMat(mat[e], ms0, maxInd+1);
				length /= 2; e++;
				for(; e<eMax && length > 0; e++) {
					if(length%2 != 0) {
						productMat(ms0, mat[e], ms1, maxInd+1);
						mstmp = ms1; ms1 = ms0; ms0 = mstmp;
					}
					length /= 2;
				}
			} else {
				for(e=0; e<eMax && length > 0; e++) {
					if(length%2 != 0) {
						productMat(ms0, mat[e], ms1, maxInd+1);
						mstmp = ms1; ms1 = ms0; ms0 = mstmp;
					}
					length /= 2;
				}
			}
		}
		prec = sorted[n];
		product(ms0, x, y, maxInd+1);
		for(o=0; o<nocc; o++) {
			cumul[o] = 0;
			for(a=0; a<set->cardinal; a++)
				cumul[o] += y[zeroInd[o][a]];
			for(p=0; p<(v-1); p++)
				cumul[o] += y[stateInd[o][p]];
		}
		for(o=0; o<nocc-1; o++)
			cumul[o+1] += y[stateInd[o][v-1]];
		sink = y[maxInd] + y[stateInd[nocc-1][v-1]];
	}
	printf("Max sink[%d] %.2lE\n", v, sink);
	for(o=0; o<nocc; o++) {
		monfree((void*)stateInd[o]);
		monfree((void*)zeroInd[o]);
	}
	for(i=0; i<=maxInd; i++) {
		monfree((void*)ms0[i]);
		monfree((void*)ms1[i]);
	}
	monfree((void*)ms0);
	monfree((void*)ms1);
	monfree((void*)x);
	monfree((void*)y);
	for(e=0; e<eMax; e++) {
		for(i=0; i<=maxInd; i++)
			monfree((void*)mat[e][i]);
		monfree((void*)mat[e]);
	}
	monfree((void*)mat);
	}
	for(v=*start; v<*end; v++) {
	for(o=0; o<nocc; o++) {
		zeroInd[o] = (int*) monmalloc(cardinal*sizeof(int*));
		stateInd[o] = (int*) monmalloc(v*sizeof(int));
		for(a=0; a<set->cardinal; a++)
			zeroInd[o][a] = maxInd++;
		for(p=0; p<v; p++)
			stateInd[o][p] = maxInd++;
	}
	eMax = ((TypePosition) floor(log((double)sorted[set->number-1])/log(2.)))+1;
	mat = (double***) monmalloc(eMax*sizeof(double**));
	for(e=0; e<eMax; e++) {
		mat[e] = (double**) monmalloc((maxInd+1)*sizeof(double*));
		for(i=0; i<=maxInd; i++)
			mat[e][i] = (double*) monmalloc((maxInd+1)*sizeof(double));
	}
	ms0 = (double**) monmalloc((maxInd+1)*sizeof(double*));
	ms1 = (double**) monmalloc((maxInd+1)*sizeof(double*));
	x = (double*) monmalloc((maxInd+1)*sizeof(double));
	y = (double*) monmalloc((maxInd+1)*sizeof(double));
	for(i=0; i<=maxInd; i++) {
		ms0[i] = (double*) monmalloc((maxInd+1)*sizeof(double));
		ms1[i] = (double*) monmalloc((maxInd+1)*sizeof(double));
	}
	for(i=0; i<=maxInd; i++)
		for(j=0; j<=maxInd; j++)
			mat[0][i][j] = 0.;
	for(o=0; o<nocc; o++) {
		mat[0][zeroInd[o][0]][zeroInd[o][0]] = min;
		mat[0][zeroInd[o][0]][zeroInd[o][1]] = min;
		mat[0][zeroInd[o][1]][zeroInd[o][0]] = 1.-min;
		mat[0][zeroInd[o][1]][zeroInd[o][1]] = 1.-min;
		for(p=0; p<(v-1); p++) {
			mat[0][stateInd[o][p+1]][stateInd[o][p]] = min;
			mat[0][zeroInd[o][1]][stateInd[o][p]] = 1.-min;
		}
	}
	for(o=0; o<nocc-1; o++) {
		mat[0][stateInd[o+1][v-2]][stateInd[o][0]] = min;
		mat[0][zeroInd[o+1][1]][stateInd[o][v-1]] = 1.-min;
	}
	mat[0][maxInd][stateInd[nocc-1][v-1]] = 1.;
	mat[0][maxInd][maxInd] = 1.;
	for(e=1; e<eMax; e++)
		productMat(mat[e-1], mat[e-1], mat[e], maxInd+1);
	cumul[0] = 1.;
	for(o=1; o<nocc; o++)
		cumul[o] = 0.;
	sink = 0.;
	prec = 0;
	for(n=0; n<set->number; n++) {
		TypePosition o, length;
		int k;
		for(o=0; o<nocc; o++) {
			for(p=1; p<v; p++)
				x[stateInd[o][p]] = 0.;
			x[zeroInd[o][0]] = 0.;
			x[zeroInd[o][1]] = (1.-min)*cumul[o];
			x[stateInd[o][0]] = min*cumul[o];
		}
		x[maxInd] = sink;
		if((length = sorted[n]-prec)>0) {
			if(prec == 0) {
				for(e=0; e<eMax && length%2 == 0; e++)
					length /= 2;
				setMat(mat[e], ms0, maxInd+1);
				length /= 2; e++;
				for(; e<eMax && length > 0; e++) {
					if(length%2 != 0) {
						productMat(ms0, mat[e], ms1, maxInd+1);
						mstmp = ms1; ms1 = ms0; ms0 = mstmp;
					}
					length /= 2;
				}
			} else {
				for(e=0; e<eMax && length > 0; e++) {
					if(length%2 != 0) {
						productMat(ms0, mat[e], ms1, maxInd+1);
						mstmp = ms1; ms1 = ms0; ms0 = mstmp;
					}
					length /= 2;
				}
			}
		}
		prec = sorted[n];
/*		printf("length %d\n", length);
		printf("\n\n ms0 %d\n", n);
		printMat(ms0, maxInd+1);
*/		product(ms0, x, y, maxInd+1);
/*		printf("\n\n x %d\n", n);
		printVect(x, maxInd+1);
		printf("\n\n y %d\n", n);
		printVect(y, maxInd+1);
		printf("\n\n");
*/		for(o=0; o<nocc; o++) {
			cumul[o] = 0;
			for(a=0; a<set->cardinal; a++)
				cumul[o] += y[zeroInd[o][a]];
			for(p=0; p<(v-1); p++)
				cumul[o] += y[stateInd[o][p]];
		}
		for(o=0; o<nocc-1; o++)
			cumul[o+1] += y[stateInd[o][v-1]];
		sink = y[maxInd] + y[stateInd[nocc-1][v-1]];
	}
	printf("Min sink[%d] %.2lE\n", v, sink);
	for(o=0; o<nocc; o++) {
		monfree((void*)stateInd[o]);
		monfree((void*)zeroInd[o]);
	}
	for(i=0; i<=maxInd; i++) {
		monfree((void*)ms0[i]);
		monfree((void*)ms1[i]);
	}
	monfree((void*)ms0);
	monfree((void*)ms1);
	monfree((void*)x);
	monfree((void*)y);
	for(e=0; e<eMax; e++) {
		for(i=0; i<=maxInd; i++)
			monfree((void*)mat[e][i]);
		monfree((void*)mat[e]);
	}
	monfree((void*)mat);
	}
	monfree((void*)stateInd);
	monfree((void*)zeroInd);
	monfree((void*)cumul);
}

double getProbApp(TypeSymbol *w, TypePosition size, TypePosition nocc, TypeSetOfSequences *set, TypeMarkovModel *m) {
	double p, *st0, *st1, *stmp;
	TypePosition i;
	unsigned int tot;
	TypeSymbol a, b;
	TypeNumber n;

	if(size == 0)
		return 1.;
	tot = 0;
	st0 = (double*) monmalloc(set->cardinal*sizeof(double));
	st1 = (double*) monmalloc(set->cardinal*sizeof(double));
	for(n=0; n<set->number; n++)
		tot += set->size[n];
	tot -= set->number*(size-1);
	for(a=0; a<set->cardinal; a++)
		st0[a] = 0.;
	if(w[0]<set->cardinal)
		st0[w[0]] = m->init[w[0]];
	else {
		TypeSymbol a;
		p = 0.;
		for(a=0; a<set->ambiguity.size[w[0]-set->cardinal]; a++)
			st0[set->ambiguity.set[w[0]-set->cardinal][a]] += m->init[set->ambiguity.set[w[0]-set->cardinal][a]];
	}
	for(i=1; i<size; i++) {
		for(a=0; a<set->cardinal; a++)
			st1[a] = 0.;
		if(w[i]<set->cardinal) {
			for(a=0; a<set->cardinal; a++)
				st1[w[i]] += st0[a]*m->trans[a][w[i]];
		} else {
			for(a=0; a<set->cardinal; a++)
				for(b=0; b<set->ambiguity.size[w[i]-set->cardinal]; b++)
					st1[b] += st0[a]*m->trans[a][set->ambiguity.set[w[i]-set->cardinal][b]];
		}
		stmp = st0; st0 = st1; st1 = stmp;
	}
	p = 0.;
	for(a=0; a<set->cardinal; a++)
		p += st0[a];
	monfree((void*)st0);
	monfree((void*)st1);
//printf("papp\t%lE\t%ld\n", p, tot);
	if(p >= 1.)
		p = 1.-0.00000000000000000000000001;
	if(p <= 0.)
		p = 0.00000000000000000000000001;
	return gsl_cdf_binomial_Q(1, p, tot);
}