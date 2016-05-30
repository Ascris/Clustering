#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "Alignment.h"
#include "Sequences.h"
#include "Distance.h"
#include "SuffixTree.h"
#include "PrefixDecode.h"
#include "PrefixCode.h"
#include "Utils.h"
#include "DistSeq.h"
#include "DistAln.h"


#define SIZE_BUFFER_CHAR 300
#define MAX_FUNC 10

#define HELPMESSAGE "\nusage: decode [options] <order> <input file> <output file>\n\nThe input file has to be in Fasta format.\n\nOptions:\n\n-h                  output help\n\n-f <f, s, m, x>     select the format of output for distances matrix \n                    (r: raw, p: phylip, n: nexus)\n\n-s <d, r, p>        select the type of alphabet (d: DNA, r: RNA, p: Protein)\n                    maybe useful to handle ambiguity characters\n\n-o <d, m>           select the type of output:\n                     * d: return the decoding of order <order> in Fasta format\n                     * m: return corresponding matrix of distances,\n			  format depending on option 'o' (raw by default)\n			      \n-d <l, s>           select the type of decoding:\n                     * l: use local decoding of order <order>\n                     * s: use sequence of sliding blocks  of length <order>\n		     \n-c <p, w>           select the mode of computing:\n                     * p: compute local decoding by pair of sequences\n                     * w: compute local decoding of the whole set of sequences\n\n-m <number>           select the type of dissimilarity:\n                     * 0: Euclidian distance\n                     * 1: Kullback-leiber discrepancy\n                     * 2: Pham\n                     * 3: Gilles 1\n                     * 4: Gilles 2\n                     * 5: Gilles 3\n                     * 6: Gilles 4\n"
#define AMBIMESSAGE "The ambiguity codes are for DNA and RNA:\n                                         * Y: T/U or C\n                                         * R: A or G\n                                         * M: A or C\n                                         * K: T/U or G\n                                         * W: T/U or A\n                                         * S: C or G\n                                         * B: T/U, C, or G\n                                         * D: T/U, A, or G\n                                         * H: T/U, C, or A\n                                         * V: C, A, or G\n                                         * N: T/U, C, A, or G\n\nand for protein:\n                                         * B: D or N\n                                         * Q or E\n                                         * X: Unidentified\n"
static double score(TypeSetOfSequences *dec);
static double scorebis(TypeSetOfSequences *dec);
static double findMode(TypeSetOfSequences *set, double tmin, double tmax, double tprec, TypeCodeScheme *scheme, TypeMarkovModel *model);

int main(int argc, char **argv) {
	TypePosition orderstart=1, orderend=10;
	char option[256], inputFileName[SIZE_BUFFER_CHAR], outputFileName[SIZE_BUFFER_CHAR], bufferOutput[SIZE_BUFFER_CHAR], *table, 
	outputFormat = 'r', typeDec = 'l', typeAlphabet = '?', typeCalc = 'g', type = 't';
	TypeSetOfSequences *set, seq;
	TypeAlignment aln, atmp;
	int fixed = 0;
	double threshold = 0.001, tmin = 1E-20, tmax=0.1, tstep = 0.00001, qmin = -25, qmax = -3, qprec = 0.5;
	double thre;
	TypeNumber n;
	TypeDistance distA, distB;
	TypePosition l, tot, lmax = 50;
	TypeSuffixTree *suffixTree;
	TypeMarkovModel *model;
	TypeCodeScheme *scheme;
/*	TypeDistFunction *distfunc[MAX_FUNC]=
	{computeProba, computeKullbackLeiber1, computePham, computeCommon, computeCommonBis, computeGillesPham, computeMatchesBis, computeMatches, computeAlex, computeAlexBis};
*/		
	FILE *fi, *fo;
	int i = 1, typeDist = 0;

	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &outputFormat) == 1)
				i++;
		}
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &typeAlphabet) == 1)
				i++;
		}
		if(option['c']) {
			option['c'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &typeCalc) == 1)
				i++;
		}
		if(option['m']) {
			option['m'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &tmin) == 1)
				i++;
			if(typeDist >= MAX_FUNC)
				typeDist = 0;
		}
		if(option['t']) {
			option['t'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &threshold) == 1)
				i++;
		}
		if(option['y']) {
			option['y'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &type) == 1)
				i++;
		}
	if(option['h']) {
			printf("%s\n", HELPMESSAGE);
			exitProg(ExitOk, NULL);
		}
	}
	if (i>=argc || sscanf(argv[i++], "%s", inputFileName) != 1) exitProg(ErrorArgument, HELPMESSAGE);
	if (i>=argc || sscanf(argv[i++], "%s", outputFileName) != 1) exitProg(ErrorArgument, HELPMESSAGE);

	switch(typeAlphabet) {
		case 'd':
			table = (char*) monmalloc((strlen(DNA)+1)*sizeof(char));
			strcpy(table, DNA);
			break;
		case 'r':
			table = (char*) monmalloc((strlen(RNA)+1)*sizeof(char));
			strcpy(table, RNA);
			break;
		case 'p':
			table = (char*) monmalloc((strlen(PRO)+1)*sizeof(char));
			strcpy(table, PRO);
			break;
		case '?':
		default:
			table = (char*) monmalloc(sizeof(char));
			table[0] = '\0';
	}
	if(fi = fopen(inputFileName, "r")) {
		aln = readAlignement(fi, table, typeAlphabet == '?');
		switch(typeAlphabet) {
			case 'd':
			case 'r':
				aln.ambiguity = getXNAAmbiguity();
				break;
			case 'p':
				aln.ambiguity = getProteinAmbiguity();
				break;
			case '?':
			default:
				aln.ambiguity.number = 0;
		}
		aln.cardinal -= aln.ambiguity.number;
		fclose(fi);
	} else {
		exitProg(ErrorReading, inputFileName);
	}
	fixAlignmentAmbiguity(&aln);
	set=toSequences(&aln);

	if(!(fo = fopen(outputFileName, "w")))
		exitProg(ErrorWriting, outputFileName);
	distA = computeWholeDistancePairAln(aln, computeNorm1Aln);
	scheme = (TypeCodeScheme*) monmalloc(sizeof(TypeCodeScheme));
	scheme->suffixTree = getSuffixTree(set);
	scheme->code = (TypePosition*) monmalloc(scheme->suffixTree->size*sizeof(TypePosition));
	scheme->buffSize = INC_SIZE_CODE;
	scheme->lengthCode = (TypePosition*) monmalloc(scheme->buffSize*sizeof(TypePosition));
	if(type == 't') {
		int l;
		model = estimateMarkovModel(set);
//		for(thre=tmin; thre<=tmax; thre *= 10.0) {
		for(l=tmin; l<=-1; l++) {
			double t;
			int k;
			thre = pow(10.0, (double) l);
			for(k=0; k<10; k++) {
//			for(t=thre; t<thre*10; t+=thre) {
				double corr, sc;
				TypeSetOfSequences *dec;
				t = ((double)k+1.)*thre;
				scheme->cardCode = 0;
				buildCodeThreshold(t, scheme->suffixTree->root, 0, 1., model, scheme);
//printLengthDistribution(stdout, scheme->lengthCode,scheme->cardCode);
				dec = getDecodedFromScheme(scheme);
//printf("cardinal dec = %ld\n", dec->cardinal);
				distB = computeWholeDistanceDec(dec);
				corr = computeCorrelation(distA, distB);
				monfree((void*)distB.table);
				sc = score(dec);
				printf("%lE\t%lf\t%.2lf\n", t, corr, sc);
				fprintf(fo, "%lE\t%lf\t%.2lf\n", t, corr, sc);
				for(n=0; n<dec->number; n++)
					monfree((void*) dec->sequence[n]);
				monfree((void*) dec->sequence);
				monfree((void*) dec->size);
				monfree((void*) dec);
			}
		}
		fprintf(stdout, "\n\n%.4lE\n\n", findMode(set, qmin, qmax, qprec, scheme, model));
		freeModel(model);
	} else {
		for(l = lmax; l>=1; l--) {
			double corr;
			TypeSetOfSequences *dec;
			scheme->cardCode = 0;
			buildCodeLength(l, scheme->suffixTree->root, 0, scheme);
//printLengthDistribution(stdout, scheme->lengthCode,scheme->cardCode);
			dec = getDecodedFromScheme(scheme);
//printf("cardinal dec = %ld\n", dec->cardinal);
			distB = computeWholeDistanceDec(dec);
			corr = computeCorrelation(distA, distB);
			monfree((void*)distB.table);
			fprintf(fo, "%ld\t%lf\n", l, corr);
			fprintf(stdout, "%ld\t%lf\n", l, corr);
			for(n=0; n<dec->number; n++)
				monfree((void*) dec->sequence[n]);
			monfree((void*) dec->sequence);
			monfree((void*) dec->size);
			monfree((void*) dec);
		}
	}
		
	freeScheme(scheme);
	monfree((void*)distA.table);
	fprintf(stdout, "\n\n%ld\n\n", totalLength(*set));
	monfree((void*)set->size);
	for(n=0; n<set->number; n++)
		monfree((void*)set->sequence[n]);
	monfree((void*)set->sequence);
	monfree((void*)set);
	fclose(fo);
/*	sprintf(bufferOutput, "%s_Ali.nex", outputFileName);
	if(!(fo = fopen(bufferOutput, "w")))
		exitProg(ErrorWriting, bufferOutput);
	printDistanceNexus(fo, distA);
	fclose(fo);
	sprintf(bufferOutput, "%s_New.nex", outputFileName);
	if(!(fo = fopen(bufferOutput, "w")))
		exitProg(ErrorWriting, bufferOutput);
	printDistanceNexus(fo, distB);
	fclose(fo);
*/
;
	exitProg(ExitOk,NULL);
	return 0;
}

double findMode(TypeSetOfSequences *set, double tmin, double tmax, double tprec, TypeCodeScheme *scheme, TypeMarkovModel *model) {
	double tmid, t, smax, tres, sc,  scl, scr;
	TypeNumber n;
//	TypeCodeScheme *scheme;
	TypeSetOfSequences *dec;
	
/*	scheme = (TypeCodeScheme*) monmalloc(sizeof(TypeCodeScheme));
	scheme->suffixTree = getSuffixTree(set);
	scheme->code = (TypePosition*) monmalloc(scheme->suffixTree->size*sizeof(TypePosition));
	scheme->buffSize = INC_SIZE_CODE;
	scheme->lengthCode = (TypePosition*) monmalloc(scheme->buffSize*sizeof(TypePosition));
*/	while((tmax-tmin)>4*tprec) {
		tmid = (tmax+tmin)/2.;
		scheme->cardCode = 0;
		buildCodeThreshold(exp(tmid-3.*tprec/2.), scheme->suffixTree->root, 0, 1., model, scheme);
		dec = getDecodedFromScheme(scheme);
		scl = score(dec);
		for(n=0; n<dec->number; n++)
			monfree((void*) dec->sequence[n]);
		monfree((void*) dec->sequence);
		monfree((void*) dec->size);
		monfree((void*) dec);
		scheme->cardCode = 0;
		buildCodeThreshold(exp(tmid+3.*tprec/2.), scheme->suffixTree->root, 0, 1., model, scheme);
		dec = getDecodedFromScheme(scheme);
		scr = score(dec);
		for(n=0; n<dec->number; n++)
			monfree((void*) dec->sequence[n]);
		monfree((void*) dec->sequence);
		monfree((void*) dec->size);
		monfree((void*) dec);
		if(scl>scr)
			tmax = tmid+3.*tprec/2.;
		else
			tmin = tmid-3.*tprec/2.;
	}
	if(scl>scr) {
		smax = scl;
		tres = exp(tmid-3.*tprec/2.);
	} else {
		smax = scr;
		tres = exp(tmid+3.*tprec/2.);
	}
	scheme->cardCode = 0;
	buildCodeThreshold(exp(tmid), scheme->suffixTree->root, 0, 1., model, scheme);
	dec = getDecodedFromScheme(scheme);
	sc = score(dec);
	for(n=0; n<dec->number; n++)
		monfree((void*) dec->sequence[n]);
	monfree((void*) dec->sequence);
	monfree((void*) dec->size);
	monfree((void*) dec);
	if(sc>smax) {
		smax = scl;
		tres = exp(tmid);
	}
//	freeScheme(scheme);
	return tres;
}


double scorebis(TypeSetOfSequences *dec) {
	TypeNumber n, *last, *present;
	TypePosition p, *total;
	double score;

	TypeSymbol c;

	last =(TypeNumber*) monmalloc(dec->cardinal*sizeof(TypeNumber));
	present =(TypeNumber*) monmalloc(dec->cardinal*sizeof(TypeNumber));
	total =(TypePosition*) monmalloc(dec->cardinal*sizeof(TypePosition));
	for(c=0; c<dec->cardinal; c++) {
		present[c] = 0;
		total[c] = 0;
		last[c] = -1;
	}
	for(n=0; n<dec->number; n++) {
		for(p=0; p<dec->size[n]; p++)
			if(dec->sequence[n][p]<dec->cardinal) {
				if(last[dec->sequence[n][p]] != n) {
					last[dec->sequence[n][p]] = n;
					present[dec->sequence[n][p]]++;
				}
				total[dec->sequence[n][p]]++;
			}
	}
	score = 0;
	for(c=0; c<dec->cardinal; c++)
		if(total[c]>1 && total[c] <3* present[c])
			score += (pow(total[c], 0.5));
//			score += pow(((double)present[c])/((double)total[c]), 20.);
//	score /= dec->cardinal;
	monfree((void*)last);
	monfree((void*)total);
	monfree((void*)present);
	return score;
}




double score(TypeSetOfSequences *dec) {
	TypeNumber n, *present, *total;
	TypePosition p;
	double score;
	int *out;
	TypeSymbol c;

	present =(TypeNumber*) monmalloc(dec->cardinal*sizeof(TypeNumber));
	total =(TypeNumber*) monmalloc(dec->cardinal*sizeof(TypeNumber));
	out =(int*) monmalloc(dec->cardinal*sizeof(int));
	for(c=0; c<dec->cardinal; c++) {
		present[c] = -1;
		total[c] = 0;
		out[c] = 0;
	}
	for(n=0; n<dec->number; n++) {
		for(p=0; p<dec->size[n]; p++)
			if(dec->sequence[n][p]<dec->cardinal && !out[dec->sequence[n][p]]) {
				if(present[dec->sequence[n][p]] == n)
					out[dec->sequence[n][p]] = 1;
				else {
					present[dec->sequence[n][p]] = n;
					total[dec->sequence[n][p]]++;
				}
			}
	}
	score = 0;
	for(c=0; c<dec->cardinal; c++)
		if(!out[c])
//			score += (total[c]*(total[c]-1))/2;
			score += (pow(total[c], 1)*(total[c]-1));
	monfree((void*)out);
	monfree((void*)total);
	monfree((void*)present);
	return score;
}