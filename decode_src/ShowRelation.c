#include <stdlib.h>
#include <string.h>
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

int main(int argc, char **argv) {
	TypePosition orderstart=1, orderend=10;
	char option[256], inputFileName[SIZE_BUFFER_CHAR], outputFileName[SIZE_BUFFER_CHAR], bufferOutput[SIZE_BUFFER_CHAR], *table, 
	outputFormat = 'r', typeDec = 'l', typeAlphabet = '?', typeCalc = 'g';
	TypeSetOfSequences *set;
	TypeAlignment aln, atmp;
	int fixed = 0;
	double threshold = 0.001, tmin = 0.000000001, tmax=0.1, tstep = 0.00001;
	double thre;
	TypeDistance distA, distB;
	TypePosition l;
	int kmin = -8;
	double *tt;
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
			if((i+1)<argc && sscanf(argv[i+1], "%d", &typeDist) == 1)
				i++;
			if(typeDist >= MAX_FUNC)
				typeDist = 0;
		}
		if(option['t']) {
			option['t'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &threshold) == 1)
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
	if(!(fo = fopen(outputFileName, "w")))
		exitProg(ErrorWriting, outputFileName);
	atmp.number = aln.number;
	atmp.name = aln.name;
	atmp.table = aln.table;
	atmp.title = aln.title;
	atmp.ambiguity = aln.ambiguity;
	atmp.empty = aln.empty;
	atmp.cardinal = aln.cardinal;
	atmp.comment = aln.comment;
	atmp.sequence = (TypeSymbol**) monmalloc(atmp.number*sizeof(TypeSymbol*));
	tt = (double*) monmalloc(10*(-kmin)*sizeof(double));
printf("size %ld\n", aln.size);
	for(l=500; l<aln.size; l *= 10){
		TypePosition h, tot;
		double max, best, t;
		atmp.size = l;
		for(h=0; h<aln.size-l; h += l) {
			TypeNumber n;
			int k, ki, ks, km;
			for(n=0; n<aln.number; n++)
				atmp.sequence[n] = &(aln.sequence[n][h]);
			distA = computeWholeDistancePairAln(atmp, computeNorm1Aln);
			set=toSequences(&atmp);
			if(set->number == aln.number) {
				double t;
				int j;
				TypeMarkovModel *model;
				TypeCodeScheme *scheme;
				model = estimateMarkovModel(set);
				scheme = (TypeCodeScheme*) monmalloc(sizeof(TypeCodeScheme));
				scheme->suffixTree = getSuffixTree(set);
				scheme->code = (TypePosition*) monmalloc(scheme->suffixTree->size*sizeof(TypePosition));
				scheme->buffSize = INC_SIZE_CODE;
				scheme->lengthCode = (TypePosition*) monmalloc(scheme->buffSize*sizeof(TypePosition));
				max = 0.;
				for(k=kmin; k<0; k++) {
					for(j=0; j<9; j++) {
						double corr;
						TypeSetOfSequences *dec;
						t = (j+1)*pow(10.,(double)k);
						scheme->cardCode = 0;
						buildCodeThreshold(t, scheme->suffixTree->root, 0, 1., model, scheme);
						dec = getDecodedFromScheme(scheme);
						distB = computeWholeDistanceDec(dec);
						corr = computeCorrelation(distA, distB);
						tt[-k*10-j-1] = corr;
						monfree((void*)distB.table);
//printf("%lE\t%lf (%d-%d)\n", t, corr, k, j);
						for(n=0; n<dec->number; n++)
							monfree((void*) dec->sequence[n]);
						monfree((void*) dec->sequence);
						monfree((void*) dec->size);
						monfree((void*) dec);
						if(corr>max) {
							max = corr;
							best = t;
							km = -k*10-j-1;
						}
					}
				}
				freeModel(model);
				freeScheme(scheme);
				for(ki=km; ki>=0 && max-tt[ki]<=(1-max)/10.; ki--);
				for(ks=km; ks<10*(-kmin) && max-tt[ks]<=(1-max)/10.; ks++);
				fprintf(fo, "%ld\t%lE\t%lE\t%lE\n", tot, best, (10-ki%10)*pow(10., -1-ki/10), (10-ks%10)*pow(10., -1-ks/10));
				fprintf(stdout, "%ld\t%lE\t%lE\t%lE\n", tot, best, (10-ki%10)*pow(10., -1-ki/10), (10-ks%10)*pow(10., -1-ks/10));
			}
			monfree((void*)distA.table);
			monfree((void*)set->size);
			for(n=0; n<set->number; n++)
				monfree((void*)set->sequence[n]);
			monfree((void*)set->sequence);
			monfree((void*)set);
		}
	}
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

	exitProg(ExitOk,NULL);
	return 0;
}
/*for(k = 0; k<10*(-kmin); k++)
	printf("%d\t%lE\t%lE\n", k, (10-k%10)*pow(10., -1-k/10), tt[k]);
printf("\n max %d\t%lE\n", km, tt[km]);
printf("\n bounds %d\t%d\n", ki, ks);
*/
//					fprintf(fo, "%lE\t%lf\t%lf\n", t, corr, ((double)dec->cardinal)/((double)tot));
