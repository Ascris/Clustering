#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "Sequences.h"
#include "SuffixTree.h"
#include "PrefixCode.h"
#include "PrefixDecode.h"
#include "Utils.h"


#define SIZE_BUFFER_CHAR 300
#define MAX_FUNC 10

#define HELPMESSAGE "\nusage: decauto [options] <input file> <output file>\n\nThe input file has to be in Fasta format.\n\nOptions:\n\n-h                  output help\n\n-s <d, r, p>        select the type of alphabet (d: DNA, r: RNA, p: Protein)\n                    maybe useful to handle ambiguity characters\n"
#define AMBIMESSAGE "The ambiguity codes are for DNA and RNA:\n                                         * Y: T/U or C\n                                         * R: A or G\n                                         * M: A or C\n                                         * K: T/U or G\n                                         * W: T/U or A\n                                         * S: C or G\n                                         * B: T/U, C, or G\n                                         * D: T/U, A, or G\n                                         * H: T/U, C, or A\n                                         * V: C, A, or G\n                                         * N: T/U, C, A, or G\n\nand for protein:\n                                         * B: D or N\n"
static double score(TypeSetOfSequences *dec);

int main(int argc, char **argv) {		
	TypePosition orderstart=1, orderend=10, length = 10;
	char option[256], inputFileName[SIZE_BUFFER_CHAR], outputFileName[SIZE_BUFFER_CHAR], bufferOutput[SIZE_BUFFER_CHAR], *table, 
	outputFormat = 'r', typeDec = 'l', typeAlphabet = 'd', typeCalc = 'g';
	TypeSetOfSequences set;
	int fixed = 0, flagThre = 1;
	double threshold = 0.001, tmin = -25, tmax = -3, tprec = 0.5;
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
			flagThre = 1;
		}
		if(option['l']) {
			option['l'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%ld", &length) == 1)
				i++;
			flagThre = 0;
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
		set = readSequencesFasta(fi, table, typeAlphabet == '?');
		switch(typeAlphabet) {
			case 'd':
			case 'r':
				set.ambiguity = getXNAAmbiguity();
				break;
			case 'p':
				set.ambiguity = getProteinAmbiguity();
				break;
			case '?':
			default:
				set.ambiguity.number = 0;
		}
		set.cardinal -= set.ambiguity.number;
		fclose(fi);
	} else {
		exitProg(ErrorReading, inputFileName);
	}
	
	if(fo = fopen(outputFileName, "w")) {
		TypeSetOfSequences *dec;
		double tmid, t, smax, tres, sc,  scl, scr;
		TypeNumber n;
		TypeCodeScheme *scheme;
		TypeMarkovModel *model;

		fixSequencesAmbiguity(&set);
		scheme = (TypeCodeScheme*) monmalloc(sizeof(TypeCodeScheme));
		scheme->suffixTree = getSuffixTree(&set);
		scheme->code = (TypePosition*) monmalloc(scheme->suffixTree->size*sizeof(TypePosition));
		scheme->buffSize = INC_SIZE_CODE;
		scheme->lengthCode = (TypePosition*) monmalloc(scheme->buffSize*sizeof(TypePosition));
		model = estimateMarkovModel(&set);
		while((tmax-tmin)>4*tprec) {
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
printf("%.4lE\t%lf\n", tres, smax);
		scheme->cardCode = 0;
		buildCodeThreshold(tres, scheme->suffixTree->root, 0, 1., model, scheme);
		dec = getDecodedFromScheme(scheme);
		freeModel(model);
		freeScheme(scheme);
		printDecodedSequencesFastaAlt(fo, &set, dec, 50);
		fclose(fo);
	} else {
		exitProg(ErrorWriting, outputFileName);
	}
	exitProg(ExitOk,NULL);
	return 0;
}
/*		switch(typeCalc) {
			case 'p':
				dist = computeWholeDistancePair(set, computeAlex);
				break;
			case 'g':
			default:
				dist = computeWholeDistanceGlobal(&set);
		}
*/


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