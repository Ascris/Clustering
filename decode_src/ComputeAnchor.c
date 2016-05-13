#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "Alignment.h"
#include "Sequences.h"
#include "Distance.h"
#include "SuffixTree.h"
#include "PrefixDecode.h"
#include "Utils.h"
#include "DistSeq.h"


#define SIZE_BUFFER_CHAR 300
#define MAX_FUNC 10
#define SPECIAL -1

#define HELPMESSAGE "\nusage: anchor [options] <input file> <output file>\n\nThe input file has to be in Fasta format.\n\nreturn a table when each line is for a pair fo sequences in which:\n column uni counts the number of decoded symbols ocurring once in each sequence\n- column homol counts the number of columns equal\n-column both count the number of column with a same (unique) decoded symbols\n\n"

static TypeSetOfSequences *computeAnchor(FILE *f, TypeAlignment *aln);
static void printDecodedAlignment(FILE *f, TypeAlignment *aln, TypeSetOfSequences *dec);

TypeSetOfSequences *computeAnchor(FILE *f, TypeAlignment *aln) {
	TypeNumber i, j;
	TypeSetOfSequences *dec, *set;
	TypePosition o, p, *occi, *occj, li, lj;
	TypeNumber n;
	TypeSymbol c, *listi, *listj, ni, nj;

	dec = getCoded(set=toSequences(aln));
	monfree((void*) set->sequence);
	monfree((void*) set->size);
	monfree((void*) set);
	listi =(TypePosition*) monmalloc(dec->cardinal*sizeof(TypeSymbol));
	listj =(TypePosition*) monmalloc(dec->cardinal*sizeof(TypeSymbol));
	occi =(TypePosition*) monmalloc(dec->cardinal*sizeof(TypePosition));
	occj =(TypePosition*) monmalloc(dec->cardinal*sizeof(TypePosition));
	for(c=0; c<dec->cardinal; c++) {
		occi[c] = 0;
		occj[c] = 0;
	}
	fprintf(f, "seq1\tseq2\tsize1\tsize2\tsizemin\tunique\thomol\tboth\t% homol\t% both among homol\t% both among unique\n");
	for(i=1; i<dec->number; i++) {
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
			TypePosition indi, indj, sum, uni, homo;
			nj=0;
			lj=0;
			for(p=0; p<dec->size[j]; p++)
				if(dec->sequence[j][p]<dec->cardinal) {
					if(occj[dec->sequence[j][p]] == 0)
						listj[nj++] = dec->sequence[j][p];
					lj++;
					occj[dec->sequence[j][p]]++;
				}
			sum = 0;
			uni = 0;
			homo = 0;
			indi = 0;
			indj = 0;

			for(c=0; c<nj; c++)
				if(occj[listj[c]] == 1 && occi[listj[c]] == 1)
					uni++;
			for(p=0; p<aln->size; p++) {
				if(aln->sequence[i][p] == aln->sequence[j][p] && aln->sequence[i][p] != aln->empty)
					homo++;
				if(aln->sequence[i][p] == aln->sequence[j][p] && aln->sequence[i][p] != aln->empty && dec->sequence[i][indi]==dec->sequence[j][indj] && occi[dec->sequence[i][indi]] == 1 && occj[dec->sequence[i][indi]] == 1)
					sum++;
				if(aln->sequence[i][p] != aln->empty)
					indi++;
				if(aln->sequence[j][p] != aln->empty)
					indj++;
			}
			for(c=0; c<nj; c++)
					occj[listj[c]] = 0;
			fprintf(f, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\n", i,j,dec->size[i], dec->size[j], MIN(dec->size[i], dec->size[j]), uni, homo, sum, (100.*homo)/((double)MIN(dec->size[i], dec->size[j])), ((double)100.*sum)/((double)homo), ((double)100.*sum)/((double)uni));
		}
		for(c=0; c<ni; c++)
			occi[listi[c]] = 0;
	}
	monfree((void*) occi);
	monfree((void*) occj);
	monfree((void*) listi);
	monfree((void*) listj);
	return dec;
}

void printDecodedAlignment(FILE *f, TypeAlignment *aln, TypeSetOfSequences *dec) {
	TypeNumber n;
	TypeSymbol *nIdent, *sIdent, s;
	TypePosition *occ;

	occ = (TypePosition*) monmalloc(dec->cardinal*sizeof(TypePosition));
	for(s=0; s<dec->cardinal; s++)
		occ[s] = 0;
	for(n=0; n<dec->number; n++) {
		TypePosition position;
		for(position=0; position<dec->size[n]; position ++)
			occ[dec->sequence[n][position]]++;
	}
	nIdent = (TypeSymbol*) monmalloc(dec->cardinal*sizeof(TypeSymbol));
	for(s=0; s<dec->cardinal; s++)
		nIdent[s] = SPECIAL;
	sIdent = (TypeSymbol*) monmalloc((aln->cardinal+aln->ambiguity.number)*sizeof(TypeSymbol));
	for(s=0; s<aln->cardinal; s++)
		sIdent[s] = 0;
	for(n=0; n<aln->number; n++) {
		TypePosition position, ind;
		if(aln->name != NULL)
			fprintf(f, "%s", aln->name[n]);
		else
			fprintf(f, "Sequence %d", n+1);
		ind = 0;
		if(aln->table != NULL) {
			for(position=0; position<aln->size; position ++) {
				if(aln->sequence[n][position] != aln->empty) {
					if(occ[dec->sequence[n][ind]]>1) {
						if(nIdent[dec->sequence[n][ind]] == SPECIAL)
							nIdent[dec->sequence[n][ind]] = sIdent[aln->sequence[n][position]]++;
						fprintf(f, "\t%c%ld", aln->table[aln->sequence[n][position]], nIdent[dec->sequence[n][ind]]);
					} else
						fprintf(f, " \t%c", tolower(aln->table[aln->sequence[n][position]]));
					ind++;
				}
			}
		} else {
			for(position=0; position<aln->size; position ++) {
				if(aln->sequence[n][position] != aln->empty) {
					if(nIdent[dec->sequence[n][ind]] == SPECIAL)
						nIdent[dec->sequence[n][ind]] = sIdent[aln->sequence[n][position]]++;
					fprintf(f, "\t%ld-%ld", aln->sequence[n][position], nIdent[dec->sequence[n][ind]]);
					ind++;
				}
			}
		}
		fprintf(f, "\n");
	}
	monfree((void*) nIdent);
	monfree((void*) sIdent);
}

int main(int argc, char **argv) {	
	TypePosition orderstart=1, orderend=10;
	char option[256], inputFileName[SIZE_BUFFER_CHAR], outputFileName[SIZE_BUFFER_CHAR], bufferOutput[SIZE_BUFFER_CHAR], *table, 
	outputFormat = 'r', typeDec = 'l', typeAlphabet = '?', typeCalc = 'g';
	TypeSetOfSequences *set, *dec;
	TypeAlignment aln;
	int fixed = 0;
	double threshold = 0.001;
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
	if(!(fo = fopen(outputFileName, "w")))
		exitProg(ErrorWriting, outputFileName);
	setThreshold(threshold);
	dec = computeAnchor(fo, &aln);
	fclose(fo);
	sprintf(bufferOutput, "%s_ali.csv", outputFileName);
	if(!(fo = fopen(bufferOutput, "w")))
		exitProg(ErrorWriting, bufferOutput);
	printDecodedAlignment(fo, &aln, dec);
	fclose(fo);
	
	exitProg(ExitOk,NULL);
	return 0;
}
