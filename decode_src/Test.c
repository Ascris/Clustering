#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "Sequences.h"
#include "MarkovEmbed.h"
#include "SuffixTree.h"
#include "Utils.h"


#define SIZE_BUFFER_CHAR 300
#define MAX_FUNC 9

#define HELPMESSAGE "\nusage: decode [options] <order> <input file> <output file>\n\nOrder has to be an integer. Input file has to be in Fasta format.\nOutput is in Fasta format.\n\nOptions:\n\n-h                  output help\n\n-s <d, r, p>        select the type of alphabet (d: DNA, r: RNA, p: Protein)\n                    maybe useful to handle ambiguity characters\n\n-d <l, s>           select the type of decoding:\n                     * l: use local decoding of order <order>\n                     * s: use sequence of sliding blocks  of length <order>\n"
#define AMBIMESSAGE "The ambiguity codes are for DNA and RNA:\n                                         * Y: T/U or C\n                                         * R: A or G\n                                         * M: A or C\n                                         * K: T/U or G\n                                         * W: T/U or A\n                                         * S: C or G\n                                         * B: T/U, C, or G\n                                         * D: T/U, A, or G\n                                         * H: T/U, C, or A\n                                         * V: C, A, or G\n                                         * N: T/U, C, A, or G\n\nand for protein:\n                                         * B: D or N\n"

int main(int argc, char **argv) {		
	TypePosition order, min;
	char option[256], inputFileName[SIZE_BUFFER_CHAR], outputFileName[SIZE_BUFFER_CHAR], bufferOutput[SIZE_BUFFER_CHAR], *table, 
	outputFormat = 'r', typeOut = 'c', typeDec = 'l', typeAlphabet = '?', typeCalc = 'p';
	TypeSetOfSequences set;
	double threshold = 0.01;
	FILE *fi, *fo;
	int i = 1;

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
		if(option['d']) {
			option['d'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &typeDec) == 1)
				i++;
		}
		if(option['m']) {
			option['m'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%ls", &min) == 1)
				i++;
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
//	if (i>=argc || sscanf(argv[i++], "%d", &order) != 1) exitProg(ErrorArgument, HELPMESSAGE);
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
//		printf("mean length %f\n", meanLength(set),  log(meanLength(set))/log(4.0));
//printSequencesFasta(stdout, set, 50);
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
printf("MUMMER\n");
getMUM(&set, min);
exitProg(ExitOk,"bye bye");

{
TypePosition **bord, p, size = 10, st = 5, en = 15;
TypeSymbol w[20], *wo, a, b, card = 2;
TypeMarkovModel *m;
clock_t start, end;

printf("sizeof(short) = %d\nsizeof(int) = %d\nsizeof(long) = %d\nsizeof(void*) = %d\nsizeof(TypeSuffixNode) = %d\n size info = %d\nMAXINT %ld\nMAXSHORT %d\n", sizeof(short), sizeof(int), sizeof(long), sizeof(void*), sizeof(TypeSuffixNode), INT_MAX, SHRT_MAX);
for(p=0; p<20; p++)
putc('.',stdout);
putc('\n',stdout);
m = estimateMarkovModel(&set);
qsort(set.size, set.number, sizeof(TypePosition), comparePosition);
getBounds(&st, &en, 2, &set, set.size, 0.001, m);

exitProg(ExitOk,"bye bye");
double pr;
/*
w[0] = 0;
w[1] = 1;
w[2] = 0;
w[3] = 0;
w[4] = 1;
w[5] = 0;
w[6] = 1;
w[7] = 1;
w[8] = 1;
w[9] = 1;

bord = getNext(w, size, card);
printf("[P]");
for(p=0; p<size; p++)
	printf("\t%d",p);
printf("\n");
printf("[S]");
for(p=0; p<size; p++)
	printf("\t%d",w[p]);
printf("\n");
for(a=0; a<card; a++){
	printf("[%d]", a);
	for(p=0; p<size; p++)
		printf("\t%d",bord[a][p]);
	printf("\n");
}*/
for(a=0; a<set.cardinal; a++) {
	printf("%lf", m->init[a]);
	for(b=0; b<set.cardinal; b++)
		printf("\t%lf", m->trans[a][b]);
	printf("\n");
}

size = 9;
wo = &(set.sequence[7][189]);
bord = getNext(wo, size, set.cardinal);
printf("[P]");
for(p=0; p<size; p++)
	printf("\t%d",p);
printf("\n");
printf("[S]");
for(p=0; p<size; p++)
	printf("\t%d",wo[p]);
printf("\n");
for(a=0; a<set.cardinal; a++){
	printf("[%d]", a);
	for(p=0; p<size; p++)
		printf("\t%d",bord[a][p]);
	printf("\n");
}
start = clock();
pr = getProb(wo, size, 2, &set, m);
end = clock();
printf("%d %d prob %lE time %.2lEs\n", 63, size, pr, ((double) end-start)/((double)CLOCKS_PER_SEC));
printf("\n");
start = clock();
pr = getProbMat(wo, size, 2, &set, m);
end = clock();
printf("%d %d prob %lE time %.2lEs\n", 63, size, pr, ((double) end-start)/((double)CLOCKS_PER_SEC));
printf("\n");
//	exitProg(ExitOk,"bye bye");
/*
printf("%d %d prob %lE\n", 63, 12, getProb(&(set.sequence[0][63]), 12, 2, &set, m));
printf("%d %d prob %lE\n", 63, 11, getProb(&(set.sequence[0][63]), 11, 2, &set, m));
printf("%d %d prob %lE\n", 64, 11, getProb(&(set.sequence[0][64]), 11, 2, &set, m));

size = 12;
wo = &(set.sequence[0][63]);
bord = getNext(wo, size, set.cardinal);
printf("[P]");
for(p=0; p<size; p++)
	printf("\t%d",p);
printf("\n");
printf("[S]");
for(p=0; p<size; p++)
	printf("\t%d",wo[p]);
printf("\n");
for(a=0; a<set.cardinal; a++){
	printf("[%d]", a);
	for(p=0; p<size; p++)
		printf("\t%d",bord[a][p]);
	printf("\n");
}
pr = getProb(wo, size, 2, &set, m);
printf("%d %d prob %lE\n", 63, size, pr);
printf("\n");

size = 11;
wo = &(set.sequence[0][64]);
bord = getNext(wo, size, set.cardinal);
printf("[P]");
for(p=0; p<size; p++)
	printf("\t%d",p);
printf("\n");
printf("[S]");
for(p=0; p<size; p++)
	printf("\t%d",wo[p]);
printf("\n");
for(a=0; a<set.cardinal; a++){
	printf("[%d]", a);
	for(p=0; p<size; p++)
		printf("\t%d",bord[a][p]);
	printf("\n");
}
pr = getProb(wo, size, 2, &set, m);
printf("%d %d prob %lE\n", 64, size, pr);
//printf("prob %lE\n", getProb(w, size, 2, &set, m));
	exitProg(ExitOk,"bye bye");
*/}
	if(fo = fopen(outputFileName, "w")) {
/*double del;
clock_t start = clock();
printf("ok read start\n");
computeSuffixTree(set);
del = ((double)clock() - start) / CLOCKS_PER_SEC;
printf("Time elapsed: %f ratio = %f Mo/s\n", del, 1048576.0*del/((double)set.size[0]));
*///exitProg(ExitOk,"Bye");
		setThreshold(threshold);
		switch(typeDec) {
			case 's':
//				printSuffixTree(fo, computeSuffixTree(set));
//				break;
			case 'l':
			default:
//				printDecodedSequencesFasta(stdout, &set, computeCoded(&set), 50);
				printDecodedSequencesFasta(fo, &set, computeCoded(&set), 50);
//				printSuffixTree(fo, computeSuffixTree(&set));
		}	
		fclose(fo);
	} else {
		exitProg(ErrorWriting, outputFileName);
	}
	exitProg(ExitOk,"bye bye");
	return 0;
}
