#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "Utils.h"

#define PROGRESS_SIZE 50

static long progress_max;
static int progress_cur, progress_flag = 0;

void setProgress(int flag) {
	progress_flag = flag;
}

void initProgress(char *message, long max) {
	int i;
if(!progress_flag)
	return;
	progress_max = max;
	progress_cur = 0;
	printf("\nComputing %s\n", message);
	for(i=0; i<=PROGRESS_SIZE; i++)
		printf("_");
	printf("\n");
}

void updateProgress(long val) {
	int m;
if(!progress_flag)
	return;
	m = (PROGRESS_SIZE*val)/progress_max;
//printf("%ld %ld\n", val, progress_max);
	for(; progress_cur<=m; progress_cur++)
		fputc('*', stdout);
	fflush(stdout);
}

void flushProgress() {
if(!progress_flag)
	return;
	for(; progress_cur<=PROGRESS_SIZE; progress_cur++)
		fputc('*', stdout);
	fflush(stdout);
	printf("\n");
}

void printPreambleNexus(FILE *f, char **name, int number) {
	int i;
	fprintf(f, "#NEXUS\n\n");
	fprintf(f, "BEGIN taxa;\n");
	fprintf(f, "    DIMENSIONS ntax=%ld;\n", number);
	fprintf(f, "TAXLABELS\n");
	for(i=0; i<number; i++)
		fprintf(f, "'%s'\n", name[i]);
	fprintf(f, ";\n");
	fprintf(f, "END;\n");
}

double log2(double x) {
	return log(x)/log(2);
}

void exitProg(TypeExit code, char *message) {
	switch(code)
	{
		case ErrorArgument:
			printf("Error when reading arguments\n");
			break;
		case ErrorInit:
			printf("Problem of initialization\n");
			break;
		case ErrorReading:
			printf("Problem when reading a file\n");
			break;
		case ErrorWriting:
			printf("Problem when writing a file\n");
			break;
		case ErrorMemory:
			printf("Not enougth memory\n");
			break;
		case ErrorExec:
			printf("Problem during execution...\n");
			break;
	}
	if(message != NULL)
		printf("%s\n", message);
	if(code == ExitOk)
		exit(EXIT_SUCCESS);
	else
		exit(EXIT_FAILURE);
}

/****************************************************/
/************** Allocations ********************/
/****************************************************/
void *monmalloc(long size) {
	void *point;
	if(size<0 || !(point = malloc(size))) {
		char tmp[200];
		sprintf(tmp, "Try to allocate %ld bytes\n", size);
		exitProg(ErrorMemory, tmp);
	}
	return point;
} 

void *monrealloc(void *in, long size) {
	void *point;
/*	if(size<0 || !(point = realloc(in, size))) {
		char tmp[200];
		sprintf(tmp, "Try to reallocate %ld bytes\n", size);
		exitProg(ErrorMemory, tmp);
	}
*/	point = realloc(in, size);
	return point;
}

void monfree(void *in) {
	if(in != 0L)
		free(in);
}
int IsSeparator(char c) {
	return (c == ' ' || c == '\t' || c == '\n' || c == '\r');
}
int IsItemSeparator(char c) {
	return (c == ' ' || c == '\t' || c == ';');
}
int IsLineSeparator(char c) {
	return (c == '\n' || c == '\r');
}

int readLine(FILE *f, char *buffer) {
	char c;
	int tot = 0;
	while((c=getc(f))!=EOF && IsLineSeparator(c));
	while(c!=EOF && !IsLineSeparator(c)) {
		buffer[tot++] = c;
		c=getc(f);
	}
	if(c!=EOF)
		buffer[tot++] = '\n';
	buffer[tot] = '\0';
	return tot;
}

int readItem(FILE *f, char *buffer) {
	char c;
	int tot = 0;
	while((c=getc(f))!=EOF && IsSeparator(c));
	while(c!=EOF && !IsSeparator(c)) {
		buffer[tot++] = c;
		c=getc(f);
	}
	buffer[tot++] = '\n';
	buffer[tot] = '\0';
	return tot;
}

char *truncFileName(char* name) {
	int i;
	for(i=strlen(name)-1; i>=0 && name[i] != '.'; i--);
	if(i>=0)
	 name[i] = '\0';
	return name;
}

char *getExtension(char* name) {
	int i, length;
	
	length = strlen(name);
	for(i=length-1; i>0 && name[i] != '.'; i--);
	if(i>0)
		return &(name[i+1]);
	else
		return &(name[length]);
}

int tokenize(char *src, char **dest) {
	int indice = 0, pos = 0;
	do {
		while(src[pos] != '\0' && IsSeparator(src[pos]))
			 pos++;
		if(src[pos] != '\0') {
			dest[indice++] = &(src[pos]);
		}
		while(src[pos] != '\0' && !IsSeparator(src[pos]))
			 pos++;
		if(src[pos] != '\0') {
			src[pos++] = '\0';
		}
	} while(src[pos] != '\0');
	dest[indice] = NULL;
	return indice;
}

int find(char *src, char **dest, int size) {
	int i;
	for(i=0; i<size && strcmp(src, dest[i]) != 0; i++);
	return i;
}

void fixSpace(char *src) {
	int i;
	for(i=0; src[i] != '\0'; i++)
		if(IsSeparator(src[i]))
			src[i] = '_';
}
