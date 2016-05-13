#ifndef UtilsF
#define UtilsF
#include <stdlib.h>
#include <stdio.h>

typedef enum TS {
	ErrorArgument=0,
	ErrorInit,
	ErrorReading,
	ErrorWriting,
	ErrorMemory,
	ErrorExec,
	ExitOk
} TypeExit;

typedef struct REGRESSION {
	double a, b, r;
} TypeCoeffReg;

void printPreambleNexus(FILE *f, char **name, int number);
void exitProg(TypeExit code, char *message);
void *monmalloc(long size);
void *monrealloc(void *in, long size);
int IsSeparator(char c);
int IsItemSeparator(char c);
int IsLineSeparator(char c);
int readLine(FILE *f, char *buffer);
int readItem(FILE *f, char *buffer);
int tokenize(char *src, char **dest);
int find(char *src, char **dest, int size);
void fixSpace(char *src);
char *truncFileName(char* name);
char *getExtension(char* name);
double log2(double x);
void initProgress(char *message, long max);
void updateProgress(long val);
void setProgress(int flag);
void flushProgress();

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define ABS(x) ((x)<(0)?(-x):(x))
#define NEG_INFTY -9.999e99
#define POS_INFTY 9.999e99
#define SPECIAL -1
//#define DEBUG -1
#endif
