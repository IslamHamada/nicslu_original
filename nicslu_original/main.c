#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nicslu.h"
#include "nicslu_util.h"
void printLU2(SNicsLU* nicslu, int i) {
	uint__t* llen = nicslu->llen;
	uint__t* ulen = nicslu->ulen;
	unsigned char* arr = (unsigned char*)(((uint__t *)nicslu->lu_array2[i]) + ulen[i]);
	printf("%d: ", i);
	for (int j = 0; j < nicslu->ulen[i]; j++) {
		printf("%d ", ((unsigned int*)arr)[j]);
	}

	for (int j = 0; j < ulen[i]; j++) {
		printf("%1.3f ", ((double*)(arr + 4 * ulen[i]))[j]);
	}
	printf("\t");
	for (int j = 0; j < llen[i]; j++) {
		printf("%d ", ((unsigned int*)(arr + 12 * ulen[i]))[j]);
	}

	for (int j = 0; j < nicslu->llen[i]; j++) {
		printf("%1.3f ", ((double*)(arr + 12 * ulen[i] + 4 * llen[i]))[j]);
	}

	printf("\n");
}

void fprintLen(SNicsLU* nicslu, int i, FILE* f) {
	uint__t ll = nicslu->llen[i];
	uint__t ul = nicslu->ulen[i];
	fprintf(f, "%d: %d %d\n", i, ul, ll);
}

void fprintLU2(SNicsLU* nicslu, int i, FILE* f) {
	uint__t* llen = nicslu->llen;
	uint__t* ulen = nicslu->ulen;
	unsigned char* arr = (unsigned char*)(((uint__t *)nicslu->lu_array2[i]) + ulen[i]);
	fprintf(f, "%d: ", i);
	for (int j = 0; j < nicslu->ulen[i]; j++) {
		fprintf(f, "%d ", ((unsigned int*)arr)[j]);
	}

	for (int j = 0; j < ulen[i]; j++) {
		fprintf(f, "%f ", ((double*)(arr + 4 * ulen[i]))[j]);
	}
	fprintf(f, "\t");
	for (int j = 0; j < llen[i]; j++) {
		fprintf(f, "%d ", ((unsigned int*)(arr + 12 * ulen[i]))[j]);
	}

	for (int j = 0; j < nicslu->llen[i]; j++) {
		fprintf(f, "%f ", ((double*)(arr + 12 * ulen[i] + 4 * llen[i]))[j]);
	}

	fprintf(f, "\n");
}

void fprintLU(SNicsLU* nicslu, int i, FILE* f) {
	uint__t* llen = nicslu->llen;
	uint__t* ulen = nicslu->ulen;
	unsigned char* arr = (unsigned char*)(((uint__t *)nicslu->lu_array2[i]));
	fprintf(f, "%d: ", i);
	for (int j = 0; j < nicslu->ulen[i]; j++) {
		fprintf(f, "%d ", ((unsigned int*)arr)[j]);
	}

	for (int j = 0; j < ulen[i]; j++) {
		fprintf(f, "%f ", ((double*)(arr + 4 * ulen[i]))[j]);
	}
	fprintf(f, "\t");
	for (int j = 0; j < llen[i]; j++) {
		fprintf(f, "%d ", ((unsigned int*)(arr + 12 * ulen[i]))[j]);
	}

	for (int j = 0; j < nicslu->llen[i]; j++) {
		fprintf(f, "%f ", ((double*)(arr + 12 * ulen[i] + 4 * llen[i]))[j]);
	}

	fprintf(f, "\n");
}

int getUInt(FILE* file, int* num) {
	char c;
	int digits = 0;
	c = fgetc(file);
	int v = 0;
	while (c >= '0'&&c <= '9') {
		v *= 10;
		v += c - '0';
		c = fgetc(file);
	}
	*num = v;
	if (c == '\n') return 1;
	if (c == ',') return 2;
	if (c == ' ') return 3;
	return 0;
}

int getDouble(FILE* file, double* num) {
	char c;
	c = fgetc(file);
	double v1 = 0;
	int sign = 1;

	if (c == '-') {
		sign = -1;
		c = fgetc(file);
	}

	while (c >= '0'&&c <= '9') {
		v1 *= 10;
		v1 += c - '0';
		c = fgetc(file);
	}

	c = fgetc(file);
	double d = 10;

	while (c >= '0'&&c <= '9') {
		v1 += (c - '0') / d;

		d *= 10;
		c = fgetc(file);
	}

	*num = sign * v1;
	if (c == '\n') return 1;
	if (c == ',') return 2;
	if (c == ' ') return 3;
	return 0;
}

int parseCSVInt(FILE* file, uint__t* ap) {
	char c;
	int v = 0;
	int i = 0;
	int out = getUInt(file, &v);
	while (out == 2) {
		ap[i++] = v;
		out = getUInt(file, &v);
	}
	ap[i] = v;
	if (out == 1)return 1;
	return 0;
}

int parseCSVDouble(FILE* file, real__t* ax) {
	char c;
	double v = 0;
	int i = 0;
	int out = getDouble(file, &v);
	while (out == 2) {
		ax[i++] = v;
		out = getDouble(file, &v);
	}
	ax[i] = v;

	if (out == 1)return 1;
	return 0;
}

int main(void)
{
	int ret;
	uint__t n, nnz, i;
	real__t *ax;
	uint__t *ai, *ap;
	SNicsLU *nicslu;
	real__t *x, *b, err;
	ax = NULL;
	ai = NULL;
	ap = NULL;
	x = NULL;
	b = NULL;

	nicslu = (SNicsLU *)malloc(sizeof(SNicsLU));
	NicsLU_Initialize(nicslu);

	//ret = NicsLU_ReadTripletColumnToSparse("ASIC_100k.mtx", &n, &nnz, &ax, &ai, &ap);
	//ret = NicsLU_ReadTripletColumnToSparse("m1.mtx", &n, &nnz, &ax, &ai, &ap);
	FILE * fp;
	fp = fopen("out.csv", "r");
	getUInt(fp, &n);
	getUInt(fp, &nnz);
	fgetc(fp); fgetc(fp); fgetc(fp);
	ap = (uint__t*)malloc(sizeof(uint__t)*(n + 1));
	ai = (uint__t*)malloc(sizeof(uint__t)*nnz);
	ax = (real__t*)malloc(sizeof(real__t)*nnz);
	parseCSVInt(fp, ap);
	fgetc(fp); fgetc(fp); fgetc(fp);
	parseCSVInt(fp, ai);
	fgetc(fp); fgetc(fp); fgetc(fp);
	parseCSVDouble(fp, ax);
	for (i = 0; i<nnz; ++i)
		ax[i] = 1.;
	fclose(fp);
	ret = NICS_OK;
	if (ret != NICS_OK) goto EXIT;

	x = (real__t *)malloc(sizeof(real__t)*(n + n));
	b = x + n;
	for (i = 0; i<n + n; ++i) x[i] = 1.;

	NicsLU_CreateMatrix(nicslu, n, nnz, ax, ai, ap);
	nicslu->cfgf[0] = 1.;

	NicsLU_Analyze(nicslu);
	printf("analysis time: %.8g\n", nicslu->stat[0]);

	NicsLU_Factorize(nicslu);
	printf("factorization time: %.8g\n", nicslu->stat[1]);

	err = _I_NicsLU_CreateETree(nicslu);

	FILE* len = fopen("olenWindowsS.txt", "w");
	FILE* luData = fopen("oluWindowsS.txt", "w");
	for (int ii = 0; ii<n; ii++) {
		int i = nicslu->aeg_data[ii];
		fprintLen(nicslu, i, len);
	}
	for (int ii = 0; ii<n; ii++) {
		int i = nicslu->aeg_data[ii];
		fprintLU(nicslu, i, luData);
	}
	fclose(len);
	fclose(luData);
	//NicsLU_ReFactorize(nicslu, ax);
	//printf("re-factorization time: %.8g\n", nicslu->stat[2]);

	NicsLU_Solve(nicslu, x);
	printf("substitution time: %.8g\n", nicslu->stat[3]);

	NicsLU_Residual(n, ax, ai, ap, x, b, &err, 1, 0);
	printf("Ax-b (1-norm): %.8g\n", err);
	NicsLU_Residual(n, ax, ai, ap, x, b, &err, 2, 0);
	printf("Ax-b (2-norm): %.8g\n", err);
	NicsLU_Residual(n, ax, ai, ap, x, b, &err, 0, 0);
	printf("Ax-b (infinite-norm): %.8g\n", err);
	exit(0);

	printf("NNZ(L+U-I): %ld\n", nicslu->lu_nnz);
	NicsLU_Flops(nicslu, NULL);
	NicsLU_Throughput(nicslu, NULL);
	NicsLU_ConditionNumber(nicslu, NULL);
	printf("flops: %.8g\n", nicslu->stat[5]);
	printf("throughput (bytes): %.8g\n", nicslu->stat[12]);
	printf("condition number: %.8g\n", nicslu->stat[6]);
	NicsLU_MemoryUsage(nicslu, NULL);
	printf("memory (Mbytes): %.8g\n", nicslu->stat[21] / 1024. / 1024.);

EXIT:
	NicsLU_Destroy(nicslu);
	free(ax);
	free(ai);
	free(ap);
	free(nicslu);
	free(x);
	return 0;
}

//int main(int argc, char *argv[]) {
//	int ret;
//	uint__t n, nnz, i;
//	real__t *ax;
//	uint__t *ai, *ap;
//	SNicsLU *nicslu;
//	real__t *x, *b, err;
//
//	if (argc == 1)
//	{
//		printf("usage: demop <#threads>\n");
//		return -1;
//	}
//
//	ax = NULL;
//	ai = NULL;
//	ap = NULL;
//	x = NULL;
//	b = NULL;
//
//	nicslu = (SNicsLU *)malloc(sizeof(SNicsLU));
//	NicsLU_Initialize(nicslu);
//
//	//ret = NicsLU_ReadTripletColumnToSparse("ASIC_100k.mtx", &n, &nnz, &ax, &ai, &ap);
//	FILE * fp;
//	fp = fopen("out.csv", "r");
//	getUInt(fp, &n);
//	getUInt(fp, &nnz);
//	fgetc(fp); fgetc(fp); fgetc(fp);
//	ap = (uint__t*)malloc(sizeof(uint__t)*(n + 1));
//	ai = (uint__t*)malloc(sizeof(uint__t)*nnz);
//	ax = (real__t*)malloc(sizeof(real__t)*nnz);
//	parseCSVInt(fp, ap);
//	fgetc(fp); fgetc(fp); fgetc(fp);
//	parseCSVInt(fp, ai);
//	fgetc(fp); fgetc(fp); fgetc(fp);
//	parseCSVDouble(fp, ax);
//	for (i = 0; i<nnz; ++i)
//		ax[i] = 1.;
//	fclose(fp);
//	ret = NICS_OK;
//	if (ret != NICS_OK) goto EXIT;
//
//	x = (real__t *)malloc(sizeof(real__t)*(n + n));
//	b = x + n;
//	for (i = 0; i<n + n; ++i) x[i] = 1.;
//
//	NicsLU_CreateMatrix(nicslu, n, nnz, ax, ai, ap);
//	nicslu->cfgf[0] = 1.;
//
//	NicsLU_Analyze(nicslu);
//	printf("analysis time: %.8g\n", nicslu->stat[0]);
//
//	ret = NicsLU_CreateScheduler(nicslu);
//	printf("time of creating scheduler: %.8g\n", nicslu->stat[4]);
//	printf("suggestion: %s.\n", ret == 0 ? "parallel" : "sequential");
//
//	NicsLU_CreateThreads(nicslu, atoi(argv[1]), TRUE);
//	printf("total cores: %d, threads created: %d\n", (int)(nicslu->stat[9]), (int)(nicslu->cfgi[5]));
//
//	NicsLU_BindThreads(nicslu, FALSE);
//
//	NicsLU_Factorize_MT(nicslu);
//	printf("factorization time: %.8g\n", nicslu->stat[1]);
//
//	//NicsLU_ReFactorize_MT(nicslu, ax);
//	//printf("re-factorization time: %.8g\n", nicslu->stat[2]);
//
//	NicsLU_Solve(nicslu, x);
//	printf("substitution time: %.8g\n", nicslu->stat[3]);
//
//	NicsLU_Residual(n, ax, ai, ap, x, b, &err, 1, 0);
//	printf("Ax-b (1-norm): %.8g\n", err);
//	NicsLU_Residual(n, ax, ai, ap, x, b, &err, 2, 0);
//	printf("Ax-b (2-norm): %.8g\n", err);
//	NicsLU_Residual(n, ax, ai, ap, x, b, &err, 0, 0);
//	printf("Ax-b (infinite-norm): %.8g\n", err);
//	exit(0);
//
//	printf("NNZ(L+U-I): %ld\n", nicslu->lu_nnz);
//
//	NicsLU_Flops(nicslu, NULL);
//	NicsLU_Throughput(nicslu, NULL);
//	NicsLU_ConditionNumber(nicslu, NULL);
//	printf("flops: %.8g\n", nicslu->stat[5]);
//	printf("throughput (bytes): %.8g\n", nicslu->stat[12]);
//	printf("condition number: %.8g\n", nicslu->stat[6]);
//	NicsLU_MemoryUsage(nicslu, NULL);
//	printf("memory (Mbytes): %.8g\n", nicslu->stat[21] / 1024. / 1024.);
//
//EXIT:
//	NicsLU_Destroy(nicslu);
//	free(ax);
//	free(ai);
//	free(ap);
//	free(nicslu);
//	free(x);
//	return 0;
//}