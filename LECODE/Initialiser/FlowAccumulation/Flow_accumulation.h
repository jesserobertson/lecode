/******************************************************************************

Sequential MFD-md algorithm

Reference:
Qin C et al. An adaptive approach to selecting a flow-partition exponent for a multiple-flow-direction algorithm.
International Journal of Geographical Information Science, 2007, 21(4): 443-458.


Code by ZHAN Li-Jun
Date: 2011.07


Input:  gridded DEM:(Format: ASCII)
Output: 1. getDegreeAccu(): flow accumulation(FA)
		2. getDegreeSCA(): SCA:(SCA = FA * cellSize )
******************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>
#include <cmath>
#include <string.h>
#include <ctype.h>

#define NUL '\0'

using namespace std;


struct HeadFile
{
	int nCols;
	int nRows;
	double xllCenter;
	double yllCenter;
	double cellSize;
	int noData;
};

struct AsciiFile
{
	HeadFile* headfile;
	double* data;
};

static int direct8X[8] = {1, 1, 1, 0, -1, -1, -1, 0 };
static int direct8Y[8] = { -1, 0, 1, 1, 1, 0, -1, -1 };

struct directNode
{
	int dirY;
	int dirX;
};

struct weight
{
	double wei[8];
};

//read the ASCII file
AsciiFile* readAscFile(ifstream& ifstr) ;

//write the ASCII file
void writeAscFile(double* pAccumulation, HeadFile* headfile, const char* filename);

//compute the matrix of flow fraction
weight* getWeightMatrix(double* pDEM, int xSize, int ySize, double cellSize);

// determine the matrix of the flow direction, recorded as a Binary number)
int* getCodeDirectMatrix(double* pDEM,int xSize, int ySize);
int* getReverseDirectMatrix(int* pDirectMatrix, int xSize, int ySize) ;


//(get the indegree matrix)
int* getDegreeMatrix(int* pReverseDirectMatrix, int xSize, int ySize);

// compute the flow accumulation
double* getDegreeFA(int* degreeMatrix, int* reverseDirectMatrix, int* directMatrix, weight* pWeightMatrix, double* pDEM, int xSize, int ySize);


// compute the SCA
double* getDegreeSCA( double* FA, double cellsize, int xSize, int ySize);


char *F2C_trimer( char *str ) ;

void flowacc( char *fdem, char *faccu );

extern "C" void flowacc_( char *fdem, char *faccu ) { flowacc( fdem, faccu );}
