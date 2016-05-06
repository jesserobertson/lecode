/******************************************************************************

Flow-transfer-matrix-based sequential MFD-md algorithm

Reference:
Qin C et al. An adaptive approach to selecting a flow-partition exponent for a multiple-flow-direction algorithm.
International Journal of Geographical Information Science, 2007, 21(4): 443-458.


Code by ZHAN Li-Jun
Date: 2011.07


Input:  gridded DEM:(Format: ASCII)
Output: flow accumulation(FA)

******************************************************************************/
#include "Flow_accumulation.h"

#include <stdio.h>

#ifndef NULL
#define NULL   ((void *) 0)
#endif

//read the ASCII file
AsciiFile* readAscFile(ifstream& ifstr)
{
	AsciiFile* file = new AsciiFile();
	file->headfile = new HeadFile();
	for (int i = 0; i < 6; i++)
	{
		string str;
		switch (i)
		{
		case 0:
			ifstr>>str;
			ifstr>>file->headfile->nCols;
			break;
		case 1:
			ifstr>>str;
			ifstr>>file->headfile->nRows;
			break;
		case 2:
			ifstr>>str;
			ifstr>>file->headfile->xllCenter;
			break;
		case 3:
			ifstr>>str;
			ifstr>>file->headfile->yllCenter;
			break;
		case 4:
			ifstr>>str;
			ifstr>>file->headfile->cellSize;
			break;
		case 5:
			ifstr>>str;
			ifstr>>file->headfile->noData;
			break;
		}
	}

	int ncol = file->headfile->nCols;
	int nrow = file->headfile->nRows;
	file->data = new double[ncol*nrow];

	for (int i = 0; i < nrow; i++) // read the body
	{
		for (int j = 0; j < ncol; j++)
		{
			ifstr>>file->data[i*ncol+j];
		}
	}

	ifstr.close();
	return file;
}

//write the ASCII file
void writeAscFile(double* pAccumulation, HeadFile* headfile, const char* filename)
{
	ofstream fw(filename);
	if (!fw)
	{
		return;
	}
	for (int i = 0; i < 6; i++)
	{
		switch (i)
		{
		case 0:
			fw<<"ncols ";
			fw<<headfile->nCols;
			fw<<endl;
			break;
		case 1:
			fw<<"nrows ";
			fw<<headfile->nRows;
			fw<<endl;
			break;
		case 2:
			fw<<"xllcorner ";
			fw<<headfile->xllCenter;
			fw<<endl;
			break;
		case 3:
			fw<<"yllcorner ";
			fw<<headfile->yllCenter;
			fw<<endl;
			break;
		case 4:
			fw<<"cellsize ";
			fw<<headfile->cellSize;
			fw<<endl;
			break;
		case 5:
			fw<<"NODATA_value ";
			fw<<headfile->noData;
			fw<<endl;
			break;
		}
	}

	for(int i = 0; i < headfile->nRows; i++)
	{
		for (int j = 0; j < headfile->nCols; j++)
		{
			fw<<pAccumulation[i*headfile->nCols + j]<<" ";
		}
		fw<<endl;
	}
	fw.close();
}

//compute the matrix of flow fraction
weight* getWeightMatrix(double* pDEM, int xSize, int ySize, double cellSize)
{
	weight* pWeightMatrix = new weight[xSize*ySize];
	double slopeArray[8];
	double Sum = 0.0;

	double maxslope = 0.0;
	double fe = 0.0;

	for (int i = 0; i < ySize; i++)
	{
		for (int j = 0; j < xSize; j++)
		{
			Sum = 0.0;
			maxslope = 0.0;
			fe = 0.0;

			if(pDEM[i*xSize + j] != -9999)//
			{
				for (int m = 0; m < 8; m++)
				{
					slopeArray[m] = 0.0;
				}
				for (int k = 0; k < 8; k++)
				{
					int Row = i + direct8Y[k];
					int Col = j + direct8X[k];
					double value;
					if (Row >= 0 && Row < ySize && Col >= 0 && Col < xSize)
					{
						value = pDEM[Row*xSize + Col];
						if (value == -9999)//
						{
							slopeArray[k] = 0;
						}
						else
						{
							if (pDEM[i*xSize + j] > value)
							{
								double slope = (pDEM[i*xSize + j] - value)/cellSize;

								if (k % 2 == 0)
								{
									slope = slope / (double)sqrt(2.0) ;
									slopeArray[k] = slope;

								}
								else
								{
									slopeArray[k] = slope;

								}
							}
						}
					}

					if (maxslope < slopeArray[k])
					{
						maxslope = slopeArray[k];
					}
				}

				fe = 8.9*((maxslope>1.0)?1.0:maxslope) + 1.1; // flow partition exponent in (Qin et al., 2007)

				for (int k = 0; k < 8; k++)
				{
					if (slopeArray[k] > 0.0)
					{
						if (k % 2 == 0)
						{
							Sum = Sum + (double)pow(slopeArray[k], fe) * (double)sqrt(2.0) / 4.0;
						}
						else
						{
							Sum = Sum + (double)pow(slopeArray[k], fe) / 2.0;
						}
					}

				}

				double weightuint = 0.0000000;
				double temp = 0.0000000;
				for (int m = 0; m < 8; m++)
				{
					if (slopeArray[m] > 0.0)
					{
						if (m % 2 == 0)
						{
							if (Sum > 0)
							{
								temp = slopeArray[m];
								weightuint = (double)pow(temp, fe)*(double)sqrt(2.0)/(4.0*Sum);
								pWeightMatrix[i*xSize + j].wei[m] = weightuint;
							}
							else
							{
								pWeightMatrix[i*xSize + j].wei[m] = 0;
							}

						}
						else
						{
							if (Sum > 0)
							{
								temp = slopeArray[m];
								weightuint = (double)pow(temp, fe) /(2.0*Sum);
								pWeightMatrix[i*xSize + j].wei[m] = weightuint;
							}
							else
							{
								pWeightMatrix[i*xSize + j].wei[m] = 0;
							}

						}
					}
					else
					{
						pWeightMatrix[i*xSize + j].wei[m] = 0.0;
					}
				}
			}
			else
			{
				for (int m = 0; m < 8; m++)
				{
					pWeightMatrix[i*xSize + j].wei[m] = 0.0;
				}
			}


		}
	}

	return pWeightMatrix;
}

/////////////////////////////////////////////////////////
// determine the matrix of the flow direction, recorded as a Binary number)
int* getCodeDirectMatrix(double* pDEM,int xSize, int ySize)
{
	int* pDirectMatrix = new int[xSize*ySize];
	double slopeArray[8];

	for (int i = 0; i < ySize; i++)
	{
		for (int j = 0; j < xSize; j++)
		{
			if (pDEM[i*xSize + j] != -9999)
			{
				for (int m = 0; m < 8; m++)
				{
					slopeArray[m] = 0.0;
				}
				for (int k = 0; k < 8; k++)
				{
					int Row = i + direct8Y[k];
					int Col = j + direct8X[k];
					double value;
					if (Row >= 0 && Row < ySize && Col >= 0 && Col < xSize)
					{
						value = pDEM[Row*xSize + Col];
						if (value == -9999)
						{
							slopeArray[k] = 0;
						}
						else
						{
							if (pDEM[i*xSize + j] > value)
							{
								double slope = pDEM[i*xSize + j] - value;

								if (k % 2 == 0)
								{
									slope = slope / (double)sqrt(2.0) ;
									slopeArray[k] = slope;
								}
								else
								{
									slopeArray[k] = slope;
								}
							}
						}

					}
				}

				pDirectMatrix[i*xSize + j] = 0;
				double maxSlope = -1.0;
				for (int m = 0; m < 8; m++)
				{
					if (slopeArray[m] > 0)
					{
						if (m == 0)
						{
							pDirectMatrix[i*xSize + j] += 128;
						}
						else if (m == 1)
						{
							pDirectMatrix[i*xSize + j] += 1;
						}
						else if (m == 2)
						{
							pDirectMatrix[i*xSize + j] += 2;
						}
						else if (m == 3)
						{
							pDirectMatrix[i*xSize + j] += 4;
						}
						else if (m == 4)
						{
							pDirectMatrix[i*xSize + j] += 8;
						}
						else if (m == 5)
						{
							pDirectMatrix[i*xSize + j] += 16;
						}
						else if (m == 6)
						{
							pDirectMatrix[i*xSize + j] += 32;
						}
						else if (m == 7)
						{
							pDirectMatrix[i*xSize + j] += 64;
						}
						else
						{
							cout<<"Error!"<<endl;
						}
					}

				}
			}
			else
			{
				pDirectMatrix[i*xSize + j] = 0;
			}


		}
	}

	return pDirectMatrix;
}

//////////////////////////////////
int* getReverseDirectMatrix(int* pDirectMatrix, int xSize, int ySize) //( get the reverse flow direction matrix )
{
	int* reverseDirectMatrix = new int[xSize*ySize];

	for (int i = 0; i < ySize; i++)
	{
		for (int j = 0; j < xSize; j++)
		{
			reverseDirectMatrix[i*xSize + j] = 0;
		}
	}
	for (int i = 0; i < ySize; i++)
	{
		for (int j = 0; j < xSize; j++)
		{
			int mask = pDirectMatrix[i*xSize + j];
			if (mask & 1)
			{
				reverseDirectMatrix[i*xSize + j+1] += 16;
			}
			if (mask & 2)
			{
				reverseDirectMatrix[(i+1)*xSize + j+1] += 32;
			}
			if (mask & 4)
			{
				reverseDirectMatrix[(i+1)*xSize + j] += 64;
			}
			if (mask & 8)
			{
				reverseDirectMatrix[(i+1)*xSize + j-1] += 128;
			}
			if (mask & 16)
			{
				reverseDirectMatrix[i*xSize + j-1] += 1;
			}
			if (mask & 32)
			{
				reverseDirectMatrix[(i-1)*xSize + j-1] += 2;
			}
			if (mask & 64)
			{
				reverseDirectMatrix[(i-1)*xSize + j] += 4;
			}
			if (mask & 128)
			{
				reverseDirectMatrix[(i-1)*xSize + j+1] += 8;
			}

		}
	}

	return reverseDirectMatrix;
}

////////////////////////////////////
//(get the indegree matrix)
int* getDegreeMatrix(int* pReverseDirectMatrix, int xSize, int ySize)
{
	int* reverseDirectMatrix = new int[xSize*ySize];

	for (int i = 0; i < ySize; i++)
	{
		for (int j = 0; j < xSize; j++)
		{
			reverseDirectMatrix[i*xSize + j] = 0;
			int mask = pReverseDirectMatrix[i*xSize + j];
			if (mask & 1)
			{
				reverseDirectMatrix[i*xSize + j]++;
			}
			if (mask & 2)
			{
				reverseDirectMatrix[i*xSize + j]++;
			}
			if (mask & 4)
			{
				reverseDirectMatrix[i*xSize + j]++;
			}
			if (mask & 8)
			{
				reverseDirectMatrix[i*xSize + j]++;
			}
			if (mask & 16)
			{
				reverseDirectMatrix[i*xSize + j]++;
			}
			if (mask & 32)
			{
				reverseDirectMatrix[i*xSize + j]++;
			}
			if (mask & 64)
			{
				reverseDirectMatrix[i*xSize + j]++;
			}
			if (mask & 128)
			{
				reverseDirectMatrix[i*xSize + j]++;
			}

		}
	}

	return reverseDirectMatrix;
}

////////////////////////////////////
// compute the flow accumulation
//
double* getDegreeFA(int* degreeMatrix, int* reverseDirectMatrix, int* directMatrix, weight* pWeightMatrix, double* pDEM, int xSize, int ySize)
{
	double* accumulation = new double[xSize*ySize];
	for (int i = 0; i < ySize; i++)
	{
		for (int j = 0; j < xSize; j++)
		{
			accumulation[i*xSize + j] = 1.0;
		}
	}

	bool flag = true;
	int count = 0;


	while(flag)
	{
		flag = false;
		//count++;
		for (int i = 0; i < ySize; i++)
		{
			for (int j = 0; j < xSize; j++)
			{
				if (pDEM[i*xSize + j] == -9999)
				{
					accumulation[i*xSize + j] = -9999;
				}
				else
				{
					if (degreeMatrix[i*xSize + j] == 0)
					{
						degreeMatrix[i*xSize + j] = -1;
						double accu = 0.0;
						int dir = reverseDirectMatrix[i*xSize + j];
						if (dir & 1)
						{
							accu += accumulation[i*xSize + j+1]*pWeightMatrix[i*xSize + j+1].wei[5];
						}
						if (dir & 2)
						{
							accu += accumulation[(i+1)*xSize + j+1]*pWeightMatrix[(i+1)*xSize + j+1].wei[6];
						}
						if (dir & 4)
						{
							accu += accumulation[(i+1)*xSize + j]*pWeightMatrix[(i+1)*xSize + j].wei[7];
						}
						if (dir & 8)
						{
							accu += accumulation[(i+1)*xSize + j-1]*pWeightMatrix[(i+1)*xSize + j-1].wei[0];
						}
						if (dir & 16)
						{
							accu += accumulation[i*xSize + j-1]*pWeightMatrix[i*xSize + j-1].wei[1];
						}
						if (dir & 32)
						{
							accu += accumulation[(i-1)*xSize + j-1]*pWeightMatrix[(i-1)*xSize + j-1].wei[2];
						}
						if (dir & 64)
						{
							accu += accumulation[(i-1)*xSize + j]*pWeightMatrix[(i-1)*xSize + j].wei[3];
						}
						if (dir & 128)
						{
							accu += accumulation[(i-1)*xSize + j+1]*pWeightMatrix[(i-1)*xSize + j+1].wei[4];
						}
						accumulation[i*xSize +j] += accu;
						/////////////////////////////////////
						/////////////////////////////////////
						/////////////////////////////////////
						int dir1 = directMatrix[i*xSize + j];
						if (dir1 & 1)
						{
							degreeMatrix[i*xSize + j+1]--;
						}
						if (dir1 & 2)
						{
							degreeMatrix[(i+1)*xSize + j+1]--;
						}
						if (dir1 & 4)
						{
							degreeMatrix[(i+1)*xSize + j]--;
						}
						if (dir1 & 8)
						{
							degreeMatrix[(i+1)*xSize + j-1]--;
						}
						if (dir1 & 16)
						{
							degreeMatrix[i*xSize + j-1]--;
						}
						if (dir1 & 32)
						{
							degreeMatrix[(i-1)*xSize + j-1]--;
						}
						if (dir1 & 64)
						{
							degreeMatrix[(i-1)*xSize + j]--;
						}
						if (dir1 & 128)
						{
							degreeMatrix[(i-1)*xSize + j+1]--;
						}
						//////////////
						//////////////
						flag = true;
					}
				}

			}
		}
	}


	//cout<<count<<endl;

	return accumulation;
}

////////////////////////////////////
// compute the SCA
//
double* getDegreeSCA( double* FA, double cellsize, int xSize, int ySize)
{
	double* SCA = new double[xSize*ySize];

	for (int i = 0; i < ySize; i++)
	{
		for (int j = 0; j < xSize; j++)
		{
			if (FA[i*xSize + j] > 0)
			{
				SCA[i*xSize + j] = FA[i*xSize + j]*cellsize;
			}
			else
			{
				SCA[i*xSize + j] = -9999;
			}
		}
	}

	return SCA;
}

/* ============================================================================
 * Function F2C_trimer()
 * When passing string from fortran to C some end character are mismatching.
 * This function ensures the consistency between string for the 2 languages
 * ============================================================================*/
char *F2C_trimer( char *str )
{
	char *ibuf = str, *obuf = str;
	int i = 0, cnt = 0;
	if (str) {
		for (ibuf = str; *ibuf && isspace(*ibuf); ++ibuf) ;
		if (str != ibuf) memmove(str, ibuf, ibuf - str);
		while (*ibuf) {
			if (isspace(*ibuf) && cnt)
				ibuf++;
			else {
				if (!isspace(*ibuf))
					cnt = 0;
				else {
					*ibuf = ' ';
					cnt = 1;
				}
				obuf[i++] = *ibuf++;
			}
		}
		obuf[i] = NUL;
		while (--i >= 0) {
			if (!isspace(obuf[i])) break;
		}
		obuf[++i] = NUL;
	}

	return str;
}

void flowacc( char *fdem, char *faccu )
{

	const char* demf ;
	const char* accum ;

	demf = F2C_trimer( fdem ) ;
	accum = F2C_trimer( faccu ) ;

	ifstream ifstr;
	ifstr.open(demf);

	 //read the file
	AsciiFile* ascfile = readAscFile(ifstr);

	double* buf = ascfile->data;
	int xSize = ascfile->headfile->nCols;
	int ySize = ascfile->headfile->nRows;
	double cellSize = ascfile->headfile->cellSize;

	// pre-processing
	clock_t starttime1, endtime1, starttime2, endtime2;

	starttime1 = clock();
	int* codeDirectMatrix = getCodeDirectMatrix(buf, xSize, ySize);
	int* reverseDirectMatrix = getReverseDirectMatrix(codeDirectMatrix, xSize, ySize);
	int* degreeMatrix= getDegreeMatrix(reverseDirectMatrix,xSize,ySize);
	weight* pWeightMatrix = getWeightMatrix(buf,xSize,ySize,(double)cellSize);
	endtime1 =clock();

	//1. calculate flow accumulation
	starttime2 = clock();
	double* FA = getDegreeFA(degreeMatrix,reverseDirectMatrix,codeDirectMatrix,pWeightMatrix,buf,xSize,ySize);
	endtime2 = clock();

	//output the result
	writeAscFile(FA,ascfile->headfile,accum);

	//clear the resources
	delete []buf;
	delete []codeDirectMatrix;
	delete []reverseDirectMatrix;
	delete []degreeMatrix;
	delete []pWeightMatrix;
	delete []FA;

	delete ascfile->headfile;
	delete ascfile;

}


