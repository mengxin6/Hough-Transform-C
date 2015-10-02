#include"HTLineDetect.h"
#include<math.h>
#include<stdlib.h>
#include<iostream>
using namespace std;


const double pi = 3.1415926;

void HTLineDetection(unsigned char* BinaryImageDst, int* LineNumber, struct LineParameter* DetectedLine, int nHeight, int nWidth)
{
	int VoteThreshold=100;

	float DiagonalLength = sqrt((float)(nHeight*nHeight+nWidth*nWidth));
	float MinimumAngle = 0;
	float MaximumAngle = 90;
	float AngleInterval_Coarse = 2;
	float AngleInterval_Fine = 0.1;
	float DistanceInterval_Coarse =2;
	float DistanceInterval_Fine =1;

	int NumberofAngles_Coarse = (int)(MaximumAngle/AngleInterval_Coarse)+1;
	int NumberofDistances_Coarse = (int)DiagonalLength/DistanceInterval_Coarse+1;

	//build and initialize the vote table for coarse angle and distance
	unsigned int** VoteTable_Coarse = (unsigned int**)malloc( NumberofAngles_Coarse*sizeof(unsigned int*));
	for(int i=0; i< NumberofAngles_Coarse;i++)
	{
		VoteTable_Coarse[i] = (unsigned int*)malloc(NumberofDistances_Coarse*sizeof(unsigned int));
	}
	for(int i=0; i<NumberofAngles_Coarse;i++)
	{
		for(int j=0; j<NumberofDistances_Coarse;j++)
		{
			VoteTable_Coarse[i][j]=0;
		}
	}
	//vote table initialization finished

	
	HoughTransform_Line(BinaryImageDst, MinimumAngle, AngleInterval_Coarse, NumberofAngles_Coarse, 0, DistanceInterval_Coarse, NumberofDistances_Coarse, VoteTable_Coarse, nHeight, nWidth); 
	
	int MaxVote=110;
	int M=0;
	int N=0; //  VoteTable_Coarse[M][N] has the maximum votes. 

	while (MaxVote>VoteThreshold)
	{
		FindMaxVote(VoteTable_Coarse,NumberofAngles_Coarse, NumberofDistances_Coarse, &M, &N);  
		MaxVote = VoteTable_Coarse[M][N];
		float Angle_withMaxVote_Coarse = MinimumAngle+M*AngleInterval_Coarse;
		float Distance_withMaxVote_Coarse = N*DistanceInterval_Coarse;
		(*LineNumber)++;
		DetectedLine[(*LineNumber)-1].angle = Angle_withMaxVote_Coarse;
		DetectedLine[(*LineNumber)-1].distance=Distance_withMaxVote_Coarse;
		cout<<Angle_withMaxVote_Coarse<<" "<<Distance_withMaxVote_Coarse<<" "<<MaxVote<<endl;

		VoteTable_Coarse[M][N]=0;

		if(N>0){
			VoteTable_Coarse[M][N-1]=0;
			if(M>0){
				VoteTable_Coarse[M-1][N-1]=0;
			}
			if(M<NumberofAngles_Coarse-1){
				VoteTable_Coarse[M+1][N-1]=0;
			}
		}

		if(N<NumberofDistances_Coarse-1){
			VoteTable_Coarse[M][N+1]=0;
			if(M>0){
				VoteTable_Coarse[M-1][N+1]=0;
			}
			if(M<NumberofAngles_Coarse-1){
				VoteTable_Coarse[M+1][N+1]=0;
			}
		}

		if(M>0)
		{
			VoteTable_Coarse[M-1][N]=0;
		}

		if(M<NumberofAngles_Coarse-1)
		{
			VoteTable_Coarse[M+1][N]=0;
		}

		
	}
	
			
	


	for(int i=0; i< NumberofAngles_Coarse;i++)
	{
		free(VoteTable_Coarse[i]);
	}

	free(VoteTable_Coarse);


}






void HoughTransform_Line(unsigned char* BinaryImageDst, float MinAngle, float AngleInterval, int No_Angles, float MinDistance,  float DistanceInterval, int No_Distances,unsigned int** VoteTable, int nHeight, int nWidth )
{
	for (int x=0; x<nWidth; x++)
	{
		for(int y=0; y<nHeight; y++)
		{
			if(*(BinaryImageDst+y*nWidth+x)==255)
			{
				for (int i=0; i<No_Angles; i++)
				{
					float Angle = MinAngle+AngleInterval*((float)i);
					float Distance = ((float)x)*cos(Angle/180*pi)+((float)y)*sin(Angle/180*pi);
					int j = (int)((Distance-MinDistance)/DistanceInterval);
					if(j < No_Distances)
						{
							VoteTable[i][j]=((VoteTable[i][j])+1);
						}

				}
			}
		}
	}
}





void FindMaxVote(unsigned int** VoteTable, int No_Angles, int No_Distances, int* M, int* N)
{
	unsigned int *MaxVote_for_Diff_Angles = (unsigned int*)malloc(No_Angles*sizeof(unsigned int));

	for(int i=0; i<No_Angles; i++)
	{
		int j = FindArrayMax(VoteTable[i], No_Distances); 
		MaxVote_for_Diff_Angles[i] = VoteTable[i][j];
	}

	*M = FindArrayMax(MaxVote_for_Diff_Angles, No_Angles);
	*N = FindArrayMax(VoteTable[*M], No_Distances);

	free(MaxVote_for_Diff_Angles);
}


int FindArrayMax(unsigned int* Array, int No_of_Element)
{
	int N=0; //N is used to indicate the order of the element int the array which has the max value: Array[N] is the maximum in Array
	for(int i=1; i<No_of_Element; i++)
	{
		N = Larger_in_Array(Array, N, i);
	}

	return N;

}




int Larger_in_Array(unsigned int* Array, int x, int y)
{
	if(Array[x]>=Array[y])
		return x;
	else 
		return y;
}
