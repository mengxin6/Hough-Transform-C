	struct LineParameter{
		float angle;
		float distance;
	};

void HTLineDetection(unsigned char* BinaryImageDst, int* LineNumber, struct LineParameter* DetectedLine, int nHeight, int nWidth);


void HoughTransform_Line(unsigned char* BinaryImageDst, float MinAngle, float AngleInterval, int No_Angles, float MinDistance,  float DistanceInterval, int No_Distances, unsigned int** VoteTable, int nHeight, int nWidth );

void FindMaxVote(unsigned int** VoteTable, int No_Angles, int No_Distances, int* M, int* N);   // To get the *M and *N,  VoteTable[*M][*N] has the maximum vote number

int FindArrayMax(unsigned int* Array, int No_of_Element); //return int n that Array[n] has the maximum value is this array

int Larger_in_Array(unsigned int* Array, int x, int y);   // Find the larger one between Array[x] and Array[y], and return the corresponding x/y.
