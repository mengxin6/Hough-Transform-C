###input parameters explaination for function: 
#####void HTLineDetection(unsigned char* BinaryImageDst, int* LineNumber, struct LineParameter* DetectedLine, int nHeight, int nWidth)
#####BinaryImageDst: 
an unsigned char pointer pointing to the location where a binary image is located. The binary image is used
as the input of hough transform, and can be get from image processing techniques like edge detection.
#####LineNumber: 
an int pointer pointing to the location where the number of detected lines are to be stored when running this function.
#####DetectedLine: 
a struct pointer pointing to the location where parameters of detected lines are to be stored when running this function.
#####nHeight: 
an int indicating the height of the input binary image (number of pixels in a column)
#####nWidth: 
an int indicating the width of the input binary image (number of pixels in a row).
