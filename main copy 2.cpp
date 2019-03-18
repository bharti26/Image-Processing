/****************************************************************************/
/*  Project Description: Demo code for Assignment 4				      		*/
/*  Jing Zhang  Nov. 2, 2007                         						*/
/****************************************************************************/


#include <stdio.h>
#include <fftw3.h>
#include<math.h>
#include "image.h"
#include "ipTool.h"
#include<string>
#define WIDTH  256
#define HEIGHT 256
#include <iostream>
#define MAXLEN 256
using namespace std;
//global varialbles;
//function declaration;

#include<stdlib.h>


void RGBToHSI(int R, int G, int B, double &H, double &S, double &I)
{

    I = (double)(R + G + B) / (double)3;
    double m = min(R, min(G, B));
    double pi = 3.141;

    if (I == 0)
    {
        S = 0;
    }
    else if (I > 0)
    {
        S = 1.0 - (m / I);
    }

    if (G >= B)
    {
        if (R == G && G == B)
        {
            H = 0;
        }
        else
        {
            double temp = acos((R - (G / 2.0) - (B / 2.0)) / (sqrt(pow(R, 2) + pow(G, 2) + pow(B, 2) - (R*G) - (R*B) - (G*B))));
            H = (temp * (180.0)) / pi;
        }
    }
    else if (B > G)
    {
        double temp = acos((R - (G / 2.0) - (B / 2.0)) / (sqrt(pow(R, 2) + pow(G, 2) + pow(B, 2) - (R*G) - (R*B) - (G*B))));

        H = 360 - ((temp * (180.0)) / pi);
    }
}

void HSIToRGB(double H, double S, double I, int &R, int &G, int &B)
{
    double pi = 3.141;

    if (H == 0)
    {
        R = I + (2 * I*S);
        G = I - (I*S);
        B = I - (I*S);
    }
    else if (H > 0 && H < 120)
    {
        R = I + I*S*(cos(H*(pi / 180.0)) / cos((pi / 3.0) - H*(pi / 180.0)));
        G = I + I*S*(1.0 - (cos(H*(pi / 180.0)) / cos((pi / 3.0) - H*(pi / 180.0))));
        B = I - I*S;
    }
    else if (H == 120)
    {
        R = I - (I*S);
        G = I + (2 * I*S);
        B = I - (I*S);
    }
    else if (H > 120 && H < 240)
    {
        R = I - I*S;
        G = I + I*S*((cos(H*(pi / 180.0) - (2.0 * pi) / 3.0)) / (cos((pi)-H*(pi / 180.0))));
        B = I + I*S*(1 - ((cos(H*(pi / 180.0) - (2.0 * pi) / 3.0)) / (cos((pi)-H*(pi / 180.0)))));
    }
    else if (H == 240)
    {
        R = I - (I*S);
        G = I - (I*S);
        B = I + (2 * I*S);
    }
    else if (H > 240 && H < 360)
    {
        R = I + (I*S)*(1 - cos(H*(pi / 180.0) - (4.0*pi / 3.0)) / cos((5.0*pi / 3.0) - H*(pi / 180.0)));
        G = I - (I*S);
        B = I + (I*S)*(cos(H*(pi / 180.0) - (4.0*pi / 3.0)) / cos((5.0*pi / 3.0) - H*(pi / 180.0)));
    }
}

void fft(int rows, int columns,Image &src, Image &magDFT,Image &centered,Image &freqDFT,Image &freqCentered,Image &tgt)
{
    //cout<<"rahul";
     fftw_plan planR, planG, planB,filterPlan;
    fftw_complex *inR, *inG, *inB, *outR,*outRfinal;

    // allocate input arrays
    inR = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * rows * columns);

    for (int i=0;i<rows;i++) {
        for (int j=0;j<columns;j++) {
        double value =src(i,j);
        inR[j+columns * i][0]=value;
        }
    }


    // allocate output arrays
    outR = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * rows * columns);
        outRfinal = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * rows * columns);


    // create plans
    planR = fftw_plan_dft_2d(rows,columns, inR, outR, FFTW_FORWARD, FFTW_ESTIMATE);


    // TODO: assign color-values to input arrays

    // perform FORWARD fft
    fftw_execute(planR);

    for (int i=0;i<rows;i++) {
        for (int j=0;j<columns;j++) {

        double realPart = outR[j+columns *i][0] / (double) (rows*columns);
        double imagPart = outR[j+columns *i][1] / (double) (rows*columns);
        double mag = sqrt((realPart*realPart) + (imagPart *imagPart));

        if(mag>1)
        {
            mag=255;

        }

        else{

            mag=mag*255;
        }
        magDFT(i,j)=(int) mag;
        }
    }
for (int i=rows/2;i<rows;i++) {
        for (int j=columns/2;j<columns;j++) {

// for Freq

				int v1 = freqDFT(i, j);

				int v2 = freqDFT(i - ((rows) / 2), j - ((columns) / 2));

				freqCentered(i, j) = v2;

				freqCentered((i - (rows / 2)), j - (columns / 2)) = v1;


v1 = magDFT(i,j);
 v2= magDFT((i-(rows/2)),(j-(columns/2)));
centered(i,j)=v2;
centered((i-rows/2),j-(columns/2)) =v1;


        }
}

for (int i=rows/2;i<rows;i++) {
        for (int j=0;j<columns/2;j++) {

int v1 = freqDFT(i, j);

				int v2 = freqDFT(i - ((rows) / 2), j - ((columns) / 2));

				freqCentered(i, j) = v2;

				freqCentered((i - (rows / 2)), j - (columns / 2)) = v1;



 v1 = magDFT(i,j);
 v2= magDFT((i-(rows/2)),(j+(columns/2)));
centered(i,j)=v2;
centered((i-rows/2),j+(columns/2)) =v1;


        }
}

/*
filterPlan = fftw_plan_dft_2d(rows, columns, outR, outRfinal, FFTW_BACKWARD, FFTW_ESTIMATE);



		fftw_execute(filterPlan);



		for (int i = 0; i < rows; i++){

			for (int j = 0; j < columns; j++)

			{

				int k = i*(columns) + j;

				// normalize values
				double realPart = outRfinal[k][0] / (double)(rows * columns);
				double imagPart = outRfinal[k][1] / (double)(rows * columns);

				// magnitude
				double magnitude = sqrt((realPart * realPart) + (imagPart * imagPart));


				invDFT(i, j) = (int)magnitude;

			}

		}
		*/


    fftw_destroy_plan(planR);

    fftw_free(inR); fftw_free(outR);
}

void filtering(int rows, int columns,Image &src, Image &magDFT,Image &centered,Image &freqDFT,Image &freqCentered,Image &filtered,Image &tgt,int filter_type, char *filteredName)
{
    //cout<<"paul";
   // int filter_type;
    int thresh1,thresh2;
    thresh1=10;
    thresh2= 25;
    int startx =100;
    int starty=300;
    int sizex= 40;
    int sizey= 140;

  //char filteredName[250];
   // strcpy(filteredName,argv[7]);
    filtered.initialize(src.NR,src.NC);

    rows = src.NR;
    columns= src.NC;

    fftw_plan planR, planG, planB,filterPlan;
    fftw_complex *inR, *inG, *inB, *outR,*outRfinal;

    // allocate input arrays
    inR = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * rows * columns);

    for (int i=0;i<rows;i++) {
        for (int j=0;j<columns;j++) {
                if(i>= startx && i<= startx+sizex && j>= starty && j<=starty+sizey){
        double value =src(i,j);
        inR[j+columns * i][0]=value;
        }
    }
    }

    // allocate output arrays
    outR = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * rows * columns);
        outRfinal = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * rows * columns);


    // create plans
    planR = fftw_plan_dft_2d(rows,columns, inR, outR, FFTW_FORWARD, FFTW_ESTIMATE);


    // TODO: assign color-values to input arrays

    // perform FORWARD fft
    fftw_execute(planR);

    for (int i=0;i<rows;i++) {
        for (int j=0;j<columns;j++) {
                if(i>= startx && i<= startx+sizex && j>= starty && j<=starty+sizey){

                freqDFT(i,j) = ( outR[j+columns *i][0] )/ 255;

        double realPart = outR[j+columns *i][0] / (double) (rows*columns);
        double imagPart = outR[j+columns *i][1] / (double) (rows*columns);
        double mag = sqrt((realPart*realPart) + (imagPart *imagPart));

        if(mag>1)
        {
            mag=255;

        }

        else{

            mag=mag*255;
        }
        magDFT(i,j)=(int) mag;
        }
    }
    }
for (int i=rows/2;i<rows;i++) {
        for (int j=columns/2;j<columns;j++) {
                if(i>= startx && i<= startx+sizex && j>= starty && j<=starty+sizey){

// for Freq

				int v1 = freqDFT(i, j);

				int v2 = freqDFT(i - ((rows) / 2), j - ((columns) / 2));

				freqCentered(i, j) = v2;

				freqCentered((i - (rows / 2)), j - (columns / 2)) = v1;


v1 = magDFT(i,j);
 v2= magDFT((i-(rows/2)),(j-(columns/2)));
centered(i,j)=v2;
centered((i-rows/2),j-(columns/2)) =v1;


        }
}

}
for (int i=rows/2;i<rows;i++) {
        for (int j=0;j<columns/2;j++) {
                if(i>= startx && i<= startx+sizex && j>= starty && j<=starty+sizey){

int v1 = freqDFT(i, j);

				int v2 = freqDFT(i - ((rows) / 2), j - ((columns) / 2));

				freqCentered(i, j) = v2;

				freqCentered((i - (rows / 2)), j - (columns / 2)) = v1;



 v1 = magDFT(i,j);
 v2= magDFT((i-(rows/2)),(j+(columns/2)));
centered(i,j)=v2;
centered((i-rows/2),j+(columns/2)) =v1;


        }
}


}



/*
		for (int i = 0; i < rows; i++){

			for (int j = 0; j < columns; j++)

			{

				int k = i*(columns) + j;

				// normalize values
				double realPart = outRfinal[k][0] / (double)(rows * columns);
				double imagPart = outRfinal[k][1] / (double)(rows * columns);

				// magnitude
				double magnitude = sqrt((realPart * realPart) + (imagPart * imagPart));


				invDFT(i, j) = (int)magnitude;

			}

		}
*/
if(filter_type==2){

        cout<< "lowpass";
    for (int i = 0; i < rows; i++){

				for (int j = 0; j < columns; j++)

				{
				    if(i>= startx && i<= startx+sizex && j>= starty && j<=starty+sizey){

					double D = sqrt(pow(rows / 2 - i, 2) + pow(columns / 2 - j, 2));//do + x and y with i and j for ROI

					if ((D) <= thresh1)

					{

						filtered(i, j) = freqCentered(i,j);

					}

					else

					{

						filtered(i, j) = 0;

					}

				}

			}

    }
        filtered.save(filteredName);

        //Do the folding to bring to the center - First and Third Quadrants

			for (int i = rows / 2; i< rows; i++){

				for (int j = columns / 2; j < columns; j++)

				{
				    if(i>= startx && i<= startx+sizex && j>= starty && j<=starty+sizey){

					// for Freq

					int v1 = filtered(i, j);//Start Get Value from center Middle till end column

					int v2 = filtered(i - ((rows) / 2), j - ((columns) / 2));//start get value from  (0,0) till (0,column)

					filtered(i, j) = v2;

					filtered((i - (rows / 2)), j - (columns / 2)) = v1;

				}

			}
			}


			//Second and Fourth Quadrants

			for (int i = (rows) / 2; i < rows; i++){

				for (int j = 0; j < (columns) / 2; j++)

				{
				    if(i>= startx && i<= startx+sizex && j>= starty && j<=starty+sizey){

					// for frquency image.

					int v1 = filtered(i, j);

					int v2 = filtered(i - ((rows) / 2), j + ((columns) / 2));

					filtered(i, j) = v2;

					filtered((i - (rows / 2)), j - (columns / 2)) = v1;

				}

			}
			}




			//Multiply the filter values with the FFT OUT

			//outDFT



			for (int i = 0; i < rows; i++){

				for (int j = 0; j < columns; j++)

				{
				    if(i>= startx && i<= startx+sizex && j>= starty && j<=starty+sizey){

					if (filtered(i, j) == 0)

					{

						outR[i*columns + j][0] = 0;

						outR[i*columns + j][1] = 0;

					}

				}

			}
			}
						filter_type=1;

}


if(filter_type==3){


        cout<< "highpass";
    for (int i = 0; i < rows; i++){

				for (int j = 0; j < columns; j++)

				{

					double D = sqrt(pow(rows / 2 - i, 2) + pow(columns / 2 - j, 2));//do + x and y with i and j for ROI

					if ((D) >= thresh1)

					{

						filtered(i, j) = freqCentered(i,j);

					}

					else

					{

						filtered(i, j) = 0;

					}

				}

			}


        filtered.save(filteredName);


			for (int i = rows / 2; i< rows; i++){

				for (int j = columns / 2; j < columns; j++)

				{

					int v1 = filtered(i, j);

					int v2 = filtered(i - ((rows) / 2), j - ((columns) / 2));

					filtered(i, j) = v2;

					filtered((i - (rows / 2)), j - (columns / 2)) = v1;

				}

			}


			for (int i = (rows) / 2; i < rows; i++){

				for (int j = 0; j < (columns) / 2; j++)

				{


					int v1 = filtered(i, j);

					int v2 = filtered(i - ((rows) / 2), j + ((columns) / 2));

					filtered(i, j) = v2;

					filtered((i - (rows / 2)), j - (columns / 2)) = v1;

				}

			}

			for (int i = 0; i < rows; i++){

				for (int j = 0; j < columns; j++)

				{

					if (filtered(i, j) == 0)

					{

						outR[i*columns + j][0] = 0;

						outR[i*columns + j][1] = 0;

					}

				}

			}
						filter_type=1;

}



if(filter_type==4)
{
cout << "\nBandpass\n";

			for (int i = 0; i < rows; i++){

				for (int j = 0; j < columns; j++)

				{

					double D = sqrt(pow(rows / 2 - i, 2) + pow(columns / 2 - j, 2));

					if ((D) >= thresh1 && (D) <= thresh2)
					{

						filtered(i, j) = freqCentered(i, j);

					}

					else

					{

						filtered(i, j) = 0;

					}

				}
			}
filtered.save(filteredName);




			for (int i = rows / 2; i< rows; i++){

				for (int j = columns / 2; j < columns; j++)

				{

					int v1 = filtered(i, j);

					int v2 = filtered(i - ((rows) / 2), j - ((columns) / 2));

					filtered(i, j) = v2;

					filtered((i - (rows / 2)), j - (columns / 2)) = v1;

				}

			}


			for (int i = (rows) / 2; i < rows; i++){

				for (int j = 0; j < (columns) / 2; j++)

				{


					int v1 = filtered(i, j);

					int v2 = filtered(i - ((rows) / 2), j + ((columns) / 2));

					filtered(i, j) = v2;

					filtered((i - (rows / 2)), j - (columns / 2)) = v1;

				}

			}

			for (int i = 0; i < rows; i++){

				for (int j = 0; j < columns; j++)

				{

					if (filtered(i, j) == 0)

					{

						outR[i*columns + j][0] = 0;

						outR[i*columns + j][1] = 0;

					}

				}

			}
						filter_type=1;

}


if(filter_type==5)
{
cout << "\nBandstop\n";

			for (int i = 0; i < rows; i++){

				for (int j = 0; j < columns; j++)

				{

					double D = sqrt(pow(rows / 2 - i, 2) + pow(columns / 2 - j, 2));

					if ((D) <= thresh1 || (D) >= thresh2)

					{

						filtered(i, j) = freqCentered(i, j);

					}

					else

					{

						filtered(i, j) = 0;

					}

				}
			}

filtered.save(filteredName);




			for (int i = rows / 2; i< rows; i++){

				for (int j = columns / 2; j < columns; j++)

				{

					int v1 = filtered(i, j);

					int v2 = filtered(i - ((rows) / 2), j - ((columns) / 2));

					filtered(i, j) = v2;

					filtered((i - (rows / 2)), j - (columns / 2)) = v1;

				}

			}


			for (int i = (rows) / 2; i < rows; i++){

				for (int j = 0; j < (columns) / 2; j++)

				{


					int v1 = filtered(i, j);

					int v2 = filtered(i - ((rows) / 2), j + ((columns) / 2));

					filtered(i, j) = v2;

					filtered((i - (rows / 2)), j - (columns / 2)) = v1;

				}

			}

			for (int i = 0; i < rows; i++){

				for (int j = 0; j < columns; j++)

				{

					if (filtered(i, j) == 0)

					{

						outR[i*columns + j][0] = 0;

						outR[i*columns + j][1] = 0;

					}

				}

			}
			filter_type=1;
}



filterPlan = fftw_plan_dft_2d(rows, columns, outR, outRfinal, FFTW_BACKWARD, FFTW_ESTIMATE);

		fftw_execute(filterPlan);

if(filter_type==1)
{

    for (int i = 0; i < rows; i++){

			for (int j = 0; j < columns; j++)

			{

				int k = i*(columns) + j;

				// normalize values
				double realPart = outRfinal[k][0] / (double)(rows * columns);
				double imagPart = outRfinal[k][1] / (double)(rows * columns);

				// magnitude
				double magnitude = sqrt((realPart * realPart) + (imagPart * imagPart));


				tgt(i, j) = (int)magnitude;

			}

		}
}


 fftw_destroy_plan(planR);

    fftw_free(inR); fftw_free(outR);
		}




void rgbimage(int rows, int columns,Image &src, Image &magDFT,Image &centered,Image &freqDFT,Image &freqCentered,Image &filtered,Image &tgt,int filter_type,char *filteredName)
{


    int Red,Green,Blue;
    double Hue,Saturation,Intensity;

int thresh1,thresh2;
    thresh1=10;
    thresh2= 25;
  //char filteredName[250];
   // strcpy(filteredName,argv[7]);
    filtered.initialize(src.NR,src.NC);

    rows = src.NR;
    columns= src.NC;

    fftw_plan planR, planG, planB,filterPlan;
    fftw_complex *inR, *inG, *inB, *outR,*outRfinal;

    // allocate input arrays
    inR = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * rows * columns);

    for (int i=0;i<rows;i++) {
        for (int j=0;j<columns;j++) {
              Red = src(i,j,0);
            Green = src(i,j,1);
            Blue = src(i,j,2);

            RGBToHSI(Red,Green,Blue,Hue,Saturation,Intensity);
        inR[j+columns * i][0]=Intensity;
        }
    }


    // allocate output arrays
    outR = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * rows * columns);
        outRfinal = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * rows * columns);


    // create plans
    planR = fftw_plan_dft_2d(rows,columns, inR, outR, FFTW_FORWARD, FFTW_ESTIMATE);


    // TODO: assign color-values to input arrays

    // perform FORWARD fft
    fftw_execute(planR);

    for (int i=0;i<rows;i++) {
        for (int j=0;j<columns;j++) {

                freqDFT(i,j) = ( outR[j+columns *i][0] )/ 255;


        double realPart = outR[j+columns *i][0] / (double) (rows*columns);
        double imagPart = outR[j+columns *i][1] / (double) (rows*columns);
        double mag = sqrt((realPart*realPart) + (imagPart *imagPart));

        if(mag>1)
        {
            mag=255;

        }

        else{

            mag=mag*255;
        }
        magDFT(i,j)=(int) mag;
        }
    }
for (int i=rows/2;i<rows;i++) {
        for (int j=columns/2;j<columns;j++) {

// for Freq

				int v1 = freqDFT(i, j);

				int v2 = freqDFT(i - ((rows) / 2), j - ((columns) / 2));

				freqCentered(i, j) = v2;

				freqCentered((i - (rows / 2)), j - (columns / 2)) = v1;


v1 = magDFT(i,j);
 v2= magDFT((i-(rows/2)),(j-(columns/2)));
centered(i,j)=v2;
centered((i-rows/2),j-(columns/2)) =v1;


        }
}

for (int i=rows/2;i<rows;i++) {
        for (int j=0;j<columns/2;j++) {

int v1 = freqDFT(i, j);

				int v2 = freqDFT(i - ((rows) / 2), j - ((columns) / 2));

				freqCentered(i, j) = v2;

				freqCentered((i - (rows / 2)), j - (columns / 2)) = v1;



 v1 = magDFT(i,j);
 v2= magDFT((i-(rows/2)),(j+(columns/2)));
centered(i,j)=v2;
centered((i-rows/2),j+(columns/2)) =v1;


        }
}






/*
		for (int i = 0; i < rows; i++){

			for (int j = 0; j < columns; j++)

			{

				int k = i*(columns) + j;

				// normalize values
				double realPart = outRfinal[k][0] / (double)(rows * columns);
				double imagPart = outRfinal[k][1] / (double)(rows * columns);

				// magnitude
				double magnitude = sqrt((realPart * realPart) + (imagPart * imagPart));


				invDFT(i, j) = (int)magnitude;

			}

		}
*/
if(filter_type==2){

        cout<< "lowpass";
    for (int i = 0; i < rows; i++){

				for (int j = 0; j < columns; j++)

				{

					double D = sqrt(pow(rows / 2 - i, 2) + pow(columns / 2 - j, 2));//do + x and y with i and j for ROI

					if ((D) <= thresh1)

					{

						filtered(i, j) = freqCentered(i,j);

					}

					else

					{

						filtered(i, j) = 0;

					}

				}

			}


        filtered.save(filteredName);

        //Do the folding to bring to the center - First and Third Quadrants

			for (int i = rows / 2; i< rows; i++){

				for (int j = columns / 2; j < columns; j++)

				{

					// for Freq

					int v1 = filtered(i, j);//Start Get Value from center Middle till end column

					int v2 = filtered(i - ((rows) / 2), j - ((columns) / 2));//start get value from  (0,0) till (0,column)

					filtered(i, j) = v2;

					filtered((i - (rows / 2)), j - (columns / 2)) = v1;

				}

			}



			//Second and Fourth Quadrants

			for (int i = (rows) / 2; i < rows; i++){

				for (int j = 0; j < (columns) / 2; j++)

				{

					// for frquency image.

					int v1 = filtered(i, j);

					int v2 = filtered(i - ((rows) / 2), j + ((columns) / 2));

					filtered(i, j) = v2;

					filtered((i - (rows / 2)), j - (columns / 2)) = v1;

				}

			}





			//Multiply the filter values with the FFT OUT

			//outDFT



			for (int i = 0; i < rows; i++){

				for (int j = 0; j < columns; j++)

				{

					if (filtered(i, j) == 0)

					{

						outR[i*columns + j][0] = 0;

						outR[i*columns + j][1] = 0;

					}

				}

			}
						filter_type=1;

}


if(filter_type==3){


        cout<< "highpass";
    for (int i = 0; i < rows; i++){

				for (int j = 0; j < columns; j++)

				{

					double D = sqrt(pow(rows / 2 - i, 2) + pow(columns / 2 - j, 2));//do + x and y with i and j for ROI

					if ((D) >= thresh1)

					{

						filtered(i, j) = freqCentered(i,j);

					}

					else

					{

						filtered(i, j) = 0;

					}

				}

			}


        filtered.save(filteredName);


			for (int i = rows / 2; i< rows; i++){

				for (int j = columns / 2; j < columns; j++)

				{

					int v1 = filtered(i, j);

					int v2 = filtered(i - ((rows) / 2), j - ((columns) / 2));

					filtered(i, j) = v2;

					filtered((i - (rows / 2)), j - (columns / 2)) = v1;

				}

			}


			for (int i = (rows) / 2; i < rows; i++){

				for (int j = 0; j < (columns) / 2; j++)

				{


					int v1 = filtered(i, j);

					int v2 = filtered(i - ((rows) / 2), j + ((columns) / 2));

					filtered(i, j) = v2;

					filtered((i - (rows / 2)), j - (columns / 2)) = v1;

				}

			}

			for (int i = 0; i < rows; i++){

				for (int j = 0; j < columns; j++)

				{

					if (filtered(i, j) == 0)

					{

						outR[i*columns + j][0] = 0;

						outR[i*columns + j][1] = 0;

					}

				}

			}
						filter_type=1;

}



if(filter_type==4)
{
cout << "\nBandpass\n";

			for (int i = 0; i < rows; i++){

				for (int j = 0; j < columns; j++)

				{

					double D = sqrt(pow(rows / 2 - i, 2) + pow(columns / 2 - j, 2));

					if ((D) >= thresh1 && (D) <= thresh2)
					{

						filtered(i, j) = freqCentered(i, j);

					}

					else

					{

						filtered(i, j) = 0;

					}

				}
			}
filtered.save(filteredName);




			for (int i = rows / 2; i< rows; i++){

				for (int j = columns / 2; j < columns; j++)

				{

					int v1 = filtered(i, j);

					int v2 = filtered(i - ((rows) / 2), j - ((columns) / 2));

					filtered(i, j) = v2;

					filtered((i - (rows / 2)), j - (columns / 2)) = v1;

				}

			}


			for (int i = (rows) / 2; i < rows; i++){

				for (int j = 0; j < (columns) / 2; j++)

				{


					int v1 = filtered(i, j);

					int v2 = filtered(i - ((rows) / 2), j + ((columns) / 2));

					filtered(i, j) = v2;

					filtered((i - (rows / 2)), j - (columns / 2)) = v1;

				}

			}

			for (int i = 0; i < rows; i++){

				for (int j = 0; j < columns; j++)

				{

					if (filtered(i, j) == 0)

					{

						outR[i*columns + j][0] = 0;

						outR[i*columns + j][1] = 0;

					}

				}

			}
						filter_type=1;

}


if(filter_type==5)
{
cout << "\nBandstop\n";

			for (int i = 0; i < rows; i++){

				for (int j = 0; j < columns; j++)

				{

					double D = sqrt(pow(rows / 2 - i, 2) + pow(columns / 2 - j, 2));

					if ((D) <= thresh1 || (D) >= thresh2)

					{

						filtered(i, j) = freqCentered(i, j);

					}

					else

					{

						filtered(i, j) = 0;

					}

				}
			}

filtered.save(filteredName);




			for (int i = rows / 2; i< rows; i++){

				for (int j = columns / 2; j < columns; j++)

				{

					int v1 = filtered(i, j);

					int v2 = filtered(i - ((rows) / 2), j - ((columns) / 2));

					filtered(i, j) = v2;

					filtered((i - (rows / 2)), j - (columns / 2)) = v1;

				}

			}


			for (int i = (rows) / 2; i < rows; i++){

				for (int j = 0; j < (columns) / 2; j++)

				{


					int v1 = filtered(i, j);

					int v2 = filtered(i - ((rows) / 2), j + ((columns) / 2));

					filtered(i, j) = v2;

					filtered((i - (rows / 2)), j - (columns / 2)) = v1;

				}

			}

			for (int i = 0; i < rows; i++){

				for (int j = 0; j < columns; j++)

				{

					if (filtered(i, j) == 0)

					{

						outR[i*columns + j][0] = 0;

						outR[i*columns + j][1] = 0;

					}

				}

			}
			filter_type=1;
}



filterPlan = fftw_plan_dft_2d(rows, columns, outR, outRfinal, FFTW_BACKWARD, FFTW_ESTIMATE);

		fftw_execute(filterPlan);

if(filter_type==1)
{

    for (int i = 0; i < rows; i++){

			for (int j = 0; j < columns; j++)

			{

				int k = i*(columns) + j;

				// normalize values
				double realPart = outRfinal[k][0] / (double)(rows * columns);
				double imagPart = outRfinal[k][1] / (double)(rows * columns);

				// magnitude
				double magnitude = sqrt((realPart * realPart) + (imagPart * imagPart));

                 Red = src(i,j,0);
                Green = src(i,j,1);
                Blue = src(i,j,2);

                RGBToHSI(Red,Green,Blue,Hue,Saturation,Intensity);

                Intensity = magnitude;

                HSIToRGB(Hue,Saturation,Intensity,Red,Green,Blue);
                tgt(i,j,0) = Red;
                tgt(i,j,1) = Green;
                tgt(i,j,2) = Blue;

				//tgt(i, j) = (int)magnitude;

			}

		}
}


 fftw_destroy_plan(planR);

    fftw_free(inR); fftw_free(outR);


}

//function body;
int main(int argc, char* argv[]){


int rows;
int columns;
Image src;
 Image magDFT;
 Image centered;
 Image freqDFT;
 Image freqCentered;
 Image filtered;
 Image invDFT;
 Image tgt;
int filter_type;


 char srcName[250];
 char magDFTName[250];
  char centeredName[250];
char freqDFTName[250];
  char freqCenteredName[250];

  char filteredName[250];
char invDFTName[250];
  char tgtName[250];

	if (argc<2)
	{
		  fprintf(stderr, "Please tell me the filter's name ^^\n");
		  exit(0);


	}
 std::string arg = argv[1];

    if(arg == "fft") {
            strcpy(srcName,argv[2]);
                strcpy(magDFTName,argv[3]);
    strcpy(centeredName,argv[4]);
                  strcpy(freqDFTName,argv[5]);
                        strcpy(freqCenteredName,argv[6]);
    //strcpy(invDFTName,argv[7]);


    strcpy(tgtName,argv[7]);
    //strcpy(tgtName,argv[5]);
    src.read(srcName);

        magDFT.initialize(src.NR,src.NC);
        centered.initialize(src.NR,src.NC);
        freqDFT.initialize(src.NR,src.NC);
                freqCentered.initialize(src.NR,src.NC);
        invDFT.initialize(src.NR,src.NC);

      tgt.initialize(src.NR,src.NC);
    rows = src.NR;
    columns= src.NC;
        // FORWARD TRANSFORM
        fft(rows, columns,src, magDFT,centered,freqDFT,freqCentered,tgt);
        magDFT.save(magDFTName);
        centered.save(centeredName);
                invDFT.save(invDFTName);


    }
    else if (arg == "filtering"){
        // BACKWARD TRANSFORM

        strcpy(srcName,argv[2]);
                strcpy(magDFTName,argv[3]);
    strcpy(centeredName,argv[4]);
                    strcpy(freqDFTName,argv[5]);
                        strcpy(freqCenteredName,argv[6]);
    strcpy(filteredName,argv[7]);

    //strcpy(invDFTName,argv[8]);


    strcpy(tgtName,argv[8]);
    filter_type=atoi(argv[9]);
    src.read(srcName);

        magDFT.initialize(src.NR,src.NC);
        centered.initialize(src.NR,src.NC);
         freqDFT.initialize(src.NR,src.NC);
                freqCentered.initialize(src.NR,src.NC);
                filtered.initialize(src.NR,src.NC);

                        invDFT.initialize(src.NR,src.NC);


      tgt.initialize(src.NR,src.NC);
    rows = src.NR;
    columns= src.NC;
        // FORWARD TRANSFORM
//fft(rows, columns,src, magDFT,centered,freqDFT,freqCentered,tgt);







               filtering(rows, columns,src, magDFT,centered,freqDFT,freqCentered,filtered,tgt,filter_type,filteredName);
                magDFT.save(magDFTName);
        centered.save(centeredName);
                //filtered.save(filteredName);
                               // invDFT.save(invDFTName);
                                tgt.save(tgtName);


    }



     else if (arg == "rgb"){
        // BACKWARD TRANSFORM

        strcpy(srcName,argv[2]);
                strcpy(magDFTName,argv[3]);
    strcpy(centeredName,argv[4]);
                    strcpy(freqDFTName,argv[5]);
                        strcpy(freqCenteredName,argv[6]);
    strcpy(filteredName,argv[7]);

    //strcpy(invDFTName,argv[8]);


    strcpy(tgtName,argv[8]);
    filter_type=atoi(argv[9]);
    src.read(srcName);

        magDFT.initialize(src.NR,src.NC);
        centered.initialize(src.NR,src.NC);
         freqDFT.initialize(src.NR,src.NC);
                freqCentered.initialize(src.NR,src.NC);
                filtered.initialize(src.NR,src.NC);

                        invDFT.initialize(src.NR,src.NC);


      tgt.initialize(src.NR,src.NC);
    rows = src.NR;
    columns= src.NC;
        // FORWARD TRANSFORM
//fft(rows, columns,src, magDFT,centered,freqDFT,freqCentered,tgt);







               rgbimage(rows, columns,src, magDFT,centered,freqDFT,freqCentered,filtered,tgt,filter_type,filteredName);
                magDFT.save(magDFTName);
        centered.save(centeredName);
                //filtered.save(filteredName);
                               // invDFT.save(invDFTName);
                                tgt.save(tgtName);


    }



	//fprintf(stderr, "Unknow Filter Name\n");
	return 0;

}
