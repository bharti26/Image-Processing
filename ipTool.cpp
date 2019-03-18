// ipTool.cpp: implementation of the ipTool class.
//
//////////////////////////////////////////////////////////////////////

#include "ipTool.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ipTool::ipTool()
{

}

ipTool::~ipTool()
{

}

void ipTool::binarize(Image &src, Image &tgt, int thresh)
{

	for(int i=0;i<tgt.NR;i++){
		for(int j=0;j<tgt.NC;j++){		
			if (src(i,j,RED)>thresh)
				tgt(i,j,RED)=255;
			else
				tgt(i,j,RED)=0;
		}
	}	

}

void ipTool::OutputROI(Image &src, ROI roi, Image &tgt)
{
	for(int i=0;i<tgt.NR;i++){
		for(int j=0;j<tgt.NC;j++){		
			if (roi.InROI(i,j))
			{
				tgt(i,j,RED)=src(i,j,RED);
			}
			else
				tgt(i,j,RED)=0;
		}
	}	

}
