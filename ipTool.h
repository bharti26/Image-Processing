// ipTool.h: interface for the ipTool class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_IPTOOL_H__65F5F541_0319_42DF_9ABC_143FAB67972C__INCLUDED_)
#define AFX_IPTOOL_H__65F5F541_0319_42DF_9ABC_143FAB67972C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "image.h"
#include "ROI.h"
class ipTool  
{
public:
	static void OutputROI(Image& src, ROI roi,Image& tgt);
	static void binarize(Image& src, Image& tgt, int thresh);
	ipTool();
	virtual ~ipTool();

};

#endif // !defined(AFX_IPTOOL_H__65F5F541_0319_42DF_9ABC_143FAB67972C__INCLUDED_)
