// Copyright (c) 2005-2012 Shi-Qing Xin (xinshiqing@163.com) and Guo-Jin Wang (wanggj@zju.edu.cn).
// NOTE: this is an old version. For the lastest version, please email us.
// The code is free for research purpose, but requires a token charge for commercial purpose. 
// Users are forbidden to reproduce, republish, redistribute, or resell the code without our permission.
//
// In this code, we implemented chen and han's algorithm [1990].
// Furthermore, we gave two techniques to improve the CH algorithm.
// If you have any ideas about improving the code, we are very grateful.
// My personal website: http://sites.google.com/site/xinshiqing/Home
// 
// We are debted to Prof. Han, who gave us many helpful ideas.
// We must thank Surazhsky and Kirsanov for their knowledge share.
//
// xin_wang_geodesic.cpp : program entry.
//

#include "stdafx.h"
#include <tchar.h>
#include "BaseModel.h"
#include "RichModel.h"
#include "ExactMethodForDGP.h"
#include "PreviousCH.h"
#include "ImprovedCHWithFilteringRule.h"
#include "XinWangImprovedCH.h"


int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		cout << "No enough input parameters! input your parameters like this: [filename]." << endl;
		return 0;
	}
	else if (argc > 2)
	{
		cout << "Useless parameters exist. Please check again!" << endl;
	}
	cout << "Please use release mode!\n";

	CRichModel model(argv[1]);
	cout <<"----------------------model info begin----------------\n";	
	cout << "File name:\t" << model.GetFileName() << endl;		
	try
	{
		model.LoadModel();
	}
	catch(const char* msg)
	{
		cout << "ERRORS happen!\n" << msg << endl;		
		return 1;
	}

	model.Preprocess();

	if (!model.HasBeenLoad() || !model.HasBeenProcessed())
	{
		cout << "The model fails to be handled." << endl;
		return 1;
	}
	

	cout << "Face number:\t" << model.GetNumOfFaces() << endl;
	cout << "Vertex number:\t" << model.GetNumOfVerts() << endl;
	cout << "Edge number:\t" << model.GetNumOfEdges() << endl;
	if (model.IsClosedModel())
		cout << "It is a closed model." << endl;
	else 
		cout << "It is an open model with " << model.GetNumOfHoles() << " holes." << endl;
	if (model.GetNumOfIsolated() > 1)
	{
		cout << "Perhaps it is composed of several components." << endl;
	}
	cout <<"----------------------model info end------------------\n";

	int indexOfSourceIndex = 0;
	cout << "Please input your source vertex index (0 ~ " << model.GetNumOfVerts() - 1 << "):  ";
	cin >> indexOfSourceIndex;
	if (indexOfSourceIndex < 0 || indexOfSourceIndex >=  model.GetNumOfVerts())
	{
		cout << "your input is wrong, so we assume the source index is zero." << endl;
		indexOfSourceIndex = 0;
	}
	CExactMethodForDGP * algorithm = new CXinWangImprovedCH(model, indexOfSourceIndex);
	
	cout << "----------------------Algorithm begins:----------------\n";
	cout << "Algorithm name:\t" << algorithm->GetAlgorithmName() << endl;
	algorithm->Execute();
	cout << "Running time (s):\t" << algorithm->GetRunTime() / 1000.0 << endl;
	cout << "Peak memory: (M)\t" << algorithm->GetMemoryCost() << endl;
	cout << "Total window number (including deleted windows):\t" << algorithm->GetWindowNum() << endl;
	cout << "The max size of queue:\t" << algorithm->GetMaxLenOfQue() << endl;
	cout << "The max depth (levels) of the sequence tree\t" <<algorithm->GetDepthOfSequenceTree() << endl;  
	cout << "----------------------Algorithm ends:------------------\n";
	cout << "----------------------NOTE:----------------\n";
	cout << "The peak memory is not equivalent to the total windows generated, since most of nodes are deleted during the execution of CH algorith. However, for the MMP algorithm, only a small portion of windows are deleted during the execution.\n";
	cout << "----------------------NOTE END----------\n";
	cout << "----------------------Backtrace a shortest path to a vertex---\n";
	cout << "Please input your favarite destination vertex index (0 ~ " << model.GetNumOfVerts() - 1 << ",  \'-1\' exit.)";
	
	int dest = 0;
	cin >> dest;
	if (!(dest < 0 || dest >=  model.GetNumOfVerts()))
	{
		cout << "Since we have scaled the model, so the output coordinate may be different from the original model.\n";
		vector<CPoint3D> resultpoints;
		algorithm->BackTrace(dest, resultpoints);	
		for (int i = resultpoints.size() -1; i >=0; --i)
		{
			cout << "(" << resultpoints[i].x << ", " << resultpoints[i].y << ", " << resultpoints[i].z << ")\n";
		}
	}	
	cout << "----------------------End backtracing-------------------------------\n";
	cout << "----------------------Backtrace a shortest path to a surface point---\n";
	cout << "NOTE: in fact, to backtrace a shortest path to a surface point, we can use a simple idea by Surazhsky et al. [2005]. ";
	cout << "After we finish this part of code, I will upload to my website. Please stay tuned.\n";	
	cout << "----------------------End backtracing-------------------------------\n";
	delete algorithm;
	algorithm = NULL;
	char ch;
	cout << "\nIf you want to exit? yes or no?\t" ;
	while (cin >> ch)
	{
		if (ch == 'y')
			break;
	}
	
	return 0;
}

