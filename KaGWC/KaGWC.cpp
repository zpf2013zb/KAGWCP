// SGKTEST.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <time.h>
#include <iostream>
#include <functional>

#include "IBRTree.h"
#include "SearchIBRTree.h"
#include "Problem1.h"
#include "Problem2.h"
using namespace std;

string basepath = "..\\dataHotel\\";

string locFile = basepath + "loc";
string treeFile = basepath + "rtree";
string textFile = basepath + "doc";
string invertedFile = basepath + "invertedfile";

string btreeFolder = basepath + "btreeindex\\";
string leafnodeFile = basepath + "leafnodes";
string indexnodeFile = basepath + "indexnodes";
string subdocFolder = basepath + "subdoc\\";
string invertedFolder = basepath + "inverted\\";

int numOfEntry = 40;		//For hotel data, the fanout of R-tree is 40. it must be consistent with NODENUM defined in data.h!
double alpha = 0.5;

int main()
{
	IBRTree irtree;
	//irtree.BuildIBRTree();
	//cout << 2 << endl;

	irtree.ReadTree();
	
	IStatistics *out;
	irtree.GetTree()->getStatistics(&out);
	int nodeNum = out->getNumberOfData();
	//cout << 1 << endl;

	Query *Q = new Query("2,5,100,52", 44, -151);

	clock_t start, finish;
	//cout << 2 << endl;

	//the following code shows how to retrieve the top-k nearest neighbor that contains some query keywords	
	
	/*
	SearchIBRTree ss(Q, nodeNum, 10);
	irtree.GetTree()->queryStrategy(ss);
	ISpatialIndex* tree= irtree.GetTree();
	
	map<double, int, greater<double>>::reverse_iterator iter;
	int num = 0;
	for(iter = ss.topk.rbegin(); iter!= ss.topk.rend(); ++iter)
	{
	cout<<++num<<" Point ID:"<<iter->second<<" Distance:"<<iter->first<<endl;
	}
	*/
	
	
	
	//the following code shows how to use the algorithms for the TYPE1 query
	start = clock();
	Problem1Appr pb1a(Q, nodeNum);
	irtree.GetTree()->queryStrategy(pb1a);
	finish = clock();
	//cout<<pb1a;	
	cout<< (finish - start) / CLOCKS_PER_SEC<<endl;

	/*
	start = clock();
	Problem1Baseline pb1e1(Q);
	pb1e1.BaseLine();
	finish = clock();
	cout<<pb1e1;
	cout<< (finish - start) / CLOCKS_PER_SEC<<endl;

	start = clock();
	Problem1IBRTree pb1e2(Q, nodeNum);
	irtree.GetTree()->queryStrategy(pb1e2);
	finish = clock();
	cout<<pb1e2;
	cout<< (finish - start) / CLOCKS_PER_SEC<<endl;
	//*/

	/*
	//the following code shows how to use the algorithms for the TYPE2 query
	start = clock();
	Problem2Appr1 pb2a1(Q, nodeNum);
	irtree.GetTree()->queryStrategy(pb2a1);
	finish = clock();
	cout<<pb2a1;
	cout<< (finish - start) / CLOCKS_PER_SEC<<endl;

	start = clock();
	Problem2Appr2 pb2a2(Q, nodeNum);
	irtree.GetTree()->queryStrategy(pb2a2);
	finish = clock();
	cout<<pb2a2;
	cout<< (finish - start) / CLOCKS_PER_SEC<<endl;

	start = clock();
	Problem2Exact1 pb2e(Q, nodeNum);
	irtree.GetTree()->queryStrategy(pb2e);
	finish = clock();
	cout<<pb2e;
	cout<< (finish - start) / CLOCKS_PER_SEC<<endl;		//*/

	return 0;
}

