#ifndef _PROBLEM2_
#define _PROBLEM2_

#include <iostream>
#include <SpatialIndex.h>
#include <map>
#include <set>
#include <vector>
#include "Tool.h"
#include "data.h"
#include "btree.h"
#include "IBRtree.h"

using namespace std;

extern string textFile;
extern string locFile;
extern string invertedFile;
extern string subdocFolder;
extern string btreeFolder;
extern string treeFile;

extern int numOfEntry;
extern double alpha;

class Problem2Appr1 : public SpatialIndex::IQueryStrategy
{
private:	
	priority_queue<RTreeNode *, deque<RTreeNode *>, NodeValueLess<RTreeNode> > queueNode;	
	Query *Q;
	int nodeNum;
	vector<int> keywords;
	KEYTYPE keywordsBM;	
	vector<int> group;
	double cost;
	map<int, Point> objLoc;
	bool usedInAPP2;
				
public:

	Query furthestObj;		//this structure is used for the second approximate algorithm
	
	Problem2Appr1(Query *Q, int nodeNum, bool usedInAPP2=false)
	{
		this->Q = Q;
		this->nodeNum = nodeNum;		

		istringstream iss(Q->text);
		int wid; char c;
		while(iss>>wid)
		{
			keywords.push_back(wid);
			iss>>c;
		}

		keywordsBM = pow(2.0, (int)keywords.size())-1;
		cost = 0;
		this->usedInAPP2 = usedInAPP2;
	}

	~Problem2Appr1()
	{
		while( !queueNode.empty())
		{
			RTreeNode *p = queueNode.top();
			queueNode.pop();
			delete p;
		}
	}

	friend ostream& operator<<(ostream& out,Problem2Appr1 & t)    
	{
		out<<t.getCost()<<":";
		vector<int>::iterator iter = t.group.begin();
		for(; iter != t.group.end(); ++iter)
			out<<*iter<<" ";
		out<<endl;
		return out;
	}

	void setQuery(Query *Q)
	{
		this->Q = Q;
	}

	double getCost(vector<int> *objs = NULL)
	{
		double maxDiam = 0;
		for(int i=0;i<group.size();i++)
		{
			int o1 = group[i];
			if (objs != NULL)	objs->push_back(o1);
			for(int j=i+1;j<group.size();j++)
			{
				int o2 = group[j];
				double x1 = objLoc[o1].m_pCoords[0], y1 = objLoc[o1].m_pCoords[1];
				double x2 = objLoc[o2].m_pCoords[0], y2 = objLoc[o2].m_pCoords[1];
				double dist = sqrt( (x1-x2) * (x1-x2) + (y1-y2) * (y1-y2));
				if(dist > maxDiam)
					maxDiam = dist;
			}
		}
		cost = alpha * cost + (1-alpha) * maxDiam;
		return cost;
	}

	vector<int> getGroup()
	{
		return group;
	}

	double getRes(Query *oriQ,vector<int> *objs)	//this function is used in the second approximate algorithm
	{
		double res=0, maxDist=0, maxDiam=0;		
		for(int i=0;i<group.size();i++)
		{
			int o1 = group[i];
			double x1 = objLoc[o1].m_pCoords[0], y1 = objLoc[o1].m_pCoords[1];
			objs->push_back(o1);
			double dist = sqrt( (x1-oriQ->x) * (x1-oriQ->x) + (y1-oriQ->y) * (y1-oriQ->y));
			if(dist > maxDist)
				maxDist = dist;

			for(int j=i+1;j<group.size();j++)
			{
				int o2 = group[j];
				double x2 = objLoc[o2].m_pCoords[0], y2 = objLoc[o2].m_pCoords[1];
				double dist = sqrt( (x1-x2) * (x1-x2) + (y1-y2) * (y1-y2));
				if(dist > maxDiam)
					maxDiam = dist;
			}
		}
		res = alpha * maxDist + (1-alpha) * maxDiam;
		return res;
	}

	void getNextEntry(const IEntry& entry, id_type& nextEntry, bool& hasNext)
	{	
		const INode* n = dynamic_cast<const INode*>(&entry);

		if (n != 0)
		{
			int pid = n->getIdentifier();
			string btFile = btreeFolder + MyTool::IntToString(pid);
			char *btfname = new char[btFile.size()+1];
			memcpy(btfname, btFile.c_str(), btFile.size());		
			btfname[btFile.size()] = '\0';
			BTree *bt = new BTree(btfname, 0);	
			
			map<int, vector<int>*> objectTexts;		//the text description of an object
			set<int> indexID;						//the id of child nodes which contains query keywords				
			map<int, vector<int>*>::iterator iterMap;
			for(unsigned int k=0;k<keywords.size();k++)
			{
				int wordID = keywords[k];
				VECTYPE *data = new VECTYPE[DIMENSION];
				bool flag = bt->search(wordID, &data);
				if(flag)
				{					
					for(int i=0;i<DIMENSION;i++)
					{
						if(data[i] > 0)							
						{
							unsigned char mask = 1;
							for(int j=0;j<8;j++)
							{
								if((data[i] & mask) > 0)
								{
									int index = i*8+j;
									indexID.insert(index);
									iterMap = objectTexts.find(index);
									if(iterMap ==  objectTexts.end())
									{
										vector<int> *p = new vector<int>();		
										p->push_back(k);
										objectTexts[index] = p;
									}
									else
									{
										vector<int> *p = iterMap->second;
										p->push_back(k);
									}
								}
								mask = mask << 1;
							}
						}
					}
				}
				delete data;
			}
			delete bt;
			delete btfname;

			set<int>::iterator si = indexID.begin();
			for(;si != indexID.end(); ++si)	
			{
				uint32_t cChild = *si;
				int cid = n->getChildIdentifier(cChild);
				
				IShape *out;
				n->getChildShape(cChild, &out);
				Region* pr = dynamic_cast<Region *>(out);	

				double dx = Q->x - pr->m_pHigh[0];
				double dy = Q->y - pr->m_pHigh[1];											

				int key = MyTool::ConvertToInt(objectTexts[cChild]);
				if( (key & keywordsBM) == 0)
					continue;

				double dist;
				RTreeNode *rtp = new RTreeNode;

				if(n->isLeaf())
				{
					dist = sqrt(dx * dx + dy * dy);
					rtp->isNode = false;
				}
				if(n->isIndex())
				{
					double coor[2];
					coor[0] = Q->x; coor[1] = Q->y;
					Point pt(coor, 2);

					dist = MyTool::ComputeMinPossible(pr, &pt);
					rtp->isNode = true;
				}
				
				rtp->identifier = n->getChildIdentifier(cChild);
				rtp->minValue = dist;
				rtp->bitmap = MyTool::ConvertToInt(objectTexts[cChild]);
				rtp->pr = pr;

				queueNode.push(rtp);
			}
			map<int, vector<int>*>::iterator oi = objectTexts.begin();
			for(; oi != objectTexts.end(); ++oi)
			{
				vector<int> *p = oi->second;
				delete p;
			}
		}

		if (!queueNode.empty())		
		{			
			while (!queueNode.empty())
			{
				RTreeNode *p = queueNode.top();
				queueNode.pop();

				if( p->isNode )		//if the node is not an object, we read the next page.
				{	
					KEYTYPE keyValue = p->bitmap;
					if( (keyValue & keywordsBM) > 0)
					{
						nextEntry = p->identifier;
						hasNext = true;
						delete p;
						break;
					}
				}
				else
				{
					KEYTYPE key = p->bitmap;
					int interS = key & keywordsBM;
					if(interS > 0)		//contains some keywords that have not been covered
					{
						group.push_back(p->identifier);
						keywordsBM = keywordsBM - interS;
						double coor[2];
						coor[0] = p->pr->m_pLow[0];
						coor[1] = p->pr->m_pLow[1];
						Point pp(coor,2);
						objLoc[p->identifier] = pp;

						if(keywordsBM == 0)
						{
							hasNext = false;
							cost = p->minValue;
							//the following code is for the second approximate algorithm
							if(usedInAPP2)
							{
								furthestObj.x = p->pr->m_pLow[0];
								furthestObj.y = p->pr->m_pLow[1];
								KEYTYPE mask = 1;
								for(int i=0;i<sizeof(KEYTYPE) * 8;i++)
								{
									if( (interS & mask) > 0)
									{
										furthestObj.text = MyTool::IntToString(keywords[i]);
										//any word is OK, used for finding next nearest object
										break;
									}
									mask = mask << 1;
								}
							}
							delete p;
							break;
						}
					}
				}
				delete p;
			}
		}
	
		else 
			hasNext = false;
	}
};



class Problem2Appr2: public SpatialIndex::IQueryStrategy
{
private:		
	priority_queue<RTreeNode *, deque<RTreeNode *>, NodeValueLess<RTreeNode> > queueNode;	
	Query *Q;
	int nodeNum;
	vector<int> *group;
	double cost;	
	IBRTree *ibrtree;
	Query fo;	

public:
	Problem2Appr2(Query *Q, int nodeNum)
	{
		this->Q = Q;
		this->nodeNum = nodeNum;		
		this->ibrtree = new IBRTree();
		ibrtree->ReadTree();

		Problem2Appr1 *pb2a1 = new Problem2Appr1(Q, nodeNum, true);
		ibrtree->GetTree()->queryStrategy(*pb2a1);
		fo = pb2a1->furthestObj;
		group = new vector<int>();
		cost = pb2a1->getCost(group);		
		delete pb2a1;
	}

	~Problem2Appr2()
	{
		delete group;
		delete ibrtree;
		while( !queueNode.empty())
		{
			RTreeNode *p = queueNode.top();
			queueNode.pop();
			delete p;
		}
	}

	double getCost()	{return cost;}
	vector<int> * getGroup()	{return group;}

	friend ostream& operator<<(ostream& out,Problem2Appr2 & t)   
	{
		out<<t.cost<<":";
		vector<int>::iterator iter = t.group->begin();
		for(; iter != t.group->end(); ++iter)
			out<<*iter<<" ";
		out<<endl;
		return out;
	}

	void getNextEntry(const IEntry& entry, id_type& nextEntry, bool& hasNext)
	{
		const INode* n = dynamic_cast<const INode*>(&entry);

		if(n!=0)
		{
			int pid = n->getIdentifier();
			string btFile = btreeFolder + MyTool::IntToString(pid);
			char *btfname = new char[btFile.size()+1];
			memcpy(btfname, btFile.c_str(), btFile.size());		
			btfname[btFile.size()] = '\0';
			BTree *bt = new BTree(btfname, 0);	

			set<int> indexID;	//the id of child nodes which contains query keywords						
			
			int wordID = atoi( (fo.text).c_str());
			VECTYPE *data = new VECTYPE[DIMENSION];
			bool flag = bt->search(wordID, &data);
			if(flag)
			{					
				for(int i=0;i<DIMENSION;i++)
				{
					if(data[i] > 0)							
					{
						unsigned char mask = 1;
						for(int j=0;j<8;j++)
						{
							if((data[i] & mask) > 0)
							{
								int index = i*8+j;
								indexID.insert(index);
							}
							mask = mask << 1;
						}
					}
				}
			}
			delete data;		
			delete bt;
			delete btfname;

			set<int>::iterator si = indexID.begin();
			for(;si != indexID.end(); ++si)	
			{
				uint32_t cChild = *si;
				int cid = n->getChildIdentifier(cChild);
				
				IShape *out;
				n->getChildShape(cChild, &out);
				Region* pr = dynamic_cast<Region *>(out);	

				double dx = Q->x - pr->m_pHigh[0];
				double dy = Q->y - pr->m_pHigh[1];											
				
				double dist;
				RTreeNode *rtp = new RTreeNode;

				if(n->isLeaf())
				{
					dist = sqrt(dx * dx + dy * dy);
					rtp->isNode = false;
				}
				if(n->isIndex())
				{
					double coor[2];
					coor[0] = Q->x; coor[1] = Q->y;
					Point pt(coor, 2);

					dist = MyTool::ComputeMinPossible(pr, &pt);
					rtp->isNode = true;
				}
				
				if(alpha * dist < cost)
				{
					rtp->identifier = n->getChildIdentifier(cChild);
					rtp->minValue = dist;			
					rtp->pr = pr;
					queueNode.push(rtp);
				}
				else 
					delete rtp;
			}
		
			while (!queueNode.empty())
			{
				RTreeNode *p = queueNode.top();
				queueNode.pop();

				if( p->isNode )		//if the node is not an object, we read the next page.
				{	
					nextEntry = p->identifier;
					hasNext = true;
					delete p;
					return;
				}
				else
				{
					if(alpha * p->minValue > cost)		//if its distance is already larger than the cost
						break;

					Query *newq = new Query();
					newq->text = Q->text;
					newq->x = p->pr->m_pLow[0];
					newq->y = p->pr->m_pLow[1];

					Problem2Appr1 pb2a1(newq, nodeNum, true);
					ibrtree->GetTree()->queryStrategy(pb2a1);
					vector<int> *objs = new vector<int>();
					double newCost = pb2a1.getRes(Q, objs);
					if(newCost < cost)
					{
						cost = newCost;
						group = objs;
					}
					else
						delete objs;
					
					delete newq;
				}
				delete p;
			}

			hasNext = false;
		}

		else
			hasNext = false;
	}	

};


class NodeSet{
public:
	vector<RTreeNode *> nodes;
	double minCost;
	double maxDiam;
	double maxDist;
	KEYTYPE bitmap;
       
	NodeSet()
	{
		minCost = maxDiam = maxDist = 0.0;
		bitmap = 0;
	}

   static bool CompareValueLess(const NodeSet * a,const  NodeSet * b)
   {
	   return a->minCost > b->minCost;
   }
};


class Problem2Exact1 : public SpatialIndex::IQueryStrategy
{
private:	
	priority_queue<NodeSet *, deque<NodeSet *>, NodeValueLess<NodeSet> > queueNodeSet;	
	Query *Q;
	int nodeNum;
	NodeSet *ns;				//the current processed node set
	NodeSet *optimal;			//final result
	vector<int> keywords;
	KEYTYPE keywordsBM;
	set<int> group;
	double cost;
	map<int, RTreeNode *> nodesBuffer;
	
	vector< vector<NodeSet *> * > setLists;		//存放每一个node的children nodes
	vector<NodeSet *> *inverted;		//the word contained in which node sets
				
public:
	
	Problem2Exact1(Query *Q, int nodeNum)
	{
		this->Q = Q;
		this->nodeNum = nodeNum;

		istringstream iss(Q->text);
		int wid; char c;
		while(iss>>wid)
		{
			keywords.push_back(wid);
			iss>>c;
		}

		keywordsBM = pow(2.0, (int)keywords.size())-1;
		ns = new NodeSet();
		optimal = NULL;

		inverted = new vector<NodeSet *>[keywords.size()];
		
		IBRTree *rt = new IBRTree();
		rt->ReadTree();
		Problem2Appr2 pb2a2(Q, nodeNum);
		rt->GetTree()->queryStrategy(pb2a2);
		cost = pb2a2.getCost();

		vector<int> *iniGroup = pb2a2.getGroup();
		vector<int>::iterator vi = iniGroup->begin();
		for(; vi != iniGroup->end(); ++vi)
		{
			group.insert(*vi);
		}
	}
	
	~Problem2Exact1()
	{
		map<int, RTreeNode *>::iterator iter;
		for(iter = nodesBuffer.begin(); iter != nodesBuffer.end(); ++iter)
		{
			delete iter->second;
		}
		while(!queueNodeSet.empty())
		{
			NodeSet *p = queueNodeSet.top();
			delete p;
			queueNodeSet.pop();
		}
		delete ns;
	}

	
	friend ostream& operator<<(ostream& out,Problem2Exact1 & t)   
	{
		out<<t.cost<<":";
		set<int>::iterator iter = t.group.begin();
		for(; iter != t.group.end(); ++iter)
			out<<*iter<<" ";
		out<<endl;
		return out;
	}

	void EnumerateOneNode(vector<NodeSet *> *combination, bool isLeaf = false);	
		//enumerate all possible nodeset from one node

	void EnumerateMultiNode(vector<NodeSet *> *combination, bool isLeaf = false);
		//enumerate all possible nodeset from multiple nodes

	void SelectNodeSet();

	void ExhaustiveSearch(NodeSet *res, vector<NodeSet *> **lists, int size);

	void getNextEntry(const IEntry& entry, id_type& nextEntry, bool& hasNext);
};


#endif