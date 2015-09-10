#include "Problem2.h"

void Problem2Exact1::EnumerateOneNode(vector<NodeSet *> *combination, bool isLeaf)
{	
	vector<NodeSet *> *L = setLists[0];
	vector<NodeSet *>::iterator ni;
	for(ni = L->begin(); ni != L->end(); ++ni)				
		combination->push_back(*ni);
	
	set<long long> mutexpairs;			//two nodes that cannot combine, key: node1.id + node2.id
	set<long long>::iterator pairIter;
	map<long long, double> pairDist;

	int maxSize = keywords.size();
	for(unsigned int i=2;i<=maxSize;i++)
	{			
		vector<NodeSet *> *Li = new vector<NodeSet *>;	//Li is a list, store all possible nodesets with size i

		KEYTYPE levelBM = 0;							//the keywords covered by this level

		for(unsigned int m=0;m<L->size();m++)
		{
			NodeSet * p1 = (*L)[m];
			for(unsigned int n=m+1;n<L->size();n++)
			{
				NodeSet * p2 = (*L)[n];

				KEYTYPE ubm = p1->bitmap | p2->bitmap;
				int num = MyTool::getNumOf1(ubm);
				if( num < i)		//E.g., three nodes ABC, but only covers 2 words, thus can be ignored
					continue;

				bool sharePrefix = true;		//check if two nodesets with size n share the same (n-1)nodes
				for(int k=0;k<i-2;k++)
				{
					if( (p1->nodes)[k] != (p2->nodes)[k])
					{
						sharePrefix = false;
						break;
					}
				}
				if(!sharePrefix)		//nodesets are ordered alphabetically, thus should be break
					break;//continue;

				RTreeNode *rtp1 = (p1->nodes)[i-2];
				RTreeNode *rtp2 = (p2->nodes)[i-2];

				if( i==2 )		//find all mutexpairs
				{
					long long key =  ( (long long)rtp1->identifier << 32 ) + rtp2->identifier;

					double diam;
					if(!isLeaf)
						diam = MyTool::ComputeMBRDist(rtp1->pr, rtp2->pr);
					else
					{
						double dx = rtp1->pr->m_pLow[0] - rtp2->pr->m_pLow[0];
						double dy = rtp1->pr->m_pLow[1] - rtp2->pr->m_pLow[1];
						diam = sqrt(dx *dx + dy * dy);
					}
					double maxDist = p1->maxDist > p2->maxDist ? p1->maxDist : p2->maxDist;
					pairDist[key] = diam;

					double minCost = (1-alpha) * diam + alpha * maxDist;
					if(minCost < cost)
					{
						NodeSet *p = new NodeSet();
						p->nodes.insert(p->nodes.begin(), p1->nodes.begin(), p1->nodes.end());
						p->nodes.push_back( rtp2 );
						p->bitmap = ubm;
						p->maxDiam = diam;
						p->maxDist = maxDist;
						p->minCost = minCost;
						Li->push_back(p);
						levelBM = levelBM | ubm;
					}
					else		//the two pairs cannot be combined
					{	
						mutexpairs.insert(key);
					}
				}

				else		//E.g., given ABC and ABD, if C and D are mutexpair, then the two nodesets cannot be combined							
							//However, even though, ABCD may still not be a candidate nodeset
							//because the maxDiam could be the distance between C and D, which may be larger
				{
					long long key =  ( (long long)rtp1->identifier << 32 ) + rtp2->identifier;
					if( (pairIter = mutexpairs.find(key)) != mutexpairs.end() )
						continue;

					double maxDist = p1->maxDist > p2->maxDist ? p1->maxDist : p2->maxDist;
					double maxDiam = p1->maxDiam > p2->maxDiam ? p1->maxDiam : p2->maxDiam;
					double diam = pairDist[key];
					maxDiam = diam > maxDiam ? diam : maxDiam;
					double minCost = (1-alpha) * maxDiam + alpha * maxDist;
					if(minCost < cost)
					{
						NodeSet *p = new NodeSet();
						p->nodes.insert(p->nodes.begin(), p1->nodes.begin(), p1->nodes.end());
						p->nodes.push_back( rtp2 );
						p->bitmap = ubm;
						p->maxDiam = maxDiam;
						p->maxDist = maxDist;
						p->minCost = minCost;
						Li->push_back(p);
						levelBM = levelBM | ubm;
					}
				}					
			}
		}

		delete L; L = Li;
		//注意delete L的时候已经把setLists中申请的空间释放了, 但是其中的nodeset, 会在后续过程释放combination时delete

		if(levelBM == keywordsBM)		//insert all possible nodesets
			combination->insert(combination->begin(), Li->begin(), Li->end());
		else				//if this level cannot cover all query keywords, we can stop		
		{
			vector<NodeSet *>::iterator iter = Li->begin();
			for(; iter != Li->end(); ++iter)
				delete *iter;
			delete Li;
			return;
		}
	}
	delete L;
}


void Problem2Exact1::EnumerateMultiNode(vector<NodeSet *> *combination, bool isLeaf)
{
	int min = 1<<31-1, rare;
	for(int i=0;i<keywords.size(); i++)			//find the word coveres by the least nodes
	{
		int temp = inverted[i].size();
		if(temp < min)
		{
			min = temp;
			rare = i;
		}
	}

	vector< vector<NodeSet *> * >::iterator iter = setLists.begin();		
	map<int, int> node2idx;
	vector<int> sizes;
	vector<NodeSet *> *L = new vector<NodeSet *>();
	int index = 0;

	for(; iter != setLists.end(); ++iter)		//get each node's child nodes, and filter out some of them
	{
		vector<NodeSet *> *p = *iter;
		vector<NodeSet *>::iterator ni;
		int count = 0;
		for(ni = p->begin(); ni!= p->end(); ++ni)
		{
			RTreeNode *rtp = (*ni)->nodes[0];

			bool satisfied = false;
			for(int i=0;i<inverted[rare].size();i++)	//compute the distance between a node and each node covering rare word
			{
				NodeSet *t = inverted[rare][i];
				if( (*ni) == t)
				{
					count++;
					satisfied = true; 
					break;
				}
				RTreeNode *rtp2 = t->nodes[0];
				double maxDist = (*ni)->maxDist > t->maxDist ? (*ni)->maxDist : t->maxDist;
				double diam;
				if(!isLeaf)
					diam = MyTool::ComputeMBRDist(rtp->pr, rtp2->pr);
				else
				{
					double dx = rtp->pr->m_pLow[0] - rtp2->pr->m_pLow[0];
					double dy = rtp->pr->m_pLow[1] - rtp2->pr->m_pLow[1];
					diam = sqrt(dx * dx + dy * dy);
				}
				if((1-alpha) * diam + alpha * maxDist < cost)		//can be combined with one node is ok
				{
					count++;
					satisfied = true;
					break;
				}
			}
			if(satisfied)		//insert this node into L
			{
				int nid = rtp->identifier;
				node2idx[nid] = index;		//map the node' ID to its parent node's position
				L->push_back(*ni);
			}
		}
		sizes.push_back(count);
		delete p;		//setLists里的vector *先释放, 其内容转到了L中
		index ++;
	}//the first filtering step, must at least can be combined with one node covering the rare word (do more?)

	vector<NodeSet *> *Li = new vector<NodeSet *>;		

	KEYTYPE levelBM = 0;

	set<long long> mutexpairs;
	set<long long>::iterator pairIter;
	map<long long, double> pairDist;
	
	for(unsigned int m=0;m<L->size();m++)
	{
		NodeSet * p1 = (*L)[m];

		for(unsigned int n=m+1;n<L->size();n++)
		{
			NodeSet * p2 = (*L)[n];

			KEYTYPE ubm = p1->bitmap | p2->bitmap;
			int num = MyTool::getNumOf1(ubm);
			if( num < 2)
				continue;

			RTreeNode *rtp1 = (p1->nodes)[0];
			RTreeNode *rtp2 = (p2->nodes)[0];
			long long key = ( (long long)rtp1->identifier << 32 ) + rtp2->identifier;

			double diam;
			if(!isLeaf)
				diam = MyTool::ComputeMBRDist(rtp1->pr, rtp2->pr);
			else
			{
				double dx = rtp1->pr->m_pLow[0] - rtp2->pr->m_pLow[0];
				double dy = rtp1->pr->m_pLow[1] - rtp2->pr->m_pLow[1];		
				diam = sqrt(dx *dx + dy * dy);
			}
			pairDist[key] = diam;
			double maxDist = p1->maxDist > p2->maxDist ? p1->maxDist : p2->maxDist;
			double minCost = (1-alpha) * diam + alpha * maxDist;
			if(minCost < cost)
			{
				NodeSet *p = new NodeSet();
				p->nodes.insert(p->nodes.begin(), p1->nodes.begin(), p1->nodes.end());
				p->nodes.push_back( (p2->nodes)[0]);
				p->bitmap = ubm;
				p->maxDiam = diam;
				p->maxDist = maxDist;
				p->minCost = minCost;
				Li->push_back(p);
				levelBM = levelBM | ubm;
			}		
			else
			{
				mutexpairs.insert(key);
			}
		}
	}
	vector<NodeSet *>::iterator ni = L->begin();
	for(; ni != L->end(); ++ni)
	{	
		delete *ni;					
	}
	delete L; L = Li;

	if(levelBM != keywordsBM)
	{
		vector<NodeSet *>::iterator ni = Li->begin();
		for(; ni != Li->end(); ++ni)
			delete *ni;
		delete Li;
		return;
	}

	for(unsigned int i=3;i<=keywords.size();i++)
	{
		vector<NodeSet *> *Li = new vector<NodeSet *>;	

		levelBM = 0;

		for(unsigned int m=0;m<L->size();m++)
		{
			NodeSet * p1 = (*L)[m];
			
			for(unsigned int n=m+1;n<L->size();n++)
			{
				NodeSet * p2 = (*L)[n];

				KEYTYPE ubm = p1->bitmap | p2->bitmap;
				int num = MyTool::getNumOf1(ubm);
				if( num < i)
					continue;

				bool sharePrefix = true;		//could be improved!
				for(int k=0;k<i-2;k++)
				{
					if( (p1->nodes)[k] != (p2->nodes)[k])
					{
						sharePrefix = false;
						break;
					}
				}
				if(!sharePrefix)
					break;

				int lastID = (p2->nodes)[i-2]->identifier;		//get ID of the last node in one nodeset to be combined

				int key_num = keywords.size(); int node_num = setLists.size(); 
				int lastPossibleNode = node_num - key_num + i - 1;
				if(node2idx[lastID] < lastPossibleNode)
					//This is to make sure that from the rest node, we can select at least one child node from each of them.
					//For example, A(A1, A2) B(B1, B2) C(C1, C2, C3) D(D1, D2), five keywords t1, t2, t3, t4, t5
					//If we are trying to combine A1A2B1 and A1A2B2, we can get A1A2B1B2, however at least 4 keywords are already
					//covered (each node must cover at least one), and thus only 1 word is left to be covered by C and D.
					//This means that either C or D is not necessary, which contradicts to our assumtion. 
					//Actually, this should be enumerated in the node set A(A1, A2), B(B1, B2), C(C1, C2, C3) {or D(D1, D2)}
					continue;

				RTreeNode *rtp1 = (p1->nodes)[i-2];
				RTreeNode *rtp2 = (p2->nodes)[i-2];
				long long key = ( (long long)rtp1->identifier << 32 ) + rtp2->identifier;

				if( !mutexpairs.empty() && (pairIter = mutexpairs.find(key)) != mutexpairs.end())
					continue;

				double diam = pairDist[key];		
				double maxDist = p1->maxDist > p2->maxDist ? p1->maxDist : p2->maxDist;
				double maxDiam = p1->maxDiam > p2->maxDiam ? p1->maxDiam : p2->maxDiam;
				maxDiam = diam > maxDiam ? diam : maxDiam;
				double minCost = (1-alpha) * maxDiam + alpha * maxDist;					
				if(minCost < cost)
				{
					NodeSet *p = new NodeSet();
					p->nodes.insert(p->nodes.begin(), p1->nodes.begin(), p1->nodes.end());
					p->nodes.push_back( (p2->nodes)[i-2]);
					p->bitmap = ubm;
					p->maxDiam = maxDiam;
					p->maxDist = maxDist;
					p->minCost = minCost;
					Li->push_back(p);
					levelBM = levelBM | ubm;
				}					
			}
		}

		vector<NodeSet *>::iterator ni = L->begin();
		for(; ni != L->end(); ++ni)
		{
			NodeSet *p = *ni;
			int lastID = p->nodes.back()->identifier;
			if(node2idx[lastID] == setLists.size()-1 && p->nodes.size() >= setLists.size())	
					//if the last node is from the last parent node:; A1A2B1B2C1 is not allowed
				combination->push_back(p);
			else
				delete p;
		}
		delete L;	L = Li;

		if( levelBM != keywordsBM)		//all possible nodesets do not cover all keywords
		{
			for(ni = Li->begin(); ni != Li->end(); ++ni)
				delete *ni;
			delete Li;
			return;
		}
	}

	//do not forget to process the last!
	ni = L->begin();
	for(; ni != L->end(); ++ni)
	{
		NodeSet *p = *ni;
		combination->push_back(p);
	}
	delete L;

}



void Problem2Exact1::SelectNodeSet()
{
	vector<NodeSet *> *p = new vector<NodeSet *>();

	if(setLists.size() == 1)		//only one node covers all the query keywords
		EnumerateOneNode(p);
	else
		EnumerateMultiNode(p);

	vector<NodeSet *>::iterator iter;
	for(iter = p->begin(); iter != p->end(); ++iter)
	{
		NodeSet *n = *iter;
		if(n->bitmap == keywordsBM)
			queueNodeSet.push(n);
		else
			delete n;
	}

	delete p;			
}

void Problem2Exact1::ExhaustiveSearch(NodeSet *res, vector<NodeSet *> **lists, int size)	//find the result recursively
{
	if(size == 0)
	{
		if(res->minCost <= cost)
		{
			cost = res->minCost;
			if(optimal != res)
			{
				delete optimal;
				optimal = res;
			}
		}
		return;
	}

	int min = 1<<31-1, rare;
	for(int i=0;i<size; i++)
	{
		int temp = lists[i]->size();
		if(temp < min)
		{
			min = temp;
			rare = i;
		}
	}

	vector<NodeSet *> ** newlists = new vector<NodeSet *>* [size-1];
	int count = 0;
	for(int i=0;i<size;i++)
	{
		if(i != rare)
			newlists[count++] = lists[i];
	}

	vector<NodeSet *>::iterator iter = lists[rare]->begin();
	for(; iter != lists[rare]->end(); ++iter)
	{
		NodeSet *p = *iter;
		double maxDist = p->maxDist > res->maxDist ? p->maxDist : res->maxDist;
		vector<RTreeNode *>::iterator ri = res->nodes.begin();
		double maxDiam = 0;
		for(; ri != res->nodes.end(); ++ri)
		{
			RTreeNode *t = *ri;
			double dx = t->pr->m_pLow[0] - p->nodes[0]->pr->m_pLow[0];
			double dy = t->pr->m_pLow[1] - p->nodes[0]->pr->m_pLow[1];
			double diam = sqrt(dx * dx + dy * dy);
			if(diam > maxDiam)
				maxDiam = diam;
		}
		if(res->maxDiam > maxDiam)
			maxDiam = maxDiam;
		double minCost = alpha * maxDist + (1-alpha) * maxDiam;
		if(minCost < cost)
		{
			double tempDiam = res->maxDiam;
			double tempDist = res->maxDist;
			double tempminCost = res->minCost;
			res->maxDiam = maxDiam;
			res->maxDist = maxDist;
			res->minCost = minCost;
			res->nodes.push_back(p->nodes[0]);
			ExhaustiveSearch(res, newlists, size-1);		//begin the recursion
			res->nodes.pop_back();
			res->maxDiam = tempDiam;
			res->maxDist = tempDist;
			res->minCost = tempminCost;
		}
		else 
			return;
	}
	delete newlists;
}


void Problem2Exact1::getNextEntry(const IEntry& entry, id_type& nextEntry, bool& hasNext)
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
		
		map<int, vector<int>*> objectTexts;		//the keywords contained in a node
		set<int> indexID;						//the id of child nodes which contain query keywords				
		map<int, vector<int>*>::iterator iterMap;
		int wsize = keywords.size();			
		for(unsigned int k=0;k<wsize;k++)
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
									vector<int> *p = new vector<int>();		// delete in RTreeNode, not a good way...
									p->push_back(k);		//the index of the word in keywords, not the real word id!
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

		vector<RTreeNode *> *candidateNodes = new vector<RTreeNode *>();
		set<int>::iterator si = indexID.begin();
		for(;si != indexID.end(); ++si)	
		{
			uint32_t cChild = *si;
			int cid = n->getChildIdentifier(cChild);
			
			IShape *out;
			n->getChildShape(cChild, &out);
			Region* pr = dynamic_cast<Region *>(out);				
			
			map<int, RTreeNode *>::iterator oi = nodesBuffer.find(cid);
			if( oi != nodesBuffer.end())	//this node is in the buffer, read it directly
			{
				RTreeNode *rtp = oi->second;
				if(alpha * rtp->minValue > cost)
					continue;
				candidateNodes->push_back(rtp);
			}
			else
			{
				double dist;
				if(n->isIndex())
				{
					double coor[2];
					coor[0] = Q->x; coor[1] = Q->y;
					Point pt(coor, 2);
					dist = MyTool::ComputeMinPossible(pr, &pt);
				}
				else
				{
					double dx = Q->x - pr->m_pLow[0];
					double dy = Q->y - pr->m_pLow[1];
					dist = sqrt ( dx * dx + dy * dy ); 
				}

				if( alpha * dist > cost)		//the node cannot contribute final results
					continue;					
				
				RTreeNode *rtp = new RTreeNode();
				rtp->identifier = n->getChildIdentifier(cChild);
				rtp->minValue = dist;
				rtp->bitmap = MyTool::ConvertToInt(objectTexts[cChild]);
				rtp->pr = pr;

				candidateNodes->push_back(rtp);
				nodesBuffer[rtp->identifier] = rtp;
			}
		}				
		
		vector<NodeSet *> *L = new vector<NodeSet *>;						//注意什么时候delete!!
		vector<RTreeNode *>::iterator iter = candidateNodes->begin();		//convert RTreeNode to NodeSet
		for(; iter != candidateNodes->end(); ++iter)
		{
			NodeSet *p = new NodeSet();
			p->nodes.push_back(*iter);
			p->maxDist = (*iter)->minValue;
			p->bitmap = (*iter)->bitmap;
			p->minCost = alpha * (*iter)->minValue;
			L->push_back(p);

			int mask = 1;
			for(int i=0; i<keywords.size(); i++, mask = mask<<1)
			{
				if( (p->bitmap & mask) > 0)
					inverted[i].push_back(p);
			}
		}
		setLists.push_back(L);
		delete candidateNodes;
		
		if( !(ns->nodes.empty()) )		//we are only retrieving a node in ns, we keep reading, until all nodes in ns are read
		{
			id_type page = ns->nodes.back()->identifier;
			ns->nodes.pop_back();
			hasNext = true;
			nextEntry = page;
		}
		else
		{					
			if(n->isIndex())			//if not reaching leave nodes, we only generate NodeSet
				SelectNodeSet();
			else						//now we begin to enumerate the results
			{
				NodeSet *res = new NodeSet();
				int size = keywords.size();
				vector<NodeSet *> ** lists = new vector<NodeSet *>* [size];
				for(int i=0;i<size;i++)
				{
					lists[i] = &(inverted[i]);
				}
				ExhaustiveSearch(res, lists, size);
				delete lists;//*/
			}

			setLists.clear();
			for(int i=0;i<keywords.size();i++)
				inverted[i].clear();

			if(!queueNodeSet.empty())		//we still have nodeset in queue, we need to read nodes for the next one
			{
				delete ns;
				ns = queueNodeSet.top();
				queueNodeSet.pop();			

				if(ns->minCost > cost)
				{
					group.clear();
					vector<RTreeNode *>::iterator nodes = (optimal->nodes).begin();
					for(; nodes != (optimal->nodes).end(); nodes++)
						group.insert( (*nodes)->identifier);
					hasNext = false;
				}
				else
				{
					id_type page = ns->nodes.back()->identifier;
					ns->nodes.pop_back();
					hasNext = true;
					nextEntry = page;
				}
			}
			else
			{	
				hasNext = false;
			}
		}
	}

	else
		hasNext = false;
}
