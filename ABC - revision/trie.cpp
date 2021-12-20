#include "trie.h"



Trie::Trie()
{
	this->key = 0; //root key is always 0
	//this->isLeaf = false; //see h file
	this->flag = false;
	//do I need to initialize the 'sons' vector?
	this->blockVec.clear();
}

Trie::Trie(int num)
{
	this->key = num; //root key is always 0
	//this->isLeaf = false; //see h file
	this->flag = false;
	//do I need to initialize the 'sons' vector?
	
}

void Trie::insert(vector<int> &block)
{
	// start from root node
	Trie* curr = this;
	for (int i = 0; i < block.size(); i++)
	{
		// create a new node if path doesn't exists
		int index = curr->findOrCreateSon(block[i]);

		// go to next node
		curr = &curr-> sons[index];
	}

	// mark current node as leaf
	if (curr->flag == false)
	{
		curr->flag = true;
		if (blockVec.size()<block.size())
		{
			for (size_t i = blockVec.size(); i < block.size(); i++)
			{
				blockVec.push_back(0);
			}
		}
		blockVec[block.size() - 1]++; //updating the number of unique blocks of the relevant size
	}
	
}

int Trie::findOrCreateSon(int sonNum) {//need to test this to see why it won't insert second block to Trie
	//binary search the vector of sons creates son if didn't exist before and return son's index 
	//Need testing
	if (sons.size() == 0)
	{
		Trie newTrie(sonNum);
		vector<Trie>::iterator it = sons.begin();
		sons.insert(it, newTrie);
		return 0;
	}
	if (sonNum > sons.back().key)
	{
		Trie newTrie(sonNum);
		sons.push_back(newTrie);
		return sons.size() - 1;
	}
	size_t left = 0;
	size_t right = sons.size() - 1;
	while (left <= right)
	{
		size_t mid = (left + right) / 2;
		if (sons[mid].key == sonNum)
		{
			return mid;
		}
		if (sons[mid].key < sonNum)
		{
			left = mid + 1;
			if (sons[left].key > sonNum)
			{
				Trie newTrie(sonNum);
				vector<Trie>::iterator it = sons.begin() + left;
				sons.insert(it, newTrie);
				return left;
			}
		}
		if (sons[mid].key > sonNum)
		{
			if (mid == 0)
			{
				Trie newTrie(sonNum);
				vector<Trie>::iterator it = sons.begin();
				sons.insert(it, newTrie);
				return 0;
			}
			if (sons[mid - 1].key < sonNum)
			{
				Trie newTrie(sonNum);
				vector<Trie>::iterator it = sons.begin() + mid;
				sons.insert(it, newTrie);
				return mid;
			}
			right = mid - 1;
		}
	}

}


Trie::~Trie()
{
}