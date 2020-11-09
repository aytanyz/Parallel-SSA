#include <cstdio>
#include <iostream>
#include <stdlib.h>



using namespace std;

class MinIndexedPQ 
{
    int NMAX;
	int N;
	int *heap;
	int *index; 
	double *keys;

    void swap(int i, int j) 
	{
        int temp = heap[i]; 
		heap[i] = heap[j]; 
		heap[j] = temp;
		
        index[heap[i]] = i; 
		index[heap[j]] = j;
    }

    void bubbleUp(int k)    
	{
        while(k > 1 && keys[heap[k/2]] > keys[heap[k]])   
		{
            swap(k, k/2);
            k = k/2;
        }
    }

    void bubbleDown(int k)  
	{
        int j;
        while(2*k <= N) 
		{
            j = 2*k;
            if(j < N && keys[heap[j]] > keys[heap[j+1]])
                j++;
            if(keys[heap[k]] <= keys[heap[j]])
                break;
            swap(k, j);
            k = j;
        }
    }

public:
    // Create an empty MinIndexedPQ which can contain atmost NMAX elements
    MinIndexedPQ(int NMAX)  
	{
        this->NMAX = NMAX;
        N = 0;
		
        keys = new double[NMAX + 1];
        heap = new int[NMAX + 1];
        index = new int[NMAX + 1];
        for(int i = 0; i <= NMAX; i++)
            index[i] = -1;
    }
	
    ~MinIndexedPQ() 
	{
        delete [] keys;
        delete [] heap;
        delete [] index;
    }

    // check if the PQ is empty
    bool isEmpty()  
	{
        return N == 0;
    }

    // check if i is an index on the PQ
    bool contains(int i)    
	{
        return index[i] != -1;
    }

    // return the number of elements in the PQ
    int size()  
	{
        return N;
    }

    // associate key with index i; 0 < i < NMAX
    void insert(int i, double key) 
	{
        N++;
        index[i] = N;
        heap[N] = i;
        keys[i] = key;
		//cout<<"inserted key: "<<keys[i]<<endl;
        bubbleUp(N);
    }

    // returns the index associated with the minimal key
    int minIndex()  
	{
        return heap[1];
    }

    // returns the minimal key
    double minKey()    
	{
		//cout<<"returning key: "<<keys[heap[1]]<<endl;
        return keys[heap[1]];
    }

    // delete the minimal key and return its associated index
    // Warning: Don't try to read from this index after calling this function
    int deleteMin() 
	{
        int min = heap[1];
        swap(1, N--);
        bubbleDown(1);
        index[min] = -1;
        heap[N+1] = -1;
        return min;
    }

    // returns the key associated with index i
    double keyOf(int i)    
	{
        return keys[i];
    }

    // change the key associated with index i to the specified value
    void changeKey(int i, double key)  
	{
		//cout<<"---------------------------"<<endl;
		//cout<<"we change "<<keys[i]<<" to "<<key<<endl;
        keys[i] = key;
        bubbleUp(index[i]);
        bubbleDown(index[i]);
		//cout<<"---------------------------"<<endl;
    }

    // decrease the key associated with index i to the specified value
    void decreaseKey(int i, double key)    
	{
        keys[i] = key;
        bubbleUp(index[i]);
    }

    // increase the key associated with index i to the specified value
    void increaseKey(int i, double key)    
	{
        keys[i] = key;
        bubbleDown(index[i]);
    }

    // delete the key associated with index i
    void deleteKey(int i)   
	{
        int ind = index[i];
        swap(ind, N--);
        bubbleUp(ind);
        bubbleDown(ind);
        index[i] = -1;
    }
};

/*
int main()  
{
    MinIndexedPQ pq(10);
    for(int i=0; i<10; i++)
        pq.insert(i, 100-i*0.3);
	
	cout<<pq.minIndex()<<" "<<pq.minKey()<<endl;
	printf("%d %d\n", pq.minIndex(), pq.minKey());	// (9, 91)
    pq.deleteKey(2);
    //printf("%d %d\n", pq.minKey(), pq.size());
    pq.deleteMin();
    printf("%d %d %d\n", pq.minIndex(), pq.minKey(), pq.size());

    return 0;
}
*/