
#include <iostream>
#include <math.h>
#include <limits.h>
#include <vector>
#include <string>
#include <omp.h>
#include <chrono>
#include <utility> 
#include <time.h> 
#include <mutex>
#include <random>
#include <algorithm>
#include <stdlib.h>
#include <thread>
#include <queue>
#include <deque>
#include <stack>
#include <atomic>
#include <pthread.h>

#include "MinIndexedPriorityQueue.cpp"

//#include <ff/ff.hpp>
//#include <ff/pipeline.hpp>


using namespace std;
//using namespace ff;

#define ENDTIME   400			// end of time
#define M 		  7				// number of reaction 
#define N		  3				// number of species 
#define tDelay	  0.1			// delay time for reactions NCD, CD reactions

//int 	x[M];	 				// population of chemical species		
double 	rate[N];		    	// reaction rates
int		with_delay[N];			// if reaaction is CD = 1 , NCD = 2, ND = 3 

vector<int>	dependency_graph[M];
vector<int>	reaction_dependency_graph[N];
vector<int>	depends_on[N];
vector<int>	affects[N];
vector<pair<int, int>>	reactant[N];
vector<pair<int, int>>	product[N];


class SubNetwork
{
	public: 
		vector<int> reactions;
		vector<int> species;
};

class StateStack //stack
{
    public:
        double time_when_changes_happend;
        double safe_time;
        int safe_value_of_x[M];
        double safe_value_of_tau[N];
};

class NotificationQueue //queue
{
    public:
        double _time;
        pair<int, int> _shared_x; 
};

void init()
{	
	// reaction rates	
	rate[0] = 0.01;
	rate[1] = 0.01;
	rate[2] = 0.001;
	
	
	// A + B --r0--> C	
	reactant[0].push_back(make_pair(0, 1)); 	//r0(A, 1)
	reactant[0].push_back(make_pair(1, 1)); 	//r0(B, 1)
	product[0].push_back(make_pair(2, 1));	 	//r0(C, 1)		
	
	// C + D --r1--> E + F + D		
	reactant[1].push_back(make_pair(2, 1)); 	//r1(C, 1)
	reactant[1].push_back(make_pair(3, 1));	    //r1(D, 1)
	product[1].push_back(make_pair(4, 1));		//r1(E, 1)
	product[1].push_back(make_pair(5, 1));		//r1(F, 1)
	product[1].push_back(make_pair(3, 1));	    //r1(D, 1)
	
	// F + F --r2--> G		
	reactant[2].push_back(make_pair(5, 2)); 	//r2(F, 2)
	product[2].push_back(make_pair(6, 1));	 	//r2(G, 1)
}	
	
double computePropensity(int id, int x[M])
{
	double p;
	
	switch(id)
	{
		case 0:		
				{ 
					p = rate[0]*x[0]*x[1]; 
					break;
				}	
		case 1:
				{ 
					p = rate[1]*x[2]*x[3]; 
					break; 
				}	
		case 2:
				{ 
					p = rate[2]*x[5]*(x[5]-1)/2; 
					break; 
				}	
	}	
	return p;
}

void reactionDependency()
{
	// A + B --r0--> C
	depends_on[0].push_back(0);		// A
	depends_on[0].push_back(1);		// B
	affects[0].push_back(3);		// C
	
	// C + D --r1--> E + F + D	
	depends_on[1].push_back(2);		// C
	depends_on[1].push_back(3);		// D
	affects[1].push_back(4);		// E
	affects[1].push_back(5);		// F
	affects[1].push_back(3);		// D
	
	// F + F --r2--> G		
	depends_on[2].push_back(5);		// F
	affects[2].push_back(6);		// G
}

void printDependencyGraph()
{
	for(int i=0; i<N; i++)
	{
		cout<<"For reaction "<<i<<" there is affected reactions: ";
		for(int j=0; j<dependency_graph[i].size(); j++)
			cout<<dependency_graph[i][j]<<" ";
		cout<<endl;
	}		
}

void generateDependencyGraph()		// Dependency between reactions
{
	for(int vi=0; vi<N; vi++) //vi
	{
		for(int vj=0; vj<N; vj++) //vj
		{		
			bool hasSameElement = false;	
			
			if(vi==vj)
			{
				hasSameElement = true;
			}
			else
			{							
				for(int i=0; i<affects[vi].size(); i++)
				{
					for(int j=0; j<depends_on[vj].size(); j++)
					{
						if(affects[vi][i] == depends_on[vj][j])
						{
							hasSameElement = true;
							break;
						}
					}
				}
			}
			
			if(hasSameElement)
				dependency_graph[vi].push_back(vj); // R.vj depends on R.vi
		}
	}
	//printDependencyGraph();
}

double generateSeed()
{
	auto thread_id = std::this_thread::get_id();
	uint64_t* ptr=(uint64_t*) &thread_id;
	double seed = *ptr + time(NULL);
	return seed;
}

double generateRandomNumber(double seed) 
{
	thread_local static std::mt19937 rng(seed);
	thread_local std::uniform_real_distribution<double> urd;
	return urd(rng, decltype(urd)::param_type{0,1});
}

double computeTau(double seed, double p)
{	
	double rand_num = (1/p)*log(1/generateRandomNumber(seed));
	return rand_num;
}

void nextReactionMethods()
{
	cout.precision(20);
    double my_tau[N];	//puntative time per reaction
	double my_p[N];		//propensity per reaction
    int my_x[M];		//abundances of all species
	double t = 0;
	int iter = 0;
 
	
	//Finding seed for RandomNumberGenerator
	double my_seed = generateSeed();

    //Initial abundance of species
    my_x[0] = 1000000;		// A
	my_x[1] = 1000000;		// B
	my_x[2] = 100;		    // C
	my_x[3] = 10;		    // D
	my_x[4] = 100;		    // E
	my_x[5] = 1000000;		// F
	my_x[6] = 1;		    // G

	//Computing propensities	
	for(int i=0; i<N; i++)
		my_p[i] = computePropensity(i, my_x);	

	//computing tau
	for(int i=0; i<N; i++)
		my_tau[i] = computeTau(my_seed, my_p[i]);
		
	//generating MinHeap
	MinIndexedPQ my_queue(N);
	for(int i=0; i<N; i++)
		my_queue.insert(i, my_tau[i]);
	
	//-------------------------------------STARTING ITERATIONS--------------------------------------//
	while(t<ENDTIME)
	{
		iter++;			
		
		//taking min value from our heap
		double min_tau = my_queue.minKey();
		int fired_reaction = my_queue.minIndex();
	
		
		int idOfElement;
		int changedAmount;
		int temp_x[M];

		//updating time
		t = min_tau;
		
		// updating abundance of reactants
		for(int i=0; i<(int)reactant[fired_reaction].size(); i++)
		{
            idOfElement   = reactant[fired_reaction][i].first;
			changedAmount = reactant[fired_reaction][i].second;
			
			my_x[idOfElement] -= changedAmount;
		}	

		// updating abundance of products
		for(int i=0; i<(int)product[fired_reaction].size(); i++)
		{
			idOfElement   = product[fired_reaction][i].first;
			changedAmount = product[fired_reaction][i].second;
			
			my_x[idOfElement] += changedAmount;
		}					
		
		// computing new propensity value and tau for each reaction that depends on fired reactionDependency
		for(int i=0; i<dependency_graph[fired_reaction].size(); i++)
		{
			int effected_reaction = dependency_graph[fired_reaction][i];
            		
			for(int j=0; j<N; j++)
			{				
				if(effected_reaction == j) 
				{
					double old_p = my_p[j];
					double new_p = computePropensity(effected_reaction, my_x);
					double new_tau;
					
					if(effected_reaction != fired_reaction)
                        new_tau = t + (old_p/new_p)*(my_tau[j]-t); 
					else
						new_tau = t + computeTau(my_seed, new_p);

					if(new_tau<0)
						cout<<"Negative tau"<<endl;

					my_tau[j] = new_tau;
					my_p[j] = new_p;
					my_queue.changeKey(j, new_tau);
				}
			}
		}
	}

    cout<<"t: "<<t<<endl;
    cout<<"iter: "<<iter<<endl<<"x: ";
    for(int i=0; i<M; i++)
        cout<<" "<<my_x[i];
    cout<<endl;
}




int main(int argc, char * argv[])
{
//Initializing starting values
	init();	
//Finding dependenceis between reactions
	reactionDependency();
//generating dependency graph
	generateDependencyGraph();

//Starting our simulation
	auto startReac = chrono::high_resolution_clock::now();
	nextReactionMethods();
	auto elapsedReac = chrono::high_resolution_clock::now() - startReac;
	cout << "Time	chrono: " << chrono::duration_cast<chrono::microseconds>(elapsedReac).count() << " microseconds"<< endl; 
	
	return 0;
}

