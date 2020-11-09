
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

#define ENDTIME   0.000000001	// end of time
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

//queue <double> rand_queue;

mutex my_mutex1;
mutex my_mutex2;
mutex my_mutex3;
mutex my_mutex4;

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
        pair<int, int> _shared_x;       //pair(id_of_shared_element, amount_of_this_shared_elemet)
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
	//reactant[1].push_back(make_pair(3, 1));	 //r1(D, 1)
	product[1].push_back(make_pair(4, 1));		//r1(E, 1)
	product[1].push_back(make_pair(5, 1));		//r1(F, 1)
	//product[1].push_back(make_pair(3, 1));	//r1(D, 1)
	
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
					//cout<<"p[0]="<<p<<endl; 
					break;
				}			// A + B --r0--> C
		case 1:
				{ 
					p = rate[1]*x[2]*x[3]; 
					//cout<<"p[1]="<<p<<endl; 
					break; 
				}			// C + D --r1--> E + F + D	
		case 2:
				{ 
					p = rate[2]*x[5]*(x[5]-1)/2; 
					//cout<<"p[2]="<<p<<endl; 
					break; 
				}		// F + F --r2--> G	
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
	printDependencyGraph();
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

void printPartitions(SubNetwork * alpha, SubNetwork *beta, SubNetwork *shared)
{
	cout<<"Alpha reactions: ";
	for(int j=0; j<alpha->reactions.size(); j++)
		cout<<alpha->reactions[j]<<" ";
	cout<<endl;
	
	cout<<"Alpha species: ";
	for(int j=0; j<alpha->species.size(); j++)
		cout<<alpha->species[j]<<" ";
	cout<<endl;
	
	cout<<"Beta reactions: ";
	for(int j=0; j<beta->reactions.size(); j++)
		cout<<beta->reactions[j]<<" ";
	cout<<endl;
	
	cout<<"Beta species: ";
	for(int j=0; j<beta->species.size(); j++)
		cout<<beta->species[j]<<" ";
	cout<<endl;
	
	cout<<"Shared species: ";
	for(int j=0; j<shared->species.size(); j++)
		cout<<shared->species[j]<<" ";
	cout<<endl;
}

void partition(SubNetwork * alpha, SubNetwork *beta, SubNetwork *shared)
{
	alpha->reactions.push_back(0);
	alpha->species.push_back(0);
	alpha->species.push_back(1);
	
	
	beta->reactions.push_back(1);
	beta->reactions.push_back(2);
	beta->species.push_back(3);
	beta->species.push_back(4);
	beta->species.push_back(5);
	beta->species.push_back(6);
	
	shared->species.push_back(2);
	
	printPartitions(alpha, beta, shared);
}

void nextReactionMethods( int i_th, SubNetwork *object, SubNetwork *shared,
		 deque <NotificationQueue> &notification_for_this_thread,  deque<NotificationQueue> &notification_for_other_thread,
		 bool &this_thread_reached_t_end, bool &other_thread_reached_t_end)
{
	
	auto startReac = chrono::high_resolution_clock::now();
	cout.precision(20);

	int the_number_of_undo_reactions = 0;

 	int m = object->species.size();

    int n = object->reactions.size();
    int n_shared_species = shared->species.size();

	double my_tau[n];	//puntative time per reaction
	double my_p[n];		//propensity per reaction
    int my_x[M];		//abundances of all species
	double t = 0;
	int iter = 0;
 
    stack<StateStack> stack_of_shared_species[n_shared_species];
	
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
	for(int i=0; i<n; i++)
		my_p[i] = computePropensity(object->reactions[i], my_x);

	//computing tau
	for(int i=0; i<n; i++)
		my_tau[i] = computeTau(my_seed, my_p[i]);
		
	//generating MinHeap
	MinIndexedPQ my_queue(n);
	for(int i=0; i<n; i++)
		my_queue.insert(i, my_tau[i]);
	
	//-------------------------------------STARTING ITERATIONS--------------------------------------//
	while(t<ENDTIME || !other_thread_reached_t_end)
	{
		iter++;			

		int t_sleep = 200 + 10 * i_th;
		this_thread::sleep_for(chrono::milliseconds(t_sleep));
		
		//taking min value from our heap
		double min_tau = my_queue.minKey();
		int min_tau_index = my_queue.minIndex();
		int fired_reaction = object->reactions[min_tau_index];
				
		int idOfElement;
		int changedAmount;
		int temp_x[M];

		//updating time
		double t_old = t;
		t = min_tau;

        bool returned_to_saved_state = false;



		bool undo_reaction = false;
//---------------------------------------------------------------------------------------------
        std::unique_lock<std::mutex> locker(my_mutex1);
        if(!notification_for_this_thread.empty())
        {
            NotificationQueue temp_notification_for_this_thread = notification_for_this_thread.front();           

            if(temp_notification_for_this_thread._time < t)
            {
                //In order to find which element has been changed regarting this notification
				int my_id;
                int id_of_shared_element = temp_notification_for_this_thread._shared_x.first; //c
                for(int i=0; i<n_shared_species; i++)
				{
                    if(shared->species[i] == id_of_shared_element)
                    {
                        my_id = i;
                        break;
                    }
				}     

				if(stack_of_shared_species[my_id].empty() || temp_notification_for_this_thread._time > t_old)  
				{					
					notification_for_this_thread.pop_front();
					int id_of_shared_element = temp_notification_for_this_thread._shared_x.first;
					int amount_of_this_shared_elemet = temp_notification_for_this_thread._shared_x.second;
					my_x[id_of_shared_element] = amount_of_this_shared_elemet;    
				}         
				else
				{
					while(!stack_of_shared_species[my_id].empty())
					{
						StateStack temp_stack_of_shared_species = stack_of_shared_species[my_id].top();
						if(!notification_for_other_thread.empty() && temp_notification_for_this_thread._time > notification_for_other_thread.front()._time)
							break;
						else if(temp_stack_of_shared_species.time_when_changes_happend > temp_notification_for_this_thread._time)
						{
							undo_reaction = true;stack_of_shared_species[my_id].pop();							
							the_number_of_undo_reactions++;

							t = temp_stack_of_shared_species.safe_time;				
							
							for(int i=0; i<M; i++)
							{
								my_x[i] = temp_stack_of_shared_species.safe_value_of_x[i];
							}
							for(int i=0; i<n; i++)
							{
								my_tau[i] = temp_stack_of_shared_species.safe_value_of_tau[i];
								my_queue.changeKey(i, my_tau[i]);
							}

							// undo all changes that has been done on shared elements before time t
							for(int i=0; i<n_shared_species; i++)
							{
								while(!stack_of_shared_species[i].empty())
								{
									temp_stack_of_shared_species = stack_of_shared_species[i].top();
									if(temp_stack_of_shared_species.safe_time > t)
									{
										stack_of_shared_species[i].pop();
										the_number_of_undo_reactions++;
										t = temp_stack_of_shared_species.safe_time;				
										
										for(int i=0; i<M; i++)
										{
											my_x[i] = temp_stack_of_shared_species.safe_value_of_x[i];
										}
										for(int i=0; i<n; i++)
										{
											my_tau[i] = temp_stack_of_shared_species.safe_value_of_tau[i];
											my_queue.changeKey(i, my_tau[i]);
										}
									}
									else
									{
										break;
									}
									
								}
							}

							// delete all notes for other thread
							while(!notification_for_other_thread.empty() && notification_for_other_thread.back()._time > t)
								notification_for_other_thread.pop_back();
						}
						else
							break;
						
					}                
            	}
				
        	}           
        }
        locker.unlock();

        if(undo_reaction)
            continue;
//---------------------------------------------------------------------------------------------


        if(t<ENDTIME && this_thread_reached_t_end)
            this_thread_reached_t_end = false;
        else if(t_old>ENDTIME && this_thread_reached_t_end)
        {
            t = t_old;
			continue;
        }

		for(int i=0; i<M; i++)
			temp_x[i] = my_x[i];

        bool bool_array_of_shared_elements[n_shared_species];
        for(int i=0; i<n_shared_species; i++)
            bool_array_of_shared_elements[i] = false;

		bool change_on_shared_species = false;
		
		// updating abundance of reactants
		for(int i=0; i<(int)reactant[fired_reaction].size(); i++)
		{
            idOfElement   = reactant[fired_reaction][i].first;
			changedAmount = reactant[fired_reaction][i].second;

			for(int j=0; j<n_shared_species; j++)
			{
				if(idOfElement == shared->species[j])
				{
                    change_on_shared_species = true;
                    bool_array_of_shared_elements[j] = true;

					temp_x[idOfElement] = my_x[idOfElement] - changedAmount;

					NotificationQueue temp_notification_for_other_thread;
					temp_notification_for_other_thread._time = t;
					temp_notification_for_other_thread._shared_x = make_pair(idOfElement, temp_x[i]);
					std::unique_lock<std::mutex> locker(my_mutex3);		
					notification_for_other_thread.push_back(temp_notification_for_other_thread);
					locker.unlock();
                    break;
				}
			}			
			temp_x[idOfElement] = my_x[idOfElement] - changedAmount;
		}	

		// updating abundance of products
		for(int i=0; i<(int)product[fired_reaction].size(); i++)
		{
			idOfElement   = product[fired_reaction][i].first;
			changedAmount = product[fired_reaction][i].second;
			
			for(int j=0; j<n_shared_species; j++)
			{
				if(idOfElement == shared->species[j])
				{
					change_on_shared_species = true;
                    bool_array_of_shared_elements[j] = true;

					temp_x[i] = my_x[idOfElement] + changedAmount;

					
					NotificationQueue temp_notification_for_other_thread;
					temp_notification_for_other_thread._time = t;
					temp_notification_for_other_thread._shared_x = make_pair(idOfElement, temp_x[i]);			
					std::unique_lock<std::mutex> locker(my_mutex4);
					notification_for_other_thread.push_back(temp_notification_for_other_thread);
					locker.unlock();
                    break;
				}
			}
			temp_x[i] = my_x[idOfElement] + changedAmount;
		}

		// push the time, x, and tau for save the safe state if there have been done changes on shared elements
		if(change_on_shared_species == true)
		{
			StateStack temp_stack_of_shared_species;
			temp_stack_of_shared_species.time_when_changes_happend = t;
            temp_stack_of_shared_species.safe_time = t_old;
			for(int i=0; i<M; i++)
				temp_stack_of_shared_species.safe_value_of_x[i] = my_x[i];
			for (int i=0; i<n; i++)
				temp_stack_of_shared_species.safe_value_of_tau[i] = my_tau[i];

            for(int i=0; i<n_shared_species; i++)
            {
                if(bool_array_of_shared_elements[i])
                    stack_of_shared_species[i].push(temp_stack_of_shared_species);
            }			
		}
		

		for(int i=0; i<M; i++)
			my_x[i] = temp_x[i];
					
		
		// computing new propensity value and tau for each reaction that depends on fired reactionDependency
		for(int i=0; i<dependency_graph[fired_reaction].size(); i++)
		{
			int effected_reaction = dependency_graph[fired_reaction][i];
			//cout<<"Looking for effected reaction: "<<effected_reaction<<endl;
					
			for(int j=0; j<object->reactions.size(); j++)
			{				
				if(effected_reaction == object->reactions[j]) //we should look dependences only from this sub-network
				{
					double old_p = my_p[j];
					double new_p = computePropensity(effected_reaction, ref(my_x));
					double new_tau;
					
					if(effected_reaction != fired_reaction)
					{
						new_tau = t + (old_p/new_p)*(my_tau[j]-t); 
						//new_tau = computeTau(my_seed, new_p);
					}
					else
					{
						new_tau = t + computeTau(my_seed, new_p);
						//new_tau = computeTau(my_seed, new_p);
					}

					if(new_tau<0)
						break;

					my_tau[j] = new_tau;
					my_p[j] = new_p;
					my_queue.changeKey(j, new_tau);
				}
			}
		}

		if(t>ENDTIME
			this_thread_reached_t_end = true;
	}
	
	auto elapsedReac = chrono::high_resolution_clock::now() - startReac;

	int t_sleep = 200 + 10 * i_th;
	this_thread::sleep_for(chrono::milliseconds(t_sleep));	


	cout<<"THREAD ["<<i_th<<"]: "<<endl;
	cout<<" Time	chrono: " << chrono::duration_cast<chrono::microseconds>(elapsedReac).count() << " microseconds"<< endl; 
	cout<<" UNDO =  "<<the_number_of_undo_reactions<<endl;
	cout<<" Iter: "<<iter<<endl;
	cout<<" my_x:";
	for(int i=0; i<M; i++)
		cout<<" "<<my_x[i];
	cout<<endl<<endl;
}




int main(int argc, char * argv[])
{
	cout.precision(20);

	SubNetwork *alpha 	= new SubNetwork;
	SubNetwork *beta	= new SubNetwork;
	SubNetwork *shared	= new SubNetwork;

    deque<NotificationQueue> notify_alpha;
    deque<NotificationQueue> notify_beta;

	bool alpha_reached_t_end = false;
	bool beta_reached_t_end = false;
	
	
//Initializing starting values
	init();	
//Finding dependenceis between reactions
	reactionDependency();
//generating dependency graph
	generateDependencyGraph();
//partitioning our species into 2 sub-networks
	partition(alpha, beta, shared);


vector<thread*> myThread(2);
cout.precision(5);

//Starting our simulation
	auto startReac = chrono::high_resolution_clock::now();
		myThread[0] = new thread(nextReactionMethods, 1, ref(alpha), ref(shared), ref(notify_alpha), ref(notify_beta), 
									ref(alpha_reached_t_end), ref(beta_reached_t_end));
		
		myThread[1] = new thread(nextReactionMethods, 2, ref(beta), ref(shared), ref(notify_beta), ref(notify_alpha),
									ref(beta_reached_t_end), ref(alpha_reached_t_end));

		for(int i_th=0; i_th<2; i_th++)
			myThread[i_th]->join();

	auto elapsedReac = chrono::high_resolution_clock::now() - startReac;
	cout << "Time	chrono: " << chrono::duration_cast<chrono::microseconds>(elapsedReac).count() << " microseconds"<< endl; 
	
	return 0;
}

