// The code below shows the Partitioning Reactions method 
// that was written based on the NRM method. 


#include <iostream>
#include <math.h>
#include <limits.h>
#include <vector>
#include <string>
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
#include <pthread.h>

#include "MinIndexedPriorityQueue.cpp"

using namespace std;

#define T_MAX   0.3		// termination time of time
#define M 	    7		// number of reactions 
#define N		3		// number of species 
	
double 	rate[N];		// reaction rates
vector<int>	dependency_graph[M];
vector<int>	reaction_dependency_graph[N];
vector<int>	depends_on[N];
vector<int>	affects[N];
vector<pair<int, int>> reactant[N];
vector<pair<int, int>> product[N];
mutex my_mutex1;

class SubNetwork {
    public: 
        vector<int> reactions;
        vector<int> species;
};

class SafeStack { //stack 
    public:
        double  time_change; 		
        double time; 			
        int value_x[M];			
        double value_tau[M];		
};

class SharedElement {		 	    //NotificationQueue
    public:
        double time_when_changed;  	// at what time element has been changed
	int ID;			                // id of the element
	int value;			            // the amount that been changed
};

void init() {	
	// reaction rates	
	rate[0] = 0.01;
	rate[1] = 0.01;
	rate[2] = 0.02;
	
	
	// A + B --r0--> C	
	reactant[0].push_back(make_pair(0, 1)); 	//r0(A, 1)
	reactant[0].push_back(make_pair(1, 1)); 	//r0(B, 1)
	product[0].push_back(make_pair(2, 1));	    //r0(C, 1)		
	
	// C + D --r1--> E + F + D		
	reactant[1].push_back(make_pair(2, 1)); 	//r1(C, 1)
	//reactant[1].push_back(make_pair(3, 1));	//r1(D, 1)
	product[1].push_back(make_pair(4, 1));	    //r1(E, 1)
	product[1].push_back(make_pair(5, 1));	    //r1(F, 1)
	//product[1].push_back(make_pair(3, 1));	//r1(D, 1)
	
	// F + F --r2--> G		
	reactant[2].push_back(make_pair(5, 2)); 	//r2(F, 2)
	product[2].push_back(make_pair(6, 1));  	//r2(G, 1)
}	
	
double computePropensity(int p_id, int x[M]) {
	double p;
	
	switch(p_id) {
		case 0:	{ 
			p = rate[0]*x[0]*x[1]; 
			break;
		}	
		case 1: { 
			p = rate[1]*x[2]*x[3]; 
			break; 
		}	
		case 2: { 
			p = rate[2]*x[5]*(x[5]-1)/2; 
			break; 
		}	
	}	
	return p;
}

void reactionDependency() {
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

void printDependencyGraph() {
	for(int i=0; i<N; i++) {
		cout<<"For reaction "<<i<<" there is affected reactions: ";
		for(int j=0; j<dependency_graph[i].size(); j++)
			cout<<dependency_graph[i][j]<<" ";
		cout<<endl;
	}		
}

void generateDependencyGraph() {		// Dependency between reactions
    for(int vi=0; vi<N; vi++) {
        for(int vj=0; vj<N; vj++) {		
	        bool hasSameElement = false;	
			
	        if(vi==vj)
	             hasSameElement = true;
	        else {							
	            for(int i=0; i<affects[vi].size(); i++) {
	                for(int j=0; j<depends_on[vj].size(); j++) {
		                if(affects[vi][i] == depends_on[vj][j])  {
		                    hasSameElement = true;
		                    break;
		                }
	                }
	            }
	        }
			
            if(hasSameElement)
                dependency_graph[vi].push_back(vj); // Rvj depends on Rvi
        }
    }
       //printDependencyGraph();
}

double generateSeed()  {
    auto thread_id = std::this_thread::get_id();
    uint64_t* ptr=(uint64_t*) &thread_id;
    double seed = *ptr + time(NULL);
    return seed;
}

double generateRandomNumber(double seed)  {
    thread_local static std::mt19937 rng(seed);
    thread_local std::uniform_real_distribution<double> urd;
    return urd(rng, decltype(urd)::param_type{0,1});
}

double computeTau(double seed, double p) {	
    double rand_num = (1/p)*log(1/generateRandomNumber(seed));
    return rand_num;
}

void partition(SubNetwork * alpha, SubNetwork *beta, SubNetwork *shared) {
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
    // printPartitions(alpha, beta, shared);
}

void nextReactionMethods( int i_th, SubNetwork *object, SubNetwork *shared, deque <SharedElement> &notification_for_this_thread,  deque<SharedElement> &notification_for_other_thread, bool &this_thread_reached_t_end, bool &other_thread_reached_t_end) {	
    auto startReac = chrono::high_resolution_clock::now();

    int count_undo_reactions = 0; 
    int m = object->species.size();
    int n = object->reactions.size();
    int m_shared = shared->species.size(); 

    double tau[n];			//puntative time per reaction
    double p[n];			//propensity per reaction
    int x[M];				//abundances of all species	
    int temp_x[M];
    double t = 0;
    double t_old;
    int elementID; 
    int used_value; 
    bool undo_reaction;
 
    stack<SafeStack> shared_species_stack[m_shared];
    double my_seed = generateSeed();

    //Initial abundance of species
    x[0] = 1000000;		// A
    x[1] = 1000000;		// B
    x[2] = 1000000;	    // C
    x[3] = 1000000;	    // D
    x[4] = 100000;		// E
    x[5] = 1000050;		// F
    x[6] = 100;		    // G

    //Computing propensities	
    for(int i=0; i<n; i++)
        p[i] = computePropensity(object->reactions[i], x);
	
    //computing tau
    for(int i=0; i<n; i++)
        tau[i] = computeTau(my_seed, p[i]);
		
    //generating MinHeap
    MinIndexedPQ minHeap(n);
    for(int i=0; i<n; i++)
        minHeap.insert(i, tau[i]);
	
    while(t<T_MAX || !other_thread_reached_t_end){
        double min_tau = minHeap.minKey();
        int min_tauID = minHeap.minIndex();
        int fired_reactionID = object->reactions[min_tauID];

        //updating time
        t_old = t;
        t = min_tau;
        undo_reaction = false;

        std::unique_lock<std::mutex> locker(my_mutex1);	
        if(!notification_for_this_thread.empty()) {
            SharedElement temp_notification_for_this_thread = notification_for_this_thread.front();        
            if(temp_notification_for_this_thread.time_when_changed < t)  {
                //In order to find which element has been changed regarding this notification
                int shared_elementID;
                for(int i=0; i<m_shared; i++) {
                    if(shared->species[i] == temp_notification_for_this_thread.ID)  {
                        shared_elementID = i;
                        break;
                    }
	            }     
                if(shared_species_stack[shared_elementID].empty()) {				
                    notification_for_this_thread.pop_front();
                    shared_elementID = temp_notification_for_this_thread.ID;
                    int amount_of_this_shared_elemet = temp_notification_for_this_thread.value;
                    x[shared_elementID] = amount_of_this_shared_elemet;    
                } 
                else if(temp_notification_for_this_thread.time_when_changed <= t && temp_notification_for_this_thread.time_when_changed >= shared_species_stack[shared_elementID].top().time_change)  {
                    notification_for_this_thread.pop_front();
                    int shared_elementID = temp_notification_for_this_thread.ID;
                    int amount_of_this_shared_elemet = temp_notification_for_this_thread.value;
                    x[shared_elementID] = amount_of_this_shared_elemet;   
                }
                else {
                    while(!shared_species_stack[shared_elementID].empty()) {
                        SafeStack temp_shared_species_stack = shared_species_stack[shared_elementID].top();
                        if(!notification_for_other_thread.empty() && temp_notification_for_this_thread.time_when_changed > notification_for_other_thread.front().time_when_changed) {
                            // other thread should do the undo
                            break;
                        }
                        else if(temp_shared_species_stack.time_change > temp_notification_for_this_thread.time_when_changed) {
                            undo_reaction = true;
                            //other thread did change before this thread. we should undo,return to the safe state		
                            shared_species_stack[shared_elementID].pop();			
                            t = temp_shared_species_stack.time;				
                            for(int i=0; i<M; i++)
                                x[i] = temp_shared_species_stack.value_x[i];
                            for(int i=0; i<n; i++) {
                                tau[i] = temp_shared_species_stack.value_tau[i];
                                minHeap.changeKey(i, tau[i]);
                            }
                            while(!notification_for_other_thread.empty() && notification_for_other_thread.back().time_when_changed > t) {
                                notification_for_other_thread.pop_back();
                            }
                        }
                        else {
                            // no need to undo any reaction, we will apply this change maybe in the future
                            break;
                        }				
                    }                
                }			
            }           
        }
        locker.unlock();

        if(undo_reaction)
            continue;
        if(t<T_MAX && this_thread_reached_t_end) 
            this_thread_reached_t_end = false;
        else if(t_old>T_MAX && this_thread_reached_t_end) {
            t = t_old;
            continue;
        }
        for(int i=0; i<M; i++)
            temp_x[i] = x[i];

        bool bool_array_of_shared_elements[m_shared];
        for(int i=0; i<m_shared; i++)
            bool_array_of_shared_elements[i] = false;

        bool change_on_shared = false;

            // updating abundance of reactants
        for(int i=0; i<(int)reactant[fired_reactionID].size(); i++) {
            elementID   = reactant[fired_reactionID][i].first;
	        used_value = reactant[fired_reactionID][i].second;
            for(int j=0; j<m_shared; j++) {
	            if(elementID == shared->species[j]) {
                    change_on_shared = true;
                    bool_array_of_shared_elements[j] = true;
                            
                    temp_x[elementID] = x[elementID] - used_value;
                    SharedElement temp_notification_for_other_thread;
                    temp_notification_for_other_thread.time_when_changed = t;
                    temp_notification_for_other_thread.ID = elementID; 
                    temp_notification_for_other_thread.value = temp_x[elementID];
                    
                    std::unique_lock<std::mutex> locker(my_mutex1);
                    notification_for_other_thread.push_back(temp_notification_for_other_thread);
		            locker.unlock();
                    break;
	            }
	        }			
	        temp_x[elementID] = x[elementID] - used_value;   
        }	

            // updating abundance of products
        for(int i=0; i<(int)product[fired_reactionID].size(); i++) {
            elementID   = product[fired_reactionID][i].first;
            used_value = product[fired_reactionID][i].second;
	    
            for(int j=0; j<m_shared; j++) {
	            if(elementID == shared->species[j]) {
                    change_on_shared = true;
                    bool_array_of_shared_elements[j] = true;

		            temp_x[elementID] = x[elementID] + used_value;		 
                    SharedElement temp_notification_for_other_thread;
		            temp_notification_for_other_thread.time_when_changed = t;
		            temp_notification_for_other_thread.ID = elementID;			
		            temp_notification_for_other_thread.value = temp_x[elementID];
		            std::unique_lock<std::mutex> locker(my_mutex1);
                    notification_for_other_thread.push_back(temp_notification_for_other_thread);
		            locker.unlock();
                    break;
	            }
	        }
	        temp_x[elementID] = x[elementID] + used_value;
        }
        // push the time, x, and tau for save the safe state if there have been done changes on shared elements
        if(change_on_shared == true) {
            SafeStack temp_shared_species_stack;
            temp_shared_species_stack.time_change = t; 
            temp_shared_species_stack.time = t_old;
	    
            for(int i=0; i<M; i++)
                temp_shared_species_stack.value_x[i] = x[i];
            for (int i=0; i<n; i++)
                temp_shared_species_stack.value_tau[i] = tau[i];
            for(int i=0; i<m_shared; i++) {
                if(bool_array_of_shared_elements[i])
                    shared_species_stack[i].push(temp_shared_species_stack);
            }			
        }
		
        for(int i=0; i<M; i++)
            x[i] = temp_x[i];

        // computing new propensity value and tau for each reaction that depends on fired     reactionDependency
        for(int i=0; i<dependency_graph[fired_reactionID].size(); i++) {
	        int effected_reaction = dependency_graph[fired_reactionID][i];
	        for(int j=0; j<object->reactions.size(); j++) {				
	            if(effected_reaction == object->reactions[j])  {
                    double old_p = p[j];
                    double new_p = computePropensity(effected_reaction, ref(x));
                    double new_tau;
					
	                if(effected_reaction != fired_reactionID)
		                new_tau = t + (old_p/new_p)*(tau[j]-t); 
	                else
		                new_tau = t + computeTau(my_seed, new_p);
                      
                    if(new_tau<0)
		                break;

                    tau[j] = new_tau;
                    p[j] = new_p;
                    minHeap.changeKey(j, new_tau);
	            }
	        }
        }

        if(t>T_MAX)
	        this_thread_reached_t_end = true;
    }
    auto elapsedReac = chrono::high_resolution_clock::now() - startReac;
}

int main(int argc, char * argv[]) {

    SubNetwork *alpha 	= new SubNetwork;
    SubNetwork *beta	= new SubNetwork;
    SubNetwork *shared	= new SubNetwork;

    deque<SharedElement> alpha_notifications; 
    deque<SharedElement> beta_notifications;

    bool alpha__is_end = false;         
    bool beta__is_end = false;
    
    init();	
    reactionDependency();
    generateDependencyGraph();
    partition(alpha, beta, shared);
    
    vector<thread*> myThread(2);
    auto startReac = chrono::high_resolution_clock::now();
    myThread[0] = new thread(nextReactionMethods, 1, ref(alpha), ref(shared), ref(alpha_notifications), ref(beta_notifications), ref(alpha__is_end), ref(beta__is_end));
    myThread[1] = new thread(nextReactionMethods, 2, ref(beta), ref(shared), ref(beta_notifications), ref(alpha_notifications), ref(beta__is_end), ref(alpha__is_end));
    
    for(int i_th=0; i_th<2; i_th++)
        myThread[i_th]->join();
    auto elapsedReac = chrono::high_resolution_clock::now() - startReac;

    return 0;
}