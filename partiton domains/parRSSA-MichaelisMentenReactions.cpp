// The code below shows the Domain Decomposion method 
// that was written based on the RSSA method. 
// This is code for experiment – Michaelis-Menten Reactions 

// g++ -O3 -pthread RSSA-MichaelisMentenReactions.cpp

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

using namespace std;

#define M 		    4 			// number of reactions 
#define N		    3 			// number of species 
#define Nthreads  	256         // max number of Threads in the machine
#define T_max    	0.5			// termination time
#define interval  	100			

class Partions {
    public:
        int x[M];       
};

int x[M];	 			        // population of chemical species		
int abundance_P[100][Nthreads];
int abundance_ES[100][Nthreads];

void init_abundance() {
    x[0]  = 360000;    		//E
    x[1]  = 360000;    		//S
    x[2]  = 0;        		//P
    x[3]  = 0;        		//ES
}

void init(double rate[N], vector<pair<int, int>> reactant[N], vector<pair<int, int>> product[N]) {	
    // reaction rates	
    rate[0] = 0.01;
    rate[1] = 250;
    rate[2] = 1;
	
    // E + S -> ES      
    reactant[0].push_back(make_pair(0, 1)); //r0(E, 1)
    reactant[0].push_back(make_pair(1, 1)); //r0(S, 1)
    product[0].push_back(make_pair(3, 1));	//r0(ES, 1)		
	
    // ES -> E + S
    reactant[1].push_back(make_pair(3, 1));	//r1(ES, 1)
    product[1].push_back(make_pair(0, 1));	//r1(E, 1)
    product[1].push_back(make_pair(1, 1));	//r1(S, 1)
	
    // ES -> E + P
    reactant[2].push_back(make_pair(3, 1)); //r2(ES, 1)
    product[2].push_back(make_pair(0, 1)); 	//r2(E, 1)
    product[2].push_back(make_pair(2, 1));	//r2(P, 1)
}

double computePropensity(int p_id, int num_partitions, int sub_domain[M], double rate[N]) {
    double p;
	
    switch(p_id){
        case 0: { 
            p = rate[0]*num_partitions*sub_domain[0]*sub_domain[1]; 
            break;
        }   
        case 1: { 
            p = rate[1]*sub_domain[3];  
            break; 
        }
        case 2: { 
            p = rate[2]*sub_domain[3];
            break; 
        } 
    }	
    return p;
}

double generateSeed(){
    auto thread_id = std::this_thread::get_id();
    uint64_t* ptr=(uint64_t*) &thread_id;
    double seed = *ptr + time(NULL);
    return seed;
}

double generateRandomNumber(double seed) {
    thread_local static std::mt19937 rng(seed);
    thread_local std::uniform_real_distribution<double> urd;
    return urd(rng, decltype(urd)::param_type{0,1});
}

void partitioning(Partions *partition, int num_partitions) {
    int more;
    int amount;

    for(int i=0; i<M; i++)   {
        more = x[i] % num_partitions;
        amount = x[i] / num_partitions;

        for(int j=0; j<num_partitions; j++) {
            if(more>0) {
                partition[j].x[i] = amount + 1;
                more--;
            }
            else
                partition[j].x[i] = amount;            
        }            
    }
}

void computeBounds(int num_partitions, int sub_domain[N], double rate[N], int *upper_x, int *lower_x, double *upper_p, double *lower_p) {
    for(int i=0; i<M; i++){
        upper_x[i] = 1.001 * sub_domain[i];
        lower_x[i] = 0.999 * sub_domain[i];  
    }   
    upper_p[0] = rate[0] * num_partitions * upper_x[0] * upper_x[1]; 
    lower_p[0] = rate[0] * num_partitions * lower_x[0] * lower_x[1]; 
	
    upper_p[1] = rate[1] * upper_x[3]; 
    lower_p[1] = rate[1] * lower_x[3];
	
    upper_p[2] = rate[2] * upper_x[3];
    lower_p[2] = rate[2] * lower_x[3];
}

void rejectionSSA(int i_th, int num_partitions, int partArr[M]) {
    auto startReac = chrono::high_resolution_clock::now();
    double rate[N]; 			// reaction rates
    int sub_domain[M];		    	        
    double upper_p[N];
    double lower_p[N];
    int upper_x[M];
    int lower_x[M];
    vector<pair<int, int>>	reactant[N];
    vector<pair<int, int>>	product[N];
    double t = 0;
    double tau;
    double my_seed;
    int reac_id;
    int elementID; 
    double upper_p_sum;
    bool recalculate_intervals = true;

    int update_value;
    double show_result = 0;
    double update_interval = (T_max*1.0)/(interval*1.0);
    int index = 0;

    init(rate, reactant, product);
    my_seed = generateSeed();

    for(int i=0; i<M; i++)
        sub_domain[i] = partArr[i];

    while(t<T_max)   {    
        if(recalculate_intervals) {
            computeBounds(num_partitions, sub_domain, rate, ref(upper_x), ref(lower_x), ref(upper_p), ref(lower_p));
            
            upper_p_sum=0;
            for(int i=0; i<N; i++)
                upper_p_sum += upper_p[i];
            
            recalculate_intervals = false;
        }

        while(t<T_max && !recalculate_intervals) {
            double u = 1;
            bool accepted = false; 

            while(!accepted) {
                double r1 = generateRandomNumber(my_seed);
                double r2 = generateRandomNumber(my_seed);
                double r3 = generateRandomNumber(my_seed);
                double sum = 0;

                for(int i=0; i<N; i++){
                    sum += upper_p[i];
                    if(sum > upper_p_sum*r1) {
                        reac_id = i;
                        break;
                    }
                }

                if(r2 <= (lower_p[reac_id]/upper_p[reac_id]))
                    accepted = true;
                else {
                    int new_p = computePropensity(reac_id, num_partitions,sub_domain, rate);
                    if(r2 <= (new_p/upper_p[reac_id]))
                        accepted = true;                        
                }
                u = u * r3;              
            } 

            if(t >= show_result) {
                abundance_P[index][i_th] = sub_domain[2];
                abundance_ES[index][i_th] = sub_domain[3];
                index++;
                show_result += update_interval;
            }
            
            tau = (-1/upper_p_sum)*log(u);  
            t += tau;   

            // updating abundance of reactants
            for(int i=0; i<(int)reactant[reac_id].size(); i++) {
                elementID   = reactant[reac_id][i].first;
                update_value = reactant[reac_id][i].second;
                sub_domain[elementID] -= update_value;
            }

            // updating abundance of products
            for(int i=0; i<(int)product[reac_id].size(); i++){
                elementID   = product[reac_id][i].first;
                update_value = product[reac_id][i].second;
                sub_domain[elementID] += update_value;
            } 
            
            for(int i=0; i<M; i++){
                if(sub_domain[i]>=upper_x[i] || sub_domain[i]<=lower_x[i]) {
                    recalculate_intervals = true;
                    break;
                }
            }
        }   
    } 
}

int main(int argc, char * argv[])
{
    if (argc<2) {
        std::cerr << "use: " << argv[0]  << " number1 as num_partitions\n";
        return -1;
    }	
    
    int num_partitions =  std::stol(argv[1]);    	
    
    Partions *partition = new Partions[num_partitions];
    init_abundance();
    partitioning(ref(partition), num_partitions);    
    
    vector<thread*> myThread(num_partitions);
    for (int i_th = 0; i_th < num_partitions; i_th++)
        myThread[i_th] = new thread(rejectionSSA, i_th,  num_partitions, partition[i_th].x);

    for(int i_th=0; i_th<num_partitions; i_th++)
        myThread[i_th]->join();	

    return 0;
}
