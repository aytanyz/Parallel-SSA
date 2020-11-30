// The code below shows the Synchronized DD method 
// that was written based on the RSSA method. 
// This is code for experiment – Lotka Reactions 

// g++ -O3 -pthread SynchRSSA-LotkaReactions.cpp

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

using namespace std;

#define M 		    4 		// number of reactions 
#define N		    3 		// number of species 
#define Nthreads  	256     // max number of threads in the machine
#define T_max     	60

class Partions {
    public:
        int x[M];       
};

class SendedCurrentState {
    public:
        vector<int> x;
        double t;
        double update_time;
};

class AverageValue {
    public:
        vector<int> old_x;        
        double t;
        double update_time;
        vector<int> average_x;
};

int 	x[M];	 				// population of chemical species		
mutex   mutex_update_and_average;
int     abundance_Y1[600];
int     index = 1;

vector<bool> check_new_average;
queue<SendedCurrentState> update[Nthreads];
AverageValue generated_average[Nthreads];
int have_terminated[Nthreads];

void init_abundance() {
    x[0]  = 100;        	//X
    x[1]  = 1000;       	//Y1
    x[2]  = 1000;       	//Y2
    x[3]  = 0;          	//Z
}

void init(double rate[N], vector<pair<int, int>> reactant[N], vector<pair<int, int>> product[N]) {	
    // reaction rates	
    rate[0] = 10;  	//c1X
    rate[1] = 0.01;
    rate[2] = 10;
	
    // X + Y1 -> 2Y1     
    //reactant[0].push_back(make_pair(0, 1)); 	//r0(X, 1)
    reactant[0].push_back(make_pair(1, 1)); 	//r0(Y1, 1)
    product[0].push_back(make_pair(1, 2));	    //r0(Y1, 2)		
        
    // Y1 + Y2 -> 2Y1 
    reactant[1].push_back(make_pair(1, 1));	    //r1(Y1, 1)
    reactant[1].push_back(make_pair(2, 1));	    //r1(Y2, 1)
    product[1].push_back(make_pair(2, 2));	    //r1(Y2, 1)
        
    // Y2 -> Z
    reactant[2].push_back(make_pair(2, 1)); 	//r2(Y2, 1)
    //product[2].push_back(make_pair(3, 1));	//r2(Z, 1)
}

double computePropensity(int p_id, int num_partitions, int sub_domain[M], double rate[N]) {
    double p;

    switch(p_id) {
        case 0: { 
        p = rate[0]*sub_domain[1]; //c0X 
        break;
        }   
        case 1: { 
        p = rate[1]*num_partitions*sub_domain[1]*sub_domain[2];  
        break; 
        }
        case 2: { 
            p = rate[2]*sub_domain[2];
            break; 
        } 
    }	
    return p;
}

double generateSeed() {
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

    for(int i=0; i<M; i++)  {
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

void computeBounds(int num_partitions, int sub_domain[N], double rate[N], int *upper_x, int *lower_x, double *upper_p, double *lower_p)  {
    for(int i=0; i<M; i++) {
        upper_x[i] = 1.001 * sub_domain[i];
        lower_x[i] = 0.999 * sub_domain[i];  
    }    

    upper_p[0] = rate[0] * upper_x[1]; 
    lower_p[0] = rate[0] * lower_x[1]; 
	
    upper_p[1] = rate[1] * num_partitions * upper_x[1] * upper_x[2]; 
    lower_p[1] = rate[1] * num_partitions * lower_x[1] * lower_x[2];
        
    upper_p[2] = rate[2] * upper_x[2];
    lower_p[2] = rate[2] * lower_x[2];
}

void redistributeAverage(int i_th, int num_partitions) {
    vector<int> temp_sum_x(M);
    AverageValue temp_average;

    for(int i=0; i<num_partitions; i++)  {
        for(int j=0; j<M; j++)
            temp_sum_x[j] += update[i].front().x[j];
    }    

    abundance_Y1[index] = temp_sum_x[1];
    index++;

    for(int i=0; i<M; i++)
        temp_sum_x[i] = temp_sum_x[i]/num_partitions;

    for(int i=0; i<num_partitions; i++) {
        temp_average.old_x = update[i].front().x;
        temp_average.t = update[i].front().t;
        temp_average.update_time  = update[i].front().update_time;
        temp_average.average_x  = temp_sum_x;
        update[i].pop();
        generated_average[i] = temp_average;
        check_new_average[i] = true;
    }
}

void sendDomainVector(int i_th, int num_partitions, vector<int> domain, double time, double update_time)  {
    // check if we get value from all threads for generating average
    bool can_redistribute = true;
    SendedCurrentState temp;
    
    temp.x = domain;
    temp.t = time;
    temp.update_time = update_time;

    std::unique_lock<std::mutex> locker(mutex_update_and_average);      
    update[i_th].push(temp);

    for(int i=0; i<num_partitions; i++)
        if(update[i].empty()) {
            can_redistribute = false;
            break;
        }
    for(int i=0; i<num_partitions; i++)
        if(check_new_average[i]) {
            can_redistribute = false;
            break;
        }
    if(can_redistribute)
        redistributeAverage(i_th, num_partitions);
    locker.unlock();
}

void deleteSentDomains(int i_th) {  
    while(!update[i_th].empty())
        update[i_th].pop();        
}

bool allTerminated(int num_partitions) {
    bool flag = true;
    for(int i=0; i<num_partitions; i++)
        if(!have_terminated[i]) {          
            flag = false;
            break;
        } 
    return flag;
}

void rejectionSSA(int i_th, int num_partitions, int num_update, int partArr[M]) {
    double rate[N];		    // reaction rates
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

    int idOfElement;
    int changedAmount;
    double interval = (T_max*1.0)/(num_update*1.0);
    double update_time = interval;
    bool recalculate_intervals = true;
    double upper_p_sum;

    init(rate, reactant, product);
    my_seed = generateSeed();

    for(int i=0; i<M; i++)
        sub_domain[i] = partArr[i];
    
    while(true) {         
        if(recalculate_intervals) {
            computeBounds(num_partitions, sub_domain, rate, ref(upper_x), ref(lower_x), ref(upper_p), ref(lower_p));
            
            upper_p_sum=0;
            for(int i=0; i<N; i++)
                upper_p_sum += upper_p[i];  
            recalculate_intervals = false;
        }

        if(allTerminated(num_partitions))
            break;        
        
        while(!recalculate_intervals) {
            double u = 1;
            bool accepted = false; 
            while(!accepted) {            
                if(upper_p_sum <= 0) {
                    have_terminated[i_th] = true;
                    recalculate_intervals = true;
                    break;
                }
                
                double r1 = generateRandomNumber(my_seed);
                double r2 = generateRandomNumber(my_seed);
                double r3 = generateRandomNumber(my_seed);
                double sum=0;
                    
                for(int i=0; i<N; i++) {
                    sum += upper_p[i];
                    if(sum > upper_p_sum*r1) {
                        reac_id = i;
                        break; 
                    }
                }

                if(r2 <= (lower_p[reac_id]/upper_p[reac_id]))
                    accepted = true;
                else
                {
                    int new_p = computePropensity(reac_id, num_partitions,sub_domain, rate);
                    if(r2 <= (new_p/upper_p[reac_id]))
                        accepted = true;                        
                }
                u = u * r3;              
            }  

             // Gets generated average 
            if(check_new_average[i_th])
            {                
                std::unique_lock<std::mutex> locker(mutex_update_and_average); 
                int diff;
                bool flag = false;
                for(int i=0; i<M; i++)
                {
                    diff = sub_domain[i] - generated_average[i_th].old_x[i];
                    if(generated_average[i_th].average_x[i] + diff > 1.001*generated_average[i_th].old_x[i] 
                        || generated_average[i_th].average_x[i] + diff < 0.999*generated_average[i_th].old_x[i]) {
                        flag = true;
                        break;
                    }
                }

                if(flag) {
                    have_terminated[i_th] = false;                        
                    t = generated_average[i_th].t;
                    update_time = generated_average[i_th].update_time;
                    for(int i=0; i<M; i++)
                        sub_domain[i] = generated_average[i_th].average_x[i]; 
                    deleteSentDomains(i_th);               
                }
                else{
                    for(int i=0; i<M; i++) {
                            diff = sub_domain[i] - generated_average[i_th].old_x[i];
                            sub_domain[i] = generated_average[i_th].average_x[i] + diff; 
                    } 
                }                  
                check_new_average[i_th] = false;
                locker.unlock();
                break;
            }         

            tau = (-1/upper_p_sum)*log(u);  
            t += tau;
            
            if(sub_domain[1]==0 && sub_domain[2]==0)
                cout<<"ERROR"<<endl;

            if(t>T_max && !have_terminated[i_th])
                have_terminated[i_th] = true;   

            // updating abundance of reactants
            for(int i=0; i<(int)reactant[reac_id].size(); i++) {
                idOfElement   = reactant[reac_id][i].first;
                changedAmount = reactant[reac_id][i].second;

                sub_domain[idOfElement] -= changedAmount;
            }
            // updating abundance of products
            for(int i=0; i<(int)product[reac_id].size(); i++) {
                idOfElement   = product[reac_id][i].first;
                changedAmount = product[reac_id][i].second;

                sub_domain[idOfElement] += changedAmount;
            } 


            // Sends domains for redistrebuting average
            if(t > update_time && t < T_max) {    
                update_time += interval;
                vector<int> domain(M);
                for(int i=0; i<M; i++)
                    domain[i] = sub_domain[i];
                sendDomainVector(i_th, num_partitions, domain, t, update_time);                
            }         
            
            for(int i=0; i<M; i++){
                if(sub_domain[i]>=upper_x[i] || sub_domain[i]<=lower_x[i]) {
                    recalculate_intervals = true;
                    break;
                }
            }
           
            if(allTerminated(num_partitions))
                break;
        }   
    }
}

int main(int argc, char * argv[]) {   
    int num_update =  600;    	
    int num_partitions = std::stol(argv[1]);	
    
    // initilization
    for(int th=0; th<num_partitions; th++) {
        check_new_average.push_back(false);
        have_terminated[th] = false;
    }
    
    Partions *partition = new Partions[num_partitions];
    init_abundance();
    abundance_Y1[0] = x[1];
    partitioning(ref(partition), num_partitions);

    auto startReac = chrono::high_resolution_clock::now(); 
    vector<thread*> myThread(num_partitions);
    for (int i_th = 0; i_th < num_partitions; i_th++)
        myThread[i_th] = new thread(rejectionSSA, i_th,  num_partitions, num_update, partition[i_th].x);

    for(int i_th=0; i_th<num_partitions; i_th++)
        myThread[i_th]->join();		 
    auto elapsedReac = chrono::high_resolution_clock::now() - startReac;
    
    return 0;
}