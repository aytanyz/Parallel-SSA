//g++ -O3 -pthread  -o parlarge parlarge.cpp

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

#define M 		  4 			// number of reaction 
#define N		  3 			// number of species 
#define Nthreads  256           // number of Threads in the machine
#define tDelay	  10			// delay time for reactions NCD, CD reactions
//#define num_partitions 4

class Partions
{
    public:
        int x[M];       
};

class SendedCurrentState
{
    public:
        vector<int> x;
        double t;
        double update_time;
        int num_iter;
};

class AverageValue
{
    public:
        vector<int> old_x;        
        double t;
        double update_time;
        vector<int> average_x;
        int num_iter;
};

int 	x[M];	 				// population of chemical species		
int		with_delay[N];			// if reaaction is CD = 1 , NCD = 2, ND = 3 
double  ENDTIME;                // end of time

mutex   mutex_overall;
mutex   mutex_update_and_average;
mutex   my_mutex1;

int     overall_x[M];
double  overall_time;
int     overall_real_iter;
int     overall_iter;
double  overall_t_bound;

queue<SendedCurrentState> update[Nthreads];
vector<bool> check_for_new_average;
AverageValue generated_average[Nthreads];
int count_redistribute = 0;
int have_terminated[Nthreads];

void init_abundance()
{
    x[0]  = 360000;     //E
    x[1]  = 360000;     //S
    x[2]  = 0;          //P
    x[3]  = 0;          //ES

    ENDTIME = 1;
}

void init(double rate[N], vector<pair<int, int>> reactant[N], vector<pair<int, int>> product[N])
{	
 	// reaction rates	
	rate[0] = 0.01;
    rate[1] = 250;
    rate[2] = 1;
	
    reactant[0].push_back(make_pair(0, 1)); 	//r0(E, 1)
    reactant[0].push_back(make_pair(1, 1)); 	//r0(S, 1)
	product[0].push_back(make_pair(3, 1));	 	//r0(ES, 1)		
		
	reactant[1].push_back(make_pair(3, 1));	    //r1(ES, 1)
	product[1].push_back(make_pair(0, 1));		//r1(E, 1)
	product[1].push_back(make_pair(1, 1));		//r1(S, 1)
	
	reactant[2].push_back(make_pair(3, 1)); 	//r2(ES, 1)
	product[2].push_back(make_pair(0, 1)); 	    //r2(E, 1)
	product[2].push_back(make_pair(2, 1));	 	//r2(P, 1)
}

double computePropensity(int id, int num_partitions, int my_x[M], double rate[N])
{
	double p;

	switch(id)
	{
		case 0:		
				{ 
					p = rate[0]*num_partitions*my_x[0]*my_x[1]; 
					break;
				}   
		case 1:
				{ 
					p = rate[1]*my_x[3];  
					break; 
				}
		case 2:
				{ 
					p = rate[2]*my_x[3];
					break; 
				} 
	}	
	return p;
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

void partitioning(Partions *partition, int num_partitions)
{
    int more;
    int amount;

    for(int i=0; i<M; i++)
    {
        more = x[i] % num_partitions;
        amount = x[i] / num_partitions;

        for(int j=0; j<num_partitions; j++)
        {
            if(more>0)
            {
                partition[j].x[i] = amount + 1;
                more--;
            }
            else
                partition[j].x[i] = amount;            
        }            
    }
}

void computeBounds(int num_partitions, int my_x[N], double rate[N], int *upper_x, int *lower_x, double *upper_p, double *lower_p, double &t_bound)
{
    auto startReac = chrono::high_resolution_clock::now();
    for(int i=0; i<M; i++)
    {
        upper_x[i] = 1.1 * my_x[i];
        lower_x[i] = 0.9 * my_x[i];  
    }    

    double p;                
    for(int i=0; i<N; i++)
    {
        p = computePropensity(i, num_partitions, my_x, rate);
        upper_p[i] = 1.1 * p;
        lower_p[i] = 0.9 * p; 
    }
    auto elapsedReac = chrono::high_resolution_clock::now() - startReac;
    t_bound += chrono::duration_cast<chrono::nanoseconds>(elapsedReac).count();
}

void rejectionSSA(int i_th, int num_partitions, int num_update, int sub_domain[M])
{
    auto startReac = chrono::high_resolution_clock::now();

    double rate[N];
    int my_x[M];		    	        // reaction rates
    double my_p[N];
    double upper_p[N];
    double lower_p[N];
    int upper_x[M];
    int lower_x[M];
    vector<pair<int, int>>	reactant[N];
    vector<pair<int, int>>	product[N];

    double t=0;
    int delta = 0.1;    // 10%
    double p_sum;
    double tau;
    double my_seed;
    double rand2;
    double sum;
    int reac_id;
    double t_bound = 0;

    int idOfElement;
    int changedAmount;
    int iter = 0;
    int real_iter = 0;
    int count_recalculate = 0;
    double interval = ENDTIME/num_update;
    double update_time = interval;
    int last_x[M];
    int last_iter;

    init(rate, reactant, product);
    my_seed = generateSeed();

    for(int i=0; i<M; i++)
        my_x[i] = sub_domain[i];


    computeBounds(num_partitions, my_x, rate, ref(upper_x), ref(lower_x), ref(upper_p), ref(lower_p), ref(t_bound));

    p_sum=0;
    for(int i=0; i<N; i++)
        p_sum += upper_p[i];

    bool recalculate_x = false;

    int which[N];
    for(int i=0; i<N; i++)
        which[i]=0;
    
    while(t<ENDTIME)
    {    
        if(recalculate_x)
        {
            count_recalculate++;
            computeBounds(num_partitions, my_x, rate, ref(upper_x), ref(lower_x), ref(upper_p), ref(lower_p), ref(t_bound));
            p_sum=0;
            for(int i=0; i<N; i++)
                p_sum += upper_p[i];
        }
          
        recalculate_x = false;

        while(t<ENDTIME && !recalculate_x)
        {
            double u = 1;
            bool accepted = false; 
            while(!accepted)
            {
                double r1 = generateRandomNumber(my_seed);
                double r2 = generateRandomNumber(my_seed);
                double r3 = generateRandomNumber(my_seed);
                double sum=0;

                for(int i=0; i<N; i++)
                {
                    sum += upper_p[i];
                    cout.precision(10);
                    if(sum > p_sum*r1)
                    {
                        reac_id = i;
                        break;
                    }
                }

                if(r2 <= (lower_p[reac_id]/upper_p[reac_id]))
                    accepted = true;
                else
                {
                    int new_p = computePropensity(reac_id, num_partitions,my_x, rate);
                    if(r2 <= (new_p/upper_p[reac_id]))
                        accepted = true;                        
                }
                u = u * r3;              
            }  

            tau = (-1/upper_p[reac_id])*log(u);  
            t += tau;            

            iter++;
            which[reac_id]++;
            // updating abundance of reactants
            for(int i=0; i<(int)reactant[reac_id].size(); i++)
            {
                idOfElement   = reactant[reac_id][i].first;
                changedAmount = reactant[reac_id][i].second;

                my_x[idOfElement] -= changedAmount;
            }
            // updating abundance of products
            for(int i=0; i<(int)product[reac_id].size(); i++)
            {
                idOfElement   = product[reac_id][i].first;
                changedAmount = product[reac_id][i].second;

                my_x[idOfElement] += changedAmount;
            } 
            
            for(int i=0; i<M; i++)
            {
                if(my_x[i]>=upper_x[i] || my_x[i]<=lower_x[i])
                {
                    recalculate_x = true;
                    break;
                }
            }
        }   
    }
    auto elapsedReac = chrono::high_resolution_clock::now() - startReac;
    std::unique_lock<std::mutex> locker(mutex_overall);        
    for(int i=0; i<M; i++)
        overall_x[i] += my_x[i];
    double temp_time = chrono::duration_cast<chrono::milliseconds>(elapsedReac).count();
    overall_time = (overall_time>temp_time) ? overall_time : temp_time; 
    overall_iter += iter;
    overall_t_bound += t_bound;
    locker.unlock(); 
}

int main(int argc, char * argv[])
{
     if (argc<3) 
	 {
        std::cerr << "use: " << argv[0]  << " number1 as num_update, number2 as num_partitions\n";
        return -1;
    }	
    
    int num_update =  std::stol(argv[1]);    	
	int num_partitions = std::stol(argv[2]);	// The number of partitions
    
    // initilization
    for(int th=0; th<num_partitions; th++)
    {
        check_for_new_average.push_back(false);
        have_terminated[th] = false;
    }

    for(int i=0; i<M; i++)
        overall_x[i] = 0;
    overall_time = 0;
    overall_t_bound =0;
    overall_real_iter = 0;
    overall_iter = 0;
    
    Partions *partition = new Partions[num_partitions];
    init_abundance();
    partitioning(ref(partition), num_partitions);
    //-----------------------------------------------------------------------------------
    auto startReac = chrono::high_resolution_clock::now();
    
    vector<thread*> myThread(num_partitions);
    for (int i_th = 0; i_th < num_partitions; i_th++)
        myThread[i_th] = new thread(rejectionSSA, i_th,  num_partitions, num_update, partition[i_th].x);

    for(int i_th=0; i_th<num_partitions; i_th++)
        myThread[i_th]->join();
		
 
    auto elapsedReac = chrono::high_resolution_clock::now() - startReac;
    //-----------------------------------------------------------------------------------

    cout<<endl<<"count_redistribute = "<<count_redistribute<<endl;
    cout<<"Total_time: " << chrono::duration_cast<chrono::milliseconds>(elapsedReac).count() << " milliseconds"<< endl<<endl; 
    cout<<"overall_iter: "<<overall_iter<<endl;
    cout<<"overall_real_iter: "<<overall_real_iter<<endl;
    cout<<"Overall_x: ";
    for(int i=0; i<M; i++)
        cout<<" "<<overall_x[i];
    cout<<endl<<"Overall_time: "<<overall_time << " milliseconds"<< endl<<endl; 
    cout<<endl<<"Overall_t_bound: "<<overall_t_bound << " nanoseconds"<< endl<<endl; 
    return 0;
}