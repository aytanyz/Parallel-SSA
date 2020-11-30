
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
int     num_partitions = 2;
double  ENDTIME;                // end of time

queue<SendedCurrentState> update[Nthreads];
vector<bool> flag(num_partitions);
AverageValue generated_average[Nthreads];
int count_redistribute = 0;
int have_terminated[Nthreads];

void init_abundance()
{
    x[0]  = 360000;   //E
    x[1]  = 360000;   //S
    x[2]  = 0;        //P
    x[3]  = 0;        //ES

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

double computePropensity(int id, double rate[N])
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
					p = rate[1]*x[3];  
					break; 
				}
		case 2:
				{ 
					p = rate[2]*x[3];
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

void computeBounds(double rate[N], int *upper_x, int *lower_x, double *upper_p, double *lower_p)
{
    for(int i=0; i<M; i++)
    {
        upper_x[i] = 1.1 * x[i];
        lower_x[i] = 0.9 * x[i];  
    }    

    double p;                
    for(int i=0; i<N; i++)
    {
        p = computePropensity(i, rate);
        upper_p[i] = 1.1 * p;
        lower_p[i] = 0.9 * p; 
    }
}

void rejectionSSA()
{
    auto startReac = chrono::high_resolution_clock::now();

    double rate[N];    	        // reaction rates
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

    int idOfElement;
    int changedAmount;
    int iter = 0;
    int real_iter = 0;
    int count_recalculate = 0;
    int last_x[M];
    int last_iter;

    init(rate, reactant, product);
    my_seed = generateSeed();

    computeBounds( rate, ref(upper_x), ref(lower_x), ref(upper_p), ref(lower_p));

    p_sum=0;
    for(int i=0; i<N; i++)
        p_sum += upper_p[i];

    bool recalculate_x = false;
    
    while(t<ENDTIME)
    {  

        if(recalculate_x)
        {
            count_recalculate++;
            computeBounds(rate, ref(upper_x), ref(lower_x), ref(upper_p), ref(lower_p));
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
                    int new_p = computePropensity(reac_id, rate);
                    if(r2 <= (new_p/upper_p[reac_id]))
                        accepted = true;                        
                }
                u = u * r3;              
            }  

            tau = (-1/upper_p[reac_id])*log(u);  
            t += tau;
            
            iter++;
                        
            for(int i=0; i<(int)reactant[reac_id].size(); i++)
            {
                idOfElement   = reactant[reac_id][i].first;
                changedAmount = reactant[reac_id][i].second;

                x[idOfElement] -= changedAmount;
            }
            // updating abundance of products
            for(int i=0; i<(int)product[reac_id].size(); i++)
            {
                idOfElement   = product[reac_id][i].first;
                changedAmount = product[reac_id][i].second;

                x[idOfElement] += changedAmount;
            } 
            
            for(int i=0; i<M; i++)
            {
                if(x[i]>=upper_x[i] || x[i]<=lower_x[i])
                {
                    recalculate_x = true;
                    break;
                }
            }
        }   
    }
    auto elapsedReac = chrono::high_resolution_clock::now() - startReac;
    cout<<"   iter = "<<iter<<endl;
    cout<<"   recalculate = "<<count_recalculate<<endl;
    cout<<"   last_x:"; 
    for(int i=0; i<M; i++)
        cout<<" "<<x[i];
    cout<<endl;
    cout<<"Time: " << chrono::duration_cast<chrono::milliseconds>(elapsedReac).count() << " milliseconds"<< endl;
}

int main(int argc, char * argv[])
{   
    init_abundance();
    
    //-----------------------------------------------------------------------------------
    auto startReac = chrono::high_resolution_clock::now();

    rejectionSSA();
 
    auto elapsedReac = chrono::high_resolution_clock::now() - startReac;
    cout<<endl<< "Total_time: " << chrono::duration_cast<chrono::milliseconds>(elapsedReac).count() << " milliseconds"<< endl; 

    //-----------------------------------------------------------------------------------
    return 0;
}