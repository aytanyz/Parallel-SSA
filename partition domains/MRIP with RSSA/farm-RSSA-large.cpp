// g++ -std=c++17 -Wall -g -I ~/fastflow -pthread -DNDEBUG -DTRACE_FASTFLOW -finline-functions -O3 -fopenmp farm.cpp -o ff

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
#include <time.h>
#include <algorithm>
#include <stdlib.h>
#include <thread>
#include <atomic>
#include <pthread.h>

#include <ff/ff.hpp>
#include <ff/pipeline.hpp>


using namespace std;
using namespace ff;

#define M 		  4 			// number of reaction 
#define N		  3 			// number of species 
#define Nthreads  256           // number of Threads in the machine
#define tDelay	  10			// delay time for reactions NCD, CD reactions
#define ENDTIME   0.5

class SaveResult
    {
        public:
            int x[M];
            double time;
            int iter;
    };

// First Stage. Emitter of the Farm
struct Emitter: ff_node_t< double >
{	
	int num_worker;
	
	Emitter(int n_worker) 
	{
		this->num_worker = n_worker;
	}
	
	
	int svc_init()
	{
		// in order to give each thread the different seed
		auto thread_id = std::this_thread::get_id();
		uint64_t* ptr=(uint64_t*) &thread_id;
		double seed = *ptr + time(NULL);		
		srand(seed);
		return 0;
	}
	
	
	double* svc(double *) 
	{		
		for(int i=0; i<num_worker; i++)
		{
			double rand_seed = (double)rand();
			ff_send_out(new double(rand_seed));
		}
			
		return EOS;
	};
};

// Second Stage. Worker of the Farm
struct Simulation: ff_node_t< double , vector<SaveResult> > 
{
    int 	num_worker;
	int		num_simulation;
	int 	x[M];	 				// population of chemical species		
	int		with_delay[N];			// if reaaction is CD = 1 , NCD = 2, ND = 3 
    
    vector<SaveResult> result;      //result[0] for min, result[1] for max
	
	Simulation(int n_worker, int n_simulation)
	{
		this->num_worker = n_worker;
		this->num_simulation = n_simulation;
	}

    int svc_init()
	{

        SaveResult temp;
        result.push_back(temp);        
        result.push_back(temp);

		for(int i=0; i<M; i++)
        {
            result[0].x[i] = 100000000;
            result[1].x[i] = -1;
        }
        result[0].time = 100000000;
        result[1].time = -1;
        result[0].iter = 100000000;
        result[1].iter = -1;

		return 0;
	}	

	vector<SaveResult>* svc(double *seed)
	{
		double 	&my_seed = *seed;
		for(int sim=0; sim<num_simulation/num_worker; sim++)
		{
			init_abundance();
			rejectionSSA();
		}
		return (new vector<SaveResult>(result));
	};

	void init_abundance()
	{
		x[0]  = 360000;   //E
		x[1]  = 360000;   //S
		x[2]  = 0;        //P
		x[3]  = 0;        //ES
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
				
				cout.precision(10);

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
		
        for(int i=0; i<M; i++)
        {
            if(x[i]<result[0].x[i])
                result[0].x[i] = x[i];
            if(x[i]>result[1].x[i])
                result[1].x[i] = x[i];
        }
        if(iter<result[0].iter)
            result[0].iter = iter;        
        if(iter>result[1].iter)
            result[1].iter = iter;

        double last_time = chrono::duration_cast<chrono::milliseconds>(elapsedReac).count();
        if(last_time<result[0].time)
             result[0].time = last_time;
        if(last_time>result[1].time)
             result[1].time = last_time;
    };

};

// Third Stage. Collector of the Farm
struct Collector: ff_node_t< vector<SaveResult>, int > 
{		
	int num_simulation;
	vector<SaveResult> result;
	
	Collector(int num_simulation):num_simulation(num_simulation){}

    int svc_init()
	{
        SaveResult temp;
        result.push_back(temp);
        result.push_back(temp);

		for(int i=0; i<M; i++)
        {
            result[0].x[i] = 100000000;
            result[1].x[i] = -1;
        }
        result[0].time = 100000000;
        result[1].time = -1;
        result[0].iter = 100000000;
        result[1].iter = -1;
        
		return 0;
	}			
	
	int* svc(vector<SaveResult> *v) 
	{				
		vector<SaveResult> &thread_result = *v;
        
		 for(int i=0; i<M; i++)
        {
            if(thread_result[0].x[i]<result[0].x[i])
                result[0].x[i] = thread_result[0].x[i];
            if(thread_result[1].x[i]>result[1].x[i])
                result[1].x[i] = thread_result[1].x[i];
        }
        if(thread_result[0].iter<result[0].iter)
            result[0].iter = thread_result[0].iter;        
        if(thread_result[1].iter>result[1].iter)
            result[1].iter = thread_result[1].iter;

        if(thread_result[0].time<result[0].time)
             result[0].time = thread_result[0].time;
        if(thread_result[1].time>result[1].time)
             result[1].time = thread_result[1].time;
		
		return GO_ON;		
    };
	
	void svc_end()
	{	
        cout<<"MIN x:";
		for(int i=0; i<M; i++)
            cout<<" "<<result[0].x[i];
        cout<<endl;
        cout<<"MIN time: "<<result[0].time<<" milliseconds"<<endl<<"MIN iter: "<<result[0].iter<<endl<<endl;

        cout<<"MAX x:";
		for(int i=0; i<M; i++)
            cout<<" "<<result[1].x[i];
        cout<<endl;
        cout<<"MAX time: "<<result[1].time<<" milliseconds"<<endl<<"MAX iter: "<<result[1].iter<<endl<<endl;
	};	

};

int main(int argc, char * argv[])
{
	 if (argc<3) 
	 {
        std::cerr << "use: " << argv[0]  << " number1 number2  as num_worker  num_simulation\n";
        return -1;
    }
	
	int num_worker 		= std::stol(argv[1]);	// The number of workers
	int num_simulation	= std::stol(argv[2]);	// The number of simulations
	
//************** Creating nodes of FastFlow *******************
	
	Emitter emitter(num_worker);	
	Collector collect(num_simulation);
	
	std::vector<std::unique_ptr<ff_node> > W; 
	for(int i=0;i<num_worker;++i) 
		W.push_back(make_unique<Simulation>(num_worker, num_simulation));
											
	ff_Farm<int> farm(std::move(W), emitter, collect); 				
	ff_Pipe<> pipe(farm);
	
//*************************************************************
	
	
//************************** Fast Flow ************************	
	auto startReac = chrono::high_resolution_clock::now();
	ffTime(START_TIME);
	
	if (pipe.run_and_wait_end()<0) 
	{
		error("running pipe");
		return -1;
	}	
	
	auto elapsedReac = chrono::high_resolution_clock::now() - startReac;
	ffTime(STOP_TIME);
	
	
	cout<<"Time	ffTime(GET_TIME): "<< ffTime(GET_TIME) <<endl;		
	cout << "Time	chrono: " << chrono::duration_cast<chrono::microseconds>(elapsedReac).count() << " microseconds"<< endl; 
	//std::cout << "Time	pipe.ffTime(): " << farm.ffTime() << "\n";	
//*****************************************************************
	
		
//pipe.ffStats(std::cout);
	return 0;
}

