#include <iostream>
#include <vector>
#include <fstream>
#include <numeric>
#include <iomanip>
#include "include/Choose_Ske.h"
#include "include/utils.h"

int Heavy_Thes=500;
int BILI=4;

using namespace std;
vector<string> items;
unordered_map<string, int>freq1,freq0;
long double Join_Ground_Truth;

int memaccess_i = 0, memaccess_q = 0, memaccess_d = 0;

int thres_mid = 16, thres_large = 0;

int item_num = 0, flow_num = 0;
vector<double> experiment_MEM,experiment_MAX,experiment_MIN,experiment_ARE, experiment_AAE,experiment_CON,experiment_YangBen,experiment_INT,experiment_JOT;
vector<vector<double> >experiment_all;
int data_delimiter;
void readFile_CAIDA(const char* filename, int length=13, int MAX_ITEM = INT32_MAX){
	ifstream inFile(filename, ios::binary);

	if (!inFile.is_open())
		cout << "File fail." << endl;

	int max_freq = 0;
	char key[length];
	for (int i = 0; i < MAX_ITEM; ++i)
	{
		inFile.read(key, length);
		if (inFile.gcount() < length)
			break;
		items.push_back(string(key, KEY_LEN));
	}
	data_delimiter=items.size()/2;
	freq0.clear();
	freq1.clear();
    unsigned int maxvalue = 0, minvalue = UINT32_MAX;
	for (int i=0; i<items.size();i++){
		if(i<data_delimiter){
    		freq1[items[i]]++;
        }
		else {
            freq0[items[i]]++;
        }
		maxvalue = max(maxvalue, ((unsigned*)items[i].c_str())[0]);
		minvalue = min(minvalue, ((unsigned*)items[i].c_str())[0]);
	}
	inFile.close();
	for (auto pr : freq0)
		max_freq = max(max_freq, pr.second);
	for (auto pr : freq1)
		max_freq = max(max_freq, pr.second);
	
	for (auto pr : freq0){
		if(freq1.count(pr.first)){
			if(pr.second>Heavy_Thes)thres_large++;
			Join_Ground_Truth+=1ll*pr.second*freq1[pr.first];
		}
	}


	item_num = items.size();
	flow_num = freq0.size()+freq1.size();
    cout << "dataset name: " << filename << endl;
	cout << flow_num << " flows, " << items.size() << " items read" << endl;;
	cout << "max freq = " << max_freq << endl;
	cout << "large flow sum = " << thres_large  << endl;
    cout << "maxDataValue = " << maxvalue << "   minDataValue = "<< minvalue << endl;
	cout <<"Join Ground Truth = "<<Join_Ground_Truth<< endl<< endl;
}


void readFile_zipf(const char* filename,int length=4, int MAX_ITEM = INT32_MAX){
	ifstream inFile(filename, ios::binary);

	if (!inFile.is_open())
		cout << "zipf File fail." << endl;

	int max_freq = 0;
	char key[length];
	for (int i = 0; i < MAX_ITEM; ++i)
	{
		inFile.read(key, length);
		if (inFile.gcount() < length)
			break;
		items.push_back(string(key, KEY_LEN));
	}
	data_delimiter=items.size()/2;
	for (int i=0; i<items.size();i++){
		if(i<data_delimiter)
		freq1[items[i]]++;
		else 
		freq0[items[i]]++;
	}
	inFile.close();


	for (auto pr : freq0)
		max_freq = max(max_freq, pr.second);
	for (auto pr : freq1)
		max_freq = max(max_freq, pr.second);
	
	for (auto pr : freq0){
		if(freq1.count(pr.first)){
			if(pr.second>Heavy_Thes)thres_large++;
			Join_Ground_Truth+=1ll*pr.second*freq1[pr.first];
		}
	}
	item_num = items.size();
	flow_num = freq0.size()+freq1.size();
	cout << "this" << items[0] << " " << items[1] << endl;
	cout << "dataset name: " << filename << endl;
	cout << flow_num << " flows, " << items.size() << " items read" << endl;;
	cout << "max freq = " << max_freq << endl;
	cout << "large flow sum = " << thres_large  << endl;
	cout <<"Join Ground Truth = "<<Join_Ground_Truth<< endl<< endl;
}

class interval{
    public:
        // [ ] 区间
        unsigned int begin;
        unsigned int end;
        interval(unsigned int a, unsigned int b){
            begin = a;
            end = b;
        }
};

vector<interval> itv;
vector<unsigned int>itv_items;
void generate_range(unsigned int min, unsigned int max, unsigned int intervalnum, unsigned int max_size){
    if(intervalnum%2 != 0) return;
    vector<unsigned int> v,p;
    for(int i = 0; i < intervalnum; ++i){
        unsigned int temp = rand();
        v.push_back(temp);
        temp = rand() % max_size;
        p.push_back(temp);
    }
    sort(v.begin(), v.end());
    for(int i = 0; i < intervalnum; i++){
        itv.push_back(interval(v[i],v[i]+p[i]));
        // cout << v[i] << "   " << v[i] + p[i] << endl;
    }
    
    unsigned int temp = 0;
    for(auto pr : itv){
        // cout << pr.begin << "  " << pr.end << endl;
        temp = pr.begin;
        while(temp < pr.end){
            itv_items.push_back(temp);
            temp++;
        }
    }
}
void generate_range_clean(){
    itv.clear();
    itv_items.clear();
}

void test_ske(int mem_in_byte,int D,int CHOOSE)
{
	int d = D;  //counts of hash function
	int w = mem_in_byte;  //   bits/counter_size/hash_counts

	long double _ARE = 0, _AAE = 0 , Mx=0,Mn=Join_Ground_Truth,X2=0,X=0,Con=0,INT=0,JOT=0;
	vector<double>all;
	all.clear();
	//generate the pseudo-random numbers for FAGMS sketches; use EH3 and CW2B
	for (int i = 0; i < testcycles; ++i)
	{
		Xi ***xi_pm1;
		Xi **xi_bucket;
		std::vector<int> counter_width;
		if(CHOOSE == 3){
			int _w = w;
			int heavysize=_w/8;
			heavysize=heavysize/64;
			_w-=heavysize*8*8;
			int midsize=_w/2;
			midsize/=4;
			int n=midsize/8;
			midsize=n*8;
			_w-=midsize*4;
			counter_width = closestPrimes(_w/d/2, d);
		}
		else counter_width = closestPrimes(w/d/2, d);

		unsigned int I1, I2;
		xi_pm1 = new Xi**[d];
		xi_bucket = new Xi*[d];
		for (int j = 0; j < d; j++)
		{	
			xi_pm1[j] = new Xi*[counter_width[j]];
			for (int k = 0; k < counter_width[j]; ++k){
				I1 = M_Random_Generate();
				I2 = M_Random_Generate();
				xi_pm1[j][k] = new Xi_EH3(I1, I2);
			}
			I1 = M_Random_Generate();
			I2 = M_Random_Generate();
			xi_bucket[j] = new Xi_CW4B(I1, I2, counter_width[j]);
		}

		Sketch *bcm1,*bcm0;
		bcm0=Choose_Sketch(w,d, i * 123,CHOOSE);
		bcm1=Choose_Sketch(w,d, i * 123,CHOOSE);
		
	    if(CHOOSE == 3){
			ClassifierSum* temp = dynamic_cast<ClassifierSum*>(bcm0);
			temp->Light->xi_pm1 = xi_pm1;
			temp = dynamic_cast<ClassifierSum*>(bcm1);
			temp->Light->xi_pm1 = xi_pm1;
		}
        
        // generate_range(0, UINT32_MAX, 4, 400000);

		clock_t st;
		st=clock();
		for (int j=0; j<items.size();j++){
		// if(i&1)
			if(j<data_delimiter)
				bcm1->Insert(items[j].c_str());
			else
				bcm0->Insert(items[j].c_str());
		}

		INT+=1.0*(clock()-st)/CLOCKS_PER_SEC;

        st = clock();
		long double my=bcm0->Join(bcm1);
		// cout << my << endl;
		JOT += 1.0*(clock()-st)/CLOCKS_PER_SEC;

		_AAE+=abs(my-Join_Ground_Truth);
		Mx=max(Mx,abs(my-Join_Ground_Truth));
		Mn=min(Mn,abs(my-Join_Ground_Truth));
		_ARE+=1.0*abs(my-Join_Ground_Truth)/Join_Ground_Truth;
		X+=my;
		X2+=my*my;
		Con+=(my-Join_Ground_Truth)*(my-Join_Ground_Truth);
		all.push_back((my-Join_Ground_Truth));
		
		for (int j = 0; j < d; j++){
			delete xi_pm1[j];
			delete xi_bucket[j];
		}
		delete [] xi_pm1;
		delete [] xi_bucket;
        // generate_range_clean();
	}
	_AAE/=testcycles;
	_ARE/=testcycles;
	Con/=testcycles;
	X/=testcycles;
	INT/=testcycles;
	JOT/=testcycles;
	X2/=testcycles;
	X=X2-X*X;
	experiment_AAE.push_back(_AAE);
	experiment_ARE.push_back(_ARE);
	experiment_MAX.push_back(Mx);
	experiment_MIN.push_back(Mn);
	experiment_MEM.push_back(mem_in_byte);
	experiment_CON.push_back(Con);
	experiment_YangBen.push_back(X);
	experiment_INT.push_back(INT);
	experiment_JOT.push_back(JOT);
	experiment_all.push_back(all);
	cout << "ARE = " << _ARE << ", AAE = " << _AAE << endl;
}


#include <boost/program_options.hpp>
using namespace boost::program_options;
void ParseArg(int argc, char *argv[],char*filename,int&len,int&sz,int&CHOOSE)
{
    options_description opts("Join Options");

    opts.add_options()
        ("memory,m", value<int>()->required(), "memory size")
        ("version,v", value<int>()->required(), "version")
        ("keylen,l", value<int>()->required(), "keylen")
        ("filename,f", value<string>()->required(), "trace file")
        ("help,h", "print help info")
        ;
    boost::program_options::variables_map vm;
    store(parse_command_line(argc, argv, opts), vm);
    

    if(vm.count("help"))
    {
        cout << opts << endl;
        exit(0);
    }
	if (vm.count("keylen"))
    {
        len=vm["keylen"].as<int>();
    }
    else
    {
        printf("please use -l to specify the keylen.\n");
        exit(0);
    }
	if (vm.count("version"))
    {
        CHOOSE=vm["version"].as<int>();
    }
    else
    {
        printf("please use -v to specify the algorithm to use.\n");
        exit(0);
    }
    if (vm.count("memory"))
    {
        sz=vm["memory"].as<int>();
    }
    else
    {
        printf("please use -m to specify the memory size.\n");
        exit(0);
    }
    if (vm.count("filename"))
    {
        strcpy(filename,vm["filename"].as<string>().c_str());
    }
    else
    {
        printf("please use -f to specify the trace file.\n");
        exit(0);
    }
}



int main(int argc,char *argv[]){
	char filename[100];
	int sz,CHOOSE,len;
    ParseArg(argc, argv,filename,len,sz,CHOOSE);
	if(len==13)readFile_CAIDA(filename,len);
	else readFile_zipf(filename,len);
	
	test_ske(sz,3,CHOOSE);
	ofstream oFile;
	char oFilename[50];
	sprintf(oFilename,"result.csv",CHOOSE);
	oFile.open(oFilename, ios::app);
	if (!oFile) return 0;
	// oFile << "Dataset:" << filename << endl;
	// oFile <<"Memory," << "AAE," << "ARE,"<<"VAR,"<<"Max,"<<"Min,"<<"Insert time,"<<"Join time,"<<"Join Ground Truth = "<<Join_Ground_Truth<< endl;
	cout<<experiment_ARE.size()<<endl;
	for (int o = 0; o < experiment_ARE.size(); o++){
		oFile << fixed << setprecision(15)<< experiment_MEM.at(o) << "," << experiment_AAE.at(o) << "," << experiment_ARE.at(o)<< "," << experiment_CON.at(o)<< "," << experiment_MAX.at(o)<< "," << experiment_MIN.at(o) << "," << experiment_INT.at(o) << "," << experiment_JOT.at(o);
		// for(auto item:experiment_all[o]){
		// 	oFile<<","<<item;
		// }
	}
    oFile << endl;
	oFile.close();
	experiment_MEM.clear();
	experiment_AAE.clear();
	experiment_ARE.clear();
	experiment_CON.clear();
	experiment_MAX.clear();
	experiment_MIN.clear();
	experiment_INT.clear();
	experiment_YangBen.clear();
	experiment_JOT.clear();
	experiment_all.clear();
	return 0;
}