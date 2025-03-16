#include "MurmurHash.h"
# include "params.h"
# include <iostream>
# include <string.h>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <immintrin.h>
#include <stdexcept>
#include "xis.h"
#include "utils.h"
#ifdef UNIX
#include <x86intrin.h>
#else
#include <emmintrin.h>
#endif
using namespace std;
extern int Heavy_Thes;
extern int BILI;
inline int SIMD_match_8_S(uint32_t *ID,uint32_t key){
	const __m128i item = _mm_set1_epi32((int)key);
	int matched;
	__m128i *keys_p = (__m128i *)ID;
	__m128i a_comp = _mm_cmpeq_epi32(item, keys_p[0]);
	__m128i b_comp = _mm_cmpeq_epi32(item, keys_p[1]);
	a_comp = _mm_packs_epi32(a_comp, b_comp);
	matched = _mm_movemask_epi8(a_comp);
	if (matched != 0) {
		int matched_index = __builtin_ctz((uint32_t)matched)>>1;
		return matched_index;
    }else return -1;
}
struct Hash_table_S{
	int n,hash_seed;const int m=8;
	unsigned **Id,**Key;
	Hash_table_S(int bucket_sum,int _hash_seed=1000){
		n=bucket_sum;
		Id=new unsigned*[n]();
		for(int i=0;i<n;i++)
			Id[i]=new unsigned[m]();
		Key=new unsigned*[n]();
		for(int i=0;i<n;i++)
			Key[i]=new unsigned[m]();
		for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
		Id[i][j]=Key[i][j]=0;
	}
	bool count(unsigned id)const{
		unsigned i=MurmurHash32(&id, 4, hash_seed)%n;
		if(SIMD_match_8_S(Id[i],id)>=0)return 1;
		i=MurmurHash32(&id, 4, hash_seed*2)%n;
		if(SIMD_match_8_S(Id[i],id)>=0)return 1;
		return 0;
	}
	int query(unsigned id)const{
		unsigned i=MurmurHash32(&id, 4, hash_seed)%n;
		int j=SIMD_match_8_S(Id[i],id);
		if(j>=0)return Key[i][j];
		i=MurmurHash32(&id, 4, hash_seed*2)%n;
		j=SIMD_match_8_S(Id[i],id);
		if(j>=0)return Key[i][j];
		return -1;
	}
	void insert(unsigned id,unsigned cnt=1){
		unsigned i=MurmurHash32(&id, 4, hash_seed)%n;
		int j;
		j=SIMD_match_8_S(Id[i],id);
		if(j>=0){
			Key[i][j]+=cnt;
			return;
		}
		i=MurmurHash32(&id, 4, hash_seed*2)%n;
		j=SIMD_match_8_S(Id[i],id);
		if(j>=0){
			Key[i][j]+=cnt;
			return;
		}
		for(int j=0;j<m;j++)
		if(!solid(i,j)){
			Id[i][j]=id;
			Key[i][j]=cnt;
			return;
		}
		i=MurmurHash32(&id, 4, hash_seed)%n;
		for(int j=0;j<m;j++)
		if(!solid(i,j)){
			Id[i][j]=id;
			Key[i][j]=cnt;
			return;
		}
		moretime:Expansion();
		i=MurmurHash32(&id, 4, hash_seed)%n;
		for(int j=0;j<m;j++)
		if(!solid(i,j)){
			Id[i][j]=id;
			Key[i][j]=cnt;
			return;
		}
		i=MurmurHash32(&id, 4, hash_seed*2)%n;
		for(int j=0;j<m;j++)
		if(!solid(i,j)){
			Id[i][j]=id;
			Key[i][j]=cnt;
			return;
		}
		goto moretime;
		exit(1);
	}
	bool solid(unsigned i,unsigned j){
		if(Key[i][j]==0)return 0;
		unsigned id=Id[i][j];
		unsigned ii=MurmurHash32(&id, 4, hash_seed)%n;
		if(ii==i)return 1;
		ii=MurmurHash32(&id, 4, hash_seed*2)%n;
		if(ii==i)return 1;
		return 0;
	}
	void Expansion(){
		unsigned **nwid;
		unsigned **nwkey;
		nwid=new unsigned*[n+n]();
		nwkey=new unsigned*[n+n]();
		for(int i=0;i<n;i++){
			nwid[i]=Id[i];
			nwkey[i]=Key[i];
		}
		for(int i=0;i<n;i++){
			nwid[n+i]=new unsigned[m]();
			nwkey[n+i]=new unsigned[m]();
			memmove(nwid[n+i],nwid[i],m*sizeof(unsigned));
			memmove(nwkey[n+i],nwkey[i],m*sizeof(unsigned));
		}
		delete[]Id;
		delete[]Key;
		Id=nwid;
		Key=nwkey;
		n=n+n;
	}
	~Hash_table_S()
	{
		for (int i = 0; i < n; i++){
			delete[]Id[i];
			delete[]Key[i];
		}
		delete[]Id;
		delete[]Key;
	}
};

class Short_SF_Sketch
{
public:
	int w, d;
	int index[MAX_HASH_NUM];
	int* counter[MAX_HASH_NUM];
	int MAX_CNT, MIN_CNT, hash_seed;
	vector<int> counter_width;
	Xi ***xi_pm1;
public:
	Short_SF_Sketch(int _w, int _d, int _hash_seed = 1000)
	{
		d = _d, w = _w/_d/2;
		hash_seed = _hash_seed;
		counter_width = closestPrimes(w, d);
		for (int i = 0; i < d; i++)
		{
			counter[i] = new int[counter_width[i]]();
			memset(counter[i], 0, sizeof(int) * counter_width[i]);
		}
		MAX_CNT = INT32_MAX;
		MIN_CNT = INT32_MIN;
	}
	void InsertRange(unsigned int begin, unsigned int end){
		
		for(int i = 0; i < d; ++i){
			unsigned int alpha = begin / counter_width[i] + 1, beta = end / counter_width[i];
			unsigned int b_bucket = begin % counter_width[i], e_bucket = end % counter_width[i];
			
			for(int j = 0; j < counter_width[i]; ++j){
				if(j == b_bucket) alpha--;
				else if(j == e_bucket) beta--;
				if(alpha > beta) continue;
				else if(alpha == beta) counter[i][j] += xi_pm1[i][j]->element(alpha);
				else counter[i][j] += xi_pm1[i][j]->interval_sum(alpha, beta);
			}
		}
	}
	void Insert(unsigned id,int k=1)
	{
		int g = 0;
		for (int i = 0; i < d; i++)
		{
			// index[i] = ((unsigned)MurmurHash32(&id, 4, hash_seed + i)) % w;
			index[i] = id % counter_width[i];
			if(counter[i][index[i]] != MAX_CNT && counter[i][index[i]] != MIN_CNT){
				counter[i][index[i]] += xi_pm1[i][0]->element(id / counter_width[i]) * k;
			}
		}
	}
	long double Join(Short_SF_Sketch*other){
		long double res[MAX_HASH_NUM];
		for (int i = 0; i < d; i++){
			long double k=0;
			for (int j = 0; j < counter_width[i]; j++)
				k+=1ll*counter[i][j]*other->counter[i][j];
			res[i]=k;
		}
		sort(res, res + d);
		if (d % 2 == 0)
			return ((res[d / 2] + res[d / 2 - 1]) / 2);
		else
			return (res[d / 2]);
	}
	long double JoinWhere(Short_SF_Sketch*other, Short_SF_Sketch*where){
		long double res[MAX_HASH_NUM];
		for (int i = 0; i < d; i++){
			long double k=0;
			for (int j = 0; j < counter_width[i]; j++)
				k+=1ll*counter[i][j]*other->counter[i][j]*where->counter[i][j];
			res[i]=k;
		}
		sort(res, res + d);
		if (d % 2 == 0)
			return ((res[d / 2] + res[d / 2 - 1]) / 2);
		else
			return (res[d / 2]);
	}
	int Query(unsigned id)const
	{
		int temp;
		int res[MAX_HASH_NUM];
		int index[MAX_HASH_NUM];
		int g;
		for (int i = 0; i < d; i++)
		{
			// index[i] = ((unsigned)MurmurHash32(&id, 4, hash_seed + i)) % w;
			index[i] = id % counter_width[i];
			temp = counter[i][index[i]];
			g = xi_pm1[i][0]->element(id / counter_width[i]);
			res[i] = (g == 1 ? temp : -temp);
		}

		sort(res, res + d);
		if (d % 2 == 0)
			return ((res[d / 2] + res[d / 2 - 1]) / 2);
		else
			return (res[d / 2]);
	}

	~Short_SF_Sketch()
	{
		for (int i = 0; i < d; i++)
			delete[]counter[i];
	}
};

class ClassifierSum:public Sketch
{
public:
	int n;const int m=8;
	struct data{
		unsigned Id,Key;
		data(unsigned a=0,unsigned b=0):Id(a),Key(b){}
		bool operator<(const data&other)const{
			if(Key!=other.Key)return Key<other.Key;
			return Id<other.Id;
		}
	}**Mid;
	int hash_seed;
	Short_SF_Sketch *Light;
	Hash_table_S *Heavy;
public:
	ClassifierSum(int _w, int _d, int _hash_seed = 1000)
	{
		int heavysize=_w/8;
		heavysize=heavysize/64;
		Heavy=new Hash_table_S(heavysize,_hash_seed*2);
		_w-=heavysize*8*8;
		int midsize=_w/2;
		midsize/=4;
		n=midsize/m;
		midsize=n*m;
		Mid=new data*[n];
		for(int i=0;i<n;i++)
		Mid[i]=new data[m]();
		_w-=midsize*4;
		hash_seed=_hash_seed*3;
		
		Light=new Short_SF_Sketch(_w,_d,_hash_seed);
	}
	void Insert(const unsigned id){
		if(Heavy->count(id)){Heavy->insert(id);return;}
		unsigned k=MurmurHash32(&id, 4, hash_seed)%n;
		data*bucket=Mid[k];
		for(int i=0;i<m;i++)
		if(bucket[i].Id==id||bucket[i].Key==0){
			bucket[i].Id=id;
			bucket[i].Key++;
			if(bucket[i].Key==Heavy_Thes){
				Heavy->insert(id,bucket[i].Key);
				bucket[i].Id=bucket[i].Key=0;
			}
			return;
		}
		int to=0,mn=bucket[0].Key;
		for(int i=1;i<m;i++)
		if(bucket[i].Key<mn){
			mn=bucket[i].Key;
			to=i;
		}
		Light->Insert(bucket[to].Id,bucket[to].Key);
		bucket[to].Key=1;
		bucket[to].Id=id;
		return;
	}
	void Insert(const char* str)
	{
		unsigned id=((unsigned*)str)[0];
		Insert(id);
		return;
	}
	bool CheckHeavy(const char* str){
		unsigned id=((unsigned*)str)[0];
		return Heavy->count(id);
	}
	long double Join(Sketch*_other){
		ClassifierSum* other=(ClassifierSum*)_other;
		long double re=Light->Join(other->Light);

		for(int i=0;i<Heavy->n;i++)
		for(int j=0;j<Heavy->m;j++)
		if(Heavy->solid(i,j)){
			re+=1.0*Heavy->Key[i][j]*other->Query(Heavy->Id[i][j]);
		}

		for(int i=0;i<other->Heavy->n;i++)
		for(int j=0;j<other->Heavy->m;j++)
		if(other->Heavy->solid(i,j)){
			re+=1.0*QueryMid(other->Heavy->Id[i][j])*other->Heavy->Key[i][j];
		}

		for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
		if(Mid[i][j].Key!=0){
			re+=1.0*Mid[i][j].Key*other->QueryMid(Mid[i][j].Id);
		}

		for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
		if(other->Mid[i][j].Key!=0){
			re+=1.0*other->Mid[i][j].Key*Light->Query(other->Mid[i][j].Id);
		}
		return re;
	}
	int QueryMid(const unsigned id)const {
		int re=Light->Query(id);
		unsigned k=MurmurHash32(&id, 4, hash_seed)%n;
		data*bucket=Mid[k];
		for(int i=0;i<m;i++)
		if(bucket[i].Id==id)re+=bucket[i].Key;
		return re;
	}
	int Query(const unsigned id)const{
		int re=Light->Query(id);
		if(Heavy->count(id))re+=Heavy->query(id); 
		unsigned k=MurmurHash32(&id, 4, hash_seed)%n;
		data*bucket=Mid[k];
		for(int i=0;i<m;i++)
		if(bucket[i].Id==id)re+=bucket[i].Key;
		return re;
	}
	int Query(const char* str)const
	{
		// unsigned id= (unsigned)MurmurHash32(str, KEY_LEN, fingerprint);
		unsigned id=((unsigned*)str)[0];
		return Query(id);
	}



	int QueryLight(const char* str)const
	{
		unsigned id=((unsigned*)str)[0];
		return Light->Query(id);
	}
	int QueryHeavy(const char* str)const
	{
		unsigned id=((unsigned*)str)[0];
		if(Heavy->count(id)) return Heavy->query(id);
		else return 0;
	}
	int QueryMidwithoutLight(const char* str)const
	{
		unsigned id=((unsigned*)str)[0];
		return QueryMid(id) - Light->Query(id);
	}
	int QueryLight(const unsigned id)const
	{
		return Light->Query(id);
	}
	int QueryHeavy(const unsigned id)const
	{
		if(Heavy->count(id)) return Heavy->query(id);
		else return 0;
	}
	int QueryMidwithoutLight(const unsigned id)const
	{
		return QueryMid(id) - Light->Query(id);
	}


	void JoinSp(Sketch *_other){
		ClassifierSum* other=(ClassifierSum*)_other;
		long double re = Light->Join(other->Light);
		cout << "light light :" <<  re << endl;

		long double hh = 0, hm = 0, hl = 0;
		for(int i=0;i<Heavy->n;i++)
		for(int j=0;j<Heavy->m;j++)
		if(Heavy->solid(i,j)){
			hh += 1.0*Heavy->Key[i][j]*other->QueryHeavy(Heavy->Id[i][j]);
			hm += 1.0*Heavy->Key[i][j]*other->QueryMidwithoutLight(Heavy->Id[i][j]);
			hl += 1.0*Heavy->Key[i][j]*other->QueryLight(Heavy->Id[i][j]);
		}
		cout << "heavy heavy :" <<  hh << endl;
		cout << "heavy mid :" <<  hm << endl;
		cout << "heavy light :" <<  hl << endl;

		long double mh = 0, lh = 0;
		for(int i=0;i<other->Heavy->n;i++)
		for(int j=0;j<other->Heavy->m;j++)
		if(other->Heavy->solid(i,j)){
			mh += 1.0*QueryMidwithoutLight(other->Heavy->Id[i][j])*other->Heavy->Key[i][j];
			lh += 1.0*QueryLight(other->Heavy->Id[i][j])*other->Heavy->Key[i][j];
		}
		cout << "mid heavy :" <<  mh << endl;
		cout << "light heavy :" <<  lh << endl;

		long double ml = 0, mm = 0;
		for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
		if(Mid[i][j].Key!=0){
			mm += 1.0*Mid[i][j].Key*other->QueryMidwithoutLight(Mid[i][j].Id);
			ml += 1.0*Mid[i][j].Key*other->QueryLight(Mid[i][j].Id);
		}
		cout << "mid mid :" <<  mm << endl;
		cout << "mid light :" <<  ml << endl;

		long double lm = 0;
		for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
		if(other->Mid[i][j].Key!=0){
			lm += 1.0*other->Mid[i][j].Key*Light->Query(other->Mid[i][j].Id);
		}
		cout << "light mid :" <<  lm << endl;
	}

	long double JoinWhere(Sketch *_other, Sketch *_where, unsigned int begin, unsigned end){
		ClassifierSum* other=(ClassifierSum*)_other;
		ClassifierSum* where=(ClassifierSum*)_where;

		long double re=Light->JoinWhere(other->Light, where->Light);

		for(int i=0;i<Heavy->n;i++)
		for(int j=0;j<Heavy->m;j++)
		if(Heavy->solid(i,j)){
			if(Heavy->Id[i][j] >= begin &&Heavy->Id[i][j] < end)
				re+=1.0*Heavy->Key[i][j]*other->Query(Heavy->Id[i][j]);
		}

		for(int i=0;i<other->Heavy->n;i++)
		for(int j=0;j<other->Heavy->m;j++)
		if(other->Heavy->solid(i,j)){
			if(other->Heavy->Id[i][j] >= begin &&other->Heavy->Id[i][j] < end)
				re+=1.0*QueryMid(other->Heavy->Id[i][j])*other->Heavy->Key[i][j];
		}

		for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
		if(Mid[i][j].Key!=0){
			if(Mid[i][j].Id >= begin &&Mid[i][j].Id < end)
				re+=1.0*Mid[i][j].Key*other->QueryMid(Mid[i][j].Id);
		}

		for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
		if(other->Mid[i][j].Key!=0){
			if(other->Mid[i][j].Id >= begin && other->Mid[i][j].Id < end)
				re+=1.0*other->Mid[i][j].Key*Light->Query(other->Mid[i][j].Id);
		}
		return re;
	}

	~ClassifierSum()
	{
		for (int i = 0; i < n; i++)
			delete[]Mid[i];
		delete[]Mid;
		Heavy->~Hash_table_S();
		Light->~Short_SF_Sketch();
	}
};