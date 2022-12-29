#ifndef UTIL_H
#define UTIL_H
#include<cstdarg>
#include<cstdint>
#include<string>
#include<sstream>
#include<sys/time.h>
#include<time.h>
#define MAXCHARS 1000
#include <cmath>
#include <vector>
#include <stdint.h>
#include "using.h"

struct HashValue {
  uint64_t h0;
  uint64_t h1;
  uint64_t& operator[](int i){ return i==0? h0 : h1;}
  bool operator==(const HashValue& h) const {
    if(h0 != h.h0 || h1 != h.h1) return false;
    return true;
  }
};

std::ostream& operator <<(std::ostream& os, HashValue& h);
std::istream& operator >>(std::istream& is, HashValue& h);

bool has_char(string& s, char c);
//const char* get_bytes(const char* key, int& nbytes);
//const char* get_bytes(const string& key, int& nbytes);

struct my_cout {
  bool doit;
  my_cout(bool d = false) : doit(d){}; 
  void set(void){doit = true;}
  void unset(void) {doit = false;}
  void operator()(int n){
    if(doit) cout << n<<": ";
  }
};

double log_2(double x, double eps = 1.0e-8);

class LFSR128 { // implements the xorshift128+ algorithm
  uint64_t s[2];
public:
  LFSR128(void) {}
  LFSR128(uint64_t seed){ 
    if(seed == 0) throw "LFSR128: Seed must be non-zero.";
    step(seed);
    for(int i = 0;i < 1000;i++) step();
  }
  LFSR128(uint64_t* ss) {
    step(ss);
    for(int i = 0;i < 1000;i++) step();
  }
  uint64_t step(void){
	uint64_t x = s[0];
	uint64_t const y = s[1];
	s[0] = y;
	x ^= x << 23; // a
	s[1] = x ^ y ^ (x >> 17) ^ (y >> 26); // b, c
	return s[1] + y;
  }
  uint64_t step(uint64_t seed){
    s[0] = seed;
    s[1] = 0;
    return step();
  }
  uint64_t step(uint64_t* ss){
    s[0] = ss[0];
    s[1] = ss[1];
    return step();
  }
  double uniform(void){ return step()/(double)UINT64_MAX;}
  LFSR128& operator=(LFSR128& x){
    s[0] = x.s[0];
    s[1] = x.s[1];
    return *this;
  }
  void reseed(uint64_t* ss){
    s[0] = ss[0];
    s[1] = ss[1];
  }
};

template<typename T1, typename T2>
struct IndexPair {
  T1 item1;
  T2 item2;
  int sort_item; // NOTE: moved from static to private on 4/7/22
  IndexPair(void) {}
  IndexPair(const T1& t1, const T2& t2, int si = 1): item1(t1),item2(t2){sort_item = si;}
  bool operator <(const IndexPair& p)const {
    return (sort_item == 1? item1 < p.item1 : item2 < p.item2);
  }
  bool operator ==(const IndexPair& p) const{
    return (sort_item == 1? item1 == p.item1 : item2 == p.item2);
  }
  IndexPair& operator=(const IndexPair& ip){
    item1 = ip.item1;
    item2 = ip.item2;
    //    sort_item = ip.sort_item;
    return *this;
  }
  void set_sort_item(int i){sort_item = i;}
};



class PopCount {
  char byte_count[256];
 public:
  PopCount(void);
  uint operator()(uint64_t word);
};

class LCG64 { // 64-bit LCG. Parameters from Wikipedia
  static const uint64_t a;
  static const uint64_t b;
  static const double uint64_max;
  uint64_t x;
public:
  LCG64(uint64_t seed = 0){// : a(2862933555777941757), b(3037000493), uint64_max((double)UINT64_MAX){
    reseed(seed);
  }
  LCG64(const LCG64& lcg){// : a(2862933555777941757), b(3037000493), uint64_max((double)UINT64_MAX){
    x = lcg.x; // copy state
  }
  uint64_t step(void){return x = a*x+b;}
  double uniform(void){ return step()/uint64_max;}
  uint64_t reseed(uint64_t seed){
    x = seed;
    for(int i = 0;i < 10;i++)step(); // avoid crappy seeds
    return step();
  }
  LCG64& operator=(const LCG64& lcg){
    x = lcg.x;
    return *this;
  }
};

class BloomFilter { // Have we seen a 64-bit key or a string yet? Remember it for next time
  uint64_t word[4];
  LCG64 lcg;
  int _density;
 public:
  BloomFilter(void) { clear();}
  bool hash(uint64_t key, bool set = true){
    bool ret = true;
    uint64_t hash_val = lcg.reseed(key);
    for(int i = 0;i < 8;i++){
      int j = hash_val & 3; // low order 2-bits
      hash_val >>= 2;
      uint64_t mask  = (1<<(hash_val & (0x3f))); // next six bits
      hash_val >>= 6;
      if(  (mask&word[j]) == 0){ // bit not set
        ret = false; 
        if(set){ // set the bit
          word[j] |= mask; 
          _density++;
        }
      }
    }
    return ret;
  }
  bool hash(std::string& s, bool set = true){
    uint64_t hval = 0;
    for(int i = 0;i < s.size();i++){
      hval ^= (s[i] << (i%8));
    }
    return hash(hval,set);
  }
  void clear(void){
    for(int i = 0;i < 4;i++){
      word[i] = 0;
      _density = 0;
    }
  }
  bool has(uint64_t key){return hash(key,false);}
  bool has(std::string& s){return hash(s,false);}
  int density(void){ return _density; }
};
      
      


#define MAXFIELDS 100000  //Not easy to make field an Array because the header references are circular

class StringSplitter { // Functor to split a C-string into fields using field separator fs.  NOTE: leading
  // blanks in a field are skipped, so the field separator  is really fs followed by any number of blanks.
  // If fs == ' ', the splitter will also recognize '\t'.  If fs == '\t', ONLY '\t' is recognized.
  int nf; // no. of fields split out
  char* field[MAXFIELDS]; // array of ptrs to fields
public:
  StringSplitter(void):  nf(0) {} 
  int operator()(char* record,char fs = ' '); /* input a record, split into fields using field separator fs, 
                                               * return nf (Note: fs is over-written with '\0', making
                                               * record a sequence of C-strings)*/
  char* operator[](int i); // subsequent access to the fields

};


void split_string(char* line, std::vector<char*>& token, char fs = ' ');
void deblank(char* string);
std::string format(const char* fstr,...);
std::string str_time(double t);
std::string str_date(double t);
double get_time(int date, std::string time);
std::string currentDateTime(void);
class InternalTime { // floating pt seconds <-> 38-bit integer microseconds (3.18 days = 2^38 microseconds)
  int a;
  char b;
  inline void convert(double delta){
    double n = delta*1000000/128;
    a = (int) n; // high-order signed 31 bits
    b = (char) ((n-a)*128); // low-order signed 7 bits
  }
public:
  InternalTime(double delta = 0){ 
    if(fabs(delta) > 259200) throw format("InternalTime: bad input %lf\n",delta); // abs(delta) >= 3 days
    convert(delta);
  }
  inline InternalTime& operator=(double delta){convert(delta); return *this;}
  inline InternalTime& operator=(const InternalTime& t){a = t.a;b = t.b; return *this;}
  inline operator double() const{return (a*128.0 + b)/1000000.0;}
};

uint32_t str_hash(const char* s, uintptr_t salt = 0x55555555); // hash a c-string

class String : public std::string { // example of a useable KEY class (see KeyIndex.h) 
public:
  static String null;
  String(void) : std::string() {}
  String(const std::string& s) : std::string(s){}
  String(std::string& s) : std::string(s){}
  String(const char* p) : std::string(p){}
  String(char* p) : std::string(p){}
  uint32_t hash(uintptr_t salt = 0) const{return str_hash(this->c_str(),salt);}
  uint32_t hash1(uintptr_t salt = 0) const; // alt hash function
  operator int() const{ return hash();}
  template<typename T>
  String(T x){ // also add type conversions from int, double, char, and other classes 
    std::stringstream iss;
    iss << x;
    *this = iss.str();
  }
};

class Welford { // online computation of mean and variance
  double S1;
  double S2;
  double W;
  double tau; // exponential decay coefficient
  std::string ser;
public:
  Welford(double t = 1.0) : W(0),S1(0),S2(0),tau(t) {}
  void reset(double t = 1.0){W = S1 = S2 = 0;tau = t;}
  std::string& serialize(void){
    ser = format("%.18f %.18f %.18f %.18f",S1,S2,W,tau); // for PersistentState
    return ser;
  }
  Welford(std::string& s){ restart(s);}
  void restart(std::string& s){ // restart from serialization
    std::istringstream iss(s);
    iss >> S1 >> S2 >> W >> tau;
  }
  void update(double x, double weight = 1.0); // add another data point
  double tot_weight(void){return W;}
  double mean(void){return W == 0? 0:S1/W;}
  double variance(void){return W == 0? 0 : S2/W;} // NOTE: for correctness, should be n-1, but WGAS.
};


template <typename T1, typename T2>
ostream& operator<<(ostream& os, const IndexPair<T1,T2>& p){
  os << p.item1 <<": "<< p.item2 << endl;
  return os;
}

#endif
