#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include <pthread.h>

#include <queue>
#include <list>

using namespace std;

#define encode(ch) code[ch-'A']

constexpr int NUM_THREADS = 4;
constexpr int LEN = 6;
constexpr int AA_NUMBER = 20;
constexpr double EPSILON = 1e-010;

static long M, M1, M2;
const short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3};

struct Bacteria
{
  int index;
  char* filename;
  unsigned int count;
  double* tv;
  long* ti;
  
  Bacteria(int i, char* name) : index(i), filename(name) {}
};

struct Comparison
{
  Bacteria* b1;
  Bacteria* b2;
  
  Comparison(Bacteria* b1, Bacteria* b2) : b1(b1), b2(b2) {}
};

static int number_bacteria;
static queue<Bacteria*> bacteria_waiting;
static list<Bacteria*> bacteria_done;
static queue<Comparison*> compare_waiting;
static char** output;

static int thread_ids[NUM_THREADS];
static pthread_t threads[NUM_THREADS];
static pthread_mutex_t bacteria_waiting_mutex = PTHREAD_RECURSIVE_MUTEX_INITIALIZER;
static pthread_mutex_t bacteria_done_mutex = PTHREAD_RECURSIVE_MUTEX_INITIALIZER;
static pthread_mutex_t compare_waiting_mutex = PTHREAD_RECURSIVE_MUTEX_INITIALIZER;

void* worker(void* data)
{
  // BACTERIA LOADING/COMPUTATION
  Bacteria* bact;
  FILE* bact_file;
  
  // main input vectors
  int* vector = new int[M];
  long* second = new long[M1];
  double* second_div = new double[M1];
  
  long one_l[AA_NUMBER];
  double one_l_div[AA_NUMBER];
  
  pthread_mutex_lock(&bacteria_waiting_mutex);
  
  while (!bacteria_waiting.empty())
  {
    bact = bacteria_waiting.front();
    bacteria_waiting.pop();
    pthread_mutex_unlock(&bacteria_waiting_mutex);
    
    memset(vector, 0, M * sizeof(long));
    memset(second, 0, M1 * sizeof(long));
    memset(one_l, 0, AA_NUMBER * sizeof(long));
    
    long indexs;
    long total = 0;
    long total_l = 0;
    long complement = 0;
    
    char ch;
    char buffer[LEN - 1];
    short enc;
    bact_file = fopen(bact->filename, "r");
    
    while ((ch = fgetc(bact_file)) != EOF)
    {
      if (ch == '>')
      {
        while (fgetc(bact_file) != '\n'); // skip rest of line
        fread(buffer, sizeof(char), LEN - 1, bact_file);
        
        complement++;
        indexs = 0;
        
        for (int i = 0; i < LEN - 1; ++i)
        {
          enc = encode(buffer[i]);
          one_l[enc]++;
          total_l++;
          indexs = indexs * AA_NUMBER + enc;
        }
        
        second[indexs]++;
      }
      else if (ch != '\n')
      {
        enc = encode(ch);
        one_l[enc]++;
        total_l++;
        vector[indexs * AA_NUMBER + enc]++;
        total++;
        indexs = (indexs % M2) * AA_NUMBER + enc;
        second[indexs]++;
      }
    }
    
    fclose(bact_file);
    
    long total_plus_complement = total + complement;
    
    for (int i = 0; i < AA_NUMBER; ++i)
    {
      one_l_div[i] = (double) one_l[i] / total_l;
    }
    
    for (int i = 0; i < M1; ++i)
    {
      second_div[i] = (double) second[i] / total_plus_complement;
    }
    
    double half_total = total * 0.5;
    int i_mod_aa_number = 0;
    int i_div_aa_number = 0;
    long i_mod_M1 = 0;
    long i_div_M1 = 0;
    
    unsigned int count = 0;
    double* tv = new double[M];
    long* ti = new long[M]; 
    
    for (int i = 0; i < M; ++i)
    {
      double stochastic = (second_div[i_div_aa_number] * one_l_div[i_mod_aa_number] +
                           second_div[i_mod_M1] * one_l_div[i_div_M1]) * half_total;
      
      if (i_mod_aa_number == AA_NUMBER - 1)
      {
        i_mod_aa_number = 0;
        i_div_aa_number++;
      }
      else
      {
        i_mod_aa_number++;
      }
      
      if (i_mod_M1 == M1 - 1)
      {
        i_mod_M1 = 0;
        i_div_M1++;
      }
      else
      {
        i_mod_M1++;
      }
      
      if (stochastic > EPSILON) 
      {
        tv[count] = (vector[i] - stochastic) / stochastic;
        ti[count] = i;
        count++;
      }
    }
    
    realloc(tv, count * sizeof(double));
    realloc(ti, count * sizeof(long));
    bact->tv = tv;
    bact->ti = ti;
    bact->count = count;
    
    pthread_mutex_lock(&bacteria_done_mutex);
    pthread_mutex_lock(&compare_waiting_mutex);
    
    for (list<Bacteria*>::iterator it = bacteria_done.begin(); it != bacteria_done.end(); ++it)
    {
      compare_waiting.push(new Comparison(*it, bact));
    }
    
    pthread_mutex_unlock(&compare_waiting_mutex);
    bacteria_done.push_back(bact);
    pthread_mutex_unlock(&bacteria_done_mutex);
    
    pthread_mutex_lock(&bacteria_waiting_mutex);
  }
  
  pthread_mutex_unlock(&bacteria_waiting_mutex);
  
  // CLEANUP
  
  delete[] vector;
  delete[] second;
  delete[] second_div;
  
  // BACTERIA COMPARISON
  
  Comparison* comp;
  Bacteria* b1;
  Bacteria* b2;
  char* output_line;
  
  pthread_mutex_lock(&compare_waiting_mutex);
  
  while (!compare_waiting.empty())
  {
    comp = compare_waiting.front();
    compare_waiting.pop();
    pthread_mutex_unlock(&compare_waiting_mutex);
    
    b1 = comp->b1;
    b2 = comp->b2;
    if (b1->index > b2->index) swap(b1, b2); // to preserve original output order
    
    double correlation = 0;
    double vec_len1 = 0;
    double vec_len2 = 0;
    unsigned int p1 = 0;
    unsigned int p2 = 0;
    long n1, n2;
    double t1, t2;
    
    while (p1 < b1->count && p2 < b2->count)
    {
      n1 = b1->ti[p1];
      n2 = b2->ti[p2];
      
      if (n1 < n2)
      {
        t1 = b1->tv[p1++];
        vec_len1 += t1 * t1;
      }
      else if (n2 < n1)
      {
        t2 = b2->tv[p2++];
        vec_len2 += t2 * t2;
      }
      else
      {
        t1 = b1->tv[p1++];
        t2 = b2->tv[p2++];
        vec_len1 += t1 * t1;
        vec_len2 += t2 * t2;
        correlation += t1 * t2;
      }
    }
    
    while (p1 < b1->count)
    {
      t1 = b1->tv[p1++];
      vec_len1 += (t1 * t1);
    }
    
    while (p2 < b2->count)
    {
      t2 = b2->tv[p2++];
      vec_len2 += (t2 * t2);
    }
    
    output_line = new char[40];
    double corr_val = correlation / (sqrt(vec_len1) * sqrt(vec_len2));
    sprintf(output_line, "%2d %2d -> %.20lf\n", b1->index, b2->index, corr_val);
    output[b1->index * number_bacteria + b2->index] = output_line;
    delete comp;
  }

  return NULL;
}

inline void init_M()
{
  M2 = 1;
  
  for (int i = 0; i < LEN - 2; ++i)
  {
    M2 *= AA_NUMBER;
  }
  
  M1 = M2 * AA_NUMBER;
  M = M1 * AA_NUMBER;
}

inline void read_input_file(const char* list_name)
{
  static char* prefix = "data/";
  FILE* list_file = fopen(list_name, "r");
  fscanf(list_file, "%d", &number_bacteria);
  char buff[20];
  char* bacteria_name;
  
  output = (char**) calloc(number_bacteria * number_bacteria, sizeof(char*));
  
  for (long i = 0; i < number_bacteria; ++i)
  {
    bacteria_name = new char[25];
    strcpy(bacteria_name, prefix);
    fscanf(list_file, "%s", buff);
    strcat(bacteria_name, buff);
    strcat(bacteria_name, ".faa");
    bacteria_waiting.push(new Bacteria(i, bacteria_name));
  }
  
  fclose(list_file);
}

inline void create_threads()
{
  for (int i = 0; i < NUM_THREADS; ++i)
  {
    thread_ids[i] = i;
    pthread_create(&threads[i], NULL, worker, NULL);
  }
  
  for (int i = 0; i < NUM_THREADS; ++i)
  {
    pthread_join(threads[i], NULL);
  }
}

inline void display_output()
{
  for (int i = 0; i < number_bacteria * number_bacteria; ++i)
  {
    if (output[i])
    {
      fputs(output[i], stdout);
      delete output[i];
    }
  }
}

inline void final_cleanup()
{
  free(output);
  
  Bacteria* bact;
  for (list<Bacteria*>::iterator it = bacteria_done.begin(); it != bacteria_done.end(); ++it)
  {
    bact = *it;
    delete bact->filename;
    delete bact->ti;
    delete bact->tv;
    delete bact;
  }
}

double time_diff(struct timespec* t1, struct timespec* t2)
{
  return (double) ((t2->tv_sec - t1->tv_sec) + (t2->tv_nsec - t1->tv_nsec) / 1000000000.0L);
}

int main(int argc, char** argv)
{
  struct timespec t1;
  struct timespec t2;
  clock_gettime(CLOCK_MONOTONIC, &t1);
  
  init_M();
  read_input_file(argv[1]);
  create_threads();
  display_output();
  final_cleanup();
  
  clock_gettime(CLOCK_MONOTONIC, &t2);
  printf("Total: %.3f seconds\n", time_diff(&t1, &t2));
  return 0;
}
