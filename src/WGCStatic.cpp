#include <assert.h>
#include <limits>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "WGCStatic.hh"

std::string WGCStatic::ds_fpath;
Itemset::item_t WGCStatic::max_item = 0;
Itemset::item_t WGCStatic::min_item = 0;
Transactionset::transaction_t WGCStatic::num_transactions = 0;
double WGCStatic::beta = -1.0;
unsigned int WGCStatic::merge_k = std::numeric_limits<unsigned int>::max();
int64_t WGCStatic::start_usec = 0;

int64_t t_usec()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return (int64_t)tv.tv_sec * pow(10, 6) + tv.tv_usec;
}

double WGCStatic::elapsed_time_sec()
{
  return (double)(t_usec() - WGCStatic::start_usec) / pow(10, 6);
}

void WGCStatic::init_timer()
{
  WGCStatic::start_usec = t_usec();
}
