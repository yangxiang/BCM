#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "WG.hh"

template <typename T>
WG<T>::WG(uint64_t num_trans, uint64_t min_item, uint64_t max_item)
    : db(NULL), num_trans(num_trans), min_item(min_item), max_item(max_item), _item_width(max_item - min_item + 1)
{
  try
  {
	  db= new T* [num_trans];
	  for (uint64_t i=1; i<=num_trans; i++)
	  {
		db[i-1]=new T[_item_width];
	  }
  }
  catch (int e)
  {
	std::cout<<"An exception of intial allocation of memory: "<<e<<std::endl;
	exit(-1);
  }
  
}


template <typename T>
WG<T>::~WG()
{
  try
  {
	  for (uint64_t i=1; i<=num_trans; i++)
	  {
		delete [] db[i-1];
	  }
	  delete [] db;
  }
  catch (int e)
  {
	std::cout<<"An exception of final release of memory: "<<e<<std::endl;
	exit(-1);
  }
}

template <typename T> WG<T>::WG(uint64_t num_trans, uint64_t min_item, uint64_t max_item);
template <typename T> WG<T>::~WG();

