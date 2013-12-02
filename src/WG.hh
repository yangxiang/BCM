#ifndef __WG
#define __WG

template <typename T>
class WG
{

	friend std::ostream& operator<< (std::ostream& os, const WG<T>& tip_db)
	{
	  for (std::size_t i=1; i<=tip_db.num_trans; ++i)
	  {
		for (std::size_t j=tip_db.min_item; j<=tip_db.max_item; ++j)
		{
		  std::cout << tip_db.get_weight(i, j)<<", ";
		}
		std::cout << "\n";
	  }
	  return os;
	}

public:
  WG(uint64_t num_trans, uint64_t min_item, uint64_t max_item);
  ~WG();

  void set_weight(int transaction, int item, T weight);
  
  T get_weight(int transaction, int item) const;

protected:

  T **db;
  uint64_t num_trans;
  uint64_t min_item;
  uint64_t max_item;
  uint64_t _item_width;
};


template <typename T>
inline
void WG<T>::set_weight(int transaction, int item, T weight)
{
	db[transaction-1][item-min_item]=weight;
}
  

template <typename T>
inline
T WG<T>::get_weight(int transaction, int item) const
{

  return db[transaction-1][item-min_item];
}


#endif
