#include <algorithm>
#include <assert.h>
#include <errno.h>
#include <ext/algorithm>
#include <ext/functional>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <math.h>
#include <list>
#include <queue>
#include <string>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include "BitDb.hh"
#include "WG.cpp"
#include "CartesianProductDb.hh"
#include "WGCStatic.hh"
#define MaxTransRatio 1.0
#define MaxItemsRatio 1.0

typedef float WeightMeasure;

float TransOverItems=1.0; 

struct weightcomp {
  bool operator() (const WeightMeasure& lhs, const WeightMeasure& rhs) const
  {return lhs>rhs;}
};

struct sizecomp {
  bool operator() (const Transactionset::size_type& lhs, const Transactionset::size_type& rhs) const
  {return lhs>rhs;}
};

struct edgeinfo{
	Transactionset::transaction_t v_a;
	Itemset::item_t v_b;
};

static void usage() {
	std::cout << "\nUsage:\n"
		"	bicomp [-h] [-g gamma] [-c C_para] [-t tau_para]  [-a absolution] input_data_file output_pattern_file\n"
		"Description:\n"
		"	-h	Print the help message.\n"
		"	-g	gamma; default value 0.7; range 0<gamma && gamma<=1.\n"
		"	-c	C_para: integer; default value 25; range C_para>(tau_para+2)^2\n"
		"	-t	tau_para: integer; default value 2; range tau_para>=0.\n"
		"	-a	abs_choice: default value 0 (do not convert); range abs_choice={0,1}.\n"
		<< std::endl;
}

int main(int argc, char **argv)
{
  if (argc == 1)
  {
	usage();
	return 1;
  }
  
  struct timeval tv;
  gettimeofday(&tv, 0);
  srand48(tv.tv_sec + tv.tv_usec);

  WGCStatic::init_timer();

  //measure redone
  //Accuracy redone (contribution_a and contribution_b are original weights)
  
  std::string output_file;
  
  //default values begin
  double gamma=0.7; //range 0<gamma && gamma<=1  
  int C_para=25; //range C_para>(tau_para+2)^2\n
  int tau_para=2; //range range tau_para>=0
  int abs_choice=0; //range abs_choice={0,1}
  //default values end	
  
  
  int input_para_counter=1;
  while(input_para_counter<argc){
	if(strcmp("-h",argv[input_para_counter])==0){
		usage();
		return 1;
	}

	if(strcmp("-g",argv[input_para_counter])==0){
		input_para_counter++;
		gamma=atof(argv[input_para_counter++]);
	}	
	else if(strcmp("-c",argv[input_para_counter])==0){
		input_para_counter++;
		C_para=atoi(argv[input_para_counter++]);
	}
	else if(strcmp("-t",argv[input_para_counter])==0){
		input_para_counter++;
		tau_para=atoi(argv[input_para_counter++]);
	}
	else if(strcmp("-a",argv[input_para_counter])==0){
		input_para_counter++;
		abs_choice=atoi(argv[input_para_counter++]);
	}
	else{
		WGCStatic::ds_fpath= argv[input_para_counter++];
		output_file=argv[input_para_counter++];
		break;
	}
  }  
  
  //range check begin
  if (!(tau_para>=0 && C_para> (tau_para+2)*(tau_para+2) && (abs_choice==0 || abs_choice==1) ))
  {
	std::cout<<"At least one input parameter is out of range"<<std::endl;
	usage();
	exit(-1);
  }
  //range check done 

  //managing input and output files
  const std::string fi_fname(WGCStatic::ds_fpath);

  FILE *fi_file;
  if ((fi_file = fopen(fi_fname.c_str(), "rb")) == NULL)
  {

    std::cout << "Dataset " << fi_fname << " does not exist. Exit.\n";

    exit(1);
  }

  std::ostringstream out_wcdb;
  out_wcdb<<output_file;
  std::ofstream wcdb_outp(out_wcdb.str().c_str());

  if(!wcdb_outp)
  {	
	std::cout<<"Fail to create output file. Exit."<<std::endl;
	exit(-1);
  }
  //managing input and output files done
  
  std::cout << "Skim dataset file to determine its size: " << fi_fname << "\n";
  
  uint64_t Number_of_Rows=0;
  uint64_t Number_of_Cols=0;
  
  size_t line_size = 0;
  char *fi_line = NULL;

  if (getline(&fi_line, &line_size, fi_file)!=-1)
  {
	Number_of_Rows++;
    int item_digit_len;
	WeightMeasure item_weight;
    char *fi_line_scan = fi_line;
	//std::cout<<"reach 0"<<std::endl;
    while (sscanf(fi_line_scan, "%f%n", &item_weight, &item_digit_len) == 1)
	{
		//std::cout<<"item_weight:"<<item_weight<<",";
		//std::cout<<"item_digit_len:"<<item_digit_len<<",";
		fi_line_scan += item_digit_len;
		Number_of_Cols++;
	}
	//std::cout<<"reach 1"<<std::endl;
	while((getline(&fi_line, &line_size, fi_file)!=-1))
		Number_of_Rows++;
	//std::cout<<"reach 2"<<std::endl;
  }

  if (NULL != fi_line)
  {
    free(fi_line);
  }

  if (fclose(fi_file) == -1)
  {
    perror("fclose fi_file");
    exit(1);
  }  
  
  if ((Number_of_Rows<=0)||(Number_of_Cols<=0))
  {
	std::cout<<"The dataset format is incorrect. First Col or Row is not detected. Exit"<<std::endl;
	exit(1);
  }
  WGCStatic::min_item=static_cast<Itemset::item_t>(1);
  WGCStatic::max_item=static_cast<Itemset::item_t>(Number_of_Cols);
  WGCStatic::num_transactions = static_cast<Transactionset::transaction_t>(Number_of_Rows);
  TransOverItems=WGCStatic::num_transactions/(WGCStatic::max_item-WGCStatic::min_item+1.0);
  
  std::cout<<"TransOverItems: "<<TransOverItems<<std::endl;
  
  //Create graph object
  WG<WeightMeasure> graph(WGCStatic::num_transactions, WGCStatic::min_item, WGCStatic::max_item);
  //Create graph object done
  
  if ((fi_file = fopen(fi_fname.c_str(), "rb")) == NULL)
  {

    std::cout << "Dataset " << fi_fname << " does not exist. Exit.\n";

    exit(1);
  }

  std::cout << "Reading dataset file " << fi_fname << "\n";  

  line_size = 0;
  fi_line = NULL;

  uint64_t row_index=1;
  while (getline(&fi_line, &line_size, fi_file) != -1)
  {

    int item_digit_len;
	WeightMeasure weight;

    char *fi_line_scan = fi_line;
	uint64_t col_index=1;
    while (sscanf(fi_line_scan, "%f%n", &weight, &item_digit_len) == 1)
    {
      if (weight < -std::numeric_limits<WeightMeasure>::max() || weight > std::numeric_limits<WeightMeasure>::max())
      {

		
		std::cout << "weight out of range: " << weight << "\n";
		exit(1);
      }
	  else
	  {
		//std::cout<<weight<<",";
		if (abs_choice==1)
		{
			weight=fabs(weight);//option for abs
		}
		
		if ((row_index<=WGCStatic::num_transactions)&&((col_index)<=WGCStatic::max_item-WGCStatic::min_item+1))
			graph.set_weight(row_index,col_index,weight);
		else
		{
			std::cout << "Dataset format is wrong, possibly not a rectangle. Exit \n";
			exit(1);			
		}
			
	  }
      fi_line_scan += item_digit_len;
	  col_index++;
    }
	//std::cout<<std::endl;
	row_index++;
  }

  if (NULL != fi_line)
  {
    free(fi_line);
  }

  if (fclose(fi_file) == -1)
  {
    perror("fclose fi_file");
    exit(1);
  }
  
  	std::cout<<std::endl;
  
  /*
  for (Transactionset::transaction_t i=1; i<=WGCStatic::num_transactions; i++)
  {
	for (Itemset::item_t j=WGCStatic::min_item; j<=WGCStatic::max_item; j++)
	{
		std::cout<<graph.get_weight(i,j)<<",";
	}
	std::cout<<std::endl;
  }
  */
   
  std::cout<<std::endl;
  //std::cout<<graph;  
      
  
  //starting approximation algorithm
  BitDb edgeCover(WGCStatic::num_transactions, WGCStatic::min_item, WGCStatic::max_item);//To set and query whether an edge is covered or not.
  
  std::multimap<WeightMeasure, edgeinfo, weightcomp> edgeRank;//multimap structure ranking edges according to weight from large to small.
  
  
  CartesianProductDb wcdb; //wcdb saves all the dense components
  
  for (Transactionset::transaction_t i=1; i<=WGCStatic::num_transactions; i++)
  {
	for (Itemset::item_t j=WGCStatic::min_item; j<=WGCStatic::max_item; j++)
	{
		
		edgeinfo tmpinfo;
		tmpinfo.v_a=i;
		tmpinfo.v_b=j;
		edgeRank.insert(std::pair<WeightMeasure, edgeinfo>(graph.get_weight(i,j),tmpinfo));
		edgeCover.set_zero(i,j);//To ensure all zeroes at the beginning
	}
  }
  
  WeightMeasure MAXWEIGHT=edgeRank.begin()->first;
  WeightMeasure Cutoff=gamma*MAXWEIGHT;//The cutoff value of low weight edges that is not worth intializing a new cluster
  
 
  
  
  for (std::multimap<WeightMeasure, edgeinfo, weightcomp>::iterator edgeit=edgeRank.begin(); edgeit!=edgeRank.end(); edgeit++)
  {
	//std::cout<<"first edge weight:"<<edgeit->first<<std::endl;
	
	if (edgeit->first<Cutoff)
		break;
	
	if (edgeCover.exists(edgeit->second.v_a, edgeit->second.v_b))
		continue;
	
	//std::cout<<"cutoff value:"<<Cutoff<<std::endl;
	
	//mark the edge as selected
	edgeCover.insert(edgeit->second.v_a, edgeit->second.v_b);

	
	Transactionset V_Trans;
	Itemset V_Items;
	
	V_Trans.push_back(edgeit->second.v_a);
	V_Items.push_back(edgeit->second.v_b);

	std::set<Transactionset::transaction_t> SelectedVertices_Trans;
	std::set<Transactionset::transaction_t> SelectedVertices_Items;	
	SelectedVertices_Trans.insert(edgeit->second.v_a);
	SelectedVertices_Items.insert(edgeit->second.v_b);
	//std::cout<<"First vertex: "<<edgeit->second.v_a<<"; Second vertex: "<<edgeit->second.v_b<<";";
	
	
	WeightMeasure cluster_weight=edgeit->first;
	
	float density=((float)(cluster_weight))/(SelectedVertices_Trans.size()*SelectedVertices_Items.size());
	
	while(true)
	{
		Transactionset::transaction_t bestnode_a=-1;
		WeightMeasure Contribution_a=-1;
		
		for(Transactionset::transaction_t i=1; i<=WGCStatic::num_transactions; i++)
		{
			if (SelectedVertices_Trans.find(i)!=SelectedVertices_Trans.end())
				continue;
				
			WeightMeasure tmp_Weight=0;
			


			for(std::set<Itemset::item_t>::iterator Selected_it=SelectedVertices_Items.begin(); Selected_it!=SelectedVertices_Items.end(); Selected_it++)
			{
				tmp_Weight+=graph.get_weight(i, *Selected_it);
			}
			
			if (tmp_Weight>Contribution_a)
			{
				Contribution_a=tmp_Weight;
				bestnode_a=i;
			}
			
		}
		double C_v_a=((double)Contribution_a)/SelectedVertices_Items.size();//nomalize to average edge
		float lambda_a=std::max(1.0, C_para/pow((SelectedVertices_Trans.size()+tau_para+1),2));
		float alpha_a=1-1/((float)(lambda_a*(SelectedVertices_Trans.size()+tau_para+1)));	
	
		Itemset::item_t bestnode_b=-1;
		WeightMeasure Contribution_b=-1;
		
		for(Itemset::item_t i=WGCStatic::min_item; i<=WGCStatic::max_item; i++)
		{
			if (SelectedVertices_Items.find(i)!=SelectedVertices_Items.end())
				continue;
				
			WeightMeasure tmp_Weight=0;
			


			for(std::set<Transactionset::transaction_t>::iterator Selected_it=SelectedVertices_Trans.begin(); Selected_it!=SelectedVertices_Trans.end(); Selected_it++)
			{
				tmp_Weight+=graph.get_weight(*Selected_it, i);
		
			}
			if (tmp_Weight>Contribution_b)
			{
				Contribution_b=tmp_Weight;
				bestnode_b=i;
			}
			
		}
		double C_v_b=((double)Contribution_b)/SelectedVertices_Trans.size();//nomalize to average edge	
		float lambda_b=std::max(1.0, C_para/pow((SelectedVertices_Items.size()+tau_para+1),2));
		float alpha_b=1-1/((float)(lambda_b*(SelectedVertices_Items.size()+tau_para+1)));	
		
		bool Choose_A=false;
		if ((C_v_a>=alpha_a*density) && (C_v_b>=alpha_b*density))
		{
			if (C_v_a>C_v_b)
				Choose_A=true;
			else
				Choose_A=false;
				
//			std::cout<<"SelectedVertices_Trans.size()/(1.0*SelectedVertices_Items.size()): "<<SelectedVertices_Trans.size()/(1.0*SelectedVertices_Items.size())<<"; TransOverItems*MaxTransRatio: "<<TransOverItems*MaxTransRatio<<std::endl;
//			std::cout<<"SelectedVertices_Items.size()/(1.0*SelectedVertices_Trans.size()): "<<SelectedVertices_Items.size()/(1.0*SelectedVertices_Trans.size())<<"; MaxItemsRatio/TransOverItems: "<<MaxItemsRatio/TransOverItems<<std::endl;
			
//			if ((SelectedVertices_Trans.size()/(1.0*SelectedVertices_Items.size()))>TransOverItems*MaxTransRatio)
//				Choose_A=false;
//			else if ((SelectedVertices_Items.size()/(1.0*SelectedVertices_Trans.size()))>MaxItemsRatio/TransOverItems)
//				Choose_A=true;
				
//			std::cout<<"SelectedVertices_Trans.size()/(1.0*SelectedVertices_Items.size()): "<<SelectedVertices_Trans.size()/(1.0*SelectedVertices_Items.size())<<std::endl;
//			std::cout<<"SelectedVertices_Items.size()/(1.0*SelectedVertices_Trans.size()): "<<SelectedVertices_Items.size()/(1.0*SelectedVertices_Trans.size())<<std::endl;

//			if (SelectedVertices_Trans.size()*WGCStatic::num_transactions>SelectedVertices_Items.size()*(WGCStatic::max_item-WGCStatic::min_item+1.0))
//				Choose_A=false;
//			else if (SelectedVertices_Items.size()*(WGCStatic::max_item-WGCStatic::min_item+1.0)>SelectedVertices_Trans.size()*WGCStatic::num_transactions)
//				Choose_A=true;
				
//			if (SelectedVertices_Trans.size()>SelectedVertices_Items.size())
//				Choose_A=false;
//			else if (SelectedVertices_Items.size()>SelectedVertices_Trans.size())
//				Choose_A=true;
				
		}
		else if (C_v_a>=alpha_a*density)
		{
			Choose_A=true;
		}
		else if (C_v_b>=alpha_b*density)
		{
			Choose_A=false;
		}
		else//None of them are above the threshold
			break;
		

		//std::cout<<"c_v_C: "<<c_v_C<<"; alpha_n: "<<alpha_n<<"; density: "<<density<<"; bestnode: "<<bestnode<<"; Contribution: "<<Contribution<<std::endl;
		
		if(Choose_A)
		{
			V_Trans.push_back(bestnode_a);

			for(std::set<Itemset::item_t>::iterator Selected_it=SelectedVertices_Items.begin(); Selected_it!=SelectedVertices_Items.end(); Selected_it++)
				edgeCover.insert(bestnode_a, *Selected_it);
			
			SelectedVertices_Trans.insert(bestnode_a);
			cluster_weight+=Contribution_a;
		}
		else
		{
			V_Items.push_back(bestnode_b);

			for(std::set<Transactionset::transaction_t>::iterator Selected_it=SelectedVertices_Trans.begin(); Selected_it!=SelectedVertices_Trans.end(); Selected_it++)
				edgeCover.insert(*Selected_it, bestnode_b);
			
			SelectedVertices_Items.insert(bestnode_b);
			cluster_weight+=Contribution_b;
		
		}
		
		if (SelectedVertices_Trans.size()==WGCStatic::num_transactions)
			break;
		
		if (SelectedVertices_Items.size()==WGCStatic::max_item-WGCStatic::min_item)
			break;

		//Update Density !!!
		density=((float)(cluster_weight))/(SelectedVertices_Trans.size()*SelectedVertices_Items.size());
			
	}
	
	//SelectedVertices size must be no less than 2.
	CartesianProduct cp(V_Items, V_Trans);
    wcdb.push_back(cp);	
	
  }
  
  
  //std::cout<<std::endl;
  
  for (CartesianProductDb::iterator w_it=wcdb.begin(); w_it!=wcdb.end(); w_it++)
  {

	wcdb_outp<<*w_it;
	
	WeightMeasure cluster_weight=0;
	
	for (Transactionset::iterator Trans_it=w_it->transactionset.begin(); Trans_it!=w_it->transactionset.end(); Trans_it++)
	{	
		for(Itemset::iterator Item_it=w_it->itemset.begin(); Item_it!=w_it->itemset.end(); Item_it++)
		{
			//std::cout<<graph.get_weight(*Trans_it, *Item_it)<<" ";
			cluster_weight+=graph.get_weight(*Trans_it, *Item_it);//each edge weight has been calculated twice
		}
		//std::cout<<std::endl;
	}
	wcdb_outp<<";Density "<<(1.0*cluster_weight)/(w_it->transactionset.size()*w_it->itemset.size())<<std::endl;
	//std::cout<<std::endl;
  }
   
  wcdb_outp.close();
  //ending approximation algorithm


  return 0;

}
