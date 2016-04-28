#ifndef OPTREE_H
#define OPTREE_H

#include<vector>
#include<new>
#include<math.h>
#include<map>
#include <inttypes.h>

using namespace std;

static const double rho0_b = 0.5; // global variable prior stopping probability
static const int prec = 1; // precision of the OPTree's
static const double alpha0_b=0.5;
static const double rho0_b_c = 0.5;
static const double alpha0_b_c=0.5;


typedef uint32_t INDEX_TYPE;
typedef pair<INDEX_TYPE, INDEX_TYPE> DATA_TYPE;


double log_exp_x_plus_exp_y(double x, double y);

class GBT;

class OPTREE_NODE { // node for OPTree

  //friend class OPTREE;
	friend class GBT;
 private:

	INDEX_TYPE I0; // indicators for whether the node covers block 0 in each dimension
	INDEX_TYPE I1; // indicators for whether the node covers block 1 in each dimension
	INDEX_TYPE I2; // indicators for whether the node covers block 2 in each dimension

	INDEX_TYPE n_gt1_b[3][2]; // indicator for whether children in each dimension has at least 1 data point
	INDEX_TYPE n_gt2_b[3][2]; // indicator for whether children in each dimension has at least 2 data point

	unsigned int level; // the location of the first divisible 7
	unsigned int k_div2; // number of dimensions divisible but not intact
	unsigned int k_div3; // number of dimensions intact

	double log_rho; // stopping probability; this is meaningless for leaf nodes
	unsigned int n_x[2]; // number of data points from the two samples

	double logphi; // log expected conditional density of the data
	double logQ;



	DATA_TYPE xlast; // when n_x == 1, x is the LAST data point added to the node

	int samplelast;

	OPTREE_NODE * gbc[7]; // an array of pointers to its GBC


	GBT *gbt;




 public:
	OPTREE_NODE(GBT *gbt, INDEX_TYPE I0, INDEX_TYPE I1, INDEX_TYPE I2, unsigned int level, unsigned int k_div2, unsigned int k_div3) :
	  I0(I0), I1(I1), I2(I2), level(level), k_div2(k_div2),k_div3(k_div3),log_rho(log(rho0_b)), logphi(0), logQ(0), gbt(gbt)
	  {
	    n_x[0] = 0;
	    n_x[1] = 0;

	    for (int j = 0; j<3; j++) {
	      n_gt1_b[j][0] = 0;
	      n_gt1_b[j][1] = 0;
	      n_gt2_b[j][0] = 0;
	      n_gt2_b[j][1] = 0;
	    }


	  }; // constructor for OPTREE_NODE


	~OPTREE_NODE() {}; // empty desctructor

	bool is_root() {return !(~(I0 & I1 & I2));}

	double get_log_size() {
	  return log(pow(2,k_div2)*pow(3,k_div3));
	}

	void free_subGBT();

	// Leaf node is a node that has reached the support threshold.
	// Leaf nodes have no (left and right) children.

	bool data_inbound(DATA_TYPE x);
	   // check whether a data point is within the bound of the node


	int add_data(DATA_TYPE & x, int sample); // add a data point to the node

	pair<pair<int,int>, double> which_dim_to_split();
	void print_coupling_subtree(int offset, int width);


};


class GBT {

  friend class OPTREE_NODE;


 private:

  OPTREE_NODE *root;

  map<uint64_t, OPTREE_NODE *> gbc_map;
  //public:

  unsigned int k;

 public:

  GBT(unsigned int k) : k(k), root(new OPTREE_NODE(this, ~(~0 << k), ~(~0<<k), ~(~0 <<k), 0, 0, k )) {};
  ~GBT();
    //map<INDEX_TYPE, OPTREE_NODE *>::iterator iter;

  void add_data(DATA_TYPE & x, int sample) {root->add_data(x, sample);};
  void add_data(DATA_TYPE *x, unsigned n_x, int sample) {for(unsigned i=0; i<n_x; i++) root->add_data(x[i], sample);}
  void add_data(vector<DATA_TYPE> x, int sample) {for (unsigned i=0; i< x.size(); i++) {root->add_data(x[i], sample);}}

  unsigned get_size() {
    return gbc_map.size();
  }


  int update();
  double get_root_log_rho() {return root->log_rho;}
  double get_root_logphi() {return root->logphi;}
  double get_root_logQ() {return root->logQ;}

  void print_coupling_tree(int width) { root->print_coupling_subtree(0, width); };


};




#endif
