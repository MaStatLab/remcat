#ifndef COND_TREE_H
#define COND_TREE_H

#include<cstdlib>
#include<stdio.h>

#include<vector>
#include<new>
#include<math.h>
#include<map>
#include<utility>
#include<inttypes.h>

#include "optree.h"

using namespace std;

static const double gamma0 = 0.5;
static const double alpha0 = 0.5;
static const double c = 1;

typedef uint32_t COV_TYPE;
typedef DATA_TYPE GENE_TYPE;

class COND_GBT;
class COND_NODE;


class COND_NODE {

  friend class OPTREE_NODE;
  friend class GBT;
  friend class COND_GBT;

 public:

  INDEX_TYPE I0; // indicators for whether the node covers block 0 in each dimension
  INDEX_TYPE I1; // indicators for whether the node covers block 1 in each dimension
  INDEX_TYPE n_gt1_b0; // indicator for whether child 1 in each dimension has at least 1 data point
  INDEX_TYPE n_gt2_b0; // indicator for whether child 1 in each dimension has at least 2 data point
  INDEX_TYPE n_gt1_b1; // indicator for whether child 1 in each dimension has at least 1 data point
  INDEX_TYPE n_gt2_b1; // indicator for whether child 1 in each dimension has at least 2 data point
  uint32_t size; // the size of the node i.e. number of blocks

  unsigned int level; // the location of the first divisible 7
  unsigned int k_div; // number of dimensions available for division; size = 2^(k_div)

  double logrho; // stopping probability; this is meaningless for leaf nodes
  double loggamma;

  //uint ng[2][3]; // counts of three genotypes

  int n_x[2]; // number of data points

  double logphi, logpsi; // log expected conditional density of the data
  double logM0, logM1, logM;

  // double m_stop;
  double m_divide,m_divide_psi;

  GENE_TYPE xglast;
  COV_TYPE xelast; // when n_x == 1, x is the LAST data point added to the node
  int samplelast;

  //DATA_TYPE xlast2; // when n_x =2, x is the second LAST data point added to the node
  // if n_x >2 then this value is not useful since the node has children

  COND_NODE * gbc[3]; // an array of pointers to its GBC

  COND_GBT *cond_gbt;
  //GBT *gbt_x1;

  GBT *gbt_all;

  COND_NODE(COND_GBT *cond_gbt, INDEX_TYPE I0, INDEX_TYPE I1, uint32_t size, unsigned int level, unsigned int k_div);
  ~COND_NODE() {};


  bool is_root() {return !(~(I0 & I1));}
  void free_subGBT();
  bool data_inbound(COV_TYPE & xe); // check whether a data point is within the bound of the node

  int add_data(GENE_TYPE xg, COV_TYPE xe, int sample); // add a data point to the node


  uint64_t get_index();


};


class COND_GBT{

  friend class COND_NODE;
  friend class OPTREE_NODE;
  friend class GBT;

 public:
  //private:

  COND_NODE *root;

  map<uint64_t, COND_NODE *> gbc_map;
  //public:

  int kg, ke;
  double rho0;
  int mode;
  // kg is number of genotypes $G$
  // ke is number of covariates $E$

 COND_GBT(int kg, int ke, double rho0, int mode): kg(kg), ke(ke), rho0(rho0), mode(mode)
  {
    root = new COND_NODE(this, ~0, ~0, 1 << ke, 0, ke);
  };
  ~COND_GBT();

  void add_data(GENE_TYPE &xg, COV_TYPE &xe, int &sample) {root->add_data(xg,xe,sample);};
  void add_data(GENE_TYPE *xg, COV_TYPE *xe, int *sample, unsigned n_x) {for(int i=0; i<n_x; i++) root->add_data(xg[i],xe[i],sample[i]);}

  unsigned get_size() {
    return gbc_map.size();
  }


  int update();

  double get_root_logrho() { return root->logrho;}
  double get_root_logpsi() { return root->logpsi;}
  double get_root_loggamma() { return root->loggamma;}
  double get_root_logphi() { return root->logphi; }

};





#endif
