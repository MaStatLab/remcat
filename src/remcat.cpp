#include<iostream>
#include<iomanip>
using std::cout;
using std::cerr;
using std::endl;
using std::setw;

#include "optree.h"
#include "remcat.h"

#include <math.h>

int get_dim(INDEX_TYPE index, int d) { //get the value of dimension d of index
	  return (index >> (d + d - 2)) & ~(~(INDEX_TYPE)0 << 2);
}



COND_NODE::COND_NODE(COND_GBT *cond_gbt, INDEX_TYPE I0, INDEX_TYPE I1, uint32_t size, unsigned int level, unsigned int k_div) :
  I0(I0), I1(I1), size(size), level(level), k_div(k_div),loggamma(log(gamma0)), logphi(0),logpsi(0),logM0(0),logM1(0),logM(0),cond_gbt(cond_gbt)
{
  n_gt1_b0 = 0;
  n_gt2_b0 = 0;
  n_gt1_b1 = 0;
  n_gt2_b1 = 0;

  n_x[0]=n_x[1]=0;

  gbt_all = new GBT(cond_gbt->kg);

  switch (cond_gbt->mode) {

  case 1: logrho = log(cond_gbt->rho0); break;// constant initial rho
  case 2: logrho = log(1 - (1 - cond_gbt->rho0)/pow((level+1),2.0)); break; // 1-rho = (1-rho0)/level^2
  case 3: logrho = log(1 - (1 - cond_gbt->rho0)*exp((-1.0)*level * c)); break; // 1-rho = (1-rho0) e^{-c*level}
  case 4: logrho = log(1 - c * pow(1-cond_gbt->rho0, level+1)); break; // rho = 1- c*(1-rho0)^{level+1}
  case 5: logrho = log(cond_gbt->rho0) - level * c; break;// rho = rho0 * e^{-level}

  }

}; // constructor for OPTREE_NODE


bool COND_NODE::data_inbound(COV_TYPE &xe) {

  //uint k;
  //k = cond_gbt->kg;
  //cout << oct << I0 << "," << oct << I1 << "," << k << endl;
  //cout << oct << x << endl;
  return !( ~((~xe & I0) | (xe & I1)) & (~(~0 << cond_gbt->ke)) );
}

int COND_NODE::add_data(GENE_TYPE xg, COV_TYPE xe, int sample) {

  if (!data_inbound(xe)) {

    // cout << "Data out of bound!" << endl;
    return 1; // data out of bound; nothing happens; returns 1

  }

  else {

    //cout << "Try: " << oct << I0 << "," << oct << I1 << "," << level << endl;


    n_x[sample]++;
    // ng[sample][xg]++; // genotype x = 0, 1, or 2

    gbt_all->add_data(xg,sample);

    //uint k = cond_gbt->kg; // number of genotypes
    //DATA_TYPE x = (x1 >> k) & (~(~0 << cond_gbt->ke));

    if (n_x[0] + n_x[1] == 1) {
      xglast = xg;
      xelast = xe;
      samplelast = sample;

      n_gt1_b0 = (~xe) & I0 & (~(~0 << cond_gbt->ke));
      n_gt1_b1 = xe & I1 & (~(~0 << cond_gbt->ke))  ;

    }

    else { // first update counters for child nodes

      n_gt2_b0 |= (~xe) & n_gt1_b0;
      n_gt2_b0 &=(~(~0 << cond_gbt->ke));
      n_gt1_b0 |= (~xe) & I0;
      n_gt1_b0 &= (~(~0 << cond_gbt->ke));

      n_gt2_b1 |= xe & n_gt1_b1;
      n_gt2_b1 &= (~(~0 << cond_gbt->ke));
      n_gt1_b1 |= xe & I1;
      n_gt1_b1 &= (~(~0 << cond_gbt->ke));
    }

    if (n_x[0] + n_x[1] == 2) {

      if ( ! (level >=1 && level < cond_gbt->ke && ((((I0 & I1) >> (level-1)) & 3) == 3 )  ) ) { // GBC 0 has been added so we skip

	cond_gbt->gbc_map.insert(pair<uint64_t, COND_NODE *> ( ((uint64_t) I0 << 32) | ((uint64_t) I1) , this   ) );

     }

    }

    if ( level < cond_gbt->ke ) { // so this node has GBC


      if ( n_x[0] + n_x[1] == 2 ) {

	/* add the node (pointer) to the map
	 the key is the last 32 bits of I0 and last 32 bits of I1 joined
	*/

	gbc[0]= new COND_NODE( cond_gbt, I0, I1, size, level+1, k_div);
	gbc[1]= new COND_NODE( cond_gbt, I0, I1 & (~(1 << level)), size >> 1, level+1, k_div-1); // half size : size >> 1
	gbc[2]= new COND_NODE( cond_gbt, I0 & (~(1 << level)), I1, size >> 1, level+1, k_div-1);


	for(int i=0; i <= 2; i++) {

	  gbc[i]->add_data(xglast,xelast,samplelast);
	  gbc[i]->add_data(xg,xe,sample);
	}

      }

      if (n_x[0] + n_x[1] > 2) { // n_x >= 3 so the GBC must already exist; just need to update them

	for(int i=0; i <= 2; i++) gbc[i]->add_data(xg,xe,sample);
      }
    }

    return 0; //success
  }


}


int COND_GBT::update() {

  map<uint64_t, COND_NODE *>::iterator iter;
  uint i;
  INDEX_TYPE div_ind;
  uint64_t index;
  double m_stop,m_stop_psi;
  double m_divide,m_divide_psi;
  double m_divide_curr,m_divide_psi_curr;
  double loglambda0;

  COND_NODE *NODE_CURR;
  COND_NODE *CHILD_0;
  COND_NODE *CHILD_1;
  int n_x_c0;
  int n_x_c1;

  for ( iter = gbc_map.begin(); iter != gbc_map.end(); iter++ ) { // every node in this map will have at least 2 data points

    NODE_CURR = iter->second;
    //cout << hex << NODE_CURR->I0

    //cout << "OK " << NODE_CURR->size << endl;
    //cout << oct << NODE_CURR->I0 << "," << NODE_CURR->I1 << endl;
    /*
    cout << NODE_CURR->gbt_x1->root->n_x << endl;
    cout << NODE_CURR->gbt_x1->k << endl;
    cout << NODE_CURR->gbt_x1->gbc_map.size() << endl;
    */

    NODE_CURR->gbt_all->update();

    // cout << NODE_CURR->gbt_all->get_root_logphi() << "," << NODE_CURR->gbt_controls->get_root_logphi() << "," << NODE_CURR->gbt_cases->get_root_logphi() << endl;

    NODE_CURR->logM0 = NODE_CURR->gbt_all->get_root_logphi();

    NODE_CURR->logM = NODE_CURR->gbt_all->get_root_logQ();

    NODE_CURR->loggamma = NODE_CURR->gbt_all->get_root_log_rho();

    //cout << "OK2" << endl;

    //NODE_CURR->gbt_x1->update();
    //cout << "OK2" << endl;

    //NODE_CURR->gbt_x2->update();
    //NODE_CURR->gbt_x12->update();



    if (NODE_CURR->size == 1) { // cannot partition in Y

      NODE_CURR->logphi = NODE_CURR->logM;
      NODE_CURR->logpsi = log(gamma0) + NODE_CURR->logM0;
      NODE_CURR->logrho = 0;

    }


    else { // divisible node


      m_divide = 0;
      m_divide_psi = 0;

      div_ind = (NODE_CURR->I0 & NODE_CURR->I1);
      index = iter->first;

      if (NODE_CURR->k_div > 1) loglambda0 = (-1.0) * log(NODE_CURR->k_div);
      else loglambda0 = 0;


      for (i=0; i < ke; i++) { // for each dimension i


	m_divide_curr = 0;
	m_divide_psi_curr = 0;

	if ( (  div_ind >> i ) & 1 ) { // if the ith dimension is divisible

	  //cout << "OK2" << endl;


	  if ( ( (NODE_CURR->n_gt2_b0) >> i ) & 1 ) { /* if child 0 in dimension i has more than 2 data points

						      then it must be in the map				       */
	    // cout << "OK3" << endl;


	    CHILD_0 = gbc_map.find( index & ( ~((uint64_t) 1 << i))  )->second;
	    n_x_c0 = CHILD_0->n_x[0] + CHILD_0->n_x[1];

	    m_divide_curr += CHILD_0->logphi;
	    m_divide_psi_curr += CHILD_0->logpsi;


	  }

	  else if ( ( (NODE_CURR->n_gt1_b0) >> i ) & 1 ) { /* if child 0 in dimension i has 1 data point */

	    n_x_c0 = 1;

	    m_divide_curr += (-1.0) * kg * log(3); //kg = 1 Phi(A) = 1/(3^kg) for atomic sets
	    m_divide_psi_curr += (-1.0) * kg * log(3) - log(2); // Psi(A) = 1/(2*3^kg) for atomic sets

	  }

	  else { /* if child 0 has no data points*/

	    n_x_c0 = 0;

	    m_divide_psi_curr += log(gamma0);
	  }

	  n_x_c1 = NODE_CURR->n_x[0]+NODE_CURR->n_x[1] - n_x_c0;

	  if ( n_x_c1 > 1 ) { /* if child 1 in dimension i has more than 2 data points
						       then it must be in the map				       */

	    CHILD_1 = gbc_map.find( index & ( ~ ((uint64_t) 1 << (i+32)) ))->second;

	    m_divide_curr += CHILD_1->logphi;
	    m_divide_psi_curr += CHILD_1->logpsi;
	  }

	  else if (n_x_c1 == 1) { /* if child 0 in dimension i has 1 data point */

	    m_divide_curr += (-1.0) * kg * log(3);
	    m_divide_psi_curr += (-1.0) * kg * log(3) - log(2);

	  }

	  else { // n_x_c1 == 0

	    m_divide_psi_curr += log(gamma0);
	  }


	  m_divide_curr += loglambda0;
	  m_divide_psi_curr += loglambda0;

	  if (m_divide==0) {
	    m_divide = m_divide_curr;
	    m_divide_psi = m_divide_psi_curr;
	  }
	  else {
	    m_divide = log_exp_x_plus_exp_y (m_divide, m_divide_curr);
	    m_divide_psi = log_exp_x_plus_exp_y (m_divide_psi, m_divide_psi_curr);
	  }

	}
      }



      m_divide += log(1-exp(NODE_CURR->logrho));
      m_divide_psi += log(1-exp(NODE_CURR->logrho));

      NODE_CURR->m_divide = m_divide; // for later figuring out which dimension to partition
      NODE_CURR->m_divide_psi = m_divide_psi;

      m_stop = NODE_CURR->logrho + NODE_CURR->logM;
      m_stop_psi = NODE_CURR->logrho + log(gamma0) + NODE_CURR->logM0;

      NODE_CURR->logphi = log_exp_x_plus_exp_y(m_divide,m_stop);
      NODE_CURR->logrho = m_stop - NODE_CURR->logphi;


      NODE_CURR->logpsi = log_exp_x_plus_exp_y(m_divide_psi,m_stop_psi);

      // cout << oct << NODE_CURR->I0 << "," << NODE_CURR->I1 << ",," << dec << NODE_CURR->size << ",," << NODE_CURR->logM0 << "," << NODE_CURR->logM1 << "," << NODE_CURR->logM << ",," << NODE_CURR->logphi << "," << NODE_CURR->logpsi << endl;

    }

    //cout << oct << NODE_CURR->I0 << "," << oct << NODE_CURR->I1 << "," << dec << NODE_CURR->n_x << "," << NODE_CURR->logrho << endl;
    //cout << NODE_CURR->size << "," << NODE_CURR->logphi << "," << NODE_CURR->gbt_x1->get_root_logphi() << endl << endl;


  }

  return 0;
}



void COND_NODE::free_subGBT() {
	if (level < cond_gbt->ke && n_x[0] + n_x[1] >= 2) {
		gbc[0]->free_subGBT();
		gbc[1]->free_subGBT();
		gbc[2]->free_subGBT();
	}

	delete gbt_all;

	delete this;
}


COND_GBT::~COND_GBT() {
  root->free_subGBT();
}





// double COND_NODE::get_split_prob(int i) { //get the splitting probability of dimension i

//   uint64_t index = get_index();
//   double m_divide_curr = 0;
//   COND_NODE * CHILD_0;
//   COND_NODE * CHILD_1;

//   int n_x_0;
//   int n_x_1;
//   double loglambda0 = (-1.0) * log(k_div);

//   INDEX_TYPE div_ind = (I0 & I1) & ~(~(uint64_t) 0 << cond_gbt->ke);

//   if (k_div == 0) {
//     cerr << "Error: Node indivisible! Can't get splitting probibilities!" << endl;
//     exit(1);
//   }


//   if ( (  div_ind >> i ) & 1 ) { // if the ith dimension is divisible


//     if (k_div > 1) {


//       if ( ( (n_gt2_b0) >> i ) & 1 ) { /* if child 0 in dimension i has more than 2 data points
// 						   then it must be in the map				       */

// 	CHILD_0 = cond_gbt->gbc_map.find( index & ( ~((uint64_t) 1 << i))  )->second;
// 	n_x_0 = CHILD_0->n_x;
// 	m_divide_curr += CHILD_0->logphi;

//       }

//       else if ( ( (n_gt1_b0) >> i ) & 1 ) { /* if child 0 in dimension i has 1 data point */

// 	n_x_0 = 1;
// 	m_divide_curr += (-1.0) * log(2) * (cond_gbt->kg);

//       }

//       else { /* if child 0 has no data points*/

// 	n_x_0 = 0;

//       }

//       n_x_1 = n_x - n_x_0;

//       if ( n_x_1 > 1 ) { /* if child 1 in dimension i has more than 2 data points
// 						   then it must be in the map				       */

// 	CHILD_1 = cond_gbt->gbc_map.find( index & ( ~ ((uint64_t) 1 << (i+32)) ))->second;
// 	m_divide_curr += CHILD_1->logphi;

//       }

//       else if ( n_x_1 == 1) { /* if child 0 in dimension i has 1 data point */

// 	m_divide_curr += (-1.0) * log(2) * (cond_gbt->kg);

//       }

//       else {}


//       m_divide_curr += loglambda0 + log(1-rho0);


//       return exp(m_divide_curr - m_divide);
//     }

//     else return 1; // if k_div == 1 then i is the only dimension to split
//   }

//   else {
//     return 0;
//   } // if dimension i is indivisible then return 0;


// }

uint64_t COND_NODE::get_index() {

  return ((uint64_t) I0 << 32) | ((uint64_t) I1);

}


// pair<int, double> COND_NODE::get_map_split_dim() {

//   int curr_max = 0;
//   double curr_max_split_prob = -1;
//   double curr_split_prob;
//   uint i;

//   if ( k_div == 0 ) {
//     cerr << "Node not divisible" << endl;
//     exit(1);
//   }


//   for ( i = 0; i < cond_gbt->ke; i++ ) {
//     if ( ((I0 & I1) >> i) & 1 ) {
//       curr_split_prob = get_split_prob(i);
//       //cout << curr_split_prob <<endl;

//       if (curr_split_prob > curr_max_split_prob) {

// 	curr_max = i;
// 	curr_max_split_prob = curr_split_prob;
//       }
//     }
//   }

//   return make_pair(curr_max, curr_max_split_prob);
// }


// void COND_GBT::print_structure ( COND_NODE *node, int l, int full) {

//   map<uint64_t, COND_NODE *>::iterator iter;
//   COND_NODE * CHILD_0 = 0;
//   COND_NODE * CHILD_1 = 0;
//   int i;

//   if ( !node ) {

//     padding ( '\t', l );

//     puts ( "~" );

//   }

//   else {
//     //cout << "Parent=" << node << endl;
//     uint64_t index = node->get_index();

//     if (node->k_div >= 1) {
//       pair<int, double> split_dim = node->get_map_split_dim();
//       i = split_dim.first;


//       iter = gbc_map.find( index & ( ~((uint64_t) 1 << i))  );
//       if (iter != gbc_map.end() ) {
// 	CHILD_0 = iter->second;

//       }
//       else {
// 	CHILD_0 = 0;
//       }

//       iter = gbc_map.find( index & ( ~ ((uint64_t) 1 << (i+32)) ));
//       if (iter != gbc_map.end() ) {
// 	CHILD_1=iter->second;

//       }

//       else {

// 	CHILD_1 = 0;

//       }

//       if ( full ||  ( node->logrho < log(0.5)) ) {
// 	print_structure ( CHILD_1, l + 1 , full);
//       }

//       padding ( '\t', l );

//       printf ( "%1.2f, (%d, %1.2f)\n", exp(node->logrho), split_dim.first, split_dim.second );

//       if ( full ||  ( node->logrho < log(0.5) ) ) {

// 	print_structure ( CHILD_0, l + 1 , full);


//       }
//     }

//     else { // if the node is not divisible

//       padding ( '\t', l );
//       printf ( "n_x=%d, Leaf, %1.2f \n", node->n_x, exp(node->logrho));

//     }

//   }

// }


// void COND_GBT::print_map_tree(int full) {

//   print_structure(root, 0, full);


// }





