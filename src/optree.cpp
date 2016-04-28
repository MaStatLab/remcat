/* 
	Class functions for OPTREE_NODE and GBT classes
	Li Ma 12.22.2009
*/

#include<iostream>
#include<iomanip>
#include<limits>
using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::numeric_limits;

#include "optree.h"
#include <math.h>
#include <float.h>


double log_exp_x_plus_exp_y(double x, double y) {

    double result;
    if ( x - y >= 100 ) result = x;
    else if ( x - y <= -100 ) result = y;
    else {
      if (x > y) {
      result = y + log( 1 + exp(x-y) );
      }
      else result = x + log( 1 + exp(y-x) );
    }
    
    return result;
};



int OPTREE_NODE::add_data(DATA_TYPE & x, int sample) {
 

  if (!data_inbound(x)) {

    return 1; // data out of bound; nothing happens; returns 1
    
  }
  
  else {
    
    n_x[sample]++;

    if (n_x[0] + n_x[1] == 1) {
      xlast=x;
      samplelast=sample;
    }

    if (n_x[sample] == 1) {
	
	n_gt1_b[0][sample] = x.first & I0 & ~(~0 << gbt->k);
	n_gt1_b[1][sample] = (~x.first) & (~x.second) & I1 & ~(~0 << gbt->k);
	n_gt1_b[2][sample] = x.second & I2 & ~(~0 << gbt->k); 
	
    }

    else { // first update counters for child nodes
	
	n_gt2_b[0][sample] |= x.first & n_gt1_b[0][sample];
	n_gt1_b[0][sample] |= x.first & I0 & ~(~0 << gbt->k);

	n_gt2_b[2][sample] |= x.second & n_gt1_b[2][sample];
	n_gt1_b[2][sample] |= x.second & I2 & ~(~0 << gbt->k);

	n_gt2_b[1][sample] |= (~x.first) & (~x.second) & n_gt1_b[1][sample]; 
	n_gt1_b[1][sample] |= (~x.first) & (~x.second) & I1 & ~(~0 << gbt->k);

      
    }

    if (n_x[0] + n_x[1] == 2) {

      if ( ! (level >=1 && level < gbt->k && ((((I0 & I1 & I2) >> (level-1)) & 3) == 3 )  ) ) { // GBC 0 has been added so we skip
	//cout << (((uint64_t) I0 << 42) | ((uint64_t) I1 << 21) | ((uint64_t) I2)) << endl;

	gbt->gbc_map.insert(pair<uint64_t, OPTREE_NODE *> 
			    ( ((uint64_t) I0 << 42) | ( (uint64_t) I1 << 21) | ((uint64_t) I2), this  ) ); 
	
      }
    }


    if ( level < gbt->k ) { // so this node has GBC
 

      if ((n_x[0] + n_x[1] == 2)) {
	 
	/* add the node (pointer) to the map
	 the key is the last 32 bits of I0 and last 32 bits of I1 joined 
	*/
	   
	gbc[0]= new OPTREE_NODE( gbt, I0, I1, I2, level+1, k_div2, k_div3);
	gbc[1]= new OPTREE_NODE( gbt, I0, I1 & (~(1 << level)), I2 & (~(1 << level)), level+1, k_div2, k_div3-1);
	gbc[2]= new OPTREE_NODE( gbt, I0 & (~(1 << level)), I1, I2 & (~(1 << level)), level+1, k_div2, k_div3-1);
	gbc[4]= new OPTREE_NODE( gbt, I0 & (~(1 << level)), I1 & (~(1 << level)), I2, level+1, k_div2, k_div3-1);
	gbc[3]= new OPTREE_NODE( gbt, I0, I1, I2 & (~(1 << level)), level+1, k_div2+1, k_div3-1);
	gbc[5]= new OPTREE_NODE( gbt, I0, I1 & (~(1 << level)), I2, level+1, k_div2+1, k_div3-1);
	gbc[6]= new OPTREE_NODE( gbt, I0 & (~(1 << level)), I1, I2, level+1, k_div2+1, k_div3-1);
	

	for(int i=0; i < 7; i++) {

	  gbc[i]->add_data(xlast,samplelast);
	  gbc[i]->add_data(x,sample);
	}
	
      }
  
      if (n_x[0] + n_x[1] > 2) { // n_x >= 3 so the GBC must already exist; just need to update them

	for(int i=0; i < 7; i++) gbc[i]->add_data(x,sample);
      }     
    }
        
    return 0; //success 
  }
}




void OPTREE_NODE::free_subGBT() {
	if (level < gbt->k && n_x[0]+n_x[1] >= 2) {
	  for (int i = 0; i < 7; i++){
	    gbc[i]->free_subGBT();
	  }
	}

	delete this;
}


GBT::~GBT() {

  root->free_subGBT();
  
};

bool OPTREE_NODE::data_inbound(DATA_TYPE x) {

  return !( ~((x.first & I0) | ((~x.first) & (~x.second) & I1) | (x.second & I2)) & (~(~0 << gbt->k)) )  ;

}; // check whether a data point is within the bound of the node



pair<pair<int,int>, double> OPTREE_NODE::which_dim_to_split() {
  
  

  if (k_div2 == 0 && k_div3 == 0) {
    cerr << "Error: No dim to split! This should not happen!" << endl;
    exit(1);
  }

  else if (n_x[0] + n_x[1] <= 1) { // if the node has no more than 1 data point
    
    // This should not happen by design
    cout << "Should not happen by design!" << endl;
    return make_pair(make_pair(-1,-1),0); // does not split such a node 
  }
      
  else {

   
    uint64_t index = ((uint64_t) I0 << (21*2) ) | ((uint64_t) I1 << 21) | (uint64_t) I2;
    // the index of the node in the map
    
  
    OPTREE_NODE *CHILD_0;
    OPTREE_NODE *CHILD_1;

    INDEX_TYPE div_ind = (I0 & I1) | (I1 & I2) | (I2 & I0);
    INDEX_TYPE div3_ind = I0 & I1 & I2;
    INDEX_TYPE div2_ind = div_ind & (~div3_ind);
    int n_x_0;
    int n_x_1;

    int n_x_0_1, n_x_1_1, n_x_0_2, n_x_1_2;


    double loglambda0;
    double alpha0_b_c_r;

    double m_divide = 0;
    double m_stop;
    double m_divide_curr;
    int i, j, j1, j2;

    double curr_max_m_divide = -DBL_MAX;
    double max_lambda;
    
    int max_divide_dim; // the current best dimension to cut
    int max_divide_way; // the current best way to cut the dimension (useful for intact dimensions)
    
    INDEX_TYPE n_gt1_b[3]; // indicator for whether child 1 in each dimension has at least 1 data point
    INDEX_TYPE n_gt2_b[3]; // indicator for whether child 1 in each dimension has at least 2 data point


    if (k_div2 == 1 &&  k_div3 == 0)  loglambda0 = 0;
    else loglambda0 = (-1.0) * log(k_div2 + 3*k_div3);
    
    for (j = 0; j < 3; j++) {
	 
      n_gt1_b[j] = this->n_gt1_b[j][0] | this->n_gt1_b[j][1];
      n_gt2_b[j] = this->n_gt2_b[j][0] | this->n_gt2_b[j][1] | (this->n_gt1_b[j][0] & this->n_gt1_b[j][1]);

    }       
  

    for (i=0; i< gbt->k; i++) { //dimension to cut

      if ( (  div3_ind >> i ) & 1 ) { // if the ith dimension is intact
	
	alpha0_b_c_r = 2*alpha0_b_c;

	for ( j = 0; j < 3 ; j++ ) {

	  m_divide_curr = 0;
	      
	  j1 = (j+1) % 3;
	  j2 = (j+2) % 3;


	  if ( ( (n_gt2_b[j]) >> i ) & 1 ) { /* if child 0 in dimension i has more than 2 data points
							 then it must be in the map				       */

	    CHILD_0 = gbt->gbc_map.find( index & ( ~((uint64_t) 1 << (i+ 21*(2-j1)) | (uint64_t) 1 << (i +21 *(2-j2))))  )->second;

	    n_x_0 = CHILD_0->n_x[0] + CHILD_0->n_x[1];
	    n_x_0_1 = CHILD_0->n_x[0];
	    n_x_0_2 = CHILD_0->n_x[1];
	    
	    m_divide_curr += CHILD_0->logQ;
	      
	  }

	  else {
	      
	    if  (( (this->n_gt1_b[j][0]) >> i ) & 1 ) {
		
	      n_x_0 = 1;
	      n_x_0_1 = 1;
	      n_x_0_2 = 0;
	      
	      m_divide_curr += (-1.0) * ( this->get_log_size() - log(3) );
	      
	      }
	    
	    else if (( (this->n_gt1_b[j][1]) >> i ) & 1 ) {
	      
	      n_x_0 = 1;
	      n_x_0_1 = 0;
	      n_x_0_2 = 1;
	      
	      m_divide_curr += (-1.0) * ( this->get_log_size() - log(3) );
	    }

	    else {
	      
	      n_x_0 = 0;
	      n_x_0_1 = 0;
	      n_x_0_2 = 0;
	    }
	    
	  }

	  n_x_1 = n_x[0] + n_x[1] - n_x_0;
	  n_x_1_1 = n_x[0] - n_x_0_1;
	  n_x_1_2 = n_x[1] - n_x_0_2;

	  if (n_x_1 > 1) {

	    CHILD_1 = gbt->gbc_map.find( index & ~((uint64_t) 1 << (i+21*(2-j))) )->second;
	    
	    m_divide_curr += CHILD_1->logQ;
	      
	  }

	  else if (n_x_1 == 1) {
	           
	    m_divide_curr += (-1.0) * ( this->get_log_size() - log(3) +log(2) );
	
	  }

	  m_divide_curr += loglambda0 
	      
	    + lgamma(n_x_0_1+alpha0_b_c) + lgamma(n_x_1_1+alpha0_b_c_r) - lgamma(n_x_0_1 + n_x_1_1 + alpha0_b_c + alpha0_b_c_r)
	    - ( lgamma(alpha0_b_c) + lgamma(alpha0_b_c_r) - lgamma( alpha0_b_c + alpha0_b_c_r ) )
	    
	    + lgamma(n_x_0_2+alpha0_b_c) + lgamma(n_x_1_2+alpha0_b_c_r) - lgamma(n_x_0_2 + n_x_1_2 + alpha0_b_c + alpha0_b_c_r)
	    - ( lgamma(alpha0_b_c) + lgamma(alpha0_b_c_r) - lgamma( alpha0_b_c + alpha0_b_c_r) );

	  
	  if (m_divide == 0) {
	    m_divide = m_divide_curr;
	  }
	  else { 
	    m_divide = log_exp_x_plus_exp_y (m_divide, m_divide_curr);
	  }
	  
	  if (m_divide_curr > curr_max_m_divide) {
	    
	    curr_max_m_divide = m_divide_curr;
	    max_divide_dim = i + 1;
	    max_divide_way = j + 1;
	    
	  }
	  
	}
	
      }

      
      if ( (div2_ind >> i) & 1 ) {
	    
	alpha0_b_c_r = alpha0_b_c;
	    
	m_divide_curr = 0;
	
	if ((I0 >> i) & 1) j=0; // block 0 is included in node i
	else j=1; // block 0 is not included in node i
	
	j1 = (j+1) % 3;
	j2 = (j+2) % 3;

	if ( ( (n_gt2_b[j]) >> i ) & 1 ) { /* if child 0 in dimension i has more than 2 data points
							 then it must be in the map				       */

	  CHILD_0 = gbt->gbc_map.find( index & ( ~((uint64_t) 1 << (i+ 21*(2-j1)) | (uint64_t) 1 << (i +21 *(2-j2))))  )->second;

	  n_x_0 = CHILD_0->n_x[0] + CHILD_0->n_x[1];
	  n_x_0_1 = CHILD_0->n_x[0];
	  n_x_0_2 = CHILD_0->n_x[1];
	  
	  m_divide_curr += CHILD_0->logQ;
	  
	}

	else {
	  
	  if  (( (this->n_gt1_b[j][0]) >> i ) & 1 ) {
		
	    n_x_0 = 1;
	    n_x_0_1 = 1;
	    n_x_0_2 = 0;
	    
	    m_divide_curr += (-1.0) * ( this->get_log_size() - log(2));
	    
	  }
	  
	  else if (( (this->n_gt1_b[j][1]) >> i ) & 1 ) {
	      
	    n_x_0 = 1;
	    n_x_0_1 = 0;
	    n_x_0_2 = 1;
	    
	    m_divide_curr += (-1.0) * ( this->get_log_size() - log(2));
	  }
	  
	  else {
	    
	    n_x_0 = 0;
	    n_x_0_1 = 0;
	    n_x_0_2 = 0;
	  }
	  
	}
	
	n_x_1 = n_x[0] + n_x[1] - n_x_0;
	n_x_1_1 = n_x[0] - n_x_0_1;
	n_x_1_2 = n_x[1] - n_x_0_2;

	if (n_x_1 > 1) {
	  
	  CHILD_1 = gbt->gbc_map.find( index & ~((uint64_t) 1 << (i+21*(2-j))) )->second;

	  m_divide_curr += CHILD_1->logQ;	  
	}

	else if (n_x_1 == 1) {

	  m_divide_curr += (-1.0) * ( this->get_log_size() - log(2) );
	}

	m_divide_curr += loglambda0 
	      
	  + lgamma(n_x_0_1+alpha0_b_c) + lgamma(n_x_1_1+alpha0_b_c_r) - lgamma(n_x_0_1 + n_x_1_1 + alpha0_b_c + alpha0_b_c_r)
	  - ( lgamma(alpha0_b_c) + lgamma(alpha0_b_c_r) - lgamma( alpha0_b_c + alpha0_b_c_r ) )
	  
	  + lgamma(n_x_0_2+alpha0_b_c) + lgamma(n_x_1_2+alpha0_b_c_r) - lgamma(n_x_0_2 + n_x_1_2 + alpha0_b_c + alpha0_b_c_r)
	  - ( lgamma(alpha0_b_c) + lgamma(alpha0_b_c_r) - lgamma( alpha0_b_c + alpha0_b_c_r) );
	
	
	if (m_divide==0) {
	  m_divide = m_divide_curr;
	}
	else { 
	  m_divide = log_exp_x_plus_exp_y (m_divide, m_divide_curr);
	}
	
	if (m_divide_curr > curr_max_m_divide) {
	  
	  curr_max_m_divide = m_divide_curr;
	  max_divide_dim = i + 1;
	  max_divide_way = j + 1;
	    
	}
	
      }


      max_lambda = exp(curr_max_m_divide - m_divide);
    }

    return make_pair(make_pair(max_divide_dim,max_divide_way),max_lambda);
  }

}
  


void padding ( char ch, int n ) {

  int i;

  for ( i = 0; i < n; i++ )

    putchar ( ch );
}

 

void OPTREE_NODE::print_coupling_subtree(int offset, int width) {
  
  int i, j1, j2, n_x_0, n_x_1, n_x_0_1, n_x_0_2, n_x_1_1, n_x_1_2;
  
  INDEX_TYPE div_ind = (I0 & I1) | (I1 & I2) | (I2 & I0);
  INDEX_TYPE div3_ind = I0 & I1 & I2;
  INDEX_TYPE div2_ind = div_ind & (~div3_ind);
 
  uint64_t print_ID = 0;
  uint64_t print_ID_left = 0;
  uint64_t print_ID_right = 0;

  OPTREE_NODE *CHILD_0;
  OPTREE_NODE *CHILD_1;


  
  //for (i=gbt->k - 1; i >= 0; i--) {
    
  for (i=0; i < gbt->k ; i++) {
    
    print_ID = (print_ID << 3) | ((I0 >> i) & 1) | (((I1 >> i) & 1) << 1) + (((I2 >> i) & 1) << 2);

  }

  if (k_div2 + k_div3 >= 1) { // node divisible

    if (log_rho >= log(0.5) || n_x[0] == 0 || n_x[1] == 0) {
      padding ( '\t', offset);
      printf ( "ID: %o\n", print_ID);
      padding ( '\t', offset);
      if (log_rho >= log(0.5)) {
	printf ( "Coupled: logrho=%1.2f\n", (-1.0)*log_rho/log(10));
      }
      else printf ( "Not coupled: logrho=%1.2f\n", (-1.0)*log_rho/log(10));
      
      padding ( '\t', offset);
      printf ("n0=%d, n1=%d\n", n_x[0], n_x[1]);
      padding ( '\t', offset);
      //printf("n1/N=%1.2f\%\n", (double) n_x[1] / (n_x[0]+n_x[1]) * 100);
      printf("(n1/N1)/(n0/N0)=%1.2f\n", (double) n_x[1] / gbt->root->n_x[1] / n_x[0] * gbt->root->n_x[0] );

    }

    else {
      
      pair<pair<int,int>, double> split_dim= which_dim_to_split();
      int i = split_dim.first.first - 1;
      int j = split_dim.first.second - 1;
      double lambda = split_dim.second;
      uint64_t index = ((uint64_t) I0 << 21*2) | ((uint64_t) I1 << 21) | ((uint64_t) I2);



      
      //if ( ( (div3_ind >> i) & 1 ) | ( (div2_ind >> i) & 1 ) ) { 
      if ( (div_ind >> i) & 1 ) { 
	
	j1 = (j+1) % 3;
	j2 = (j+2) % 3;
	  
	print_ID_left = print_ID & ~((uint64_t) 1 << 3*(gbt->k -1 - i)+j1) & ~((uint64_t) 1 << 3*(gbt->k - 1 - i) +j2);
	print_ID_right = print_ID & ~((uint64_t) 1 << 3*(gbt->k - 1 - i)+j);


// 	if ( ( n_gt2_b[j][1] >> i ) & 1 ) { /* if child 0 in dimension i has more than 2 case data points
// 					       then it must be in the map				       */

	if ( (( n_gt2_b[j][0] | n_gt2_b[j][1] | (n_gt1_b[j][0] & n_gt1_b[j][1] )) >> i) & 1)  { // If Child 0 has at least 2 data points

	  CHILD_0=gbt->gbc_map.find( index & ( ~((uint64_t) 1 << (i+ 21*(2-j1)) | (uint64_t) 1 << i +21 *(2-j2))) )->second;
	  CHILD_0->print_coupling_subtree(offset + width, width);
	  

	  n_x_0 = CHILD_0->n_x[0] + CHILD_0->n_x[1];
	  n_x_0_1 = CHILD_0->n_x[0];
	  n_x_0_2 = CHILD_0->n_x[1];

	}

	else {

	  if  (( (this->n_gt1_b[j][0]) >> i ) & 1 ) {
		
	    n_x_0 = 1;
	    n_x_0_1 = 1;
	    n_x_0_2 = 0;
    
	  }
	    
	  else if (( (this->n_gt1_b[j][1]) >> i ) & 1 ) {
	      
	    n_x_0 = 1;
	    n_x_0_1 = 0;
	    n_x_0_2 = 1;
	  }

	  else {
	      
	    n_x_0 = 0;
	    n_x_0_1 = 0;
	    n_x_0_2 = 0;
	  }

	  
	  padding ( '\t', offset+width);
	  printf ( "ID: %o\n", print_ID_left);
	  padding ( '\t', offset+width);
	  printf ( "Stop: logrho=%1.2f\n", (-1.0)*log(rho0_b_c)/log(10) );
	  padding ( '\t', offset+width);
	  printf ( "n0=%d, n1 =%d\n", n_x_0_1, n_x_0_2);
	  //padding ( '\t', offset+width);
	  //printf("(n1/N1)/(n0/N0)=Inf\n");
	}

	padding( '\t', offset);
	printf( "ID: %o\n", print_ID);
	padding( '\t', offset);
	printf( "not stop: logrho=%1.2f\n", (-1.0)*log_rho/log(10));
	padding ( '\t', offset);
	printf ( "Cut: %d, Way: %d, lambda = %1.2f\n", i+1, j+1, lambda);
	padding( '\t', offset);
	printf ("n0=%d, n1=%d\n", n_x[0], n_x[1]);
	padding ( '\t', offset);
	//printf("n1/N=%1.2f\%\n", (double) n_x[1] / (n_x[0]+n_x[1]) * 100);
	printf("(n1/N1)/(n0/N0)=%1.2f\n", (double) n_x[1] / gbt->root->n_x[1] / n_x[0] * gbt->root->n_x[0] );
	


	n_x_1 = n_x[0] + n_x[1] - n_x_0; // number of data points in right child
	n_x_1_1 = n_x[0] - n_x_0_1;
	n_x_1_2 = n_x[1] - n_x_0_2;


	if (n_x_1 >= 2) { 
		
	  CHILD_1 = gbt->gbc_map.find( index & ~((uint64_t) 1 << (i+21*(2-j))) )->second;
	  CHILD_1->print_coupling_subtree(offset + width, width);
	  
	}
	
	else {

	  padding ( '\t', offset+width);
	  printf ( "ID: %o\n", print_ID_right);
	  padding ( '\t', offset+width);
	  printf ( "Stop: logrho=%1.2f\n", (-1.0)*log(rho0_b_c)/log(10) );
	  padding ( '\t', offset+width);
	  printf ( "n0=%d, n1 =%d\n",n_x_1_1, n_x_1_2);
	  //padding ( '\t', offset+width);
	  //printf("n1/N=100\%\n");
	  //printf("(n1/N1)/(n0/N0)=Inf\n");
	  
	    
	}
	
      }
    }

  }

  else {

    padding ( '\t', offset);
    printf ( "ID: %o\n", print_ID);
    padding ( '\t', offset);
    printf ( "Stop: Leaf!\n");
    padding ( '\t', offset);
    printf ("n0=%d, n1=%d\n", n_x[0], n_x[1]);
    padding ( '\t', offset);
    //printf("n1/N=%1.2f\%\n", (double) n_x[1] / (n_x[0]+n_x[1]) * 100);
    printf("(n1/N1)/(n0/N0)=%1.2f\n", (double) n_x[1] / gbt->root->n_x[1] / n_x[0] * gbt->root->n_x[0] );

  }


}


int GBT::update() {

  map<uint64_t, OPTREE_NODE *>::iterator iter;
  int i, j, j1, j2;
  INDEX_TYPE div_ind;
  INDEX_TYPE div2_ind;
  INDEX_TYPE div3_ind;
  
  uint64_t index;

  double alpha0_b_r, alpha0_b_c_r;


  double m_divide;
  double m_stop;
  double m_divide_curr;

  double m_divide0;
  double m_stop0;
  double m_divide0_curr;

  
  double loglambda0;
 
  INDEX_TYPE n_gt1_b[3]; // indicator for whether child 1 in each dimension has at least 1 data point
  INDEX_TYPE n_gt2_b[3]; // indicator for whether child 1 in each dimension has at least 2 data point
 	

  OPTREE_NODE *NODE_CURR;
  OPTREE_NODE *CHILD_0;
  OPTREE_NODE *CHILD_1;
  int n_x_0; // number of data in left child
  int n_x_1; // number of data in right child
  int n_x_0_1, n_x_1_1, n_x_0_2, n_x_1_2;

  int n_x;

  for ( iter = gbc_map.begin(); iter != gbc_map.end(); iter++ ) { // every node in this map will have at least 2 data points
    
    NODE_CURR = iter->second;
    n_x = NODE_CURR->n_x[0] + NODE_CURR->n_x[1]; // total number of data points in the node
    

    if (NODE_CURR->k_div2 == 0 && NODE_CURR->k_div3 == 0) { // atomic node get answers directly
      
      NODE_CURR->logphi = 0;
      NODE_CURR->logQ = 0;
      NODE_CURR->log_rho = 0;
   
      //NODE_CURR->m_stop = log(rho0_b);
      //NODE_CURR->m_divide = log(1-rho0_b);

    }
    
    else if (n_x == 1) {
      
      NODE_CURR->logQ = NODE_CURR->logphi 
	= (-1.0) * (NODE_CURR->k_div2*log(2.0) + NODE_CURR->k_div3*log(3.0));
      NODE_CURR->log_rho = log(rho0_b_c);
      
      
      cout << "This should not occur by design! " << endl; 
    }

    else { // divisible node with more than one data point
 
      m_divide = 0; 
      m_divide0 = 0;

      div3_ind = (NODE_CURR->I0 & NODE_CURR->I1 & NODE_CURR->I2);
      div_ind = (NODE_CURR->I0 & NODE_CURR->I1) |(NODE_CURR->I1 & NODE_CURR->I2) |(NODE_CURR->I2 & NODE_CURR->I0)  ;
      div2_ind = div_ind & (~div3_ind);

      index = iter->first;
    
      if (NODE_CURR->k_div2 == 1 &&  NODE_CURR->k_div3 == 0)  loglambda0 = 0;
      else loglambda0 = (-1.0) * log(NODE_CURR->k_div2 + 3*NODE_CURR->k_div3);

      
	
      for (j = 0; j < 3; j++) {
	 
	n_gt1_b[j] = NODE_CURR->n_gt1_b[j][0] | NODE_CURR->n_gt1_b[j][1];
	n_gt2_b[j] = NODE_CURR->n_gt2_b[j][0] | NODE_CURR->n_gt2_b[j][1] | (NODE_CURR->n_gt1_b[j][0] & NODE_CURR->n_gt1_b[j][1]);

      }       
      

      for (i=0; i < k; i++) { // for each dimension i


	if ( (  div3_ind >> i ) & 1 ) { // if the ith dimension is intact
	  
	  alpha0_b_r = 2*alpha0_b;
	  alpha0_b_c_r = 2*alpha0_b_c;

	  for ( j = 0; j < 3 ; j++ ) {
	    
	    m_divide_curr = 0;
	    m_divide0_curr = 0;

	    j1 = (j+1) % 3;
	    j2 = (j+2) % 3;

	    if ( ( (n_gt2_b[j]) >> i ) & 1 ) { /* if child 0 in dimension i has at least 2 data points
							 then it must be in the map				       */

	      CHILD_0 = gbc_map.find( index & ( ~((uint64_t) 1 << (i+ 21*(2-j1)) | (uint64_t) 1 << i +21 *(2-j2)))  )->second;

	      n_x_0_1 = CHILD_0->n_x[0];
	      n_x_0_2 = CHILD_0->n_x[1];
	      n_x_0 = n_x_0_1 + n_x_0_2;
	      
	      m_divide0_curr += CHILD_0->logphi;
	      m_divide_curr += CHILD_0->logQ;

	    }
	
	    else {
	      
	      if  (( (NODE_CURR->n_gt1_b[j][0]) >> i ) & 1 ) {
		
		n_x_0 = 1;
		n_x_0_1 = 1;
		n_x_0_2 = 0;

		m_divide0_curr += (-1.0) * ( NODE_CURR->get_log_size() - log(3) );
		m_divide_curr += (-1.0) * ( NODE_CURR->get_log_size() - log(3) );
		
	      }

	      else if (( (NODE_CURR->n_gt1_b[j][1]) >> i ) & 1 ) {

		n_x_0 = 1;
		n_x_0_1 = 0;
		n_x_0_2 = 1;

		m_divide0_curr += (-1.0) * ( NODE_CURR->get_log_size() - log(3) );
		m_divide_curr += (-1.0) * ( NODE_CURR->get_log_size() - log(3) );
	      }

	      else {
		
		n_x_0 = 0;
		n_x_0_1 = 0;
		n_x_0_2 = 0;
	      }

	    }

	  
	    n_x_1 = n_x - n_x_0;
	    n_x_1_1 = NODE_CURR->n_x[0] - n_x_0_1;
	    n_x_1_2 = NODE_CURR->n_x[1] - n_x_0_2;

	    if (n_x_1 > 1) {

	      CHILD_1 = gbc_map.find( index & ~((uint64_t) 1 << (i+21*(2-j))) )->second;
	      	      
	      //if (n_x_1 != n_x_1_1 + n_x_1_2) cout << "Error! This should not occur!" << endl;

	      m_divide0_curr += CHILD_1->logphi;
	      m_divide_curr += CHILD_1->logQ;
	      
	    }

	    else if (n_x_1 == 1) {

	      m_divide0_curr += (-1.0) * ( NODE_CURR->get_log_size() - log(3) + log(2) );
	      m_divide_curr += (-1.0) * ( NODE_CURR->get_log_size() - log(3) +log(2) );
		
	    }

	    else {} // if n_x_1 == 0, do nothing


	    m_divide0_curr += loglambda0 + lgamma(n_x_0+alpha0_b) + lgamma(n_x_1+alpha0_b_r) 
	      - lgamma(n_x_0 + n_x_1 + alpha0_b + alpha0_b_r)
	      - ( lgamma(alpha0_b) + lgamma(alpha0_b_r) - lgamma(alpha0_b+alpha0_b_r)  ); 
	    
	    m_divide_curr += loglambda0 
	      
	      + lgamma(n_x_0_1+alpha0_b_c) + lgamma(n_x_1_1+alpha0_b_c_r) - lgamma(n_x_0_1 + n_x_1_1 + alpha0_b_c + alpha0_b_c_r)
	      - ( lgamma(alpha0_b_c) + lgamma(alpha0_b_c_r) - lgamma( alpha0_b_c + alpha0_b_c_r ) )
	      
	      + lgamma(n_x_0_2+alpha0_b_c) + lgamma(n_x_1_2+alpha0_b_c_r) - lgamma(n_x_0_2 + n_x_1_2 + alpha0_b_c + alpha0_b_c_r)
	      - ( lgamma(alpha0_b_c) + lgamma(alpha0_b_c_r) - lgamma( alpha0_b_c + alpha0_b_c_r) );
	    
	    if (m_divide0 == 0) m_divide0 = m_divide0_curr; 
	    else m_divide0 = log_exp_x_plus_exp_y (m_divide0, m_divide0_curr);
	    
	    if (m_divide == 0) m_divide = m_divide_curr;
	    else m_divide = log_exp_x_plus_exp_y (m_divide, m_divide_curr);
	     	    
	  }
	}
	  
	
	  
	if ( (div2_ind >> i) & 1 ) {
	  
	  alpha0_b_r = alpha0_b;
	  alpha0_b_c_r = alpha0_b_c;
	  
	  m_divide_curr = 0;
	  m_divide0_curr = 0;
	  
	  if ((NODE_CURR->I0 >> i) & 1) j=0;
	  else j=1;
	  
	  j1 = (j+1) % 3;
	  j2 = (j+2) % 3;
	  
	  if ( ( (n_gt2_b[j]) >> i ) & 1 ) { /* if child 0 in dimension i has more than 2 data points
						then it must be in the map				       */
	    
	    CHILD_0 = gbc_map.find( index & ( ~((uint64_t) 1 << (i+ 21*(2-j1)) | (uint64_t) 1 << i +21 *(2-j2)))  )->second;
	      
	    n_x_0 = CHILD_0->n_x[0] + CHILD_0->n_x[1];
	    n_x_0_1 = CHILD_0->n_x[0];
	    n_x_0_2 = CHILD_0->n_x[1];
	    
	    m_divide0_curr += CHILD_0->logphi;
	    m_divide_curr += CHILD_0->logQ;
	      
	  }
	  
	  else {
	    
	    if  (( (NODE_CURR->n_gt1_b[j][0]) >> i ) & 1 ) {
	      
	      n_x_0 = 1;
	      n_x_0_1 = 1;
	      n_x_0_2 = 0;
	      
	      m_divide0_curr += (-1.0) * ( NODE_CURR->get_log_size() - log(2) );
	      m_divide_curr += (-1.0) * ( NODE_CURR->get_log_size() - log(2) );
		
	    }
	    
	    else if (( (NODE_CURR->n_gt1_b[j][1]) >> i ) & 1 ) {
	      
	      n_x_0 = 1;
	      n_x_0_1 = 0;
	      n_x_0_2 = 1;
	      
	      m_divide0_curr += (-1.0) * ( NODE_CURR->get_log_size() - log(2) );
	      m_divide_curr += (-1.0) * ( NODE_CURR->get_log_size() - log(2) );
	    }
	    
	    else {
	      
	      n_x_0 = 0;
	      n_x_0_1 = 0;
	      n_x_0_2 = 0;
	    }
	    
	  }
	  
	  
	  n_x_1 = n_x - n_x_0;
	  n_x_1_1 = NODE_CURR->n_x[0] - n_x_0_1;
	  n_x_1_2 = NODE_CURR->n_x[1] - n_x_0_2;

	  if (n_x_1 > 1) {

	    CHILD_1 = gbc_map.find( index & ~((uint64_t) 1 << (i+21*(2-j))) )->second;
	    
	    m_divide0_curr += CHILD_1->logphi;
	    m_divide_curr += CHILD_1->logQ;
	    
	  }
	  
	  else if (n_x_1 == 1) {
	    
	    m_divide0_curr += (-1.0) * ( NODE_CURR->get_log_size() - log(2) );
	    m_divide_curr += (-1.0) * ( NODE_CURR->get_log_size() - log(2) );
	      
	  }
	  
	  else {} // if n_x_1 == 0, do nothing


	  m_divide0_curr += loglambda0 + lgamma(n_x_0+alpha0_b) + lgamma(n_x_1+alpha0_b_r) 
	    - lgamma(n_x_0 + n_x_1 + alpha0_b + alpha0_b_r)
	    - ( lgamma(alpha0_b) + lgamma(alpha0_b_r) - lgamma(alpha0_b+alpha0_b_r)  ); 
	  
	  m_divide_curr += loglambda0 
	      
	    + lgamma(n_x_0_1+alpha0_b_c) + lgamma(n_x_1_1+alpha0_b_c_r) - lgamma(n_x_0_1 + n_x_1_1 + alpha0_b_c + alpha0_b_c_r)
	    - ( lgamma(alpha0_b_c) + lgamma(alpha0_b_c_r) - lgamma( alpha0_b_c + alpha0_b_c_r ) )
	    
	    + lgamma(n_x_0_2+alpha0_b_c) + lgamma(n_x_1_2+alpha0_b_c_r) - lgamma(n_x_0_2 + n_x_1_2 + alpha0_b_c + alpha0_b_c_r)
	    - ( lgamma(alpha0_b_c) + lgamma(alpha0_b_c_r) - lgamma( alpha0_b_c + alpha0_b_c_r) );
	  
	  
	  if (m_divide0 == 0) m_divide0 = m_divide0_curr; 
	  else m_divide0 = log_exp_x_plus_exp_y (m_divide0, m_divide0_curr);
	  
	  if (m_divide == 0) m_divide = m_divide_curr;
	  else m_divide = log_exp_x_plus_exp_y (m_divide, m_divide_curr);
	  	  
	}

      }
      
      
      m_divide0 += log( 1 - rho0_b ); 
      m_stop0 = log(rho0_b) + n_x * (-1.0) * NODE_CURR->get_log_size();
      NODE_CURR->logphi = log_exp_x_plus_exp_y(m_stop0, m_divide0);
      
      m_divide += log( 1 - rho0_b_c );
      m_stop = log(rho0_b_c) + NODE_CURR->logphi;

      NODE_CURR->logQ = log_exp_x_plus_exp_y(m_stop, m_divide);
      NODE_CURR->log_rho = m_stop - NODE_CURR->logQ;      
    }
  }
  
  return 0;
}





