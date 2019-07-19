#include "CurvilinearTransformation.h"

#include "kernels/KernelUtils.h"


// void Linear::CurvilinearTransformation::test(){
//   return;
// }



void getBoundaryCurves3D(int num_points,
			 double offset_x, double offset_y, double offset_z,
			 double width_x, double width_y , double width_z ,
			 double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
			 double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
			 double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
			 double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
			 double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
			 double* back_bnd_x, double* back_bnd_y, double* back_bnd_z){

  double dy= 1.0/(num_points-1);
  double dx= 1.0/(num_points-1);
  
  double pi = 3.14159265359;



  int nx = 1/width_x*num_points;
  int ny = 1/width_y*num_points;
  int nz = 1/width_z*num_points;  

  double dx= 1.0/(nx-1);
  double dy= 1.0/(ny-1);
  double dz= 1.0/(nz-1);


  kernels::idx2 id_xy(ny,nx);// back front
  kernels::idx2 id_xz(nz,nx);// botton top
  kernels::idx2 id_yz(nz,ny);//left right


  
  // int i_0 =  (offset_x-0.0)/width_x;
  // int j_0 =  (offset_y-0.0)/width_y;
  

  //std::cout<<i_0<<std::endl;
  //std::cout<<j_0<<std::endl;
  //std::cout<<std::endl;

  
  // for(int i = 0 ; i< num_points; i++){
  //   left_bnd_x[i] =  offset_x;
  //   right_bnd_x[i] = width_x+offset_x;
  //   bottom_bnd_x[i] = width_x*dx*i + offset_x;
  //   top_bnd_x[i] = width_x*dx*i + offset_x;

  //   left_bnd_y[i] = width_y*dy*i + offset_y;
  //   right_bnd_y[i] = width_y*dy*i + offset_y;
  //   bottom_bnd_y[i] =offset_y;
  //   top_bnd_y[i] = width_y+offset_y + 0.0*std::sin(2*pi*top_bnd_x[i]);
  // }


  //top bottom
  for(int i = 0 ; i< nx; i++){
    for(int k = 0 ; k< nz; k++){      
      bottom_bnd_x[id_xz(k,i)] = 0+dx*i;      
      bottom_bnd_y[id_xz(k,i)] = 0;
      bottom_bnd_z[id_xz(k,i)] = 0+dz*k;      

      top_bnd_x[id_xz(k,i)] = 0+dx*i;      
      top_bnd_y[id_xz(k,i)] = 1;
      top_bnd_z[id_xz(k,i)] = 0+dz*k;      
    }
  }


  //left right
  for(int j = 0 ; j< ny; j++){
    for(int k = 0 ; k< nz; k++){
      left_bnd_x[id_yz(k,j)]  = 0;
      left_bnd_y[id_yz(k,j)]  = dy*j;
      left_bnd_z[id_yz(k,j)]  = dz*k;

      right_bnd_x[id_yz(k,j)] = 1;
      right_bnd_y[id_yz(k,j)] = dy*j;
      right_bnd_z[id_yz(k,j)] = dz*k;
    }
  }


  //front back
  for(int i = 0 ; i< nx; i++){
    for(int j = 0 ; j< ny; j++){

      front_bnd_x[id_xy(j,i)]  = dx*i;
      front_bnd_y[id_xy(j,i)]  = dy*j;
      front_bnd_z[id_xy(j,i)]  = 0;      


      back_bnd_x[id_xy(j,i)] = dx*i;
      back_bnd_y[id_xy(j,i)] = dy*j;
      back_bnd_z[id_xy(j,i)] = 1;      
    }
  }

  
  // double h0y = (top_bnd_y[0] - bottom_bnd_y[0])/(ny-1);
  // double hny = (top_bnd_y[nx-1] - bottom_bnd_y[nx-1])/(ny-1);

  // for(int i = 0 ; i< ny; i++){

  //   left_bnd_y[i]  = bottom_bnd_y[0] + h0y*i;
  //   right_bnd_y[i] = bottom_bnd_y[nx-1] + hny*i;
    
  // }
  
}




void getBoundaryCurves(int num_points,double offset_x, double offset_y,double width_x, double width_y ,double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y){

 
  
  double pi = 3.14159265359;

  int nx = 1/width_x*num_points;
  int ny = 1/width_y*num_points;

  double dx= 1.0/(nx-1);
  double dy= 1.0/(ny-1);

  // int i_0 =  (offset_x-0.0)/width_x;
  // int j_0 =  (offset_y-0.0)/width_y;
  

  //std::cout<<i_0<<std::endl;
  //std::cout<<j_0<<std::endl;
  //std::cout<<std::endl;

  
  // for(int i = 0 ; i< num_points; i++){
  //   left_bnd_x[i] =  offset_x;
  //   right_bnd_x[i] = width_x+offset_x;
  //   bottom_bnd_x[i] = width_x*dx*i + offset_x;
  //   top_bnd_x[i] = width_x*dx*i + offset_x;

  //   left_bnd_y[i] = width_y*dy*i + offset_y;
  //   right_bnd_y[i] = width_y*dy*i + offset_y;
  //   bottom_bnd_y[i] =offset_y;
  //   top_bnd_y[i] = width_y+offset_y + 0.0*std::sin(2*pi*top_bnd_x[i]);
  // }


  for(int i = 0 ; i< nx; i++){
    
    bottom_bnd_x[i] = 0 + dx*i;
    top_bnd_x[i]    = 0 + dx*i;
    
    bottom_bnd_y[i] = 0.0;
    top_bnd_y[i]    = 1 + 0.0*std::sin(2*pi*top_bnd_x[i]);

    //std::cout<<top_bnd_x[i]<<std::endl;
    //std::cout<<i<<std::endl;
    //std::cout<<std::endl;
  }

  //std::exit(-1);
  
  for(int i = 0 ; i< ny; i++){
    
    left_bnd_x[i]  = 0;
    right_bnd_x[i] = 1;
    bottom_bnd_x[i] = dx*i;
    top_bnd_x[i] = dx*i;

    left_bnd_y[i] = dy*i;
    right_bnd_y[i] = dy*i;
    bottom_bnd_y[i] = 0;
    top_bnd_y[i] = 1;
  }
}


void transFiniteInterpolation3D(int mx, int my, int mz,
				int k_m, int k_p ,
				int j_m, int j_p ,
				int i_m, int i_p ,
				int num_points,
				double* left_bnd,
				double* right_bnd,
				double* bottom_bnd,
				double* top_bnd,
				double* front_bnd,
				double* back_bnd,
				double* curvilinear, int dim
				){
  kernels::idx3 id_xyz(num_points,num_points,num_points);
  kernels::idx2 id_xy(my,mx);// back front
  kernels::idx2 id_xz(mz,mx);// bottom top
  kernels::idx2 id_yz(mz,my);//left right

  double mesh_size_x = 1.0/(mx-1);
  double mesh_size_y = 1.0/(my-1);
  double mesh_size_z = 1.0/(mz-1);

  int i_0;
  int j_0;
  int k_0;    

  double u,v,w;
  double uv,vw,uw,uvw1,uvw2;

  double q,r,s;


  for(int i=i_m ; i< i_p ; i++){
    i_0=i-i_m;
    q=mesh_size_x * i;	
    for(int j=j_m ; j< j_p ; j++){
      j_0=j-j_m;
      r=mesh_size_y * j;
      for(int k=k_m ; k< k_p ; k++){
	k_0=k-k_m;
	s=mesh_size_z * k;
	
	u=(1-q)*left_bnd[id_yz(k,j)]+q*right_bnd[id_yz(k,j)];

	v=(1-r)*top_bnd[id_xz(k,i)]+r*bottom_bnd[id_xz(k,i)];

	w=(1-s)*front_bnd[id_xy(j,i)]+s*back_bnd[id_xy(j,i)];

	uw=(1-q)* ((1 -r)*left_bnd[id_yz(k,1)] 
		   +r *left_bnd[id_yz(k,my)])
	  +q*((1-r)*right_bnd[id_yz(k,1)] 
	      +r *right_bnd[id_yz(k,my)]);

	uv=(1-r)*((1-s)*top_bnd[id_xz(1,i)]
		  + s*top_bnd[id_xz(mz,i)])
	  + r*((1-s)*bottom_bnd[id_xz(1,i)]
	       + s*bottom_bnd[id_xz(mz,i)]);

	vw = (1-s)*((1-q)*front_bnd[id_xy(j,1)] 
		    + q*front_bnd[id_xy(j,mx)])
	       + s*((1-q)*back_bnd[id_xy(j,1)] 
	              + q*back_bnd[id_xy(j,mx)]);
	
	uvw1=(1-q)*((1-r)*((1-s)*top_bnd[id_xz(1,1)] 
		              +s*top_bnd[id_xz(mz,1)])
		       +r*((1-s)*bottom_bnd[id_xz(1,1)]
			   +s*bottom_bnd[id_xz(mz,1)]));

	uvw2=q*((1-r)*((1-s)*top_bnd[id_xz(1,mx)] 
			 + s*top_bnd[id_xz(mz,mx)])
		  + r*((1-s)*bottom_bnd[id_xz(1,mx)]
		       + s*bottom_bnd[id_xz(mz,mx)]));

	curvilinear[id_xyz(k_0,j_0,i_0)] = u + v + w - uv - uw - vw + uvw1 + uvw2;

      }
    }
  }
}





void transFiniteInterpolation(int mx, int my, int j_m, int j_p, int i_m, int i_p, int num_points, double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y, double* curvilinear_x , double* curvilinear_y ){

   double mesh_size_x = 1.0/(num_points-1);
   double mesh_size_y = 1.0/(num_points-1);

   double r;
   double q;
   kernels::idx2 id_xy(num_points,num_points);


   for(int j = 0 ; j < num_points ; j ++){
     //printf("%f \n",top_bnd_y[j]);
} 
   //printf("\n");
   for(int j =0 ; j < num_points ; j++) {
   for(int i =0 ; i < num_points ; i++) {

      	q = (i)*mesh_size_x;
        r = (j)*mesh_size_y;
        
        curvilinear_x[id_xy(j,i)] = (1-q)*left_bnd_x[j]+q*right_bnd_x[j]+(1-r)*bottom_bnd_x[i]+r*top_bnd_x[i]-
               (1-q)*(1-r)*left_bnd_x[0]-q*(1-r)*right_bnd_x[0]-r*(1-q)*top_bnd_x[0]-
               (r*q)*top_bnd_x[num_points-1];


        curvilinear_y[id_xy(j,i)] = (1-q)*left_bnd_y[j]+q*right_bnd_y[j]+(1-r)*bottom_bnd_y[i]+r*top_bnd_y[i]-
               (1-q)*(1-r)*left_bnd_y[0]-q*(1-r)*right_bnd_y[0]-r*(1-q)*top_bnd_y[0]-
               (r*q)*top_bnd_y[num_points-1];


   }
   }

   //  for(int j = 0 ; j < num_points_y ; j ++){
    // printf("%f \n",curvilinear_y[id_xy(j,num_points_y-1)]);
   //} 

}


double lagrangeBasis(double x,double* points,int i,int num_points){
  double result=1;
  for (int j = 0 ; j< num_points ; j ++){

    if (j != i) {
      result *= (x-points[j])/(points[i]-points[j]);
    }
  }
   
  return result;
}


void interpolate(double x, double y, double* orig_mesh_x , double* orig_mesh_y, double* dest_mesh, int num_nodes,double& result){

  double a_x=0;
  double a_y=0;
  
  result=0;
  
  kernels::idx2 id_xy(num_nodes,num_nodes);

  
  for (int j = 0 ; j< num_nodes ; j ++){
    for (int i = 0 ; i< num_nodes ; i ++){
      a_x=lagrangeBasis(x,orig_mesh_x,i,num_nodes);
      a_y=lagrangeBasis(y,orig_mesh_y,j,num_nodes);
      result += dest_mesh[id_xy(j,i)] * a_x*a_y;
    }
  }
}

void getValuesAtQuadNodes(double* orig_mesh_x , double* orig_mesh_y, double* dest_mesh, int num_nodes, double* results){

  kernels::idx2 id_xy(num_nodes,num_nodes);

  for (int j = 0 ; j< num_nodes ; j ++){
    for (int i = 0 ; i< num_nodes ; i ++){
      interpolate(kernels::legendre::nodes[num_nodes-1][i],kernels::legendre::nodes[num_nodes-1][j],orig_mesh_x,orig_mesh_y,dest_mesh,num_nodes,results[id_xy(j,i)]);
    }
  }
}

void computeDerivatives_x(int i, int j , double* values , int num_nodes, double& der_x){
  //  double* vals_x;
  //  double* vals_y

    //  vals_x = malloc(sizeof(double)*num_nodes*num_nodes);
    //  vals_y = malloc(sizeof(double)*num_nodes*num_nodes);


  //  free(vals_x);
  //  free(vals_y);
  
   kernels::idx2 id_xy(num_nodes,num_nodes);

   der_x = 0.0;

  // for (int j = 0 ; j< num_nodes ; j ++){
  //   for (int i = 0 ; i< num_nodes ; i ++){
   for (int n = 0 ; n< num_nodes ; n ++){
     der_x += kernels::dudx[num_nodes-1][i][n] * values[id_xy(j,n)];
   }
  //   }
  // }

}

void computeDerivatives_y (int i, int j , double* values , int num_nodes, double& der_y){
  //  double* vals_x;
  //  double* vals_y
    //  vals_x = malloc(sizeof(double)*num_nodes*num_nodes);
    //  vals_y = malloc(sizeof(double)*num_nodes*num_nodes);
  //  free(vals_x);
  //  free(vals_y);
  
  kernels::idx2 id_xy(num_nodes,num_nodes); 

  der_y = 0.0;
  // for (int j = 0 ; j< num_nodes ; j ++){
  //   for (int i = 0 ; i< num_nodes ; i ++){
      for (int n = 0 ; n< num_nodes ; n ++){
         der_y += kernels::dudx[num_nodes-1][j][n] * values[id_xy(n,i)];
      }
  //   }
  // }
}



void metricDerivativesAndJacobian(int num_nodes, double* curvilinear_x, double* curvilinear_y,double* gl_vals_x,double* gl_vals_y,double* q_x,double* q_y,double* r_x,double* r_y,double* jacobian){

  double* unif_mesh= new double[num_nodes];
  for(int i = 0; i< num_nodes ; i++){
    unif_mesh[i]=i*1.0/(num_nodes-1);
  }

  getValuesAtQuadNodes(unif_mesh,unif_mesh,curvilinear_x,num_nodes,gl_vals_x);
  getValuesAtQuadNodes(unif_mesh,unif_mesh,curvilinear_y,num_nodes,gl_vals_y);

  double x_der_x; 
  double x_der_y;
  double y_der_x; 
  double y_der_y;
  
  kernels::idx2 id_xy(num_nodes,num_nodes);
  
  for(int j = 0 ; j < num_nodes ; j ++){
    for(int i = 0 ; i< num_nodes ; i ++){

     computeDerivatives_x(j,i,gl_vals_x,num_nodes,x_der_x);
     computeDerivatives_y(j,i,gl_vals_x,num_nodes,x_der_y);

     computeDerivatives_x(j,i,gl_vals_y,num_nodes,y_der_x);
     computeDerivatives_y(j,i,gl_vals_y,num_nodes,y_der_y);

// J = x_q.*y_r - y_q.*x_r; % Jacobian (determinant of metric)
// q_x = y_r./J;   % 1/J y_r
// r_x = - y_q./J; % -1/J y_q
// q_y = - x_r./J; % -1/J x_r
// r_y = x_q./J;   %  1/J x_q
		     
     jacobian[id_xy(j,i)]=x_der_x*y_der_y-x_der_y*y_der_x;
     
     q_x[id_xy(j,i)]=y_der_y/jacobian[id_xy(j,i)];
     q_y[id_xy(j,i)]=-x_der_y/jacobian[id_xy(j,i)];
     r_x[id_xy(j,i)]=-y_der_x/jacobian[id_xy(j,i)];
     r_y[id_xy(j,i)]=x_der_x/jacobian[id_xy(j,i)];

     //      printf("%f\n",der_x);
     //      printf("%f\n",der_y);
    }
 }

  


}


// void computeCurvilinearTransformationCoefficients(double* dest_mesh_x,double* dest_mesh_y, int num_nodes){

// double* dest_mesh_x_der_x;
// double* dest_mesh_x_der_y;
// double* dest_mesh_y_der_x;
// double* dest_mesh_y_der_y;

// dest_mesh_x_der_x= malloc(sizeof(double)*num_nodes*num_nodes);
// dest_mesh_x_der_y= malloc(sizeof(double)*num_nodes*num_nodes);
// dest_mesh_x_der_y= malloc(sizeof(double)*num_nodes*num_nodes);


// }
