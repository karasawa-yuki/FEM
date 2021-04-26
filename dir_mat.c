#include"./functions.c"

int main(int argc,char *argv[])
{
  int N = atoi(argv[1]);
  int np = N-1, ne = N,nl = N+1;
  double h = 1./N,z0=0,z1=0;
  int i,j,k;
  double **npxy,**mat_a,**mat_m,**npf;
  int **elnp,**npnp;
  double *vec_b,*vec_x;

  // make matrix and vector
  elnp = imatrix(1,N,0,2); init_mat_int(elnp,1,N,0,2,0);
  npnp = imatrix(1,N+1,0,3); init_mat_int(npnp,1,N+1,0,3,0);
  npxy = dmatrix(1,N+1,1,DIM); init_mat_double(npxy,1,N+1,1,DIM,0);
  npf = dmatrix(1,N+1,1,DIM); init_mat_double(npf,1,N+1,1,DIM,0);
  mat_a = dmatrix(1,N+1,1,N+1); init_mat_double(mat_a,1,N+1,1,N+1,0);
  mat_m = dmatrix(1,N+1,1,N+1); init_mat_double(mat_m,1,N+1,1,N+1,0);
  vec_b = dvector(1,N+1); init_vec_double(vec_b,1,N+1,0);
  vec_x = dvector(1,N+1); init_vec_double(vec_x,1,N+1,0);

  
  // initializatioon
  make_elnp(elnp,1,N,DIM+1);
  initail_set_npxy_and_npf(npxy,npf,1,N+1,DIM,h);
  // initail setting
  turn_by_element_mat(mat_a,mat_m,N,npxy,elnp);
  make_vec_b(mat_a,vec_b,mat_m,npf,h,z0,z1,np,N);
  vec_x[1]=z0;      vec_x[N+1]=z1;

  cg_method_matrix(mat_a,vec_b,vec_x,1,N+1);

  // print to answer
  print_answer(1,N+1,npxy,vec_x);
  
  free_imatrix(elnp,1,N,0,2);       free_imatrix(npnp,1,N+1,0,2);
  free_dmatrix(npxy,1,N+1,1,DIM);     free_dmatrix(npf,1,N+1,1,DIM);
  free_dmatrix(mat_a,1,N+1,1,N+1);    free_dmatrix(mat_m,1,N+1,1,N+1);
  free_dvector(vec_b,1);              free_dvector(vec_x,1);
}
