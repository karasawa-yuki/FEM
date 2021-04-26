void make_vec_b(double **a,double *vec_b,double **mat_m,double **npf,double h,double z0,double z1,int np,int N)
{
  int i,j;
  for(i=2;i<=np+1;i++){
    for(j=1;j<=N+1;j++){
      vec_b[i] += mat_m[i][j]*npf[j][1];
    }
  }
  a[1][1] = 1;
  a[N+1][1] = 0;
  for(j=2;j<=N;j++){
    a[1][j] = 0;
    a[N+1][j] = 0;
  }
  a[1][N+1] = 0;
  a[N+1][N+1] = 1;
  for(i=2;i<=np+1;i++){
    vec_b[i] -= z0*a[i][1] + z1*a[i][N+1];
    a[i][1] = 0;  a[i][N+1] = 0;
  }
  vec_b[1]=z0;
  vec_b[N+1]=z1;
}

void cg_method_matrix(double **mat_a,double *vec_b,double *vec_x,int start,int end)
{
  int i,j,k;
  double eps=1.0e-10,r_inner,alpha,beta;
  double *vec_r,*vec_ap,*vec_p;

  // make vector
  vec_r=dvector(start,end);    vec_p=dvector(start,end);     vec_ap=dvector(start,end);
  init_vec_double(vec_r,start,end,0); init_vec_double(vec_p,start,end,0);
  init_vec_double(vec_ap,start,end,0);

  // ready of cg method
  for(i=start;i<=end;i++){
    vec_r[i]=vec_b[i];
    for(j=start;j<=end;j++){
      vec_r[i] -= mat_a[i][j]*vec_x[i];
    }
    vec_p[i]=vec_r[i];
   }

  // cg method
  for(k=0; ;k++){
    r_inner = inner_product(start,end,vec_r,vec_r);
    matrix_vector_product(mat_a,vec_p,vec_ap,end);
    alpha = inner_product(start,end,vec_r,vec_p)/ inner_product(start,end,vec_ap,vec_p);
    for(i=start;i<=end;i++){
      vec_x[i] += alpha * vec_p[i];
      vec_r[i] -= alpha * vec_ap[i];
    }
    if(sqrt(inner_product(start,end,vec_r,vec_r)) < eps) break;
    beta = inner_product(start,end,vec_r,vec_r)/r_inner;
    for(i=start;i<=end;i++){
      vec_p[i] = vec_r[i] + beta*vec_p[i];
    }
  }

  // free matrix and vector
  free_dvector(vec_r,start);     free_dvector(vec_p,start);      free_dvector(vec_ap,start);
}

void turn_by_element_mat(double **mat_a,double **mat_m,int nn,double **npxy,int **elnp)
{
  int i,j,k;
  double **a,**m,meas;
  a = dmatrix(1,DIM+1,1,DIM+1);    m = dmatrix(1,DIM+1,1,DIM+1);
  init_mat_double(a,1,DIM+1,1,DIM+1,0);  init_mat_double(m,1,DIM+1,1,DIM+1,0);
  make_m(m);
  
  for(k=1;k<=nn;k++){
    meas = ek_mat(npxy,a,elnp,k);
    for(i=1;i<=DIM+1;i++){
      for(j=1;j<=DIM+1;j++){
	mat_a[elnp[k][i]][elnp[k][j]] += a[i][j];
	mat_m[elnp[k][i]][elnp[k][j]] += meas*m[i][j];
      }
    }
  }
  free_dmatrix(a,1,DIM+1,1,DIM+1);  free_dmatrix(m,1,DIM+1,1,DIM+1);  
}
