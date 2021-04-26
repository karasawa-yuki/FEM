void init_mat_int(int **a,int i1,int i2,int j1,int j2,int c){
  int i,j;
  for(i=i1;i<=i2;i++){
    for(j=j1;j<=j2;j++){
      a[i][j] = c;
    }
  }
}

void init_mat_double(double **a,int i1,int i2,int j1,int j2,double c){
  int i,j;
  for(i=i1;i<=i2;i++){
    for(j=j1;j<=j2;j++){
      a[i][j] = c;
    }
  }
}

void init_vec_int(int *a,int i1,int i2,int c){
  int i;
  for(i=i1;i<=i2;i++){
    a[i] = c;
  }
}

void init_vec_double(double *a,int i1,int i2,double c){
  int i;
  for(i=i1;i<=i2;i++){
    a[i] = c;
  }
}

void initail_set_npxy_and_npf(double **npxy,double **npf,int start,int end,int dim,double h)
{
  int i;
  for(i=start;i<=end;i++){
    npxy[i][1]=(i-1)*h;  npf[i][1]=func(npxy[i][1]);
  }
}
