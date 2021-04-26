void print_mat(double *a,int *ia,int *ja,int N)
{
  int i,j;
  for(i=0;i<=N;i++){
    for(j=0;j<=N;j++){
      if(matrix_crs_find(ia,ja,i,j) != -1){
	printf("%10.10f ",a[matrix_crs_find(ia,ja,i,j)]);
      } else {
	printf("%10.10f ",0.0);
      }
    }
    printf("\n");
  }
}

void print_mat1(double *a,int *ia,int *ja,int N)
{
  int i,j;
  for(i=1;i<=N+1;i++){
    for(j=1;j<=N+1;j++){
      if(matrix_crs_find(ia,ja,i,j) != -1){
	printf("%10.10f ",a[matrix_crs_find(ia,ja,i,j)]);
      } else {
	printf("%10.10f ",0.0);
      }
    }
    printf("\n");
  }
}

void print_matrix(double **a,int N)
{
  int i,j;
  for(i=0;i<=N;i++){
    for(j=0;j<=N;j++){
      printf("%f ",a[i][j]);
    }
    printf("\n");
  }
}

void print_matrix1(double **a,int N)
{
  int i,j;
  for(i=1;i<=N+1;i++){
    for(j=1;j<=N+1;j++){
      printf("%f ",a[i][j]);
    }
    printf("\n");
  }
}

void print_vec(double *b,int N)
{
  int i;
  for(i=0;i<=N;i++){
    printf("%f\n",b[i]);
  }
}

void print_vec1(double *b,int N)
{
  int i;
  for(i=1;i<=N+1;i++){
    printf("%f\n",b[i]);
  }
}
 
void print_answer(int start,int end,double **npxy,double *vec_x)
{
  int i,j;
  for(i=start;i<=end;i++){
    for(j=1;j<=DIM;j++){
      printf("%10.10f ",npxy[i][j]);
    }
    printf("%10.10f\n",vec_x[i]);
  }
}
