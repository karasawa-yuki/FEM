#include<stdio.h>
#include<stdlib.h>
#include<math.h>

/*
   function f
*/
double func(double x)
{
  return 1;
  //return pi*pi*sin(pi*x);
}

/*
  Conversion from elnp to npnp
*/
void from_elnp_to_npnp(int **elnp, int **npnp, int nn,int mm)
{
  int i,j,k,l,count;

  for(i=1;i<=nn+1;i++){
    for(j=1;j<=mm;j++){
      npnp[i][j] = -1;
    }
  }

  for(k=1;k<=nn;k++){
    for(i=1;i<=elnp[k][0];i++){
      for(j=1;j<=elnp[k][0];j++){
	count = 1;
	while(npnp[elnp[k][i]][count] >= 0){
	  if(npnp[elnp[k][i]][count] == elnp[k][j]){
	    count = -1;
	    break;
	  }
	  count++;
	}
	if(count != -1){
	  npnp[elnp[k][i]][count] = elnp[k][j];
	  npnp[elnp[k][i]][0] = count;
	}
      }
    }
  }
}




void make_elnp(int **elnp,int start, int end, int end1)
{
  int i,k;
  for(k=start;k<=end;k++){
    elnp[k][0] = end1;
    for(i=1;i<=end1;i++){
      elnp[k][i] = k+i-1;
    }
  }
}



void make_npel(int **elnp,int **npel,int start,int end,int end1)
{
  int i,j,k,count;
  
  for(i=start;i<=end;i++){
    npel[i][0] = DIM+1;
    for(j=1;j<=end1;j++){
      npel[i][j] = -1;
    }
  }

  for(k=start;k<=end;k++){
    for(i=1;i<=elnp[k][0];i++){
      count = 1;
      while(npel[elnp[k][i]][count] >= 0){
	count++;
      }
      npel[elnp[k][i]][count] = k;
    }
  }
}


void sort(int *sort_vec)
{
  int i,j,dummy;
  
  for(i=1;i<DIM;i++){
    for(j=i+1;j<=DIM;j++){
      if(sort_vec[i] > sort_vec[j]){
        dummy = sort_vec[i];
	sort_vec[i] = sort_vec[j];
	sort_vec[j] = dummy;
      }
    }
  }
}

int diff_make_elel(int *diff1,int *diff2)
{
  int i;
  for(i=1;i<=DIM;i++){
    if(diff1[i] != diff2[i]){
      return -1;
    }
  }
  return 1;
}

double ek_mat(double **npxy,double **mat_a,int **elnp,int k)
{
  int i,j,l;
  double h=1.,lsv,dummy;
  double **mat,**phi;

  mat = dmatrix(1,DIM+1,1,DIM+1);    phi = dmatrix(1,DIM+1,1,DIM+1);
  init_mat_double(mat,1,DIM+1,1,DIM+1,0);  init_mat_double(phi,1,DIM+1,1,DIM+1,0); 
  for(i=1;i<=DIM+1;i++){
    for(j=1;j<=DIM+1;j++){
      mat_a[i][j] = 0.0;
    }
  }

  for(i=1;i<=elnp[k][0];i++){
    mat[i][1] = 1;
    for(j=2;j<=DIM+1;j++){
      mat[i][j] = npxy[elnp[k][i]][j-1];
    }
  }

   
  lsv = det_mat(mat);
  make_phi(phi,npxy,elnp,k);
  for(i=1;i<=DIM;i++){
    h = h/(double)i;
  }
  for(i=1;i<=DIM+1;i++){
    for(j=i;j<=DIM+1;j++){
      dummy = 0;
      for(l=2;l<=DIM+1;l++){
	dummy += phi[i][l]*phi[j][l];
      }
      mat_a[i][j] = dummy*h*lsv;
      if(i!=j){
	mat_a[j][i] = mat_a[i][j];
      }
    }
  }
  free_dmatrix(mat,1,DIM+1,1,DIM+1);     free_dmatrix(phi,1,DIM+1,1,DIM+1);
  return lsv;
}


double det_mat(double **a)
{
  int i,j;
  double sum=0;
  if(DIM==1){
    return a[1][1]*a[2][2] - a[1][2]*a[2][1];
  } else if(DIM==2){
    sum += a[1][1]*(a[2][2]*a[3][3] - a[2][3]*a[3][2]);
    sum += a[2][1]*(a[3][2]*a[1][3] - a[2][1]*a[3][3]);
    sum += a[3][1]*(a[1][2]*a[2][3] - a[1][3]*a[2][2]);
    return sum;
  } else if(DIM==3){
    for(i=2;i<=DIM+1;i++){
      for(j=2;j<=DIM+1;j++){
	a[i][j] -= a[1][j];
      }
    }
    sum += a[2][2]*(a[3][3]*a[4][4] - a[3][4]*a[4][3]);
    sum += a[3][2]*(a[4][3]*a[2][4] - a[3][2]*a[4][4]);
    sum += a[4][2]*(a[2][3]*a[3][4] - a[2][4]*a[3][3]);
    return sum;
  }
}

void make_phi(double **phi,double **npxy,int **elnp,int k)
{
  int i,j,l,m;
  double *b,**a,**a_cp,a_det,ab_det;
  b = dvector(1,DIM+1);  a = dmatrix(1,DIM+1,1,DIM+1);  a_cp = dmatrix(1,DIM+1,1,DIM+1);
  init_vec_double(b,1,DIM+1,0); init_mat_double(a,1,DIM+1,1,DIM+1,0);
  init_mat_double(a_cp,1,DIM+1,1,DIM+1,0);
  
  for(i=1;i<=DIM+1;i++){
    a[i][1] = 1.0;
    for(j=2;j<=DIM+1;j++){
      a[i][j] = npxy[elnp[k][i]][j-1];
    }
  }
  a_det = det_mat(a);

  for(m=1;m<=DIM+1;m++){
    for(j=1;j<=DIM+1;j++){
      b[j] = 0.0;
    }
    b[m] = 1.0;
    for(l=1;l<=elnp[k][0];l++){
      for(i=1;i<=DIM+1;i++){
	for(j=1;j<=DIM+1;j++){
	  if(j==l){
	    a_cp[i][j] = b[i];
	  } else {
	    a_cp[i][j] = a[i][j];
	  }
	}
      }
      ab_det = det_mat(a_cp);
      phi[m][l] = ab_det/a_det;
    }
  }
  free_dvector(b,1);    free_dmatrix(a,1,DIM+1,1,DIM+1);    free_dmatrix(a_cp,1,DIM+1,1,DIM+1);
}

void make_m(double **m)
{
  int i,j;
  double h=1.,d=1.;
  for(i=1;i<=DIM+1;i++){
    for(j=1;j<=DIM+1;j++){
      m[i][j] = 0.0;
    }
  }
  for(i=1;i<=DIM+2;i++){
    h = h/(double)i;
  }
  for(i=1;i<=DIM;i++){
    d = d*(double)i;
  }
  for(i=1;i<=DIM+1;i++){
    for(j=i;j<=DIM+1;j++){
      if(i==j){
        m[i][j] = 2.*d*h;
      }else{
        m[i][j] = 1.*d*h;
	m[j][i] = m[i][j];
      }
    }
  }
}

void from_elnp_and_npel_to_elel(int **npel,int **elnp,int **elel,int start,int end)
{
  int i,j,ii,jj,k,l,count1,count2,damy;
  int *hm,*hm1;
  hm = ivector(1,DIM);  hm1 = ivector(1,DIM);
  
  for(k=start;k<=end;k++){
    elel[k][0] = DIM+1;
    for(j=1;j<=DIM+1;j++){
      elel[k][j] = -1;
    }
  free_ivector(hm,1);     free_ivector(hm1,1);
}
