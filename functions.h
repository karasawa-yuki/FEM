#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <math.h>

#define pi 4.0*atan(1.0)
#define DIM 1
#define Max(a,b) (((a)>(b))?(a):(b))
#define Min(a,b) (((a)<(b))?(a):(b))
#define Sgn(a) (((a)==0.0)?0.0:((a)/(fabs(a))))

double arg( double *x );
double *dvector(int i, int j);
void free_dvector(double *a, int i);
int *ivector(int i, int j);
void free_ivector(int *a, int i);
double **dmatrix(int nr1, int nr2, int nl1, int nl2);
void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2);
int **imatrix(int nr1, int nr2, int nl1, int nl2);
void free_imatrix(int **a, int nr1, int nr2, int nl1, int nl2);
void matrix_vector_product( double **a, double *b, double *c, int n );
double vector_norm1( double *a, int m, int n );
double inner_product( int m, int n, double *a, double *b);
void cg( double **a, double *b, double *x, int n, double EPS, int KMAX );


//my_func.c
double func(double x);
void from_elnp_and_npel_to_elel(int **npel,int **elnp,int **elel,int start,int end);
void make_elnp(int **elnp,int start, int end, int end1);
void make_npel(int **elnp,int **npel,int start,int end,int end1);
void sort(int *sort_vec);
int diff_make_elel(int *diff1,int *diff2);
double ek_mat(double **npxy,double **mat_a,int **elnp,int k);
void from_elnp_to_npnp(int **elnp, int **npnp, int nn,int mm);
double det_mat(double **a);
void make_phi(double **phi,double **npxy,int **elnp,int k);
void make_m(double **m);
void initail_set_npxy_and_npf(double **npxy,double **npf,int start,int end,int dim,double h);


//crs_fanc.c
void a_and_b_dirichlet_crs1(double *a,int *ia,int *ja,double *b,double z0,double z1,int N);
int matrix_crs_find(int *ia, int *ja, int i, int j);
double matrix_value_crs(double *a,int *ia,int *ja,int N,int i,int j);
void make_b_crs(double *vec_b,double *mat_m,int *mat_im,int *mat_jm,double **npf,int N);
void turn_by_element_a_and_m_crs(double *mat_a,int *mat_ia,int *mat_ja,double *mat_m,int **npnp,int **elnp,double **npxy,int count,int n,int N,double h);
void cg_method_crs_0_N(int *mat_i, int *mat_j, double *mat, double *vec_b, double *vec_x, int start, int end);


//mat_fanc.c
void make_vec_b(double **a,double *vec_b,double **mat_m,double **npf,double h,double z0,double z1,int np,int N);
void cg_method_matrix(double **mat_a,double *vec_b,double *vec_x,int start,int end);
void turn_by_element_mat(double **a,double **m,int nn,double **npxy,int **elnp);


//initial.c
void init_mat_int(int **a,int i1,int i2,int j1,int j2,int c);
void init_mat_double(double **a,int i1,int i2,int j1,int j2,double c);
void init_vec_int(int *a,int i1,int i2,int c);
void init_vec_double(double *a,int i1,int i2,double c);


//func_print.c
void print_mat(double *a,int *ia,int *ja,int N);
void print_mat1(double *a,int *ia,int *ja,int N);
void print_matrix(double **a,int N);
void print_matrix1(double **a,int N);
void print_vec(double *b,int N);
void print_vec1(double *b,int N);
void print_answer(int start,int end,double **npxy,double *vec_x);

#endif
