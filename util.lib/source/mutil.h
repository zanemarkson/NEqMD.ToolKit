//-----------------  Macros  -------------------//
#ifndef PARALENGTH
#define PARALENGTH 30
#endif

#ifndef HBAR
#define HBAR 1
#endif

#ifndef AMU2AU
#define AMU2AU 1822.88839
#endif


#ifndef PI
#define PI 3.1415926535
#endif



#define IA 16807 

#define IM 2147483647 

#define AM (1.0/IM) 

#define IQ 127773 

#define IR 2836 

#define MASK 123459876


//------------------------------------------------------------//








// Matrix Utility Routines ...


void dtri2sym( char uplo, int lda, double * tri, double * sym);

void ztri2sym( char uplo, int lda, double _Complex * tri, double _Complex * sym);

void dmatprint(int dimrow, int dimcol, double *p);

void double_allocate( int size, double * head);

void dnummat( double coeff, double * inp, int rlda, int clda, int rdim, int cdim, double * out );

void dmatplus( int rdim, int cdim, double * left, double * right, int leftlda, int rightlda, double * res);







void dkron( double *M1, int dimrow1, int dimcol1, 
            double *M2, int dimrow2, int dimcol2, 
            double *R );
void dmatmul( int dimrowleft, int dimcolright, int dimcommon, 
              double *left, double *right, double *result );

void dmatassgn(int dimrow, int dimcol,  double  *porigin, double  *ptarget);

void dmatrixprint(int dimrow, int dimcol, double *p);

void dzeros(int dimrow, int dimcol, double *p);

void dones(int dimrow, int dimcol, double *p);

int  double_factorial(int x);

double dtrace(int dimofmat, double *p);

void dtranspose( int dimofmat, double *porigin, double *ptarget );

void dvec2mat(int dimrow, int dimcol, double *pvec, double *pmat);

int  factorial(int x);

void zkron( double _Complex *M1, int dimrow1, int dimcol1, 
            double _Complex *M2, int dimrow2, int dimcol2, 
            double _Complex *R );
            
void zmatmul( int dimrowleft, int dimcolright, int dimcommon, 
              double _Complex *left, double _Complex *right, double _Complex *result );

void zmatassgn(int dimrow, int dimcol,  double _Complex *porigin, double _Complex *ptarget);

void zmatprint(int dimrow, int dimcol, double _Complex *p);

void zones(int dimrow, int dimcol, double _Complex *p);

double _Complex  ztrace(int dimofmat, double _Complex *p);

void ztranspose( int dimofmat, double _Complex *porigin, double _Complex *ptarget );

void zvec2mat(int dimrow, int dimcol, double _Complex *pvec, double _Complex *pmat);

void zzeros(int dimrow, int dimcol, double _Complex *p);

void tril(int nofline, int dimofmat, double *porigin, double *ptarget);

void triu(int nofline, int dimofmat, double *porigin, double *ptarget);

void dLiouvillian(int dim, double *hilbert, double *L);

void zLiouvillian(int dim, double *hilbert, double _Complex *L);

void imatrixprint(int dimrow, int dimcol, int *p);

void ioutput(FILE *pout, int dimrow, int dimcol, int *M);

void loutput(FILE *pout, int dimrow, int dimcol, long *M);

void doutput(FILE *pout, int dimrow, int dimcol, double *M);

void zoutput(FILE *pout, int dimrow, int dimcol, double _Complex *M);

void zdmatmul( int dimrowleft, int dimcolright, int dimcommon, 
               double _Complex *left, double *right, double _Complex *result );

void dzmatmul( int dimrowleft, int dimcolright, int dimcommon, 
               double  *left, double _Complex *right, double _Complex *result );

void iones(int dimrow, int dimcol, int *p);

void izeros(int dimrow, int dimcol, int *p);


void mark(int choice);

double dmax( int length , double * sample );

double dmaxAbs( int length , double * sample ) ;

int dmaxID( int length , double * sample ) ;

int dmaxAbsID( int length , double * sample ) ;

double dmin( int length , double * sample );

double dminAbs( int length , double * sample );

int dminID( int length , double * sample ) ;

int dminAbsID( int length , double * sample ) ;

int imax( int length , int * sample );

int imin( int length , int * sample );

long ranzm( long * idum );

void dtranspose_nonsquare( int nrow , int ncol , double *porigin, double *ptarget ) ;

double hermite( int n , double x ) ;

double gaussianDistribution( double sigma , double mu , double x ) ;

double hcPsi( int n , double xp , double w , double A );






// Fortran-to-C interfaces 

void dsyev_f2c( int dimension, double * sym, double * eigvectors,  double * eigvalues  );






// Fortran Routines: BLAS.F, LAPACK.F, LINPACK.F

void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *, int * );

void zgemm_(char *, char *, int *, int *, int *, double _Complex *,
            double _Complex *, int *, double _Complex *, int *,
            double _Complex *, double _Complex *, int *);

void dgemm_(char *, char *, int *, int *, int *, double  *,
            double  *, int *, double  *, int *,
            double  *, double  *, int *);

void dgemv_(char *, int *, int *, double *, double *, int *, 
		double *, int *, double *, double *, int *);
		
// void nomalisevector_( int * , double * ) ;

// void matrixdiagonal_( int * , double * , double * , double * ) ;

// void dotproduct_( int * , double * , double * , double * ) ;

// schmidtorthogonal_( int * , int * , double * ) ;












