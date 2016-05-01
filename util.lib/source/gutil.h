#ifndef PARALENGTH
#define PARALENGTH 30
#endif

#ifndef HBAR
#define HBAR 1
#endif

#ifndef AMU2AU
#define AMU2AU 1822.88839
#endif

#ifndef CM2HARTREE
#define CM2HARTREE 1.000/219474.63
#endif

#ifndef A2BOHR
#define A2BOHR 0.529177249
#endif



typedef struct PARA
{
  char field[50];
  char value[90];

} para;



/*
void pickdfdr( FILE * p_of_dfdr_file, double * p_of_dfdr_triangle, int len_dfdr_file,
               int len_dfdr, char * kw, int mode_number);
*/

void pick( FILE * p_of_info_file, char * kw, int nskip, double * database );

void dx2dr( double * u, int ulda, double * h, int hlda, double * hnorm);

void masswt( int na, double * mass, double * unwtmat, double * wtmat);

void massdivide( char direction, int na, double * mass, double * undivided, double * divided);

void mo2nbo(int natom, int nbasis, double * u, double * pre, double * post);

void readgauss( char * flag, FILE * gaussoutput, int * ires, double * dres );

double boltzmannfactor( double freq, int quanta, double T );

double boltzmannAvgEnergy( double freq , double T );

int dropball( double rnd , double freq , double temperature );

void effquanta( int nmode, int irmode, double temperature, 
                double * freqlist, double cutoff, double * quanta);

void thetasequencegen( int nmode, int ntraj, int * seed, double * theta);

void thetagen( int nmode, int * seed, double * theta);

void reducedmasscalc( int natom, double * mass, double * hesscart, char mwcornot, 
                      double * reducedmass);

void checkDistribution( int ntraj , int nbin , double * sample , 
                        double *xx , int *yy );


int tellIR_range( double position , double front , double tail ) ;

int tellIR_list( int place , int * list , int length_of_list ) ;

int tellActive_list( int place , int * list , int length_of_list ) ;

void loadCNDOdipole( char direction , FILE * pcndoLogFile , double * buffer ) ;


// Fortran subroutines ...

void gausvib_( int * , double * , double * , double * , double * ) ;
//            ( natom, mass, coordxyzinp, hessianinp, Dmatrix )








