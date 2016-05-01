#ifndef PARALENGTH
#define PARALENGTH 30
#endif

#ifndef HBAR
#define HBAR 1
#endif

#ifndef AMU2AU
#define AMU2AU 1822.88839
#endif




void fload( FILE * pfile, double * buffer);

void fnload( FILE * pfile, double * buffer , int _nnumber_ ) ;

void int_fload( FILE * pfile, int * buffer) ;

void int_fnload( FILE * pfile, double * buffer , int _nnumber_ ) ;

int flength( FILE * pfile);

int fwc( FILE * pfile );

void fcopy(FILE * ctrlc, FILE * ctrlv);

void fsub(FILE * ctrlc, char ctrlf, char ctrls, FILE * ctrlv);

int fsearch( FILE * pfile, char kw[]  );

void fdel( char filename[] );

int stellblank( char * line ) ;

int freadline( char * s , int n , register FILE * iop , char symbol ) ;

int strwordcount( char * string ) ;

void strpickword( char * string , int position , char * cache ) ;

int preLoadEntry( FILE * pfile , char entryName[] , char commentSymbol ) ;

int preLoadGRO( FILE * pgroinput ) ;

int fgrepwc( FILE * pfile , char keyword[] ) ;

