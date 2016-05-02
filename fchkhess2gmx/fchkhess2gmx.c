#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <math.h>
#include <time.h>
#include "gutil.h"
#include "mutil.h"
#include "futil.h"



#ifndef YES
#define YES 1
#endif

#ifndef NO
#define NO 0
#endif

#ifndef A2BOHR
#define A2BOHR 0.52917721092
#endif

#ifndef NM2BOHR
#define NM2BOHR 0.052917721092
#endif


#ifndef HESSAU2GMX
#define HESSAU2GMX 9.37582E5
#endif


#ifndef AUT2PS
#define AUT2PS 2.418884326505000E-05
#endif


#ifndef MAXLINE
#define MAXLINE 1000
#endif

#ifndef MAXCHARINLINE
#define MAXCHARINLINE 1500
#endif

#ifndef CM2HARTREE
#define CM2HARTREE 1.000/219474.63
#endif


#ifndef HBAR
#define HBAR 1
#endif

#ifndef AMU2AU
#define AMU2AU 1822.88839
#endif



int main( int argc, char * argv[] )
{

/*-o-o-o-o-o-o-o some constants... o-o-o-o-o-o-o-o-*/

double done = 1.0000;

double dzero = 0.0000;


/*-o-o-o-o-o-o-o Variables... o-o-o-o-o-o-o-o-*/

register FILE * pg09hess , *pgmxhess ;

char g09hessname[ MAXCHARINLINE ] , gmxhessname[ MAXCHARINLINE ];

double * g09hess , * g09hess_tri , * gmxhess ;

int len_g09hessname , len_g09hessname_ext , len_gmxhessname ;

int len_g09hessFile ;


int natom , ncart ;

int icart , iatom , iconstant ;

int itmp; double dtmp; char chartmp[ 100 ];

char * keyword;

char buffer[ MAXCHARINLINE ] ;

char cache[ MAXCHARINLINE ] ;

char tmpString[ MAXCHARINLINE ] ;



char ** pcmd ; pcmd = argv ;

int icmd = 1 ;
  




/*-o-o-o-o-o-o-o Debugging... o-o-o-o-o-o-o-o-*/

FILE * debug;

//-o-o-o-o-o-o-o-o Recording cmd-line and time stamp o-o-o-o-o-o-o-o-o-//

time_t current_time;

time( &current_time );

char now[ 300 ] ;

strcpy( now , ctime( &current_time ) );

int lennow = strlen( now ) ;

*( now + lennow - 1 ) = ' ';

  
  
  printf("\n**********************************************************************\n");
    printf("* G_FCHKHESS2GMX_D : Stand-Alone Utility to Convert Hessian          *\n");
    printf("* From G09 .fchk file into gmx format and unit.                      *\n");
    printf("*                                                                    *\n");
    printf("*  ");
for( icmd = 0 ; icmd < argc ; icmd ++ )
{
  printf("%s " , *( pcmd + icmd ) );
}
printf("\n");
    printf("*                                                                    *\n");
    printf("*                                                                    *\n");
    printf("* Current Time : %s                           *\n" , now );
    printf("*                                                                    *\n");
    printf("*                                                                    *\n");
    printf("**********************************************************************\n");


/*-o-o-o-o-o-o-o-o Parsing the Command-Line Arguments o-o-o-o-o-o-o-o-o-*/

//char ** pv ; pv = argv ;

if( argc == 1 )
{
  printf("\nUsage : %s [ input g09 Hessian file name (unit = Hartree/Bohr^2 )] [ output Hessian file name (unit = kJ/mol/nm^2 ) ] \n\n\n" , *argv  );

  exit( 1 );
}
else if( argc == 2 )
{
  strcpy( g09hessname , *( pcmd + 1 ) );

  if( strcmp( g09hessname , "-h" ) == 0 )
  {
    printf("\nUsage : %s [ input g09 Formatted Checkpoint (.fchk .fch ) files containing hessian matrix (unit = Hartree/Bohr^2 )] [ output Hessian file name (unit = kJ/mol/nm^2 ) ] \n\n\n" , *argv  );
    
    exit( 3 ) ;
  
  }
  else
  {
    len_g09hessname = strlen( g09hessname ) ;

    if( *( g09hessname + len_g09hessname - 4 ) == '.' )
    
      len_g09hessname_ext = 3 ;
 
    else if( *( g09hessname + len_g09hessname - 5 ) == '.' )
    
      len_g09hessname_ext = 4 ;
 
    else if( *( g09hessname + len_g09hessname - 3 ) == '.' )

      len_g09hessname_ext = 2 ;

    else
    {	   
      printf("\nCan't you just name your file regularly???\n");

      exit( 5 );
    }

    strncpy( gmxhessname , g09hessname , len_g09hessname - len_g09hessname_ext );
  
    *( gmxhessname + len_g09hessname - len_g09hessname_ext ) = '\0' ;

    strcat( gmxhessname , "cmw" );
  
  }

}
else if( argc == 3 )
{
  printf("\nGood ... you have both input and output file names ... \n");

  strcpy( g09hessname , *( pcmd + 1 ) );

  strcpy( gmxhessname , *( pcmd + 2 ) );
}
else
{
  printf("\nToo many command-line arguments ...\n\n");

  exit( 69 );

}


if( ( pg09hess = fopen( g09hessname , "r" ) ) == NULL )
{
  printf("\nUser specified input Gaussian Hessian file does not exist ... \n");

  exit( 6 );
}



fsearch( pg09hess , "atoms" ) ; fsearch( pg09hess , "I" ) ; 

fscanf( pg09hess , "%s" , cache ) ; 

natom = atoi( cache ) ;

ncart = 3 * natom ; 

 
len_g09hessFile = ncart * ( ncart + 1 ) / 2 ;

rewind( pg09hess ) ;



/*-o-o-o-o-o-o-o Allocating Space ... o-o-o-o-o-o-o-o-*/

g09hess_tri = ( double * ) calloc( len_g09hessFile , sizeof( double ) );

g09hess = ( double * ) calloc( ncart * ncart , sizeof( double ) );

gmxhess = ( double * ) calloc( ncart * ncart , sizeof( double ) );



/*-o-o-o-o-o-o-o ... Loading : Hessian, mass, cart_q0 ... o-o-o-o-o-o-o-o-*/

/* Loading ... */



rewind( pg09hess );

fsearch( pg09hess , "Constants" ) ; 

fsearch( pg09hess , "N=" ) ;

fscanf( pg09hess , "%s" , cache ) ; 

itmp = atoi( cache ) ;

if( itmp != len_g09hessFile )
{
  printf( "\nAccording to .fchk file, NAtom = [ %d ] , so there should be [ %d ] numbers in \"Cartesian Force Constants\" section\n" , natom , len_g09hessFile ) ;
  
  printf( "\nHowever, there are in total [ %d ] numbers in this section ... Something is wrong ...\n\n" , itmp ) ;
  
  exit( 37 ) ;

}

//fscanf( pg09hess , "%s" , chartmp );

//printf("\nNow we are at %s ... \n" , chartmp );

dtmp = 0.0000 ;

for( itmp = 0 ; itmp < len_g09hessFile ; itmp ++ )
{
  fscanf( pg09hess , "%lf" , &dtmp ) ;
  
  *( g09hess_tri + itmp ) = dtmp ;

}



/* Change triangular Hessian(Cart) into symmetric matrix  */

dtri2sym( 'L' , ncart , g09hess_tri , g09hess );


/*
debug = fopen( "g09hess.hess" , "wb+" );

doutput( debug , ncart , ncart , g09hess );

fclose( debug );
*/



/* Converting in progress ... */

for( itmp = 0 ; itmp < ncart * ncart ; itmp ++ )
{
  *( gmxhess + itmp ) = *( g09hess + itmp ) * HESSAU2GMX ;
}


pgmxhess = fopen( gmxhessname , "wb+" );

fprintf( pgmxhess , "GROMACS Hessian File , unit = kJ/mol/(nm^2) , Generated by %s : \n %d %d \n" , *argv , ncart , ncart );

doutput( pgmxhess , ncart , ncart , gmxhess );

fclose( pgmxhess );


/*-o-o-o-o-o-o-o ... Em...What else? ... o-o-o-o-o-o-o-o-*/



return( 0 );

}






