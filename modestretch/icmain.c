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

#ifndef KBAU
#define KBAU 3.1668107400994983e-06
#endif


int nbin , atompick ;

typedef struct groform
{
  //char resnumber[10] ;        
  int resnumber;
  char resname[10] ;
  char atomname[10] ;
  int atomnumber ;
  double cx ;
  double cy ;
  double cz ;
  double vx ;
  double vy ;
  double vz ;

} GRO ;


/*
double hcPsi( int n , double xp , double w , double A )
{
  double wavefunc = A * hermite( n , xp ) * exp( -0.50 * xp * xp ) ;
  printf("% 10.8E * hermite( %d , %lf ) * exp( -0.50 * %lf * %lf ) = % 10.8E\n\n" , A , n , xp , xp , xp , wavefunc ) ;
  printf("\nhc = % 12.8E\n" , wavefunc ) ;
  return( wavefunc ) ;

}


double gaussianDistribution( double sigma , double mu , double x )
{
  double g = 1.00 / sigma / sqrt( 2.00 * PI ) * exp( -1.00 * ( x - mu ) * ( x - mu ) / 2.00 / ( sigma * sigma ) ) ;
  
  return g ;

}

*/
 
double wignerPositionNormalize( int n , double w , double ZPEscalingFactor )
{
  double scaledOmega = ( n + 0.50 * ZPEscalingFactor ) / ( n + 0.50 ) * w ;
  double leftCutOff = -1.00 * sqrt( ( 2.0 * n + 1.00 ) / scaledOmega );
  int nstepOneSide = 1200 ;
  int nsteps = 2 * nstepOneSide ;
  double stepsize = -1.00 * leftCutOff / nstepOneSide ;
  
  int istep = 0 ;
  double accumulatedP = 0.00 ;
  double wavefunc = 0.00 ;
  
  //double scaledOmega = w ;
  
  for( istep = 0 ; istep < nsteps ; istep ++ )
  {
    wavefunc = hcPsi( n , ( leftCutOff + stepsize * istep ) * sqrt( scaledOmega ) , scaledOmega , 1.00 ) ;
    accumulatedP = accumulatedP + stepsize * wavefunc * wavefunc ;
    //printf("\n% 12.8f * % 12.8f * % 12.8f = % 12.8f\n" , stepsize , wavefunc , wavefunc , stepsize * wavefunc * wavefunc ) ;
  }
  
  return( 1.00 / accumulatedP ) ;

}

double wignerVelocityNormalize( int n , double w , double ZPEscalingFactor )
{
  double scaledOmega = ( n + 0.50 * ZPEscalingFactor ) / ( n + 0.50 ) * w ;
  double leftCutOff = -1.00 * sqrt( ( 2.0 * n + 1.0 ) * scaledOmega ) ;
  int nstepOneSide = 400 ;
  int nsteps = 2 * nstepOneSide ;
  double stepsize = -1.00 * leftCutOff / nstepOneSide ;
  
  
  int istep = 0 ;
  double accumulatedP = 0.00 ;
  double wavefunc = 0.00 ;
  
  //double scaledOmega = w ;
  
  for( istep = 0 ; istep < nsteps ; istep ++ )
  {
    wavefunc = hcPsi( n , ( leftCutOff + stepsize * istep ) / sqrt( scaledOmega ) , scaledOmega , 1.00 ) ;
    accumulatedP = accumulatedP + stepsize * wavefunc * wavefunc ;
  }
  
  return( 1.00 / accumulatedP ) ;

}




double canonicalPositionNormalize( int n , double w , double ZPEscalingFactor , double T )
{
  
  double scaledOmega = ( n + 0.50 * ZPEscalingFactor ) / ( n + 0.50 ) * w ;
  double ctp = sqrt( ( 2.0 * n + 1.00 ) / scaledOmega );
  double leftCutOff = -1.00 * ctp;

  int nstepOneSide = 1200 ;
  int nsteps = 2 * nstepOneSide ;
  double stepsize = -1.00 * leftCutOff / nstepOneSide ;
  
  double beta = 1.00 / ( KBAU * T ) ;
  double u = tanh( scaledOmega * beta / 2.00 ) ;
  double mu = 0.00 ;
  double sigma = sqrt( ( 1.00 / scaledOmega ) / ( 2.00 * u ) ) ;
  
  int istep = 0 ;
  double accumulatedP = 1E-16 ;
  double wavefunc = 0.00 ;
  
  for( istep = 0 ; istep < nsteps ; istep ++ )
  {
    wavefunc = gaussianDistribution( sigma , mu , ( leftCutOff + stepsize * istep ) ) ;
    //printf("\nAdded % 12.8E ... \n\n" , hc ) ;
    accumulatedP = accumulatedP + stepsize * wavefunc ;
    //printf("\nUp to istep = %d , accumulatedP = % 10.8E , difference towards half is % 16.12E... \n" , istep , accumulatedP , 0.50 - accumulatedP ) ;
  }
  
  return( 1.00 / accumulatedP ) ;
 
}

double canonicalVelocityNormalize( int n , double w , double ZPEscalingFactor , double T )
{
  double scaledOmega = ( n + 0.50 * ZPEscalingFactor ) / ( n + 0.50 ) * w ;
  double ctp = sqrt( ( 2.0 * n + 1.0 ) * scaledOmega ) ; 
  double leftCutOff = -1.00 * ctp ;

  int nstepOneSide = 400 ;
  int nsteps = 2 * nstepOneSide ;
  double stepsize = -1.00 * leftCutOff / nstepOneSide ;
  
  double beta = 1.00 / ( KBAU * T ) ;
  double u = tanh( scaledOmega * beta / 2.00 ) ;
  double mu = 0.00 ;
  double sigma = sqrt( ( 1.00 * scaledOmega ) / ( 2.00 * u ) ) ;
  
  int istep = 0 ;
  double accumulatedP = 1E-16  ;
  double wavefunc = 0.00 ;
  
  for( istep = 0 ; istep < nsteps ; istep ++ )
  {
    wavefunc = gaussianDistribution( sigma , mu , ( leftCutOff + stepsize * istep ) ) ;
    
    accumulatedP = accumulatedP + stepsize * wavefunc ;
    //printf("\nUp to istep = %d , accumulatedP = % 10.8E , difference towards half is % 16.12E... \n" , istep , accumulatedP , 0.50 - accumulatedP ) ;
  }
  
  return( 1.00 / accumulatedP ) ;

}





int main( int argc, char * argv[] )
{

  // =====> Declaring functions ... 
  void icgen( int natomselect , int * irstatus_list , long * seeda , long * seedn , double fixedPhaseAngle ,
            double nquanta , double irfront , double irtail , int asinitial , int asfinal , 
            double T , double ZPEscalingFactor , double * mass , double * vibfreq, 
            double * cart_q0, double * U , double * cart_q,  double * cart_p , 
            double * wignerPositionNormArray , double * wignerVelocityNormArray , 
            double * canonicalPositionNormArray , double * canonicalVelocityNormArray ,
            int debuggingMode ) ;
            
  // =====> Some Constants ...

  double done = 1.0000;

  double dzero = 0.0000;

  time_t current_time;

  time( &current_time );

  char now[ 300 ] ;

  strcpy( now , ctime( &current_time ) );

  int lennow = strlen( now ) ;

  *( now + lennow - 1 ) = ' ';

  // =====> Declaring variables ... 

  int max_cmd_args = 40 ;

  if( argc > max_cmd_args )
  {
    printf("\nToo many command-line arguments ... \n");

    exit( 6 );
  }

  int icmd ;

  char ** pcmd ; pcmd = argv ;
  
  char ** pcmd_tmp ;

  printf("\n****************************************************************\n");
    printf("*  G_NEQICGEN_MODE : Non-EQuilibrium Initial Condition GENerator. *\n");
    printf("*                                                              *\n");
  printf("*  ");
  for( icmd = 0 ; icmd < argc ; icmd ++ )
  {
    printf("%s " , *( pcmd + icmd ) );
  }
  printf("\n");
    printf("*                                                              *\n");
    printf("*                                                              *\n");
    printf("* Current Time : %s                     *\n" , now );
    printf("*                                                              *\n");
    printf("* Note : 1) icmain takes 3N modes and will only pass the 3N-6  *\n");
    printf("*                                                              *\n");
    printf("*           modes into icgen ;                                 *\n");
    printf("*                                                              *\n");
    printf("*        2) Inside of this program , everything works in AU.   *\n");
    printf("*                                                              *\n");
    printf("****************************************************************\n");




  FILE * pgro , * pfreq , * pdxdr , * pmass , * phess ;
  
  FILE * debug ;
  
  double * freq , * vibfreq , * dxdr , * vibdxdr , * mass , * hess ;
  
  //char * inpgroname , * freqname , * dxdrname , * massname , * hessname ;
  
  char inpgroname[ MAXCHARINLINE ] , freqname[ MAXCHARINLINE ] , dxdrname[ MAXCHARINLINE ] , massname[ MAXCHARINLINE ] , hessname[ MAXCHARINLINE ] , outnameprefix[ MAXCHARINLINE ] ;
  
  char frequnit[ 10 ] ;
  
  char hessunit[ 10 ] ;
  
  char massunit[ 10 ] ;
  
  char crdunit[ 10 ] ;
  
  char outunit[ 10 ] ;
  
  int outUnitType = 19 ; // 19 : gmx ; 21 : g09 ;
  
  int exgro , exfreq , exdxdr , exmass , exselect , exhess ;
  
  int exas , exoutatom ;
  
  int exPhaseAngle ; int opPhaseAngle = YES ;
  
  int debuggingMode = NO ;
  
  // opPhaseAngle = YES means the phase angle sampling will be performed and it is default ;
  // Also, when opPhaseAngle = YES , fixedPhaseAngle = 9.00 which is meaningless for sampling ;
  // The value parsed in from command-line argument for fixedPhaseAngle should be [0,2), but 
  // if other number is parsed in, it is also OK because the sin and cos function will take care
  // of this simple periodic issue ...
  
  double fixedPhaseAngle = 9.00 ;

  int natom , ncart , nmode , ntraj ;
  
  int natomoutput = 0 ;
  
  int natomselect , ncartselect , nmodeselect ;
  
  int randomtype , randomseed ;
  
  int startingFrame ;

  double T;
  
  double ZPEscalingFactor ;
  
  int excitationType ; // 0 for real-life broad-band ir excitation ; 1 for specific ir-excitation ;
 
  double nquanta ;

  char irmodelistfilename[ MAXCHARINLINE ];
  
  char activemodelistfilename[ MAXCHARINLINE ] ;
  
  double irfront , irtail ;
  
  int asinitial , asfinal ;

  int softrange ;

  double * cart_q, * cart_p, * cart_q0;


  int itmp , iatom , icart , imode , itraj ; 
  
  double dtmp ;
  
  char ctmp ;
  
  char buffer[ MAXCHARINLINE ] ;
  
  char cache[ MAXCHARINLINE ] ;
  
  char tmpString[ MAXCHARINLINE ] ;
  

  // =====> Parsing command line input arguments ...

  // -----> Defaults ...
  
  T = 300.00 ;
  
  randomtype = 9 ;
  
  startingFrame = 0 ;
  
  ntraj = 1000 ;
  
  nquanta = 1.00 ;
  
  irfront = 2000.000 * CM2HARTREE ;
  
  irtail  = 2100.000 * CM2HARTREE ;
  
  softrange = -100 ;
  
  ZPEscalingFactor = 1.0000 ;

  
  strcpy( dxdrname , "dxdr.deb" );
  
  strcpy( freqname , "vibfrequency.deb" );
  
  strcpy( hessname , "system.gess" ) ;
  
  strcpy( massname , "system.mass" ) ;
  
  strcpy( outnameprefix , "proic" );

  strcpy( frequnit , "au" ) ;
  
  strcpy( massunit , "amu") ;
  
  strcpy( hessunit , "gmx" ) ;
  
  strcpy( outunit , "gmx" ) ;

  // -----> Command-Line arguments ...
  
  exgro = 18 ; exdxdr = 20 ; exfreq = 22 ; exmass = 24 ; exselect = 26 , exhess = 28 ;
  
  exas = 32 ; exoutatom = 66 ;

  if( argc == 1 )
  {
    printf("\n\nNo command-line arguments provided ... Mission aborting ...\n\n");
    
    printf("\nPlease refer to the usage by typing ' %s -h '\n\n" , * argv );
    
    exit(1); 

  }

  char * flag ;
  
  icmd = 1 ;

  while( icmd < argc )
  {  
    pcmd ++ ; 

    flag = * pcmd ;

    printf("\nNo.%d argument , Currently @ flag = %s ...\n\n" , icmd , flag );

    if( ( * flag == '-' ) && ( strlen( flag ) == 2 ) )
    {
      switch ( *( flag + 1 ) )
      {

	      case 'c' : strcpy( tmpString , *( ++ pcmd ) );
	      
	                 strcpy( crdunit , *( ++ pcmd ) );
	                 
	                 strcpy( inpgroname , tmpString ) ;
	                   
	                 if( * crdunit == '-' )
	                 {
	                   printf("\nHey, you did not specify the unit of you Coordinate file ... !\n");
	                   
	                   exit( 11 ) ;
	                 
	                 }
	                   
	                 printf("\nCommand-line argument indicates : Input Eq. Coordinate File name : %s . Format is %s ...\n" , inpgroname , crdunit ); 
	                 
	                 	           
	                 exgro = 19 ;	           
	                 	                 
	                 icmd = icmd + 3 ;
	                 
	                 break ;
	               
	                 
	                 
	      case 'f' : strcpy( dxdrname , *( ++ pcmd ) ) ; 
			 
			         printf("\nCommand-line argument indicates : Input dxdr File name : %s ...\n" , dxdrname ); 
	      
	                 icmd = icmd + 2 ; 
	                 
	                 exdxdr = 21 ;
	                 
	                 break ;

	      case 'w' : strcpy( freqname , *( ++ pcmd ) );
	      
	                 strcpy( frequnit , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Input frequency File name : %s in the unit of %s ...\n" , freqname , frequnit ); 
	                 
	                 if( * frequnit == '-' )
	                 {
	                   printf("\nHey, you did not specify the unit of you frequency file ... !\n");
	                   
	                   exit( 11 );
	                 }
	                 
	                 icmd = icmd + 3 ;
	                 
	                 exfreq = 23 ;
	                 
	                 break ;

	      case 'm' : strcpy( massname , *( ++ pcmd ) );
	      
	                 strcpy( massunit , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Input mass File name : %s in the unit of %s ...\n" , massname , massunit ); 

	                 if( * massunit == '-' )
	                 {
	                   printf("\nHey, you did not specify the unit of you mass file ... !\n");
	                   
	                   exit( 11 );
	                 }

	                 icmd = icmd + 3 ;
	                 
	                 exmass = 25 ;
	                 
	                 break ;
	      
	      case 'N' : strcpy( outnameprefix , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Output file names will have a prefix : %s ...\n" , outnameprefix ); 
	                 
	                 strcpy( outunit , *( ++ pcmd ) ) ;
	                 
	                 if( * outunit == '-' )
	                 {
	                   pcmd -- ;
	                   
	                   printf("\nSince you did not specify the unit of you Output File , by default we use \"nm\" as the length unit !\n");
	                   
	                   strcpy( outunit , "gmx" ) ;
	                   
	                   outUnitType = 19 ;
	                   
	                   icmd = icmd + 2 ;
	                 }
	                 else
	                 {
	                   printf("\nCommand-line argument indicates : Output file will have unit type : %s ...\n" , outunit ); 
	                   
	                   icmd = icmd + 3 ;
	                 }
	                 
	                 break ;

	      case 'H' : strcpy( hessname , *( ++ pcmd ) );
	      
	                 strcpy( hessunit , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Input Hessian File name : %s . Unit is %s ...\n" , hessname , hessunit ); 
	                 
	                 if( * hessunit == '-' )
	                 {
	                   printf("\nHey, you did not specify the unit of you Hessian file ... !\n");
	                   
	                   exit( 11 );
	                 }
	                 	                 
	                 icmd = icmd + 3 ;
	                 
	                 exhess = 29 ;
	                 
	                 break ;
	                 
	      case 'T' : T = atof( *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Simulation Temperature is %lf ... \n" , T ); 
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;
	      
	      case 'z' : ZPEscalingFactor = atof( *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : ZPE will be scaled by %lf ... \n" , ZPEscalingFactor ); 
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;
	      
	      case 'R' : randomtype = atoi( *( ++ pcmd ) );
	                 
	                 printf("\nCommand-line argument indicates : Random number seeds are of type %d ... ( 9 for fixed number , 11 for time_NULL)\n" , randomtype ); 
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;


	      case 'S' : startingFrame = atoi( *( ++ pcmd ) );
	                 
	                 printf("\nCommand-line argument indicates : Starting Frame is %d ... \n" , startingFrame ); 
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;


	      case 'i' : pcmd_tmp = pcmd ;
	      
	                 itmp = strcmp( *( pcmd_tmp + 1 ) , "list" ) ;
	      
	                 if( itmp == 0 )
	                 {
	                   excitationType = 1 ;
	                   
	                   pcmd ++ ;
	                   
	                   strcpy( irmodelistfilename , * ( ++ pcmd ) );
	                   
	                   printf("\nCommand-line argument indicates : File containing IR-Excitation list is name ===> %s <=== \n" , irmodelistfilename ) ;
	                   
	                 }
	                 else
	                 {
	                   excitationType = 0 ;
	                   
	                   irfront = atof( *( ++ pcmd ) ) * CM2HARTREE ;
	      
	                   irtail = atof( *( ++ pcmd ) ) * CM2HARTREE ;
	      
	                   printf("\nCommand-line argument indicates : Modes within the range %lf - %lf cm^-1 are excited by IR pulse ... \n" , irfront / CM2HARTREE , irtail / CM2HARTREE ); 
	                 
	                 }
	                 
	                 icmd = icmd + 3 ;
	                 
	                 break ;



	      case 'r' : strcpy( cache , *( ++ pcmd ) ) ;
	      
	                 printf("\nCache is [ %s ] \n" , cache ) ;
	      
	                 icmd = icmd + 2 ;
	      
	                 if( strcmp( cache , "start" ) == 0 )
	                 {
	                   asinitial = 1 ;
	                 }
	                 else if( strcmp( cache , "ir" ) == 0 )
	                 {
	                   asinitial = -2 ;
	                 }
	                 else if( strcmp( cache , "list" ) == 0 )
	                 {
	                   asinitial = -3 ;
	                   
	                   asfinal = -3 ;
	                   
	                   exas = 46 ;
	                   
	                 }
	                 else if( *( cache + 0 ) == '-' || icmd > argc || *( cache + 0 ) > '9' || *( cache + 0 ) < '0' )
	                 {
	                   printf("\nError : Wrong Normal Mode Active Space Format. When specifying the active NM space, use integers only !\n");
	                   
	                   exit( 124 );
	                 }
	                 else
	                 {
	                   asinitial = atoi( cache ) ;
	                 }
	      
	                 strcpy( cache , *( ++ pcmd ) ) ;
	                 
	                 printf("\nCache is [ %s ] \n" , cache ) ;
	                 
	                 icmd = icmd + 1 ;
	                 
	                 if( strcmp( cache , "end" ) == 0 )
	                 {
	                   exas = 34 ;
	                   
	                   asfinal = 1 ; // Temporarily , because there is the "if( asinitial > asfinal )" below ... 
	                 }
	                 else if( strcmp( cache , "ir" ) == 0 )
	                 {
	                   exas = 38 ;
	                   
	                   asfinal = -2 ;
	                 }
	                 else if( *( cache + 0 ) == '-' || icmd > argc || *( cache + 0 ) > '9' || *( cache + 0 ) < '0' )
	                 {
	                   if( asinitial != -3 )
	                   {
	                     //asfinal = atoi( cache ) ;
	                     //exas = 33 ;
	                     //printf("\nCommand-line argument indicates : Active Normal Mode Space is between [ # %d ] and [ # %d ] \n" , asinitial , asfinal );
	                   
	                     printf("\nError : Command Line Error. Both \"Starting-From\" and \"End-At\" normal mode sequence number need to be provided!!!\n\n");
	                     
	                     exit( 125 ) ;
	                   }
	                   else
	                   {
	                     strcpy( activemodelistfilename , cache );
	                     
	                     asfinal = -3 ;
	                   }

	                 }
	                 else
	                 {
	                   asfinal = atoi( cache ) ;
	                   
	                   exas = 33 ;
	                   
	                   printf("\nCommand-line argument indicates : Active Normal Mode Space is between [ # %d ] and [ # %d ] \n" , asinitial , asfinal );
	                 }
	                 
	                 
	                 
	                 if( asinitial == -2 && asfinal != -2 )
	                 {
	                   exas = 37 ;
	                 }
	                 else if( asinitial != -2 && asfinal == -2 )
	                 {
	                   exas = 38 ;
	                 }
	                 else if( asinitial == -2 && asfinal == -2 )
	                 {
	                   exas = 39 ;
	                 }
	                 else if( asinitial == -3 && asfinal == -3 )
	                 {
	                   exas = 46 ;
	                   
	                   printf("\nActive space will be specified by file %s ...\n\n" , activemodelistfilename ) ;
	                 }
	                 else
	                 {
	                   if( asinitial > asfinal )
	                   {
	                     itmp = asinitial ;
	                   
	                     asinitial = asfinal ; 
	                   
	                     asfinal = itmp ;
	                   }
                 
	                 }

	                 
	                 //pcmd ++ ;
	                 
	                 break ;


	      case 't' : ntraj = atoi( *( ++ pcmd ) ) ;
	      
	                 printf("\nCommand-line argument indicates : # of trajectory of NE-MD is %d ... \n" , ntraj ); 
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;
	     
	     
	      case 'a' : strcpy( tmpString , *( ++ pcmd ) );
	      
                     if( strcmp( tmpString , "YES" ) == 0 || strcmp( tmpString , "yes" ) == 0 || strcmp( tmpString , "Yes" ) == 0 )
                     {
                       opPhaseAngle = YES ;
                       
                       printf("\nCommand-line argument indicates : Sampling of phase angle will be performed ...\n" );

                       fixedPhaseAngle = 9.00 ;

                     }
                     else if( strcmp( tmpString , "NO" ) == 0 || strcmp( tmpString , "no" ) == 0 || strcmp( tmpString , "No" ) == 0 )
                     {
                       opPhaseAngle = NO ;
                       
                       printf("\nCommand-line argument indicates : Sampling of phase angle will NOT be performed ...\n" );
                       
                       printf("\nSince no specific partition method is provided, energy will be equally distributed into Ep and Ek ...\n");

                       fixedPhaseAngle = 0.25 ;
                          
                     }
                     else if( *( tmpString ) <= '9' && *( tmpString ) >= '0' )
                     {
                       opPhaseAngle = NO ;

                       fixedPhaseAngle = atof( tmpString ) ;

                       printf("\nCommand-line argument indicates : Sampling of phase angle will NOT be performed ...\n" );
                       
                       printf("\nAlso, the partition of total energy will be based on user-provided number [ %lf ] ...\n\n" , fixedPhaseAngle ) ;

                     }
                     else
                     {
                       printf("\nInvalid choice of phase angle operation mode ... Aborting ...\n") ;

                       exit( 409 ) ;

                     }
	                 
	                 	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;
	      
	      
	      case 'q' : nquanta = atof( *( ++ pcmd ) ) ;
	      
	                 printf("\nCommand-line argument indicates : The vib-state IR excite modes to in NE-MD is %lf ... \n" , nquanta ); 
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;


	      case 's' : strcpy( cache , *( ++ pcmd ) ) ;
	      
	                 if( strcmp( cache , "all" ) == 0 || strcmp( cache , "ALL" ) == 0 || strcmp( cache , "All" ) == 0 )
	                 {
	                   printf("\nCommand-line argument indicates : All available atoms are selected for generating NEIC ... \n" ); 
	                   
	                   exselect = 28 ;
	                 }
	                 else
	                 {
	                   natomselect = atoi( cache ) ;
	      
	                   printf("\nCommand-line argument indicates : The first %d atoms are selected for generating NEIC ... \n" , natomselect ); 
	                   
	                   exselect = 27 ;
	                   
	                 }
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;


	      case 'x' : strcpy( cache , *( ++ pcmd ) ) ;
	      
	                 if( strcmp( cache , "all" ) == 0 || strcmp( cache , "ALL" ) == 0 || strcmp( cache , "All" ) == 0 )
	                 {
	                   printf("\nCommand-line argument indicates : All available atoms are selected for output ... \n" ); 
	                   
	                   exoutatom = 68 ;
	                 }
	                 else
	                 {
	                   natomoutput = atoi( cache ) ;
	      
	                   printf("\nCommand-line argument indicates : The first %d atoms are selected for output... \n" , natomoutput ); 
	                   
	                   exoutatom = 67 ;
	                   
	                 }
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;


          case 'g' : strcpy( cache , *( ++ pcmd ) ) ;

                     printf("\nCommand-line argument indicates : Debugging Mode Invoking : %s ...\n" , cache );

                     if( strcmp( cache , "YES" ) == 0 || strcmp( cache , "Yes" ) == 0 || strcmp( cache , "yes" ) == 0 || strcmp( cache , "Y" ) == 0 || strcmp( cache , "y" ) == 0 )
                     {
                       debuggingMode = YES ;
                     }
                     else if( strcmp( cache , "NO" ) == 0 || strcmp( cache , "No" ) == 0 || strcmp( cache , "no" ) == 0 || strcmp( cache , "N" ) == 0 || strcmp( cache , "n" ) == 0 )
                     {
                       debuggingMode = NO ;
                     }
                     else
                     {
                       printf("\nInvalid choice of debugging mode ... Mission Aborted ...\n\n") ;

                       exit( 21 ) ;
                     }

                     icmd = icmd + 2 ;

                     break ;


	      case 'h' : printf("\nUsage:  %s [ -f 'input dxdr file name' ] [ -c 'input coordinate file name' 'file type ( gro/g09inp/crd.a/crd.nm)' ] \n  \
	              [ -w 'input vibrational frequency file' 'unit : au or cm' ] [ -m 'input mass file' 'unit of mass = au / amu'] \n  \
	              [ -N 'Output file name prefix'] [ -H 'Hessian file name' ' unit : au or gmx' ] [ -T Simulation temperature ] [ -z Scaling Factor of ZPE ] [ -t # of Trajectory ] \n\
	              [ -s # of atom selected ] [ -x # of atom selected for output ] [ -q to which level IR excite active modes to ]\n\
	              [ -i ir-excitation range ( See Explanation below ) ] [ -R random seed type ( 9 for fixed number , 11 for time_NULL )] [ -r initial and final mode # of active space ]\n\
	              [ -a option for whether sample the phase angle. Use floating number to specify the fixed phase angle. ][ -g debugging mode ]\n\n" , * argv ); 
	                 
	                 printf("\nFor -i flag, available input format are : \n\t(1) Two floating number seperated by space key to indicate IR range front and tail ...\n");
	                 
	                 printf("\n\t(2) -i list [ listname ] indicate all the IR-Excited modes are listed in the file \"listname\" ...Empty file indicates EqMD situation.\n");
	                 
	                 printf("\n\t(3) -i -1.00 -1.00 also indicates EqMD situation.\n") ;
	                 
	                 printf("\nFor -r flag , \"ir\" \"ir\" will cause the active space to include only ir-Excited modes ...\n\n" ) ;
	                 
	                 printf("\nFor -r flag , \"list\" \"[ active list file name ]\" will cause the active space to be specified by a stand-along file ...\n\n" ) ;
	                 
	                 printf("\nFor -s flag, \"All\" or \"ALL\" or \"all\" will select all available atoms. Default value is also choosing all available atoms ...\n\n") ;
	                 
	                 printf("\nFor -x flag, \"All\" or \"ALL\" or \"all\" will select all available atoms for output. Default value is also choosing all atoms which were selected for NEqIC generation ...\n\n") ;
	      
	                 printf("\nFor -a flag, YES/Yes/yes will cause the phase-angle sampling to be performed (default); NO/No/no will cause the phase-angle to be fixed at 0.25, meaning Ep=Ek; Floating number can be used to specify the partition. \ne.g. \"0.3\" being specified will make the fixed phase angle 0.3*pi ...\n\n") ;
	                 
	                 
	                 

	                 icmd = icmd + 1 ; 
	                 
	                 exit(1) ;

	      default : printf("\n\nInvalid option ' %s ' ... Please refer to the usage by typing ' %s -h '\n\n" , flag , * argv ); 
	      
	                icmd = argc ; 
	                
	                exit(1);

      
      
      
      }
    
    }
    else
    {
        printf("\n\nInvalid option ' %s ' ... Please refer to the usage by typing ' %s -h '\n\n" , flag , * argv );

	    exit(1);
      
      
    }
    
 
  } 
  

  // =====> Checking up on somethings that can cause mission aborted immediately ...
  /*
  if( softrange >= irfront )
  {
    printf("\nReally? While modes with frequency below % 12.8f will be treated as soft-modes per your request, you also said modes between % 12.8f and % 12.8f will be IR-Excited ...\n" , softrange , irfront , irtail );
  
    printf("\nWe can procede for now, but please make sure your input parameters are correct after the run ... \n\n");
    
    //exit( 72 );
  
  } 
  */
  
  
  
  
  


  // =====> pre-Loading .gro file to get NAtom info ...  
  
  char grotitlestring[MAXLINE];
  
  int iline = 3 ;
  
  int iload = 0 ;
  
  int blank_signal , groinfo , info ;

  char tmp_char ;
  
  int natomgroline , natomgrotitle ;
  
  int natomInput ; // only for the g09inp type use ...
  
  int unittype = 1 ;
  
  // 1 ... .gro ; 3 ... G09 input ; 5 ... Generic crd in angstrom ( crd.a ); 7 ... Generic crd in nm ( crd.nm ) ;
  

  if( ( pgro = fopen( inpgroname , "r" ) ) == NULL )
  {
    printf("\nUser defined .gro file does not exist ... \n");
   
    exit( 3 );
  }

  
  
  if( strcmp( outunit , "gmx" ) == 0 )
  {
    outUnitType = 19 ; 
  }
  else if( strcmp( outunit , "g09" ) == 0 )
  {
    outUnitType = 21 ;
  }
  else
  {
    printf("\nERROR : Unknown Unit Type for Output File ... Exiting ...\n\n") ;
    
    exit( 21 ) ;
  }
  
  
  if( strcmp( crdunit , "gro" ) == 0 )
  {
    unittype = 1 ;
    
    rewind( pgro );
    
    printf("\nCurrent character is %c ... \n" , fgetc( pgro ) );
    
    fskip( pgro , 1 );

    fscanf( pgro , "%d" , &natomgrotitle );
    
    fskip( pgro , 1 );
    
    printf("\n Second line of .gro file says it is describing %d atoms ... \n\n" , natomgrotitle );
    
  
    
 
    printf("\nNow let's read the actual .gro file ans see how many atoms it is describing ... \n");

    while( ( groinfo = freadline( buffer , MAXCHARINLINE , pgro , ';' ) ) != 0 )
    { 
      //printf("\n//--------------> WORKING ON NO. %d LINE ... <-------------//\n" , iline );
    
      //info = freadline( buffer , MAXCHARINLINE , pinputfile , ';' );
    
      //printf("\nNow we are at : %s\n" , buffer );
    
      //printf("\n INFO is %d ...\n" , info );
    
      blank_signal = stellblank( buffer ) ;
    
      if( blank_signal == 0 )
      {
        //printf("\nNo.%d line is a blank line ... Moving on ...\n" , iline ) ;
      }
      else if( blank_signal == 1 )
      {  
        //printf("\nNo.%d line is NOT a blank line ... loading ...\n" , iline );
        
        if( ( tmp_char = getfirst( buffer ) ) == ';' )
        {
          //printf("\nThis is a comment line ... So nothing will be loaded ...\n");
          
          //fskip( pinputfile , 1 );
        }
        else
        {
          //printf("\nLine reads : %s ...\n" , buffer );
          
          iload ++ ;
        
        }
        //printf("\n%s\n" , buffer );
      }
      else
      {
        printf("\nSomething is wrong with the reading file part ...\n");
      
        exit(1);
      }
    
      iline ++ ;
  
    }//while( info != 0 );

    natomgroline = iload - 1 ;
  
    printf("\nIt can be seen that this .gro file is describing %d atoms ...\n" , natomgroline );
  
    if( natomgrotitle > natomgroline )
    {
      printf("\nYour .gro file is self-contradictory ... While the second line of your .gro file says there will be %d atoms, there are actually only %d atoms being described ... \n" ,  natomgrotitle , natomgroline );
  
      printf("\nWe will take all the atoms we can to procede ... \n");
    
      natom = natomgroline ;
    }
    else if( natomgrotitle == natomgroline )
    {
      printf("\nOkay ... Your .gro file is fine ... NAtom will be %d ... \n" , natomgrotitle );
    
      natom = natomgrotitle ;
    }
    else if( natomgrotitle < natomgroline )
    {
      printf("\nThe second line of your .gro file indicates there are %d atoms in this file but there are more atoms ( %d atoms ) being described insided ... We will take the first %d atoms ... \n" , natomgrotitle , natomgroline , natomgrotitle );
  
      natom = natomgrotitle ;
    }
    else
    {
      printf("\nSomething is wrong with checking the .gro file ... Please take a look at it ... \n");
    
      exit( 81 );
    }
  
    rewind( pgro );

  }
  else if( strcmp( crdunit , "g09inp" ) == 0 )
  {
    unittype = 3 ;
    
    rewind( pgro );
    
    while( ( info = freadline( buffer , MAXCHARINLINE , pgro , '!' ) ) != 0 )
    {
      blank_signal = stellblank( buffer ) ;
 
      if( blank_signal != 0 )
      {
        tmp_char = getfirst( buffer ) ;
 
        printf("\nLine reads : %s \n\n" , buffer ) ;
 
        if( tmp_char == '#' )
        {
          break ;
        }
 
      }
 
    } // Searching for the 1st line of #Route section ...
 
 
    while( ( info = freadline( buffer , MAXCHARINLINE , pgro , '!' ) ) != 0 )
    {
      blank_signal = stellblank( buffer ) ;
 
      if( blank_signal == 0 )  break ;
 
    }  // Searching for the end of #Route section ...
 
 
    while( ( info = freadline( buffer , MAXCHARINLINE , pgro , '!' ) ) != 0 )
    {
      blank_signal = stellblank( buffer ) ;
 
      if( blank_signal == 0 )  break ;
 
    }  // Searching for the end of Comment card section ...
   
   
    fskip( pgro , 1 ) ;
    
    iload = 0 ; iline = 0 ;
 
    while( ( info = freadline( buffer , MAXCHARINLINE , pgro , '!' ) ) != 0 )
    {
      //printf("\n//--------------> WORKING ON NO. %d LINE ... <-------------//\n" , iline );
 
      //info = freadline( buffer , MAXCHARINLINE , pinputfile , ';' );
 
      //printf("\nNow we are at : %s\n" , buffer );
 
      //printf("\n INFO is %d ...\n" , info );
 
      blank_signal = stellblank( buffer ) ;
 
	  if( blank_signal == 0 )
      {
        //printf("\nNo.%d line is a blank line ... Moving on ...\n" , iline ) ;
 
        continue ;
      }
      else if( blank_signal == 1 )
      {
        //printf("\nNo.%d line is NOT a blank line ... loading ...\n" , iline );
 
        if( ( tmp_char = getfirst( buffer ) ) == '!' )
        {
          //printf("\nThis is a comment line ... So nothing will be loaded ...\n");
 
          //fskip( pinputfile , 1 );
 
          continue ;
        }
        else
        {
          iload ++ ;
 
        }
 
        //printf("\n%s\n" , buffer );
      }
      else
      {
        printf("\nSomething is wrong with the reading file part ...\n");
 
        exit( 77 );
      }
 
    }
 
    natomInput = iload ;
    
    natom = natomInput ;
 
    printf("\nYour input G09 .inp file provides geometry for %d atoms ...\n" , natomInput ) ;
  
  
  }
  else if( strcmp( crdunit , "crd.a" ) == 0 )
  {
    unittype = 5 ; 
    
    rewind( pgro );
    
    itmp = flength( pgro ) ;
    
    if( itmp % 3 != 0 )
    {
      printf( "\nYour file contains [ %d ] numbers and is not a multiple of 3 ... It cannot be a crd.a type coordinate file ...\n\n" , itmp ) ;
      
      exit( 79 ) ;
    
    }
    else
    {
      natom = itmp / 3 ; 
      
      printf("\nYour input crd file provides geometry for [ %d ] atoms ...\n" , natom ) ;
    
    }
  
  
  
  }
  else if( strcmp( crdunit , "crd.nm" ) == 0 )
  {
    unittype = 7 ; 
    
    rewind( pgro );
    
    itmp = flength( pgro ) ;
    
    if( itmp % 3 != 0 )
    {
      printf( "\nYour file contains [ %d ] numbers and is not a multiple of 3 ... It cannot be a crd.a type coordinate file ...\n\n" , itmp ) ;
      
      exit( 79 ) ;
    
    }
    else
    {
      natom = itmp / 3 ; 
      
      printf("\nYour input crd file provides geometry for [ %d ] atoms ...\n" , natom ) ;
    
    }
  
  
  
  }
  
  rewind( pgro );
  


  ncart = 3 * natom ;
 
  switch ( natom )
  {
    case 1 : nmode = 0 ; break;

    case 2 : nmode = 1 ; break;

    default : nmode = ncart - 0 ; break;
 
  }
  
  
  if( exselect == 28 )
  {
    natomselect = natom ;
  }
  else if( exselect == 26 )
  {
    natomselect = natom ;
  }
  else if( exselect == 27 && natomselect > natom )
  {
    printf("\nThere are only %d atoms in this system ... you cannot select more than that ... \n" , natom );
    
    /*
    if( natomgrotitle > natomgroline && natomselect <= natomgrotitle )
    {
      printf("\nAlthough ... the second line of your initial .gro file did indicate there were supposed to be %d atoms in system ... So go back and make sure what you are trying to do ... \n" , natomgrotitle );
    }
    else if( natomgrotitle < natomgroline && natomselect <= natomgroline )
    {
      printf("\nAlthough ... your initial .gro file did describe %d atoms in system ... So go back and make sure what you are trying to do ... \n" , natomgroline );
    }
    */
    
    exit( 78 );
  }
  else if( exselect == 27 && natomselect <= natom )
  {
    printf("\nYou have selected %d atoms for NEIC ... There are %d atoms en toto in this system ... \n" , natomselect , natom );
  
  }
  else
  {
    printf("\nSomething is wrong with the atom selection process ... NAtom = %d , NAtomSelect = %d ... \n" , natom , natomselect );
    
    exit( 78 );
  }
  
  
  ncartselect = 3 * natomselect ;
  
  switch ( natomselect )
  {
    case 1 : nmodeselect = 0 ; break;

    case 2 : nmodeselect = 1 ; break;

    default : nmodeselect = ncartselect - 6 ; break;
 
  }
  
  
  
  if( exoutatom == 66 )
  {
    natomoutput = natomselect ;
    
    printf("\nBy default , all selected atoms will be output into configuration files ...\n\n") ;
  
  }
  else if( exoutatom == 68 )
  {
    natomoutput = natomselect ;
    
    printf("\nPer user's request , all selected atoms will be output into configuration files ...\n\n") ;
    
  }
  else if( exoutatom == 67 )
  {
    if( natomoutput > natomselect ) // Will be dead ...
    {
      printf("\nERROR : We have only [ %d ] atoms selected for NEqIC-generation , hence it is not possible to output [ %d ] atoms ...\n\n" , natomselect , natomoutput ) ;
      
      exit( 92 ) ;
      
    }
    else
    {
      printf("\nConfirmation : [ %d ] atoms have been selected for output ...\n\n\n" , natomoutput ) ;
    }
  
  }
  
  
  

  switch( exas )
  {
    case 32 : asinitial = 1 ; // Human Label ...
    
              asfinal = nmodeselect ; // Human Label ...
              
              printf("\nBy default, all the normal modes [ # 1 ] - [ # %d ] are considered as \"active\" for sampling" , asfinal ) ;
              
              break ;
              
              
    case 34 : asfinal = nmodeselect ; // Human Label ...
    
              printf("\nPer user's request, normal modes [ # %d ] - [ # %d ] are considered as \"active\" for sampling" , asinitial , asfinal ) ;
              
              break ;
              
              
    case 33 : printf("\nAgain , Per user's request, normal modes [ # %d ] - [ # %d ] are considered as \"active\" for sampling" , asinitial , asfinal ) ;
    
              break ;
              
    
    case 37 : printf( "\nPer user's request, active space will start from the first IR-Mode ... and will end at [ # %d ] ...\n" , asfinal ) ;
    
              break ;
              
              
    case 38 : printf( "\nPer user's request, active space will start from the [ # %d ] ... and will end at the last IR-Mode ...\n" , asinitial ) ;
    
              break ;
              
    case 39 : printf( "\nPer user's request, active space will start from the first IR-Mode ... and will end at the last IR-Mode ...\n" ) ;
    
              printf( "\nAlso, only IR-Mode will be marked as \"active\" ... the range of \"active\" space will be determined later after IR range is determined ...\n\n" ) ;
              
              break ;
              
    case 46 : printf( "\nPer user's request, active space is specified by list file [ %s ] ...\n" , activemodelistfilename ) ;
    
              break ;
              
    default :  printf("\nUNKOWN Error : Active Space Flag Corrupted ... \n");
    
               exit( 126 ) ;
  
  
  }
  

  
  
  
  
  
  
  
  
  
  
 
  // =====> Allocating memories ... 

  mass = calloc( natomselect , sizeof( double ) ); dzeros( natomselect , 1 , mass );


  cart_q0 = calloc( ncartselect , sizeof( double ) ); dzeros( ncartselect , 1 , cart_q0 );


  freq = calloc( ncartselect , sizeof( double ) );  dzeros( ncartselect , 1 , freq );

  
  vibfreq = calloc( nmodeselect , sizeof( double ) );  dzeros( nmodeselect , 1 , vibfreq );


  dxdr = calloc( ncartselect * ncartselect , sizeof( double ) );  dzeros( ncartselect * ncartselect , 1 , dxdr );


  vibdxdr = calloc( ncartselect * nmodeselect , sizeof( double ) );  dzeros( ncartselect * nmodeselect , 1 , vibdxdr );


  cart_q = calloc( ncartselect * ntraj , sizeof(double) ); dzeros( ntraj * ncartselect , 1 , cart_q );

  
  cart_p = calloc( ncartselect * ntraj , sizeof(double) ); dzeros( ntraj * ncartselect , 1 , cart_p );



  // =====> Loading Mass File ... 
  
  if( ( pmass = fopen( massname ,"r" ) ) == NULL )
  {
    printf("\nUser defined mass info file name is %s ...\n\n" , massname );
    
    printf("\nFile containing mass information does not exist...\n");
  
    exit(1);
  }
  else
  {
    printf("\nUser defined mass info file name is %s ...\n\n" , massname );
  
  }
  
  itmp = flength( pmass ) ;
  
  rewind( pmass );


  if( itmp != natom && itmp != natomselect )
  {
    printf("\nSomething is wrong with the mass file ... There are %d atomic mass in the file while the total NAtom in this system is %d and %d atoms are selected ... \n" , itmp , natom , natomselect );
  
    exit(1);
  }
  else if( itmp == natomselect && natom == natomselect )
  {
    printf("\nOKay ... I see you provided the mass info for the selected atoms ... Although you are selecting the whole system ... \n");
  }
  else if( itmp == natomselect && natom != natomselect )
  {
    printf("\nOKay ... I see you provided the mass info for just the selected atoms ... \n");
  }
  else if( itmp == natom && natom != natomselect )
  {
    printf("\nOKay ... I see you provided the mass for all the atom in this system ... \n"); 
  }
  else if( itmp > natomselect && itmp < natom )
  {
    printf("\nOKay ... you provided %d numbers for mass. The total NAtom in system is %d while %d atoms are selected ...\n" , itmp , natom , natomselect );

    printf("\nSo ... I will assume the first %d numbers in your mass file correspond to the selected atoms ...\n" , natomselect );
  }
  else
  {
    printf("\nSomething is wrong with the mass file ... There are %d atomic mass in the file while the total NAtom in this system is %d and %d atoms are selected ... \n" , itmp , natom , natomselect );

    exit(1);

  }


  // =====> Deciding whether we need to perform the diagonalization of Hessian ... 

  int diag_type = 0 ;
  
  switch( exdxdr * exfreq * exhess )
  {
    case 12760 : printf("\nNo eigenvectors or frequencies provided ... \n");
    
                 printf("\nFound User-Provided Hessian file : %s  ... Will perform Hessian diagonalization to obtain eigenvectors and frequencies ... \n" , hessname );
                 
                 diag_type = 1 ;
                 
                 break ;
                 
    case 12320 : printf("\nNo eigenvectors or frequencies provided ... \n") ;
    
                 printf("\nUser-Provided Hessian file NOT FOUND ... Will NOT perform Hessian diagonalization ... Mission Aborting ... \n");
                 
                 exit( 11 ) ;
                 
                 break ;
                 
    case 12880 : printf("\nNo eigenvector file name specified ... will use the default name %s ... \n" , dxdrname );
    
                 printf("\nFrequency file specified : %s ... \n" , freqname );
                 
                 diag_type = 0 ;
                 
                 break ;
                 
    case 12936 : printf("\nNo frequency file name specified ... will use the default name %s ... \n" , freqname );
  
                 printf("\nEigenvector file specified : %s ... \n" , dxdrname );
                 
                 diag_type = 0 ;
                 
                 break ;
                 
    case 13524 : printf("\nBoth eigenvector file name and frequency name are specified : %s and %s ... \n" ,  dxdrname , freqname );
    
                 diag_type = 0 ;
                 
                 break ;
                 
    case 14007 : printf("\nBoth eigenvector file name and frequency name are specified : %s and %s ... \n" ,  dxdrname , freqname );
    
                 printf("\nYou also provided the name of Hessian file : %s ... \n" , hessname ) ;
                 
                 if( strcmp( dxdrname , "none" ) == 0  && strcmp( freqname , "none" ) == 0 )
                 {
                   printf("\nI see the name of both eigenvector file and frequency file are NONE ... Hessian diagonalization will be lauched ... \n");
                   
                   diag_type = 1 ;
                 }
                 else
                 {
                   printf("\nI am confused ... What are you trying to do ? For now, I will neglect the hessian file and use the eigenvector file and frequency file you specified ... \n");
                 
                   diag_type = 0 ;
                 }
                 
                 break ;
                 
    default    : printf("\nEither you only provided the eigenvector file name or only the frequency name ... \n");
    
                 printf("\nHowever, since you provided the name of Hessian file : %s ... I will try to diagonalize the Hessian ... \n" , hessname );
                 
                 diag_type = 1 ;
                 
                 break ;

  }
  
  
  int len_dxdr , len_freq ;
  
  double massconvert ;
  
  if( strcmp( massunit , "au" ) == 0 )
  {
    massconvert = 1.000 ;
  }
  else if( strcmp( massunit , "amu" ) == 0 )
  {
    massconvert = AMU2AU ;
  }
  else
  {
    printf("\nUNKNOWN ERROR : The unit of mass appears as [ %s ] ...\n\n" , massunit ) ;
    
    exit( 29 ) ;
  }
  
  
  if( diag_type == 1 )
  {
    printf("\n============== Launching the diagonalization machinary ... =============\n");
    
    dxdrfreq( massname , hessname , hessunit , natom , natomselect , dxdr , freq );
    
    printf("\n============== Diagonalization Done ... =============\n");
    
    rewind( pmass );
    
    //debug = fopen( "mass.dat" , "wb+" );
    
    for( iatom = 0 ; iatom < natomselect ; iatom ++ )
    {
      fscanf( pmass , "%lf" , mass + iatom );
      
      *( mass + iatom ) = *( mass + iatom ) * massconvert ;
      
      //fprintf( debug , "% 12.8f\n" , *( mass + iatom ) );
      
    }
    
    printf("\n============== Done with loading Mass and Generating Eigenvectors and Frequencies ... =============\n");
  
  }
  else
  {
    printf("\n============== Loading Mass , Eigenvectors and Frequencies ... =============\n");
    
    // -----> mass ...
    
    rewind( pmass );
    
    //debug = fopen( "mass.dat" , "wb+" );
    
    for( iatom = 0 ; iatom < natomselect ; iatom ++ )
    {
      fscanf( pmass , "%lf" , mass + iatom );
      
      *( mass + iatom ) = *( mass + iatom ) * massconvert ;
      
      //fprintf( debug , "% 12.8f\n" , *( mass + iatom ) );
      
    }
    
    //fclose( debug );
    
    // -----> dxdr ... 
    
    if( ( pdxdr = fopen( dxdrname ,"r" ) ) == NULL )
    {
      printf("\nUser defined eigenvector file %s was NOT found ...\n\n" , dxdrname );
  
      exit( 9 );
    }
    else
    {
      len_dxdr = flength( pdxdr );
      
      if( len_dxdr != ncartselect * ncartselect )
      {
        printf("\nThere are %d number in eigenvector file ... %d != %d * %d ... Mission aborting ... \n" , len_dxdr , len_dxdr , ncartselect , ncartselect );
      
        exit( 11 );
      }
      else
      {
        rewind( pdxdr );
      
        fload( pdxdr , dxdr );
              
      }
      
    }
  


    // -----> frequency ...
    
    double freqconvert ;

    if( strcmp( frequnit , "au" ) == 0 )
    {
      printf("\nThe unit these frequencies are in is ATOMIC UNIT ... \n");
      
      freqconvert = 1.0000 ;
    }
    else if( strcmp( frequnit , "cm" ) == 0 )
    {
      printf("\nThe unit these frequencies are in is CM^-1 ... \n");
      
      freqconvert = CM2HARTREE ;
    }
    else if( strcmp( frequnit , "none" ) == 0 )
    {
      printf("\nWhat do you mean 'none'? I assume it is 'au' ... \n");

      freqconvert = 1.0000 ;  
    }
    else
    {
      printf("\nWhat did you say about the unit of you frequency again ??? \n");

      exit( 107 );
    } 
    
    
    if( ( pfreq = fopen( freqname ,"r" ) ) == NULL )
    {
      printf("\nUser defined frequency file %s was NOT found ...\n\n" , freqname );
  
      exit( 9 );
    }
    else
    {
      len_freq = flength( pfreq );
      
      if( len_freq != ncartselect )
      {
        printf("\nThere are %d number in frequency file ... while there are %d modes ( 6-mode included ) ... Mission aborting ... \n" , len_freq , ncartselect );
      
        exit( 11 );
      }
      else
      {
        rewind( pfreq );
      
        fload( pfreq , freq );
        
        for( icart = 0 ; icart < ncartselect ; icart ++ )
        {
          *( freq + icart ) = *( freq + icart ) * freqconvert ;
          //*( freq + icart ) = *( freq + icart ) * CM2HARTREE ; // Converting cm into Hartree ; 
        }
        
      
      }
     
    }

    

  

    printf("\n============== Done with loading Mass , Eigenvectors and Frequencies ... =============\n");

  
  }

  // -----> Removing translational and rotational motion ... 
  
  
  for( icart = 0 ; icart < ncartselect ; icart ++ )
  {
    for( imode = 0 ; imode < nmodeselect ; imode ++ )
    {
      *( vibdxdr + icart * nmodeselect + imode ) = *( dxdr + icart * ncartselect + imode + 6 );
    }
  }

  
  if( debuggingMode == YES )
  {
    debug = fopen( "vibdxdr.deb" , "wb+" );
  
    doutput( debug , ncartselect , nmodeselect , vibdxdr ) ;
  
    fclose( debug );
  }
  
  for( icart = 6 ; icart < ncartselect ; icart ++ )
  {
    *( vibfreq + icart - 6 ) = *( freq + icart ) ;
  }
  
  if( debuggingMode == YES )
  {
    debug = fopen( "freq_vib.deb" , "wb+" );
  
    doutput( debug , nmodeselect , 1 , vibfreq ) ;
  
    fclose( debug );
  
  }

  // =====> Loading Equilibrium Configuration for solute , i.e. Loading .gro file ... Yikes ... 

  GRO atomlist[ natomselect ] ; 
  
  double * array_buffer = calloc( ncart , sizeof( double ) ) ;
  
  rewind( pgro );
  
  if( unittype == 1 )
  {
    fskip( pgro , 2 );
  
    for( iatom = 0 ; iatom < natomselect ; iatom ++ )
    {
      fscanf( pgro , "%5d%5s" , &atomlist[ iatom ].resnumber , atomlist[ iatom ].resname );

      //printf( "%s\t" , atomlist[ iatom ].resname );

      fscanf( pgro , "%s" , atomlist[ iatom ].atomname );

      //printf( "%s" , atomlist[ iatom ].atomname );

      fscanf( pgro , "%d" , &atomlist[ iatom ].atomnumber );

      //printf( "\nWorking on No. %d atom ...\n" , atomlist[ iatom ].atomnumber );

      fscanf( pgro , "%lf" , &atomlist[ iatom ].cx ); //printf("\n Cx is %lf ...\t" , atomlist[ iatom ].cx);

      fscanf( pgro , "%lf" , &atomlist[ iatom ].cy ); //printf("\n Cy is %lf ...\t" , atomlist[ iatom ].cy);
  
      fscanf( pgro , "%lf" , &atomlist[ iatom ].cz ); //printf("\n Cz is %lf ...\n\n" , atomlist[ iatom ].cz);
    
      //fscanf( pgro , "%lf" , &atomlist[ iatom ].vx ); //printf("\n Vx is %lf ...\t" , atomlist[ iatom ].cx);

      //fscanf( pgro , "%lf" , &atomlist[ iatom ].vy ); //printf("\n Vy is %lf ...\t" , atomlist[ iatom ].cy);
  
      //fscanf( pgro , "%lf" , &atomlist[ iatom ].vz ); //printf("\n Vz is %lf ...\n\n" , atomlist[ iatom ].cz);

  
    }
  
    rewind( pgro );
  
    fskip( pgro , natomgroline + 1 );

    double boxvector[3];

    fscanf( pgro , "%lf" , boxvector + 0 );

    fscanf( pgro , "%lf" , boxvector + 1 );

    fscanf( pgro , "%lf" , boxvector + 2 );

    rewind( pgro );
    
      
    for( iatom = 0 ; iatom < natomselect ; iatom ++ )
    {
      *( cart_q0 + iatom * 3 + 0 ) = atomlist[ iatom ].cx / NM2BOHR ; // Convert into Bohr ;
    
      *( cart_q0 + iatom * 3 + 1 ) = atomlist[ iatom ].cy / NM2BOHR ;
    
      *( cart_q0 + iatom * 3 + 2 ) = atomlist[ iatom ].cz / NM2BOHR ;
    
    }
  

  }
  else if( unittype == 3 )
  {
    while( ( info = freadline( buffer , MAXCHARINLINE , pgro , '!' ) ) != 0 )
    {
      blank_signal = stellblank( buffer ) ;
 
      if( blank_signal != 0 )
      {
        tmp_char = getfirst( buffer ) ;
 
        printf("\nLine reads : %s \n\n" , buffer ) ;
 
        if( tmp_char == '#' )
        {
          break ;
        }
 
      }
 
    } // Searching for the 1st line of #Route section ...
 
 
    while( ( info = freadline( buffer , MAXCHARINLINE , pgro , '!' ) ) != 0 )
    {
      blank_signal = stellblank( buffer ) ;
 
      if( blank_signal == 0 )  break ;
 
    }  // Searching for the end of #Route section ...
 
 
    while( ( info = freadline( buffer , MAXCHARINLINE , pgro , '!' ) ) != 0 )
    {
      blank_signal = stellblank( buffer ) ;
 
      if( blank_signal == 0 )  break ;
 
    }  // Searching for the end of Comment card section ...
   
   
    fskip( pgro , 1 ) ;
    
    
    iload = 0 ; iline = 0 ;
    
    while( ( info = freadline( buffer , MAXCHARINLINE , pgro , '!' ) ) != 0 )
    {
      //printf("\n//--------------> WORKING ON NO. %d LINE ... <-------------//\n" , iline );
 
      //info = freadline( buffer , MAXCHARINLINE , pinputfile , ';' );
 
      //printf("\nNow we are at : %s\n" , buffer );
 
      //printf("\n INFO is %d ...\n" , info );
 
      blank_signal = stellblank( buffer ) ;
 
      if( blank_signal == 0 )
      {
        //printf("\nNo.%d line is a blank line ... Moving on ...\n" , iline ) ;
 
        continue ;
      }
      else if( blank_signal == 1 )
      {
        //printf("\nNo.%d line is NOT a blank line ... loading ...\n" , iline );
  
        if( ( tmp_char = getfirst( buffer ) ) == '!' )
        {
          //printf("\nThis is a comment line ... So nothing will be loaded ...\n");
  
          //fskip( pinputfile , 1 );
  
          continue ;
        }
        else
        {
          strpickword( buffer , 1 , cache ) ;
  
          strcpy( atomlist[ iload ].atomname , cache ) ;
  
          printf("\nWe have an atom [ %s ] ...\n" , cache ) ;
  
          strpickword( buffer , 2 , cache ) ;
  
          atomlist[ iload ].cx = atof( cache ) ;
  
          strpickword( buffer , 3 , cache ) ;
  
          atomlist[ iload ].cy = atof( cache ) ;
  
          strpickword( buffer , 4 , cache ) ;
  
          atomlist[ iload ].cz = atof( cache ) ;
  
  
          atomlist[ iload ].atomnumber = 0 ;
  
  
          iload ++ ;
          
  
        }
  
        //printf("\n%s\n" , buffer );
      }
      else
      {
        printf("\nSomething is wrong with the reading file part ...\n");
  
        exit( 73 );
      }
 
 
      if( iload == natom ) break ;
  
      iline ++ ;
  
    }
    
    for( iatom = 0 ; iatom < natomselect ; iatom ++ )
    {
      *( cart_q0 + iatom * 3 + 0 ) = atomlist[ iatom ].cx / A2BOHR ; // Convert from Angstrom into Bohr ;
    
      *( cart_q0 + iatom * 3 + 1 ) = atomlist[ iatom ].cy / A2BOHR ;
    
      *( cart_q0 + iatom * 3 + 2 ) = atomlist[ iatom ].cz / A2BOHR ;
    
    }
  

 
  
  }
  else if( unittype == 5 )
  {
    rewind( pgro ) ;
    
    fload( pgro , array_buffer ) ;
    
    
    for( iatom = 0 ; iatom < natomselect ; iatom ++ )
    {
      *( cart_q0 + iatom * 3 + 0 ) = *( array_buffer + iatom * 3 + 0 ) / A2BOHR ; // Convert from Angstrom into Bohr ;
    
      *( cart_q0 + iatom * 3 + 1 ) = *( array_buffer + iatom * 3 + 1 ) / A2BOHR ;
    
      *( cart_q0 + iatom * 3 + 2 ) = *( array_buffer + iatom * 3 + 2 ) / A2BOHR ;
    
    }
  

  
  }
  else if( unittype == 7 )
  {
    rewind( pgro ) ;
    
    fload( pgro , array_buffer ) ;
    
    
    for( iatom = 0 ; iatom < natomselect ; iatom ++ )
    {
      *( cart_q0 + iatom * 3 + 0 ) = *( array_buffer + iatom * 3 + 0 ) / NM2BOHR ; // Convert from nm into Bohr ;
    
      *( cart_q0 + iatom * 3 + 1 ) = *( array_buffer + iatom * 3 + 1 ) / NM2BOHR ;
    
      *( cart_q0 + iatom * 3 + 2 ) = *( array_buffer + iatom * 3 + 2 ) / NM2BOHR ;
    
    }
  

  
  }
  
  
  
  
  if( debuggingMode == YES )
  {
    debug = fopen("equilibrium_geom_au.deb", "wb+");
    
    for( iatom = 0 ; iatom < natomselect ; iatom ++ ) 
    {
      fprintf( debug , "% 16.12f\t% 16.12f\t% 16.12f\n\n" , *( cart_q0 + iatom * 3 + 0 ) , *( cart_q0 + iatom * 3 + 1 ) , *( cart_q0 + iatom * 3 + 2 ) );
    }
  
    fclose(debug);
  
  
  
    debug = fopen("equilibrium_geom_gmx.deb", "wb+");
  
    for( iatom = 0 ; iatom < natomselect ; iatom ++ )
    {
      fprintf( debug , "% 16.12f\t% 16.12f\t% 16.12f\n\n" , *( cart_q0 + iatom * 3 + 0 ) / 10.00 , *( cart_q0 + iatom * 3 + 1 ) / 10.00 , *( cart_q0 + iatom * 3 + 2 ) / 10.00 );
    }

    fclose(debug);

  }



  // =====> Loading Done... Now Random Numbers ... 

  long * seeda, * seedn ;

  seeda = calloc( nmodeselect , sizeof( long ) ); 

  seedn = calloc( nmodeselect , sizeof( long ) ); 

  pid_t pid ;
  
  pid = getpid() ;
  
  long tmp_seeda , tmp_seedn ;
  
  if( randomtype == 9 )
  {
    for( imode = 0 ; imode < nmodeselect ; imode ++ )
    {
      *( seeda + imode ) = 19880708 + imode * 99999 ; 
    
      * (seedn + imode ) = 19880123 + ( nmodeselect - imode ) * 99999 ; 
    }
  
  }
  else if( randomtype == 11 )
  {
    tmp_seeda = time( NULL ) + pid ;
    
    tmp_seedn = 11 * pid ;
    
    for( imode = 0 ; imode < nmodeselect ; imode ++ )
    {
      *( seeda + imode ) = ranzm( &tmp_seeda );
      
      *( seedn + imode ) = ranzm( &tmp_seedn );
      
    }
  
  }
  else
  {
    randomseed = randomtype ;
    
    tmp_seeda = randomseed ;
    
    tmp_seedn = ranzm( &tmp_seeda ) ;
    
    tmp_seeda = randomseed ;
    
    for( imode = 0 ; imode < nmodeselect ; imode ++ )
    {
      *( seeda + imode ) = ranzm( &tmp_seeda );
      
      *( seedn + imode ) = ranzm( &tmp_seedn );
      
    }
  
  }
  /*
  else
  {
    printf("\nInvalid choice of random seed type ... Mission aborting ...\n");
    
    exit( 11 );
  }
  */

  
  
  
  if( debuggingMode == YES )
  {
    debug = fopen( "seed.theta.0.deb" , "wb+" );
  
    loutput( debug , nmodeselect , 1 , seeda );
  
    fclose( debug );
  
  
    debug = fopen( "seed.level.0.deb" , "wb+" );
  
    loutput( debug , nmodeselect , 1 , seedn );
  
    fclose(debug);
  
  }
  
  
  
  // =====> Making IR-Mode-List ...
  
  int * irmodelist ;
  
  int NofIRModes = 0 ;
  
  FILE * pirmodelistfile ; 
  
  switch ( excitationType )
  {
    case 0 : printf("\nIR-Excitation range is from %lf to %lf ... \n" , irfront , irtail );
    
             break ;
    
    case 1 : if( ( pirmodelistfile = fopen( irmodelistfilename , "r" ) ) == NULL )
  	         {
  	           printf("\nUser specified IR-Excitation List not found ... Mission Aborting ... \n") ;
  	                     
  	           exit( 57 );
  	       
  	         }
  	       
  	         NofIRModes = flength( pirmodelistfile );
  	       
  	         if( NofIRModes != 0 )
  	         {
  	           rewind( pirmodelistfile );
  	       
  	           irmodelist = calloc( NofIRModes , sizeof( int ) );
  	       
  	           int_fload( pirmodelistfile , irmodelist );
  	       
  	           printf("\nAccording to user-provided file there are %d mode(s) which will be IR-Excited ...\n" , NofIRModes );
  	       
  	           printf("\nThey are ... : \n");
  	       
  	           for( itmp = 0 ; itmp < NofIRModes ; itmp ++ )
  	           {
  	             printf("\nMode # %d ...\n" , *( irmodelist + itmp ) );
  	           }
  	         } 
  	         else
  	         {
  	           printf("\nAlllll right ... The ir-excitation list you provided is empty ...  so I assume this is EqMD situation ...\n");
  	         
  	           printf("\nPlease make sure this is correct ... Proceeding now ...\n");
  	       
  	         }
  	       
  	         break ;
  
    
  }
  
  
  //printf("\n NAtomSelect = %d \n\n", natomselect);
  
  



  // =====> Making IR-Status-List Based on IR-Mode-List ...
 
  int * irstatus_list = calloc( nmodeselect , sizeof( int ) );
  
  izeros( nmodeselect , 1 , irstatus_list ) ;
  
  int irstatus ;
  
  int firstIRMode = -1 , lastIRMode = -1 ;

  
  for( imode = 0 ; imode < nmodeselect ; imode ++ )
  {
    if( excitationType == 0 )
    {
      irstatus = tellIR_range( *( vibfreq + imode ) , irfront , irtail );
          
    }
    else if( excitationType == 1 )
    {
      if( NofIRModes == 0 )
      {
        irstatus = 0 ;
      }
      else
      {
        irstatus = tellIR_list( imode + 1 , irmodelist , NofIRModes );
      }  
    
    }
    else
    {
      printf("\nWTF is going on? The type of IR excitation is %d ...\n" , excitationType );
        
      exit( 9 );
    }
    
    
    *( irstatus_list + imode ) = irstatus ;
    
    
    if( irstatus == 1 && firstIRMode == -1 )
    {
      firstIRMode = imode + 1 ; // Human Label ...
    }
    
    if( irstatus == 1 )
    {
      lastIRMode = imode + 1 ; // Human Label ... 
    }
    
    
    
  }
  
  
  int * activemodelist ;
  
  int NofActiveModes = 0 ;
  
  FILE * pactivemodelistfile ; 
  
  
  
  
  switch( exas )
  {
    case 37 : asinitial = firstIRMode ;
    
              break ;
              
    case 38 : asfinal = lastIRMode ;
    
              break ;
              
    case 39 : asinitial = firstIRMode ;
    
              asfinal = lastIRMode ;
              
              break ;
    
    case 32 : printf("\nAgain ... By default, all the normal modes [ # 1 ] - [ # %d ] are considered as \"active\" for sampling" , asfinal ) ;
              
              break ;
              
              
    case 34 : printf("\nPer user's request, normal modes [ # %d ] - [ # %d ] are considered as \"active\" for sampling" , asinitial , asfinal ) ;
              
              break ;
              
              
    case 33 : printf("\nAgain , Per user's request, normal modes [ # %d ] - [ # %d ] are considered as \"active\" for sampling" , asinitial , asfinal ) ;
    
              break ;
              
    case 46 : if( ( pactivemodelistfile = fopen( activemodelistfilename , "r" ) ) == NULL )
  	          {
  	            printf("\nUser specified Active Mode List not found ... Mission Aborting ... \n") ;
  	                     
  	            exit( 57 );
  	       
  	          }
  	       
  	          NofActiveModes = flength( pactivemodelistfile );
  	       
  	          if( NofActiveModes != 0 )
  	          {
  	            rewind( pactivemodelistfile );
  	        
  	            activemodelist = calloc( NofActiveModes , sizeof( int ) );
  	       
  	            int_fload( pactivemodelistfile , activemodelist );
  	       
  	            printf("\nAccording to user-provided file there are %d mode(s) which will be included in Active Space ...\n" , NofActiveModes );
  	       
  	            printf("\nThey are ... : \n");
  	       
  	            for( itmp = 0 ; itmp < NofActiveModes ; itmp ++ )
  	            {
  	              printf("\nMode # %d ...\n" , *( activemodelist + itmp ) );
  	            }
  	          
  	          } 
  	          else
  	          {
  	            printf("\nAlllll right ... The ir-excitation list you provided is empty ...  so I assume all modes are Active ...\n");
  	         
  	            printf("\nPlease make sure this is correct ... Proceeding now ...\n");
  	            
  	            NofActiveModes = nmodeselect ;
  	            
  	            activemodelist = calloc( NofActiveModes , sizeof( int ) );
  	            
  	            for( itmp = 0 ; itmp < nmodeselect ; itmp ++ )
  	            {
  	              *( activemodelist + itmp ) = itmp + 1 ;
  	            }
  	       
  	          }
  	          
  	          asinitial = imin( NofActiveModes , activemodelist ) ;
  	          
  	          asfinal = imax( NofActiveModes , activemodelist ) ;
  	       
  	          break ;
              
    
    default :  printf("\nUNKOWN Error : Active Space Flag Corrupted ... \n");
    
               exit( 126 ) ;
  }
  
  
  
  
  for( imode = 0 ; imode < nmodeselect ; imode ++ )
  {
    irstatus = *( irstatus_list + imode ) ;
    
    if( exas == 46 ) 
    {
      irstatus = tellActive_list( imode + 1 , activemodelist , NofActiveModes );
      
    }
    else if( exas == 39 )
    {
      if( irstatus != 1 )
      {
        irstatus = -1 ;
      }
    
    }
    else
    {
      if( imode + 1 > asfinal || imode + 1 < asinitial )
      {
        irstatus = -1 ;
      }
    
    }
    

    *( irstatus_list + imode ) = irstatus ;
    
  }
  
  
  


  
  //icgen( natomselect , ntraj , randomtype , excitationType , irmodelistfilename , nquanta , irfront , irtail , asinitial , asfinal , T , ZPEscalingFactor , mass , vibfreq , cart_q0 , vibdxdr , cart_q , cart_p );
     
  /*
  void icgen( int natomselect , int ntraj , 
              int excitationType , char irmodelistfilename[] , 
              double irfront , double irtail , double softrange , 
              double T , double * mass , double * vibfreq, double * cart_q0,
              double * U ,   double * cart_q,  double * cart_p)
  */
  
  // ======> Adjusting the units ... from ( Bohr , au_time ) to ( nm , ps )
  
  for( itmp = 0 ; itmp < ncartselect * ntraj ; itmp ++ )
  {
    *( cart_q + itmp ) = *( cart_q + itmp ) * NM2BOHR ;
    
    *( cart_p + itmp ) = *( cart_p + itmp ) * ( NM2BOHR / AUT2PS ) ;
  
  }
  
  
  
  
    //=================> For Wigner Distribution Normalization Factors <============//
  
  double * wignerPositionNormArray , * wignerVelocityNormArray ;
  
  int ilevel ;
    
  if( fixedPhaseAngle == 9.00 )
  {
    wignerPositionNormArray = calloc( nmodeselect * 3 , sizeof( double ) ) ; dzeros( nmodeselect , 3 , wignerPositionNormArray ) ; 
    
    wignerVelocityNormArray = calloc( nmodeselect * 3 , sizeof( double ) ) ; dzeros( nmodeselect , 3 , wignerVelocityNormArray ) ; 
    
    for( imode = 0 ; imode < nmodeselect ; imode ++ )
    {
      for( ilevel = 0 ; ilevel < 3 ; ilevel ++ )
        {
          *( wignerPositionNormArray + 3 * imode + ilevel ) = wignerPositionNormalize( ilevel , *( vibfreq + imode ) , ZPEscalingFactor ) ; 
          *( wignerVelocityNormArray + 3 * imode + ilevel ) = wignerVelocityNormalize( ilevel , *( vibfreq + imode ) , ZPEscalingFactor ) ; 
        }
  
    }
  
  
    
    

  }
  

//=================> For Canonical Distribution Normalization Factors <============//
 
  double * canonicalPositionNormArray , * canonicalVelocityNormArray ;
  
  //int ilevel ;
  /*  
  if( fixedPhaseAngle == 9.00 )
  {
    canonicalPositionNormArray = calloc( nmodeselect * 3 , sizeof( double ) ) ; dzeros( nmodeselect , 3 , canonicalPositionNormArray ) ; 
    
    canonicalVelocityNormArray = calloc( nmodeselect * 3 , sizeof( double ) ) ; dzeros( nmodeselect , 3 , canonicalVelocityNormArray ) ; 
    
    for( imode = 0 ; imode < nmodeselect ; imode ++ )
    {
      for( ilevel = 0 ; ilevel < 3 ; ilevel ++ )
        {
          *( canonicalPositionNormArray + 3 * imode + ilevel ) = canonicalPositionNormalize( ilevel , *( vibfreq + imode ) , ZPEscalingFactor , T ) ; 
          *( canonicalVelocityNormArray + 3 * imode + ilevel ) = canonicalVelocityNormalize( ilevel , *( vibfreq + imode ) , ZPEscalingFactor , T ) ; 
        }
  
    }
  
  
    

  }
  */
  


  
  
  
  
  // =====> Output the data into dat files ... 
  
  
  char coordinateName[ MAXCHARINLINE ] , velocityName[ MAXCHARINLINE ] ;
  
  double * tmp_crd = calloc( ncartselect , sizeof( double ) );
  
  double * tmp_vel = calloc( ncartselect , sizeof( double ) );
  
  register FILE * pcrd , * pvel ;
 
  for( itraj = 0 ; itraj < ntraj ; itraj ++ )
  {
    printf("\n-o-o-o-o-o-o-o-o Working on Traj# %d -o-o-o-o-o-o-o\n" , itraj );
    
    dzeros( ncartselect , 1 , tmp_crd ) ; dzeros( ncartselect , 1 , tmp_vel );
    
    icgen( natomselect , irstatus_list , seeda , seedn , fixedPhaseAngle , nquanta , irfront , irtail , asinitial , asfinal , 
           T , ZPEscalingFactor , mass , vibfreq , cart_q0 , vibdxdr , tmp_crd , tmp_vel , 
           wignerPositionNormArray , wignerVelocityNormArray , canonicalPositionNormArray , canonicalVelocityNormArray , debuggingMode );
    
    
    
    if( outUnitType == 19 ) // gmx 
    {
      for( icart = 0 ; icart < ncartselect ; icart ++ )
      {
        *( tmp_crd + icart ) = *( tmp_crd + icart ) * NM2BOHR ;
      
        *( tmp_vel + icart ) = *( tmp_vel + icart ) * NM2BOHR / AUT2PS ;
      }
    }
    else if( outUnitType == 21 ) // G09
    {
      for( icart = 0 ; icart < ncartselect ; icart ++ )
      {
        *( tmp_crd + icart ) = *( tmp_crd + icart ) * NM2BOHR * 10.00 ;
      
        *( tmp_vel + icart ) = *( tmp_vel + icart ) * NM2BOHR / AUT2PS * 10.00 ;
      }
    }

    
    sprintf( coordinateName , "%s.%05d.crd" , outnameprefix , itraj + startingFrame );
    
    sprintf( velocityName , "%s.%05d.vel" , outnameprefix , itraj + startingFrame );
  
    printf("\n[\tnatomoutput\t]=[\t%d\t]" , natomoutput ) ;

    pcrd = fopen( coordinateName , "wb+" );
    
    pvel = fopen( velocityName , "wb+" );
    
    doutput( pcrd , natomoutput , 3 , tmp_crd );
    
    doutput( pvel , natomoutput , 3 , tmp_vel );
    
    fclose( pcrd ) ; fclose( pvel );
  
  
  } 
  
  
  
  
  
  
  
  //-o-o-o-o-o-o-o ... All right... Post Calcs ... Boltzmann Energy ... o-o-o-o-o-o-o-o-//
  
  printf("\n\nALL RIGHT THEN ... NOW IT'S THE TIME TO DEBUG ...\n\n");
  
  printf("\n\n...orz...orz...orz...orz...orz...orz...orz...orz...orz\n\n");
  
  double * boltenergy;

  boltenergy = calloc( nmodeselect , sizeof( double ) );
  
  double fb , omegacurr;
  
  int istate ;
  
  for( imode = 0 ; imode < nmodeselect ; imode ++ )
  {
    irstatus = *( irstatus_list + imode ) ;
    
    dtmp = *( vibfreq + imode );
    
    omegacurr = *( vibfreq + imode );
    
    printf("\nFrequency of %d mode is % 12.8f ... \n" , imode , omegacurr );
    
    if ( irstatus == 0 )
    {
      dtmp = boltzmannAvgEnergy( omegacurr , T );
  
      *( boltenergy + imode ) = dtmp ;
    }
    /*
    for ( istate = 0 ; istate < 400 ; istate ++ )
    {
      fb = boltzmannfactor( omegacurr, istate, T );
      
      dtmp = dtmp + fb * omegacurr * ( istate + 0.500 );
    
    }
    */
    else
    {
      *( boltenergy + imode ) = dtmp ;
    }

  }
  
  
  if(  debuggingMode == YES )
  {
    debug = fopen( "boltenergy.deb", "wb+" );
  
    doutput( debug , nmodeselect , 1 , boltenergy );

    fclose( debug );

  }





  return( 0 );






}






