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
  double atommass ;
  

} GRO ;



int main(int argc, char * argv[])
{

  //-o-o-o-o-o-o-o-o Declaring Variables o-o-o-o-o-o-o-o-o-//
  
  FILE * pmass , * phess , * pEqCrd , * poutDXDR , * poutFreq ; //, *pmo2ao;

  char massFileName[ 100 ] , hessFileName[ 100 ] , EqCrdFileName[ 100 ] ;
  
  char outDXDRFileName[ 100 ] , outFreqFileName[ 100 ] ;

  char massunit[ 100 ] , hessunit[ 100 ] , crdunit[ 100 ] ;
  
  int sdouble = sizeof(double);

  int natom , ncart , nmode ;
  
  int iatom , icart , imode ;

  int natomselect , natomprovide , ncartselect , ncartprovide ;

  int len_hess_cart, len_tri_hess_cart;

  int j,k,m;

  int itmp; double dtmp;

  char * keyword;

  int debuggingMode = NO ;

  //char gessname[ 100 ] ;

  double * hess_cart , * hess_cart_select ;
  
  double * freq , * vib_freq , * dxdr , * vib_dxdr ;
  
  double hessconvert , crdconvert , massconvert ;

  double * mass_provide , * mass; 

  char ** pcmd ; pcmd = argv ;

  int icmd = 1 ;
  
  
  int iline , iload ;
  
  int irow , icol ;
  
  int blank_signal , groinfo ;

  char buffer[ MAXCHARINLINE ] ;
  
  char cache[ MAXCHARINLINE ] ;
  
  char tmpString[ MAXCHARINLINE ] ;

  


  
  // ---> For debugging output 

  FILE * debug;



  //-o-o-o-o-o-o-o-o Recording cmd-line and time stamp o-o-o-o-o-o-o-o-o-//
  
  time_t current_time;

  time( &current_time );

  char now[ 300 ] ;

  strcpy( now , ctime( &current_time ) );

  int lennow = strlen( now ) ;

  *( now + lennow - 1 ) = ' ';

    
    
    printf("\n**********************************************************************\n");
      printf("* G_GMXFREQ_D : Stand-Alone Utility to Calculate Normal Modes        *\n");
      printf("* From Hessian File and Mass File.                                   *\n");
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

 


  //-o-o-o-o-o-o-o Read input command-line arguments and input files while deciding some parameters ...o-o-o-o-o-o-o-o-o-o-//

  // -------> Parsing the Command Line Arguments ... 
  
  pcmd = argv ; 
  
  
  
  //int exn = 10 ; int exr  = 16 ; int exR = 18 ; int exH  = 22 ; int exM = 28 ; int exL = 30 ; 

  int exc = 20 ; 
  
  int exs = 18 ; int exm = 88 ;
  
  int exH = 70 ; 
  
  int exo = 22 , exw = 26 ;

  char * flag ;
  
  
  int internalOrNot = NO ;
  
  
  icmd = 1 ;
  
  printf("\n%d command-line arguments provided ...\n" , argc );
  
  if( argc == 1 )
  {
    printf("\n\nNo command-line arguments provided ... Mission aborting ...\n\n");
    
    printf("\nPlease refer to the usage by typing ' %s -h '\n\n" , * argv );
    
    exit(1); 
  
  
  }

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
	      
	                 strcpy( buffer , *( ++ pcmd ) );
	                 
	                 if( strcmp( tmpString , "none" ) == 0 )
	                 {
	                   strcpy( EqCrdFileName , "dirtroad.anthem" ) ;
	                   
	                   internalOrNot = NO ;
	                   
	                   printf("\nCommand-line argument indicates : NO INPUT Eq. Coordinate ... Hence NO Internal Coordinate Transformation \n" ); 
	                 
	                 }
	                 else
	                 {
	                   strcpy( EqCrdFileName , tmpString ) ;
	                   
	                   strcpy( crdunit , buffer ) ;
	                   
	                   if( * crdunit == '-' )
	                   {
	                     printf("\nHey, you did not specify the unit of you Coordinate file ... !\n");
	                     
	                     exit( 11 ) ;
	                   
	                   }
	                 
	                   internalOrNot = YES ;
	                   
	                   printf("\nCommand-line argument indicates : Input Eq. Coordinate File name : %s . Format is %s ...\n" , EqCrdFileName , crdunit ); 
	                 
	                 
	                 }
	                 	                 	           
	                 exc = 21 ;	           
	                 	                 
	                 icmd = icmd + 3 ;
	                 
	                 break ;
	               
	      

	      case 'H' : strcpy( hessFileName , *( ++ pcmd ) );
	      
	                 strcpy( hessunit , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Input Hessian File name : %s . Unit is %s ...\n" , hessFileName , hessunit ); 
	                 
	                 if( * hessunit == '-' )
	                 {
	                   printf("\nHey, you did not specify the unit of you Hessian file ... !\n");
	                   
	                   exit( 11 );
	                 }
	                 	           
	                 exH = 71 ;	           
	                 	                 
	                 icmd = icmd + 3 ;
	                 
	                 break ;


	      case 'm' : strcpy( massFileName , *( ++ pcmd ) );
	      
	                 strcpy( massunit , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Input Mass File name : %s . Unit is %s ...\n" , massFileName , massunit ); 
	                 
	                 if( * crdunit == '-' )
	                 {
	                   printf("\nHey, you did not specify the unit of you Coordinate file ... !\n");
	                   
	                   exit( 11 );
	                 }
	                 	           
	                 exm = 89 ;	           
	                 	                 
	                 icmd = icmd + 3 ;
	                 
	                 break ;
	                 
               

	      case 's' : strcpy( tmpString , *( ++ pcmd ) ) ;
	      
	                 if( strcmp( tmpString , "all" ) == 0 || strcmp( tmpString , "All" ) == 0 )
	                 {
	                   exs = 20 ;
	                   
	                   printf("\nCommand-line argument indicates : All atoms will be chosen as solute ...\n" );
	                 }
	                 else
	                 {
	                   printf("\nReceived information : %s ...\n" , tmpString ) ;
	                   
	                   natomselect = atoi( tmpString ); 
	                  
	                   exs = 19 ;
	                   
	                   printf("\nCommand-line argument indicates : First %d atoms will be chosen as solute ...\n" , natomselect );
	                 }
	                 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ;
	      

	      case 'o' : strcpy( outDXDRFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Output Normal Mode file name : %s ...\n" , outDXDRFileName ); 
	                 
	                 exo = 23 ;
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ; 
	                 
	                 
	      case 'w' : strcpy( outFreqFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Output frequency file name : %s ...\n" , outFreqFileName ); 
	                 
	                 exw = 27 ;
	                 
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




		                 
	      case 'h' : printf("\nUsage:  % 15s [ -c 'input Equilibrium Geometry' ( gro / gmxcrd / g09crd )] [ -H 'Hessian file name' ' unit : au or gmx' ] [ -s # of atoms chosen as the solute ]" , * argv) ;
	                 printf("\n                        [ -m Mass file name & unit ( au or amu ) ] [ -w output frequency file name ] [ -o output dxdr file name ]\n\n");
	                 //printf("\n           [ -c 'input EqMD gro file of solvent' ][ -p P Group Name ][ -q Q Group Name ][ -N atom number of L-Shape reference atom ]");
	                 //printf("\n           [ -w whether to perform van der Waals contacting check ( YES / Yes / yes ) or ( NO / No / no ) or floating point number to indicate scaling factor ]");
	                 //printf("\n           [ -s # of atoms in solute molecule ] [ -t distance threshold to accept one alignment , unit = Angstrom ]\n\n" ); 
	                 
	                 //printf("\nNote : 1) For \"-w\" option, YES/Yes/yes will cause the vdw-check to perform with default scaling 1.00 while NO/No/no will shut down the vdw-check.\n");
	                 //printf("\n          Specifying a floating number will also cause the vdw-check to perform but the floating number will be the user-defined scaling factor (vdwFactor).\n");

	                 //printf("\n       2) For \"-t\" option, YES/Yes/yes will cause the universal-distance-check to perform with default threshold 1.20Å while NO/No/no will shut down the unidist-check.\n");
	                 //printf("\n          Specifying a floating number will also cause the unidist-check to perform but the floating number will be the user-defined distance threshold (vdwFactor).\n");

	                 printf("\nNote : 1) For \"-s\" option, [ -s all ] or [ -s All ] indicates all atoms chosen as solute;\n");
	                 printf("\n          Default for -s is all atoms when nresidue = 1  or natom in 1st residue when nresidue != 1 \n");
	                 
	                 printf("\n       2) For \"-c\" option, [ -c none none ] will turn off the internal coordinate transformation and perform regular Cartesian coordinate normal mode analysis ... \n\n\n");


	                 icmd = icmd + 1 ; 
	                 
	                 exit( 1 ) ;
	                 
	                 break ;



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

 
  // -------> Default File Names
  

  if( internalOrNot == YES && exc == 20 )
  {
    strcpy( crdunit , "gro" ) ;
    
    strcpy( EqCrdFileName , "system.gro" ) ;
    
    printf("\nNo input .gro file provided, default \"system.gro\" in play ...\n") ;
  
  }

  
  int lenEqGROFileName = strlen( EqCrdFileName ) ;
  
  if( exo == 22 )
  {
    strncpy( outDXDRFileName , EqCrdFileName , lenEqGROFileName - 4 ) ;
    
    *( outDXDRFileName + lenEqGROFileName - 4 ) = '\0' ;
    
    strcat( outDXDRFileName , ".dxdr" ) ;
    
    printf("\nBy default , output DXDR file name will be [ %s ] ...\n\n" , outDXDRFileName ) ;
  
  }
  
  
  
  if( exw == 26 )
  {
    strncpy( outFreqFileName , EqCrdFileName , lenEqGROFileName - 4 ) ;
    
    *( outFreqFileName + lenEqGROFileName - 4 ) = '\0' ;
    
    strcat( outFreqFileName , ".freq" ) ;
    
    printf("\nBy default , output DXDR file name will be [ %s ] ...\n\n" , outFreqFileName ) ;
  
  }

  

  if( exH == 70 )
  {
    strncpy( hessFileName , EqCrdFileName , lenEqGROFileName - 4 ) ;
    
    *( hessFileName + lenEqGROFileName - 4 ) = '\0' ;
    
    strcat( hessFileName , ".gess" ) ;
    
    printf("\nBy default , searching for Hessian file with the name %s ...\n\n" , hessFileName ) ;
    
    strcpy( hessunit , "gmx" ) ;
  
  }
  

  if( exm == 88 )
  {
    strncpy( massFileName , EqCrdFileName , lenEqGROFileName - 4 ) ;
    
    *( massFileName + lenEqGROFileName - 4 ) = '\0' ;
    
    strcat( massFileName , ".mass" ) ;
    
    printf("\nBy default , searching for Mass file with the name %s ...\n\n" , massFileName ) ;
    
    strcpy( massunit , "amu" ) ;
  
  }


 
  if( strcmp( hessunit , "au" ) == 0 )
  {
    printf("\nThe unit this Hessian file in is ATOMIC UNIT : Hartree/(Bohr^2) ... \n");
  
    hessconvert = 1.0000 ;
  }
  else if( strcmp( hessunit , "gmx" ) == 0 )
  {
    printf("\nThe unit this Hessian file in is GMX-Unit : kJ/mol/(nm^2) ... \n");
  
    hessconvert = HESSAU2GMX ;
  }
  else
  {
    printf("\nWhat did you say about the unit of you frequency again ??? \n");
  
    exit( 107 );
  }




  printf("\n===> Done Defining Default Fila Names <===\n\n\n") ;
  
  
  
  // -------> Confirming File Access ... 

  
  if( ( pmass = fopen( massFileName ,"r" ) ) == NULL )
  {  
    printf("\nUser defined mass-file [ %s ] does not exist...\n" , massFileName );

    exit( 63 );
  }
  
  
  
  if( ( pEqCrd = fopen( EqCrdFileName ,"r" ) ) == NULL && internalOrNot == YES )
  {  
    printf("\nUser defined Equilibrium Coordinate File [ %s ] does not exist...\n" , EqCrdFileName );

    exit( 63 );
  }
  


  if( ( phess = fopen( hessFileName ,"r" ) ) == NULL )
  {  
    printf("\nUser defined Hessian-file [ %s ] does not exist...\n" , hessFileName );

    exit( 63 );
  }
  


  printf("\n===> Done Confirming File Access <===\n\n\n") ;

  
  
  //-o-o-o-o-o-o-o Loading Eq. Coordinate , Hessian & Mass ...o-o-o-o-o-o-o-o-o-o-//
  
  //-----> Pre-Loading Equilibrium Geometry ...
  
  char grotitlestring[MAXLINE];
  
  iline = 3 ;
  
  iload = 0 ;
  
  int natomgroline , natomgrotitle ;
  
  int exVelocity = NO ;
  

  if( internalOrNot == YES && ( strcmp( crdunit , "gro" ) ) == 0 )
  {
    rewind( pEqCrd );
    
    //printf("\nCurrent character is %c ... \n" , fgetc( pEqGRO ) );
    
    fskip( pEqCrd , 1 );

    fscanf( pEqCrd , "%d" , &natomgrotitle );
    
    fskip( pEqCrd , 1 );
    
    printf("\n Second line of .gro file says it is describing %d atoms ... \n\n" , natomgrotitle );
    
    
    while( ( groinfo = freadline( buffer , MAXCHARINLINE , pEqCrd , ';' ) ) != 0 )
    {
      itmp = inLineWC( buffer ) ;
      
      break ;
    }
        
    if( itmp == 6 )
    {
      exVelocity = NO ;
      
      printf("\nI see there is no velocity information in .gro file ...\n\n") ;
      
    }
    else if( itmp == 9 )
    {
      exVelocity = YES ;
      
      printf("\nI see velocity information is also included in .gro file ...\n\n") ;
    
    }
    else
    {
      printf("\nPlease check the format of you .gro file ... There are %d words in one line of your molecular specification ... \n" , itmp ) ;
      
      exit( 456 );
    }
    
  
    printf("\nNow let's pre-load the .gro file ans see how many atoms it is describing ... \n");
    
    
    
    natomgroline = preLoadGRO( pEqCrd ) ;
  
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
  
    ncart = 3 * natom ;

    
    itmp = 0 ; 
   
    if( exs == 18 )
    {
      natomselect = natom ;
    }
    else if( exs == 19 && natomselect > natom ) // Will be dead ... 
    {
      printf("\nThere are only %d atoms in this system ... you cannot select more than that ... \n" , natom );
      
      if( natomgrotitle > natomgroline && natomselect <= natomgrotitle )
      {
        printf("\nAlthough ... the second line of your initial .gro file did indicate there were supposed to be %d atoms in system ... So go back and make sure what you are trying to do ... \n" , natomgrotitle );
      }
      else if( natomgrotitle < natomgroline && natomselect <= natomgroline )
      {
        printf("\nAlthough ... your initial .gro file did describe %d atoms in system ... So go back and make sure what you are trying to do ... \n" , natomgroline );
      }
      
      exit( 78 );
    }
    else if( exs == 19 && natomselect <= natom )
    {
      printf("\nYou have selected %d atoms as solute ... There are %d atoms en toto in this system ... \n" , natomselect , natom );    

    }
    else if( exs == 20 )
    {
      natomselect = natom ;
    }
    else
    {
      printf("\nSomething is wrong with the atom selection process ... NAtom = %d , NAtomSelect = %d ... \n" , natom , natomselect );
      
      exit( 78 );
    }
  
    ncartselect = 3 * natomselect ;
  
  
  
  }
  else if( internalOrNot == YES && ( ( strcmp( crdunit , "gmxcrd" ) ) == 0 || ( strcmp( crdunit , "g09crd" ) ) == 0 ) )
  {
    rewind( pEqCrd );
    
    itmp = flength( pEqCrd ) ;
    
    natom = itmp / 3 ;
    
    ncart = itmp ;
  
    if( exs == 18 )
    {
      natomselect = natom ;
    }
    else if( exs == 19 && natomselect > natom ) // Will be dead ... 
    {
      printf("\nThere are only %d atoms in this system ... you cannot select more than that ... \n" , natom );
            
      exit( 78 );
    }
    else if( exs == 19 && natomselect <= natom )
    {
      printf("\nYou have selected %d atoms as solute ... There are %d atoms en toto in this system ... \n" , natomselect , natom );    

    }
    else if( exs == 20 )
    {
      natomselect = natom ;
    }
    else
    {
      printf("\nSomething is wrong with the atom selection process ... NAtom = %d , NAtomSelect = %d ... \n" , natom , natomselect );
      
      exit( 78 );
    }
  
    ncartselect = 3 * natomselect ;
    
  }
  
  
  //-----> Pre-Loading Mass ...

  itmp = flength( pmass ) ;

   
  rewind( pmass ) ;

  
  if( internalOrNot == NO )
  {
    natom = itmp ;
    
    ncart = 3 * natom ;
    
    if( exs == 18 )
    {
      natomselect = natom ;
      
      printf("\nBy default, all available atoms will be chosen as solute ...\n\n") ;
    }
    else if( exs == 20 )
    {
      natomselect = natom ;
      
      printf("\nPer user's request, all available atoms will be chosen as solute ...\n\n") ;

    }
    else if( exs == 19 )
    {
      if( natomselect == natom )
      {
        printf("\nOKay ... I see you provided the mass for all the atom in this system ... \n"); 
      }
      else if( natomselect > natom ) // will be dead
      {
        printf("\nERROR : There are only %d atoms according to mass file, I cannot select that many atoms as solute ...\n\n" , natomselect ) ;
        
        exit( 89 ) ;
      }
      else if( natomselect < natom )
      {
        printf("\nOKay ... you provided %d floating-point numbers for mass so the total NAtom in system is %d while %d atoms are selected ...\n" , itmp , natom , natomselect );

        printf("\nSo ... I will assume the first %d numbers in your mass file correspond to the selected atoms ...\n" , natomselect );
      }
      else
      {
        printf("\nSomething is wrong with the mass file ... There are %d atomic mass in the file while the total NAtom in this system is %d and %d atoms are selected ... \n" , itmp , natom , natomselect );

        exit( 81 );

      }
    
    }
    
    ncartselect = natomselect * 3 ;
  
  
  }
  else if( internalOrNot == YES )
  {
    if( itmp == natomselect )
    {
      printf("\nOKay ... I see you provided the mass info for just the selected atoms ... \n");
    }
    else if( itmp == natom )
    {
      printf("\nOKay ... I see you provided the mass for all the atom in this system ... \n"); 
    }
    else if( itmp > natom )
    {
      printf("\nOKay ... you provided %d numbers for mass. The total NAtom in system is %d while %d atoms are selected ...\n" , itmp , natom , natomselect );

      printf("\nSo ... I will assume the first %d numbers in your mass file correspond to the selected atoms ...\n" , natomselect );
    }
    else
    {
      printf("\nSomething is wrong with the mass file ... There are %d atomic mass in the file while the total NAtom in this system is %d and %d atoms are selected ... \n" , itmp , natom , natomselect );

      exit( 81 );

    }
  
  }
  
  // ---> Loading Equilibrium Geometry ...
  
  GRO EqAtomList[ natomselect ] ; 
  
  double * EqCrd = calloc( ncartselect , sizeof( double ) ) ;
  
  dzeros( ncartselect , 1 , EqCrd ) ;
  
  double boxvector[ 3 ]; dzeros( 3 , 1 , boxvector ) ;
  
  char tmp_char ;
  
  if( internalOrNot == YES && ( strcmp( crdunit , "gro" ) ) == 0 ) // We want every coordinate to be in BOHR unit
  {
    crdconvert = NM2BOHR ;
    
    rewind( pEqCrd );
      
    fskip( pEqCrd , 1 );

    fskip( pEqCrd , 1 );
      
    printf("\nNow let's read the actual .gro file  ... \n");
    
    iload = 0 ; iline = 0 ;
    
    while( ( groinfo = freadline( buffer , MAXCHARINLINE , pEqCrd , ';' ) ) != 0 )
    { 
      //printf("\n//--------------> WORKING ON NO. %d LINE ... <-------------//\n" , iline );
      
      blank_signal = stellblank( buffer ) ;
      
      if( blank_signal == 0 )
      {
        //printf("\nNo.%d line is a blank line ... Moving on ...\n" , iline ) ;
        
        continue ;
      }
      else if( blank_signal == 1 )
      {  
        //printf("\nNo.%d line is NOT a blank line ... loading ...\n" , iline );
        
        if( ( tmp_char = getfirst( buffer ) ) == ';' )
        {
          //printf("\nThis is a comment line ... So nothing will be loaded ...\n");
          
          //fskip( pinputfile , 1 );
          
          continue ;
        }
        else
        {
          //printf("\nLine reads : %s ...\n" , buffer );
                
          sscanf( buffer , "%5d%5s" , &EqAtomList[ iload ].resnumber , EqAtomList[ iload ].resname );

          //printf( "%s\t" , EqAtomList[ iatom ].resname );

          //sscanf( pEqGRO , "%s" , EqAtomList[ iload ].atomname );
          
          strpickword( buffer , 2 , cache ) ;
          
          strcpy( EqAtomList[ iload ].atomname , cache ) ;

          //printf( "%s" , EqAtomList[ iatom ].atomname );

          //sscanf( pEqGRO , "%d" , &EqAtomList[ iload ].atomnumber );
          
          strpickword( buffer , 3 , cache ) ;
          
          EqAtomList[ iload ].atomnumber = atoi( cache ) ;

          //printf( "\nWorking on No. %d atom ...\n" , EqAtomList[ iatom ].atomnumber );

          //sscanf( pEqGRO , "%lf" , &EqAtomList[ iload ].cx ); //printf("\n Cx is %lf ...\t" , EqAtomList[ iatom ].cx);
          
          strpickword( buffer , 4 , cache ) ; EqAtomList[ iload ].cx = atof( cache ) ;
          
          *( EqCrd + 3 * iload + 0 ) = EqAtomList[ iload ].cx / crdconvert ;
          
          //sscanf( pEqGRO , "%lf" , &EqAtomList[ iload ].cy ); //printf("\n Cy is %lf ...\t" , EqAtomList[ iatom ].cy);
          
          strpickword( buffer , 5 , cache ) ; EqAtomList[ iload ].cy = atof( cache ) ;
          
          *( EqCrd + 3 * iload + 1 ) = EqAtomList[ iload ].cy / crdconvert ;
          
          //sscanf( pEqGRO , "%lf" , &EqAtomList[ iload ].cz ); //printf("\n Cz is %lf ...\n\n" , EqAtomList[ iatom ].cz);
          
          strpickword( buffer , 6 , cache ) ; EqAtomList[ iload ].cz = atof( cache ) ;
          
          *( EqCrd + 3 * iload + 2 ) = EqAtomList[ iload ].cz / crdconvert ;
          
          
          if( exVelocity == YES )
          {
            strpickword( buffer , 7 , cache ) ; EqAtomList[ iload ].vx = atof( cache ) ;
            
            strpickword( buffer , 8 , cache ) ; EqAtomList[ iload ].vy = atof( cache ) ;
            
            strpickword( buffer , 9 , cache ) ; EqAtomList[ iload ].vz = atof( cache ) ;
          
          }
          
          //fscanf( pgroinput , "%lf" , &EqAtomList[ iatom ].vx ); //printf("\n Vx is %lf ...\t" , EqAtomList[ iatom ].cx);

          //fscanf( pgroinput , "%lf" , &EqAtomList[ iatom ].vy ); //printf("\n Vy is %lf ...\t" , EqAtomList[ iatom ].cy);
    
          //fscanf( pgroinput , "%lf" , &EqAtomList[ iatom ].vz ); //printf("\n Vz is %lf ...\n\n" , EqAtomList[ iatom ].cz);

          
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
      
      if( iload == natomselect ) break ;

    }
    
    
    /*
    if( iline < natomgroline )
    {
      fskip( pEqCrd , natomgroline - iline ) ;
    }
    
     
    
    fscanf( pgroinput , "%lf" , boxvector + 0 );

    fscanf( pgroinput , "%lf" , boxvector + 1 );

    fscanf( pgroinput , "%lf" , boxvector + 2 );
    */
    
    
  
  }
  else if( internalOrNot == YES && ( strcmp( crdunit , "grocrd" ) ) == 0 )
  {
    crdconvert = NM2BOHR ;
    
    rewind( pEqCrd ) ;
    
    for( icart = 0 ; icart < ncartselect ; icart ++ )
    {
      fscanf( pEqCrd , "%lf" , &dtmp ) ;
    
      *( EqCrd + icart ) = dtmp / crdconvert ;
    
    }
  
  }
  else if( internalOrNot == YES && ( strcmp( crdunit , "g09crd" ) ) == 0 )
  {
    crdconvert = 1.0000 ;
    
    rewind( pEqCrd ) ;
    
    for( icart = 0 ; icart < ncartselect ; icart ++ )
    {
      fscanf( pEqCrd , "%lf" , &dtmp ) ;
    
      *( EqCrd + icart ) = dtmp / crdconvert ;
    
    }
  
  }
  else if( internalOrNot == NO )
  {
    printf("\nAlready told you ... NO INTERNAL COORDINATE BUSINESS ...\n\n") ;
  }
  else
  {
    printf("\nUNKOWN ERROR : INVALID COORDINATE FILE FORMAT ...\n\n");
    
    exit( 73 ) ;
  }
  
  printf("\n===> Finished Loading Equilibrium Geometry <===\n\n\n") ;
  
  //dtranspose_nonsquare( natomselect , 3 , EqCrd , EqCrd ) ;
  
  
  /*
  debug = fopen("transposed_eqgeom.deb" , "wb+") ;
  
  doutput( debug , 3 , natomselect , EqCrd ) ;
  
  fclose( debug ) ;
  
  
  printf("\n===> Done Transposing Equilibrium Geometry <===\n\n\n") ;
  */
  
  // -------> Loading Mass ...  
  
  mass = calloc( natomselect , sizeof( double ) );
  
  dzeros( natomselect , 1 , mass );
  
  if( strcmp( massunit , "au" ) == 0 )
  {
    massconvert = AMU2AU ;
  }
  else if( strcmp( massunit , "amu" ) == 0 )
  {
    massconvert = 1.000 ;
  }
  else
  {
    printf("\nUNKOWN ERROR : INVALID MASS FILE FORMAT ...\n\n");
    
    exit( 75 ) ;
  }
  
  
  for( iatom = 0 ; iatom < natomselect ; iatom ++ )
  {
    fscanf( pmass , "%lf" , &dtmp ) ;
    
    *( mass + iatom ) = dtmp / massconvert ;
  }
   


  
  printf("\n===> Done Loading Mass ( in %s ) <===\n\n\n" , massunit ) ;

  
  //------> Hessian File :

  fskip( phess , 2 );

  len_hess_cart = flength( phess );

  if( len_hess_cart == ncart * ncart )
  {
    printf("\nOKay ... I see you provided the Hessian for all %d atom in system ... \n" , natom );
  }
  else if( len_hess_cart == ncartselect * ncartselect )
  {
    printf("\nOKay ... I see you only provided the Hessian for the %d selected atoms ... \n" , natomselect );
  }
  else if( len_hess_cart < ncart * ncart && len_hess_cart > ncartselect * ncartselect )
  {
    printf("\nThere are %d numbers in the Hessian file which is less than the square of total %d Cartesian coordinates in system but more than that of %d Cartesian coordinates of selected atoms ... \n" , len_hess_cart , ncart , ncartselect );

    printf("\nSo I am taking the first %d number ( %d * %d ) as the Hessian for selected atoms ...\n" , ncartselect * ncartselect , ncartselect , ncartselect );

  }
  else
  {
    printf("\nSomething is wrong with the Hessian file ... There are %d numbers in the Hessian file ... \n" , len_hess_cart );

    exit( 77 );

  }

  rewind( phess );

  //-------> Allocating space for Hessian and Diagonalization Process ...

  ncartprovide = sqrt( len_hess_cart ); 
  
  printf("\nWell ... the Hessian file you provided is a %d * %d square matrix ...\n" , ncartprovide , ncartprovide );

  hess_cart = calloc( len_hess_cart , sizeof(double));  dzeros( len_hess_cart , 1, hess_cart );

  hess_cart_select = calloc( ncartselect * ncartselect , sizeof(double) ) ; 
  
  dzeros( ncartselect , ncartselect , hess_cart_select );
  
  double * DMatrix = calloc( ncartselect * ncartselect , sizeof( double ) ) ;
  
  dzeros( ncartselect , ncartselect , DMatrix ) ;
  
  double * tmp_hessian = calloc( ncartselect * ncartselect , sizeof(double) ) ; 
  
  dzeros( ncartselect , ncartselect , tmp_hessian );


  //tri_hess_cart = calloc( ncart*(ncart+1)/2 , sizeof(double)); 
  //dzeros(ncart*(ncart+1)/2, 1, tri_hess_cart);



  //-------> Loading Hessian Matrix ...

  rewind( phess );

  fskip( phess , 2 );

  fload( phess , hess_cart );

  for( j = 0 ; j < len_hess_cart ; j ++ )
  {
    *( hess_cart + j ) = *( hess_cart + j ) / hessconvert ;
  }


  printf("\n===> Done Loading Whole Provided Hessian <===\n\n\n") ;


  //-------> Selecting the desired part of Hessian matrix 

  for( j = 0 ; j < ncartselect ; j ++ )
  {
    for( k = 0 ; k < ncartselect ; k ++ )
    {
      *( hess_cart_select + j * ncartselect + k ) = *( hess_cart + j * ncartprovide + k );
    }

  }

  printf("\n===> Done Picking The Selected Hessian <===\n\n\n") ;


  if( debuggingMode == YES )
  {
    debug = fopen("hess_cart.deb", "wb+");

    doutput( debug , ncart , ncart , hess_cart);

    fclose(debug);


    debug = fopen("hess_cart_select.deb", "wb+");

    doutput( debug , ncartselect , ncartselect , hess_cart_select );

    fclose(debug);
  }



  //-------> Transpose Hessian_cart_select to pass into gausvib_ ...
  
  
  
  dtranspose( ncartselect , hess_cart_select , hess_cart_select ) ;


  printf("\n===> Done Transposing Selected Hessian <===\n\n\n") ;


  //void gausvib_( int * , double * , double * , double * , double * , double * ) ;
  //            ( natom, mass, coordxyzinp, hessianinp, Dmatrix )




  //-------> Generate D-Matrix and transform into internal coord ...

  if( internalOrNot == YES )
  {
    gausvib_( &natomselect , mass , EqCrd , hess_cart_select , DMatrix ) ;
  
    dtranspose( ncartselect , DMatrix , DMatrix ) ;  
  
    printf("\n===> Done Generating D-Matrix ... Currently in C-Fashion <===\n\n\n") ;
  
  }
  else if( internalOrNot == NO )
  {
    for( icart = 0 ; icart < ncartselect ; icart ++ )
    {
      *( DMatrix + icart * ncartselect + icart ) = 1.000 ;
    }
  
  }
  else
  {
    printf("\nUNKNOWN ERROR : [ internalOrNot ] = %d \n\n" , internalOrNot ) ;
    
    exit( 15 ) ;
  
  }
  
  
  if( debuggingMode == YES )
  {
    debug = fopen( "DMatrix.deb" , "wb+") ;
  
    doutput( debug , ncartselect , ncartselect , DMatrix ) ;
  
    fclose( debug ) ;
  }

  


  for( j = 0 ; j < natomselect ; j ++ )   *( mass + j ) = ( *( mass + j ) ) * AMU2AU;

  
  printf("\n===> Done Putting Mass In AU <===\n\n\n") ;


  if( debuggingMode == YES )
  {
    debug = fopen("mass_au.deb", "wb+");

    doutput( debug , natomselect , 1 , mass );

    fclose(debug);
  }

  //-------> Performing mass-weighting for the Force constant matrix  

  dtranspose( ncartselect , hess_cart_select , hess_cart_select ) ; 
  
  //transpose back to C-Fashion
  

  masswt( natomselect , mass , hess_cart_select , hess_cart_select );

  
  
  printf("\n===> Done Mass-Weighting Selected Hessian <===\n\n\n") ;



  if( debuggingMode == YES )
  {
    debug = fopen("hess_masswt_select.deb", "wb+");

    doutput( debug , ncartselect , ncartselect , hess_cart_select );

    fclose(debug);
  }


  //-------> Calculating f_int = (D.') * f_mwc * D ;
  
  double done = 1.0000 ;
  
  double dzero = 0.0000 ;
  
  int nmode_trans_rot , nvibmodes ;
  
  if( natom == 2 )
  {
    nmode_trans_rot = 5 ;
  }
  else if( natom >= 3 )
  {
    nmode_trans_rot = 6 ;
  }
  
  nvibmodes = ncartselect - nmode_trans_rot ;
  
  double * hess_int_mwc = calloc( ncartselect * ncartselect , sizeof( double ) ) ;
  
  dzeros( ncartselect , ncartselect , hess_int_mwc ) ;
  
  double * vib_hess_int_mwc = calloc( nvibmodes * nvibmodes , sizeof( double ) ) ;
  
  dzeros( nvibmodes , nvibmodes , vib_hess_int_mwc ) ;
  
  
  
  dgemm_( "N" , "T" , &ncartselect , &ncartselect , &ncartselect , &done , DMatrix , &ncartselect , hess_cart_select , &ncartselect , &dzero , tmp_hessian , &ncartselect ) ;

  dgemm_( "N" , "T" , &ncartselect , &ncartselect , &ncartselect , &done , tmp_hessian , &ncartselect , DMatrix , &ncartselect , &dzero , hess_int_mwc , &ncartselect ) ;
  
  //dtransopose

  //    SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )



  printf("\n===> Done Calculating f_int = (D.') * f_mwc * D <===\n\n\n") ;
  
  
  if( debuggingMode == YES )
  {
    debug = fopen( "hess_int_mwc.deb" , "wb+") ;
  
    doutput( debug , ncartselect , ncartselect , hess_int_mwc ) ;
  
    fclose( debug ) ;
  }
  
  
  
  
  //-------> Performing matrix-diagonalization for Hess_Cart and l=D*L ( for internalOrNot == YES )
  
  
  freq = calloc( ncartselect , sizeof( double ) ) ; dzeros( ncartselect , 1 , freq );

  vib_freq = calloc( nvibmodes , sizeof(double));  dzeros( nvibmodes , 1 , vib_freq );

  dxdr = calloc( ncartselect * ncartselect , sizeof(double));  
  
  dzeros( ncartselect * ncartselect , 1 , dxdr );
  
  vib_dxdr = calloc( ncartselect * nvibmodes , sizeof( double ) );  
  
  dzeros( ncartselect * nvibmodes , 1 , vib_dxdr ) ;
  
  double * tmp_dxdr = calloc( nvibmodes * nvibmodes , sizeof( double ) ) ;
  
  dzeros( nvibmodes , nvibmodes , tmp_dxdr ) ;
  
  double * tmp_dxdr_2 = calloc( nvibmodes * nvibmodes , sizeof( double ) ) ;
  
  dzeros( nvibmodes , nvibmodes , tmp_dxdr_2 ) ;
  

  if( internalOrNot == YES )
  {
    
    // ---> Taking the vibrational ( 3N - 6 )-by-( 3N - 6 ) block 
    
    for( icart = nmode_trans_rot ; icart < ncartselect ; icart ++ )
    {
      for( itmp = nmode_trans_rot ; itmp < ncartselect ; itmp ++ ) 
      {
        *( vib_hess_int_mwc + ( icart - nmode_trans_rot ) * nvibmodes + ( itmp - nmode_trans_rot ) ) = *( hess_int_mwc + icart * ncartselect + itmp ) ;
    
      }
  
    }

    if( debuggingMode == YES )
    {
      debug = fopen( "vib_hess_int_mwc.deb" , "wb+" ) ;
  
      doutput( debug , nvibmodes , nvibmodes , vib_hess_int_mwc ) ;
  
      fclose( debug ) ;
    }
  
    // ---> Diagonalization
  
    dsyev_f2c( nvibmodes , vib_hess_int_mwc , tmp_dxdr, vib_freq  );

    dtranspose( nvibmodes , tmp_dxdr , tmp_dxdr );
  
    
    
    if( debuggingMode == YES )
    {
      
      /*
      debug = fopen( "vib_dxdr_internal.deb" , "wb+" ) ;
  
      doutput( debug , nvibmodes , nvibmodes , tmp_dxdr ) ;
  
      fclose( debug ) ;
      */
      
      debug = fopen( "vib_freq_internal.deb" , "wb+" ) ;
      
      doutput( debug , nvibmodes , 1 , vib_freq ) ;
      
      fclose( debug ) ;
      
    }



  
    printf("\n===> Done Diagonalizing ( 3N - 6 )-by-( 3N - 6 ) internal coordinate mass-weighted Hessian and transpose 3N-6 normal modes back to C <===\n\n\n") ;
    
  
    // ---> Putting ( 3N - 6 )-by-( 3N - 6 ) dxdr into vib_dxdr which is 3N-by-( 3N - 6 )

    for( icart = nmode_trans_rot ; icart < ncartselect ; icart ++ )
    {
      for( imode = 0 ; imode < nvibmodes ; imode ++ ) 
      {
        *( vib_dxdr + icart * nvibmodes + imode ) = *( tmp_dxdr + ( icart - nmode_trans_rot ) * nvibmodes + imode ) ;
    
      }
  
    }
  
  
    printf("\n===> Done Putting ( 3N - 6 )-by-( 3N - 6 ) dxdr into vib_dxdr which is 3N-by-( 3N - 6 ) <===\n\n\n") ;
  
  
    if( debuggingMode == YES )
    {
      debug = fopen( "vib_dxdr_internal.deb" , "wb+" ) ;
  
      doutput( debug , ncartselect , nvibmodes , vib_dxdr ) ;
      
      printf("\n[===> Debug <===] [\t%lf\t%lf\t%lf\t]" , *( vib_dxdr + 10 ) , *( vib_dxdr + 19 ) , *( vib_dxdr + 38 ) ) ;
  
      fclose( debug ) ;
    }
  
    //---> Performing l = D*L
  
  
    dtranspose_nonsquare( ncartselect , nvibmodes , vib_dxdr , vib_dxdr ) ;

    dgemm_( "T" , "N" , &ncartselect , &nvibmodes , &ncartselect , &done , DMatrix , &ncartselect , vib_dxdr , &ncartselect , &dzero , tmp_dxdr_2 , &ncartselect ) ;
  
    dtranspose_nonsquare( nvibmodes , ncartselect , tmp_dxdr_2 , tmp_dxdr_2 ) ;


    if( debuggingMode == YES )
    {
      debug = fopen( "vib_dxdr_Cartesian.deb" , "wb+" ) ;
  
      doutput( debug , ncartselect , nvibmodes , tmp_dxdr_2 ) ;
  
      fclose( debug ) ;
    }
  
   // ---> Putting l into the 3N-by-3N dxdr matrix ...


  
    for( icart = 0 ; icart < ncartselect ; icart ++ )
    {
      for( imode = nmode_trans_rot ; imode < ncartselect ; imode ++ ) 
      {
        *( dxdr + icart * ncartselect + imode ) = *( tmp_dxdr_2 + icart * nvibmodes + ( imode - nmode_trans_rot ) ) ;
    
      }
  
    }
  
  
    printf("\n===> Done Calculating l = D * L <===\n\n\n") ;
    
    
    // ---> Assemble variable "freq" from "vib_freq"
    
    for( imode = nmode_trans_rot ; imode < ncartselect ; imode ++ ) 
    {
      *( freq + imode ) = *( vib_freq + imode - nmode_trans_rot ) ;
    
    }
  
  
  
  } 
  else if( internalOrNot == NO )
  {
    dsyev_f2c( ncartselect , hess_int_mwc , dxdr, freq  );

    dtranspose( ncartselect , dxdr , dxdr );
  
  
  
  } ////////// ________________________  //////////
  
  
  if( debuggingMode == YES )
  {
    /*
    debug = fopen( "dxdr.deb" , "wb+" ) ;
  
    doutput( debug , ncartselect , ncartselect , dxdr ) ;
  
    fclose( debug ) ;
  
    */


    debug = fopen("allfrequency.deb", "wb+");

    doutput( debug , ncartselect , 1 , freq );

    fclose( debug );
  
  }
  
  
  
  
  
  
  //-------> Arranging frequencies and w=sqrt(lambda) ... 

  if( natomselect == 1 )
  {
    printf("\nSeriously? Only one atom ? NO WAYYYYYY ... \n");

    exit( 90 );
  }
  else if( natomselect == 2  )
  {
    *( freq + 5 ) = sqrt( *( freq + 5 ) ) * 219474.6313705 ;
  }

  else
  {
    for( j = 0 ; j < ncartselect ; j ++ )
    {
      if( *( freq + j ) >= 0.000 )
	    
	      *( freq + j ) = sqrt( *( freq + j ) ) * 219474.6313705 ;
      else
	      //*( freq + j ) = -1.000 * sqrt( -1.0000 * ( *( freq + j ) ) ) * 219474.6313705 ;
	      *( freq + j ) = 0.0000 ;
  
    }

  }

  
  // Here, no matter it is internal or not, after unit change, we will need to assign freq into vib_freq again ...
  
  for( imode = 0 ; imode < nvibmodes ; imode ++ )
  {
    *( vib_freq + imode ) = *( freq + imode + nmode_trans_rot ) ;
  }

  printf("\n===> Done Putting Vib-Frequencies Together <===\n\n\n") ;

  


  poutDXDR = fopen( outDXDRFileName , "wb+" ) ;
  
  doutput( poutDXDR , ncartselect , ncartselect , dxdr ) ;

  fclose( poutDXDR ) ;


  poutFreq = fopen( outFreqFileName, "wb+");

  doutput( poutFreq , ncartselect , 1 , freq );

  fclose( poutFreq );




  printf("\n\nI assume it's Allllllll Done   ...\n\n");




/*

  FILE * pmass , * phess , * pEqCrd , * poutDXDR , * poutFreq ; //, *pmo2ao;

  char massFileName[ 100 ] , hessFileName[ 100 ] , EqCrdFileName[ 100 ] ;
  
  char outDXDRFileName[ 100 ] , outFreqFileName[ 100 ] ;




*/


















  return( 0 ) ;

}



