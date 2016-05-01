#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "mutil.h"
#include "futil.h"
#include "gutil.h"

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

#ifndef EV2HARTREE
#define EV2HARTREE 1.000/27.211
#endif

#ifndef CM2HARTREE
#define CM2HARTREE 1.000/219474.63
#endif

#ifndef DEBYE2AU
#define DEBYE2AU 0.393430307
#endif


#ifndef YES
#define YES 1
#endif


#ifndef NO
#define NO 0
#endif






typedef struct dipoleMoments
{   
  int initialOrbital ;
  int finalOrbital ;
  double dipX ;
  double dipY ;
  double dipZ ;
  double dipAbs ;
  
} DIP ;



int main( int argc , char * argv[ ] )
{

  char ** pcmd = argv ; 
  
  int icmd ;
  
  int moType ; // 9 = cndo ; 11 = deb ;
  
  char cndoLogFileName[ 100 ] ,  outputFileName[ 100 ] , energyFileName[ 100 ] , fchkFileName[ 100 ];
  
  FILE * pcndoLogFile , * penergyFile , * poutputFile , * pfchkFile ;
  
  double * cndoMO  , * MOEnergy ;
  
  int HOMO , LUMO ;
  
  int nbasis , natom , nelectron , ncistate ;
  
  int ibasis , iatom , icistate ;
  
  
  
  
  
  
  
  int itmp , itmp2 ; 
  
  double tmp ; char ctmp ; 
  
  char stmp[ 150 ] , stmp2[ 150 ] , tmpString[ 150 ];
  
  int info , signal , blank_signal ;
  
  char buffer[ MAXCHARINLINE ] , cache[ MAXCHARINLINE ] ;
  
  FILE * debug ; 
  



 
  time_t current_time;

  time( &current_time );

  char now[ 300 ] ;

  strcpy( now , ctime( &current_time ) );

  int lennow = strlen( now ) ;

  *( now + lennow - 1 ) = ' ';

    
    
    printf("\n**********************************************************************\n");
      printf("* G_FAKEG09MO_D : Format CNDO or Generic MO Coefficient in G09 form  *\n");
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
      printf("*                                                                    *\n");
      printf("**********************************************************************\n");

 



  // =====> Defaults ...
  
  strcpy( cndoLogFileName , "system.snap.log" ) ;
    
  //strcpy( outputFileName , "fakeG09mo.deb" ) ;
  
  
  

  
  
  // =====> Parsing cmd-line arguments ...
  
  int exE = NO ;
  
  int exFCHK = NO ;
  
  int exo = NO ;
    
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
	      
	      case 't' : strcpy( tmpString , *( ++ pcmd ) ) ;
	      
	                 if( strcmp( tmpString , "cndo" ) == 0 || strcmp( tmpString , "CNDO" ) == 0 )
	                 {
	                   moType = 9 ;
	                   
	                   printf("\nInput file type is : CNDO Log File ...\n");
	                 
	                 }
	                 else if( strcmp( tmpString , "deb" ) == 0 || strcmp( tmpString , "gen" ) == 0 )
	                 {
	                   moType = 11 ;
	                   
	                   printf("\nInput file type is : Generic ASCII File ...\n");
	                 
	                 }
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;
	            
	      case 'l' : strcpy( cndoLogFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Input file name : %s ...\n" , cndoLogFileName ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ; 
	      
	      
	      case 'E' : strcpy( energyFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Input energy file name : %s ...\n" , energyFileName ); 
	                 
	                 exE = YES ;
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ; 

	      
	      case 'f' : strcpy( fchkFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Input file name : %s ...\n" , fchkFileName ); 
	                 
	                 exFCHK = YES ;
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ; 
	      

	      case 'o' : strcpy( outputFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Output File name : %s ...\n" , outputFileName ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 exo = YES ;
	                 
	                 break ; 

	      case 'h' : printf("\nUsage:  %s [ -l 'Input file name (cndo log file OR Generic ASCII txt file)' ] [ -t input file type ( CNDO/cndo or gen )] [ -E energy file ] [ -f Input Gaussian .fchk File Name ] [ -o ' Output FAKE .fchk file name ' ]  \n\n" , * argv ); 
	                 
	                 //printf("\n\n==> NOTE : 1) -t \n\n");
	                 
	                 //printf("\nUsage:  %s [ -t G09 calculation type : 1=ONIOM ; 2=Point Charge ] [ -f 'input gro file name' ] [(optional) -o 'output g09 file name' ] [ -n # of layers (integer) ] [ (optional) -r radius of middle layer (real) ] [-R radius of lower layer (real) ] [ -H method for Highest layer (string) ] [ (optional) -M method for Middle layer (string) ] [ -L method for Lower layer (string) ] [ -x input GMX .itp file ]\n\n" , * argv ); 
	      
	                 //exh = 9 ;
	                 
	                 icmd = icmd + 1 ; 
	                 
	                 exit( 1 ) ;
	                 
	      

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
  
  
  // =====> Verify File Access <===== //
  
  
  if( ( pcndoLogFile = fopen( cndoLogFileName , "r" ) ) == NULL )
  {
    printf("\nUser specified CNDO Log File NOT FOUND ...\n");
    
    exit( 63 ) ;
  
  }
  
  if( moType == 11 && exE == NO )
  {
    printf("\nWhen generic MO file is provided , corresponding energy file is mandatory!\n\n") ;
    
    exit( 238 );
    
  }
  
  if( moType == 11 && ( penergyFile = fopen( energyFileName , "r" ) ) == NULL )
  {
    printf("\nUser specified Energy File NOT FOUND ...\n");
    
    exit( 63 ) ;
  
  }
  
  if( exFCHK == YES && ( pfchkFile = fopen( fchkFileName , "r" ) ) == NULL )
  {
    printf("\nUser specified Gaussian Formatted-Checkpoint File NOT FOUND! \n\n") ;
    
    exit( 63 ) ;
  
  }
  
  
  if( exo == NO )
  {
    if( exFCHK == YES )
    {
      strcpy( outputFileName , "fakeG09.fchk" ) ;
    }
    else if( exFCHK == NO )
    {
      strcpy( outputFileName , "fakeG09mo.deb" ) ;
    }
  
  }
  
  
  
  
  // =====> Loading Reference MO Coefficients ... <===== //
  
  int lmo ;
  
  if( moType == 9 )
  {
    // ---> NBasis
 
    rewind( pcndoLogFile ) ;

    fsearch( pcndoLogFile , "Basis" ) ;
  
    fscanf( pcndoLogFile , "%s" , stmp );
  
    if( strcmp( stmp , "fns=" ) == 0 )
    {
      fscanf( pcndoLogFile , "%d" , &nbasis );
    }
    else
    {
      fsearch( pcndoLogFile , "SHELL" );
    
      fsearch( pcndoLogFile , "TOTAL" ) ;
    
      fscanf( pcndoLogFile , "%d" , &nbasis ) ;
  
    }
  
    printf("\n%d Basis functions are being used ...\n" , nbasis ) ;
    
    lmo = nbasis * nbasis ;

  }
  else if( moType == 11 )
  {
    rewind( pcndoLogFile ) ;
    
    lmo = flength( pcndoLogFile ) ;
    
    nbasis = sqrt( lmo ) ;
  
  
    rewind( penergyFile ) ;
    
    itmp = flength( penergyFile ) ;
    
    if( itmp != nbasis )
    {
      printf("\nThere are %d numbers in your energy file while there are %d orbitals ...\n" , itmp , nbasis ) ;
      
      exit( 91 ) ;
    }
  
  
  }
  
  
  
  // ---> MO Coefficient @ Current Geometry
  
  int irow , icol , iblock , idMO ;

  int moPerBlock , moLeftOver , nblock ; // = 15 ;
  
  cndoMO = ( double * ) calloc( lmo , sizeof( double ) ) ; dzeros( lmo , 1 , cndoMO ) ;
  
  MOEnergy = ( double * ) calloc( nbasis , sizeof( double ) ) ; // in eV unit ...
  
  dzeros( nbasis , 1 , MOEnergy ) ;




  if( moType == 9 )
  {
    
    moPerBlock = 15 ;
  
    moLeftOver = nbasis % moPerBlock ;
  
    nblock = ( nbasis - moLeftOver ) / moPerBlock ;
    

    if( moLeftOver != 0 )
    {
      nblock = nblock + 1 ;
    }      
    else
    {
      moLeftOver = moPerBlock ;
    }
  
    printf("\nThere are %d blocks and %d orbitals leftover ... \n" , nblock , moLeftOver ) ;
    
    
    
    rewind( pcndoLogFile ) ;
    
    fsearch( pcndoLogFile , "H(D)OMO=" ) ;
  
    for( iblock = 0 ; iblock < nblock - 1 ; iblock ++ )
    {
      //printf("\n# %d Block ... \n" , iblock ) ;
    
      fsearch( pcndoLogFile , "eigvals(Ev)" ) ;
    
      for( icol = 0 ; icol < moPerBlock ; icol ++ )
      {
        fscanf( pcndoLogFile , "%lf" , MOEnergy + iblock * moPerBlock + icol ) ;
    
      }
    
      fsearch( pcndoLogFile , "shell" ) ;
    
      fskip( pcndoLogFile , 3 ) ;
    
      for( irow = 0 ; irow < nbasis ; irow ++ )
      {
        fscanf( pcndoLogFile , "%s" , stmp ) ; // AO#
      
        fscanf( pcndoLogFile , "%s" , stmp ) ; // Atom#
      
        fscanf( pcndoLogFile , "%s" , stmp ) ; // AtomName
      
        fscanf( pcndoLogFile , "%s" , stmp ) ; // AOName
      
        for( icol = 0 ; icol < moPerBlock ; icol ++ )
        {
          fscanf( pcndoLogFile , "%lf" , cndoMO + irow * nbasis + ( iblock * moPerBlock + icol ) ) ;
 
        }

      }
  
  
    }
  
    fsearch( pcndoLogFile , "eigvals(Ev)" ) ;
  
    for( icol = 0 ; icol < moLeftOver ; icol ++ )
    {
      fscanf( pcndoLogFile , "%lf" , MOEnergy + ( nblock - 1 ) * moPerBlock + icol ) ;
    }
  
    fsearch( pcndoLogFile , "shell" ) ;
  
    fskip( pcndoLogFile , 3 ) ;
  
    for( irow = 0 ; irow < nbasis ; irow ++ )
    {
      fscanf( pcndoLogFile , "%s" , stmp ) ; // AO#

      fscanf( pcndoLogFile , "%s" , stmp ) ; // Atom#

      fscanf( pcndoLogFile , "%s" , stmp ) ; // AtomName

      fscanf( pcndoLogFile , "%s" , stmp ) ; // AOName

      for( icol = 0 ; icol < moLeftOver ; icol ++ )
      { 
        fscanf( pcndoLogFile , "%lf" , cndoMO + irow * nbasis + ( ( nblock - 1 ) * moPerBlock + icol ) ) ;
      }
  
    }
  
  /*
  debug = fopen( "cndoMO.deb" , "wb+" ) ;
  
  doutput( debug , nbasis , nbasis , cndoMO ) ;
  
  fclose( debug ) ;
  */
  
  
  }
  else if( moType == 11 )
  {
    rewind( pcndoLogFile ) ;
    
    fload( pcndoLogFile , cndoMO ) ;
  
    rewind( penergyFile ) ;
    
    fload( penergyFile , MOEnergy ) ;
  
  }


  dtranspose( nbasis , cndoMO , cndoMO ) ;
  
  // ===> Assemble into G09 .fchk format
  
  
  int g09MoPerBlock = 5 ;
  
  int g09NBlockMO , g09NBlockEnergy ;
  
  int g09LeftoverMO , g09LeftoverEnergy ;
    
  int iline , iload ;
  
  
  
  
  g09LeftoverMO = lmo % g09MoPerBlock ;
  
  if( g09LeftoverMO == 0 )
  {
    g09LeftoverMO = g09MoPerBlock ;
  }
  
  g09NBlockMO = ( lmo - g09LeftoverMO ) / g09MoPerBlock + 1 ;
  
  g09LeftoverEnergy = nbasis % g09MoPerBlock ;
  
  if( g09LeftoverEnergy == 0 )
  {
    g09LeftoverEnergy = g09MoPerBlock ;
  }
  
  g09NBlockEnergy = ( nbasis - g09LeftoverEnergy ) / g09MoPerBlock + 1 ;
  

  poutputFile = fopen( outputFileName , "wb+" ) ;
  
  
  
  if( exFCHK == NO )
  {
  
    // => MO ...
    
    for( iblock = 0 ; iblock < g09NBlockMO - 1 ; iblock ++ )
    {
    
      for( icol = 0 ; icol < g09MoPerBlock ; icol ++ )
      {
        fprintf( poutputFile , "% 16.8E" , *( cndoMO + iblock * g09MoPerBlock + icol ) ) ;
      }
    
      fprintf( poutputFile , "\n");
  
    }
    
  
  
    for( icol = 0 ; icol < g09LeftoverMO ; icol ++ )
    {
      fprintf( poutputFile , "% 16.8E" , *( cndoMO + ( g09NBlockMO - 1 ) * g09MoPerBlock + icol ) ) ;
    }
  
    fclose( poutputFile ) ;
  
  
  
    // => Energy ... 

  
    debug = fopen( "fakeG09Energy.deb" , "wb+" ) ;
  
    for( iblock = 0 ; iblock < g09NBlockEnergy - 1 ; iblock ++ )
    {
    
      for( icol = 0 ; icol < g09MoPerBlock ; icol ++ )
      {
        fprintf( debug , "% 16.8E" , *( MOEnergy + iblock * g09MoPerBlock + icol ) ) ;
      }
    
      fprintf( debug , "\n");
  
    }


    for( icol = 0 ; icol < g09LeftoverEnergy ; icol ++ )
    {
      fprintf( debug , "% 16.8E" , *( MOEnergy + ( g09NBlockEnergy - 1 ) * g09MoPerBlock + icol ) ) ;
    }
  
    fclose( debug ) ;
  
  
  }
  else if( exFCHK == YES )
  {
    rewind( pfchkFile ) ;
    
    iload = 0 ; 

    // => Energy ... 
    
    while( ( info = freadline( buffer , MAXCHARINLINE , pfchkFile , '!' ) ) != 0 )
    {
      strpickword( buffer , 1 , cache ) ;
      
      if( strcmp( cache , "Alpha" ) != 0 )
      {
        fprintf( poutputFile , "%s\n" , buffer ) ;
        
        //iload ++ ;
      
      } 
      else
      {
        fprintf( poutputFile , "%s\n" , buffer ) ;
        
        printf("\nFound The Alpha!!!\n\n") ;
        
        for( iblock = 0 ; iblock < g09NBlockEnergy - 1 ; iblock ++ )
        {
    
          for( icol = 0 ; icol < g09MoPerBlock ; icol ++ )
          {
            fprintf( poutputFile , "% 16.8E" , *( MOEnergy + iblock * g09MoPerBlock + icol ) ) ;
          }
    
          fprintf( poutputFile , "\n");
  
        }
        
        for( icol = 0 ; icol < g09LeftoverEnergy ; icol ++ )
        {
          fprintf( poutputFile , "% 16.8E" , *( MOEnergy + ( g09NBlockEnergy - 1 ) * g09MoPerBlock + icol ) ) ;
        }
        
        fprintf( poutputFile , "\n");
        
        break ;
      
      
      
      }

    } 
    
    fskip( pfchkFile , g09NBlockEnergy ) ;
    
    // => MO ...
  
    while( ( info = freadline( buffer , MAXCHARINLINE , pfchkFile , '!' ) ) != 0 )
    {
      strpickword( buffer , 1 , cache ) ;
      
      if( strcmp( cache , "Alpha" ) != 0 )
      {
        fprintf( poutputFile , "%s\n" , buffer ) ;
        
        //iload ++ ;
      
      } 
      else
      {
        fprintf( poutputFile , "%s\n" , buffer ) ;
        
        printf("\nFound The Alpha!!!\n\n") ;
        
        for( iblock = 0 ; iblock < g09NBlockMO - 1 ; iblock ++ )
        {
    
          for( icol = 0 ; icol < g09MoPerBlock ; icol ++ )
          {
            fprintf( poutputFile , "% 16.8E" , *( cndoMO + iblock * g09MoPerBlock + icol ) ) ;
          }
    
          fprintf( poutputFile , "\n");
  
        }
    

        for( icol = 0 ; icol < g09LeftoverMO ; icol ++ )
        {
          fprintf( poutputFile , "% 16.8E" , *( cndoMO + ( g09NBlockMO - 1 ) * g09MoPerBlock + icol ) ) ;
        }
        
        fprintf( poutputFile , "\n");
  
        break ;
      }

    } 
    
    fskip( pfchkFile , g09NBlockMO ) ;
    
    // => The Rest of FCHK File ...
  
    while( ( info = freadline( buffer , MAXCHARINLINE , pfchkFile , '!' ) ) != 0 )
    {
      //strpickword( buffer , 1 , cache ) ;
      
      fprintf( poutputFile , "%s\n" , buffer ) ;
    
    } 
    
  
  

  
  
  
  }
  
  
  
  
  
  
  
  
  
  
  return( 0 ) ;






}



