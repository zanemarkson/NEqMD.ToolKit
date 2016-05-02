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
  
  char refMOFileName[ 100 ] , refCIListName[ 100 ] , cndoLogFileName[ 100 ] , cndoTRDFileName[ 100 ] , outputFileName[ 100 ] , outputGMHName[ 100 ] ;

  int len_cndoLogFileName ;
  
  FILE * prefMOFile , * prefCIListFile , * pcndoLogFile , * pcndoTRDFile , * poutputFile ;
  
  double * refMO ; 
  
  int * refCIList ; 
  
  double * cndoMO  , * MOEnergy ;
  
  int refMOFileLength , refCIListFileLength , lmo ; // lmo is just for convenience = nbasis^2
  
  int HOMO , LUMO ;
  
  int nbasis , natom , nelectron , ncistate ;
  
  int ibasis , iatom , icistate , icidet ;
  
  int iload = 0 ; int iline = 0 ; int ieqv = 0 ;
  
  int nword , iword ;
  
  
  int debuggingMode = NO ;
  

  char buffer[ MAXCHARINLINE ] ;
  
  char cache[ MAXCHARINLINE ] ;

  
  int itmp , itmp2 ; 
  
  double dtmp , dtmp2 ; 
  
  char ctmp , tmp_char ; 
  
  char stmp[ 150 ] , tmpString[ 150 ] ;
  
  double dtmpArray[ 150 ] ;
  
  int info , signal , blank_signal ;
  
  FILE * debug ; 
  
  // ----------------------------------> Recording Command-Line Arguments ... <---------------------------------- //
  
  time_t current_time;

  time( &current_time );

  char now[ 300 ] ;

  strcpy( now , ctime( &current_time ) );

  int lennow = strlen( now ) ;

  *( now + lennow - 1 ) = ' ';

    
    
    printf("\n**********************************************************************\n");
      printf("* G_CNDO2GMH_D : Calculate GMH Type HDA for Snapshots in MD Calcs.   *\n");
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
      printf("* Note : snapCIList is labeled in snap MO number.                    *\n");
      printf("*                                                                    *\n");
      printf("*                                                                    *\n");
      printf("**********************************************************************\n");

 
  
  // =====> Defaults ...
  
  strcpy( refMOFileName , "ref.MO" ) ;
  
  strcpy( refCIListName , "ref.CI" ) ;
  
  strcpy( cndoLogFileName , "system.snap.log" ) ;
  
  strcpy( cndoTRDFileName , "system.snap.trd" ) ;
  
  //strcpy( outputGMHName , "system.gmh" ) ;
  
  
  

  
  
  // =====> Parsing cmd-line arguments ...
  
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
	      
	      case 'm' : strcpy( refMOFileName , *( ++ pcmd ) ) ; 
			 
			         printf("\nCommand-line argument indicates : Input ref. MO coefficient File name : %s ...\n" , refMOFileName ); 
	      
	                 icmd = icmd + 2 ; 
	                 
	                 break ;
  
          case 'L' : strcpy( refCIListName , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Input ref. CI excitation File name : %s ...\n" , refCIListName ); 
	         
	                 icmd = icmd + 2 ;
	         
	                 break ;

	      case 'l' : strcpy( cndoLogFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Input cndo Log file name : %s ...\n" , cndoLogFileName ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ; 
	      
	      case 't' : strcpy( cndoTRDFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Input cndo transition dipole file name : %s ...\n" , cndoTRDFileName ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ; 

	      case 'o' : strcpy( outputFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Output File name : %s ...\n" , outputFileName ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 //exgro = 19 ;
	                 
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
	    
	      
	      /*           
	      case 'p' : strcpy( pGroupName , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Group P Name is : %s ...\n" , pGroupName ); 
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;
	                 
	      case 'q' : strcpy( qGroupName , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Group Q Name is : %s ...\n" , qGroupName ); 
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;
	                 
	      case 'N' : referenceAtomNumber = atoi( *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : # %d atom will be used as the L-shape reference ...\n" , referenceAtomNumber );

                     icmd = icmd + 2 ;
                     
                     break;

	      case 's' : natomsoluteSelect = atoi( *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : The first %d atom will be treated as solute ...\n" , natomsoluteSelect );

                     icmd = icmd + 2 ;
                     
                     exs = 19 ;
                     
                     break;
          */


	      case 'h' : printf("\nUsage:  %s [ -m 'input ref. MO Coefficient file name' ] [ -L 'Ref. CI list file name' ] [ -l 'Input cndo log file name' ] [ -t 'input cndo-Calculated transition dipole file name ' ]  \n\n" , * argv ); 
	                 
	                 //printf("\n\n==> NOTE : In order to be fool-proof , if there is mis-match on the NAtom between dxdr file and any other file , this code will set NAtom to be the number from dxdr automatically ...\n\n");
	                 
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
  
  if( ( prefMOFile = fopen( refMOFileName , "r" ) ) == NULL )
  {
    printf("\nUser specified reference MO coefficient NOT FOUND ...\n");
    
    exit( 63 ) ;
  
  }
  
  if( ( prefCIListFile = fopen( refCIListName , "r" ) ) == NULL )
  {
    printf("\nUser specified reference CI List NOT FOUND ...\n");
    
    exit( 63 ) ;
  
  }
  
  
  if( ( pcndoLogFile = fopen( cndoLogFileName , "r" ) ) == NULL )
  {
    printf("\nUser specified CNDO Log File NOT FOUND ...\n");
    
    exit( 63 ) ;
  
  }
  
  if( ( pcndoTRDFile = fopen( cndoTRDFileName , "r" ) ) == NULL )
  {
    printf("\nUser specified CNDO TRD File NOT FOUND ...\n");
    
    exit( 63 ) ;
  
  }
  
  
  
  
  // =====> Loading Reference MO Coefficients ... <===== //
  
  
  refMOFileLength = flength( prefMOFile ) ;
  
  rewind( prefMOFile ) ;
  
  refMO = ( double * ) calloc( refMOFileLength , sizeof( double ) ) ;
  
  fload( prefMOFile , refMO );
  
  
  
  // =====> Loading Reference CI List ... <===== //
  
  // const int leRefCINumber , leRefInitialOrbital , leRefFinalOrbital , ctRefCINumber , ctRefInitialOrbital , ctRefFinalOrbital ;
  
  int degenerateCase ;
  
  int leRefCINumber , leRefInitialOrbital , leRefFinalOrbital , ctRefCINumber , ctRefInitialOrbital , ctRefFinalOrbital ;
  
  int * degenerateStatesList[ 4 ] ;
  
  
  int  * leRefInitialOrbitalEqv , * leRefFinalOrbitalEqv , * ctRefInitialOrbitalEqv , * ctRefFinalOrbitalEqv ;
  
  int  EXleRefInitialOrbitalEqv ,  EXleRefFinalOrbitalEqv ,  EXctRefInitialOrbitalEqv ,  EXctRefFinalOrbitalEqv ;
  
  int  lenLERefInitialOrbitalEqv , lenLERefFinalOrbitalEqv , lenCTRefInitialOrbitalEqv , lenCTRefFinalOrbitalEqv ;
  

  EXleRefInitialOrbitalEqv = 0 ;
  
  EXleRefFinalOrbitalEqv = 0 ;

  EXctRefInitialOrbitalEqv = 0 ;

  EXctRefFinalOrbitalEqv = 0 ;

  
  refCIListFileLength = flength( prefCIListFile ) ;
  
  refCIList = ( int * ) calloc( 6 , sizeof( int ) ) ;
  
  for( itmp = 0 ; itmp < 6 ; itmp ++ ) *( refCIList + itmp ) = 0 ;
  
  
  rewind( prefCIListFile ) ;

  if( refCIListFileLength == 6 )
  {
    degenerateCase = NO ;
    
    printf("\nOKay ... Regular CI List File ... \n") ;
  
    for( itmp = 0 ; itmp < 6 ; itmp ++ )
    {
      fscanf( prefCIListFile , "%d" , refCIList + itmp ) ;
    }
    
    leRefCINumber = *( refCIList + 0 ) ;
  
    leRefInitialOrbital = *( refCIList + 1 ) ;
  
    leRefFinalOrbital = *( refCIList + 2 ) ;
  
    ctRefCINumber = *( refCIList + 3 ) ;
  
    ctRefInitialOrbital = *( refCIList + 4 ) ;
  
    ctRefFinalOrbital = *( refCIList + 5 ) ;
    
    for( iload = 0 ; iload < 4 ; iload ++ )
    {
      *( degenerateStatesList + iload ) = ( int * ) calloc( 1 , sizeof( int ) ) ;
    }
    
    leRefInitialOrbitalEqv = *( degenerateStatesList + 0 ) ; *( leRefInitialOrbitalEqv + 0 ) = leRefInitialOrbital ; lenLERefInitialOrbitalEqv = 1 ; // EXleRefInitialOrbitalEqv = YES ;

    leRefFinalOrbitalEqv = *( degenerateStatesList + 1 ) ; *( leRefFinalOrbitalEqv + 0 ) = leRefFinalOrbital ; lenLERefFinalOrbitalEqv = 1 ; 

    ctRefInitialOrbitalEqv = *( degenerateStatesList + 2 ) ; *( ctRefInitialOrbitalEqv + 0 ) = ctRefInitialOrbital ; lenCTRefInitialOrbitalEqv = 1 ;

    ctRefFinalOrbitalEqv = *( degenerateStatesList + 3 ) ; *( ctRefFinalOrbitalEqv + 0 ) = ctRefFinalOrbital ; lenCTRefFinalOrbitalEqv = 1 ;

  }
  else
  {
    printf("\nAll right ... degenerate case then ... \n") ;

    degenerateCase = YES ;
    
    for( itmp = 0 ; itmp < 6 ; itmp ++ )
    {
      fscanf( prefCIListFile , "%d" , refCIList + itmp ) ;
    }
    
    
    leRefCINumber = *( refCIList + 0 ) ;
  
    leRefInitialOrbital = *( refCIList + 1 ) ;
  
    leRefFinalOrbital = *( refCIList + 2 ) ;
  
    ctRefCINumber = *( refCIList + 3 ) ;
  
    ctRefInitialOrbital = *( refCIList + 4 ) ;
  
    ctRefFinalOrbital = *( refCIList + 5 ) ;
  
    
    printf("\nIn Reference , \n\n# %d CI : LE : %d -> %d \n\n # %d CI : CT : %d -> %d \n\n" , leRefCINumber , leRefInitialOrbital , leRefFinalOrbital , ctRefCINumber , ctRefInitialOrbital , ctRefFinalOrbital ) ;
    
    
    iload = 0 ; iline = 0 ;
  
    while( ( info = freadline( buffer , MAXCHARINLINE , prefCIListFile , '!' ) ) != 0 )
    {
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
          printf("\nThis is a comment line ... So nothing will be loaded ...\n");
        
          //fskip( pinputfile , 1 );
        
          continue ;
        }
        else
        {
          printf("\nLine reads : %s ...\n" , buffer );
              
          nword = inLineWC( buffer ) ;
          
          *( degenerateStatesList + iload ) = ( int * ) calloc( nword - 1 , sizeof( int ) ) ;
          
          strpickword( buffer , 1 , cache ) ;
          
          itmp = atoi( cache ) ;
          
          strpickword( buffer , 2 , cache ) ;
          
          if( strcmp( cache , "=" ) != 0 )
          {
            printf("\nCI List File format wrong ...\n");
            
            exit( 19 ) ;
          }  
          
          if( itmp == leRefInitialOrbital && EXleRefInitialOrbitalEqv == NO )
          {
            printf("\nLE , initial , equivalent to : \n") ;
            
            leRefInitialOrbitalEqv = *( degenerateStatesList + iload ) ;
            
            lenLERefInitialOrbitalEqv = nword - 1 ;
            
            *( leRefInitialOrbitalEqv + 0 ) = leRefInitialOrbital ;
            
            for( iword = 2 ; iword < nword ; iword ++ ) 
            {
              strpickword( buffer , iword + 1 , cache ) ;
              
              *( leRefInitialOrbitalEqv + iword - 1 ) = atoi( cache ) ;
              
              printf( "\n# %d\n" , atoi( cache ) ) ;
            }
          
            EXleRefInitialOrbitalEqv = YES ;
            
            printf("\nFour indicators are : EXleRefInitialOrbitalEqv = %d , EXleRefFinalOrbitalEqv = %d , EXctRefInitialOrbitalEqv = %d , EXctRefFinalOrbitalEqv = %d " , EXleRefInitialOrbitalEqv , \
                 EXleRefFinalOrbitalEqv , EXctRefInitialOrbitalEqv , EXctRefFinalOrbitalEqv ) ;
          
          }
          else if( itmp == leRefFinalOrbital && EXleRefFinalOrbitalEqv == NO )
          {
            printf("\nLE , final , equivalent to : \n") ;
            
            leRefFinalOrbitalEqv = *( degenerateStatesList + iload ) ;
            
            lenLERefFinalOrbitalEqv = nword - 1 ;
            
            *( leRefFinalOrbitalEqv + 0 ) = leRefFinalOrbital ;
            
            for( iword = 2 ; iword < nword ; iword ++ ) 
            {
              strpickword( buffer , iword + 1 , cache ) ;
              
              *( leRefFinalOrbitalEqv + iword - 1 ) = atoi( cache ) ;
              
              printf( "\n# %d\n" , atoi( cache ) ) ;
            }
            
            EXleRefFinalOrbitalEqv = YES ;
            
            printf("\nFour indicators are : EXleRefInitialOrbitalEqv = %d , EXleRefFinalOrbitalEqv = %d , EXctRefInitialOrbitalEqv = %d , EXctRefFinalOrbitalEqv = %d " , EXleRefInitialOrbitalEqv , \
                 EXleRefFinalOrbitalEqv , EXctRefInitialOrbitalEqv , EXctRefFinalOrbitalEqv ) ;
          
          } 
          else if( itmp == ctRefInitialOrbital && EXctRefInitialOrbitalEqv == NO )
          {
            printf("\nCT , initial , equivalent to : \n") ;
            
            ctRefInitialOrbitalEqv = *( degenerateStatesList + iload ) ;
            
            lenCTRefInitialOrbitalEqv = nword - 1 ;
            
            *( ctRefInitialOrbitalEqv + 0 ) = ctRefInitialOrbital ;
            
            for( iword = 2 ; iword < nword ; iword ++ ) 
            {
              strpickword( buffer , iword + 1 , cache ) ;
              
              *( ctRefInitialOrbitalEqv + iword - 1 ) = atoi( cache ) ;
              
              printf( "\n# %d\n" , atoi( cache ) ) ;
            }
              
            EXctRefInitialOrbitalEqv = YES ;
          
            printf("\nFour indicators are : EXleRefInitialOrbitalEqv = %d , EXleRefFinalOrbitalEqv = %d , EXctRefInitialOrbitalEqv = %d , EXctRefFinalOrbitalEqv = %d " , EXleRefInitialOrbitalEqv , \
                 EXleRefFinalOrbitalEqv , EXctRefInitialOrbitalEqv , EXctRefFinalOrbitalEqv ) ;
          
          }
          else if( itmp == ctRefFinalOrbital && EXctRefFinalOrbitalEqv == NO )
          {
            printf("\nCT , final , equivalent to : \n") ;
            
            ctRefFinalOrbitalEqv = *( degenerateStatesList + iload ) ;
            
            lenCTRefFinalOrbitalEqv = nword - 1 ;
            
            *( ctRefFinalOrbitalEqv + 0 ) = ctRefFinalOrbital ;
            
            for( iword = 2 ; iword < nword ; iword ++ ) 
            {
              strpickword( buffer , iword + 1 , cache ) ;
              
              *( ctRefFinalOrbitalEqv + iword - 1 ) = atoi( cache ) ;
              
              printf( "\n# %d\n" , atoi( cache ) ) ;
            }
            
            EXctRefFinalOrbitalEqv = YES ;
            
            printf("\nFour indicators are : EXleRefInitialOrbitalEqv = %d , EXleRefFinalOrbitalEqv = %d , EXctRefInitialOrbitalEqv = %d , EXctRefFinalOrbitalEqv = %d " , EXleRefInitialOrbitalEqv , \
                 EXleRefFinalOrbitalEqv , EXctRefInitialOrbitalEqv , EXctRefFinalOrbitalEqv ) ;
          
          }
          else
          {
            printf("\nCI List File format wrong . Degeneracy specification should be of format : [orbital #] = [orbital #] [orbital #] ... \n\n");
            
            exit( 602 ) ;
          
          }

          iload ++ ;
        
        }
      
      //printf("\n%s\n" , buffer );
      }
      else
      {
        printf("\nSomething is wrong with the Telling-Blank-Line part ...\n");
      
        exit( 618 );
      }
    
    
    }
  
    if( iload != 4 )
    {
      printf("\nWrong Format in CI List File . All 4 Relevant Orbitals Must Appear In CI List Degneracy Specification Section ...\n") ;
      
      exit( 646 ) ;
    }
  
  }
  


  if( debuggingMode == YES )
  {
    debug = fopen( "eqv.deb" , "wb+" ) ;
  
    fprintf( debug , "\n\n===> LE , initial <===\n\n") ;
  
    for( ieqv = 0 ; ieqv < lenLERefInitialOrbitalEqv ; ieqv ++ )
    {
      fprintf( debug , "%d\t" , *( leRefInitialOrbitalEqv + ieqv ) ) ;
    
    }

    fprintf( debug , "\n\n===> LE , final <===\n\n") ;
  
    for( ieqv = 0 ; ieqv < lenLERefFinalOrbitalEqv ; ieqv ++ )
    {
      fprintf( debug , "%d\t" , *( leRefFinalOrbitalEqv + ieqv ) ) ;
    
    }

    fprintf( debug , "\n\n===> CT , initial <===\n\n") ;
  
    for( ieqv = 0 ; ieqv < lenCTRefInitialOrbitalEqv ; ieqv ++ )
    {
      fprintf( debug , "%d\t" , *( ctRefInitialOrbitalEqv + ieqv ) ) ;
    
    }
    
    
    fprintf( debug , "\n\n===> CT , final <===\n\n") ;
    
    for( ieqv = 0 ; ieqv < lenCTRefFinalOrbitalEqv ; ieqv ++ )
    {
      fprintf( debug , "%d\t" , *( ctRefFinalOrbitalEqv + ieqv ) ) ;
    
    }
    
    
    fclose( debug ) ;
  }
  
  
  
  
  
  // =====> Loading Information from CNDO Log File and TRD File ... <===== //
  
  // ---> # of Electrons
  
  rewind( pcndoLogFile ) ;
  
  fsearch( pcndoLogFile , "electrons=" ) ;
  
  fscanf( pcndoLogFile , "%d" , &nelectron );
  
  if( ( nelectron % 2 ) == 1 )
  {
    HOMO = ( nelectron - 1 ) / 2 + 1 ;
    
    LUMO = HOMO + 1 ;
  }
  else
  {
    HOMO = nelectron / 2 ;
    
    LUMO = HOMO + 1 ;
  
  }
  

  
  // ---> NAtom
  
  /*
  fsearch( pcndoLogFile , "Net" ) ;
  
  fskip( pcndoLogFile , 1 );
  
  fscanf( pcndoLogFile , "%s" , stmp );
  
  itmp2 = strlen( stmp );
  
  for( itmp = 5 ; itmp < itmp2 ; itmp ++ )
  {
    *( tmpString + itmp - 5 ) = *( stmp + itmp ) ;
  
  }
  
  *( stmp + itmp ) = '\0' ;
  
  natom = atof( tmpString );
  */
  /*
  rewind( pcndoLogFile ) ;
  
  fsearch( pcndoLogFile , "RealAtoms=" ) ;
  
  fscanf( pcndoLogFile , "%d" , &natom );
  */
  
  // ---> NBasis
  
  rewind( pcndoLogFile ) ;
  
  fsearch( pcndoLogFile , "Basis" ) ;
  
  fscanf( pcndoLogFile , "%s" , stmp );
  
  if( strcmp( stmp , "fns" ) == 0 )
  {
    fscanf( pcndoLogFile , "%s" , stmp ) ; 
   
    fscanf( pcndoLogFile , "%d" , &nbasis );
  }
  else
  {
    fsearch( pcndoLogFile , "SHELL" );
    
    fsearch( pcndoLogFile , "TOTAL" ) ;
    
    fscanf( pcndoLogFile , "%d" , &nbasis ) ;
  
  }
  
  
  lmo = nbasis * nbasis ;
  
  if( refMOFileLength != nbasis * nbasis )
  {
    printf("\nSomething is wrong with the # of basis ... Log file say NBasis is %d but there are %d numbers in reference MO Coefficients...\n" , nbasis , refMOFileLength );
    
    printf("\nIt is also possible that this is related to the length of reference MO coefficients ...\n");
    
    exit( 145 );
  
  }
  
  
  
  
    // ---> NCIStates
    
    rewind( pcndoLogFile ) ;
    
    
    /*
    fsearch( pcndoLogFile , "max" );
    
    fscanf( pcndoLogFile , "%s" , stmp ) ;  fscanf( pcndoLogFile , "%s" , stmp ) ;
    
    fscanf( pcndoLogFile , "%s" , stmp ) ;
    
    if( strcmp( stmp , "states=" ) == 0 )
    {
      fscanf( pcndoLogFile , "%d" , &ncistate );
    }
    else
    {
      rewind( pcndoLogFile ) ;
    
      fsearch( pcndoLogFile , "states=" );
      
      fscanf( pcndoLogFile , "%d" , &ncistate );
    }
    */

  if( fsearch( pcndoLogFile , "\"MAX_CI=") == 1 )
  {
    fscanf( pcndoLogFile , "%d" , &ncistate );
    
    printf("\nPer user's request, en toto %d CI states are calculated ...\n\n" , ncistate ) ;
  }
  else
  {
    ncistate = 60 ;
    
    printf("\nBy default , 60 CI roots are solved ...\n");
  } 
  
        
  rewind( pcndoLogFile ) ;
  
    
    // ---> MO Coefficient @ Current Geometry
    
    int irow , icol , iblock , idMO ;

    int moPerBlock = 15 ;
    
    int moLeftOver = nbasis % moPerBlock ;
    
    int nblock = ( nbasis - moLeftOver ) / moPerBlock ;
    
    if( moLeftOver != 0 )
    {
      nblock = nblock + 1 ;
    }      
    else
    {
      moLeftOver = moPerBlock ;
    }
    
    

    rewind( pcndoLogFile ) ;
    
    cndoMO = ( double * ) calloc( lmo , sizeof( double ) ) ; dzeros( lmo , 1 , cndoMO ) ;
    
    MOEnergy = ( double * ) calloc( nbasis , sizeof( double ) ) ; // in eV unit ...
    

    
    fsearch( pcndoLogFile , "H(D)OMO=" ) ;
    
    for( iblock = 0 ; iblock < nblock - 1 ; iblock ++ )
    {
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



    if( debuggingMode == YES )
    {
      debug = fopen( "cndoMO.deb" , "wb+" ) ;
    
      doutput( debug , nbasis , nbasis , cndoMO ) ;
    
      fclose( debug ) ;
    }


    // ---> CI Excitation Determinates List @ current GEOM ( snapCIList ) ;
    

    
    double * snapDetList = calloc( 4 * ncistate , sizeof( double ) ) ;
    
    rewind( pcndoLogFile ) ;    
    
    fsearch( pcndoLogFile , "lowest" ) ;
    
    fscanf( pcndoLogFile , "%d" , &itmp ) ;
    
    if( itmp != ncistate )
    {
      printf("\ncndo output is not consistent with your choice of MAX_CI . Please check ...\n") ;
      
      exit( 486 ) ;
    
    }
    
    fsearch( pcndoLogFile , "strength" ) ;
    
    *( snapDetList + 0 * 4 + 0 ) = 1 ;
    
    *( snapDetList + 0 * 4 + 1 ) = HOMO ;
    
    *( snapDetList + 0 * 4 + 2 ) = HOMO ; 
    
    *( snapDetList + 0 * 4 + 3 ) = 0 ; 
    
    fsearch( pcndoLogFile , "1" ) ;
    
    fskip( pcndoLogFile , 1 ) ;
    
    
    for( icidet = 1 ; icidet < ncistate ; icidet ++ )
    {
      *( snapDetList + icidet * 4 + 0 ) = icidet + 1 ; 
      
      fsearch( pcndoLogFile , "(" ) ;
      
      fscanf( pcndoLogFile , "%lf" , snapDetList + icidet * 4 + 1 ) ; // Read-in ... so human label...
      
      fscanf( pcndoLogFile , "%s" , stmp ) ;
      
      fscanf( pcndoLogFile , "%lf" , snapDetList + icidet * 4 + 2 ) ;
    
    }


    if( debuggingMode == YES )
    {
      debug = fopen( "snapDetList.deb" , "wb+" ) ;
      
      
      for( icistate = 0 ; icistate < ncistate ; icistate ++ )
      {
        fprintf( debug , "\n%d\t%d\t--->\t%d\t\n" , ( int ) *( snapDetList + icistate * 4 + 0 ) , ( int ) ( *( snapDetList + icistate * 4 + 1 ) ) , ( int ) ( *( snapDetList + icistate * 4 + 2 ) ) ) ;
      }
      
      
      //doutput( debug , ncistate , 4 , snapDetList ) ;
    
      fclose( debug ) ;
    }



    // ---> CI Coefficients @ current GEOM ( snapCIList ) ;
    
    int nnout = 12 ; // following the CNDO code variable name ;

    double * ciCoefficients = calloc( ncistate * nnout , sizeof( double ) ) ;
    
    dzeros( ncistate , nnout , ciCoefficients ) ;
    
    int * assignmentCIStates = calloc( nnout , sizeof( int ) ) ;
    
    for( icistate = 0 ; icistate < nnout ; icistate ++ ) *( assignmentCIStates + icistate ) = NO ;
        

    //rewind( pcndoLogFile ) ;
    
    fsearch( pcndoLogFile , "determ-" ) ;
    
    fsearch( pcndoLogFile , "1" ) ;
    
    fskip( pcndoLogFile , 1 ) ;
    
    fsearch( pcndoLogFile , "1" ) ;
    
    fskip( pcndoLogFile , 1 ) ;
    
    *( ciCoefficients + 0 ) = 1.0000 ;
    

    for( icidet = 1 ; icidet < ncistate ; icidet ++ )
    {
      fscanf( pcndoLogFile , "%d" , &itmp ) ;
      
      if( itmp != icidet + 1 )
      {
        printf("\nSomething is wrong with reading the CI Coefficients ... \n");
        
        printf("\nWhile log file indicates the # %d CI config , we are filling the # %d ...\n" , itmp , icidet + 1 ) ;
        
        exit( 660 ) ;
        
      }
      
      for( icistate = 0 ; icistate < nnout ; icistate ++ )
      {
        fscanf( pcndoLogFile , "%lf" , ciCoefficients + icidet * nnout + icistate ) ;      
      }
    
   
    }
   

    if( debuggingMode == YES )
    {
      debug = fopen( "ciCoefficients.deb" , "wb+" ) ;
    
      doutput( debug , ncistate , nnout , ciCoefficients ) ;
    
      fclose( debug ) ;
    }


    // ---> CI Excitation List @ current GEOM ( snapCIList ) ;
    
    
    double * snapCIList = calloc( 4 * nnout , sizeof( double ) ) ;
    
    double * cieigvec = calloc( ncistate , sizeof( double ) ) ; 
    
    dzeros( ncistate , 1 , cieigvec ) ;
    
    int * dominDetList = calloc( nnout , sizeof( int ) ) ;
    
    izeros( nnout , 1 , dominDetList ) ;
    
    int detID ;
    

    *( snapCIList + 0 * 4 + 0 ) = 1 ;
    
    *( snapCIList + 0 * 4 + 1 ) = HOMO ;
    
    *( snapCIList + 0 * 4 + 2 ) = HOMO ;
    
    *( snapCIList + 0 * 4 + 3 ) = 0.00 ;
    
    *( dominDetList + 0 ) = 1 ;
    


    for( icistate = 1 ; icistate < nnout ; icistate ++ )
    {
      for( icidet = 0 ; icidet < ncistate ; icidet ++ )
      {
        //*( cieigvec + icidet ) = abs( *( ciCoefficients + icidet * nnout + icistate ) );
        *( cieigvec + icidet ) = fabs( *( ciCoefficients + icidet * nnout + icistate ) ) ;
        
        //printf("%lf\t" , abs( 0.500 ) ) ;

      }
      
      printf("\n\n\n") ;
      
      detID = dmaxID( ncistate , cieigvec ) + 1 ;  // Human Label ...
      
      *( snapCIList + icistate * 4 + 0 ) = icistate + 1 ;
      
      *( snapCIList + icistate * 4 + 1 ) = *( snapDetList + ( detID - 1 ) * 4 + 1 ) ;
      
      *( snapCIList + icistate * 4 + 2 ) = *( snapDetList + ( detID - 1 ) * 4 + 2 ) ;
      
      *( dominDetList + icistate ) = detID ; // See? Human Label ...
      
      
     
      printf("\nFor #%d CI state , dominant transition is ( using current snapshot determinate label ) %d ...\n" , icistate + 1 , detID ) ;
    
    }


    // ---> CI State Dipole vector and scalars ... & Transition Dipole between CI States 
    
    //DIP dipole[ ncistate * ncistate ] ;
    
    int ntrans = ncistate * ncistate ;
    
    int itrans ;
    
    int dipoleinfo ;

    
    double * dipole = calloc( 4 * ntrans , sizeof( double ) ) ;    
    
    int lengthTRDFile = flength( pcndoTRDFile ) ;
 
    if( lengthTRDFile != 6 * ncistate * ( ncistate - 1 ) ) 
    {
      printf("\nThere are %d transitions recorded in TRD file ... while ncistate  = %d ...\n\n" , lengthTRDFile / 6 , ncistate ) ;
      
      exit( 564 ) ;
    }
    
    // ---> Transition Dipoles ...
    
    rewind( pcndoTRDFile ) ;
    
    for( itrans = 0 ; itrans < ncistate * ncistate ; itrans ++ )
    {
      fscanf( pcndoTRDFile , "%d" , &itmp ) ;

      fscanf( pcndoTRDFile , "%d" , &itmp2 ) ;
      
      fscanf( pcndoTRDFile , "%lf" ,  dipole + 0 * ntrans + ( itmp - 1 ) * ncistate + ( itmp2 - 1 )  ) ;
      
      fscanf( pcndoTRDFile , "%lf" ,  dipole + 1 * ntrans + ( itmp - 1 ) * ncistate + ( itmp2 - 1 )  ) ;
      
      fscanf( pcndoTRDFile , "%lf" ,  dipole + 2 * ntrans + ( itmp - 1 ) * ncistate + ( itmp2 - 1 )  ) ;
      
      fscanf( pcndoTRDFile , "%lf" ,  dipole + 3 * ntrans + ( itmp - 1 ) * ncistate + ( itmp2 - 1 )  ) ;
      
      /*
      printf("\nBetween CI# %d and CI# %d , miuX = %lf , miuY = %lf , miuZ = %lf , |miu| = %lf ...\n" , \
      itmp , itmp2 , *( dipole + 0 * ntrans + ( itmp - 1 ) * ncistate + ( itmp2 - 1 ) ) , \
                     *( dipole + 1 * ntrans + ( itmp - 1 ) * ncistate + ( itmp2 - 1 ) ) , \
                     *( dipole + 2 * ntrans + ( itmp - 1 ) * ncistate + ( itmp2 - 1 ) ) , \
                     *( dipole + 3 * ntrans + ( itmp - 1 ) * ncistate + ( itmp2 - 1 ) ) ) ;
      */
      
      //printf("\nBetween CI# %d and CI# %d ...\n" , itmp , itmp2 ) ;
      
    }
    
    
    // ---> Ground State Dipole 
    
    rewind( pcndoTRDFile ) ;
    
    rewind( pcndoLogFile ) ;
    
    fsearch( pcndoLogFile , "hybridization" ) ;
    
    fsearch( pcndoLogFile , "total" ) ;
    
    fscanf( pcndoLogFile , "%lf" , dipole + 3 * ntrans + 0 * ncistate + 0 ) ;
    
    fscanf( pcndoLogFile , "%lf" , dipole + 0 * ntrans + 0 * ncistate + 0 ) ;
    
    fscanf( pcndoLogFile , "%lf" , dipole + 1 * ntrans + 0 * ncistate + 0 ) ;
    
    fscanf( pcndoLogFile , "%lf" , dipole + 2 * ntrans + 0 * ncistate + 0 ) ;
    
    
    
    // ---> CI States Dipoles 
    
    rewind( pcndoLogFile ) ;
    
    itmp = 0 ; itmp2 = 0 ; icistate = 0 ;
    
    //fsearch( pcndoLogFile , "---------polarization---------" ) ;
    
    fsearch( pcndoLogFile , "Emiss" ) ;
    
    while( ( dipoleinfo = freadline( buffer , MAXCHARINLINE , pcndoLogFile , ';' ) ) != 0 )
    {
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
          icistate ++ ;
          
          nword = inLineWC( buffer ) ;
          
          strpickword( buffer , 1 , cache ) ;
          
          itmp = atoi( cache ) ;
          
          if( itmp != icistate + 1 )
          {
            printf("\nATTENTION!!! Something is wrong in the CI State Dipole part in cndo Log file ...\n");
            
            printf("\nLooks like it hits => %s <=  \n" , cache ) ;
            
            exit( 655 ) ;
          }
          
          strpickword( buffer , 2 , cache ) ;
          
          dtmp = atof( cache ) ;
          
          *( dtmpArray + icistate ) = dtmp * EV2HARTREE ; //EV2HARTREE = (1.00/27.211) ;
          
          
          strpickword( buffer , nword - 3 , cache ) ;
          
          *( dipole + 3 * ntrans + icistate * ncistate + icistate ) = atof( cache ) ;
          
          
          strpickword( buffer , nword - 2 , cache ) ;
          
          *( dipole + 0 * ntrans + icistate * ncistate + icistate ) = atof( cache ) ;
          
          
          strpickword( buffer , nword - 1 , cache ) ;
          
          *( dipole + 1 * ntrans + icistate * ncistate + icistate ) = atof( cache ) ;


          strpickword( buffer , nword - 0 , cache ) ;
          
          *( dipole + 2 * ntrans + icistate * ncistate + icistate ) = atof( cache ) ;
          
          
          if( icistate == ncistate - 1 ) break ;


        }
    
    
      }
     
    }
    
    for( icistate = 1 ; icistate < nnout ; icistate ++ )
    {
      *( snapCIList + 4 * icistate + 3 ) = *( dtmpArray + icistate ) ;
    }

    if( debuggingMode == YES )
    {
      debug = fopen( "snapCIList.deb" , "wb+" ) ;
        
      for( icistate = 0 ; icistate < nnout ; icistate ++ )
      {
        fprintf( debug , "\n%d\t%d\t--->\t%d\t\t\tE = % 10.8E\n" , ( int ) ( *( snapCIList + icistate * 4 + 0 ) ) , ( int ) ( *( snapCIList + icistate * 4 + 1 ) ) , ( int ) ( *( snapCIList + icistate * 4 + 2 ) ) , *( snapCIList + icistate * 4 + 3 ) ) ;
      }

      fclose( debug ) ;
    
    }
 
    
    int idebug1 , idebug2 ;
   
    if( debuggingMode == YES )
    {
      debug = fopen( "transDipole.deb" , "wb+" ) ;
      
      for( idebug1 = 0 ; idebug1 < ncistate ; idebug1 ++ )
      {
        for( idebug2 = 0 ; idebug2 < ncistate ; idebug2 ++ )
        {
          if( idebug2 != idebug1 )
          {
            fprintf( debug , "%d\t%d\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n" , idebug1 + 1 , idebug2 + 1 , \
            *( dipole + 0 * ntrans + idebug1 * ncistate + idebug2 ) , *( dipole + 1 * ntrans + idebug1 * ncistate + idebug2 ) , \
            *( dipole + 2 * ntrans + idebug1 * ncistate + idebug2 ) , *( dipole + 3 * ntrans + idebug1 * ncistate + idebug2 ) ) ;
          }
        }
    
      }
    
      fclose( debug ) ;
    }

    if( debuggingMode == YES )
    {
      debug = fopen( "ciDipole.deb" , "wb+" ) ;
    
      for( idebug1 = 0 ; idebug1 < ncistate ; idebug1 ++ )
      {
        fprintf( debug , "%d\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n" , idebug1 + 1 , \
        *( dipole + 0 * ntrans + idebug1 * ncistate + idebug1 ) , *( dipole + 1 * ntrans + idebug1 * ncistate + idebug1 ) ,\
        *( dipole + 2 * ntrans + idebug1 * ncistate + idebug1 ) , *( dipole + 3 * ntrans + idebug1 * ncistate + idebug1 ) ) ;
    
      }
  
      fclose( debug ) ;
    
    }

    
    for( itmp = 0 ; itmp < 4 * ntrans ; itmp ++ ) // Change to AU
    {
      *( dipole + itmp ) = ( *( dipole + itmp ) ) * DEBYE2AU ;
    }
    
    

    

    
    // ---> Execute the transpose( refMO ) * cndoMO and making the mapping ... 
    
    double * orbitalOverlap = calloc( nbasis * nbasis , sizeof( double) ) ;
    
    dzeros( nbasis , nbasis , orbitalOverlap ) ;
    
    itmp2 = nbasis ;
    
    double done = 1.0000 ; double dzero = 0.0000 ;
    
    dgemm_( "N" , "T" , &itmp2 , &itmp2 , &itmp2 , &done , refMO , &itmp2 , cndoMO , &itmp2 , &dzero , orbitalOverlap , &itmp2 ) ;
    
    dtranspose( itmp2 , orbitalOverlap , orbitalOverlap ) ;
    
    // each row of orbitalOverlap corresponds to MO# in refMO , column corresponds to MO# in cndoMO 
   

    
    if( debuggingMode == YES )
    {
      debug = fopen( "orbitalOverlap.deb" , "wb+" ) ;
    
      doutput( debug , nbasis, nbasis, orbitalOverlap ) ;
    
      fclose( debug ) ; 
    }
    
    
    
    int * snapMOMapping = calloc( nbasis , sizeof( int ) ) ;
    
      
    double * moEigvec = calloc( nbasis , sizeof( double ) ) ; dzeros( nbasis , 1 , moEigvec ) ;
    
    int irefOrbital , isnapOrbital ;
    
    unsigned int currentSnapOrbitalInRef = -1 ;
    
    int neighborhood = 30 ;
    

    for( isnapOrbital = 0 ; isnapOrbital < HOMO ; isnapOrbital ++ )
    {
      dzeros( nbasis , 1 , moEigvec ) ;
      
      for( irefOrbital = 0 ; irefOrbital < HOMO ; irefOrbital ++ )
      {
        *( moEigvec + irefOrbital ) = fabs( *( orbitalOverlap + irefOrbital * nbasis + isnapOrbital ) ) ;
      }
     
      //printf("\nMO # %d ...\n" , isnapOrbital + 1 ) ;      
 
      currentSnapOrbitalInRef = dmaxID( nbasis , moEigvec ) ;
      
      *( snapMOMapping + isnapOrbital ) = currentSnapOrbitalInRef ;
      
    }
    
    
    
    for( isnapOrbital = HOMO ; isnapOrbital < nbasis ; isnapOrbital ++ )
    {
      dzeros( nbasis , 1 , moEigvec ) ;
      
      for( irefOrbital = HOMO ; irefOrbital < nbasis ; irefOrbital ++ )
      {
        *( moEigvec + irefOrbital ) = fabs( *( orbitalOverlap + irefOrbital * nbasis + isnapOrbital ) ) ;
      }
     
      //printf("\nMO # %d ...\n" , isnapOrbital + 1 ) ;      
 
      currentSnapOrbitalInRef = dmaxID( nbasis , moEigvec ) ;
      
      *( snapMOMapping + isnapOrbital ) = currentSnapOrbitalInRef ;
      
    }
    

    
    if( debuggingMode == YES )
    {
      debug = fopen( "orbitalMapping.deb" , "wb+" ) ;
    
      fprintf( debug , "\nsnapOrbital  ->  refOrbital\n" ) ;

      for( isnapOrbital = 0 ; isnapOrbital < nbasis ; isnapOrbital ++ )
      {
        fprintf( debug , "%d\t->\t%d\n" , isnapOrbital + 1 , *( snapMOMapping + isnapOrbital ) + 1 ) ;
      }
    
      fclose( debug ) ; 
    
    }
    
    
    // ----------------------> Making Hot List ... <------------------//
    
    int * leHotDetList , * ctHotDetList ;
    
    int nLEHotDet = 0 , nCTHotDet = 0 ;

    int iLEHotDet = 0 , iCTHotDet = 0 ;
    
    
    
    int currentDetInitialInRef , currentDetFinalInRef ;
    

     
    
    for( icidet = 1 ; icidet < ncistate ; icidet ++ )
    {
      currentDetInitialInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + icidet * 4 + 1 ) ) - 1 ) ) + 1 ;
      
      currentDetFinalInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + icidet * 4 + 2 ) ) - 1 ) ) + 1 ;
      
      //printf("\nIn process of building Hot Det. List. Currently checking up on No. %d Det ( %d -> %d ) \n" , icidet + 1 , currentDetInitialInRef , currentDetFinalInRef ) ;
      
      itmp = 1 ; itmp2 = 1 ;
      
      //int  lenLERefInitialOrbitalEqv , lenLERefFinalOrbitalEqv , lenCTRefInitialOrbitalEqv , lenCTRefFinalOrbitalEqv ;
      
      for( ieqv = 0 ; ieqv < lenLERefInitialOrbitalEqv ; ieqv ++ )
      {
        itmp = itmp * ( currentDetInitialInRef - *( leRefInitialOrbitalEqv + ieqv ) ) ;
      }
      
      for( ieqv = 0 ; ieqv < lenLERefFinalOrbitalEqv ; ieqv ++ )
      {
        itmp2 = itmp2 * ( currentDetFinalInRef - *( leRefFinalOrbitalEqv + ieqv ) ) ;
      }
      
      if( itmp == 0 && itmp2 == 0 && abs( ( int ) ( *( snapDetList + icidet * 4 + 1 ) ) - leRefInitialOrbital ) < neighborhood && abs( ( int ) ( *( snapDetList + icidet * 4 + 2 ) ) - leRefFinalOrbital ) <  neighborhood )
      {
        printf("\nWe found that # %d determinate is a HOT-Det for LE ...\n" , icidet + 1 ) ;
        
        iLEHotDet ++ ;
        
      }
    
    }
        
    nLEHotDet = iLEHotDet ;
        
    leHotDetList = ( int * ) calloc( nLEHotDet , sizeof( int ) ) ;
    
    iLEHotDet = 0 ;
    
    for( icidet = 1 ; icidet < ncistate ; icidet ++ )
    {
      currentDetInitialInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + icidet * 4 + 1 ) ) - 1 ) ) + 1 ;
      
      currentDetFinalInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + icidet * 4 + 2 ) ) - 1 ) ) + 1 ;
      
      // printf("\nCurrently revisiting No. %d Det ( %d -> %d ) \n" , icidet + 1 , currentDetInitialInRef , currentDetFinalInRef ) ;
      
      itmp = 1 ; itmp2 = 1 ;
        
      //int  lenLERefInitialOrbitalEqv , lenLERefFinalOrbitalEqv , lenCTRefInitialOrbitalEqv , lenCTRefFinalOrbitalEqv ;
      
      for( ieqv = 0 ; ieqv < lenLERefInitialOrbitalEqv ; ieqv ++ )
      {
        itmp = itmp * ( currentDetInitialInRef - *( leRefInitialOrbitalEqv + ieqv ) ) ;
      }
      
      for( ieqv = 0 ; ieqv < lenLERefFinalOrbitalEqv ; ieqv ++ )
      {
        itmp2 = itmp2 * ( currentDetFinalInRef - *( leRefFinalOrbitalEqv + ieqv ) ) ;
      }
      
      if( itmp == 0 && itmp2 == 0 && abs( ( int ) ( *( snapDetList + icidet * 4 + 1 ) ) - leRefInitialOrbital ) < neighborhood && abs( ( int ) ( *( snapDetList + icidet * 4 + 2 ) ) - leRefFinalOrbital ) <  neighborhood )
      { 
        *( leHotDetList + iLEHotDet ) = icidet ; // C-Labeling
        
        printf("\nFor LE Hot Det. List , now # %d det is in !!!\n" , icidet + 1 ) ;
        
        iLEHotDet ++ ;
        
      }
    
    }
    
    if( iLEHotDet != nLEHotDet )
    {
      printf("\nError making LE Hot Det List. At first glance we had %d Hot Det. But now we have only % ...\n" , nLEHotDet + 1 , iLEHotDet + 1 ) ;
    }

    
    
    iCTHotDet = 0 ; 
    
    for( icidet = 1 ; icidet < ncistate ; icidet ++ )
    {
      currentDetInitialInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + icidet * 4 + 1 ) ) - 1 ) ) + 1 ;
      
      currentDetFinalInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + icidet * 4 + 2 ) ) - 1 ) ) + 1 ;
      
      //printf("\nIn process of building Hot Det. List. Currently checking up on No. %d Det ( %d -> %d ) \n" , icidet + 1 , currentDetInitialInRef , currentDetFinalInRef ) ;
      
      itmp = 1 ; itmp2 = 1 ;
      
      //int  lenLERefInitialOrbitalEqv , lenLERefFinalOrbitalEqv , lenCTRefInitialOrbitalEqv , lenCTRefFinalOrbitalEqv ;
      
      for( ieqv = 0 ; ieqv < lenCTRefInitialOrbitalEqv ; ieqv ++ )
      {
        itmp = itmp * ( currentDetInitialInRef - *( ctRefInitialOrbitalEqv + ieqv ) ) ;
      }
      
      for( ieqv = 0 ; ieqv < lenCTRefFinalOrbitalEqv ; ieqv ++ )
      {
        itmp2 = itmp2 * ( currentDetFinalInRef - *( ctRefFinalOrbitalEqv + ieqv ) ) ;
      }
      
      if( itmp == 0 && itmp2 == 0 && abs( ( int ) ( *( snapDetList + icidet * 4 + 1 ) ) - ctRefInitialOrbital ) < neighborhood && abs( ( int ) ( *( snapDetList + icidet * 4 + 2 ) ) - ctRefFinalOrbital ) <  neighborhood )
      {
        printf("\nWe found that # %d determinate is a HOT-Det for CT ...\n" , icidet + 1 ) ;
        
        iCTHotDet ++ ;
        
      }
    
    }
        
    nCTHotDet = iCTHotDet ;
        
    ctHotDetList = ( int * ) calloc( nLEHotDet , sizeof( int ) ) ;
    
    
    iCTHotDet = 0 ;
    
    for( icidet = 1 ; icidet < ncistate ; icidet ++ )
    {
      currentDetInitialInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + icidet * 4 + 1 ) ) - 1 ) ) + 1 ;
      
      currentDetFinalInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + icidet * 4 + 2 ) ) - 1 ) ) + 1 ;
      
      // printf("\nCurrently revisiting No. %d Det ( %d -> %d ) \n" , icidet + 1 , currentDetInitialInRef , currentDetFinalInRef ) ;
      
      itmp = 1 ; itmp2 = 1 ;
        
      //int  lenLERefInitialOrbitalEqv , lenLERefFinalOrbitalEqv , lenCTRefInitialOrbitalEqv , lenCTRefFinalOrbitalEqv ;
      
      for( ieqv = 0 ; ieqv < lenCTRefInitialOrbitalEqv ; ieqv ++ )
      {
        itmp = itmp * ( currentDetInitialInRef - *( ctRefInitialOrbitalEqv + ieqv ) ) ;
      }
      
      for( ieqv = 0 ; ieqv < lenCTRefFinalOrbitalEqv ; ieqv ++ )
      {
        itmp2 = itmp2 * ( currentDetFinalInRef - *( ctRefFinalOrbitalEqv + ieqv ) ) ;
      }
      
      if( itmp == 0 && itmp2 == 0 && abs( ( int ) ( *( snapDetList + icidet * 4 + 1 ) ) - ctRefInitialOrbital ) < neighborhood && abs( ( int ) ( *( snapDetList + icidet * 4 + 2 ) ) - ctRefFinalOrbital ) <  neighborhood )
      { 
        *( ctHotDetList + iCTHotDet ) = icidet ; // C-Labeling
        
        printf("\nFor CT Hot Det. List , now # %d det is in !!!\n" , icidet + 1 ) ;
        
        iCTHotDet ++ ;
        
      }
    
    }
    
    if( iCTHotDet != nCTHotDet )
    {
      printf("\nError making CT Hot Det List. At first glance we had %d Hot Det. But now we have only % ...\n" , nCTHotDet + 1 , iCTHotDet + 1 ) ;
    }

    
    
    
    // ----------------------> Searching for the LE and CT state ... <------------------//
    
    /*
       int leRefCINumber , leRefInitialOrbital , leRefFinalOrbital , ctRefCINumber , ctRefInitialOrbital , ctRefFinalOrbital ;
       
       int * degenerateStatesList[ 4 ] ;
       
       int  * leRefInitialOrbitalEqv , * leRefFinalOrbitalEqv , * ctRefInitialOrbitalEqv , * ctRefFinalOrbitalEqv ;
       
       int  EXleRefInitialOrbitalEqv ,  EXleRefFinalOrbitalEqv ,  EXctRefInitialOrbitalEqv ,  EXctRefFinalOrbitalEqv ;
       
       int  lenLERefInitialOrbitalEqv , lenLERefFinalOrbitalEqv , lenCTRefInitialOrbitalEqv , lenCTRefFinalOrbitalEqv ;
  
    */
    
    
    int leSnapCINumber , leSnapInitialOrbital , leSnapFinalOrbital , ctSnapCINumber , ctSnapInitialOrbital , ctSnapFinalOrbital ;
    
    int leSnapCIAnotherNumber , ctSnapCIAnotherNumber;
    
    //int aliasLERefInitialOrbital , aliasLERefFinalOrbital , aliasCTRefInitialOrbital , aliasCTRefFinalOrbital ; 
    
    int currentDominateDet , currentSnapInitialInRef , currentSnapFinalInRef ;
    
    double LEFactor , CTFactor ; 
    
    int exSnapLE = 0 ; 
    
    int exSnapCT = 0 ;
    
    int canonicalIntersection = 0 ;
    
    double energyDiffOneAnother = 0.000 , energyDiffOneOne = 0.000 ;

    double CICoeffOnHotDet , preDominCICoeff ;
    
    double prevLEFactor , prevCTFactor ;
    
    double coefficientVicinity = 0.20 ;
    
    double dominateFactorThreshold = 0.71 ;
    
    double energyVicinity = 0.20 * EV2HARTREE ;
    
    
    //----------------------> LE Eigenstate Searching ... <------------------//
    
    exSnapLE = 0 ;
    
    /*
    aliasLERefInitialOrbital = *( leRefInitialOrbitalEqv + 0 ) ;
    
    aliasLERefFinalOrbital = *( leRefFinalOrbitalEqv + 0 ) ;
    
    aliasCTRefInitialOrbital = *( ctRefInitialOrbitalEqv + 0 ) ;
    
    aliasCTRefFinalOrbital = *( ctRefFinalOrbitalEqv + 0 ) ;
    */

    printf("\n----------> BEGIN Searching For The 1st Local-Excitation CI State ...\n\n") ;

    /* 1st Round : Simple Searching     

    for( icistate = 1 ; icistate < nnout ; icistate ++ )
    {
      if( *( assignmentCIStates + icistate ) == YES )
      {
        printf("\n# %d CI State has already been assigned ... So it cannot be LE state again ...\n" , icistate + 1 ) ;
        
        continue ;
      }
      
      
      currentSnapInitialInRef = *( snapMOMapping + ( int ) ( *( snapCIList + icistate * 4 + 1 ) - 1 ) ) + 1 ; // ATTENTION : here currentSnapInitialInRef is in human-label
      
      currentSnapFinalInRef = *( snapMOMapping + ( int ) ( *( snapCIList + icistate * 4 + 2 ) - 1 ) ) + 1 ;
      
      printf("\n---> Now it is the # %d CI Excitation --- In SnapBasis , %d -> %d --- same as In RefBasis %d -> %d <---\n" , icistate + 1 , ( int ) *( snapCIList + icistate * 4 + 1 ) , ( int ) *( snapCIList + icistate * 4 + 2 ) , currentSnapInitialInRef , currentSnapFinalInRef ) ;
      
      if( currentSnapInitialInRef == aliasLERefInitialOrbital && currentSnapFinalInRef == aliasLERefFinalOrbital  )
      {
        leSnapInitialOrbital = currentSnapInitialInRef ;
        
        leSnapFinalOrbital = currentSnapFinalInRef ;
        
        leSnapCINumber = icistate ;
        
        *( assignmentCIStates + icistate ) = YES ;
        
        LEFactor = fabs( *( ciCoefficients + ( *( dominDetList + icistate ) - 1 ) * nnout + icistate ) ) ;  
        
        printf("\n-----> Found It ! In snapshot , %d CI is LE <-----\n" , leSnapCINumber + 1 ) ;
        
        exSnapLE = 3 ;

        break ;
      }
      else
      {
        printf("\nEm...Not this one ... I mean %d CI is not LE ...\n" , icistate + 1 ) ;
        
        continue ;
      }
    
    }
    */
    
    
    /* 2nd Round : Scattered Multiple LE initial or final orbitals ...  
    
    itmp2 = 0 ; dtmp = 0.000 ; dtmp2 = 0.000 ;
    
    double preDominCICoeff = 0.000 ;
    
    CICoeffOnHotDet = 0.000 ;
    
    
    if( exSnapLE == 0 ) // && lenRevLERefFinalMapping > 1 ) // 2nd Round : If multiple orbital has leRefFinalOrbital nature , add them together ; 
    { 
      printf("\nLooks like we did not find the LE state from 1st Round searching ... Let's try 2nd Round ...\n") ; 
      
      for( icistate = 1 ; icistate < nnout && *( assignmentCIStates + icistate ) == NO ; icistate ++ )
      {
        
        if( *( assignmentCIStates + icistate ) == YES )
        {
          printf("\n# %d CI State has already been assigned ... So it cannot be LE state again ...\n" , icistate + 1 ) ;
          
          continue ;
        }
      

        preDominCICoeff = *( ciCoefficients + ( *( dominDetList + icistate ) - 1 ) * nnout + icistate ) ;    
        
        preDominCICoeff = fabs( preDominCICoeff ) ;
        
        dtmp = 0.000 ; 

        for( iLEHotDet = 0 ; iLEHotDet < nLEHotDet ; iLEHotDet ++ )
        {
          currentDetInitialInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( leHotDetList + iLEHotDet ) ) * 4 + 1 ) ) - 1 ) ) + 1 ;
      
          currentDetFinalInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( leHotDetList + iLEHotDet ) ) * 4 + 2 ) ) - 1 ) ) + 1 ;
          
          if( currentDetInitialInRef == aliasLERefInitialOrbital && currentDetFinalInRef == aliasLERefFinalOrbital )
          {
            CICoeffOnHotDet = *( ciCoefficients + ( *( leHotDetList + iLEHotDet ) ) * nnout + icistate ) ;
            
            dtmp = dtmp + CICoeffOnHotDet * CICoeffOnHotDet ; 
          }
        
        }

        dtmp = sqrt( dtmp ) ;
        
        
        dtmp2 = 0.000; 
        
        for( iCTHotDet = 0 ; iCTHotDet < nCTHotDet ; iCTHotDet ++ )
        {
          currentDetInitialInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( ctHotDetList + iCTHotDet ) ) * 4 + 1 ) ) - 1 ) ) + 1 ;
      
          currentDetFinalInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( ctHotDetList + iCTHotDet ) ) * 4 + 2 ) ) - 1 ) ) + 1 ;
          
          if( currentDetInitialInRef == aliasCTRefInitialOrbital && currentDetFinalInRef == aliasCTRefFinalOrbital )
          {
            CICoeffOnHotDet = *( ciCoefficients + ( *( ctHotDetList + iCTHotDet ) ) * nnout + icistate ) ;
            
            dtmp2 = dtmp2 + CICoeffOnHotDet * CICoeffOnHotDet ; 
          }
        
        }
        
        dtmp2 = sqrt( dtmp2 ) ;

        
        printf("\nPreviously , for # %d CI , dominate CI-Coefficient is on # %d Det , value is % 10.6f ...\n" , icistate + 1 , *( dominDetList + icistate ) , preDominCICoeff ) ;
        
        
        
        if( dtmp > preDominCICoeff && dtmp > dtmp2 )
        {
          printf("\n\n-----> Found It ! In snapshot , %d CI is LE <-----\n Looks like we found out that # %d CI has dominate combined LE nature ... with scaled coefficient value % 10.6f ... \n" , icistate + 1 , icistate + 1 , dtmp ) ;
          
          printf("\n@ this excited state , LE has % 10.6f combined nature and CT has % 10.6f combined nature ...\n" , dtmp , dtmp2 ) ;
          
          leSnapCINumber = icistate ;
          
          *( assignmentCIStates + icistate ) = YES ;
          
          LEFactor = dtmp ;
          
          exSnapLE = 3 ;
          
          break ;
        
        }
      
      } 
      
    }
    */
    
        
    /* 3rd Round : Add all LE nature together ...  */
        
    if( exSnapLE == 0 ) //&& lenRevLERefFinalMapping == 1 )
    {
      //printf("\nSo ... I guess either the 2nd Round died too , meaning we are really closed to canonical intersection ...\n") ;
      
      //printf("\nSo ... 3rd Round , Add all LE nature together !!!\n") ;
      
      for( icistate = 1 ; icistate < nnout && *( assignmentCIStates + icistate ) == NO ; icistate ++ )
      {
  
        if( *( assignmentCIStates + icistate ) == YES )
        {
          printf("\n# %d CI State has already been assigned ... So it cannot be LE state again ...\n" , icistate + 1 ) ;
          
          continue ;
        }
      
        currentDominateDet = *( dominDetList + icistate ) - 1 ;
        
        preDominCICoeff = *( ciCoefficients + currentDominateDet * nnout + icistate ) ;     
        
        preDominCICoeff = fabs( preDominCICoeff ) ; 
      
        dtmp = 0.000 ; 

        for( iLEHotDet = 0 ; iLEHotDet < nLEHotDet ; iLEHotDet ++ )
        {
          currentDetInitialInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( leHotDetList + iLEHotDet ) ) * 4 + 1 ) ) - 1 ) ) + 1 ;
      
          currentDetFinalInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( leHotDetList + iLEHotDet ) ) * 4 + 2 ) ) - 1 ) ) + 1 ;
          
          CICoeffOnHotDet = *( ciCoefficients + ( *( leHotDetList + iLEHotDet ) ) * nnout + icistate ) ;
          
          dtmp = dtmp + CICoeffOnHotDet * CICoeffOnHotDet ; 
          
        }

        dtmp = sqrt( dtmp ) ;
        
        
        dtmp2 = 0.000; 
        
        for( iCTHotDet = 0 ; iCTHotDet < nCTHotDet ; iCTHotDet ++ )
        {
          currentDetInitialInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( ctHotDetList + iCTHotDet ) ) * 4 + 1 ) ) - 1 ) ) + 1 ;
      
          currentDetFinalInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( ctHotDetList + iCTHotDet ) ) * 4 + 2 ) ) - 1 ) ) + 1 ;
          
          CICoeffOnHotDet = *( ciCoefficients + ( *( ctHotDetList + iCTHotDet ) ) * nnout + icistate ) ;
          
          dtmp2 = dtmp2 + CICoeffOnHotDet * CICoeffOnHotDet ; 
          
        
        }
        
        dtmp2 = sqrt( dtmp2 ) ;

        
        printf("\nPreviously , for # %d CI , dominate CI-Coefficient is on # %d Det , value is % 10.6f ...\n" , icistate + 1 , *( dominDetList + icistate ) , preDominCICoeff ) ;
        
        
        
        if( dtmp >= preDominCICoeff && dtmp > dtmp2 )
        {
          printf("\n\n-----> Found It ! In snapshot , %d CI is 1st LE <-----\n Looks like we found out that # %d CI has dominate combined LE nature ... with scaled coefficient value % 10.6f ... \n" , icistate + 1 , icistate + 1 , dtmp ) ;
          
          printf("\n@ this excited state , LE has % 10.6f combined nature and CT has % 10.6f combined nature ...\n" , dtmp , dtmp2 ) ;
          
          leSnapCINumber = icistate ;
          
          *( assignmentCIStates + icistate ) = YES ;
          
          LEFactor = dtmp ;
          
          exSnapLE = 3 ;
          
          itmp = 1 ;
          
          for( iLEHotDet = 0 ; iLEHotDet < nLEHotDet ; iLEHotDet ++ )
          {
            itmp = itmp * ( *( leHotDetList + iLEHotDet ) - currentDominateDet ) ;
          }
          
          if( itmp != 0 )
          {
            //canonicalIntersection = YES ;
            
            printf("\nHowever, it appears no single LE nature determinate stood out for this CI .  ...\n\n") ;
          }
          else if( itmp == 0 && preDominCICoeff <= dominateFactorThreshold )
          {
            //canonicalIntersection = YES ;
            
            printf("\nHowever, it appears the single LE nature determinate stood out ( Det # %d ) for this CI has only % lf component. ...\n\n" , currentDominateDet + 1 , preDominCICoeff ) ;
          }
          else
          {
            //canonicalIntersection = canonicalIntersection ;  // Nothing happened here ...
          }          

          if( fabs( dtmp * dtmp - dtmp2 * dtmp2 ) <= coefficientVicinity )
          {
            printf("\n! ---> BESIDES ... It looks like LE and CT are so close in terms of coefficient @ this CI state ... <--- ! \n") ;
            
            canonicalIntersection ++ ;
          }
          
          
          break ;
        
        }
      
      }       
            
    }
    
    
    
    printf("\n----------> Searching For The 1st Local-Excitation CI State DONE...\n\n") ;
  
  
  
  
  
    //----------------------> Another ... LE Eigenstate Searching ... <------------------//

    if( degenerateCase == YES )
    {
      exSnapLE = 0 ;
      
      /*
      aliasLERefInitialOrbital = *( leRefInitialOrbitalEqv + 1 ) ;
      
      aliasLERefFinalOrbital = *( leRefFinalOrbitalEqv + 1 ) ;
      
      aliasCTRefInitialOrbital = *( ctRefInitialOrbitalEqv + 0 ) ;
      
      aliasCTRefFinalOrbital = *( ctRefFinalOrbitalEqv + 0 ) ;
      */

      printf("\n----------> BEGIN Searching For The 2nd Local-Excitation CI State ...\n\n") ;
      
      /* 1st Round : Simple Searching     

      for( icistate = 1 ; icistate < nnout ; icistate ++ )
      {
        if( *( assignmentCIStates + icistate ) == YES )
        {
          printf("\n# %d CI State has already been assigned ... So it cannot be LE state again ...\n" , icistate + 1 ) ;
          
          continue ;
        }
        
        
        currentSnapInitialInRef = *( snapMOMapping + ( int ) ( *( snapCIList + icistate * 4 + 1 ) - 1 ) ) + 1 ; // ATTENTION : here currentSnapInitialInRef is in human-label
        
        currentSnapFinalInRef = *( snapMOMapping + ( int ) ( *( snapCIList + icistate * 4 + 2 ) - 1 ) ) + 1 ;
        
        printf("\n---> Now it is the # %d CI Excitation --- In SnapBasis , %d -> %d --- same as In RefBasis %d -> %d <---\n" , icistate + 1 , ( int ) *( snapCIList + icistate * 4 + 1 ) , ( int ) *( snapCIList + icistate * 4 + 2 ) , currentSnapInitialInRef , currentSnapFinalInRef ) ;
        
        if( currentSnapInitialInRef == aliasLERefInitialOrbital && currentSnapFinalInRef == aliasLERefFinalOrbital )
        {
          leSnapInitialOrbital = currentSnapInitialInRef ;
          
          leSnapFinalOrbital = currentSnapFinalInRef ;
          
          leSnapCIAnotherNumber = icistate ;
          
          *( assignmentCIStates + icistate ) = YES ;
          
          LEFactor = fabs( *( ciCoefficients + ( *( dominDetList + icistate ) - 1 ) * nnout + icistate ) ) ;  
          
          printf("\n-----> Found It ! In snapshot , %d CI is LE <-----\n" , leSnapCIAnotherNumber + 1 ) ;
          
          exSnapLE = 3 ;

          break ;
        }
        else
        {
          printf("\nEm...Not this one ... I mean %d CI is not LE ...\n" , icistate + 1 ) ;
          
          continue ;
        }
      
      }
      */
      
      
      /* 2nd Round : Scattered Multiple LE initial or final orbitals ...  
      
      itmp2 = 0 ; dtmp = 0.000 ; dtmp2 = 0.000 ;
      
      preDominCICoeff = 0.000 ;
      
      CICoeffOnHotDet = 0.000 ;
      
      
      if( exSnapLE == 0 ) // && lenRevLERefFinalMapping > 1 ) // 2nd Round : If multiple orbital has leRefFinalOrbital nature , add them together ; 
      { 
        printf("\nLooks like we did not find the LE state from 1st Round searching ... Let's try 2nd Round ...\n") ; 
        
        for( icistate = 1 ; icistate < nnout ; icistate ++ )
        {
          
          if( *( assignmentCIStates + icistate ) == YES )
          {
            printf("\n# %d CI State has already been assigned ... So it cannot be LE state again ...\n" , icistate + 1 ) ;
            
            continue ;
          }
        
          preDominCICoeff = *( ciCoefficients + ( *( dominDetList + icistate ) - 1 ) * nnout + icistate ) ;    
          
          preDominCICoeff = fabs( preDominCICoeff ) ;
          
          dtmp = 0.000 ; 

          for( iLEHotDet = 0 ; iLEHotDet < nLEHotDet ; iLEHotDet ++ )
          {
            currentDetInitialInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( leHotDetList + iLEHotDet ) ) * 4 + 1 ) ) - 1 ) ) + 1 ;
        
            currentDetFinalInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( leHotDetList + iLEHotDet ) ) * 4 + 2 ) ) - 1 ) ) + 1 ;
            
            if( currentDetInitialInRef == aliasLERefInitialOrbital && currentDetFinalInRef == aliasLERefFinalOrbital )
            {
              CICoeffOnHotDet = *( ciCoefficients + ( *( leHotDetList + iLEHotDet ) ) * nnout + icistate ) ;
              
              dtmp = dtmp + CICoeffOnHotDet * CICoeffOnHotDet ; 
            }
          
          }

          dtmp = sqrt( dtmp ) ;
          
          
          dtmp2 = 0.000; 
          
          for( iCTHotDet = 0 ; iCTHotDet < nCTHotDet ; iCTHotDet ++ )
          {
            currentDetInitialInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( ctHotDetList + iCTHotDet ) ) * 4 + 1 ) ) - 1 ) ) + 1 ;
        
            currentDetFinalInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( ctHotDetList + iCTHotDet ) ) * 4 + 2 ) ) - 1 ) ) + 1 ;
            
            if( currentDetInitialInRef == aliasCTRefInitialOrbital && currentDetFinalInRef == aliasCTRefFinalOrbital )
            {
              CICoeffOnHotDet = *( ciCoefficients + ( *( ctHotDetList + iCTHotDet ) ) * nnout + icistate ) ;
              
              dtmp2 = dtmp2 + CICoeffOnHotDet * CICoeffOnHotDet ; 
            }
          
          }
          
          dtmp2 = sqrt( dtmp2 ) ;

          
          printf("\nPreviously , for # %d CI , dominate CI-Coefficient is on # %d Det , value is % 10.6f ...\n" , icistate + 1 , *( dominDetList + icistate ) , preDominCICoeff ) ;
          
          
          
          if( dtmp > preDominCICoeff && dtmp > dtmp2 )
          {
            printf("\n\n-----> Found It ! In snapshot , %d CI is LE <-----\n Looks like we found out that # %d CI has dominate combined LE nature ... with scaled coefficient value % 10.6f ... \n" , icistate + 1 , icistate + 1 , dtmp ) ;
            
            printf("\n@ this excited state , LE has % 10.6f combined nature and CT has % 10.6f combined nature ...\n" , dtmp , dtmp2 ) ;
            
            leSnapCIAnotherNumber = icistate ;
            
            *( assignmentCIStates + icistate ) = YES ;
            
            LEFactor = dtmp ;
            
            exSnapLE = 3 ;
            
            break ;
          
          }
          else if( dtmp > preDominCICoeff && dtmp > dtmp2 && *( assignmentCIStates + icistate ) == YES )
          {
            printf("\nThis is AWKWARD ... You know , this CI state [ %d ] has already been assigned , so it cannot be LE any more ... \n" , icistate + 1 ) ;
          
            continue ;
          }
        
        
        } 
        
      }
      
      */
      
      /* 3rd Round : Add all LE nature together ...  */
      
      
      if( exSnapLE == 0 ) //&& lenRevLERefFinalMapping == 1 )
      {
        //printf("\nSo ... I guess either the 2nd Round died too , meaning we are really closed to canonical intersection ...\n") ;
        
        //printf("\nSo ... 3rd Round , Add all LE nature together !!!\n") ;
        
        for( icistate = 1 ; icistate < nnout ; icistate ++ )
        {
          
          if( *( assignmentCIStates + icistate ) == YES )
          {
            printf("\n# %d CI State has already been assigned ... So it cannot be LE state again ...\n" , icistate + 1 ) ;
            
            continue ;
          }

          currentDominateDet = *( dominDetList + icistate ) - 1 ;
          
          preDominCICoeff = *( ciCoefficients + currentDominateDet * nnout + icistate ) ;     
                  
          preDominCICoeff = fabs( preDominCICoeff ) ; 
        
          dtmp = 0.000 ; 

          for( iLEHotDet = 0 ; iLEHotDet < nLEHotDet ; iLEHotDet ++ )
          {
            currentDetInitialInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( leHotDetList + iLEHotDet ) ) * 4 + 1 ) ) - 1 ) ) + 1 ;
        
            currentDetFinalInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( leHotDetList + iLEHotDet ) ) * 4 + 2 ) ) - 1 ) ) + 1 ;
            
            CICoeffOnHotDet = *( ciCoefficients + ( *( leHotDetList + iLEHotDet ) ) * nnout + icistate ) ;
            
            dtmp = dtmp + CICoeffOnHotDet * CICoeffOnHotDet ; 
            
          }

          dtmp = sqrt( dtmp ) ;
          
          
          dtmp2 = 0.000; 
          
          for( iCTHotDet = 0 ; iCTHotDet < nCTHotDet ; iCTHotDet ++ )
          {
            currentDetInitialInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( ctHotDetList + iCTHotDet ) ) * 4 + 1 ) ) - 1 ) ) + 1 ;
        
            currentDetFinalInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( ctHotDetList + iCTHotDet ) ) * 4 + 2 ) ) - 1 ) ) + 1 ;
            
            CICoeffOnHotDet = *( ciCoefficients + ( *( ctHotDetList + iCTHotDet ) ) * nnout + icistate ) ;
            
            dtmp2 = dtmp2 + CICoeffOnHotDet * CICoeffOnHotDet ; 
            
          
          }
          
          dtmp2 = sqrt( dtmp2 ) ;

          
          printf("\nPreviously , for # %d CI , dominate CI-Coefficient is on # %d Det , value is % 10.6f ...\n" , icistate + 1 , *( dominDetList + icistate ) , preDominCICoeff ) ;
          
          
          
          if( dtmp > preDominCICoeff && dtmp > dtmp2 )
          {
            printf("\n\n-----> Found It ! In snapshot , %d CI is 2nd LE <-----\n Looks like we found out that # %d CI has dominate combined LE nature ... with scaled coefficient value % 10.6f ... \n" , icistate + 1 , icistate + 1 , dtmp ) ;
            
            printf("\n@ this excited state , LE has % 10.6f combined nature and CT has % 10.6f combined nature ...\n" , dtmp , dtmp2 ) ;
                        
            leSnapCIAnotherNumber = icistate ;
            
            *( assignmentCIStates + icistate ) = YES ;
            
            LEFactor = dtmp ;
            
            exSnapLE = 3 ;
            
            //energyDiffOneAnother = fabs( *( snapCIList + 4 * leSnapCINumber + 3 ) ) - ( *( snapCIList + 4 * leSnapCIAnotherNumber + 3 ) ) ;
            
            itmp = 1 ;
          
            for( iLEHotDet = 0 ; iLEHotDet < nLEHotDet ; iLEHotDet ++ )
            {
              itmp = itmp * ( *( leHotDetList + iLEHotDet ) - currentDominateDet ) ;
            }
          
            if( itmp != 0 )
            {
              //canonicalIntersection = YES ;
              
              printf("\nHowever, it appears no single LE nature determinate stood out for this CI .  ...\n\n") ;
            }
            else if( itmp == 0 && preDominCICoeff <= dominateFactorThreshold ) //&& energyDiffOneAnother <= energyVicinity )
            {
              //canonicalIntersection = YES ;
              
              printf("\nHowever, it appears the single LE nature determinate stood out ( Det # %d ) for this CI has only % lf component. ...\n\n" , currentDominateDet + 1 , preDominCICoeff ) ;
            }
            else
            {
              //canonicalIntersection = canonicalIntersection ;  // Nothing happened here ...
            }          
            
            
            if( fabs( dtmp * dtmp - dtmp2 * dtmp2 ) <= coefficientVicinity )
            {
              printf("\n! ---> BESIDES ... It looks like LE and CT are so close in terms of coefficient @ this CI state ... <--- ! \n") ;
              
              canonicalIntersection ++ ;
            }
            
            break ;
          
          }
        
        }       
              
      }
    
      printf("\n----------> Searching For The 2nd Local-Excitation CI State DONE...\n\n") ;
    
    
    } // if( degenerateCase == YES )
  
  
  
  
  
  
  
  
  
  
  
  
  
    
    //----------------------> CT Eigenstate Searching ... <------------------//
    
    exSnapCT = 0 ;
    
    /*
    aliasLERefInitialOrbital = *( leRefInitialOrbitalEqv + 0 ) ;
    
    aliasLERefFinalOrbital = *( leRefFinalOrbitalEqv + 0 ) ;
    
    aliasCTRefInitialOrbital = *( ctRefInitialOrbitalEqv + 0 ) ;
    
    aliasCTRefFinalOrbital = *( ctRefFinalOrbitalEqv + 0 ) ;
    */

    printf("\n----------> BEGIN Searching For The 1st Charge-Transfer CI State ...\n\n") ;

    /* 1st Round : Simple Searching   
    
    for( icistate = 1 ; icistate < nnout ; icistate ++ )
    {
      
      if( *( assignmentCIStates + icistate ) == YES )
      {
        printf("\n# %d CI State has already been assigned ... So it cannot be CT state again ...\n" , icistate + 1 ) ;
        
        continue ;
      }
      
      currentSnapInitialInRef = *( snapMOMapping + ( int ) ( *( snapCIList + icistate * 4 + 1 ) - 1 ) ) + 1 ;
      
      currentSnapFinalInRef = *( snapMOMapping + ( int ) ( *( snapCIList + icistate * 4 + 2 ) - 1 ) ) + 1 ;
     
      printf("\n---> Now it is the # %d CI Excitation --- In SnapBasis , %d -> %d --- same as In RefBasis %d -> %d <---\n" , icistate + 1 , ( int ) *( snapCIList + icistate * 4 + 1 ) , ( int ) *( snapCIList + icistate * 4 + 2 ) , currentSnapInitialInRef , currentSnapFinalInRef ) ;
 
      if( currentSnapInitialInRef == aliasCTRefInitialOrbital && currentSnapFinalInRef == aliasCTRefFinalOrbital )
      {
        ctSnapInitialOrbital = currentSnapInitialInRef ;
        
        ctSnapFinalOrbital = currentSnapFinalInRef ;
        
        ctSnapCINumber = icistate ;
        
        *( assignmentCIStates + icistate ) = YES ;
        
        CTFactor = fabs( *( ciCoefficients + ( *( dominDetList + icistate ) - 1 ) * nnout + icistate ) ) ;   
        
        printf("\n-----> Found It ! In snapshot , %d CI is CT <-----\n" , ctSnapCINumber + 1 ) ;
        
        exSnapCT = 9 ;
        
        break ;
      }
      else
      {
        printf("\nEm...Not this one ... I mean %d CI is not CT ...\n" , icistate + 1 ) ;
        
        continue ;
      }
    
    }
    */
    
    /* 2nd Round : Scattered Multiple CT initial or final orbitals ...  
    
    itmp2 = 0 ; dtmp = 0.000 ; dtmp2 = 0.000 ;
    
    preDominCICoeff = 0.000 ;
    
    CICoeffOnHotDet = 0.000 ;
    
    
    if( exSnapCT == 0 ) // && lenRevLERefFinalMapping > 1 ) // 2nd Round : If multiple orbital has leRefFinalOrbital nature , add them together ; 
    { 
      printf("\nLooks like we did not find the CT state from 1st Round searching ... Let's try 2nd Round ...\n") ; 
      
      for( icistate = 1 ; icistate < nnout ; icistate ++ )
      {
        
        if( *( assignmentCIStates + icistate ) == YES )
        {
          printf("\n# %d CI State has already been assigned ... So it cannot be CT state again ...\n" , icistate + 1 ) ;
          
          continue ;
        }
        
        preDominCICoeff = *( ciCoefficients + ( *( dominDetList + icistate ) - 1 ) * nnout + icistate ) ;    
        
        preDominCICoeff = fabs( preDominCICoeff ) ;
        
        dtmp = 0.000 ; 

        for( iLEHotDet = 0 ; iLEHotDet < nLEHotDet ; iLEHotDet ++ )
        {
          currentDetInitialInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( leHotDetList + iLEHotDet ) ) * 4 + 1 ) ) - 1 ) ) + 1 ;
      
          currentDetFinalInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( leHotDetList + iLEHotDet ) ) * 4 + 2 ) ) - 1 ) ) + 1 ;
          
          if( currentDetInitialInRef == aliasLERefInitialOrbital && currentDetFinalInRef == aliasLERefInitialOrbital )
          {
            CICoeffOnHotDet = *( ciCoefficients + ( *( leHotDetList + iLEHotDet ) ) * nnout + icistate ) ;
            
            dtmp = dtmp + CICoeffOnHotDet * CICoeffOnHotDet ; 
          }
        
        }

        dtmp = sqrt( dtmp ) ;
        
        
        dtmp2 = 0.000; 
        
        for( iCTHotDet = 0 ; iCTHotDet < nCTHotDet ; iCTHotDet ++ )
        {
          currentDetInitialInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( ctHotDetList + iCTHotDet ) ) * 4 + 1 ) ) - 1 ) ) + 1 ;
      
          currentDetFinalInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( ctHotDetList + iCTHotDet ) ) * 4 + 2 ) ) - 1 ) ) + 1 ;
          
          if( currentDetInitialInRef == aliasCTRefInitialOrbital && currentDetFinalInRef == aliasCTRefFinalOrbital )
          {
            CICoeffOnHotDet = *( ciCoefficients + ( *( ctHotDetList + iCTHotDet ) ) * nnout + icistate ) ;
            
            dtmp2 = dtmp2 + CICoeffOnHotDet * CICoeffOnHotDet ; 
          }
        
        }
        
        dtmp2 = sqrt( dtmp2 ) ;

        
        printf("\nPreviously , for # %d CI , dominate CI-Coefficient is on # %d Det , value is % 10.6f ...\n" , icistate + 1 , *( dominDetList + icistate ) , preDominCICoeff ) ;
        
        
        
        if( dtmp2 > preDominCICoeff && dtmp2 > dtmp )
        {
          printf("\n\n-----> Found It ! In snapshot , %d CI is CT <-----\n Looks like we found out that # %d CI has dominate combined LE nature ... with scaled coefficient value % 10.6f ... \n" , icistate + 1 , icistate + 1 , dtmp ) ;
          
          printf("\n@ this excited state , LE has % 10.6f combined nature and CT has % 10.6f combined nature ...\n" , dtmp , dtmp2 ) ;
          
          ctSnapCINumber = icistate ;
          
          *( assignmentCIStates + icistate ) = YES ;
          
          CTFactor = dtmp2 ;
          
          exSnapCT = 9 ;
          
          break ;
        
        }
      
      } 
      
    }
    */
    
    
    /* 3rd Round : Add all CT nature together ...  */
    
    
    if( exSnapCT == 0 ) //&& lenRevLERefFinalMapping == 1 )
    {
      //printf("\nSo ... I guess either the 2nd Round died too , meaning we are really closed to canonical intersection ...\n") ;
      
      //printf("\nSo ... 3rd Round , Add all CT nature together !!!\n") ;
      
      for( icistate = 1 ; icistate < nnout ; icistate ++ )
      {
        if( *( assignmentCIStates + icistate ) == YES )
        {
          printf("\n# %d CI State has already been assigned ... So it cannot be CT state again ...\n" , icistate + 1 ) ;
          
          continue ;
        }
        
        currentDominateDet = *( dominDetList + icistate ) - 1 ;
        
        preDominCICoeff = *( ciCoefficients + currentDominateDet * nnout + icistate ) ;     

        preDominCICoeff = fabs( preDominCICoeff ) ; 
      
        dtmp = 0.000 ; 

        for( iLEHotDet = 0 ; iLEHotDet < nLEHotDet ; iLEHotDet ++ )
        {
          currentDetInitialInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( leHotDetList + iLEHotDet ) ) * 4 + 1 ) ) - 1 ) ) + 1 ;
      
          currentDetFinalInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( leHotDetList + iLEHotDet ) ) * 4 + 2 ) ) - 1 ) ) + 1 ;
          
          CICoeffOnHotDet = *( ciCoefficients + ( *( leHotDetList + iLEHotDet ) ) * nnout + icistate ) ;
          
          dtmp = dtmp + CICoeffOnHotDet * CICoeffOnHotDet ; 
          
        }

        dtmp = sqrt( dtmp ) ;
        
        
        dtmp2 = 0.000; 
        
        for( iCTHotDet = 0 ; iCTHotDet < nCTHotDet ; iCTHotDet ++ )
        {
          currentDetInitialInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( ctHotDetList + iCTHotDet ) ) * 4 + 1 ) ) - 1 ) ) + 1 ;
      
          currentDetFinalInRef = *( snapMOMapping +  ( ( int ) ( *( snapDetList + ( *( ctHotDetList + iCTHotDet ) ) * 4 + 2 ) ) - 1 ) ) + 1 ;
          
          CICoeffOnHotDet = *( ciCoefficients + ( *( ctHotDetList + iCTHotDet ) ) * nnout + icistate ) ;
          
          dtmp2 = dtmp2 + CICoeffOnHotDet * CICoeffOnHotDet ; 
          
        
        }
        
        dtmp2 = sqrt( dtmp2 ) ;

        
        printf("\nPreviously , for # %d CI , dominate CI-Coefficient is on # %d Det , value is % 10.6f ...\n" , icistate + 1 , *( dominDetList + icistate ) , preDominCICoeff ) ;
        
        
        
        if( dtmp2 > preDominCICoeff && dtmp2 > dtmp )
        {
          printf("\n\n-----> Found It ! In snapshot , %d CI is CT <-----\n Looks like we found out that # %d CI has dominate combined LE nature ... with scaled coefficient value % 10.6f ... \n" , icistate + 1 , icistate + 1 , dtmp ) ;
          
          printf("\n@ this excited state , LE has % 10.6f combined nature and CT has % 10.6f combined nature ...\n" , dtmp , dtmp2 ) ;
          
          ctSnapCINumber = icistate ;
          
          *( assignmentCIStates + icistate ) = YES ;
          
          CTFactor = dtmp2 ;
          
          exSnapCT = 9 ;
          
          itmp = 1 ;
          
          for( iCTHotDet = 0 ; iCTHotDet < nCTHotDet ; iCTHotDet ++ )
          {
            itmp = itmp * ( *( ctHotDetList + iCTHotDet ) - currentDominateDet ) ;
          }
          
          if( itmp != 0 )
          {
            //canonicalIntersection = YES ;
            
            printf("\nHowever, it appears no single CT nature determinate stood out for this CI .  ...\n\n") ;
          }
          else if( itmp == 0 && preDominCICoeff <= dominateFactorThreshold )
          {
            //canonicalIntersection = YES ;
            
            printf("\nHowever, it appears the single CT nature determinate stood out ( Det # %d ) for this CI has only % lf component ...\n\n" , currentDominateDet + 1 , preDominCICoeff ) ;
          }
          else
          {
            canonicalIntersection = canonicalIntersection ;  // Nothing happened here ...
          }          

          if( fabs( dtmp - dtmp2 ) <= coefficientVicinity )
          {
            printf("\n! ---> BESIDES ... It looks like LE and CT are so close in terms of coefficient @ this CI state ... <--- ! \n") ;
            
            canonicalIntersection ++ ;
          }

          break ;
        
        }
      
      }       
            
    }
    
    printf("\n----------> Searching For The 1st Charge-Transfer CI DONE ...\n\n") ;
    
    
    //----------------------> Checking Again for the Canonical Intersection <------------------//    
    
    energyDiffOneOne = fabs( *( snapCIList + 4 * leSnapCINumber + 3 ) ) - ( *( snapCIList + 4 * ctSnapCINumber + 3 ) ) ;
    

    if( degenerateCase == YES )
    {
      energyDiffOneAnother = fabs( *( snapCIList + 4 * leSnapCIAnotherNumber + 3 ) ) - ( *( snapCIList + 4 * ctSnapCINumber + 3 ) ) ;
      
      if( energyDiffOneOne <= energyVicinity && energyDiffOneAnother <= energyVicinity )
      {
        printf("\n! ---> BESIDES ... It looks like LE and CT are so close in terms of energy @ this CI state ... <--- ! \n") ;
        
        canonicalIntersection ++ ;
      }
    
    }
    else
    {
      if( energyDiffOneOne <= energyVicinity )
      {
        printf("\n! ---> BESIDES ... It looks like LE and CT are so close in terms of energy @ this CI state ... <--- ! \n") ;
        
        canonicalIntersection ++ ;
      }
    
    }
    
    
    //----------------------> LE and CT Eigenstate Searching Completed <------------------//
    
    /*
    if( leSnapCINumber == ctSnapCINumber )
    {
      printf("\n=> CodeBlack <=\n") ;
   
      exit( 195 ) ;
    }
    */
    
    
    
    // -----> Do the actual GMH ...
    
    double * MU11 , * MU22 , * MU12 , * V0 ;
    
    double mu11 , mu22 , mu12 ;
    
    double * tmpMU ;
    
    double projectMU11 , projectMU22 , projectMU12 , v0 ;
    
    double dE ;
    
    
    double HDA_GSLE = 0.000 ; double HDA_GSCT = 0.000 ; double HDA_LECT = 0.000 ;
    
    double HDA_GSLE_LINEAR = 0.000 ; double HDA_GSCT_LINEAR = 0.000 ; double HDA_LECT_LINEAR = 0.000 ;
    
    
    double HDA_GSLE_ANOTHER = 0.000 ; double HDA_LECT_ANOTHER = 0.000 ; 
    
    double HDA_GSLE_LINEAR_ANOTHER = 0.000 ; double HDA_LECT_LINEAR_ANOTHER = 0.000 ; 
    


    int flagGSCT = 0 ; int flagLECT = 0 ;


    MU11 = calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , MU11 ) ;
    
    MU22 = calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , MU22 ) ;
    
    MU12 = calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , MU12 ) ;
    
    tmpMU = calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , tmpMU ) ;
    
    V0 = calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , V0 ) ;
    



    if( exSnapCT != 0 )
    {
      flagGSCT = 1 ;     

      dzeros( 3 , 1 , MU11 ) ;
    
      dzeros( 3 , 1 , MU22 ) ;
    
      dzeros( 3 , 1 , MU12 ) ;
    
      dzeros( 3 , 1 , tmpMU ) ;
    
      dzeros( 3 , 1 , V0 ) ;
    
      // Ground State Dipole
      
      *( MU11 + 0 ) = *( dipole + 0 * ntrans + 0 * ncistate + 0 ) ;
      
      *( MU11 + 1 ) = *( dipole + 1 * ntrans + 0 * ncistate + 0 ) ;
      
      *( MU11 + 2 ) = *( dipole + 2 * ntrans + 0 * ncistate + 0 ) ;
      
      printf("\nMU11 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( MU11 + 0 ) ,  *( MU11 + 1 )  , *( MU11 + 2 ) ) ;
      
      mu11 = sqrt( (*( MU11 + 0 )) * (*( MU11 + 0 )) +  (*( MU11 + 1 )) * (*( MU11 + 1 )) + (*( MU11 + 2 )) * (*( MU11 + 2 ))  ) ;
      
      // CT State Dipole
      
      *( MU22 + 0 ) = *( dipole + 0 * ntrans + ctSnapCINumber * ncistate + ctSnapCINumber ) ;
      
      *( MU22 + 1 ) = *( dipole + 1 * ntrans + ctSnapCINumber * ncistate + ctSnapCINumber ) ;
      
      *( MU22 + 2 ) = *( dipole + 2 * ntrans + ctSnapCINumber * ncistate + ctSnapCINumber ) ;
      
      printf("\nMU22 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( MU22 + 0 ) ,  *( MU22 + 1 )  , *( MU22 + 2 ) ) ;
      
      mu22 = sqrt( (*( MU22 + 0 )) * (*( MU22 + 0 )) +  (*( MU22 + 1 )) * (*( MU22 + 1 )) + (*( MU22 + 2 )) * (*( MU22 + 2 ))  ) ;
      
      // Transition Dipole Between Two States ...
      
      *( MU12 + 0 ) = *( dipole + 0 * ntrans + 0 * ncistate + ctSnapCINumber ) ;
      
      *( MU12 + 1 ) = *( dipole + 1 * ntrans + 0 * ncistate + ctSnapCINumber ) ;
      
      *( MU12 + 2 ) = *( dipole + 2 * ntrans + 0 * ncistate + ctSnapCINumber ) ;
      
      printf("\nMU12 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( MU12 + 0 ) ,  *( MU12 + 1 )  , *( MU12 + 2 ) ) ;
      
      mu12 = sqrt( (*( MU12 + 0 )) * (*( MU12 + 0 )) +  (*( MU12 + 1 )) * (*( MU12 + 1 )) + (*( MU12 + 2 )) * (*( MU12 + 2 ))  ) ;
      
      *( V0 + 0 ) = *( MU11 + 0 ) - *( MU22 + 0 ) ;
      
      *( V0 + 1 ) = *( MU11 + 1 ) - *( MU22 + 1 ) ;
      
      *( V0 + 2 ) = *( MU11 + 2 ) - *( MU22 + 2 ) ;
      
      v0 = sqrt( (*( V0 + 0 )) * (*( V0 + 0 )) +  (*( V0 + 1 )) * (*( V0 + 1 )) + (*( V0 + 2 )) * (*( V0 + 2 ))  ) ;
      
      projectMU11 = ( ( *( MU11 + 0 ) ) * ( *( V0 + 0 ) ) + ( *( MU11 + 1 ) ) * ( *( V0 + 1 ) ) + ( *( MU11 + 2 ) ) * ( *( V0 + 2 ) ) ) / v0 ;
      
      projectMU22 = ( ( *( MU22 + 0 ) ) * ( *( V0 + 0 ) ) + ( *( MU22 + 1 ) ) * ( *( V0 + 1 ) ) + ( *( MU22 + 2 ) ) * ( *( V0 + 2 ) ) ) / v0 ;
      
      projectMU12 = ( ( *( MU12 + 0 ) ) * ( *( V0 + 0 ) ) + ( *( MU12 + 1 ) ) * ( *( V0 + 1 ) ) + ( *( MU12 + 2 ) ) * ( *( V0 + 2 ) ) ) / v0 ;
    
      
      dE = *( snapCIList + 4 * ctSnapCINumber + 3 ) ;

      dE = dE / ( EV2HARTREE ) ;
      
      printf("\ndE is % 10.8f...\n" , dE ) ;
      
      HDA_GSCT = fabs( projectMU12 ) * fabs( dE ) / sqrt( v0 * v0 + 4.00 * projectMU12 * projectMU12 ) ;
      
      HDA_GSCT_LINEAR = fabs( mu12 ) * fabs( dE ) / sqrt( v0 * v0 + 4.00 * mu12 * mu12 ) ;
      
      dtmp = fabs( mu12 ) * fabs( dE ) / sqrt( ( mu11 - mu22 ) * ( mu11 - mu22 ) + 4.00 * mu12 * mu12 ) ;
      
      //HDA_GSLE = HDA_GSLE * HDA_GSLE ;
      
      printf("\nHDA between GS and CT is % 10.8E eV...\n" , HDA_GSCT ) ;
      
      printf("\nHDA between GS and CT is % 10.8E eV ( linear ) ...\n" , HDA_GSCT_LINEAR ) ;
      
      printf("\nHDA between GS and CT is % 10.8E eV ( crappy ) ...\n" , dtmp ) ;
      
      
    }  
    else
    {
      printf("\nUNFORTUNATELY ... we did not find All corresponding CT state ... Mission Aborting \n=> CodeRed <=\n\n") ;      
    }
    
   

    if( exSnapCT != 0 && exSnapLE != 0 )
    {
      printf("\n-----> HDA between LE and CT now ... <-----\n") ;
      
      flagLECT = 1 ;

      dzeros( 3 , 1 , MU11 ) ;
    
      dzeros( 3 , 1 , MU22 ) ;
    
      dzeros( 3 , 1 , MU12 ) ;
    
      dzeros( 3 , 1 , tmpMU ) ;
    
      dzeros( 3 , 1 , V0 ) ;
    
      
      
      *( MU11 + 0 ) = *( dipole + 0 * ntrans + leSnapCINumber * ncistate + leSnapCINumber ) ;
      
      *( MU11 + 1 ) = *( dipole + 1 * ntrans + leSnapCINumber * ncistate + leSnapCINumber ) ;
      
      *( MU11 + 2 ) = *( dipole + 2 * ntrans + leSnapCINumber * ncistate + leSnapCINumber ) ;
      
      printf("\nMU11 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( MU11 + 0 ) ,  *( MU11 + 1 )  , *( MU11 + 2 ) ) ;
      
      mu11 = sqrt( (*( MU11 + 0 )) * (*( MU11 + 0 )) +  (*( MU11 + 1 )) * (*( MU11 + 1 )) + (*( MU11 + 2 )) * (*( MU11 + 2 ))  ) ;
      
      *( MU22 + 0 ) = *( dipole + 0 * ntrans + ctSnapCINumber * ncistate + ctSnapCINumber ) ;
      
      *( MU22 + 1 ) = *( dipole + 1 * ntrans + ctSnapCINumber * ncistate + ctSnapCINumber ) ;
      
      *( MU22 + 2 ) = *( dipole + 2 * ntrans + ctSnapCINumber * ncistate + ctSnapCINumber ) ;
      
      printf("\nMU22 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( MU22 + 0 ) ,  *( MU22 + 1 )  , *( MU22 + 2 ) ) ;
      
      mu22 = sqrt( (*( MU22 + 0 )) * (*( MU22 + 0 )) +  (*( MU22 + 1 )) * (*( MU22 + 1 )) + (*( MU22 + 2 )) * (*( MU22 + 2 ))  ) ;
      
      *( MU12 + 0 ) = *( dipole + 0 * ntrans + leSnapCINumber * ncistate + ctSnapCINumber ) ;
      
      *( MU12 + 1 ) = *( dipole + 1 * ntrans + leSnapCINumber * ncistate + ctSnapCINumber ) ;
      
      *( MU12 + 2 ) = *( dipole + 2 * ntrans + leSnapCINumber * ncistate + ctSnapCINumber ) ;
      
      printf("\nMU12 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( MU12 + 0 ) ,  *( MU12 + 1 )  , *( MU12 + 2 ) ) ;
      
      mu12 = sqrt( (*( MU12 + 0 )) * (*( MU12 + 0 )) +  (*( MU12 + 1 )) * (*( MU12 + 1 )) + (*( MU12 + 2 )) * (*( MU12 + 2 ))  ) ;
      
      *( V0 + 0 ) = *( MU11 + 0 ) - *( MU22 + 0 ) ;
      
      *( V0 + 1 ) = *( MU11 + 1 ) - *( MU22 + 1 ) ;
      
      *( V0 + 2 ) = *( MU11 + 2 ) - *( MU22 + 2 ) ;
      
      v0 = sqrt( (*( V0 + 0 )) * (*( V0 + 0 )) +  (*( V0 + 1 )) * (*( V0 + 1 )) + (*( V0 + 2 )) * (*( V0 + 2 ))  ) ;
      
      projectMU11 = ( ( *( MU11 + 0 ) ) * ( *( V0 + 0 ) ) + ( *( MU11 + 1 ) ) * ( *( V0 + 1 ) ) + ( *( MU11 + 2 ) ) * ( *( V0 + 2 ) ) ) / v0 ;
      
      projectMU22 = ( ( *( MU22 + 0 ) ) * ( *( V0 + 0 ) ) + ( *( MU22 + 1 ) ) * ( *( V0 + 1 ) ) + ( *( MU22 + 2 ) ) * ( *( V0 + 2 ) ) ) / v0 ;
      
      projectMU12 = ( ( *( MU12 + 0 ) ) * ( *( V0 + 0 ) ) + ( *( MU12 + 1 ) ) * ( *( V0 + 1 ) ) + ( *( MU12 + 2 ) ) * ( *( V0 + 2 ) ) ) / v0 ;
    
      
      dE = ( *( snapCIList + 4 * leSnapCINumber + 3 ) ) - ( *( snapCIList + 4 * ctSnapCINumber + 3 ) ) ;

      dE = dE / ( EV2HARTREE ) ;
      
      printf("\ndE is % 10.8f...\n" , dE ) ;
      
      HDA_LECT = fabs( projectMU12 ) * fabs( dE ) / sqrt( v0 * v0 + 4.00 * projectMU12 * projectMU12 ) ;
      
      HDA_LECT_LINEAR = fabs( mu12 ) * fabs( dE ) / sqrt( v0 * v0 + 4.00 * mu12 * mu12 ) ;
      
      dtmp = fabs( mu12 ) * fabs( dE ) / sqrt( ( mu11 - mu22 ) * ( mu11 - mu22 ) + 4.00 * mu12 * mu12 ) ;
      
      //HDA_GSLE = HDA_GSLE * HDA_GSLE ;
      
      printf("\nHDA between LE and CT is % 10.8E eV ...\n" , HDA_LECT ) ;
      
      printf("\nHDA between LE and CT is % 10.8E eV ( linear ) ...\n" , HDA_LECT_LINEAR ) ;
      
      printf("\nHDA between LE and CT is % 10.8E eV ( crappy ) ...\n" , dtmp ) ;
      
      
    }  
    else
    {
      if( exSnapCT != 0 && exSnapLE == 0 )
      {
        printf("\nUNFORTUNATELY ... we did not find corresponding LE state ... Mission Aborting \n=> CodeYellow <=\n\n") ;
      }
      else if( exSnapCT == 0 && exSnapLE == 0 )
      {
        printf("\nUNFORTUNATELY ... we did not find corresponding LE and CT state ... Mission Aborting \n=> CodeOrange <=\n\n") ;
      }
    }


    if( exSnapCT != 0 && exSnapLE != 0 && degenerateCase == YES )
    {
      
      // ANOTHER
      
      printf("\n-----> HDA between LE ( ANOTHER ) and CT now ... <-----\n") ;
      
      flagLECT = 1 ;

      dzeros( 3 , 1 , MU11 ) ;
    
      dzeros( 3 , 1 , MU22 ) ;
    
      dzeros( 3 , 1 , MU12 ) ;
    
      dzeros( 3 , 1 , tmpMU ) ;
    
      dzeros( 3 , 1 , V0 ) ;
    
      
      
      *( MU11 + 0 ) = *( dipole + 0 * ntrans + leSnapCIAnotherNumber * ncistate + leSnapCIAnotherNumber ) ;
      
      *( MU11 + 1 ) = *( dipole + 1 * ntrans + leSnapCIAnotherNumber * ncistate + leSnapCIAnotherNumber ) ;
      
      *( MU11 + 2 ) = *( dipole + 2 * ntrans + leSnapCIAnotherNumber * ncistate + leSnapCIAnotherNumber ) ;
      
      printf("\nMU11 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( MU11 + 0 ) ,  *( MU11 + 1 )  , *( MU11 + 2 ) ) ;
      
      mu11 = sqrt( (*( MU11 + 0 )) * (*( MU11 + 0 )) +  (*( MU11 + 1 )) * (*( MU11 + 1 )) + (*( MU11 + 2 )) * (*( MU11 + 2 ))  ) ;
      
      *( MU22 + 0 ) = *( dipole + 0 * ntrans + ctSnapCINumber * ncistate + ctSnapCINumber ) ;
      
      *( MU22 + 1 ) = *( dipole + 1 * ntrans + ctSnapCINumber * ncistate + ctSnapCINumber ) ;
      
      *( MU22 + 2 ) = *( dipole + 2 * ntrans + ctSnapCINumber * ncistate + ctSnapCINumber ) ;
      
      printf("\nMU22 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( MU22 + 0 ) ,  *( MU22 + 1 )  , *( MU22 + 2 ) ) ;
      
      mu22 = sqrt( (*( MU22 + 0 )) * (*( MU22 + 0 )) +  (*( MU22 + 1 )) * (*( MU22 + 1 )) + (*( MU22 + 2 )) * (*( MU22 + 2 ))  ) ;
      
      *( MU12 + 0 ) = *( dipole + 0 * ntrans + leSnapCIAnotherNumber * ncistate + ctSnapCINumber ) ;
      
      *( MU12 + 1 ) = *( dipole + 1 * ntrans + leSnapCIAnotherNumber * ncistate + ctSnapCINumber ) ;
      
      *( MU12 + 2 ) = *( dipole + 2 * ntrans + leSnapCIAnotherNumber * ncistate + ctSnapCINumber ) ;
      
      printf("\nMU12 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( MU12 + 0 ) ,  *( MU12 + 1 )  , *( MU12 + 2 ) ) ;
      
      mu12 = sqrt( (*( MU12 + 0 )) * (*( MU12 + 0 )) +  (*( MU12 + 1 )) * (*( MU12 + 1 )) + (*( MU12 + 2 )) * (*( MU12 + 2 ))  ) ;
      
      *( V0 + 0 ) = *( MU11 + 0 ) - *( MU22 + 0 ) ;
      
      *( V0 + 1 ) = *( MU11 + 1 ) - *( MU22 + 1 ) ;
      
      *( V0 + 2 ) = *( MU11 + 2 ) - *( MU22 + 2 ) ;
      
      v0 = sqrt( (*( V0 + 0 )) * (*( V0 + 0 )) +  (*( V0 + 1 )) * (*( V0 + 1 )) + (*( V0 + 2 )) * (*( V0 + 2 ))  ) ;
      
      projectMU11 = ( ( *( MU11 + 0 ) ) * ( *( V0 + 0 ) ) + ( *( MU11 + 1 ) ) * ( *( V0 + 1 ) ) + ( *( MU11 + 2 ) ) * ( *( V0 + 2 ) ) ) / v0 ;
      
      projectMU22 = ( ( *( MU22 + 0 ) ) * ( *( V0 + 0 ) ) + ( *( MU22 + 1 ) ) * ( *( V0 + 1 ) ) + ( *( MU22 + 2 ) ) * ( *( V0 + 2 ) ) ) / v0 ;
      
      projectMU12 = ( ( *( MU12 + 0 ) ) * ( *( V0 + 0 ) ) + ( *( MU12 + 1 ) ) * ( *( V0 + 1 ) ) + ( *( MU12 + 2 ) ) * ( *( V0 + 2 ) ) ) / v0 ;
    
      
      dE = ( *( snapCIList + 4 * leSnapCIAnotherNumber + 3 ) ) - ( *( snapCIList + 4 * ctSnapCINumber + 3 ) ) ;

      dE = dE / ( EV2HARTREE ) ;
      
      printf("\ndE is % 10.8f...\n" , dE ) ;
      
      HDA_LECT_ANOTHER = fabs( projectMU12 ) * fabs( dE ) / sqrt( v0 * v0 + 4.00 * projectMU12 * projectMU12 ) ;
      
      HDA_LECT_LINEAR_ANOTHER = fabs( mu12 ) * fabs( dE ) / sqrt( v0 * v0 + 4.00 * mu12 * mu12 ) ;
      
      dtmp = fabs( mu12 ) * fabs( dE ) / sqrt( ( mu11 - mu22 ) * ( mu11 - mu22 ) + 4.00 * mu12 * mu12 ) ;
      
      //HDA_GSLE = HDA_GSLE * HDA_GSLE ;
      
      printf("\nHDA between LE ( ANOTHER ) and CT is % 10.8E eV ...\n" , HDA_LECT_ANOTHER ) ;
      
      printf("\nHDA between LE ( ANOTHER ) and CT is % 10.8E eV ( linear ) ...\n" , HDA_LECT_LINEAR_ANOTHER ) ;
      
      printf("\nHDA between LE ( ANOTHER ) and CT is % 10.8E eV ( crappy ) ...\n" , dtmp ) ;
      
      
    

    
    }  


   
  // ---> Final Output for future analysis ... 


  len_cndoLogFileName = strlen( cndoLogFileName ) ;

  if( flagGSCT == 0 )
  {
    printf("\nGMH Calculation FAILED For This Case Because No CT State Found ...\n") ;
    
    exit( 199 ) ;
  }
  else if( flagGSCT != 0 && flagLECT == 0 )
  {
    printf("\nGMH Calculation FAILED For This Case Because No LE State Found ...\n") ;
    
    exit( 197 ) ;
  
  } 
  else
  {
    printf("\nGMH Calculation SUCCESSFUL !\n") ;
    
    strncpy( outputGMHName , cndoLogFileName , len_cndoLogFileName - 4 ) ;

    *( outputGMHName + len_cndoLogFileName - 4 ) = '\0' ;
    
    printf( "\n%s\n" , outputGMHName ) ;

    strcat( outputGMHName , ".gmh" ) ;

    debug = fopen( outputGMHName , "wb+" ) ;
    
    if( canonicalIntersection >= 3 )
    {
      strcpy( tmpString , "CanonicalIntersection" ) ;
    }
    else if( canonicalIntersection <= 2 )
    {
      strcpy( tmpString , "validRegionForGMHCalc" ) ;
    }

    fprintf( debug , "\n%s\t% 10.8E\t% 10.8E\t% 10.8E\n" , tmpString , HDA_GSCT , HDA_LECT , HDA_LECT_ANOTHER ) ;

    fclose( debug ) ;


  }



  return( 0 ) ;

}

