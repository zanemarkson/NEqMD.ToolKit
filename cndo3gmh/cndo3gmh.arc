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






int main( int argc , char * argv[ ] )
{

  char ** pcmd = argv ; 
  
  int icmd ;
  
  char refMOFileName[ MAXCHARINLINE ] , refCIListName[MAXCHARINLINE] , cndoLogFileName[ MAXCHARINLINE ] ,  outputGMHName[ MAXCHARINLINE ] ;

  int len_cndoLogFileName ;
  
  FILE * prefMOFile , * prefCIListFile , * pcndoLogFile , * poutputGMHFile ;
  
  double * refMO ; 
  
  int * refCIList ; 
  
  double * cndoMO , * MOEnergy ;
  
  int refMOFileLength , refCIListFileLength , lmo ; // lmo is just for convenience = nbasis^2
  
  int noReferenceMode = NO ;
  
  int HOMO , LUMO ;
  
  int nbasis , natom , nelectron , ncistate ;
  
  int ibasis , iatom , icistate , icidet ;
  
  int iload = 0 ; int iline = 0 ; int ieqv = 0 ;
  
  int nword , iword ;
  
  int exout = 18  ;
  
  int debuggingMode = NO ;
  

  char buffer[ MAXCHARINLINE ] ;
  
  char cache[ MAXCHARINLINE ] ;

  
  int itmp , itmp2 ; double dtmp , dtmp2 ; char ctmp , tmp_char ; char stmp[ 150 ] , stmp2[ 150 ] ;
  
  double dtmpArray[ MAXCHARINLINE ] ;
  
  int info , signal , blank_signal ;
  
  FILE * debug ; 
  
  
  double done = 1.0000 ; double dzero = 0.0000 ;
  
  
  
  // ----------------------------------> Recording Command-Line Arguments ... <---------------------------------- //
  
  time_t current_time;

  time( &current_time );

  char now[ 300 ] ;

  strcpy( now , ctime( &current_time ) );

  int lennow = strlen( now ) ;

  *( now + lennow - 1 ) = ' ';

    
    
    printf("\n**********************************************************************\n");
      printf("* G_CNDO3GMH_D : Calculate GMH Type HDA for Snapshots in MD Calcs.   *\n");
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
  
  strcpy( refCIListName , "ref.ORB" ) ;
  
  strcpy( cndoLogFileName , "system.snap.log" ) ;
  
  //strcpy( cndoTRDFileName , "system.snap.trd" ) ;
  
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
	      
	      case 'o' : strcpy( outputGMHName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Output File name : %s ...\n" , outputGMHName ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 exout = 19 ;
	                 
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


	      case 'h' : printf("\nUsage:  %s [ -m 'input ref. MO Coefficient file name' ] [ -L 'Ref. CI list file name' ] [ -l 'Input cndo log file name' ] [ -o 'output Calculated HDA file name ' ] [ -g 'debugging mode'] \n\n" , * argv ); 
	                 
	                 //printf("\nOutput sequence : [ HDA(GS-CT) HDA(LE-CT) E(GS) E(CT) E(LE) DipoleX(GS) DipoleY(GS) DipoleZ(GS) DipoleX(CT) DipoleY(CT) DipoleZ(CT) DipoleX(LE) DipoleY(LE) DipoleZ(LE) \
	                         TransDipoleX(GS-CT) TransDipoleY(GS-CT) TransDipoleZ(GS-CT) TransDipoleX(LE-CT) TransDipoleY(LE-CT) TransDipoleZ(LE-CT)\n\n") ;

	                 printf("\n1) Output sequence : [ HDA(GS-CT) HDA(LE-CT) E(GS) E(CT) E(LE) DipoleX(GS) DipoleY(GS) DipoleZ(GS) DipoleX(CT) DipoleY(CT) DipoleZ(CT) DipoleX(LE) DipoleY(LE) DipoleZ(LE) TransDipole(GS-CT) TransDipole(LE-CT) ]\n\n") ;
	                 printf("\n2) \"none\" or \"None\" or \"NONE\" will specify no-reference-MO-mode to be invoked. Orbital label will be directly taken from CI file.\n\n") ;
	                 printf("\n3) For debugging mode , YES/Yes/yes will turn on the debugging print while NO/No/no will turn it off. Default is no debugging printing.\n\n") ;


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
  
  if( strcmp( refMOFileName , "none" ) == 0 || strcmp( refMOFileName , "NONE" ) == 0 || strcmp( refMOFileName , "None" ) == 0 )
  {
    noReferenceMode = YES ;
    
    printf("\nNo reference MO specified ... will use absolute orbital label in CI file ...\n\n") ;
  }
  else
  {
    if( ( prefMOFile = fopen( refMOFileName , "r" ) ) == NULL )
    {
      printf("\nUser specified reference MO coefficient NOT FOUND ...\n");
    
      exit( 63 ) ;
  
    }
  
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
  
  
  
  
  // =====> Loading Reference MO Coefficients ... <===== //
  
  refMO = ( double * ) calloc( refMOFileLength , sizeof( double ) ) ;
  
  if( noReferenceMode == NO )
  {
    refMOFileLength = flength( prefMOFile ) ;
  
    rewind( prefMOFile ) ;
  
    fload( prefMOFile , refMO );
  }
  
  
  // =====> Loading Reference CI List ... <===== //
  
  // const int leRefCINumber , leRefInitialOrbital , leRefFinalOrbital , ctRefCINumber , ctRefInitialOrbital , ctRefFinalOrbital ;
  
  int gsRefOrbital , leRefOrbital , ctRefOrbital ;
  
  int * degenerateStatesList[ 3 ] ;
  
  
  int  * gsRefOrbitalEqv , * leRefOrbitalEqv , * ctRefOrbitalEqv ;
  
  int  EXgsRefOrbitalEqv ,  EXleRefOrbitalEqv ,  EXctRefOrbitalEqv ;
  
  int  lenGSRefOrbitalEqv , lenLERefOrbitalEqv , lenCTRefOrbitalEqv ;
  

  EXgsRefOrbitalEqv = NO ;
  
  EXleRefOrbitalEqv = NO ;

  EXctRefOrbitalEqv = NO ;


  
  refCIListFileLength = fwc( prefCIListFile ) ;
  
  refCIList = ( int * ) calloc( 3 , sizeof( int ) ) ;
  
  for( itmp = 0 ; itmp < 3 ; itmp ++ ) *( refCIList + itmp ) = 0 ;
  
  
  rewind( prefCIListFile ) ;

  if( refCIListFileLength == 3 )
  {
    printf("\nOKay ... Regular Orbital List File ... \n") ;
  
    for( itmp = 0 ; itmp < 3 ; itmp ++ )
    {
      fscanf( prefCIListFile , "%d" , refCIList + itmp ) ;
    }
    
    gsRefOrbital = *( refCIList + 0 ) ;
  
    leRefOrbital = *( refCIList + 1 ) ;
  
    ctRefOrbital = *( refCIList + 2 ) ;
  
    
    for( iload = 0 ; iload < 4 ; iload ++ )
    {
      *( degenerateStatesList + iload ) = ( int * ) calloc( 1 , sizeof( int ) ) ;
    }
    
    gsRefOrbitalEqv = *( degenerateStatesList + 0 ) ; *( gsRefOrbitalEqv + 0 ) = gsRefOrbital ; lenGSRefOrbitalEqv = 1 ; // EXleRefInitialOrbitalEqv = YES ;

    leRefOrbitalEqv = *( degenerateStatesList + 1 ) ; *( leRefOrbitalEqv + 0 ) = leRefOrbital ; lenLERefOrbitalEqv = 1 ; 

    ctRefOrbitalEqv = *( degenerateStatesList + 2 ) ; *( ctRefOrbitalEqv + 0 ) = ctRefOrbital ; lenCTRefOrbitalEqv = 1 ;

  }
  else
  {
    printf("\nAll right ... degenerate case then ... \n") ;

    for( itmp = 0 ; itmp < 3 ; itmp ++ )
    {
      fscanf( prefCIListFile , "%d" , refCIList + itmp ) ;
    }
    
    
    gsRefOrbital = *( refCIList + 0 ) ;
  
    leRefOrbital = *( refCIList + 1 ) ;
  
    ctRefOrbital = *( refCIList + 2 ) ;
    
    
    printf("\nIn Reference , the 3 involved orbitals are # %d , # %d and # %d ... \n\n" , gsRefOrbital , leRefOrbital , ctRefOrbital ) ;
    
    
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
          
          
          if( itmp == gsRefOrbital && EXgsRefOrbitalEqv == NO )
          {
            printf("\nGS , equivalent to : \n") ;
            
            gsRefOrbitalEqv = *( degenerateStatesList + iload ) ;
            
            lenGSRefOrbitalEqv = nword - 1 ;
            
            *( gsRefOrbitalEqv + 0 ) = leRefOrbital ;
            
            for( iword = 2 ; iword < nword ; iword ++ ) 
            {
              strpickword( buffer , iword + 1 , cache ) ;
              
              *( gsRefOrbitalEqv + iword - 1 ) = atoi( cache ) ;
              
              printf( "\n# %d\n" , atoi( cache ) ) ;
            }
            
            EXgsRefOrbitalEqv = YES ;
            
            printf("\nFour indicators are : EXgsRefOrbitalEqv = %d , EXleRefOrbitalEqv = %d , EXctRefOrbitalEqv = %d " , EXgsRefOrbitalEqv , EXleRefOrbitalEqv , EXctRefOrbitalEqv ) ;
          
          }
          else if( itmp == leRefOrbital && EXleRefOrbitalEqv == NO )
          {
            printf("\nLE , equivalent to : \n") ;
            
            leRefOrbitalEqv = *( degenerateStatesList + iload ) ;
            
            lenLERefOrbitalEqv = nword - 1 ;
            
            *( leRefOrbitalEqv + 0 ) = leRefOrbital ;
            
            for( iword = 2 ; iword < nword ; iword ++ ) 
            {
              strpickword( buffer , iword + 1 , cache ) ;
              
              *( leRefOrbitalEqv + iword - 1 ) = atoi( cache ) ;
              
              printf( "\n# %d\n" , atoi( cache ) ) ;
            }
          
            EXleRefOrbitalEqv = YES ;
            
            printf("\nFour indicators are : EXgsRefOrbitalEqv = %d , EXleRefOrbitalEqv = %d , EXctRefOrbitalEqv = %d " , EXgsRefOrbitalEqv , EXleRefOrbitalEqv , EXctRefOrbitalEqv ) ;
          
          }
          else if( itmp == ctRefOrbital && EXctRefOrbitalEqv == NO )
          {
            printf("\nCT , equivalent to : \n") ;
            
            ctRefOrbitalEqv = *( degenerateStatesList + iload ) ;
            
            lenCTRefOrbitalEqv = nword - 1 ;
            
            *( ctRefOrbitalEqv + 0 ) = ctRefOrbital ;
            
            for( iword = 2 ; iword < nword ; iword ++ ) 
            {
              strpickword( buffer , iword + 1 , cache ) ;
              
              *( ctRefOrbitalEqv + iword - 1 ) = atoi( cache ) ;
              
              printf( "\n# %d\n" , atoi( cache ) ) ;
            }
              
            EXctRefOrbitalEqv = YES ;
          
            printf("\nFour indicators are : EXgsRefOrbitalEqv = %d , EXleRefOrbitalEqv = %d , EXctRefOrbitalEqv = %d " , EXgsRefOrbitalEqv , EXleRefOrbitalEqv , EXctRefOrbitalEqv ) ;
          
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
  
    if( iload != 3 )
    {
      printf("\nWrong Format in CI List File . All 3 Relevant Orbitals Must Appear In CI List Degneracy Specification Section ...\n") ;
      
      exit( 646 ) ;
    }
  
  }
  

  if( debuggingMode == YES )
  {
    debug = fopen( "eqv.deb" , "wb+" ) ;
  
    fprintf( debug , "\n\n===> GS <===\n\n") ;
  
    for( ieqv = 0 ; ieqv < lenGSRefOrbitalEqv ; ieqv ++ )
    {
      fprintf( debug , "%d\t" , *( gsRefOrbitalEqv + ieqv ) ) ;
  
    }

    fprintf( debug , "\n\n===> LE <===\n\n") ;
  
    for( ieqv = 0 ; ieqv < lenLERefOrbitalEqv ; ieqv ++ )
    {
      fprintf( debug , "%d\t" , *( leRefOrbitalEqv + ieqv ) ) ;
  
    }

    fprintf( debug , "\n\n===> CT <===\n\n") ;
  
    for( ieqv = 0 ; ieqv < lenCTRefOrbitalEqv ; ieqv ++ )
    {
      fprintf( debug , "%d\t" , *( ctRefOrbitalEqv + ieqv ) ) ;
  
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
    *( stmp2 + itmp - 5 ) = *( stmp + itmp ) ;
  
  }
  
  *( stmp + itmp ) = '\0' ;
  
  natom = atof( stmp2 );
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
  
  if( ( refMOFileLength != nbasis * nbasis ) && ( noReferenceMode == NO ) )
  {
    printf("\nSomething is wrong with the # of basis ... Log file say NBasis is %d but there are %d numbers in reference MO Coefficients...\n" , nbasis , refMOFileLength );
    
    printf("\nIt is also possible that this is related to the length of reference MO coefficients ...\n");
    
    exit( 145 );
  
  }
  else if( noReferenceMode == YES ) 
  {
    printf("\nThe no-reference-mode has been invoked. Temporarily setting refMOFileLength = nbasis^2 ...\n\n") ; 
    
    refMOFileLength = nbasis * nbasis ;
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


  
  // ---> Execute the transpose( refMO ) * cndoMO and making the mapping ... 
  
  double * orbitalOverlap = calloc( nbasis * nbasis , sizeof( double) ) ;
  
  dzeros( nbasis , nbasis , orbitalOverlap ) ;
  
  itmp2 = nbasis ;
  
  // each row of orbitalOverlap corresponds to MO# in refMO , column corresponds to MO# in cndoMO 
  
  if( noReferenceMode == NO )
  {
    dgemm_( "N" , "T" , &itmp2 , &itmp2 , &itmp2 , &done , refMO , &itmp2 , cndoMO , &itmp2 , &dzero , orbitalOverlap , &itmp2 ) ;
  
    dtranspose( itmp2 , orbitalOverlap , orbitalOverlap ) ;
  
  }
  
  

  /* 
  debug = fopen( "orbitalOverlap.deb" , "wb+" ) ;
  
  doutput( debug , nbasis, nbasis, orbitalOverlap ) ;
  
  fclose( debug ) ; 
  */ 
  
  
  
  int * snapMOMapping = calloc( nbasis , sizeof( int ) ) ;
  
    
  double * moEigvec = calloc( nbasis , sizeof( double ) ) ; dzeros( nbasis , 1 , moEigvec ) ;
  
  int irefOrbital , isnapOrbital ;
  
  unsigned int currentSnapOrbitalInRef = -1 ;
  
  int neighborhood = 10 ;
  

  if( noReferenceMode == NO )
  {

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
  
  }
  else
  {
    for( isnapOrbital = 0 ; isnapOrbital < nbasis ; isnapOrbital ++ )
    {
      *( snapMOMapping + isnapOrbital ) = isnapOrbital ;
    }
  
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
  
  // ----------------------> Read-in Dipole Matrices ... <------------------//

  lmo = nbasis * nbasis ;

  double * dipoleX = calloc( lmo , sizeof( double ) ) ; dzeros( lmo , 1 , dipoleX ) ;

  double * dipoleY = calloc( lmo , sizeof( double ) ) ; dzeros( lmo , 1 , dipoleY ) ;

  double * dipoleZ = calloc( lmo , sizeof( double ) ) ; dzeros( lmo , 1 , dipoleZ ) ;


  loadCNDOdipole( 'x' , pcndoLogFile , dipoleX ) ;
  
  loadCNDOdipole( 'y' , pcndoLogFile , dipoleY ) ;
  
  loadCNDOdipole( 'z' , pcndoLogFile , dipoleZ ) ;


  if( debuggingMode == YES )
  {
    debug = fopen( "dipoleX.deb" , "wb+") ;
    
    doutput( debug , nbasis , nbasis , dipoleX ) ;
    
    fclose( debug ) ;
  
  
    debug = fopen( "dipoleY.deb" , "wb+") ;
    
    doutput( debug , nbasis , nbasis , dipoleY ) ;
    
    fclose( debug ) ;
  
  
    debug = fopen( "dipoleZ.deb" , "wb+") ;
    
    doutput( debug , nbasis , nbasis , dipoleZ ) ;
    
    fclose( debug ) ;
    
  }





  // ----------------------> Convert into MO-Dipole-Integral ... <------------------//

  double * MOdipoleX = calloc( lmo , sizeof( double ) ) ; dzeros( lmo , 1 , MOdipoleX ) ;

  double * MOdipoleY = calloc( lmo , sizeof( double ) ) ; dzeros( lmo , 1 , MOdipoleY ) ;

  double * MOdipoleZ = calloc( lmo , sizeof( double ) ) ; dzeros( lmo , 1 , MOdipoleZ ) ;
  
  double * tmp_dipole = calloc( lmo , sizeof( double ) ) ; dzeros( lmo , 1 , tmp_dipole ) ;
  
  itmp2 = nbasis ;
  
  
  dgemm_( "N" , "T" , &itmp2 , &itmp2 , &itmp2 , &done , cndoMO , &itmp2 , dipoleX , &itmp2 , &dzero , tmp_dipole , &itmp2 ) ;
  
  dgemm_( "N" , "T" , &itmp2 , &itmp2 , &itmp2 , &done , tmp_dipole , &itmp2 , cndoMO , &itmp2 , &dzero , MOdipoleX , &itmp2 ) ;
  
  dtranspose( itmp2 , MOdipoleX , MOdipoleX ) ;
  

  dgemm_( "N" , "T" , &itmp2 , &itmp2 , &itmp2 , &done , cndoMO , &itmp2 , dipoleY , &itmp2 , &dzero , tmp_dipole , &itmp2 ) ;
  
  dgemm_( "N" , "T" , &itmp2 , &itmp2 , &itmp2 , &done , tmp_dipole , &itmp2 , cndoMO , &itmp2 , &dzero , MOdipoleY , &itmp2 ) ;
  
  dtranspose( itmp2 , MOdipoleY , MOdipoleY ) ;
  

  dgemm_( "N" , "T" , &itmp2 , &itmp2 , &itmp2 , &done , cndoMO , &itmp2 , dipoleZ , &itmp2 , &dzero , tmp_dipole , &itmp2 ) ;
  
  dgemm_( "N" , "T" , &itmp2 , &itmp2 , &itmp2 , &done , tmp_dipole , &itmp2 , cndoMO , &itmp2 , &dzero , MOdipoleZ , &itmp2 ) ;
  
  dtranspose( itmp2 , MOdipoleZ , MOdipoleZ ) ;
  
  

  if( debuggingMode == YES )
  {
    debug = fopen( "MOdipoleX.deb" , "wb+") ;
    
    doutput( debug , nbasis , nbasis , MOdipoleX ) ;
    
    fclose( debug ) ;
  
  
    debug = fopen( "MOdipoleY.deb" , "wb+") ;
    
    doutput( debug , nbasis , nbasis , MOdipoleY ) ;
    
    fclose( debug ) ;
  
  
    debug = fopen( "MOdipoleZ.deb" , "wb+") ;
    
    doutput( debug , nbasis , nbasis , MOdipoleZ ) ;
    
    fclose( debug ) ;

  }




  // ----------------------> Searching for the LE and CT state ... <------------------//

  /*
  int gsRefOrbital , leRefOrbital , ctRefOrbital ;
  
  int * degenerateStatesList[ 3 ] ;
  
  
  int  * gsRefOrbitalEqv , * leRefOrbitalEqv , * ctRefOrbitalEqv ;
  
  int  EXgsRefOrbitalEqv ,  EXleRefOrbitalEqv ,  EXctRefOrbitalEqv ;
  
  int  lenGSRefOrbitalEqv , lenLERefOrbitalEqv , lenCTRefOrbitalEqv ;
  

  EXgsRefOrbitalEqv = NO ;
  
  EXleRefOrbitalEqv = NO ;

  EXctRefOrbitalEqv = NO ;


  
  */

  //int casscfWidth = 10 ;
  
  int refGSinSnap , refLEinSnap , refCTinSnap ;
  
  
  int EXrefGSinSnap = NO , EXrefLEinSnap = NO , EXrefCTinSnap = NO ;
  
  
  
  for( isnapOrbital = 0 ; isnapOrbital < nbasis ; isnapOrbital ++ )
  {
    itmp = *( snapMOMapping + isnapOrbital )  ;
    
    if( itmp == gsRefOrbital - 1 )
    {
      refGSinSnap = isnapOrbital ;
      
      EXrefGSinSnap = YES ;
      
      printf("\n---> Found <--- @ This snapshot , # %d orbital is the GS ( # %d ) in reference ...\n" , isnapOrbital + 1 , gsRefOrbital ) ;
    }
    else if( itmp == leRefOrbital - 1 )
    {
      refLEinSnap = isnapOrbital ;
      
      EXrefLEinSnap = YES ;
      
      printf("\n---> Found <--- @ This snapshot , # %d orbital is the LE ( # %d ) in reference ...\n" , isnapOrbital + 1 , leRefOrbital ) ;
    }
    else if( itmp == ctRefOrbital - 1 )
    {
      refCTinSnap = isnapOrbital ;
      
      EXrefCTinSnap = YES ;
      
      printf("\n---> Found <--- @ This snapshot , # %d orbital is the CT ( # %d ) in reference ...\n" , isnapOrbital + 1 , ctRefOrbital ) ;
    }
    else
    {
      continue ;
    }
  
  
  }

  /*
  if( EXrefGSinSnap = NO || EXrefLEinSnap = NO || EXrefCTinSnap = NO )
  {
    printf("\nOrbital Matching Failed ... GMH Calculation Aborted ...\n") ;
    
    exit( 103 ) ;
  
  }
  */



  
  
  // -----> Do the actual GMH ...
  
  double * MU11 , * MU22 , * MU12 , * V0 ;
  
  double mu11 , mu22 , mu12 ;
  
  double * tmpMU ;
  
  double projectMU11 , projectMU22 , projectMU12 , v0 ;
  
  double dE ;
  
  double E_GS = 0.000 , E_CT = 0.000 , E_LE = 0.000 ;
  
  //double * TransDipole_GSCT , * TransDipole_LECT ;
  
  double TransDipole_GSCT , TransDipole_LECT ;
  
  double * Dipole_GS , * Dipole_CT , * Dipole_LE ;
  
  double HDA_GSLE = 0.000 ; double HDA_GSCT = 0.000 ; double HDA_LECT = 0.000 ;
  
  double HDA_GSLE_LINEAR = 0.000 ; double HDA_GSCT_LINEAR = 0.000 ; double HDA_LECT_LINEAR = 0.000 ;
  
  
  int flagGSCT = 0 ; int flagLECT = 0 ;


  MU11 = calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , MU11 ) ;
  
  MU22 = calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , MU22 ) ;
  
  MU12 = calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , MU12 ) ;
  
  tmpMU = calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , tmpMU ) ;
  
  V0 = calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , V0 ) ;
  
  

  if( EXrefGSinSnap == YES && EXrefCTinSnap == YES )  // Charge Re-Combination : CT -> GS 
  { 
    flagGSCT = 1 ;     

    dzeros( 3 , 1 , MU11 ) ;
  
    dzeros( 3 , 1 , MU22 ) ;
  
    dzeros( 3 , 1 , MU12 ) ;
  
    dzeros( 3 , 1 , tmpMU ) ;
  
    dzeros( 3 , 1 , V0 ) ;
  
    E_GS = ( *( MOEnergy + refGSinSnap ) ) ;
    
    E_CT = ( *( MOEnergy + refCTinSnap ) ) ;
    
    //TransDipole_GSCT = ( double * ) calloc( 3 , sizeof( double ) ) ;
    
    Dipole_GS = ( double * ) calloc( 3 , sizeof( double ) ) ;
    
    Dipole_CT = ( double * ) calloc( 3 , sizeof( double ) ) ;
    
    
    *( MU11 + 0 ) = *( MOdipoleX + refGSinSnap * nbasis + refGSinSnap ) ;
    
    *( Dipole_GS + 0 ) = *( MU11 + 0 ) ;
    
    *( MU11 + 1 ) = *( MOdipoleY + refGSinSnap * nbasis + refGSinSnap ) ;
    
    *( Dipole_GS + 1 ) = *( MU11 + 1 ) ;
    
    *( MU11 + 2 ) = *( MOdipoleZ + refGSinSnap * nbasis + refGSinSnap ) ;
    
    *( Dipole_GS + 2 ) = *( MU11 + 2 ) ;
    

    printf("\n[=============================================> GS - CT <=============================================]\n") ;

    printf("\nMU11 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( MU11 + 0 ) ,  *( MU11 + 1 )  , *( MU11 + 2 ) ) ;
    
    mu11 = sqrt( (*( MU11 + 0 )) * (*( MU11 + 0 )) +  (*( MU11 + 1 )) * (*( MU11 + 1 )) + (*( MU11 + 2 )) * (*( MU11 + 2 ))  ) ;
    
    *( MU22 + 0 ) = *( MOdipoleX + refCTinSnap * nbasis + refCTinSnap ) ;
    
    *( Dipole_CT + 0 ) = *( MU22 + 0 ) ; 
    
    *( MU22 + 1 ) = *( MOdipoleY + refCTinSnap * nbasis + refCTinSnap ) ;
    
    *( Dipole_CT + 1 ) = *( MU22 + 1 ) ; 
    
    *( MU22 + 2 ) = *( MOdipoleZ + refCTinSnap * nbasis + refCTinSnap ) ;
    
    *( Dipole_CT + 2 ) = *( MU22 + 2 ) ; 
    
    
    printf("\nMU22 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( MU22 + 0 ) ,  *( MU22 + 1 )  , *( MU22 + 2 ) ) ;
    
    mu22 = sqrt( (*( MU22 + 0 )) * (*( MU22 + 0 )) +  (*( MU22 + 1 )) * (*( MU22 + 1 )) + (*( MU22 + 2 )) * (*( MU22 + 2 ))  ) ;
    
    *( MU12 + 0 ) = *( MOdipoleX + refGSinSnap * nbasis + refCTinSnap ) ;
    
    *( MU12 + 1 ) = *( MOdipoleY + refGSinSnap * nbasis + refCTinSnap ) ;
    
    *( MU12 + 2 ) = *( MOdipoleZ + refGSinSnap * nbasis + refCTinSnap ) ;
    
    printf("\nMU12 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( MU12 + 0 ) ,  *( MU12 + 1 )  , *( MU12 + 2 ) ) ;
    
    mu12 = sqrt( (*( MU12 + 0 )) * (*( MU12 + 0 )) +  (*( MU12 + 1 )) * (*( MU12 + 1 )) + (*( MU12 + 2 )) * (*( MU12 + 2 ))  ) ;
    
    *( V0 + 0 ) = *( MU11 + 0 ) - *( MU22 + 0 ) ;
    
    *( V0 + 1 ) = *( MU11 + 1 ) - *( MU22 + 1 ) ;
    
    *( V0 + 2 ) = *( MU11 + 2 ) - *( MU22 + 2 ) ;
    
    printf("\nV0 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( V0 + 0 ) ,  *( V0 + 1 )  , *( V0 + 2 ) ) ;
    
    v0 = sqrt( (*( V0 + 0 )) * (*( V0 + 0 )) +  (*( V0 + 1 )) * (*( V0 + 1 )) + (*( V0 + 2 )) * (*( V0 + 2 ))  ) ;
    
    projectMU11 = ( ( *( MU11 + 0 ) ) * ( *( V0 + 0 ) ) + ( *( MU11 + 1 ) ) * ( *( V0 + 1 ) ) + ( *( MU11 + 2 ) ) * ( *( V0 + 2 ) ) ) / v0 ;
    
    projectMU22 = ( ( *( MU22 + 0 ) ) * ( *( V0 + 0 ) ) + ( *( MU22 + 1 ) ) * ( *( V0 + 1 ) ) + ( *( MU22 + 2 ) ) * ( *( V0 + 2 ) ) ) / v0 ;
    
    projectMU12 = ( ( *( MU12 + 0 ) ) * ( *( V0 + 0 ) ) + ( *( MU12 + 1 ) ) * ( *( V0 + 1 ) ) + ( *( MU12 + 2 ) ) * ( *( V0 + 2 ) ) ) / v0 ;
  
  
    TransDipole_GSCT = fabs( projectMU12 ) ;
    
    dE = E_CT - E_GS ;
    
    printf("\ndE is % 10.8f...\n" , dE ) ;
    
    printf("\nprojectMU12 = % 10.8f ...\n" , projectMU12 ) ;
    
    HDA_GSCT = fabs( projectMU12 ) * fabs( dE ) / sqrt( v0 * v0 + 4.00 * projectMU12 * projectMU12 ) ;
    
    HDA_GSCT_LINEAR = fabs( mu12 ) * fabs( dE ) / sqrt( v0 * v0 + 4.00 * mu12 * mu12 ) ;
    
    dtmp = fabs( mu12 ) * fabs( dE ) / sqrt( ( mu11 - mu22 ) * ( mu11 - mu22 ) + 4.00 * mu12 * mu12 ) ;
    
    //HDA_GSLE = HDA_GSLE * HDA_GSLE ;
    
    printf("\nHDA between GS and CT is % 10.8E ...\n" , HDA_GSCT ) ;
    
    printf("\nHDA between GS and CT is % 10.8E ( linear ) ...\n" , HDA_GSCT_LINEAR ) ;
    
    printf("\nHDA between GS and CT is % 10.8E ( crappy ) ...\n" , dtmp ) ;
    
    
  }  


  
  
  if( EXrefLEinSnap == YES && EXrefCTinSnap == YES ) // Charge Separation : LE -> CT
  {
    flagLECT = 1 ;     

    dzeros( 3 , 1 , MU11 ) ;
  
    dzeros( 3 , 1 , MU22 ) ;
  
    dzeros( 3 , 1 , MU12 ) ;
  
    dzeros( 3 , 1 , tmpMU ) ;
  
    dzeros( 3 , 1 , V0 ) ;
  
    E_LE = ( *( MOEnergy + refLEinSnap ) ) ;
    
    E_CT = ( *( MOEnergy + refCTinSnap ) ) ;
    
    
    Dipole_LE = ( double * ) calloc( 3 , sizeof( double ) ) ;
    
    //Dipole_CT = ( double * ) calloc( 3 , sizeof( double ) ) ;
    
    
    *( MU11 + 0 ) = *( MOdipoleX + refLEinSnap * nbasis + refLEinSnap ) ;
    
    *( Dipole_LE + 0 ) = *( MU11 + 0 ) ;
    
    *( MU11 + 1 ) = *( MOdipoleY + refLEinSnap * nbasis + refLEinSnap ) ;
    
    *( Dipole_LE + 1 ) = *( MU11 + 1 ) ;
    
    *( MU11 + 2 ) = *( MOdipoleZ + refLEinSnap * nbasis + refLEinSnap ) ;
    
    *( Dipole_LE + 2 ) = *( MU11 + 2 ) ;
    
    printf("\n[=============================================> LE - CT <=============================================]\n") ;

    printf("\nMU11 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( MU11 + 0 ) ,  *( MU11 + 1 )  , *( MU11 + 2 ) ) ;
    
    mu11 = sqrt( (*( MU11 + 0 )) * (*( MU11 + 0 )) +  (*( MU11 + 1 )) * (*( MU11 + 1 )) + (*( MU11 + 2 )) * (*( MU11 + 2 ))  ) ;
    
    *( MU22 + 0 ) = *( MOdipoleX + refCTinSnap * nbasis + refCTinSnap ) ;
    
    *( Dipole_CT + 0 ) = *( MU22 + 0 ) ;
    
    *( MU22 + 1 ) = *( MOdipoleY + refCTinSnap * nbasis + refCTinSnap ) ;
    
    *( Dipole_CT + 1 ) = *( MU22 + 1 ) ;
    
    *( MU22 + 2 ) = *( MOdipoleZ + refCTinSnap * nbasis + refCTinSnap ) ;
    
    *( Dipole_CT + 2 ) = *( MU22 + 2 ) ;
    
    
    printf("\nMU22 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( MU22 + 0 ) ,  *( MU22 + 1 )  , *( MU22 + 2 ) ) ;
    
    mu22 = sqrt( (*( MU22 + 0 )) * (*( MU22 + 0 )) +  (*( MU22 + 1 )) * (*( MU22 + 1 )) + (*( MU22 + 2 )) * (*( MU22 + 2 ))  ) ;
    
    *( MU12 + 0 ) = *( MOdipoleX + refCTinSnap * nbasis + refLEinSnap ) ;
    
    *( MU12 + 1 ) = *( MOdipoleY + refCTinSnap * nbasis + refLEinSnap ) ;
    
    *( MU12 + 2 ) = *( MOdipoleZ + refCTinSnap * nbasis + refLEinSnap ) ;
    
    printf("\nMU12 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( MU12 + 0 ) ,  *( MU12 + 1 )  , *( MU12 + 2 ) ) ;
    
    mu12 = sqrt( (*( MU12 + 0 )) * (*( MU12 + 0 )) +  (*( MU12 + 1 )) * (*( MU12 + 1 )) + (*( MU12 + 2 )) * (*( MU12 + 2 ))  ) ;
    
    *( V0 + 0 ) = *( MU11 + 0 ) - *( MU22 + 0 ) ;
    
    *( V0 + 1 ) = *( MU11 + 1 ) - *( MU22 + 1 ) ;
    
    *( V0 + 2 ) = *( MU11 + 2 ) - *( MU22 + 2 ) ;
    
    printf("\nV0 is % 10.8E  % 10.8E  % 10.8E ... \n" ,  *( V0 + 0 ) ,  *( V0 + 1 )  , *( V0 + 2 ) ) ;
    
    v0 = sqrt( (*( V0 + 0 )) * (*( V0 + 0 )) +  (*( V0 + 1 )) * (*( V0 + 1 )) + (*( V0 + 2 )) * (*( V0 + 2 ))  ) ;
    
    projectMU11 = ( ( *( MU11 + 0 ) ) * ( *( V0 + 0 ) ) + ( *( MU11 + 1 ) ) * ( *( V0 + 1 ) ) + ( *( MU11 + 2 ) ) * ( *( V0 + 2 ) ) ) / v0 ;
    
    projectMU22 = ( ( *( MU22 + 0 ) ) * ( *( V0 + 0 ) ) + ( *( MU22 + 1 ) ) * ( *( V0 + 1 ) ) + ( *( MU22 + 2 ) ) * ( *( V0 + 2 ) ) ) / v0 ;
    
    projectMU12 = ( ( *( MU12 + 0 ) ) * ( *( V0 + 0 ) ) + ( *( MU12 + 1 ) ) * ( *( V0 + 1 ) ) + ( *( MU12 + 2 ) ) * ( *( V0 + 2 ) ) ) / v0 ;
  
    
    TransDipole_LECT = fabs( projectMU12 ) ;
    
    dE = E_LE - E_CT ;
    
    printf("\ndE is % 10.8f...\n" , dE ) ;
    
    printf("\nprojectMU12 = % 10.8f ...\n" , projectMU12 ) ;
    
    HDA_LECT = fabs( projectMU12 ) * fabs( dE ) / sqrt( v0 * v0 + 4.00 * projectMU12 * projectMU12 ) ;
    
    HDA_LECT_LINEAR = fabs( mu12 ) * fabs( dE ) / sqrt( v0 * v0 + 4.00 * mu12 * mu12 ) ;
    
    dtmp = fabs( mu12 ) * fabs( dE ) / sqrt( ( mu11 - mu22 ) * ( mu11 - mu22 ) + 4.00 * mu12 * mu12 ) ;
    
    //HDA_GSLE = HDA_GSLE * HDA_GSLE ;
    
    printf("\nHDA between LE and CT is % 10.8E ...\n" , HDA_LECT ) ;
    
    printf("\nHDA between LE and CT is % 10.8E ( linear ) ...\n" , HDA_LECT_LINEAR ) ;
    
    printf("\nHDA between LE and CT is % 10.8E ( crappy ) ...\n" , dtmp ) ;
    
    
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
    
    if( exout != 19 )
    {
      strncpy( outputGMHName , cndoLogFileName , len_cndoLogFileName - 4 ) ;

      *( outputGMHName + len_cndoLogFileName - 4 ) = '\0' ;
    
      printf( "\n%s\n" , outputGMHName ) ;

      strcat( outputGMHName , ".gmh" ) ;
    }

    poutputGMHFile = fopen( outputGMHName , "wb+" ) ;

    fprintf( poutputGMHFile , "\n% 10.8E\t% 10.8E" , HDA_GSCT , HDA_LECT ) ;
    
    fprintf( poutputGMHFile , "\t% 10.8E\t% 10.8E\t% 10.8E" , E_GS , E_CT , E_LE ) ;
    
    fprintf( poutputGMHFile , "\t% 10.8E\t% 10.8E\t% 10.8E" , *( Dipole_GS + 0 ) , *( Dipole_GS + 1 ) , *( Dipole_GS + 2 ) ) ;
    
    fprintf( poutputGMHFile , "\t% 10.8E\t% 10.8E\t% 10.8E" , *( Dipole_CT + 0 ) , *( Dipole_CT + 1 ) , *( Dipole_CT + 2 ) ) ;
    
    fprintf( poutputGMHFile , "\t% 10.8E\t% 10.8E\t% 10.8E" , *( Dipole_LE + 0 ) , *( Dipole_LE + 1 ) , *( Dipole_LE + 2 ) ) ;

    fprintf( poutputGMHFile , "\t% 10.8E\t% 10.8E\n" , TransDipole_GSCT , TransDipole_LECT ) ;

    fclose( poutputGMHFile ) ;


  }




  return( 0 ) ;


}



