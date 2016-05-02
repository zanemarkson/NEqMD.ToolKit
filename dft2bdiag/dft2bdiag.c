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

#ifndef HARTREE2EV
#define HARTREE2EV 27.211
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

  // =====> Utility Functions Declaration ...

  double calcLocalization( double * orbital , int nBasis , int * basisList , int nComponent ) ;
  
  double calcSquaredLocalization( double * orbital , int nBasis , int * basisList , int nComponent );


  char ** pcmd = argv ; 
  
  int icmd ;
  
  char refMOFileName[ 300 ] , refCIListName[ 300 ] , cndoLogFileName[ 300 ] , indexFileName[ 300 ] , outputHDAName[ 300 ] ;

  int len_cndoLogFileName ;
  
  char basisMapFileName[ 300 ] ;
  
  int len_basisMapFileName ;
  
  FILE * prefMOFile , * prefCIListFile , * pcndoLogFile , * pIndexFile , * poutHDAFile ;
  
  FILE * pbasisMapFile ;
  
  double * refMO ; 
  
  int orbitalID ;
  
  int debuggingMode = NO ;

  char pGroupName[ 150 ] , qGroupName[ 150 ] , rGroupName[ 150 ] , xGroupName[ 150 ] , zGroupName[ 150 ];
  
  char dGroupName[ 150 ] , aGroupName[ 150 ] ;
  
  char d2GroupName[ 150 ] , a2GroupName[ 150 ] ;
  
  double * cndoFAO , * MOEnergy ;
  
  int refMOFileLength , refCIListFileLength , lmo ; // lmo is just for convenience = nbasis^2
  
  int HOMO , LUMO ;
  
  int nbasis , natom , nelectron , ncistate ;
  
  int ibasis , iatom , iao ;
  
  double done = 1.0000 ; double dzero = 0.0000 ; 
  
  int ithree = 3 ; int ione = 1 ; int izero = 0 ;
  


  char buffer[ MAXCHARINLINE ] ;
  
  int lenBuff ;
  
  char cache[ MAXCHARINLINE ] ;

  
  int itmp , itmp2 ; double dtmp ; char ctmp , tmp_char ; 
  
  char stmp[ 300 ] , stmp2[ 300 ] , tmpString[ 300 ];
  
  double dtmpArray[ 150 ] ;
  
  int info , signal , blank_signal ;
  
  register FILE * debug ; 
  
  
  int noRefSituation = NO ;
  
  // ----------------------------------> Recording Command-Line Arguments ... <---------------------------------- //
 
  time_t current_time;

  time( &current_time );

  char now[ 300 ] ;

  strcpy( now , ctime( &current_time ) );

  int lennow = strlen( now ) ;

  *( now + lennow - 1 ) = ' ';

    
    
    printf("\n**********************************************************************\n");
      printf("* G_DFT2BDIAG_D : Calculate B-Diag  HDA for Snapshots in MD Calcs.  *\n");
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
  
  strcpy( refMOFileName , "ref.MO" ) ;
  
  strcpy( refCIListName , "ref.CI" ) ;
  
  strcpy( basisMapFileName , "system.basis.map" ) ;
  
  strcpy( cndoLogFileName , "system.snap.log" ) ;
  
  strcpy( outputHDAName , "system.hda" ) ;
  
  
  strcpy( pGroupName , "Donor-Block" ) ;
  
  strcpy( qGroupName , "Acceptor-Block" ) ;
  
  strcpy( rGroupName , "Bridge-Block" ) ;
  
  strcpy( xGroupName , "Irrelevant-Block" ) ;
  
  strcpy( zGroupName , "Irrelevant-Block-2" ) ;
  
  
  // =====> Parsing cmd-line arguments ...
  
  int exR = NO , exX = NO , exZ = NO ;
  
  int exD = NO , exA = NO , exM = NO ;
  
  int exD2 = NO , exA2 = NO ;
  
  if( argc == 1 )
  {
    printf("\n\nNo command-line arguments provided ... Mission aborting ...\n\n");
    
    printf("\nPlease refer to the usage by typing ' %s -h '\n\n" , * argv );
    
    exit( 1 ); 

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
	      
	      case 'm' : strcpy( tmpString , *( ++ pcmd ) ) ; 
                     
                     if( strcmp( tmpString , "None" ) == 0 || strcmp( tmpString , "NONE" ) == 0 || strcmp( tmpString , "none" ) == 0 )
                     {
                       noRefSituation = YES ;
                       
                       printf("\nCommand-line arguments indicates : No reference will be checked. # of orbitals will be directly draged from CI List File...\n") ;
                     
                     }
                     else
                     {
	      			    strcpy( refMOFileName , tmpString ) ;
	      			    
	      			    printf("\nCommand-line argument indicates : Input ref. MO coefficient File name : %s ...\n" , refMOFileName ); 
	                 }
	                 
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
	      
	      case 'M' : strcpy( basisMapFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Input basis map file name : %s ...\n" , basisMapFileName ); 
	                 
	                 exM = YES ;
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ; 
	      
	      case 'o' : strcpy( outputHDAName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Output File name : %s ...\n" , outputHDAName ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 //exgro = 19 ;
	                 
	                 break ; 

  
          case 'n' : strcpy( indexFileName , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Input index File name : %s ...\n" , indexFileName ); 
	         
	                 icmd = icmd + 2 ;
	         
	                 break ;
	     
	      case 'p' : strcpy( pGroupName , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Group P Name is : %s ...\n" , pGroupName ); 
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;
	                 
	      case 'q' : strcpy( qGroupName , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Group Q Name is : %s ...\n" , qGroupName ); 
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;
	                 
	      case 'r' : strcpy( tmpString , *( ++ pcmd ) );
	      
	                 if( strcmp( tmpString , "None" ) == 0 || strcmp( tmpString , "NONE" ) == 0 || strcmp( tmpString , "none" ) == 0 )
	                 {
	                   exR = NO ; 
	                   
	                   printf("\nCommand-line argument indicates : There is NO Group X ...\n") ;
	                 }
	                 else
	                 {
	                   exR = YES ;
	                   
	                   strcpy( rGroupName , tmpString ) ;
	                   
	                   printf("\nCommand-line argument indicates : Group X Name is : %s ...\n" , rGroupName ); 
	                 }

	                 icmd = icmd + 2 ;
	                 
	                 break ;
	      
	      case 'x' : strcpy( tmpString , *( ++ pcmd ) );
	      
	                 if( strcmp( tmpString , "None" ) == 0 || strcmp( tmpString , "NONE" ) == 0 || strcmp( tmpString , "none" ) == 0 )
	                 {
	                   exX = NO ; 
	                   
	                   printf("\nCommand-line argument indicates : There is NO Group X ...\n") ;
	                 }
	                 else
	                 {
	                   exX = YES ;
	                   
	                   strcpy( xGroupName , tmpString ) ;
	                   
	                   printf("\nCommand-line argument indicates : Group X Name is : %s ...\n" , xGroupName ); 
	                 }
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;

	      case 'z' : strcpy( tmpString , *( ++ pcmd ) );
	      
	                 if( strcmp( tmpString , "None" ) == 0 || strcmp( tmpString , "NONE" ) == 0 || strcmp( tmpString , "none" ) == 0 )
	                 {
	                   exZ = NO ; 
	                   
	                   printf("\nCommand-line argument indicates : There is NO Group X ...\n") ;
	                 }
	                 else
	                 {
	                   exZ = YES ;
	                   
	                   strcpy( zGroupName , tmpString ) ;
	                   
	                   printf("\nCommand-line argument indicates : Group X Name is : %s ...\n" , zGroupName ); 
	                 }
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;


	      case 'D' : strcpy( tmpString , *( ++ pcmd ) );
	      
	                 if( strcmp( tmpString , "None" ) == 0 || strcmp( tmpString , "NONE" ) == 0 || strcmp( tmpString , "none" ) == 0 )
	                 {
	                   exD = NO ; 
	                   
	                   printf("\nCommand-line argument indicates : There is NO \"Donor\" fragment specification ...\n") ;
	                 }
	                 else
	                 {
	                   exD = YES ;
	                   
	                   strcpy( dGroupName , tmpString ) ;
	                   
	                   printf("\nCommand-line argument indicates : \"Donor\" fragment is called : %s ...\n" , dGroupName ); 
	                 }
	                 
	                 strcpy( tmpString , *( ++ pcmd ) );
	      
	                 if( strcmp( tmpString , "None" ) == 0 || strcmp( tmpString , "NONE" ) == 0 || strcmp( tmpString , "none" ) == 0 )
	                 {
	                   exD2 = NO ; 
	                   
	                   printf("\nCommand-line argument indicates : There is NO 2nd \"Donor\" fragment specification ...\n") ;
	                 }
	                 else
	                 {
	                   exD2 = YES ;
	                   
	                   strcpy( d2GroupName , tmpString ) ;
	                   
	                   printf("\nCommand-line argument indicates : 2nd \"Donor\" fragment is called : %s ...\n" , d2GroupName ); 
	                 }

	                 icmd = icmd + 3 ;
	                 
	                 break ;
	      

	      case 'A' : strcpy( tmpString , *( ++ pcmd ) );
	      
	                 if( strcmp( tmpString , "None" ) == 0 || strcmp( tmpString , "NONE" ) == 0 || strcmp( tmpString , "none" ) == 0 )
	                 {
	                   exA = NO ; 
	                   
	                   printf("\nCommand-line argument indicates : There is NO \"Acceptor\" fragment specification ...\n") ;
	                 }
	                 else
	                 {
	                   exA = YES ;
	                   
	                   strcpy( aGroupName , tmpString ) ;
	                   
	                   printf("\nCommand-line argument indicates : \"Acceptor\" fragment is called : %s ...\n" , aGroupName ); 
	                 }
	                 
	                 strcpy( tmpString , *( ++ pcmd ) );
	      
	                 if( strcmp( tmpString , "None" ) == 0 || strcmp( tmpString , "NONE" ) == 0 || strcmp( tmpString , "none" ) == 0 )
	                 {
	                   exA2 = NO ; 
	                   
	                   printf("\nCommand-line argument indicates : There is NO 2nd \"Acceptor\" fragment specification ...\n") ;
	                 }
	                 else
	                 {
	                   exA2 = YES ;
	                   
	                   strcpy( a2GroupName , tmpString ) ;
	                   
	                   printf("\nCommand-line argument indicates : 2nd \"Acceptor\" fragment is called : %s ...\n" , a2GroupName ); 
	                 }

	                 icmd = icmd + 3 ;
	                 
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
	    
	      case 'h' : printf("\nUsage:  %s [ -m 'input ref. MO Coefficient file name' ] [ -L 'Ref. CI list file name' ] \n\
	                                      [ -M 'input basis map file name'] [ -l 'Input cndo log file name' ]\n\
	                                      [ -o output file name ] [ -n GROMACS .ndx file name ] [ -g debugging mode (yes/y/YES/Yes or no/n/NO/No) ]\n\
	                                      [ -p Group P ] [ -q Group Q ] [ -r Group R ] [ -x Group X ] [ -z Group Z ]\n\
	                                      [ -D Name of Donor Fragment ] [ -A Name of Acceptor Fragment ]\n" , * argv ); 
	                 
                     printf("\n                Note : 1) when -m is specified as \"None\" or \"none\" or \"NONE\", no reference will be checked and orbital numbers provided in CI file will be directly used.\n\n");
                     printf("\n                Note : 2) Default group selection : [ -p \"Donor-Block\"] [ -q \"Acceptor-Block\"] [ -r \"Bridge-Block\"] [ -x \"Irrelevant-Block\"] [ -z \"Irrelevant-Block-2\"]\n\n\n") ; 
                     printf("\n                Note : 3) Default fragment names selection : [ -D \"Donor\"] [ -A \"Acceptor\"]\n\n\n") ; 
                     printf("\n                Note : 4) You must provide exactly two arguments following -D and -A. They could be \"none\" but they must be provided\n\n\n") ; 
                     printf("\n                Note : 5) If one of the partition is not necessary, \"None/none/None\" has to be specified ...\n\n\n") ; 


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
  
  if( ( prefMOFile = fopen( refMOFileName , "r" ) ) == NULL && noRefSituation == NO )
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
  
  if( ( pbasisMapFile = fopen( basisMapFileName , "r" ) ) == NULL )
  {
    printf("\nUser specified basis map File NOT FOUND ...\n");
    
    exit( 63 ) ;
  
  }
  
  
  if( ( pIndexFile = fopen( indexFileName , "r" ) ) == NULL )
  {
    printf("\nUser specified index File NOT FOUND ...\n");
    
    exit( 63 ) ;
  
  }
  
  
  // =====> Print out some information <===== //
  
  if( exD == NO )
  {
    printf("\nIn this Block-Diagonalization Analysis, MO will NOT be judged by localization for Donor ...\n\n") ;
  
  }
  
  if( exA == NO )
  {
    printf("\nIn this Block-Diagonalization Analysis, MO will NOT be judged by localization for Acceptor ...\n\n") ;
  
  }
  
  
  // =====> Loading Information from CNDO Log File and TRD File ... <===== //
  
  /*
  
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
  
  printf("\nThere are %d electrons ... So HOMO is # %d orbital ...\n\n" , nelectron , HOMO ) ;
  
  // ---> NAtom
  
  
  rewind( pcndoLogFile ) ;
  
  fsearch( pcndoLogFile , "RealAtoms=" ) ;
  
  fscanf( pcndoLogFile , "%d" , &natom );
  
  printf("\nThere are %d atoms in system ...\n" , natom ) ;
  
  
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

  */
  
  // ---> Fock in AO
  
  itmp = flength( pcndoLogFile ) ; 
  
  nbasis = sqrt( itmp ) ; 
  
  lmo = itmp ; 
  
  rewind( pcndoLogFile ) ;
  
  cndoFAO = ( double * ) calloc( lmo , sizeof( double ) ) ; dzeros( lmo , 1 , cndoFAO ) ;
  
  MOEnergy = ( double * ) calloc( nbasis , sizeof( double ) ) ; dzeros( nbasis , 1 , MOEnergy ) ;// in eV unit ...
  
  fload( pcndoLogFile , cndoFAO ) ;
  
  rewind( pcndoLogFile ) ;
  
  
  if ( debuggingMode == YES )
  {
    debug = fopen( "cndoFAO.deb" , "wb+" ) ;
  
    doutput( debug , nbasis , nbasis , cndoFAO ) ;
  
    fclose( debug ) ;
  }
  
  
  
  
  
  // ---> Basis function location information ... 
  
  itmp = flength( pbasisMapFile ) ; 
  
  natom = itmp / 5 ;
  
  rewind( pbasisMapFile ) ;
  
  
  int * tmpMap = calloc( 5 * natom , sizeof( int ) ) ;
  
  izeros( 5 * natom , 1 , tmpMap ) ;
  
  int * basisMap = calloc( 5 * natom , sizeof( int ) ) ;
  
  izeros( 5 * natom , 1 , basisMap ) ;
  
  /*
    For the "Basis Map" ==> 
  
    1st Column : N of Basis Function for Current Atom ;
    
    2nd Column : @ Input , "first bf" location ;
    
    3rd Column : After re-organization , "first bf" location ;
    
    4th Column : Atomic label ( number )
    
    
    
    For the original info from file ==>
    
    1st Column : 1 ~ natom
    
    2nd Column : Atomic label ( number )
    
    3rd Column : N of Basis Function for Current Atom ;
    
    4th Column : @ Input , "first bf" location ;
    
    5th Column : NElectron
    
  */
  
  
  int_fload( pbasisMapFile , tmpMap ) ; 
  
  rewind( pbasisMapFile ) ;
  
  int iline = 0 ;
  
  int iload = 0 ;

  
  for( iline = 0 ; iline < natom ; iline ++ )
  {
    *( basisMap + iline * 5 + 0 ) = *( tmpMap + iline * 5 + 2 ) ;
    
    *( basisMap + iline * 5 + 1 ) = *( tmpMap + iline * 5 + 3 ) ;
  
    // ---> To be determined : *( basisMap + iline * 5 + 2 )
    
    *( basisMap + iline * 5 + 3 ) = *( tmpMap + iline * 5 + 1 ) ;
    
    *( basisMap + iline * 5 + 4 ) = *( tmpMap + iline * 5 + 4 ) ;
    
    
  }  
  
  

  
  // -------------------------------> Read index file and acquire info of P, Q, R and X ... And Z ??? <---------------------------------- //
  
  
  int next_atom_index , current_atom_index ;
  
  int natomPGroup , natomQGroup , natomRGroup , natomXGroup , natomZGroup ;
  
  int ndxinfo ;
  
  // =====> P Group , Pre-Loading 
  
  rewind( pIndexFile ) ;
  
  info = fsearch( pIndexFile , pGroupName ) ;
  
  if( info == 0 )
  {
    printf("\nYour index file does not contain entry with name [ %s ] \n\n" , pGroupName ) ;
    
    exit( 61 ) ;
  }
  
  fskip( pIndexFile , 1 ) ;
  
  if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
  {
    printf("\nYour index file has hit the bottom with no specification for [ %s ] \n\n" , pGroupName ) ;
    
    exit( 62 ) ;
  }

  
  //printf("\nFor this group , the 1st grabbed tmpstring is %s ... \n" , tmpstring );
    
  //while( strcmp( tmpstring , "[") != 0 )
   
  for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
  {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
    if( strcmp( tmpString , "-" ) == 0 )
    {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
      fscanf( pIndexFile , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
        
      itmp = itmp + ( next_atom_index - current_atom_index );
        
      current_atom_index = next_atom_index ;
        
      if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
      {
         break ;  
      }

    }
    else
    {
      itmp ++ ;
        
      current_atom_index = atoi( tmpString );
        
      //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
      if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
      {
         break ;  
      }

    }

    

  }

  natomPGroup = itmp ;
  
  printf("\nBased on index file , P group %s has %d atoms ...\n" , pGroupName , natomPGroup );

  
  
  // =====> Q Group , Pre-Loading 
  
  rewind( pIndexFile ) ;
  
  info = fsearch( pIndexFile , qGroupName ) ;
  
  if( info == 0 )
  {
    printf("\nYour index file does not contain entry with name [ %s ] \n\n" , qGroupName ) ;
    
    exit( 61 ) ;
  }
  
  fskip( pIndexFile , 1 ) ;
  
  if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
  {
    printf("\nYour index file has hit the bottom with no specification for [ %s ] \n\n" , qGroupName ) ;
    
    exit( 62 ) ;
  }
  
  //printf("\nFor this group , the 1st grabbed tmpstring is %s ... \n" , tmpstring );
    
  //while( strcmp( tmpstring , "[") != 0 )
   
  for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
  {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
    if( strcmp( tmpString , "-" ) == 0 )
    {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
      fscanf( pIndexFile , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
        
      itmp = itmp + ( next_atom_index - current_atom_index );
        
      current_atom_index = next_atom_index ;
        
      if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
      {
         break ;  
      }

    }
    else
    {
      itmp ++ ;
        
      current_atom_index = atoi( tmpString );
        
      //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
      if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
      {
         break ;  
      }

    }

    

  }

  natomQGroup = itmp ;
  
  printf("\nBased on index file , Q group %s has %d atoms ...\n" , qGroupName , natomQGroup );

  
  // =====> R Group , Pre-Loading 
  
  rewind( pIndexFile ) ;
 
  info = fsearch( pIndexFile , rGroupName ) ;
  
  if( exR == NO )
  {
    if( info == 1 )  exR = YES ;    
  }
  else if( exR == YES )
  {
    if( info == 0 )
    {
      printf("\nYour index file does not contain entry with name [ %s ] \n\n" , rGroupName ) ;
      
      exit( 61 ) ;
    }

  }
  
  
  if( exR == YES )
  {
    rewind( pIndexFile ) ;
  
    info = fsearch( pIndexFile , rGroupName ) ;
    
    fskip( pIndexFile , 1 ) ;
  
    if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
    { 
      printf("\nYour index file has hit the bottom with no specification for [ %s ] \n\n" , rGroupName ) ;
    
      exit( 62 ) ;
    }
  
    //printf("\nFor this group , the 1st grabbed tmpstring is %s ... \n" , tmpstring );
    
    //while( strcmp( tmpstring , "[") != 0 )
   
    for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
    {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
      if( strcmp( tmpString , "-" ) == 0 )
      {
          //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
        
        itmp = itmp + ( next_atom_index - current_atom_index );
        
        current_atom_index = next_atom_index ;
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }
      else
      {
        itmp ++ ;
        
        current_atom_index = atoi( tmpString );
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }

    }

    natomRGroup = itmp ;
  
    printf("\nBased on index file , R group %s has %d atoms ...\n" , rGroupName , natomRGroup );
  
  }
  else
  {
    natomRGroup = 0 ;
  }

  
  // =====> X Group , Pre-Loading 
  
  rewind( pIndexFile ) ;
 
  info = fsearch( pIndexFile , xGroupName ) ;
  
  if( exX == NO )
  {
    if( info == 1 )  exX = YES ;    
  }
  else if( exX == YES )
  {
    if( info == 0 )
    {
      printf("\nYour index file does not contain entry with name [ %s ] \n\n" , xGroupName ) ;
      
      exit( 61 ) ;
    }

  }
  
  if( exX == YES )
  {
    rewind( pIndexFile ) ;
  
    info = fsearch( pIndexFile , xGroupName ) ;

    fskip( pIndexFile , 1 ) ;
  
    if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
    {
      printf("\nYour index file has hit the bottom with no specification for [ %s ] \n\n" , xGroupName ) ;
    
      exit( 62 ) ;
    }
  
    //printf("\nFor this group , the 1st grabbed tmpstring is %s ... \n" , tmpstring );
    
    //while( strcmp( tmpstring , "[") != 0 )
   
    for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
    {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
      if( strcmp( tmpString , "-" ) == 0 )
      {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
        
        itmp = itmp + ( next_atom_index - current_atom_index );
        
        current_atom_index = next_atom_index ;
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }
      else
      {
        itmp ++ ;
        
        current_atom_index = atoi( tmpString );
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }

    

    }

    natomXGroup = itmp ;
  
    printf("\nBased on index file , X group %s has %d atoms ...\n" , xGroupName , natomXGroup );
  
  }
  else
  {
    natomXGroup = 0 ;
  }

  

  // =====> Z Group , Pre-Loading 
  
  rewind( pIndexFile ) ;
 
  info = fsearch( pIndexFile , zGroupName ) ;
  
  if( exZ == NO )
  {
    if( info == 1 )  exZ = YES ;    
  }
  else if( exZ == YES )
  {
    if( info == 0 )
    {
      printf("\nYour index file does not contain entry with name [ %s ] \n\n" , zGroupName ) ;
      
      exit( 61 ) ;
    }

  }
  
  if( exZ == YES )
  {
    rewind( pIndexFile ) ;
  
    info = fsearch( pIndexFile , zGroupName ) ;

    fskip( pIndexFile , 1 ) ;
  
    if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
    {
      printf("\nYour index file has hit the bottom with no specification for [ %s ] \n\n" , zGroupName ) ;
    
      exit( 62 ) ;
    }
  
    //printf("\nFor this group , the 1st grabbed tmpstring is %s ... \n" , tmpstring );
    
    //while( strcmp( tmpstring , "[") != 0 )
   
    for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
    {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
      if( strcmp( tmpString , "-" ) == 0 )
      {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
        
        itmp = itmp + ( next_atom_index - current_atom_index );
        
        current_atom_index = next_atom_index ;
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }
      else
      {
        itmp ++ ;
        
        current_atom_index = atoi( tmpString );
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }

    

    }

    natomZGroup = itmp ;
  
    printf("\nBased on index file , Z group %s has %d atoms ...\n" , zGroupName , natomZGroup );
  
  }
  else
  {
    natomZGroup = 0 ;
  }












  if( natomPGroup + natomQGroup + natomRGroup + natomXGroup + natomZGroup != natom )
  {
    printf("\nAfter adding all available groups together ... we did not get a whole system ... Please check you index file ...\n");
    
    exit( 179 ) ;
  
  }
  
  
  

  
  // =====> P Group , Loading 
  
  //int * pGroupAtoms = calloc( natomPGroup , sizeof( int ) ) ;
  
  int * pGroupAtoms = calloc( natomPGroup , sizeof( int ) ) ;
  
  rewind( pIndexFile ) ;
  
  fsearch( pIndexFile , pGroupName ) ;
  
  fskip( pIndexFile , 1 ) ;
  
  fscanf( pIndexFile , "%s" , tmpString );

  for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
  {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
    if( strcmp( tmpString , "-" ) == 0 )
    {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
      fscanf( pIndexFile , "%d" , &next_atom_index );
        
      //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
      
      for( iatom = current_atom_index + 1 ; iatom <= next_atom_index ; iatom ++ )
      {
        *( pGroupAtoms + itmp + iatom - current_atom_index - 1 ) = iatom ;

      }
        
      itmp = itmp + ( next_atom_index - current_atom_index );
        
      current_atom_index = next_atom_index ;
        
      if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
      {
         break ;  
      }

    }
    else
    {   
      current_atom_index = atoi( tmpString );
      
      *( pGroupAtoms + itmp ) = current_atom_index ;
      
      itmp ++ ;
        
      //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
      if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
      {
         break ;  
      }

    }

    

  }


  
  // =====> Q Group , Loading 
  
  
  int * qGroupAtoms = calloc( natomQGroup , sizeof( int ) ) ;
  
  rewind( pIndexFile ) ;
  
  fsearch( pIndexFile , qGroupName ) ;
  
  fskip( pIndexFile , 1 ) ;
  
  fscanf( pIndexFile , "%s" , tmpString );

  for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
  {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
    if( strcmp( tmpString , "-" ) == 0 )
    {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
      fscanf( pIndexFile , "%d" , &next_atom_index );
        
      //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
      
      for( iatom = current_atom_index + 1 ; iatom <= next_atom_index ; iatom ++ )
      {
        *( qGroupAtoms + itmp + iatom - current_atom_index - 1 ) = iatom ;

      }
        
      itmp = itmp + ( next_atom_index - current_atom_index );
        
      current_atom_index = next_atom_index ;
        
      if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
      {
         break ;  
      }

    }
    else
    {   
      current_atom_index = atoi( tmpString );
      
      *( qGroupAtoms + itmp ) = current_atom_index ;
      
      itmp ++ ;
        
      //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
      if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
      {
         break ;  
      }

    }

    

  }



 
  // =====> R Group , Loading 

  int * rGroupAtoms = calloc( natomRGroup , sizeof( int ) ) ;
  
  if( exR == YES )
  {
    rewind( pIndexFile ) ;
  
    fsearch( pIndexFile , rGroupName ) ;
  
    fskip( pIndexFile , 1 ) ;
  
    fscanf( pIndexFile , "%s" , tmpString );

    for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
    {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
      if( strcmp( tmpString , "-" ) == 0 )
      {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
      
        for( iatom = current_atom_index + 1 ; iatom <= next_atom_index ; iatom ++ )
        {
          *( rGroupAtoms + itmp + iatom - current_atom_index - 1 ) = iatom ;

        }
        
        itmp = itmp + ( next_atom_index - current_atom_index );
        
        current_atom_index = next_atom_index ;
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }
      else
      {     
        current_atom_index = atoi( tmpString );
      
        *( rGroupAtoms + itmp ) = current_atom_index ;
      
        itmp ++ ;
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }
      }

    }
  
  }

  // =====> X Group , Loading 

  int * xGroupAtoms = calloc( natomXGroup , sizeof( int ) ) ;
  
  if( exX == YES )
  {
    rewind( pIndexFile ) ;
  
    fsearch( pIndexFile , xGroupName ) ;
  
    fskip( pIndexFile , 1 ) ;
  
    fscanf( pIndexFile , "%s" , tmpString );

    for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
    {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
      if( strcmp( tmpString , "-" ) == 0 )
      {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
      
        for( iatom = current_atom_index + 1 ; iatom <= next_atom_index ; iatom ++ )
        {
          *( xGroupAtoms + itmp + iatom - current_atom_index - 1 ) = iatom ;

        }
        
        itmp = itmp + ( next_atom_index - current_atom_index );
        
        current_atom_index = next_atom_index ;
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }
      else
      {   
        current_atom_index = atoi( tmpString );
      
        *( xGroupAtoms + itmp ) = current_atom_index ;
      
        itmp ++ ;
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }

    }
  
  }



  // =====> Z Group , Loading 

  int * zGroupAtoms = calloc( natomZGroup , sizeof( int ) ) ;
  
  if( exZ == YES )
  {
    rewind( pIndexFile ) ;
  
    fsearch( pIndexFile , zGroupName ) ;
  
    fskip( pIndexFile , 1 ) ;
  
    fscanf( pIndexFile , "%s" , tmpString );

    for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
    {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
      if( strcmp( tmpString , "-" ) == 0 )
      {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
      
        for( iatom = current_atom_index + 1 ; iatom <= next_atom_index ; iatom ++ )
        {
          *( zGroupAtoms + itmp + iatom - current_atom_index - 1 ) = iatom ;

        }
        
        itmp = itmp + ( next_atom_index - current_atom_index );
        
        current_atom_index = next_atom_index ;
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }
      else
      {   
        current_atom_index = atoi( tmpString );
      
        *( zGroupAtoms + itmp ) = current_atom_index ;
      
        itmp ++ ;
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }

    }
  
  }






// -------------------------------> Deal with reference MO ... <---------------------------------- //

  refMO = ( double * ) calloc( nbasis * nbasis , sizeof( double ) ) ;
  
  dzeros( nbasis , nbasis , refMO ) ;
  

  if( noRefSituation == NO )
  {
    refMOFileLength = flength( prefMOFile ) ;
  
    rewind( prefMOFile ) ;
  
    if( refMOFileLength != lmo )
    {
      printf("\nEm... it is a little strange ... there are %d numbers in refMO file but NBasis is %d ...\n" , refMOFileLength , nbasis ) ;
    
      exit( 12 ) ;
  
    }
    else
    {
      fload( prefMOFile , refMO ) ;
    
      rewind( prefMOFile ) ;
  
    }
  }
  
  
  
// -------------------------------> Start to Move the Fock Matrix ... <---------------------------------- //
  
  
  double * orgFAO = calloc( lmo , sizeof( double) ) ; dzeros( nbasis , nbasis , orgFAO ) ;
  
  double * tmpFAO = calloc( lmo , sizeof( double) ) ; dzeros( nbasis , nbasis , tmpFAO ) ;
  
  int atomID , startFrom , howManyFn , newStart ;
  
  int nbasisP , nbasisQ , nbasisR , nbasisX , nbasisZ ;
  
  int nelectronP , nelectronQ , nelectronR , nelectronX , nelectronZ ;
  
  int lumoP = 0 , homoP = 0 , lumoQ = 0 , homoQ = 0 , lumoR = 0 , homoR = 0 , lumoX = 0 , homoX = 0 , lumoZ = 0 , homoZ = 0 ; 

  
  int accu = 1 ;
  
  
  nbasisP = 0 ; nelectronP = 0 ;
  
  for( iatom = 0 ; iatom < natomPGroup ; iatom ++ )
  {
    atomID = *( pGroupAtoms + iatom ) - 1 ; // C-Labeling ... 
  
    howManyFn = *( basisMap + 5 * atomID + 0 ) ; // But startFrom and newStart are in human labeling...
    
    startFrom = *( basisMap + 5 * atomID + 1 ) ;
    
    newStart = accu ;
    
    *( basisMap + 5 * atomID + 2 ) = newStart ;
    
    *( basisMap + 5 * atomID + 3 ) = newStart + howManyFn - 1 ;
    
    accu = accu + howManyFn ; 
    
    nbasisP = nbasisP + howManyFn ;
    
    nelectronP = nelectronP + ( *( basisMap + 5 * atomID + 4 ) ) ;
  
  
  }
  
  homoP = ( nelectronP + ( nelectronP % 2 ) ) / 2 ; // Human-Label
  
  lumoP = homoP + 1 ; // Human-Label
  

  

  nbasisQ = 0 ; nelectronQ = 0 ;
  
  for( iatom = 0 ; iatom < natomQGroup ; iatom ++ )
  {
    atomID = *( qGroupAtoms + iatom ) - 1 ; // C-Labeling ... 
  
    howManyFn = *( basisMap + 5 * atomID + 0 ) ; // But startFrom and newStart are in human labeling...
    
    startFrom = *( basisMap + 5 * atomID + 1 ) ;
    
    newStart = accu ;
    
    *( basisMap + 5 * atomID + 2 ) = newStart ;
    
    *( basisMap + 5 * atomID + 3 ) = newStart + howManyFn - 1 ;
    
    accu = accu + howManyFn ; 
    
    nbasisQ = nbasisQ + howManyFn ;
    
    nelectronQ = nelectronQ + ( *( basisMap + 5 * atomID + 4 ) ) ;
  
  
  }
  
  homoQ = ( nelectronQ + ( nelectronQ % 2 ) ) / 2 ; // Human-Label
  
  lumoQ = homoQ + 1 ; // Human-Label



  
  nbasisR = 0 ; nelectronR = 0 ;
  
  for( iatom = 0 ; iatom < natomRGroup ; iatom ++ )
  {
    atomID = *( rGroupAtoms + iatom ) - 1 ; // C-Labeling ... 
  
    howManyFn = *( basisMap + 5 * atomID + 0 ) ; // But startFrom and newStart are in human labeling...
    
    startFrom = *( basisMap + 5 * atomID + 1 ) ;
    
    newStart = accu ;
    
    *( basisMap + 5 * atomID + 2 ) = newStart ;
    
    *( basisMap + 5 * atomID + 3 ) = newStart + howManyFn - 1 ;
    
    accu = accu + howManyFn ; 
    
    nbasisR = nbasisR + howManyFn ;
  
    nelectronR = nelectronR + ( *( basisMap + 5 * atomID + 4 ) ) ;
  
  }
  

  homoR = ( nelectronR + ( nelectronR % 2 ) ) / 2 ; // Human-Label
  
  lumoR = homoR + 1 ; // Human-Label

  
  
  nbasisX = 0 ; nelectronX = 0 ;
  
  for( iatom = 0 ; iatom < natomXGroup ; iatom ++ )
  {
    atomID = *( xGroupAtoms + iatom ) - 1 ; // C-Labeling ... 
  
    howManyFn = *( basisMap + 5 * atomID + 0 ) ; // But startFrom and newStart are in human labeling...
    
    startFrom = *( basisMap + 5 * atomID + 1 ) ;
    
    newStart = accu ;
    
    *( basisMap + 5 * atomID + 2 ) = newStart ;
    
    *( basisMap + 5 * atomID + 3 ) = newStart + howManyFn - 1 ;
    
    accu = accu + howManyFn ; 
    
    nbasisX = nbasisX + howManyFn ;
    
    nelectronX = nelectronX + ( *( basisMap + 5 * atomID + 4 ) ) ;
  
  }
  
  homoX = ( nelectronX + ( nelectronX % 2 ) ) / 2 ; // Human-Label
  
  lumoX = homoX + 1 ; // Human-Label
  


  nbasisZ = 0 ; nelectronZ = 0 ; 
  
  for( iatom = 0 ; iatom < natomZGroup ; iatom ++ )
  {
    atomID = *( zGroupAtoms + iatom ) - 1 ; // C-Labeling ... 
  
    howManyFn = *( basisMap + 5 * atomID + 0 ) ; // But startFrom and newStart are in human labeling...
    
    startFrom = *( basisMap + 5 * atomID + 1 ) ;
    
    newStart = accu ;
    
    *( basisMap + 5 * atomID + 2 ) = newStart ;
    
    *( basisMap + 5 * atomID + 3 ) = newStart + howManyFn - 1 ;
    
    accu = accu + howManyFn ; 
    
    nbasisZ = nbasisZ + howManyFn ;
    
    nelectronZ = nelectronZ + ( *( basisMap + 5 * atomID + 4 ) ) ;
  
  }

  homoZ = ( nelectronZ + ( nelectronZ % 2 ) ) / 2 ; // Human-Label
  
  lumoZ = homoZ + 1 ; // Human-Label


  
  
  if ( debuggingMode == YES )
  {
    debug = fopen( "basisMap.deb" , "wb+" ) ;
  
    ioutput( debug , natom , 5 , basisMap ) ;
  
    fclose( debug ) ;
  }  
  
  printf("\n[ NBASIS ] : P = %d ; Q = %d ; R = %d ; X = %d ; Z = %d : [ NBASIS ]\n\n" , nbasisP , nbasisQ , nbasisR , nbasisX , nbasisZ ) ;
  
  printf("\n[ NElectron ] : P = %d ; Q = %d ; R = %d ; X = %d ; Z = %d : [ NElectron ]\n\n" , nelectronP , nelectronQ , nelectronR , nelectronX , nelectronZ ) ;
  
  
  
  
  for( iatom = 0 ; iatom < natom ; iatom ++ ) // Moving Columns ...
  {
    howManyFn = *( basisMap + 5 * iatom + 0 ) ;
    
    startFrom = *( basisMap + 5 * iatom + 1 ) ;
    
    newStart = *( basisMap + 5 * iatom + 2 ) ;
    
    
    for( ibasis = 0 ; ibasis < howManyFn ; ibasis ++ )
    {
      for( iao = 0 ; iao < nbasis ; iao ++ )
      {
        *( tmpFAO + iao * nbasis + ( newStart - 1 + ibasis ) ) = *( cndoFAO + iao * nbasis + ( startFrom - 1 + ibasis ) ) ;
      }
    
    }
    
  
  }
  
  /*
  debug = fopen( "tmpFAO.deb" , "wb+" );
  
  doutput( debug , nbasis , nbasis , tmpFAO ) ;
  
  fclose( debug ) ;
  */
    
  
  
  for( iatom = 0 ; iatom < natom ; iatom ++ ) // Moving Rows ... 
  {
    howManyFn = *( basisMap + 5 * iatom + 0 ) ;
    
    startFrom = *( basisMap + 5 * iatom + 1 ) ;
    
    newStart = *( basisMap + 5 * iatom + 2 ) ;
    
    
    for( ibasis = 0 ; ibasis < howManyFn ; ibasis ++ )
    {
      for( iao = 0 ; iao < nbasis ; iao ++ )
      {
        *( orgFAO + ( newStart - 1 + ibasis ) * nbasis + iao ) = *( tmpFAO + ( startFrom - 1 + ibasis ) * nbasis + iao ) ;
      
      }
    
    }
  
  
  }
  
  
  
  if ( debuggingMode == YES )
  {
    debug = fopen( "orgFAO.deb" , "wb+" );
  
    doutput( debug , nbasis , nbasis , orgFAO ) ;
  
    fclose( debug ) ;
  }
  
  
  
  double * orgMO = calloc( nbasis * nbasis , sizeof( double ) ) ; dzeros( nbasis , nbasis , orgMO ) ;
  
  double * orgEnergy = calloc( nbasis , sizeof( double) ) ; dzeros( nbasis , 1 , orgEnergy ) ;
  
  dsyev_f2c( nbasis , orgFAO , orgMO , orgEnergy ) ;
  
  dtranspose( nbasis , orgMO , orgMO ) ;
  
  
  
  if ( debuggingMode == YES )
  {
    debug = fopen( "orgEnergy.deb" , "wb+" );
  
    doutput( debug , nbasis , 1 , orgEnergy ) ;
  
    fclose( debug ) ;

  
    debug = fopen( "orgMO.deb" , "wb+" );
  
    doutput( debug , nbasis , nbasis , orgMO ) ;
  
    fclose( debug ) ;
  
  }
  
  // -------------------------------> Divide into Blocks ... <---------------------------------- //
  
  double * reorgEnergy = calloc( nbasis , sizeof( double ) ) ; dzeros( nbasis , 1 , reorgEnergy ) ;
  
  //=> P
  
  double * faoP = calloc( nbasisP * nbasisP , sizeof( double ) ) ; dzeros( nbasisP , nbasisP , faoP ) ;
  
  double * moP = calloc( nbasisP * nbasisP , sizeof( double ) ) ; dzeros( nbasisP , nbasisP , moP ) ;
  
  double * energyP = calloc( nbasisP , sizeof( double ) ) ; dzeros( nbasisP , 1 , energyP ) ;
 
  
  for( ibasis = 0 ; ibasis < nbasisP ; ibasis ++ )
  {
    for( iao = 0 ; iao < nbasisP ; iao ++ )
    {
      *( faoP + ibasis * nbasisP + iao ) = *( orgFAO + ibasis * nbasis + iao ) ;
    }
  
  }
  
  //void dsyev_f2c( int dimension, double * sym, double * eigvectors,  double * eigvalues  );
  
  dsyev_f2c( nbasisP , faoP , moP , energyP ) ;
  
  dtranspose( nbasisP , moP , moP ) ;
  
  for( ibasis = 0 ; ibasis < nbasisP ; ibasis ++ )
  {
    *( reorgEnergy + ibasis ) = *( energyP + ibasis ) ;
  
  }
  
  
  /*
  debug = fopen( "moP.deb" , "wb+") ;
  
  doutput( debug , nbasisP , nbasisP , moP ) ;
  
  fclose( debug ) ;
  
  debug = fopen( "energyP.deb" , "wb+" ) ;
  
  doutput( debug , nbasisP , 1 , energyP ) ;
  
  fclose( debug ) ;
  */


  //=> Q

  double * faoQ = calloc( nbasisQ * nbasisQ , sizeof( double ) ) ; dzeros( nbasisQ , nbasisQ , faoQ ) ;

  double * moQ = calloc( nbasisQ * nbasisQ , sizeof( double ) ) ; dzeros( nbasisQ , nbasisQ , moQ ) ;
  
  double * energyQ = calloc( nbasisQ , sizeof( double ) ) ; dzeros( nbasisQ , 1 , energyQ ) ;

  
  for( ibasis = 0 ; ibasis < nbasisQ ; ibasis ++ )
  {
    for( iao = 0 ; iao < nbasisQ ; iao ++ )
    {
      *( faoQ + ibasis * nbasisQ + iao ) = *( orgFAO + ( nbasisP + ibasis ) * nbasis + ( nbasisP + iao ) ) ;
    }
  }
  
  //void dsyev_f2c( int dimension, double * sym, double * eigvectors,  double * eigvalues  );
  
  dsyev_f2c( nbasisQ , faoQ , moQ , energyQ ) ;
  
  dtranspose( nbasisQ , moQ , moQ ) ;
  
  for( ibasis = 0 ; ibasis < nbasisQ ; ibasis ++ )
  {
    *( reorgEnergy + nbasisP + ibasis ) = *( energyQ + ibasis ) ;
  
  }


  /*
  debug = fopen( "moQ.deb" , "wb+") ;
  
  doutput( debug , nbasisQ , nbasisQ , moQ ) ;
  
  fclose( debug ) ;
  
  debug = fopen( "energyQ.deb" , "wb+" ) ;
  
  doutput( debug , nbasisQ , 1 , energyQ ) ;
  
  fclose( debug ) ;
  */
  
  

  //=> R

  double * faoR , * moR , * energyR ;
 
  if( exR == YES )
  {
   faoR = ( double * ) calloc( nbasisR * nbasisR , sizeof( double ) ) ; dzeros( nbasisR , nbasisR , faoR ) ;
   
   moR = ( double * ) calloc( nbasisR * nbasisR , sizeof( double ) ) ; dzeros( nbasisR , nbasisR , moR ) ;
   
   energyR = ( double * ) calloc( nbasisR , sizeof( double ) ) ; dzeros( nbasisR , 1 , energyR ) ;
   
   for( ibasis = 0 ; ibasis < nbasisR ; ibasis ++ )
    {
      for( iao = 0 ; iao < nbasisR ; iao ++ )
      {
        *( faoR + ibasis * nbasisR + iao ) = *( orgFAO + ( nbasisP + nbasisQ + ibasis ) * nbasis + ( nbasisP + nbasisQ + iao ) ) ;
      }
    }
  
    //void dsyev_f2c( int dimension, double * sym, double * eigvectors,  double * eigvalues  );
  
    dsyev_f2c( nbasisR , faoR , moR , energyR ) ;
  
    dtranspose( nbasisR , moR , moR ) ;
  
    for( ibasis = 0 ; ibasis < nbasisR ; ibasis ++ )
    {
      *( reorgEnergy + nbasisP + nbasisQ + ibasis ) = *( energyR + ibasis ) ;
  
    }

    /*
    if( nbasisR != 0 )
    {
      debug = fopen( "moR.deb" , "wb+") ;
  
      doutput( debug , nbasisR , nbasisR , moR ) ;
  
      fclose( debug ) ;
  
      debug = fopen( "energyR.deb" , "wb+" ) ;
  
      doutput( debug , nbasisR , 1 , energyR ) ;
  
     fclose( debug ) ;
    }
  */
  
  }

  //=> X
    
  double * faoX , * moX , * energyX ;  

  if( exX == YES )
  {
    faoX = ( double * ) calloc( nbasisX * nbasisX , sizeof( double ) ) ; dzeros( nbasisX , nbasisX , faoX ) ;
    
    moX = ( double * ) calloc( nbasisX * nbasisX , sizeof( double ) ) ; dzeros( nbasisX , nbasisX , moX ) ;

    energyX = ( double * ) calloc( nbasisX , sizeof( double ) ) ; dzeros( nbasisX , 1 , energyX ) ;
 
    for( ibasis = 0 ; ibasis < nbasisX ; ibasis ++ )
    {
      for( iao = 0 ; iao < nbasisX ; iao ++ )
      {
        *( faoX + ibasis * nbasisX + iao ) = *( orgFAO + ( nbasisP + nbasisR + nbasisQ + ibasis ) * nbasis + ( nbasisP + nbasisR + + nbasisQ + iao ) ) ;
      }
    }
    
    
    
    dsyev_f2c( nbasisX , faoX , moX , energyX ) ;
    
    dtranspose( nbasisX , moX , moX ) ;
    
    for( ibasis = 0 ; ibasis < nbasisX ; ibasis ++ )
    {
      *( reorgEnergy + nbasisP + nbasisQ + nbasisR + ibasis ) = *( energyX + ibasis ) ;
    
    }
    
  /*
    if( nbasisX != 0 )
    {
      debug = fopen( "moX.deb" , "wb+") ;
    
      doutput( debug , nbasisX , nbasisX , moX ) ;
    
      fclose( debug ) ;
    
      debug = fopen( "energyX.deb" , "wb+" ) ;
    
      doutput( debug , nbasisX , 1 , energyX ) ;
    
      fclose( debug ) ;
    }
    */
    
  }

  //=> Z
  
  double * faoZ , * moZ , * energyZ ;

  if( exZ == YES )
  {
    faoZ = ( double * ) calloc( nbasisZ * nbasisZ , sizeof( double ) ) ; dzeros( nbasisZ , nbasisZ , faoZ ) ;
    
    moZ = ( double * ) calloc( nbasisZ * nbasisZ , sizeof( double ) ) ; dzeros( nbasisZ , nbasisZ , moZ ) ;
    
    energyZ = ( double * ) calloc( nbasisZ , sizeof( double ) ) ; dzeros( nbasisZ , 1 , energyZ ) ;


    for( ibasis = 0 ; ibasis < nbasisZ ; ibasis ++ )
    {
      for( iao = 0 ; iao < nbasisZ ; iao ++ )
      {
        *( faoZ + ibasis * nbasisZ + iao ) = *( orgFAO + ( nbasisP + nbasisR + nbasisQ + nbasisX + ibasis ) * nbasis + ( nbasisP + nbasisR + + nbasisQ + nbasisX + iao ) ) ;
      }
    }
  
  
  
    dsyev_f2c( nbasisZ , faoZ , moZ , energyZ ) ;
  
    dtranspose( nbasisZ , moZ , moZ ) ;
  
    for( ibasis = 0 ; ibasis < nbasisZ ; ibasis ++ )
    {
      *( reorgEnergy + nbasisP + nbasisQ + nbasisR + nbasisX + ibasis ) = *( energyZ + ibasis ) ;
  
    }
  
  }
  
  
  
  
  double * MO = calloc( nbasis * nbasis , sizeof( double ) ) ; dzeros( nbasis , nbasis , MO ) ;
  
  for( ibasis = 0 ; ibasis < nbasisP ; ibasis ++ )
  {
    for( iao = 0 ; iao < nbasisP ; iao ++ )
    {
      *( MO + ibasis * nbasis + iao ) = *( moP + ibasis * nbasisP + iao ) ;
    }
  }
  
  for( ibasis = 0 ; ibasis < nbasisQ ; ibasis ++ )
  {
    for( iao = 0 ; iao < nbasisQ ; iao ++ )
    {
      *( MO + ( nbasisP + ibasis ) * nbasis + ( nbasisP + iao ) ) = *( moQ + ibasis * nbasisQ + iao ) ;
    }
  }
  

  for( ibasis = 0 ; ibasis < nbasisR ; ibasis ++ )
  {
    for( iao = 0 ; iao < nbasisR ; iao ++ )
    {
      *( MO + ( nbasisP + nbasisQ + ibasis ) * nbasis + ( nbasisP + nbasisQ + iao ) ) = *( moR + ibasis * nbasisR + iao ) ;
    }
  }
  
  
  for( ibasis = 0 ; ibasis < nbasisX ; ibasis ++ )
  {
    for( iao = 0 ; iao < nbasisX ; iao ++ )
    {
      *( MO + ( nbasisP + nbasisQ + nbasisR + ibasis ) * nbasis + ( nbasisP + nbasisQ + nbasisR + iao ) ) = *( moX + ibasis * nbasisX + iao ) ;
    }
  }
  
  for( ibasis = 0 ; ibasis < nbasisZ ; ibasis ++ )
  {
    for( iao = 0 ; iao < nbasisZ ; iao ++ )
    {
      *( MO + ( nbasisP + nbasisQ + nbasisR + nbasisX + ibasis ) * nbasis + ( nbasisP + nbasisQ + nbasisR + nbasisX + iao ) ) = *( moZ + ibasis * nbasisZ + iao ) ;
    }
  }
  
  
  double * finalFock = calloc( nbasis * nbasis , sizeof( double ) ) ; dzeros( nbasis , nbasis , finalFock ) ;
  
  dzeros( nbasis , nbasis , tmpFAO ) ;
  
  dgemm_( "N" , "T" , &nbasis , &nbasis , &nbasis , &done , MO , &nbasis , orgFAO , &nbasis , &dzero , tmpFAO , &nbasis ) ;
  
  dgemm_( "N" , "T" , &nbasis , &nbasis , &nbasis , &done , tmpFAO , &nbasis , MO , &nbasis , &dzero , finalFock , &nbasis ) ;
  
  dtranspose( nbasis , finalFock , finalFock ) ;
  





  for( ibasis = 0 ; ibasis < nbasis ; ibasis ++ )
  {
    *( reorgEnergy + ibasis ) = ( *( reorgEnergy + ibasis ) ) * HARTREE2EV ;
  }



  if ( debuggingMode == YES )
  {
    debug = fopen( "ReOrgMO.deb" , "wb+" );
  
    doutput( debug , nbasis , nbasis , MO ) ;
  
    fclose( debug ) ;
  
    debug = fopen( "ReOrgFAO.deb" , "wb+" );
  
    doutput( debug , nbasis , nbasis , finalFock ) ;
  
    fclose( debug ) ;



    debug = fopen( "ReOrgEnergy.deb" , "wb+" );
  
    doutput( debug , nbasis , 1 , reorgEnergy ) ;
  
    fclose( debug ) ;
  
  }




  // ------------------> Define Fragments for MO Localization Analysis ... <---------------------- //
  
  int natomDGroup , natomAGroup ;
  
  int natomD2Group , natomA2Group ;
  
  
  
  // =====> D Fragment , Pre-Loading 
  
  rewind( pIndexFile ) ;
 
  /*
 
  info = fsearch( pIndexFile , dGroupName ) ;
  
  rewind( pIndexFile ) ;
 
  if( exD == NO )
  {
    if( info == 1 )  exD = YES ;    
  }
  else */if( exD == YES )
  {
    info = fsearch( pIndexFile , dGroupName ) ;
    
    if( info == 0 )
    {
      printf("\nYour index file does not contain entry with name [ %s ] \n\n" , dGroupName ) ;
      
      exit( 61 ) ;
    }

    rewind( pIndexFile ) ;
  
    info = fsearch( pIndexFile , dGroupName ) ;
    
    fskip( pIndexFile , 1 ) ;
  
    if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
    { 
      printf("\nYour index file has hit the bottom with no specification for [ %s ] \n\n" , dGroupName ) ;
    
      exit( 62 ) ;
    }
  
    //printf("\nFor this group , the 1st grabbed tmpstring is %s ... \n" , tmpstring );
    
    //while( strcmp( tmpstring , "[") != 0 )
   
    for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
    {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
      if( strcmp( tmpString , "-" ) == 0 )
      {
          //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
        
        itmp = itmp + ( next_atom_index - current_atom_index );
        
        current_atom_index = next_atom_index ;
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }
      else
      {
        itmp ++ ;
        
        current_atom_index = atoi( tmpString );
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }

    }

    natomDGroup = itmp ;
  
    printf("\nBased on index file , Donor fragment %s has %d atoms ...\n" , dGroupName , natomDGroup );
  
  }
  else
  {
    natomDGroup = 0 ;
  }

  
  
  // =====> D Fragment , Loading ...
  
  int * dGroupAtoms = calloc( natomDGroup , sizeof( int ) ) ;
  
  if( exD == YES )
  {
    rewind( pIndexFile ) ;
  
    fsearch( pIndexFile , dGroupName ) ;
  
    fskip( pIndexFile , 1 ) ;
  
    fscanf( pIndexFile , "%s" , tmpString );

    for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
    {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
      if( strcmp( tmpString , "-" ) == 0 )
      {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
      
        for( iatom = current_atom_index + 1 ; iatom <= next_atom_index ; iatom ++ )
        {
          *( dGroupAtoms + itmp + iatom - current_atom_index - 1 ) = iatom ;

        }
        
        itmp = itmp + ( next_atom_index - current_atom_index );
        
        current_atom_index = next_atom_index ;
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }
      else
      {   
        current_atom_index = atoi( tmpString );
      
        *( dGroupAtoms + itmp ) = current_atom_index ;
      
        itmp ++ ;
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }

    }
  
  }


  // =====> D Fragment , Generating Basis List ...
  
  int nbasisD = 0 ; 
  
  if( exD == YES )
  {
    for( iatom = 0 ; iatom < natomDGroup ; iatom ++ )
    {
      atomID = *( dGroupAtoms + iatom ) - 1 ; // C-Labeling ... 
  
      howManyFn = *( basisMap + 5 * atomID + 0 ) ; // But startFrom and newStart are in human labeling...
      
      nbasisD = nbasisD + howManyFn ;

    }
  
  }
  
  int * dBasisList = calloc( nbasisD , sizeof( int ) ) ; izeros( nbasisD , 1 , dBasisList ) ;
  
  accu = 0 ;

  if( exD == YES )
  {
    for( iatom = 0 ; iatom < natomDGroup ; iatom ++ )
    {
      atomID = *( dGroupAtoms + iatom ) - 1 ; // C-Labeling ... 
      
      newStart = *( basisMap + 5 * atomID + 2 ) ; // Human-Labeling ...
      
      howManyFn = *( basisMap + 5 * atomID + 0 ) ;
      
      for( ibasis = 0 ; ibasis < howManyFn ; ibasis ++ )
      {
        *( dBasisList + accu ) = newStart + ibasis ;
      
        accu ++ ;
      
      }
      
      
    }
  
    if ( debuggingMode == YES )
    {
      debug = fopen( "donorBasisSet.deb" , "wb+" );
  
      ioutput( debug , nbasisD , 1 , dBasisList ) ;
  
      fclose( debug ) ;
    }



  
  }




  // =====> A Fragment , Pre-Loading 
  
  rewind( pIndexFile ) ;
 
  /*
  
  info = fsearch( pIndexFile , aGroupName ) ; 
  
  rewind( pIndexFile ) ;
 
  if( exA == NO )
  {
    if( info == 1 )  exA = YES ;    
  }
  else */if( exA == YES )
  {
    info = fsearch( pIndexFile , aGroupName ) ;
    
    if( info == 0 )
    {
      printf("\nYour index file does not contain entry with name [ %s ] \n\n" , aGroupName ) ;
      
      exit( 61 ) ;
    }

    rewind( pIndexFile ) ;
  
    info = fsearch( pIndexFile , aGroupName ) ;
    
    fskip( pIndexFile , 1 ) ;
  
    if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
    { 
      printf("\nYour index file has hit the bottom with no specification for [ %s ] \n\n" , aGroupName ) ;
    
      exit( 62 ) ;
    }
  
    //printf("\nFor this group , the 1st grabbed tmpstring is %s ... \n" , tmpstring );
    
    //while( strcmp( tmpstring , "[") != 0 )
   
    for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
    {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
      if( strcmp( tmpString , "-" ) == 0 )
      {
          //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
        
        itmp = itmp + ( next_atom_index - current_atom_index );
        
        current_atom_index = next_atom_index ;
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }
      else
      {
        itmp ++ ;
        
        current_atom_index = atoi( tmpString );
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }

    }

    natomAGroup = itmp ;
  
    printf("\nBased on index file , Acceptor fragment %s has %d atoms ...\n" , aGroupName , natomAGroup );
  
  }
  else
  {
    natomAGroup = 0 ;
  }

  
  
  // =====> A Fragment , Loading ...
  
  int * aGroupAtoms = calloc( natomAGroup , sizeof( int ) ) ;
  
  if( exA == YES )
  {
    rewind( pIndexFile ) ;
  
    fsearch( pIndexFile , aGroupName ) ;
  
    fskip( pIndexFile , 1 ) ;
  
    fscanf( pIndexFile , "%s" , tmpString );

    for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
    {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
      if( strcmp( tmpString , "-" ) == 0 )
      {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
      
        for( iatom = current_atom_index + 1 ; iatom <= next_atom_index ; iatom ++ )
        {
          *( aGroupAtoms + itmp + iatom - current_atom_index - 1 ) = iatom ;

        }
        
        itmp = itmp + ( next_atom_index - current_atom_index );
        
        current_atom_index = next_atom_index ;
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }
      else
      {   
        current_atom_index = atoi( tmpString );
      
        *( aGroupAtoms + itmp ) = current_atom_index ;
      
        itmp ++ ;
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }

    }
  
  }



  // =====> A Fragment , Generating Basis List ...
  
  int nbasisA = 0 ; 
  
  if( exA == YES )
  {
    for( iatom = 0 ; iatom < natomAGroup ; iatom ++ )
    {
      atomID = *( aGroupAtoms + iatom ) - 1 ; // C-Labeling ... 
  
      howManyFn = *( basisMap + 5 * atomID + 0 ) ; // But startFrom and newStart are in human labeling...
      
      nbasisA = nbasisA + howManyFn ;

    }
  
  }
  
  int * aBasisList = calloc( nbasisA , sizeof( int ) ) ; izeros( nbasisA , 1 , aBasisList ) ;
  
  accu = 0 ;

  if( exA == YES )
  {
    for( iatom = 0 ; iatom < natomAGroup ; iatom ++ )
    {
      atomID = *( aGroupAtoms + iatom ) - 1 ; // C-Labeling ... 
      
      newStart = *( basisMap + 5 * atomID + 2 ) ; // Human-Labeling ...
      
      howManyFn = *( basisMap + 5 * atomID + 0 ) ;
      
      for( ibasis = 0 ; ibasis < howManyFn ; ibasis ++ )
      {
        *( aBasisList + accu ) = newStart + ibasis ;
      
        accu ++ ;
      
      }
      
      
    }
  
    if ( debuggingMode == YES )
    {
      debug = fopen( "acceptorBasisSet.deb" , "wb+" );
  
      ioutput( debug , nbasisA , 1 , aBasisList ) ;
  
      fclose( debug ) ;
    }


  }


  
  
  // =====> 2nd D Fragment , Pre-Loading 
  
  rewind( pIndexFile ) ;
 
  /*
 
  info = fsearch( pIndexFile , dGroupName ) ;
  
  rewind( pIndexFile ) ;
 
  if( exD == NO )
  {
    if( info == 1 )  exD = YES ;    
  }
  else */
  if( exD2 == YES )
  {
    info = fsearch( pIndexFile , d2GroupName ) ;
    
    if( info == 0 )
    {
      printf("\nYour index file does not contain entry with name [ %s ] \n\n" , d2GroupName ) ;
      
      exit( 61 ) ;
    }

    rewind( pIndexFile ) ;
  
    info = fsearch( pIndexFile , d2GroupName ) ;
    
    fskip( pIndexFile , 1 ) ;
  
    if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
    { 
      printf("\nYour index file has hit the bottom with no specification for [ %s ] \n\n" , d2GroupName ) ;
    
      exit( 62 ) ;
    }
  
    //printf("\nFor this group , the 1st grabbed tmpstring is %s ... \n" , tmpstring );
    
    //while( strcmp( tmpstring , "[") != 0 )
   
    for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
    {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
      if( strcmp( tmpString , "-" ) == 0 )
      {
          //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
        
        itmp = itmp + ( next_atom_index - current_atom_index );
        
        current_atom_index = next_atom_index ;
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }
      else
      {
        itmp ++ ;
        
        current_atom_index = atoi( tmpString );
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }

    }

    natomD2Group = itmp ;
  
    printf("\nBased on index file , 2nd Donor fragment %s has %d atoms ...\n" , d2GroupName , natomD2Group );
  
  }
  else
  {
    natomD2Group = 0 ;
  }

  
  
  // =====> 2nd D Fragment , Loading ...
  
  int * d2GroupAtoms = calloc( natomD2Group , sizeof( int ) ) ;
  
  if( exD2 == YES )
  {
    rewind( pIndexFile ) ;
  
    fsearch( pIndexFile , d2GroupName ) ;
  
    fskip( pIndexFile , 1 ) ;
  
    fscanf( pIndexFile , "%s" , tmpString );

    for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
    {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
      if( strcmp( tmpString , "-" ) == 0 )
      {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
      
        for( iatom = current_atom_index + 1 ; iatom <= next_atom_index ; iatom ++ )
        {
          *( d2GroupAtoms + itmp + iatom - current_atom_index - 1 ) = iatom ;

        }
        
        itmp = itmp + ( next_atom_index - current_atom_index );
        
        current_atom_index = next_atom_index ;
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }
      else
      {   
        current_atom_index = atoi( tmpString );
      
        *( d2GroupAtoms + itmp ) = current_atom_index ;
      
        itmp ++ ;
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }

    }
  
  }


  // =====> 2nd D Fragment , Generating Basis List ...
  
  int nbasisD2 = 0 ; 
  
  if( exD2 == YES )
  {
    for( iatom = 0 ; iatom < natomD2Group ; iatom ++ )
    {
      atomID = *( d2GroupAtoms + iatom ) - 1 ; // C-Labeling ... 
  
      howManyFn = *( basisMap + 5 * atomID + 0 ) ; // But startFrom and newStart are in human labeling...
      
      nbasisD2 = nbasisD2 + howManyFn ;

    }
  
  }
  
  int * d2BasisList = calloc( nbasisD2 , sizeof( int ) ) ; izeros( nbasisD2 , 1 , d2BasisList ) ;
  
  accu = 0 ;

  if( exD2 == YES )
  {
    for( iatom = 0 ; iatom < natomD2Group ; iatom ++ )
    {
      atomID = *( d2GroupAtoms + iatom ) - 1 ; // C-Labeling ... 
      
      newStart = *( basisMap + 5 * atomID + 2 ) ; // Human-Labeling ...
      
      howManyFn = *( basisMap + 5 * atomID + 0 ) ;
      
      for( ibasis = 0 ; ibasis < howManyFn ; ibasis ++ )
      {
        *( d2BasisList + accu ) = newStart + ibasis ;
      
        accu ++ ;
      
      }
      
      
    }
  
    if ( debuggingMode == YES )
    {
      debug = fopen( "2ndDonorBasisSet.deb" , "wb+" );
  
      ioutput( debug , nbasisD2 , 1 , d2BasisList ) ;
  
      fclose( debug ) ;
    }



  
  }






  // =====> 2nd A Fragment , Pre-Loading 
  
  rewind( pIndexFile ) ;
 
  /*
  
  info = fsearch( pIndexFile , aGroupName ) ; 
  
  rewind( pIndexFile ) ;
 
  if( exA == NO )
  {
    if( info == 1 )  exA = YES ;    
  }
  else */
  if( exA2 == YES )
  {
    info = fsearch( pIndexFile , a2GroupName ) ;
    
    if( info == 0 )
    {
      printf("\nYour index file does not contain entry with name [ %s ] \n\n" , a2GroupName ) ;
      
      exit( 61 ) ;
    }

    rewind( pIndexFile ) ;
  
    info = fsearch( pIndexFile , a2GroupName ) ;
    
    fskip( pIndexFile , 1 ) ;
  
    if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
    { 
      printf("\nYour index file has hit the bottom with no specification for [ %s ] \n\n" , a2GroupName ) ;
    
      exit( 62 ) ;
    }
  
    //printf("\nFor this group , the 1st grabbed tmpstring is %s ... \n" , tmpstring );
    
    //while( strcmp( tmpstring , "[") != 0 )
   
    for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
    {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
      if( strcmp( tmpString , "-" ) == 0 )
      {
          //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
        
        itmp = itmp + ( next_atom_index - current_atom_index );
        
        current_atom_index = next_atom_index ;
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }
      else
      {
        itmp ++ ;
        
        current_atom_index = atoi( tmpString );
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }

    }

    natomA2Group = itmp ;
  
    printf("\nBased on index file , Acceptor fragment %s has %d atoms ...\n" , aGroupName , natomAGroup );
  
  }
  else
  {
    natomA2Group = 0 ;
  }

  
  
  // =====> 2nd A Fragment , Loading ...
  
  int * a2GroupAtoms = calloc( natomA2Group , sizeof( int ) ) ;
  
  if( exA2 == YES )
  {
    rewind( pIndexFile ) ;
  
    fsearch( pIndexFile , a2GroupName ) ;
  
    fskip( pIndexFile , 1 ) ;
  
    fscanf( pIndexFile , "%s" , tmpString );

    for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
    {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
      if( strcmp( tmpString , "-" ) == 0 )
      {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
      
        for( iatom = current_atom_index + 1 ; iatom <= next_atom_index ; iatom ++ )
        {
          *( a2GroupAtoms + itmp + iatom - current_atom_index - 1 ) = iatom ;

        }
        
        itmp = itmp + ( next_atom_index - current_atom_index );
        
        current_atom_index = next_atom_index ;
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }
      else
      {   
        current_atom_index = atoi( tmpString );
      
        *( a2GroupAtoms + itmp ) = current_atom_index ;
      
        itmp ++ ;
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        if( ( ndxinfo = fscanf( pIndexFile , "%s" , tmpString ) ) == EOF )
        {
           break ;  
        }

      }

    }
  
  }



  // =====> 2nd A Fragment , Generating Basis List ...
  
  int nbasisA2 = 0 ; 
  
  if( exA2 == YES )
  {
    for( iatom = 0 ; iatom < natomA2Group ; iatom ++ )
    {
      atomID = *( a2GroupAtoms + iatom ) - 1 ; // C-Labeling ... 
  
      howManyFn = *( basisMap + 5 * atomID + 0 ) ; // But startFrom and newStart are in human labeling...
      
      nbasisA2 = nbasisA2 + howManyFn ;

    }
  
  }
  
  int * a2BasisList = calloc( nbasisA2 , sizeof( int ) ) ; izeros( nbasisA2 , 1 , a2BasisList ) ;
  
  accu = 0 ;

  if( exA2 == YES )
  {
    for( iatom = 0 ; iatom < natomA2Group ; iatom ++ )
    {
      atomID = *( a2GroupAtoms + iatom ) - 1 ; // C-Labeling ... 
      
      newStart = *( basisMap + 5 * atomID + 2 ) ; // Human-Labeling ...
      
      howManyFn = *( basisMap + 5 * atomID + 0 ) ;
      
      for( ibasis = 0 ; ibasis < howManyFn ; ibasis ++ )
      {
        *( a2BasisList + accu ) = newStart + ibasis ;
      
        accu ++ ;
      
      }
      
      
    }
  
    if ( debuggingMode == YES )
    {
      debug = fopen( "2ndAcceptorBasisSet.deb" , "wb+" );
  
      ioutput( debug , nbasisA2 , 1 , a2BasisList ) ;
  
      fclose( debug ) ;
    }


  }


  

  
  // -------------------------------> Reading the CI file ... <---------------------------------- //
  // ------ For : 1) Related Orbitals ;                       <---------------------------------- //
  // ------       2) Degenerate situation ;                   <---------------------------------- //
  // ------       3) # of electrons in each block ;  ???      <---------------------------------- //
  // -------------------------------------------------------------------------------------------- //
  
  int gsRef , leRef , ctRef ;
  
  int * gsRefEqv , * leRefEqv , * ctRefEqv ;
  int n_gsRefEqv , n_leRefEqv , n_ctRefEqv ;
  
  n_gsRefEqv = 0 ; 
  n_leRefEqv = 0 ; 
  n_ctRefEqv = 0 ;
  
  //int gsRefEqv , leRefEqv , ctRefEqv ;
  //gsRefEqv = leRefEqv = ctRefEqv = -1 ;
  
  //int * gsSnap , * leSnap , * ctSnap ;
  //int gsSnap , leSnap , ctSnap ;
  
  
  int degenerateOrNot = NO ;
  
  int throughBridgeOrNot = NO ; 
  
  int nRelatedOrbitals = 0 ;
  
  
  // =====> Getting the # of electrons ... 
  

  printf("\nAccording to the CI file ...\n\n");
  
  printf("\n[ %s ] Group has % 5d electrons , orbitals ranging from % 5d to % 5d , HOMO is the [ # % 5d ] orbital\n\n\n" , pGroupName , nelectronP , 1 , nbasisP , homoP ) ;

  printf("\n[ %s ] Group has % 5d electrons , orbitals ranging from % 5d to % 5d , HOMO is the [ # % 5d ] orbital\n\n\n" , qGroupName , nelectronQ , nbasisP + 1 , nbasisP + nbasisQ , homoQ ) ;
  
  if( exR == YES )
  {
    printf("\n[ %s ] Group has % 5d electrons , orbitals ranging from % 5d to % 5d , HOMO is the [ # % 5d ] orbital\n\n\n" , rGroupName , nelectronR , nbasisP + nbasisQ + 1 , nbasisP + nbasisQ + nbasisR , homoR ) ;
  }

  if( exX == YES )
  {
    printf("\n[ %s ] Group has % 5d electrons , orbitals ranging from % 5d to % 5d , HOMO is the [ # % 5d ] orbital\n\n\n" , xGroupName , nelectronX , nbasisP + nbasisQ + nbasisR + 1 , nbasisP + nbasisQ + nbasisR + nbasisX , homoX ) ;
  }

  if( exZ == YES )
  {
    printf("\n[ %s ] Group has % 5d electrons , orbitals ranging from % 5d to % 5d , HOMO is the [ # % 5d ] orbital\n\n\n" , zGroupName , nelectronZ , nbasisP + nbasisQ + nbasisR + nbasisX + 1 , nbasisP + nbasisQ + nbasisR + nbasisX + nbasisZ , homoZ ) ;
  }
  
  
  // =====> Dealing with degeneracy ... 
  
  rewind( prefCIListFile ) ;
  
  

  
  if( ( info = fsearch( prefCIListFile , "Degenerate" ) ) == 1 )
  {
    degenerateOrNot = YES ;
  
  }
  
  rewind( prefCIListFile ) ;
  
  
  

  // =====> Pre-Loading ... GS , LE , CT ... 
  
  rewind( prefCIListFile ) ;
  
  printf("\nNow let's pre-load the .CI file and see whether enough orbital numbers are provided ... \n");
  
  iline = 1 ;
  
  iload = 0 ;
  
  info = fsearch( prefCIListFile , "orbitals" ) ;
  
  if( info == 0 )
  {
    printf("\nWrong CIList format. Related orbitals should be listed in [ orbitals ] entry ... \n") ;
    
    exit( 63 ) ;
  }
  
  fskip( prefCIListFile , 1 );

  while( ( info = freadline( buffer , MAXCHARINLINE , prefCIListFile , '!' ) ) != 0 )
  { 
    //printf("\n//--------------> WORKING ON NO. %d LINE ... <-------------//\n" , iline );
    
    //info = freadline( buffer , MAXCHARINLINE , pinputfile , ';' );
    
    //printf("\nNow we are at : %s\n" , buffer );
    
    //printf("\n INFO is %d ...\n" , info );
    
    blank_signal = stellblank( buffer ) ;
    
    if( blank_signal == 0 )
    {
      printf("\nNo.%d line is a blank line ... Moving on ...\n" , iline ) ;
    }
    else if( blank_signal == 1 )
    {  
      //printf("\nNo.%d line is NOT a blank line ... loading ...\n" , iline );
      
      tmp_char = getfirst( buffer ) ;
      
      if( tmp_char == '!' )
      {
        printf("\nThis is a comment line ... So nothing will be loaded ...\n");
        
        //fskip( pinputfile , 1 );
      }
      else if( tmp_char == '[' )
      {
        printf("\nSTARTING OF NEW DIRECTIVE ... END READING ... \n\n");
        
        break ;
      }
      else
      {
        //printf("\nLine reads : %s ...\n" , buffer );
        
        //iload = iload + inLineWC( buffer ) ;
        
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

  nRelatedOrbitals = iload ;
  
  printf("\nThere are %d orbitals defined as related in CI file ...\n" , nRelatedOrbitals );

  if( nRelatedOrbitals > 3 )
  {
    printf("\nWARNING : There are %d orbitals mentioned in CI file while we only needed 3 ...\n" , nRelatedOrbitals ) ;
    
    printf("\nOnly first 3 orbitals are taken ... \n") ;

  }
  else if( nRelatedOrbitals < 3 )
  {
    printf("\nAt least 3 orbitals need to be defined in CIList file as \"related\" , but you have only %d ...\n" , nRelatedOrbitals ) ;
    
    exit( 137 ) ;
  }
  


  // =====> Pre-Loading ... Bridge Orbitals ... 
  
  int nbridge = 0 ; int * bridgeOrbitals ; double * bridgeOrbitalsEnergy ; 
  
  rewind( prefCIListFile ) ;
  
  printf("\nNow let's pre-load the .CI file to check out the bridge orbitals ... \n");
  
  iline = 1 ;
  
  iload = 0 ;
  
  info = fsearch( prefCIListFile , "bridges" ) ;
  
  if( info == 0 )
  {
    printf("\nIt appears that we do not need to do through-bridge calculations ...  \n") ;
    
    throughBridgeOrNot = NO ;
    
    nbridge = 0 ;
    
  }
  else
  {
    throughBridgeOrNot = YES ;
  
    fskip( prefCIListFile , 1 );

    while( ( info = freadline( buffer , MAXCHARINLINE , prefCIListFile , '!' ) ) != 0 )
    { 
      //printf("\n//--------------> WORKING ON NO. %d LINE ... <-------------//\n" , iline );
    
      //info = freadline( buffer , MAXCHARINLINE , pinputfile , ';' );
    
      //printf("\nNow we are at : %s\n" , buffer );
      
      //printf("\n INFO is %d ...\n" , info );
    
      blank_signal = stellblank( buffer ) ;
    
      if( blank_signal == 0 )
      {
        printf("\nNo.%d line is a blank line ... Moving on ...\n" , iline ) ;
      }
      else if( blank_signal == 1 )
      {  
        //printf("\nNo.%d line is NOT a blank line ... loading ...\n" , iline );
      
        tmp_char = getfirst( buffer ) ;
      
        if( tmp_char == '!' )
        {
          printf("\nThis is a comment line ... So nothing will be loaded ...\n");
        
          //fskip( pinputfile , 1 );
        }
        else if( tmp_char == '[' )
        {
          printf("\nSTARTING OF NEW DIRECTIVE ... END READING ... \n\n");
        
          break ;
        }
        else
        {
          //printf("\nLine reads : %s ...\n" , buffer );
        
          iload = iload + inLineWC( buffer ) ;
        
          //iload ++ ;
        
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

    
    if( iload == 0 )
    {
      throughBridgeOrNot = NO ; 
      
      nbridge = 0 ; 
    }
    else if( iload > 0 )
    {
      nbridge = iload ;
  
      printf("\nThere are %d orbitals defined as Bridge-Orbital in CI file ...\n" , nbridge );
    
      bridgeOrbitals = ( int * ) calloc( nbridge , sizeof( int ) ) ; 
    
      izeros( nbridge , 1 , bridgeOrbitals ) ;
      
      bridgeOrbitalsEnergy = ( double * ) calloc( nbridge , sizeof( double ) ) ;
      
      dzeros( nbridge , 1 , bridgeOrbitalsEnergy ) ;
      
      
    }
  
  
  }

  // =====> Loading ... GS , LE , CT ... 
  
  printf("\nNow let's actually load the .CI file  ... \n");
  
  rewind( prefCIListFile ) ;
  
  fsearch( prefCIListFile , "orbitals" ) ;
  
  if( ( info = fsearch( prefCIListFile , "GS" ) ) == 1 )
  {
    fscanf( prefCIListFile , "%s" , tmpString ) ;
    
    if( strcmp( tmpString , "=" ) == 0 )
    {
      fscanf( prefCIListFile , "%s" , cache ) ;
      
      gsRef = atof( cache ) ;
    }
    else
    {
      printf("\nWrong format of CIList file ... No Trailing \"=\"\n") ;
      
      exit( 67 ) ;
    }
  
  }
  else
  {
    printf("\nThere is no \"GS\" in your [ orbitals ] entry ...\n") ;
    
    exit( 69 ) ;
  }
  
  
  
  
  rewind( prefCIListFile ) ;
  
  fsearch( prefCIListFile , "orbitals" ) ;
  
  if( ( info = fsearch( prefCIListFile , "LE" ) ) == 1 )
  {
    fscanf( prefCIListFile , "%s" , tmpString ) ;
    
    if( strcmp( tmpString , "=") == 0 )
    {
      fscanf( prefCIListFile , "%s" , cache ) ;
      
      leRef = atof( cache ) ;
    }
    else
    {
      printf("\nWrong format of CIList file ... No Trailing \"=\"\n") ;
      
      exit( 67 ) ;
    }
  
  }
  else
  {
    printf("\nThere is no \"LE\" in your [ orbitals ] entry ...\n") ;
    
    exit( 69 ) ;
  }
  
  
  
  
  rewind( prefCIListFile ) ;
  
  fsearch( prefCIListFile , "orbitals" ) ;
  
  if( ( info = fsearch( prefCIListFile , "CT" ) ) == 1 )
  {
    fscanf( prefCIListFile , "%s" , tmpString ) ;
    
    if( strcmp( tmpString , "=") == 0 )
    {
      fscanf( prefCIListFile , "%s" , cache ) ;
      
      ctRef = atof( cache ) ;
    }
    else
    {
      printf("\nWrong format of CIList file ... No Trailing \"=\"\n") ;
      
      exit( 67 ) ;
    }
  
  }
  else
  {
    printf("\nThere is no \"CT\" in your [ orbitals ] entry ...\n") ;
    
    exit( 69 ) ;
  }
  
  
  printf("\nSo ... [ GS = %d ] , [ LE = %d ] , [ CT = %d ] ... \n" , gsRef , leRef , ctRef ) ;
  
  
  // =====> Loading ... Degeneracy ... 
  
  
  rewind( prefCIListFile ) ;
  
  iload = 0 ; iline = 1 ;
  
  int ieqv = 0 ;
  
  int nwords = 0 ;
  
  if( degenerateOrNot == YES )
  { 
    itmp = preLoadEntry( prefCIListFile , "Degenerate" , '!' ) ;
    
    rewind( prefCIListFile ) ;
    
    fsearch( prefCIListFile , "Degenerate" ) ;
    
    fskip( prefCIListFile , 1 ) ;
    

    while( ( info = freadline( buffer , MAXCHARINLINE , prefCIListFile , '!' ) ) != 0 )
    { 
      printf("\n//--------------> WORKING ON NO. %d LINE ... <-------------//\n" , iline );
    
      printf("\nLine Reads : %s \n" , buffer ) ;
      
      blank_signal = stellblank( buffer ) ;
    
      if( blank_signal == 0 )
      {
        printf("\nNo.%d line is a blank line ... Moving on ...\n" , iline ) ;
      
        continue ;
      }
      else if( blank_signal == 1 )
      {    
        //printf("\nNo.%d line is NOT a blank line ... loading ...\n" , iline );
      
        tmp_char = getfirst( buffer ) ;
      
        if( tmp_char == '!' )
        {
          //printf("\nThis is a comment line ... So nothing will be loaded ...\n");
        
          //fskip( pinputfile , 1 );
        
          continue ;
        }
        else if( tmp_char == '[' )
        {
          printf("\nSTARTING OF NEW DIRECTIVE ... END READING ... \n\n");
        
          break ;
        }
        else
        {
          //printf("\nLine Reads : %s \n" , buffer ) ;
          
          nwords = inLineWC( buffer ) ;
          
          if( nwords >= 2 )
          {
            strpickword( buffer , 1 , cache ) ;
          
            itmp2 = atoi( cache ) ;
            
            //strpickword( buffer , 3 , cache ) ;
          
            if( itmp2 == gsRef )
            {
              //printf("\nLine Reads : %s \n" , buffer ) ;
              printf("\nDegeneracy Info on GS :\n\n") ;
              
              n_gsRefEqv = nwords - 1 ;
              
              gsRefEqv = calloc( nwords - 1 , sizeof( int ) ) ;
              
              *( gsRefEqv + 0 ) = gsRef ;
              
              for( ieqv = 2 ; ieqv < nwords ; ieqv ++  )
              {
                strpickword( buffer , ieqv + 1 , cache ) ;
                
                *( gsRefEqv + ieqv - 1 ) = atoi( cache ) ;  
                
                printf("\n # %d orbital is equivalent to [ gsRef = %d ] ...\n" , *( gsRefEqv + ieqv - 1 ) , gsRef ) ;   
              }
              
                 
            }
            else if( itmp2 == leRef )
            {
              //printf("\nLine Reads : %s \n" , buffer ) ;
              printf("\nDegeneracy Info on LE :\n\n") ;
              
              n_leRefEqv = nwords - 1 ;
              
              leRefEqv = calloc( nwords - 1 , sizeof( int ) ) ;
              
              *( leRefEqv + 0 ) = leRef ;
              
              for( ieqv = 2 ; ieqv < nwords ; ieqv ++  )
              {
                strpickword( buffer , ieqv + 1 , cache ) ;
                
                *( leRefEqv + ieqv - 1 ) = atoi( cache ) ;  
                
                printf("\n # %d orbital is equivalent to [ leRef = %d ] ...\n" , *( leRefEqv + ieqv - 1 ) , leRef ) ;   
              }
              
            }
            else if( itmp2 == ctRef )
            {
              //printf("\nLine Reads : %s \n" , buffer ) ;
              printf("\nDegeneracy Info on CT : \n\n") ;
              
              n_ctRefEqv = nwords - 1 ;
              
              ctRefEqv = calloc( nwords - 1 , sizeof( int ) ) ;
              
              *( ctRefEqv + 0 ) = ctRef ;
              
              for( ieqv = 2 ; ieqv < nwords ; ieqv ++  )
              {
                strpickword( buffer , ieqv + 1 , cache ) ;
                
                *( ctRefEqv + ieqv - 1 ) = atoi( cache ) ;  
                
                printf("\n # %d orbital is equivalent to [ ctRef = %d ] ...\n" , *( ctRefEqv + ieqv - 1 ) , ctRef ) ;   
              }

            }
            else
            {
              //printf("\nLine Reads : %s \n" , buffer ) ;
              
              printf("\nWe actually do not need other degeneracy information ...\n") ;
            }
          
            iload ++ ;
          
          }

          
        
        }
      
      //printf("\n%s\n" , buffer );
      }
      else
      {
        printf("\nSomething is wrong with the reading file part ...\n");
      
        exit( 92 );
      }
    
      if( iload == itmp ) 
      {
        printf( "\nen toto [ %d ] degeneracy entry needed ... All degeneracy information loaded ...\n\n" , itmp ) ;
      
        break ;
      }
    
    
      iline ++ ;

    }
  
  
  
  
  
  }
  else
  {
    n_leRefEqv = 1 ;
              
    leRefEqv = calloc( 1 , sizeof( int ) ) ;
    
    *leRefEqv = leRef ;
    
    n_ctRefEqv = 1 ;
              
    ctRefEqv = calloc( 1 , sizeof( int ) ) ;
    
    *ctRefEqv = ctRef ;
    
    n_gsRefEqv = 1 ;
              
    gsRefEqv = calloc( 1 , sizeof( int ) ) ;
    
    *gsRefEqv = gsRef ;
    
  }
  
      
      
  // =====> Loading ... Bridge Orbitals ...      
  
  
  rewind( prefCIListFile ) ;
  
  if( throughBridgeOrNot == YES )
  {
    iline = 1 ;
  
    iload = 0 ;
    
    itmp = 0 ; 
  
    info = fsearch( prefCIListFile , "bridges" ) ;
  
    fskip( prefCIListFile , 1 );
    
    printf( "\n\n[ Bridge-Orbitals ]\n" ) ;
    
    for( iload = 0 ; iload < nbridge ; iload ++ )
    {
      fscanf( prefCIListFile , "%s" , cache ) ;
      
      *( bridgeOrbitals + iload ) = atoi( cache ) ;
      
      printf( "\n%d\n" , *( bridgeOrbitals + iload ) ) ;
    
    }
    
    rewind( prefCIListFile ) ;
    
    /*
    while( ( info = freadline( buffer , MAXCHARINLINE , prefCIListFile , '!' ) ) != 0 )
    { 
      //printf("\n//--------------> WORKING ON NO. %d LINE ... <-------------//\n" , iline );
    
      //info = freadline( buffer , MAXCHARINLINE , pinputfile , ';' );
    
      //printf("\nNow we are at : %s\n" , buffer );
      
      //printf("\n INFO is %d ...\n" , info );
    
      blank_signal = stellblank( buffer ) ;
    
      if( blank_signal == 0 )
      {
        printf("\nNo.%d line is a blank line ... Moving on ...\n" , iline ) ;
      }
      else if( blank_signal == 1 )
      {  
        //printf("\nNo.%d line is NOT a blank line ... loading ...\n" , iline );
      
        tmp_char = getfirst( buffer ) ;
      
        if( tmp_char == '!' )
        {
          printf("\nThis is a comment line ... So nothing will be loaded ...\n");
        
          //fskip( pinputfile , 1 );
        }
        else if( tmp_char == '[' )
        {
          printf("\nSTARTING OF NEW DIRECTIVE ... END READING ... \n\n");
        
          break ;
        }
        else
        {
          //printf("\nLine reads : %s ...\n" , buffer );
        
          itmp = inLineWC( buffer ) ;
        
          //iload ++ ;
        
        }
        //printf("\n%s\n" , buffer );
      }
      else
      {
        printf("\nSomething is wrong with the reading file part ...\n");
      
        exit(1);
      }
    
      iline ++ ;
  
    }
    */
  
  
  }

  
  // ------------------> Calculating MO overlap and Mapping ... <---------------------- //
  
  double * overlap = calloc( nbasis * nbasis , sizeof( double ) ) ;
  
  dzeros( nbasis , nbasis , overlap ) ;
  
  int * moMapping = calloc( nbasis , sizeof( int ) ) ; izeros( nbasis , 1 , moMapping ) ;
  
  int * smallOverlapFlag = calloc( nbasis , sizeof( int ) ) ; izeros( nbasis , 1 , smallOverlapFlag ) ;
  
  double dominateOverlap = 0.000 ;
    
  
  if( noRefSituation == NO )
  {
    dgemm_( "N" , "T" , &nbasis , &nbasis , &nbasis , &done , refMO , &nbasis , MO , &nbasis , &dzero , overlap , &nbasis ) ;
    
    dtranspose( nbasis , overlap , overlap ) ;
    
    
  if ( debuggingMode == YES )
  {
    debug = fopen( "overlap.deb" , "wb+" ) ;
    
    doutput( debug , nbasis , nbasis , overlap ) ;
    
    fclose( debug ) ;

  }    
    
    // ---> P Block : Occupied Orbitals ...
    
    for( ibasis = 0 ; ibasis < homoP ; ibasis ++ )
    {
      //printf("\nWe are @ %d ...\n" , ibasis ) ;
      
      dominateOverlap = dmaxAbs( homoP , overlap + ibasis * nbasis ) ;
      
      itmp = dmaxAbsID( homoP , overlap + ibasis * nbasis ) ;
    
      if( dominateOverlap <= 0.65 )
      {
        printf("\n===> WARNING : Maximum overlap for Ref-MO # %d is only % 10.6f <===\n" , ibasis + 1 , dominateOverlap ) ;
      
        *( smallOverlapFlag + ibasis ) = YES ;
      }
      
      /*
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      */
      
    
      *( moMapping + ibasis ) = itmp ; 
      
      
  
    }  printf( "\nDone for P Block - Occupied Orbitals ... \n\n") ;
    
    
    
    // ---> P Block : Virtual Orbitals ...
    
    for( ibasis = homoP ; ibasis < nbasisP ; ibasis ++ )
    {
      //printf("\nWe are @ %d ...\n" , ibasis ) ;
      
      dominateOverlap = dmaxAbs( nbasisP - homoP , overlap + ibasis * nbasis + homoP ) ;
      
      itmp = dmaxAbsID( nbasisP - homoP , overlap + ibasis * nbasis + homoP ) ;
    
      if( dominateOverlap <= 0.65 )
      {
        printf("\n===> WARNING : Maximum overlap for Ref-MO # %d is only % 10.6f <===\n" , ibasis + 1 , dominateOverlap ) ;
        
        *( smallOverlapFlag + ibasis ) = YES ;
        
      }
      
      /*
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      */
      
    
      *( moMapping + ibasis ) = homoP + itmp ; 
      
      
  
    }  printf( "\nDone for P Block - Virtual Orbitals ... \n\n") ;
    
    
    // ---> Q Block : Occupied Orbitals ...
    
    for( ibasis = nbasisP ; ibasis < nbasisP + homoQ ; ibasis ++ )
    {
      //printf("\nWe are @ %d ...\n" , ibasis ) ;
      
      dominateOverlap = dmaxAbs( homoQ , overlap + ibasis * nbasis + nbasisP ) ;
      
      itmp = dmaxAbsID( homoQ , overlap + ibasis * nbasis + nbasisP ) ;
    
      if( dominateOverlap <= 0.65 )
      {
        printf("\n===> WARNING : Maximum overlap for Ref-MO # %d is only % 10.6f <===\n" , ibasis + 1 , dominateOverlap ) ;
      
        *( smallOverlapFlag + ibasis ) = YES ;
      }
      
      /*
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      */
      
    
      *( moMapping + ibasis ) = nbasisP + itmp ; 
      
      
  
    }  printf( "\nDone for Q Block - Occupied Orbitals ... \n\n") ;
    
    // ---> Q Block : Virtual Orbitals ...
    
    for( ibasis = nbasisP + homoQ ; ibasis < nbasisP + nbasisQ ; ibasis ++ )
    {
      //printf("\nWe are @ %d ...\n" , ibasis ) ;
      
      dominateOverlap = dmaxAbs( nbasisQ - homoQ , overlap + ibasis * nbasis + nbasisP + homoQ ) ;
      
      itmp = dmaxAbsID( nbasisQ - homoQ , overlap + ibasis * nbasis + nbasisP + homoQ ) ;
    
      if( dominateOverlap <= 0.65 )
      {
        printf("\n===> WARNING : Maximum overlap for Ref-MO # %d is only % 10.6f <===\n" , ibasis + 1 , dominateOverlap ) ;
        
        *( smallOverlapFlag + ibasis ) = YES ;
        
      }
      
      /*
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      */
      
    
      *( moMapping + ibasis ) = nbasisP + homoQ + itmp ; 
      
      
  
    }  printf( "\nDone for Q Block - Virtual Orbitals ... \n\n") ;
    
    
    // ---> R Block : Occupied Orbitals ...
    
    for( ibasis = nbasisP + nbasisQ  ; ibasis < nbasisP + nbasisQ + homoR ; ibasis ++ )
    {
      //printf("\nWe are @ %d ...\n" , ibasis ) ;
      
      dominateOverlap = dmaxAbs( homoR , overlap + ibasis * nbasis + nbasisP + nbasisQ ) ;
    
      itmp = dmaxAbsID( homoR , overlap + ibasis * nbasis + nbasisP + nbasisQ ) ;
    
      if( dominateOverlap <= 0.65 )
      {
        printf("\n===> WARNING : Maximum overlap for Ref-MO # %d is only % 10.6f <===\n" , ibasis + 1 , dominateOverlap ) ;
      
        *( smallOverlapFlag + ibasis ) = YES ;
      
      }
      
      /*
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      */
      
    
      *( moMapping + ibasis ) = nbasisP + nbasisQ + itmp ; 
      
      
  
    }  printf( "\nDone for R Block - Occupied Orbitals ... \n\n") ;
  
  
    // ---> R Block : Virtual Orbitals ...
    
    for( ibasis = nbasisP + nbasisQ + homoR ; ibasis < nbasisP + nbasisQ + nbasisR ; ibasis ++ )
    {
      //printf("\nWe are @ %d ...\n" , ibasis ) ;
      
      dominateOverlap = dmaxAbs( nbasisR - homoR , overlap + ibasis * nbasis + nbasisP + nbasisQ + homoR ) ;
      
      itmp = dmaxAbsID( nbasisR - homoR , overlap + ibasis * nbasis + nbasisP + nbasisQ + homoR ) ;
    
      if( dominateOverlap <= 0.65 )
      {
        printf("\n===> WARNING : Maximum overlap for Ref-MO # %d is only % 10.6f <===\n" , ibasis + 1 , dominateOverlap ) ;
        
        *( smallOverlapFlag + ibasis ) = YES ;
        
      }
      
      /*
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      */
      
    
      *( moMapping + ibasis ) = nbasisP + nbasisQ + homoR + itmp ; 
      
      
  
    }  printf( "\nDone for R Block - Virtual Orbitals ... \n\n") ;
  
  
  
    // ---> X Block : Occupied Orbitals ...
    
    for( ibasis = nbasisP + nbasisQ + nbasisR ; ibasis < nbasisP + nbasisQ + nbasisR + homoX ; ibasis ++ )
    {
      //printf("\nWe are @ %d ...\n" , ibasis ) ;
      
      dominateOverlap = dmaxAbs( homoX , overlap + ibasis * nbasis + nbasisP + nbasisQ + nbasisR ) ;
    
      itmp = dmaxAbsID( homoX , overlap + ibasis * nbasis + nbasisP + nbasisQ + nbasisR ) ;
    
      if( dominateOverlap <= 0.65 )
      {
        printf("\n===> WARNING : Maximum overlap for Ref-MO # %d is only % 10.6f <===\n" , ibasis + 1 , dominateOverlap ) ;
      
        *( smallOverlapFlag + ibasis ) = YES ;
      
      }
      
      /*
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      */
      
    
      *( moMapping + ibasis ) = nbasisP + nbasisQ + nbasisR + itmp ; 
      
      
  
    }  printf( "\nDone for X Block - Occupied Orbitals ... \n\n") ;
  
  
    // ---> X Block : Virtual Orbitals ...
    
    for( ibasis = nbasisP + nbasisQ + nbasisR + homoX ; ibasis < nbasisP + nbasisQ + nbasisR + nbasisX ; ibasis ++ )
    {
      //printf("\nWe are @ %d ...\n" , ibasis ) ;
      
      dominateOverlap = dmaxAbs( nbasisX - homoX , overlap + ibasis * nbasis + nbasisP + nbasisQ + nbasisR + homoX ) ;
      
      itmp = dmaxAbsID( nbasisX - homoX , overlap + ibasis * nbasis + nbasisP + nbasisQ + nbasisR + homoX ) ;
    
      if( dominateOverlap <= 0.65 )
      {
        printf("\n===> WARNING : Maximum overlap for Ref-MO # %d is only % 10.6f <===\n" , ibasis + 1 , dominateOverlap ) ;
        
        *( smallOverlapFlag + ibasis ) = YES ;
        
      }
      
      /*
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      */
      
    
      *( moMapping + ibasis ) = nbasisP + nbasisQ + nbasisR + homoX + itmp ; 
      
      
  
    }  printf( "\nDone for X Block - Virtual Orbitals ... \n\n") ;
  
    // ---> Z Block : Occupied Orbitals ...
    
    printf( "\n[ NBasis ] nbasis is %d [ NBasis ] nbasisX is %d [ NBasis ] \n\n" , nbasis , nbasisX ) ;
    
    for( ibasis = nbasisP + nbasisQ + nbasisR + nbasisX ; ibasis < nbasisP + nbasisQ + nbasisR + nbasisX + homoZ ; ibasis ++ )
    {
      //printf("\nWe are @ %d ...\n" , ibasis ) ; 
      
      dominateOverlap = dmaxAbs( homoZ , overlap + ibasis * nbasis + nbasisP + nbasisQ + nbasisR + nbasisX ) ;
    
      itmp = dmaxAbsID( homoZ , overlap + ibasis * nbasis + nbasisP + nbasisQ + nbasisR + nbasisX ) ;
    
      if( dominateOverlap <= 0.65 )
      {
        printf("\n===> WARNING : Maximum overlap for Ref-MO # %d is only % 10.6f <===\n" , ibasis + 1 , dominateOverlap ) ;
      
        *( smallOverlapFlag + ibasis ) = YES ;
      
      }
      
      /*
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      */
      
    
      *( moMapping + ibasis ) = nbasisP + nbasisQ + nbasisR + nbasisX + itmp ; 
      
      
  
    }  printf( "\nDone for Z Block - Occupied Orbitals ... \n\n") ;
  
  
    // ---> Z Block : Virtual Orbitals ...
    
    for( ibasis = nbasisP + nbasisQ + nbasisR + nbasisX + homoZ ; ibasis < nbasisP + nbasisQ + nbasisR + nbasisX + nbasisZ ; ibasis ++ )
    {
      //printf("\nWe are @ %d ...\n" , ibasis ) ;
      
      dominateOverlap = dmaxAbs( nbasisZ - homoZ , overlap + ibasis * nbasis + nbasisP + nbasisQ + nbasisR + nbasisX + homoZ ) ;
      
      itmp = dmaxAbsID( nbasisZ - homoZ , overlap + ibasis * nbasis + nbasisP + nbasisQ + nbasisR + nbasisX + homoZ ) ;
    
      if( dominateOverlap <= 0.65 )
      {
        printf("\n===> WARNING : Maximum overlap for Ref-MO # %d is only % 10.6f <===\n" , ibasis + 1 , dominateOverlap ) ;
        
        *( smallOverlapFlag + ibasis ) = YES ;
        
      }
      
      /*
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      */
      
    
      *( moMapping + ibasis ) = nbasisP + nbasisQ + nbasisR + nbasisX + homoZ + itmp ; 
      
      
  
    }  printf( "\nDone for Z Block - Virtual Orbitals ... \n\n") ;


  
  }
  else if( noRefSituation == YES )
  {
    printf("\n===> WARNING : NO REFERENCE PROTOCOL INVOKED <===\n\n" ) ;
    
    for( ibasis = 0 ; ibasis < nbasis ; ibasis ++ )
    {
      *( moMapping + ibasis ) = ibasis ; 
    }

  }
  else
  {
    printf("\nSomething is wrong about the \"noRefSituation\" ... Please check !\n") ;
    
    exit( 127 ) ;
  
  }
  
  
  if ( debuggingMode == YES )
  {

    debug = fopen( "moMapping.deb" , "wb+" ) ; // ATTENTION : ALL C-LABEL HERE ...
  
    ioutput( debug , nbasis , 1 , moMapping ) ;
  
    fclose( debug ) ;
  
  }
  
  

  
  // -----------------> Presenting the HDAs ... <-------------//
  
 
  double * tmp_hda = calloc( 8 , sizeof( double ) ) ;
  
  dzeros( 8 , 1 , tmp_hda ) ;
  
  int nLECTtransition  , nGSCTtransition , nLEorbital , nGSorbital , nCTorbital ;
  
  nLECTtransition = 1 ; nGSCTtransition = 1 ; nLEorbital = 1 ; nGSorbital = 1 ; nCTorbital = 1 ;
  
  int nOrbitalInvolved = 0 ;
  
  
  itmp = 1 ; itmp2 = 1 ;
  
  int itmp3 = 1 ;
  
  //int nOrbitalInvolved = 0 ;
  
  int nTransitionInvolved = 2 ;
  
  
  nLEorbital = n_leRefEqv ;
  
  nGSorbital = n_gsRefEqv ;
  
  nCTorbital = n_ctRefEqv ;
  
  
  nLECTtransition = nLEorbital * nCTorbital ;
  
  nGSCTtransition = nGSorbital * nCTorbital ;
  
  
  
  nTransitionInvolved = nLECTtransition + nGSCTtransition ;
  
  nOrbitalInvolved = nLEorbital + nGSorbital + nCTorbital ;
  

  double * hda = calloc( nTransitionInvolved , sizeof( double ) ) ;
  
  dzeros( nTransitionInvolved , 1 , hda ) ;
  
  
  double * involvedOrbitalEnergies = calloc( nOrbitalInvolved , sizeof( double ) ) ;
  
  dzeros( nOrbitalInvolved , 1 , involvedOrbitalEnergies ) ;
  
  
  double * energyGapLECT = calloc( nLECTtransition , sizeof( double ) ) ; 
  
  dzeros( nLECTtransition , 1 , energyGapLECT ) ;
  
  double * energyGapGSCT = calloc( nGSCTtransition , sizeof( double ) ) ; 
  
  dzeros( nGSCTtransition , 1 , energyGapGSCT ) ;

  
  int * overlap_warning = calloc( nOrbitalInvolved , sizeof( int ) ) ;
  
  izeros( nOrbitalInvolved , 1 , overlap_warning ) ;




  double presentedLECTcoupling , presentedGSCTcoupling ;
  
  int ile , igs , ict , itransition , iorbital ;
  
  int cur_leRef , cur_ctRef , cur_gsRef ;
  
  int lectIndicator = 0 ; int gsctIndicator = 0 ; // C-Labeling for now ... 
  
  double energyLE , energyGS , energyCT ;
  
  double minimumLECTenergyGap = 0.00 ; double minimumGSCTenergyGap = 0.00 ; 

  
  int presentedFlag1 , presentedFlag2 , presentedFlag3 , presentedFlag4 ;
  

  // ---> Orbital Energy Info ...
  
  iorbital = 0 ; 
  
  for( ile = 0 ; ile < nLEorbital ; ile ++ )
  {
    cur_leRef = *( leRefEqv + ile ) ;
    
    *( involvedOrbitalEnergies + iorbital ) = *( reorgEnergy + ( *( moMapping + cur_leRef - 1 ) ) ) ;
    
    iorbital ++ ;
  }

  for( ict = 0 ; ict < nCTorbital ; ict ++ )
  {
    cur_ctRef = *( ctRefEqv + ict ) ;
    
    *( involvedOrbitalEnergies + iorbital ) = *( reorgEnergy + ( *( moMapping + cur_ctRef - 1 ) ) ) ;
    
    iorbital ++ ;
  }

  for( igs = 0 ; igs < nGSorbital ; igs ++ )
  {
    cur_gsRef = *( gsRefEqv + igs ) ;
    
    *( involvedOrbitalEnergies + iorbital ) = *( reorgEnergy + ( *( moMapping + cur_gsRef - 1 ) ) ) ;
    
    iorbital ++ ;
  }
  
  if( iorbital != nOrbitalInvolved )
  {
    printf("\nSomething is wrong when enumerating orbital energies ...\n\n");
    
    exit( 93 ) ;
  }



  // ---> Let's issue the warnings ...
  
  /*

  if( leRefEqv != -1 )
  {
    *( overlap_warning + 0 ) = *( smallOverlapFlag + ( *( moMapping + leRef - 1 ) ) ) ;
    
    *( overlap_warning + 1 ) = *( smallOverlapFlag + ( *( moMapping + leRefEqv - 1 ) ) ) ;
  }
  else
  {
    *( overlap_warning + 0 ) = *( smallOverlapFlag + ( *( moMapping + leRef - 1 ) ) ) ;
  }
  
  if( gsRefEqv != -1 )
  {
    *( overlap_warning + itmp + 0 ) = *( smallOverlapFlag + ( *( moMapping + gsRef - 1 ) ) ) ;
    
    *( overlap_warning + itmp + 1 ) = *( smallOverlapFlag + ( *( moMapping + gsRefEqv - 1 ) ) ) ;
  }
  else
  {
    *( overlap_warning + itmp + 0 ) = *( smallOverlapFlag + ( *( moMapping + gsRef - 1 ) ) ) ;
  }
  
  
  if( ctRefEqv != -1 )
  {
    *( overlap_warning + itmp + itmp2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    *( overlap_warning + itmp + itmp2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
  }
  else
  {
    *( overlap_warning + itmp + itmp2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + ctRef - 1 ) ) ) ;
  }
  
  */




  // -------> We Need to figure out the LE and CT that have the smallest gap ! Oct. 2nd !
  // -------> But we will output all the transitions ... ! Feb. 25th, 2015 !
  // -------> There will be flags indicating which transition has the smallest gap ...
  
    
  //	ATTENTION : LE-CT FIRST  //
  
  itransition = 0 ;
  
  for( ict = 0 ; ict < n_ctRefEqv ; ict ++ )
  {
    cur_ctRef = *( ctRefEqv + ict ) ;
    
    energyCT = *( reorgEnergy + ( *( moMapping + cur_ctRef - 1 ) ) ) ;
    
    for( ile = 0 ; ile < n_leRefEqv ; ile ++ )
    {
      cur_leRef = *( leRefEqv + ile ) ;
      
      energyLE = *( reorgEnergy + ( *( moMapping + cur_leRef - 1 ) ) ) ;
      
      *( hda + itransition ) = *( finalFock + ( *( moMapping + cur_leRef - 1 ) ) * nbasis + ( *( moMapping + cur_ctRef - 1 ) )  ) ;
      
      printf("\n[%d] <-> [%d] : HDA = % 12.6E \n" , cur_leRef , cur_ctRef , *( hda + itransition ) ) ;
      
      *( energyGapLECT + itransition ) = fabs( energyCT - energyLE ) ;
      
      printf("\n[%d] <-> [%d] : Energy Gap = % 12.6E \n" , cur_leRef , cur_ctRef , *( energyGapLECT + itransition ) ) ;
      
      itransition ++ ;
      
    }
  
  }
  
  lectIndicator = dminID( nLECTtransition , energyGapLECT ) ;
  
  minimumLECTenergyGap = fabs( *( energyGapLECT + lectIndicator ) ) ;
  
  presentedLECTcoupling = *( hda + lectIndicator ) ;
  
  

  //	ATTENTION : FOLLOWED BY GS-CT  //
  
  //itransition = 0 ;
  
  for( ict = 0 ; ict < n_ctRefEqv ; ict ++ )
  {
    cur_ctRef = *( ctRefEqv + ict ) ;
    
    energyCT = *( reorgEnergy + ( *( moMapping + cur_ctRef - 1 ) ) ) ;
    
    for( igs = 0 ; igs < n_gsRefEqv ; igs ++ )
    {
      cur_gsRef = *( gsRefEqv + igs ) ;
      
      energyGS = *( reorgEnergy + ( *( moMapping + cur_gsRef - 1 ) ) ) ;
      
      *( hda + itransition ) = *( finalFock + ( *( moMapping + cur_gsRef - 1 ) ) * nbasis + ( *( moMapping + cur_ctRef - 1 ) )  ) ;
      
      printf("\n[%d] <-> [%d] : HDA = % 12.6E \n" , cur_gsRef , cur_ctRef , *( hda + itransition ) ) ;
      
      *( energyGapGSCT + itransition - nLECTtransition ) = fabs( energyCT - energyGS ) ;
      
      printf("\n[%d] <-> [%d] : Energy Gap = % 12.6E \n" , cur_gsRef , cur_ctRef , *( energyGapGSCT + itransition- nLECTtransition ) ) ;
      
      itransition ++ ;
      
    }
  
  }
  
  gsctIndicator = dminID( nGSCTtransition , energyGapGSCT ) ;
  
  minimumGSCTenergyGap = fabs( *( energyGapGSCT + gsctIndicator ) ) ;
  
  presentedGSCTcoupling = *( hda + gsctIndicator + nLECTtransition ) ;
  
  

 


  // -----------------> Recording the Bridge Orbital Energies ... and nearest-neighbour coupling ... <-------------//

  int firstBridge , lastBridge , currentBridge , nextBridge , ibridge ; 
  
  double * vGSBridge , * vLEBridge , * vCTBridge , * vBB ; 
  
  int nGSBridge , nLEBridge , nCTBridge , nBB , iBB ;   
  
  double * nearestNeighbourCoupling ; 
  
  int nNearestNeighbourCoupling ;
  
  
  if( throughBridgeOrNot == YES )
  {
    firstBridge = *( bridgeOrbitals + 0 ) ;
  
    lastBridge = *( bridgeOrbitals + nbridge - 1 ) ;  itmp = 0 ;
  
  
    for( ibridge = 0 ; ibridge < nbridge ; ibridge ++ )
    {
      currentBridge = *( bridgeOrbitals + ibridge ) ; 
    
      *( bridgeOrbitalsEnergy + ibridge ) = *( reorgEnergy + ( *( moMapping + currentBridge - 1 ) ) ) ; 
  
    }

    
    // ---> < LE | H | Bridge >
    
    nLEBridge = nbridge * n_leRefEqv ;
    
    nCTBridge = nbridge * n_ctRefEqv ;
    
    nGSBridge = nbridge * n_gsRefEqv ;
    
    vLEBridge = calloc( nLEBridge , sizeof( double ) ) ;
    
    vCTBridge = calloc( nCTBridge , sizeof( double ) ) ;
    
    vGSBridge = calloc( nGSBridge , sizeof( double ) ) ;
    
    
    
    itransition = 0 ; 
    
    for( ile = 0 ; ile < n_leRefEqv ; ile ++ )
    {
      cur_leRef = *( leRefEqv + ile ) ;
      
      for( ibridge = 0 ; ibridge < nbridge ; ibridge ++ )
      {
        currentBridge = *( bridgeOrbitals + ibridge ) ; 
        
        *( vLEBridge + itransition ) = *( finalFock + ( *( moMapping + cur_leRef - 1 ) ) * nbasis + ( *( moMapping + currentBridge - 1 ) ) ) ;
      
        itransition ++ ;
        
      }
    
    }
    
    itransition = 0 ; 
    
    for( igs = 0 ; igs < n_gsRefEqv ; igs ++ )
    {
      cur_gsRef = *( gsRefEqv + igs ) ;
      
      for( ibridge = 0 ; ibridge < nbridge ; ibridge ++ )
      {
        currentBridge = *( bridgeOrbitals + ibridge ) ; 
        
        *( vGSBridge + itransition ) = *( finalFock + ( *( moMapping + cur_gsRef - 1 ) ) * nbasis + ( *( moMapping + currentBridge - 1 ) ) ) ;
      
        itransition ++ ;
        
      }
    
    }
    
    itransition = 0 ; 
    
    for( ict = 0 ; ict < n_ctRefEqv ; ict ++ )
    {
      cur_ctRef = *( ctRefEqv + ict ) ;
      
      for( ibridge = 0 ; ibridge < nbridge ; ibridge ++ )
      {
        currentBridge = *( bridgeOrbitals + ibridge ) ; 
        
        *( vCTBridge + itransition ) = *( finalFock + ( *( moMapping + cur_ctRef - 1 ) ) * nbasis + ( *( moMapping + currentBridge - 1 ) ) ) ;
      
        itransition ++ ;
        
      }
    
    }
    
    
    
  
    // ---> < Bridge_{k-1} | H | Bridge_k >

    nBB = nbridge - 1 ;
  
    vBB = ( double * ) calloc( nBB , sizeof( double ) ) ; 
    
    dzeros( nBB , 1 , vBB ) ;
   
    for( ibridge = 0 ; ibridge < nbridge - 1 ; ibridge ++ )
    {
      iBB = ibridge ;
    
      currentBridge = ( *( bridgeOrbitals + ibridge ) ) ; 
    
      nextBridge = ( *( bridgeOrbitals + ibridge + 1 ) ) ; 
    
      *( vBB + iBB ) = *( finalFock + ( *( moMapping + currentBridge - 1 ) ) * nbasis + ( *( moMapping + nextBridge - 1 ) ) ) ;

    }



    // ---> Put together ...

    //nNearestNeighbourCoupling = nGSBridge + nLEBridge + nBridgeCT + nBB ;
    
    //nearestNeighbourCoupling = ( double * ) calloc( nNearestNeighbourCoupling , sizeof( double ) ) ;
    



  }
  
  
  
  // ------------------> Calculating the MO Localization ... <---------------------- //
  
  double * transReOrgMO = calloc( nbasis * nbasis , sizeof( double ) ) ;
  
  dzeros( nbasis , nbasis , transReOrgMO ) ;
  
  dtranspose( nbasis , MO , transReOrgMO ) ; // Transpose the Re-Organized MO for convenience ...
  
  int le_loc , ct_loc , gs_loc ;
  
  
  //double leLocalization = 0.00 , leEqvLocalization = 0.00 , ctLocalization = 0.00 , ctEqvLocalization = 0.00 , gsLocalization = 0.00 , gsEqvLocalization = 0.00 ;
  
  double * leLocalization , * ctLocalization , * gsLocalization ;
  
  leLocalization = calloc( n_leRefEqv * ( exD + exD2 ) , sizeof( double ) ) ; dzeros( n_leRefEqv , ( exD + exD2 ) , leLocalization ) ;
  
  ctLocalization = calloc( n_ctRefEqv * ( exA + exA2 ) , sizeof( double ) ) ; dzeros( n_ctRefEqv , ( exA + exA2 ) , ctLocalization ) ;
  
  gsLocalization = calloc( n_gsRefEqv * ( exD + exD2 ) , sizeof( double ) ) ; dzeros( n_gsRefEqv , ( exD + exD2 ) , gsLocalization ) ;
  
  
  /* double calcSquaredLocalization( double * orbital , int nBasis , int * basisList , int nComponent ) */
  
  printf("\nRelevant MO Localization Circumstances ...\n\n" ) ;
  
  //--> LE <--//

  printf("\n");
  
  if( exD == YES )
  {
    for( ile = 0 ; ile < n_leRefEqv ; ile ++ )
    {
      cur_leRef = *( leRefEqv + ile ) ;
      
      orbitalID = *( moMapping + cur_leRef - 1 ) ;
      
      *( leLocalization + ile ) = calcSquaredLocalization( ( transReOrgMO + nbasis * orbitalID ) , nbasis , dBasisList , nbasisD ) ;
      
      printf("[%d] on [ Donor : %s ] : % 10.6E\t\n" , *( leRefEqv + ile ) , dGroupName , *( leLocalization + ile ) ) ;
    
    }

  }
  
  if( exD2 == YES )
  {
    for( ile = 0 ; ile < n_leRefEqv ; ile ++ )
    {
      cur_leRef = *( leRefEqv + ile ) ;
      
      orbitalID = *( moMapping + cur_leRef - 1 ) ;
      
      *( leLocalization + n_leRefEqv + ile ) = calcSquaredLocalization( ( transReOrgMO + nbasis * orbitalID ) , nbasis , d2BasisList , nbasisD2 ) ;
      
      printf("[%d] on [ 2nd Donor : %s ] : % 10.6E\t\n" , *( leRefEqv + ile ) , d2GroupName , *( leLocalization + n_leRefEqv + ile ) ) ;
    
    }

  }

  

  //--> CT <--//
  
  if( exA == YES )
  {
    for( ict = 0 ; ict < n_ctRefEqv ; ict ++ )
    {
      cur_ctRef = *( ctRefEqv + ict ) ;
      
      orbitalID = *( moMapping + cur_ctRef - 1 ) ;
      
      *( ctLocalization + ict ) = calcSquaredLocalization( ( transReOrgMO + nbasis * orbitalID ) , nbasis , aBasisList , nbasisA ) ;
      
      printf("[%d] on [ Acceptor : %s ] : % 10.6E\t\n" , *( ctRefEqv + ict ) , aGroupName , *( ctLocalization + ict ) ) ;
    
    }

  }


  if( exA2 == YES )
  {
    for( ict = 0 ; ict < n_ctRefEqv ; ict ++ )
    {
      cur_ctRef = *( ctRefEqv + ict ) ;
      
      orbitalID = *( moMapping + cur_ctRef - 1 ) ;
      
      *( ctLocalization + n_ctRefEqv + ict ) = calcSquaredLocalization( ( transReOrgMO + nbasis * orbitalID ) , nbasis , a2BasisList , nbasisA2 ) ;
      
      printf("[%d] on [ 2nd Acceptor : %s ] : % 10.6E\t\n" , *( ctRefEqv + ict ) , a2GroupName , *( ctLocalization + n_ctRefEqv + ict ) ) ;
    
    }

  }
  
  //--> GS <--//
  
  if( exD == YES )
  {
    for( igs = 0 ; igs < n_gsRefEqv ; igs ++ )
    {
      cur_gsRef = *( gsRefEqv + igs ) ;
      
      orbitalID = *( moMapping + cur_gsRef - 1 ) ;
      
      *( gsLocalization + igs ) = calcSquaredLocalization( ( transReOrgMO + nbasis * orbitalID ) , nbasis , dBasisList , nbasisD ) ;
      
      printf("[%d] on [ Donor : %s ] : % 10.6E\t\n" , *( gsRefEqv + igs ) , dGroupName , *( gsLocalization + igs ) ) ;
    
    }

  }
  
  
  if( exD2 == YES )
  {
    for( igs = 0 ; igs < n_gsRefEqv ; igs ++ )
    {
      cur_gsRef = *( gsRefEqv + igs ) ;
      
      orbitalID = *( moMapping + cur_gsRef - 1 ) ;
      
      *( gsLocalization + n_gsRefEqv + igs ) = calcSquaredLocalization( ( transReOrgMO + nbasis * orbitalID ) , nbasis , d2BasisList , nbasisD2 ) ;
      
      printf("[%d] on [ 2nd Donor : %s ] : % 10.6E\t\n" , *( gsRefEqv + igs ) , d2GroupName , *( gsLocalization + n_gsRefEqv + igs ) ) ;
    
    }

  }
  
  
  
  printf("\n\n\n");
    
  
  double presentLECT_localization = 0.00 , presentGSCT_localization = 0.00 ;
  
  itmp = 0 ; itmp2 = 0 ;
  
  
  itmp2 = dmaxAbsID( n_leRefEqv , leLocalization ) ;
  
  itmp = itmp2 ; // % n_leRefEqv ; 
  
  le_loc = *( leRefEqv + itmp ) ;
  
  
  itmp2 = dmaxAbsID( n_ctRefEqv , ctLocalization ) ;
  
  itmp = itmp2 ; // % n_ctRefEqv ; 
  
  ct_loc = *( ctRefEqv + itmp ) ;
  
  
  itmp2 = dmaxAbsID( n_gsRefEqv , gsLocalization ) ;
  
  itmp = itmp2 ; // % n_gsRefEqv ; 
  
  gs_loc = *( gsRefEqv + itmp ) ;
  
  
  presentLECT_localization = *( finalFock + ( *( moMapping + le_loc - 1 ) ) * nbasis + ( *( moMapping + ct_loc - 1 ) ) ) ;
  
  presentGSCT_localization = *( finalFock + ( *( moMapping + gs_loc - 1 ) ) * nbasis + ( *( moMapping + ct_loc - 1 ) ) ) ;
  
  
  
    
  
  
  // ------------------> Outputing ... <---------------------- //
  
  
  // ATTENTION : Here we output fabs( presented*coupling ) instead of themselves ...
  
  
  poutHDAFile = fopen( outputHDAName , "wb+" ) ;
  
  fprintf( poutHDAFile , "\n" ) ;
  
  

  for( iload = 0 ; iload < nTransitionInvolved ; iload ++ )
  {
    fprintf( poutHDAFile , "% 10.6E\t" , fabs( *( hda + iload ) ) ) ;
  
  }
  
  fprintf( poutHDAFile , "% 3d\t% 3d\t" , lectIndicator + 1 , gsctIndicator + 1 ) ; // easier for matlab to use, here is human labeling
  // These indicators are based on minimum energy gap ...ZZZZZZZZZ
  
  
  for( iload = 0 ; iload < nOrbitalInvolved ; iload ++ )
  {
    fprintf( poutHDAFile , "% 10.6E\t" , *( involvedOrbitalEnergies + iload ) ) ;
  }
  
  
  if( throughBridgeOrNot == YES )
  {
    for( ibridge = 0 ; ibridge < nbridge ; ibridge ++ )
    {
      fprintf( poutHDAFile , "% 10.6E\t" , *( bridgeOrbitalsEnergy + ibridge ) ) ;
    
    }
  
    
    for( itmp = 0 ; itmp < nLEBridge ; itmp ++ )
    {
      fprintf( poutHDAFile , "% 10.6E\t" , *( vLEBridge + itmp ) ) ;
    }
  
    for( itmp = 0 ; itmp < nCTBridge ; itmp ++ )
    {
      fprintf( poutHDAFile , "% 10.6E\t" , *( vCTBridge + itmp ) ) ;
    }
  
    for( itmp = 0 ; itmp < nGSBridge ; itmp ++ )
    {
      fprintf( poutHDAFile , "% 10.6E\t" , *( vGSBridge + itmp ) ) ;
    }
  
  
    for( ibridge = 0 ; ibridge < nbridge - 1 ; ibridge ++ )
    {
      iBB = ibridge ;
      
      fprintf( poutHDAFile , "% 10.6E\t" , *( vBB + iBB ) ) ;
      
    }
  
  
  
  
  }
  
  
  fprintf( poutHDAFile , "% 10.6E\t% 10.6E\t% 10.6E\t% 10.6E\t" , fabs( presentedLECTcoupling ) , fabs( presentedGSCTcoupling ) , fabs( presentLECT_localization ) , fabs( presentGSCT_localization ) ) ;
  fprintf( poutHDAFile , "%d\t%d\t%d\t" , le_loc , ct_loc , gs_loc ) ;
  
  // MO Localization Information ... 
  
  for( itmp = 0 ; itmp < n_leRefEqv * ( exD + exD2 ) ; itmp ++ )
  {
    fprintf( poutHDAFile , "% 10.6E\t" , *( leLocalization + itmp ) ) ;
  }
  
  for( itmp = 0 ; itmp < n_ctRefEqv * ( exA + exA2 ) ; itmp ++ )
  {
    fprintf( poutHDAFile , "% 10.6E\t" , *( ctLocalization + itmp ) ) ;
  }
  
  for( itmp = 0 ; itmp < n_gsRefEqv * ( exD + exD2 ) ; itmp ++ )
  {
    fprintf( poutHDAFile , "% 10.6E\t" , *( gsLocalization + itmp ) ) ;
  }
  
  

  // The Requested Coupling :
  
  cur_ctRef = *( ctRefEqv + 0 ) ;
  double presentGSCT_requested = *( finalFock + ( *( moMapping + gs_loc - 1 ) ) * nbasis + ( *( moMapping + cur_ctRef - 1 ) ) ) ;
  fprintf( poutHDAFile , "% 10.6E\t\n" , presentGSCT_requested ) ;
  
 

  return( 0 ) ;






}








double calcLocalization( double * orbital , int nBasis , int * basisList , int nComponent )
{
  if( nBasis < nComponent )
  {
    printf("\nERROR : There are only %d basis in total in the orbital , you are asking %d components ...\n\n" , nBasis , nComponent ) ;
    
    exit( 79 ) ;
  }
  
  int iComponent , basisID ;
  
  double sumComponent = 0.00 ;
  
  
  for( iComponent = 0 ; iComponent < nComponent ; iComponent ++ )
  {
    basisID = ( *( basisList + iComponent ) - 1 ) ;
    
    sumComponent = sumComponent + ( *( orbital + basisID ) ) ;
  
  
  }
  
  return( sumComponent ) ;



}






double calcSquaredLocalization( double * orbital , int nBasis , int * basisList , int nComponent )
{
  if( nBasis < nComponent )
  {
    printf("\nERROR : There are only %d basis in total in the orbital , you are asking %d components ...\n\n" , nBasis , nComponent ) ;
    
    exit( 79 ) ;
  }
  
  int iComponent , basisID ;
  
  double sumComponent = 0.00 ;
  
  
  for( iComponent = 0 ; iComponent < nComponent ; iComponent ++ )
  {
    basisID = ( *( basisList + iComponent ) - 1 ) ;
    
    sumComponent = sumComponent + ( *( orbital + basisID ) ) * ( *( orbital + basisID ) ) ;
  
  
  }
  
  return( sumComponent ) ;



}





