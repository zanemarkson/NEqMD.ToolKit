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

  // =====> Utility Functions Declaration ...

  double calcLocalization( double * orbital , int nBasis , int * basisList , int nComponent ) ;
  
  double calcSquaredLocalization( double * orbital , int nBasis , int * basisList , int nComponent );


  char ** pcmd = argv ; 
  
  int icmd ;
  
  char refMOFileName[ 300 ] , refCIListName[ 300 ] , cndoLogFileName[ 300 ] , indexFileName[ 300 ] , outputHDAName[ 300 ] ;

  int len_cndoLogFileName ;
  
  FILE * prefMOFile , * prefCIListFile , * pcndoLogFile , * pIndexFile , * poutHDAFile ;
  
  double * refMO ; 
  
  int orbitalID ;
  
  int debuggingMode = NO ;

  char pGroupName[ 150 ] , qGroupName[ 150 ] , rGroupName[ 150 ] , xGroupName[ 150 ] , zGroupName[ 150 ];
  
  char dGroupName[ 150 ] , aGroupName[ 150 ] ;
  
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
      printf("* G_CNDO2BDIAG_D : Calculate B-Diag  HDA for Snapshots in MD Calcs.  *\n");
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
  
  strcpy( cndoLogFileName , "system.snap.log" ) ;
  
  strcpy( outputHDAName , "system.hda" ) ;
  
  
  strcpy( pGroupName , "Donor-Block" ) ;
  
  strcpy( qGroupName , "Acceptor-Block" ) ;
  
  strcpy( rGroupName , "Bridge-Block" ) ;
  
  strcpy( xGroupName , "Irrelevant-Block" ) ;
  
  strcpy( zGroupName , "Irrelevant-Block-2" ) ;
  
  
  // =====> Parsing cmd-line arguments ...
  
  int exR = NO , exX = NO , exZ = NO ;
  
  int exD = NO , exA = NO ;
  
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

	                 icmd = icmd + 2 ;
	                 
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
	    
	      case 'h' : printf("\nUsage:  %s [ -m 'input ref. MO Coefficient file name' ] [ -L 'Ref. CI list file name' ] [ -l 'Input cndo log file name' ]\n\
	                                      [ -o output file name ] [ -n GROMACS .ndx file name ] [ -g debugging mode (yes/y/YES/Yes or no/n/NO/No) ]\n\
	                                      [ -p Group P ] [ -q Group Q ] [ -r Group R ] [ -x Group X ] [ -z Group Z ]\n\
	                                      [ -D Name of Donor Fragment ] [ -A Name of Acceptor Fragment ]\n" , * argv ); 
	                 
                     printf("\n                Note : 1) when -m is specified as \"None\" or \"none\" or \"NONE\", no reference will be checked and orbital numbers provided in CI file will be directly used.\n\n");
                     printf("\n                Note : 2) Default group selection : [ -p \"Donor-Block\"] [ -q \"Acceptor-Block\"] [ -r \"Bridge-Block\"] [ -x \"Irrelevant-Block\"] [ -z \"Irrelevant-Block-2\"]\n\n\n") ; 
                     printf("\n                Note : 3) Default fragment names selection : [ -D \"Donor\"] [ -A \"Acceptor\"]\n\n\n") ; 
                     printf("\n                Note : 4) If one of the partition is not necessary, \"None/none/None\" has to be specified ...\n\n\n") ; 


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

  
  // ---> NCIStates
    
  rewind( pcndoLogFile ) ;


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
  
  
  
  // ---> Fock in AO
  
  int irow , icol , iblock , idMO ;
  
  int rowid , colid ;

  int moPerBlock = 15 ;
  
  int moLeftOver = nbasis % moPerBlock ;
  
  int nblock = ( nbasis - moLeftOver ) / moPerBlock ;
  
  lmo = nbasis * nbasis ;
  

  
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
  
  cndoFAO = ( double * ) calloc( lmo , sizeof( double ) ) ; dzeros( lmo , 1 , cndoFAO ) ;
  
  MOEnergy = ( double * ) calloc( nbasis , sizeof( double ) ) ; dzeros( nbasis , 1 , MOEnergy ) ;// in eV unit ...
  

  
  fsearch( pcndoLogFile , "Final" ) ;
  
  fskip( pcndoLogFile , 3 ) ;
  
  for( iblock = 0 ; iblock < nblock - 1 ; iblock ++ )
  {
    //printf("\n# %d Block ... \n" , iblock ) ;
    
    /* First 15 Lines */
    
    for( irow = 0 ; irow < moPerBlock ; irow ++ )
    {
      rowid = iblock * moPerBlock + irow ;
      
      //printf("\nReading # %d line of FAO ...\n" , rowid + 1 ) ;
      
      fscanf( pcndoLogFile , "%d" , &itmp ) ;
      
      if( rowid != itmp - 1 )
      {
        printf("\nSomething is wrong with reading the 15-triangle in %d ( Human label) block ...\n" , iblock + 1 ) ;
        
        printf("\nLog file says it is the No. %d line ... you are actually @ No. %d line ...\n" , itmp , rowid + 1 ) ;
        
        exit( 382 ) ;
      }
      
      fscanf( pcndoLogFile , "%s" , stmp ) ;  fscanf( pcndoLogFile , "%s" , stmp ) ;  fscanf( pcndoLogFile , "%s" , stmp ) ;
      
      for( icol = 0 ; icol < irow + 1 ; icol ++ )
      {
        colid = iblock * moPerBlock + icol ;
        
        fscanf( pcndoLogFile , "%s" , cache ) ;
        
        *( cndoFAO + nbasis * rowid + colid ) = atof( cache ) ;
        
      }
      
      //printf("\nThe last number of this line is %lf ... \n\n" , atof( cache ) ) ;
    

    }
    
    /* The Rest */
    
    for( irow = rowid + 1 ; irow < nbasis ; irow ++ )
    {
      fscanf( pcndoLogFile , "%d" , &itmp ) ;
      
      //printf("\nReading # %d line of FAO ...\n" , irow + 1 ) ;
      
      if( irow != itmp - 1 )
      {
        printf("\nSomething is wrong with reading the big-rectangle in %d ( Human label) block ...\n" , iblock + 1 ) ;
        
        exit( 382 ) ;
      }
      
      fscanf( pcndoLogFile , "%s" , stmp ) ;  fscanf( pcndoLogFile , "%s" , stmp ) ;  fscanf( pcndoLogFile , "%s" , stmp ) ;
      

      for( icol = iblock * moPerBlock ; icol < ( iblock + 1 ) * moPerBlock ; icol ++ )
      {
        fscanf( pcndoLogFile , "%s" , cache ) ;
        
        *( cndoFAO + nbasis * irow + icol ) = atof( cache ) ;
      }
    
      //printf("\nThe last number of this line is %lf ... " , atof( cache ) ) ;
    
    }
  
    
    fskip( pcndoLogFile , 4 ) ;
    
  }
  
  /* Now only a moLeftOver-by-moLeftOver triangle left */
  
  //fskip( pcndoLogFile , 4 ) ;
  
  //printf("\nTHE FINAL BLOCK\n" ) ;
  
  for( irow = 0 ; irow < moLeftOver ; irow ++ )
  {
    rowid = ( nblock - 1 ) * moPerBlock + irow ;
    
    fscanf( pcndoLogFile , "%d" , &itmp ) ;
    
    if( rowid != itmp - 1 )
    {
      printf("\nSomething is wrong with reading the 15-triangle in the final block ...\n" ) ;
      
      printf("\nLog file says it is the No. %d line ... you are actually @ No. %d line ...\n" , itmp , rowid + 1 ) ;
      
      exit( 382 ) ;
    }
    
    fscanf( pcndoLogFile , "%s" , stmp ) ;  fscanf( pcndoLogFile , "%s" , stmp ) ;  fscanf( pcndoLogFile , "%s" , stmp ) ;
  
    for( icol = 0 ; icol < irow + 1 ; icol ++ )
    {
      colid = ( nblock - 1 ) * moPerBlock + icol ;
      
      fscanf( pcndoLogFile , "%s" , cache ) ;
      
      *( cndoFAO + nbasis * rowid + colid ) = atof( cache ) ;
      
    }
  
    //printf("\nThe last number of this line is %lf ... " , atof( cache ) ) ;
  }
  
  //void dtri2sym( char uplo, int lda, double * tri, double * sym)
  
  //dtri2sym( 'L' , nbasis , fao , fao ) ;
  
  for( irow = 0 ; irow < nbasis ; irow ++ )
  {
    for( icol = irow ; icol < nbasis ; icol ++ )
    {
      *( cndoFAO + irow * nbasis + icol ) = *( cndoFAO + icol * nbasis + irow ) ;
    
    
    }
  
  
  }
  
  
  if ( debuggingMode == YES )
  {
    debug = fopen( "cndoFAO.deb" , "wb+" ) ;
  
    doutput( debug , nbasis , nbasis , cndoFAO ) ;
  
    fclose( debug ) ;
  }
  
  
  
  
  
  // ---> Basis function location information ... 
  
  int * basisMap = calloc( 5 * natom , sizeof( int ) ) ;
  
  izeros( 5 * natom , 1 , basisMap ) ;
  
  /*
    1st Column : N of Basis Function for Current Atom ;
    
    2nd Column : @ Input , "first bf" location ;
    
    3rd Column : After re-organization , "first bf" location ;
  
  
  */
  
  int iline = 0 ;
  
  int iload = 0 ;

  
  rewind( pcndoLogFile ) ;
  
  fsearch( pcndoLogFile , "---------coordinates--------" ) ;
  
  fsearch( pcndoLogFile , "assign" ) ;
  
  fskip( pcndoLogFile , 1 ) ;
  
  while( ( info = freadline( buffer , MAXCHARINLINE , pcndoLogFile , ';' ) ) != 0 )
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
        //printf( "\nWorking on No. %d atom ...\n" , EqAtomList[ iatom ].atomnumber );

        strpickword( buffer , 7 , cache ) ; *( basisMap + 5 * iload + 0 ) = atoi( cache ) ;
        
        //sscanf( pEqGRO , "%lf" , &EqAtomList[ iload ].cy ); //printf("\n Cy is %lf ...\t" , EqAtomList[ iatom ].cy);
        
        strpickword( buffer , 8 , cache ) ; *( basisMap + 5 * iload + 1 ) = atoi( cache ) ;
        
        
        
        
        strpickword( buffer , 9 , cache ) ; *( basisMap + 5 * iload + 4 ) = atoi( cache ) ;
        
        
        iload ++ ;
        
      }
      
      //printf("\n%s\n" , buffer );
    }
    else
    {
      printf("\nSomething is wrong with the reading file part ...\n");
      
      exit(1);
    }
    
    if( iload == natom ) break ;
    
    
    iline ++ ;

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
  
  
  
  // =====> D Fragment , Pre-Loading 
  
  rewind( pIndexFile ) ;
 
  info = fsearch( pIndexFile , dGroupName ) ;
  
  rewind( pIndexFile ) ;
 
  if( exD == NO )
  {
    if( info == 1 )  exD = YES ;    
  }
  else if( exD == YES )
  {
    info = fsearch( pIndexFile , dGroupName ) ;
    
    if( info == 0 )
    {
      printf("\nYour index file does not contain entry with name [ %s ] \n\n" , dGroupName ) ;
      
      exit( 61 ) ;
    }

  }
  
  
  if( exD == YES )
  {
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
  

  
  }




  // =====> A Fragment , Pre-Loading 
  
  rewind( pIndexFile ) ;
 
  info = fsearch( pIndexFile , aGroupName ) ; 
  
  rewind( pIndexFile ) ;
 
  if( exA == NO )
  {
    if( info == 1 )  exA = YES ;    
  }
  else if( exA == YES )
  {
    info = fsearch( pIndexFile , aGroupName ) ;
    
    if( info == 0 )
    {
      printf("\nYour index file does not contain entry with name [ %s ] \n\n" , aGroupName ) ;
      
      exit( 61 ) ;
    }

  }
  
  
  if( exA == YES )
  {
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
  
  
  }


  
  
  





  
  // -------------------------------> Reading the CI file ... <---------------------------------- //
  // ------ For : 1) Related Orbitals ;                       <---------------------------------- //
  // ------       2) Degenerate situation ;                   <---------------------------------- //
  // ------       3) # of electrons in each block ;  ???      <---------------------------------- //
  // -------------------------------------------------------------------------------------------- //
  
  int gsRef , leRef , ctRef ;
  
  int gsRefEqv , leRefEqv , ctRefEqv ;
  
  gsRefEqv = leRefEqv = ctRefEqv = -1 ;
  
  int gsSnap , leSnap , ctSnap ;
  
  
  int degenerateOrNot = NO ;
  
  int throughBridgeOrNot = NO ; 
  
  int nRelatedOrbitals = 0 ;
  
  
  // =====> Getting the # of electrons ... 
  
  // int nelectronP = 0 , nelectronQ = 0 , nelectronR = 0 , nelectronX = 0 , nelectronZ = 0 ;
  
  // int lumoP = 0 , homoP = 0 , lumoQ = 0 , homoQ = 0 , lumoR = 0 , homoR = 0 , lumoX = 0 , homoX = 0 , lumoZ = 0 , homoZ = 0 ; 
  

  
  //rewind( prefCIListFile ) ;
  
  //info = fsearch( prefCIListFile , "electrons" ) ;
  
  
  /*
  
  if( info == 0 )
  {
    printf("\nThe number of electrons must be provided for each block.\n\n") ;
    
    exit( 135 ) ;
  
  }
  
  */
  

  // ---> P Group 
  
  /*
  
  if( ( info = fsearch( prefCIListFile , pGroupName ) ) == 0 )
  {
    printf("\nThe number of electrons of group [ %s ] must be provided ... \n\n" , pGroupName ) ;
    
    exit( 135 ) ;
  
  }
  
  
  
  
  fscanf( prefCIListFile , "%s" , tmpString ) ;
  
  fscanf( prefCIListFile , "%s" , cache ) ;
  
  //nelectronP = atoi( cache ) ;
  
  
  
  homoP = ( nelectronP + ( nelectronP % 2 ) ) / 2 ; // Human-Label
  
  lumoP = homoP + 1 ; // Human-Label
  
  
  
  
  rewind( prefCIListFile ) ;
  
  info = fsearch( prefCIListFile , "electrons" ) ;

  */
  

  // ---> Q Group 
  
  /*
  
  if( ( info = fsearch( prefCIListFile , qGroupName ) ) == 0 )
  {
    printf("\nThe number of electrons of group [ %s ] must be provided ... \n\n" , qGroupName ) ;
    
    exit( 135 ) ;
  
  }
  
  
  
  fscanf( prefCIListFile , "%s" , tmpString ) ;
  
  fscanf( prefCIListFile , "%s" , cache ) ;
  
  //nelectronQ = atoi( cache ) ;
  
  //homoQ = nbasisP + ( nelectron_q + ( nelectron_q % 2 ) ) / 2 ; // Human-Label
  
  homoQ = ( nelectronQ + ( nelectronQ % 2 ) ) / 2 ; // Human-Label
  
  lumoQ = homoQ + 1 ; // Human-Label
  
  
  rewind( prefCIListFile ) ;
  
  info = fsearch( prefCIListFile , "electrons" ) ;

  */

  // ---> R Group  
  
  
  /*
  if( exR == YES )
  {
    
    
    if( ( info = fsearch( prefCIListFile , rGroupName ) ) == 0 )
    {
      printf("\nThe number of electrons of group [ %s ] must be provided ... \n\n" , rGroupName ) ;
    
      exit( 135 ) ;
  
    }
    
    
    
  
    fscanf( prefCIListFile , "%s" , tmpString ) ;
    
    fscanf( prefCIListFile , "%s" , cache ) ;
    
    //nelectronR = atoi( cache ) ;
    
    //homoR = nbasisP + nbasisQ + ( nelectron_r + ( nelectron_r % 2 ) ) / 2 ; // Human-Label
    
    homoR = ( nelectronR + ( nelectronR % 2 ) ) / 2 ; // Human-Label
  
    lumoR = homoR + 1 ; // Human-Label
  
  
    rewind( prefCIListFile ) ;
  
    info = fsearch( prefCIListFile , "electrons" ) ;
  
  
  }
  */

  // ---> X Group
  
  
  /*
  if( exX == YES )
  {
   
    
    if( ( info = fsearch( prefCIListFile , xGroupName ) ) == 0 )
    {
      printf("\nThe number of electrons of group [ %s ] must be provided ... \n\n" , xGroupName ) ;
    
      exit( 135 ) ;
  
    }
    
    
    
  
    fscanf( prefCIListFile , "%s" , tmpString ) ;
    
    fscanf( prefCIListFile , "%s" , cache ) ;
    
    //nelectronX = atoi( cache ) ;
    
    //homoX = nbasisP + nbasisQ + nbasisR + ( nelectron_x + ( nelectron_x % 2 ) ) / 2 ; // Human-Label
    
    homoX = ( nelectronX + ( nelectronX % 2 ) ) / 2 ; // Human-Label
  
    lumoX = homoX + 1 ; // Human-Label
  
  
    rewind( prefCIListFile ) ;
  
    // info = fsearch( prefCIListFile , "electrons" ) ;
  
  
  }
  */
  

  // ---> Z Group
  
  
  /*
  if( exZ == YES )
  {
    
    
    if( ( info = fsearch( prefCIListFile , zGroupName ) ) == 0 )
    {
      printf("\nThe number of electrons of group [ %s ] must be provided ... \n\n" , zGroupName ) ;
    
      exit( 135 ) ;
  
    }
    
    
  
    fscanf( prefCIListFile , "%s" , tmpString ) ;
    
    fscanf( prefCIListFile , "%s" , cache ) ;
    
    //nelectronZ = atoi( cache ) ;
    
    //homoX = nbasisP + nbasisQ + nbasisR + ( nelectron_x + ( nelectron_x % 2 ) ) / 2 ; // Human-Label
    
    homoZ = ( nelectronZ + ( nelectronZ % 2 ) ) / 2 ; // Human-Label
  
    lumoZ = homoZ + 1 ; // Human-Label
  
  
    rewind( prefCIListFile ) ;
  
    // info = fsearch( prefCIListFile , "electrons" ) ;
  
  
  }
  */




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
    printf("\nAt lease 3 orbitals need to be defined in CIList file as \"related\" , but you have only %d ...\n" , nRelatedOrbitals ) ;
    
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
          
          if( inLineWC( buffer ) == 3 )
          {
            strpickword( buffer , 1 , cache ) ;
          
            itmp2 = atoi( cache ) ;
            
            strpickword( buffer , 3 , cache ) ;
          
            if( itmp2 == gsRef )
            {
              //printf("\nLine Reads : %s \n" , buffer ) ;
              
              gsRefEqv = atoi( cache ) ;  
              
              printf("\n # %d orbital is equivalent to [ gsRef = %d ] ...\n" , gsRefEqv , gsRef ) ;      
            }
            else if( itmp2 == leRef )
            {
              //printf("\nLine Reads : %s \n" , buffer ) ;
              
              leRefEqv = atoi( cache ) ;
              
              printf("\n # %d orbital is equivalent to [ leRef = %d ] ...\n" , leRefEqv , leRef ) ; 
            }
            else if( itmp2 == ctRef )
            {
              //printf("\nLine Reads : %s \n" , buffer ) ;
              
              ctRefEqv = atoi( cache ) ;
              
              printf("\n # %d orbital is equivalent to [ ctRef = %d ] ...\n" , ctRefEqv , ctRef ) ; 
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
      
      
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      
      
    
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
      
      
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      
      
    
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
      
      
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      
      
    
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
      
      
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      
      
    
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
      
      
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      
      
    
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
      
      
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      
      
    
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
      
      
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      
      
    
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
      
      
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      
      
    
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
      
      
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      
      
    
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
      
      
      else
      {
        printf("\n===> Maximum overlap for Ref-MO # %d is % 10.6f @ # %d snap orbital <===\n" , ibasis + 1 , dominateOverlap , itmp + 1 ) ;
      }
      
      
    
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
  
  
  if( leRefEqv != -1 )
  {
    nLECTtransition = nLECTtransition * 2 ;
    
    nLEorbital = nLEorbital * 2 ;
  }
  
  if( gsRefEqv != -1 )
  {
    nGSCTtransition = nGSCTtransition * 2 ;
    
    nGSorbital = nGSorbital * 2 ;
  }
  
  if( ctRefEqv != -1 )
  {
    nLECTtransition = nLECTtransition * 2 ;
  
    nGSCTtransition = nGSCTtransition * 2 ;
    
    nCTorbital = nCTorbital * 2 ;
    
  }
  
  
  nTransitionInvolved = nLECTtransition + nGSCTtransition ;
  
  nOrbitalInvolved = nLEorbital + nGSorbital + nCTorbital ;
  

  double * hda = calloc( nTransitionInvolved , sizeof( double ) ) ;
  
  dzeros( nTransitionInvolved , 1 , hda ) ;
  
  double * involvedOrbitalEnergies = calloc( nOrbitalInvolved , sizeof( double ) ) ;
  
  dzeros( nOrbitalInvolved , 1 , involvedOrbitalEnergies ) ;
  
  int * overlap_warning = calloc( nOrbitalInvolved , sizeof( int ) ) ;
  
  izeros( nOrbitalInvolved , 1 , overlap_warning ) ;




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
  
  double presentedLECTcoupling , presentedGSCTcoupling ;
  
  int lectIndicator = 0 ; int gsctIndicator = 0 ; // C-Labeling for now ... 
  
  double * energyGap = calloc( 4 , sizeof( double ) ) ; dzeros( 4 , 1 , energyGap ) ;
  
  double energyLE , energyLEeqv , energyGS , energyGSeqv , energyCT , energyCTeqv ;
  
  double minimumLECTenergyGap = 0.00 ; double minimumGSCTenergyGap = 0.00 ; 
  
  //int * overlap_warning = calloc( 8 , sizeof( int ) ) ; izeros( 8 , 1 , overlap_warning ) ;
  
  int presentedFlag1 , presentedFlag2 , presentedFlag3 , presentedFlag4 ;
  
  
  //	ATTENTION : CHARGE SEPARATION PROCESS FIRST  //
  
  if( leRefEqv != -1 && ctRefEqv != -1 )
  {
    *( hda + 0 ) = *( finalFock + ( *( moMapping + leRef - 1 ) ) * nbasis + ( *( moMapping + ctRef - 1 ) )  ) ;
    
    *( overlap_warning + 0 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + leRef - 1 ) ) ) ;
      
    *( overlap_warning + 0 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRef - 1 ) ) ) ;
        
    //*( hda + 0 ) = *( tmp_hda + 0 ) ;
    
    printf("\nLE -> CT = % 10.6E \n" , *( hda + 0 ) ) ;
    
    
    *( hda + 1 ) = *( finalFock + ( *( moMapping + leRef - 1 ) ) * nbasis + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
    
    *( overlap_warning + 1 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + leRef - 1 ) ) ) ;
    
    *( overlap_warning + 1 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
    
    //*( hda + 1 ) = *( tmp_hda + 1 ) ;

    printf("\nLE -> CT.Eqv = % 10.6E \n" , *( hda + 1 ) ) ;
    
    
    *( hda + 2 ) = *( finalFock + ( *( moMapping + leRefEqv - 1 ) ) * nbasis + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    *( overlap_warning + 2 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + leRefEqv - 1 ) ) ) ;
    
    *( overlap_warning + 2 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    //*( hda + 2 ) = *( tmp_hda + 2 ) ;

    printf("\nLE.Eqv -> CT = % 10.6E \n" , *( hda + 2 ) ) ;
    
    
    *( hda + 3 ) = *( finalFock + ( *( moMapping + leRefEqv - 1 ) ) * nbasis + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
    
    *( overlap_warning + 3 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + leRefEqv - 1 ) ) ) ;
    
    *( overlap_warning + 3 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
    
    //*( hda + 3 ) = *( tmp_hda + 3 ) ;

    printf("\nLE.Eqv -> CT.Eqv = % 10.6E \n" , *( hda + 3 ) ) ;
  
  
    
    energyCT = *( reorgEnergy + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    energyCTeqv = *( reorgEnergy + ( *( moMapping + ctRefEqv - 1 ) ) ) ;

    
    energyLE = *( reorgEnergy + ( *( moMapping + leRef - 1 ) ) ) ;
    
    energyLEeqv = *( reorgEnergy + ( *( moMapping + leRefEqv - 1 ) ) ) ;
    
    
    *( involvedOrbitalEnergies + 0 ) = energyLE ;
    
    *( involvedOrbitalEnergies + 1 ) = energyLEeqv ;
    
    *( involvedOrbitalEnergies + nLEorbital + 0 ) = energyCT ; 
    
    *( involvedOrbitalEnergies + nLEorbital + 1 ) = energyCTeqv ; 
    
  
    
    *( energyGap + 0 ) = fabs( energyCT - energyLE ) ; 
    
    *( energyGap + 1 ) = fabs( energyCTeqv - energyLE ) ;
    
    *( energyGap + 2 ) = fabs( energyCT - energyLEeqv ) ;
    
    *( energyGap + 3 ) = fabs( energyCTeqv - energyLEeqv ) ;
    
    printf("\nEnergy Gap :\nLE -> CT = % 10.6f ; LE -> CT.Eqv = % 10.6f ; LE.Eqv -> CT = % 10.6f ; LE.Eqv -> CT.Eqv = % 10.6f\n\n" , *( energyGap + 0 ) , *( energyGap + 1 ) , *( energyGap + 2 ) , *( energyGap + 3 ) ) ;
  
  
    lectIndicator = dminID( 4 , energyGap ) ;
    
    minimumLECTenergyGap = fabs( *( energyGap + itmp ) ) ;
    
    presentedLECTcoupling = *( hda + itmp ) ;
    
    presentedFlag1 = *( overlap_warning + itmp * 2 + 0 ) ;
    
    presentedFlag2 = *( overlap_warning + itmp * 2 + 1 ) ;  
  
  
  
  }
  else if( leRefEqv != -1 && ctRefEqv == -1 ) // Our Case ...
  {
    *( hda + 0 ) = *( finalFock + ( *( moMapping + leRef - 1 ) ) * nbasis + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    *( overlap_warning + 0 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + leRef - 1 ) ) ) ;
    
    *( overlap_warning + 0 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRef - 1 ) ) ) ;
        
    printf("\nLE -> CT = % 10.6E \n" , *( hda + 0 ) ) ;
    
    
    *( hda + 1 ) = *( finalFock + ( *( moMapping + leRefEqv - 1 ) ) * nbasis + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    *( overlap_warning + 1 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + leRefEqv - 1 ) ) ) ;
    
    *( overlap_warning + 1 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRef - 1 ) ) ) ;

    printf("\nLE.Eqv -> CT = % 10.6E \n" , *( hda + 1 ) ) ;
  
  
    energyCT = *( reorgEnergy + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    
    energyLE = *( reorgEnergy + ( *( moMapping + leRef - 1 ) ) ) ;
    
    energyLEeqv = *( reorgEnergy + ( *( moMapping + leRefEqv - 1 ) ) ) ;
    


    *( involvedOrbitalEnergies + 0 ) = energyLE ;
    
    *( involvedOrbitalEnergies + 1 ) = energyLEeqv ;
    
    *( involvedOrbitalEnergies + nLEorbital + 0 ) = energyCT ; 
    



    *( energyGap + 0 ) = fabs( energyCT - energyLE ) ;
    
    *( energyGap + 1 ) = fabs( energyCT - energyLEeqv ) ;
  
    
    printf("\nEnergy Gap :\nLE -> CT = % 10.6f ; LE.Eqv -> CT = % 10.6f\n\n" , *( energyGap + 0 ) , *( energyGap + 1 ) ) ;
    
    lectIndicator = dminID( 2 , energyGap ) ;  
    
    minimumLECTenergyGap = fabs( *( energyGap + itmp ) ) ;
  
    presentedLECTcoupling = *( hda + itmp ) ;
    
    presentedFlag1 = *( overlap_warning + itmp * 2 + 0 ) ;
    
    presentedFlag2 = *( overlap_warning + itmp * 2 + 1 ) ;  
  
  
  }
  else if( leRefEqv == -1 && ctRefEqv != -1 ) 
  {
    *( hda + 0 ) = *( finalFock + ( *( moMapping + leRef - 1 ) ) * nbasis + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    *( overlap_warning + 0 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + leRef - 1 ) ) ) ;
    
    *( overlap_warning + 0 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRef - 1 ) ) ) ;
        
    printf("\nLE -> CT = % 10.6E \n" , *( hda + 0 ) ) ;

    
    *( hda + 1 ) = *( finalFock + ( *( moMapping + leRef - 1 ) ) * nbasis + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
    
    *( overlap_warning + 1 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + leRef - 1 ) ) ) ;
    
    *( overlap_warning + 1 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
        
    printf("\nLE -> CT.Eqv = % 10.6E \n" , *( hda + 1 ) ) ;
  
  
  
    energyCT = *( reorgEnergy + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    energyCTeqv = *( reorgEnergy + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
    
    
    energyLE = *( reorgEnergy + ( *( moMapping + leRef - 1 ) ) ) ;
    
    

    *( involvedOrbitalEnergies + 0 ) = energyLE ;
    
    *( involvedOrbitalEnergies + nLEorbital + 0 ) = energyCT ; 
    
    *( involvedOrbitalEnergies + nLEorbital + 1 ) = energyCTeqv ; 



    *( energyGap + 0 ) = fabs( energyCT - energyLE ) ;
    
    *( energyGap + 1 ) = fabs( energyCTeqv - energyLE ) ;
  
    
    
    printf("\nEnergy Gap :\nLE -> CT = % 10.6f ; LE -> CT.Eqv = % 10.6f\n\n" , *( energyGap + 0 ) , *( energyGap + 1 ) ) ;
  
    lectIndicator = dminID( 2 , energyGap ) ;  
    
    minimumLECTenergyGap = fabs( *( energyGap + itmp ) ) ;
  
    presentedLECTcoupling = *( hda + itmp ) ;
    
    presentedFlag1 = *( overlap_warning + itmp * 2 + 0 ) ;
    
    presentedFlag2 = *( overlap_warning + itmp * 2 + 1 ) ;  
  
  
  
  }
  else if( leRefEqv == -1 && ctRefEqv == -1 ) 
  {
    *( hda + 0 ) = *( finalFock + ( *( moMapping + leRef - 1 ) ) * nbasis + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    *( overlap_warning + 0 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + leRef - 1 ) ) ) ;
    
    *( overlap_warning + 0 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRef - 1 ) ) ) ;
        
    printf("\nLE -> CT = % 10.6E \n" , *( hda + 0 ) ) ;

    energyCT = *( reorgEnergy + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    energyLE = *( reorgEnergy + ( *( moMapping + leRef - 1 ) ) ) ;
    



    *( involvedOrbitalEnergies + 0 ) = energyLE ;
    
    *( involvedOrbitalEnergies + nLEorbital + 0 ) = energyCT ; 



    *( energyGap + 0 ) = fabs( energyCT - energyLE ) ;

    printf("\nEnergy Gap :\nLE -> CT = % 10.6f\n\n" , *( energyGap + 0 ) ) ;
    
  
    lectIndicator = 0 ;
  
    minimumLECTenergyGap = fabs( *( energyGap + 0 ) ) ;
  
    presentedLECTcoupling = *( hda + 0 ) ;
    
    presentedFlag1 = *( overlap_warning + 0 * 2 + 0 ) ;
    
    presentedFlag2 = *( overlap_warning + 0 * 2 + 1 ) ;  

  }
  else
  {
    printf("\nSomething is wrong about whether there are degeneracy @ Charge Separation !!!\n") ;
    
    exit( 218 ) ;
  }
  
  
  //	ATTENTION : FOLLOWED BY CHARGE RECOMBINATION  //
  
  izeros( nOrbitalInvolved , 1 , overlap_warning ); dzeros( 4 , 1 , energyGap ) ;
  
  
  
  
  if( gsRefEqv != -1 && ctRefEqv != -1 )
  {
    *( hda + nLECTtransition + 0 ) = *( finalFock + ( *( moMapping + gsRef - 1 ) ) * nbasis + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    *( overlap_warning + 0 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + gsRef - 1 ) ) ) ;
    
    *( overlap_warning + 0 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRef - 1 ) ) ) ;
        
    printf("\nCT -> GS = % 10.6E eV\n" , *( hda + nLECTtransition + 0 ) ) ;
    
    
    *( hda + nLECTtransition + 1 ) = *( finalFock + ( *( moMapping + gsRef - 1 ) ) * nbasis + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
    
    *( overlap_warning + 1 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + gsRef - 1 ) ) ) ;
    
    *( overlap_warning + 1 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
        
    printf("\nCT.Eqv -> GS = % 10.6E eV\n" , *( hda + nLECTtransition + 1 ) ) ;
    
    
    *( hda + nLECTtransition + 2 ) = *( finalFock + ( *( moMapping + gsRefEqv - 1 ) ) * nbasis + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    *( overlap_warning + 2 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + gsRefEqv - 1 ) ) ) ;
    
    *( overlap_warning + 2 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRef - 1 ) ) ) ;
        
    printf("\nCT -> GS.Eqv = % 10.6E eV\n" , *( hda + nLECTtransition + 2 ) ) ;
    
    
    *( hda + nLECTtransition + 3 ) = *( finalFock + ( *( moMapping + gsRefEqv - 1 ) ) * nbasis + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
    
    *( overlap_warning + 3 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + gsRefEqv - 1 ) ) ) ;
    
    *( overlap_warning + 3 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
        
    printf("\nCT.Eqv -> GS.Eqv = % 10.6E eV\n" , *( hda + nLECTtransition + 3 ) ) ;
  
  
    
    energyCT = *( reorgEnergy + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    energyCTeqv = *( reorgEnergy + ( *( moMapping + ctRefEqv - 1 ) ) ) ;

    
    energyGS = *( reorgEnergy + ( *( moMapping + gsRef - 1 ) ) ) ;
    
    energyGSeqv = *( reorgEnergy + ( *( moMapping + gsRefEqv - 1 ) ) ) ;
  
    
    *( involvedOrbitalEnergies + nLEorbital + nCTorbital + 0 ) = energyGS ;
    
    *( involvedOrbitalEnergies + nLEorbital + nCTorbital + 1 ) = energyGSeqv ;
    
    
    
    *( energyGap + 0 ) = fabs( energyCT - energyGS ) ;
    
    *( energyGap + 1 ) = fabs( energyCTeqv - energyGS ) ;
    
    *( energyGap + 2 ) = fabs( energyCT - energyGSeqv ) ;
    
    *( energyGap + 3 ) = fabs( energyCTeqv - energyGSeqv ) ;
  
    printf("\nEnergy Gap :\nGS -> CT = % 10.6f ; GS -> CT.Eqv = % 10.6f ; GS.Eqv -> CT = % 10.6f ; GS.Eqv -> CT.Eqv = % 10.6f\n\n" , *( energyGap + 0 ) , *( energyGap + 1 ) , *( energyGap + 2 ) , *( energyGap + 3 ) ) ;



    gsctIndicator = dminID( 4 , energyGap ) ;
    
    minimumGSCTenergyGap = fabs( *( energyGap + itmp ) ) ;
    
    presentedGSCTcoupling = *( hda + nLECTtransition + itmp ) ;
    
    presentedFlag3 = *( overlap_warning + itmp * 2 + 0 ) ;
    
    presentedFlag4 = *( overlap_warning + itmp * 2 + 1 ) ;  
  
  
  
  
  
  }
  else if( gsRefEqv != -1 && ctRefEqv == -1 ) // Our Case ...
  {
    *( hda + nLECTtransition + 0 ) = *( finalFock + ( *( moMapping + gsRef - 1 ) ) * nbasis + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    *( overlap_warning + 0 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + gsRef - 1 ) ) ) ;
    
    *( overlap_warning + 0 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRef - 1 ) ) ) ;
        
    printf("\nCT -> GS = % 10.6E eV\n" , *( hda + nLECTtransition + 0 ) ) ;
    
    
    *( hda + nLECTtransition + 1 ) = *( finalFock + ( *( moMapping + gsRefEqv - 1 ) ) * nbasis + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    *( overlap_warning + 1 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + gsRefEqv - 1 ) ) ) ;
    
    *( overlap_warning + 1 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRef - 1 ) ) ) ;
        
    printf("\nCT -> GS.Eqv = % 10.6E eV\n" , *( hda + nLECTtransition + 1 ) ) ;
  
  
  
    
    energyCT = *( reorgEnergy + ( *( moMapping + ctRef - 1 ) ) ) ;

    
    energyGS = *( reorgEnergy + ( *( moMapping + gsRef - 1 ) ) ) ;
    
    energyGSeqv = *( reorgEnergy + ( *( moMapping + gsRefEqv - 1 ) ) ) ;
  
    
    
    *( involvedOrbitalEnergies + nLEorbital + nCTorbital + 0 ) = energyGS ;
    
    *( involvedOrbitalEnergies + nLEorbital + nCTorbital + 1 ) = energyGSeqv ;
    
    
    *( energyGap + 0 ) = fabs( energyCT - energyGS ) ;
    
    *( energyGap + 1 ) = fabs( energyCT - energyGSeqv ) ;

    printf("\nEnergy Gap :\nGS -> CT = % 10.6f ; GS.Eqv -> CT = % 10.6f\n\n" , *( energyGap + 0 ) , *( energyGap + 1 ) ) ;




    gsctIndicator = dminID( 2 , energyGap ) ;
    
    minimumGSCTenergyGap = fabs( *( energyGap + itmp ) ) ;
    
    presentedGSCTcoupling = *( hda + nLECTtransition + itmp ) ;
    
    presentedFlag3 = *( overlap_warning + itmp * 2 + 0 ) ;
    
    presentedFlag4 = *( overlap_warning + itmp * 2 + 1 ) ;  
  
  
  
  
  }
  else if( gsRefEqv == -1 && ctRefEqv != -1 ) 
  {
    *( hda + nLECTtransition + 0 ) = *( finalFock + ( *( moMapping + gsRef - 1 ) ) * nbasis + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    *( overlap_warning + 0 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + gsRef - 1 ) ) ) ;
    
    *( overlap_warning + 0 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRef - 1 ) ) ) ;
        
    printf("\nCT -> GS = % 10.6E eV\n" , *( hda + nLECTtransition + 0 ) ) ;
    
    
    *( hda + nLECTtransition + 1 ) = *( finalFock + ( *( moMapping + gsRef - 1 ) ) * nbasis + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
    
    *( overlap_warning + 1 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + gsRef - 1 ) ) ) ;
    
    *( overlap_warning + 1 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
        
    printf("\nCT.Eqv -> GS = % 10.6E eV\n" , *( hda + nLECTtransition + 1 ) ) ;
  
  
    
    energyCT = *( reorgEnergy + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    energyCTeqv = *( reorgEnergy + ( *( moMapping + ctRefEqv - 1 ) ) ) ;

    
    energyGS = *( reorgEnergy + ( *( moMapping + gsRef - 1 ) ) ) ;
    

    
    *( involvedOrbitalEnergies + nLEorbital + nCTorbital + 0 ) = energyGS ;
    
    
    
    *( energyGap + 0 ) = fabs( energyCT - energyGS ) ;
    
    *( energyGap + 1 ) = fabs( energyCTeqv - energyGS ) ;
    
    printf("\nEnergy Gap :\nGS -> CT = % 10.6f ; GS -> CT.Eqv = % 10.6f\n\n" , *( energyGap + 0 ) , *( energyGap + 1 ) ) ;

  
    gsctIndicator = dminID( 2 , energyGap ) ;
    
    minimumGSCTenergyGap = fabs( *( energyGap + itmp ) ) ;
    
    presentedGSCTcoupling = *( hda + nLECTtransition + itmp ) ;
    
    presentedFlag3 = *( overlap_warning + itmp * 2 + 0 ) ;
    
    presentedFlag4 = *( overlap_warning + itmp * 2 + 1 ) ;  
  
  
  
  }
  else if( gsRefEqv == -1 && ctRefEqv == -1 ) 
  {
    *( hda + nLECTtransition + 0 ) = *( finalFock + ( *( moMapping + gsRef - 1 ) ) * nbasis + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    *( overlap_warning + 0 * 2 + 0 ) = *( smallOverlapFlag + ( *( moMapping + gsRef - 1 ) ) ) ;
    
    *( overlap_warning + 0 * 2 + 1 ) = *( smallOverlapFlag + ( *( moMapping + ctRef - 1 ) ) ) ;
        
    printf("\nCT -> GS = % 10.6E eV\n" , *( hda + nLECTtransition + 0 ) ) ;
  
    
    energyCT = *( reorgEnergy + ( *( moMapping + ctRef - 1 ) ) ) ;
    
    energyGS = *( reorgEnergy + ( *( moMapping + gsRef - 1 ) ) ) ;
    
   
    *( involvedOrbitalEnergies + nLEorbital + nCTorbital + 0 ) = energyGS ;
    
    
    
    *( energyGap + 0 ) = fabs( energyCT - energyGS ) ;
    
    printf("\nEnergy Gap :\nGS -> CT = % 10.6f\n\n" , *( energyGap + 0 ) ) ;

    
    
    gsctIndicator = 0 ; 
    
    minimumGSCTenergyGap = fabs( *( energyGap + 0 ) ) ;
    
    presentedGSCTcoupling = *( hda + nLECTtransition + 0 ) ;
    
    presentedFlag3 = *( overlap_warning + 0 * 2 + 0 ) ;
    
    presentedFlag4 = *( overlap_warning + 0 * 2 + 1 ) ;  
  
    
  
  }
  else
  {
    printf("\nSomething is wrong about whether there are degeneracy @ Charge Re-Combination !!!\n") ;
    
    exit( 218 ) ;
  }
  
  
  


  // -----------------> Recording the Bridge Orbital Energies ... and nearest-neighbour coupling ... <-------------//

  int firstBridge , lastBridge , currentBridge , nextBridge , ibridge ; 
  
  double * vGSBridge , * vLEBridge , * vBridgeCT , * vBB ; 
  
  int nGSBridge , nLEBridge , nBridgeCT , nBB , iBB ;   
  
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

    // ---> < LE | H | First Bridge >

    if( ( leRefEqv != -1 ) && ( nbridge != 1 ) )
    {
      nLEBridge = 4 ;
    
      vLEBridge = calloc( nLEBridge , sizeof( double ) ) ;
    
      *( vLEBridge + 0 ) = *( finalFock + ( *( moMapping + leRef - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
    
      *( vLEBridge + 1 ) = *( finalFock + ( *( moMapping + leRefEqv - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
      
      *( vLEBridge + 2 ) = *( finalFock + ( *( moMapping + leRef - 1 ) ) * nbasis + ( *( moMapping + lastBridge - 1 ) ) ) ;
    
      *( vLEBridge + 3 ) = *( finalFock + ( *( moMapping + leRefEqv - 1 ) ) * nbasis + ( *( moMapping + lastBridge - 1 ) ) ) ;
    
    
    }
    else if( ( leRefEqv == -1 ) && ( nbridge != 1 ) )
    {
      nLEBridge = 2 ;
    
      vLEBridge = calloc( nLEBridge , sizeof( double ) ) ;
    
      *( vLEBridge + 0 ) = *( finalFock + ( *( moMapping + leRef - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
      
      *( vLEBridge + 1 ) = *( finalFock + ( *( moMapping + leRef - 1 ) ) * nbasis + ( *( moMapping + lastBridge - 1 ) ) ) ;
  
    }
    else if( ( leRefEqv != -1 ) && ( nbridge == 1 ) )
    {
      nLEBridge = 2 ;
    
      vLEBridge = calloc( nLEBridge , sizeof( double ) ) ;
    
      *( vLEBridge + 0 ) = *( finalFock + ( *( moMapping + leRef - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
      
      *( vLEBridge + 1 ) = *( finalFock + ( *( moMapping + leRefEqv - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
  
    }
    else
    {
      nLEBridge = 1 ; 
      
      vLEBridge = calloc( nLEBridge , sizeof( double ) ) ;
      
      *( vLEBridge + 0 ) = *( finalFock + ( *( moMapping + leRef - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
    
    
    }



    // ---> < Last Bridge | H | CT >

    if( ( ctRefEqv != -1 ) && ( nbridge != 1 ) )
    {
      nBridgeCT = 4 ;
    
      vBridgeCT = calloc( nBridgeCT , sizeof( double ) ) ;
    
      *( vBridgeCT + 0 ) = *( finalFock + ( *( moMapping + ctRef - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
    
      *( vBridgeCT + 1 ) = *( finalFock + ( *( moMapping + ctRefEqv - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
      
      *( vBridgeCT + 2 ) = *( finalFock + ( *( moMapping + ctRef - 1 ) ) * nbasis + ( *( moMapping + lastBridge - 1 ) ) ) ;
    
      *( vBridgeCT + 3 ) = *( finalFock + ( *( moMapping + ctRefEqv - 1 ) ) * nbasis + ( *( moMapping + lastBridge - 1 ) ) ) ;
      
    }
    else if( ( ctRefEqv == -1 ) && ( nbridge != 1 ) )
    {
      nBridgeCT = 2 ;
    
      vBridgeCT = calloc( nBridgeCT , sizeof( double ) ) ;
    
      *( vBridgeCT + 0 ) = *( finalFock + ( *( moMapping + ctRef - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
  
      *( vBridgeCT + 1 ) = *( finalFock + ( *( moMapping + ctRef - 1 ) ) * nbasis + ( *( moMapping + lastBridge - 1 ) ) ) ;
  
    }
    else if( ( ctRefEqv != -1 ) && ( nbridge == 1 ) )
    {
      nBridgeCT = 2 ;
    
      vBridgeCT = calloc( nBridgeCT , sizeof( double ) ) ;
    
      *( vBridgeCT + 0 ) = *( finalFock + ( *( moMapping + ctRef - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
  
      *( vBridgeCT + 1 ) = *( finalFock + ( *( moMapping + ctRefEqv - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
  
    }
    else
    {
      nBridgeCT = 1 ;
    
      vBridgeCT = calloc( nBridgeCT , sizeof( double ) ) ;
    
      *( vBridgeCT + 0 ) = *( finalFock + ( *( moMapping + ctRef - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
  
    }

  
    // ---> < GS | H | First Bridge >
  
    if( ( gsRefEqv != -1 ) && ( nbridge != 1 ) )
    {
      nGSBridge = 4 ;
    
      vGSBridge = calloc( nGSBridge , sizeof( double ) ) ;
    
      *( vGSBridge + 0 ) = *( finalFock + ( *( moMapping + gsRef - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
    
      *( vGSBridge + 1 ) = *( finalFock + ( *( moMapping + gsRefEqv - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
    
      *( vGSBridge + 2 ) = *( finalFock + ( *( moMapping + gsRef - 1 ) ) * nbasis + ( *( moMapping + lastBridge - 1 ) ) ) ;
    
      *( vGSBridge + 3 ) = *( finalFock + ( *( moMapping + gsRefEqv - 1 ) ) * nbasis + ( *( moMapping + lastBridge - 1 ) ) ) ;
    
    }
    else if( ( gsRefEqv == -1 ) && ( nbridge != 1 ) )
    {
      nGSBridge = 2 ;
    
      vGSBridge = calloc( nGSBridge , sizeof( double ) ) ;
    
      *( vGSBridge + 0 ) = *( finalFock + ( *( moMapping + gsRef - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
  
      *( vGSBridge + 1 ) = *( finalFock + ( *( moMapping + gsRef - 1 ) ) * nbasis + ( *( moMapping + lastBridge - 1 ) ) ) ;
  
    }
    else if( ( gsRefEqv != -1 ) && ( nbridge == 1 ) )
    {
      nGSBridge = 2 ;
    
      vGSBridge = calloc( nGSBridge , sizeof( double ) ) ;
    
      *( vGSBridge + 0 ) = *( finalFock + ( *( moMapping + gsRef - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
  
      *( vGSBridge + 1 ) = *( finalFock + ( *( moMapping + gsRefEqv - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;
  
    }
    else
    {
      nGSBridge = 1 ;
    
      vGSBridge = calloc( nGSBridge , sizeof( double ) ) ;
    
      *( vGSBridge + 0 ) = *( finalFock + ( *( moMapping + gsRef - 1 ) ) * nbasis + ( *( moMapping + firstBridge - 1 ) ) ) ;

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
  
  
  double leLocalization = 0.00 , leEqvLocalization = 0.00 , ctLocalization = 0.00 , ctEqvLocalization = 0.00 , gsLocalization = 0.00 , gsEqvLocalization = 0.00 ;
  
  
  /* double calcSquaredLocalization( double * orbital , int nBasis , int * basisList , int nComponent ) */
  
  
  //--> LE <--//
  
  orbitalID = *( moMapping + leRef - 1 ) ;
  
  if( exD == YES )
  {
    leLocalization = calcSquaredLocalization( ( transReOrgMO + nbasis * orbitalID ) , nbasis , dBasisList , nbasisD ) ;
  }
  
  //--> LE.Eqv <--//
  
  if( leRefEqv != -1 )
  {
    orbitalID = *( moMapping + leRefEqv - 1 ) ;
  
    if( exD == YES )
    {
      leEqvLocalization = calcSquaredLocalization( ( transReOrgMO + nbasis * orbitalID ) , nbasis , dBasisList , nbasisD ) ;
    }
  
  }
  
  
  //--> CT <--//
  
  orbitalID = *( moMapping + ctRef - 1 ) ;
  
  if( exA == YES )
  {
    ctLocalization = calcSquaredLocalization( ( transReOrgMO + nbasis * orbitalID ) , nbasis , aBasisList , nbasisA ) ;
  }
  
  //--> CT.Eqv <--//
  
  if( ctRefEqv != -1 )
  {
    orbitalID = *( moMapping + ctRefEqv - 1 ) ;
  
    if( exA == YES )
    {
      ctEqvLocalization = calcSquaredLocalization( ( transReOrgMO + nbasis * orbitalID ) , nbasis , aBasisList , nbasisA ) ;
    }
  
  }
  
  
  
  //--> GS <--//
  
  orbitalID = *( moMapping + gsRef - 1 ) ;
  
  if( exD == YES )
  {
    gsLocalization = calcSquaredLocalization( ( transReOrgMO + nbasis * orbitalID ) , nbasis , dBasisList , nbasisD ) ;
  }
  
  //--> GS.Eqv <--//
  
  if( gsRefEqv != -1 )
  {
    orbitalID = *( moMapping + gsRefEqv - 1 ) ;
  
    if( exD == YES )
    {
      gsEqvLocalization = calcSquaredLocalization( ( transReOrgMO + nbasis * orbitalID ) , nbasis , dBasisList , nbasisD ) ;
    }
  
  }
  
  
  
  printf("\nRelevant MO Localization Circumstances ...\n\n" ) ;
  
  printf("\n% 10.6E\t% 10.6E\t% 10.6E\t% 10.6E\t% 10.6E\t% 10.6E\n\n" , leLocalization , leEqvLocalization , ctLocalization , ctEqvLocalization , gsLocalization , gsEqvLocalization ) ;
  
  
  
  double presentLECT_localization = 0.00 , presentGSCT_localization = 0.00 ;
  
  if( leLocalization >= leEqvLocalization )
  {
    if( ctLocalization >= ctEqvLocalization )
    {
      presentLECT_localization = *( finalFock + ( *( moMapping + leRef - 1 ) ) * nbasis + ( *( moMapping + ctRef - 1 ) ) ) ;
    }
    else
    {
      presentLECT_localization = *( finalFock + ( *( moMapping + leRef - 1 ) ) * nbasis + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
    }

  }
  else
  {
    if( ctLocalization >= ctEqvLocalization )
    {
      presentLECT_localization = *( finalFock + ( *( moMapping + leRefEqv - 1 ) ) * nbasis + ( *( moMapping + ctRef - 1 ) ) ) ;
    }
    else
    {
      presentLECT_localization = *( finalFock + ( *( moMapping + leRefEqv - 1 ) ) * nbasis + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
    }

  }
  
  
  
  if( gsLocalization >= gsEqvLocalization )
  {
    if( ctLocalization >= ctEqvLocalization )
    {
      presentGSCT_localization = *( finalFock + ( *( moMapping + gsRef - 1 ) ) * nbasis + ( *( moMapping + ctRef - 1 ) ) ) ;
    }
    else
    {
      presentGSCT_localization = *( finalFock + ( *( moMapping + gsRef - 1 ) ) * nbasis + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
    }

  }
  else
  {
    if( ctLocalization >= ctEqvLocalization )
    {
      presentGSCT_localization = *( finalFock + ( *( moMapping + gsRefEqv - 1 ) ) * nbasis + ( *( moMapping + ctRef - 1 ) ) ) ;
    }
    else
    {
      presentGSCT_localization = *( finalFock + ( *( moMapping + gsRefEqv - 1 ) ) * nbasis + ( *( moMapping + ctRefEqv - 1 ) ) ) ;
    }

  }
  
  
  
  
  
  // ------------------> Outputing ... <---------------------- //
  
  
  // ATTENTION : Here we output fabs( presented*coupling ) instead of themselves ...
  
  
  poutHDAFile = fopen( outputHDAName , "wb+" ) ;
  
  fprintf( poutHDAFile , "\n" ) ;
  
  

  for( iload = 0 ; iload < nTransitionInvolved ; iload ++ )
  {
    fprintf( poutHDAFile , "% 10.6E\t" , fabs( *( hda + iload ) ) ) ;
  
  }
  
  fprintf( poutHDAFile , "% 3d\t% 3d\t" , lectIndicator + 1 , gsctIndicator + 1 ) ; // easier for matlab to use, here is human labeling
  
  
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
  
    for( itmp = 0 ; itmp < nBridgeCT ; itmp ++ )
    {
      fprintf( poutHDAFile , "% 10.6E\t" , *( vBridgeCT + itmp ) ) ;
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
  
  
  
  
  
  //fprintf( poutHDAFile , "\n" ) ;
  

  
  /*
  fprintf( poutHDAFile , "% 10.6E\t% 10.6E\t% 10.6E\t% 10.6E\t% 3d\t% 3d\t% 3d\t% 3d\n\n\n" , 
           fabs( presentedLECTcoupling ) , fabs( presentedGSCTcoupling ) , minimumLECTenergyGap , minimumGSCTenergyGap , presentedFlag1 , presentedFlag2 , presentedFlag3 , presentedFlag4 ) ;
  
  */
  

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





