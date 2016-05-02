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

typedef struct itpform
{
  int nr ;
  char type[ 4 ] ;
  int resnr ;
  char residue[ 10 ] ;
  char atomname[ 10 ] ;
  int cgnr ;
  double charge ;
  double mass ;

} ITP ;


typedef struct topform
{
  char typename[10];
  char bondtype[10];
  double mass; 
  double charge ;
  char ptype[3] ;
  double sigma ;
  double epsilon ;

} TOP ;



int main( int argc, char * argv[] )
{
  // -----------------------------------------> Variables ... <--------------------------------------------- //
  
  char ** pcmd ; pcmd = argv ;
  
  FILE * psoluteITP , * psolventITP , * pindex , * pEqGRO , * pinpcrd , * pinpvel ;
  
  FILE * debug ;
  
  char soluteITPName[ MAXCHARINLINE ] , solventITPName[ MAXCHARINLINE ] , indexFileName[ MAXCHARINLINE ] , eqGROFileName[ MAXCHARINLINE ] , inpCRDFileName[ MAXCHARINLINE ] , inpVELFileName[ MAXCHARINLINE ] ;
  
  char outCRDFileName[ MAXCHARINLINE ] , outVELFileName[ MAXCHARINLINE ] ;
  
  char outGROFileName[ MAXCHARINLINE ] ;
  
  char pGroupName[ 50 ] , qGroupName[ 50 ] ;
  
  int referenceAtomNumber ;
  
  int natom , natomsoluteGRO , natomsoluteSelect , natomITP , natomtypes ;
  
  int ncart ;
  
  int iatom , icart ;
  
  double * comEqP , * comEqQ , * coordinateEqN , * coordinateNEqN , * comNEqP , * comNEqQ ;
  
  double * comVelocityEqP , * comVelocityNEqP ;
  
  double * vectorEqPQ , * vectorNEqPQ , * vectorEqPN , * vectorNEqPN ;
  
  double distThreshold ;
 
  double vdwFactor = 1.000 ; 

  
  double done = 1.0000 ; double dzero = 0.0000 ; 
  
  int ithree = 3 ; int ione = 1 ; int izero = 0 ;



  
  double dtmp , dtmp2 ; double dtmpArray[ MAXCHARINLINE ] ;
  
  int itmp ;
  
  char ctmp ;
  
  char tmpString[ MAXCHARINLINE ] ;

  
  int icmd ; int info ;
  
  
  // ----------------------------------> Utility Functions ... <---------------------------------- //
  
  void tellsymbol( char * atomMDname , char * atomSymbol) ;

  //double tellvdwradius( TOP * p , int howmanytypes , char * whattype ) ;

  double tellvdwradius( char * whatsymbol ) ;

  int tellatomtype( ITP * p , int howmanytypes , char * whatatom , char * at ) ; 




  // ----------------------------------> Recording Command-Line Arguments ... <---------------------------------- //
  
  time_t current_time;

  time( &current_time );

  char now[ 300 ] ;

  strcpy( now , ctime( &current_time ) );

  int lennow = strlen( now ) ;

  *( now + lennow - 1 ) = ' ';

    
    
    printf("\n**********************************************************************\n");
      printf("* G_NEQALIGN_D : Align mol. to reference orientation to avoid clash. *\n");
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

 
  
  // ------------------------------------------> Defaults ... <------------------------------------------ //
  
  referenceAtomNumber = 3 ;
  
  strcpy( pGroupName , "Donor" ) ;
  
  strcpy( qGroupName , "Acceptor" ) ;
  
  strcpy( soluteITPName , "solute.itp" );
  
  strcpy( solventITPName , "dcm.itp" );
  
  strcpy( indexFileName , "system.index" );
  
  strcpy( inpCRDFileName , "solute.crd" );
  
  strcpy( inpVELFileName , "dummy" ); // Default : No operations on velocities ;
  
  
  
  
  
  // -------------------------------> Parsing command-line arguments ... <---------------------------------- //
  
  int exs = 18 ;
  
  int exv = 80 ; 
  
  int exVelocity = NO ;
  
  int velDump = NO ;
  
  int vdwCheck = YES ;
  
  int exo = 22 ;
  
  int ext = 88 ;
  
  
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
	      
	      case 'i' : strcpy( soluteITPName , *( ++ pcmd ) ) ; 
			 
			         printf("\nCommand-line argument indicates : Input solute itp File name : %s ...\n" , soluteITPName ); 
	      
	                 icmd = icmd + 2 ; 
	                 
	                 break ;
 
 	      case 'j' : strcpy( solventITPName , *( ++ pcmd ) ) ; 
			 
			         printf("\nCommand-line argument indicates : Input solvent itp File name : %s ...\n" , solventITPName ); 
	      
	                 icmd = icmd + 2 ; 
	                 
	                 break ;
   
              case 'n' : strcpy( indexFileName , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Input index File name : %s ...\n" , indexFileName ); 
	         
	                 icmd = icmd + 2 ;
	         
	                 break ;

	      case 'r' : strcpy( inpCRDFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Input crd file name : %s ...\n" , inpCRDFileName ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ; 
	      
	      case 'v' : strcpy( inpVELFileName , *( ++ pcmd ) ); 
	      
	                 exv = 81 ;
	      
	                 if( strcmp( inpVELFileName , "dummy" ) != 0 )  
	                 {
	                   velDump = YES ;
	                   
	                   printf("\nCommand-line argument indicates : Input vel file name : %s ...\n" , inpVELFileName ); 
	                 }
	                 else
	                 {
	                   velDump = NO ;
	                   
	                   printf("\nSo you chose to not operate on velocities ... that's ... interesting ...\n");
	                 
	                 }
	               
	                 icmd = icmd + 2 ; 
	                 
	                 break ; 
	                 
	      case 'o' : strcpy( outGROFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Output GRO file name : %s ...\n" , outGROFileName ); 
	                 
	                 exo = 23 ;
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ; 
	      
	      case 'c' : strcpy( eqGROFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Input .gro File name : %s ...\n" , eqGROFileName ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 //exgro = 19 ;
	                 
	                 break ; 
	                 
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

	      case 't' : distThreshold = atof( *( ++ pcmd ) ) ;
	      
	                 printf("\nCommand-line argument indicates : The threshold for minimum distance is % 10.6f Angstrom ...\n" , distThreshold );

                         distThreshold = distThreshold / 10.00 ; // Change to nm unit since all the calculations here are in nm unit ... I know, it IS weird ...
                         
                         icmd = icmd + 2 ;
                         
                         ext = 89 ;
                         
                         break;

	      case 'w' : strcpy( tmpString , *( ++ pcmd ) ); 
	      
	                 if( strcmp( tmpString , "YES" ) == 0 || strcmp( tmpString , "yes" ) == 0 || strcmp( tmpString , "Yes" ) == 0 )
	                 {
	                   printf("\nCommand-line argument indicates : van der Waals contacting check WILL be performed with default scaling factor = 1.00\n" );

                           vdwFactor = 1.000 ;
			   
			   vdwCheck = YES ;
                         }
                         else if( strcmp( tmpString , "NO" ) == 0 || strcmp( tmpString , "no" ) == 0 || strcmp( tmpString , "No" ) == 0 )
                         {
	                   printf("\nCommand-line argument indicates : van der Waals contacting check will NOT be performed" );

                           vdwCheck = NO ;
                         }
			 else if( *( tmpString ) <= '9' && *( tmpString ) >= '0' )
			 {
			   vdwCheck = YES ;

			   vdwFactor = atof( tmpString ) ;
			 
	                   printf("\nCommand-line argument indicates : van der Waals contacting check WILL be performed with user-defined scaling factor = % 10.6f\n" , vdwFactor );
			 
			 }
			 else
			 {
			   printf("\nInvalid choice of van der Waals radius checking ... Aborting ...\n") ;

			   exit( 409 ) ;
			 
			 }

                         icmd = icmd + 2 ;
                     
                         break;


	      case 'h' : printf("\nUsage:  %s [ -i 'input solute itp file name' ] [ -j 'input solvent itp file name' ] [ -n 'input index file name (in GROMACS .ndx format)' ]" , * argv);
	                 printf("\n           [ -r 'name of input .crd file of solute' ] [ -v name input .vel file of solute ][ -o ouput .gro file name ]");
	                 printf("\n           [ -c 'input EqMD gro file of solvent' ][ -p P Group Name ][ -q Q Group Name ][ -N atom number of L-Shape reference atom ]");
	                 printf("\n           [ -w whether to perform van der Waals contacting check ( YES / Yes / yes ) or ( NO / No / no ) or floating point number to indicate scaling factor ]");
	                 printf("\n           [ -s # of atoms in solute molecule ] [ -t distance threshold to accept one alignment , unit = Angstrom ]\n\n" ); 
	                 
	                 printf("\nNote : 1) If solute velocity file name is set as \"dummy\", then no operation on velocity will be done and velocity info will not be inserted into output file.\n");
	                 
	                 printf("\n       2) If .gro file of solvent is provided without velocity info, then no velocity could be put into the output .gro file. Hence no velocity operation will be done.\n");
	                 
			 printf("\n       3) For \"-w\" option, YES/Yes/yes will cause the vdw-check to perform with default scaling 1.00 while NO/No/no will shut down the vdw-check.\n");
			 printf("\n          Specifying a floating number will also cause the vdw-check to perform but the floating number will be the user-defined scaling factor (vdwFactor).\n\n\n");

			 printf("\n       4) For \"-p\" and \"-q\" options, default value are [ Donor ] and [ Acceptor ], respectively.\n\n\n");

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
  
  if( exv == 80 )
  {
    printf("\nBy default ... no operations will be applied on velocity information ...\n") ;
    
    velDump = NO ;
  
  }
  
  
  
  int lenEqGROFileName = strlen( eqGROFileName ) ;
  
  if( exo == 22 )
  {
    strncpy( outGROFileName , eqGROFileName , lenEqGROFileName - 4 ) ;
    
    *( outGROFileName + lenEqGROFileName - 4 ) = '\0' ;
    
    strcat( outGROFileName , ".fit.gro" ) ;
    
    printf("\nBy default , output GRO file name will be %s ...\n\n" , outGROFileName ) ;
  
  }
  
  
  if( ext == 88 )
  {
    distThreshold = -0.10 ;
  }
  
  
  
  
  
  
  
  
  
  
  
  // -------------------------------> Verify File Access ... <---------------------------------- //
  
/*
FILE * psoluteITP , * pindex , * pEqGRO , * pinpcrd , * pinpvel ;
char soluteITPName[ MAXCHARINLINE ] , indexFileName[ MAXCHARINLINE ] , eqGROFileName[ MAXCHARINLINE ] , inpCRDFileName[ MAXCHARINLINE ] , inpVELFileName[ MAXCHARINLINE ] ;
*/

  if( ( psoluteITP = fopen( soluteITPName , "r" ) ) == NULL )
  {
    printf("\nUser defined .itp file does not exist ... \n");
   
    exit( 3 );
  }

  if( ( psolventITP = fopen( solventITPName , "r" ) ) == NULL )
  {
    printf("\nUser defined solvent.itp file does not exist ... \n");
   
    exit( 3 );
  }

  if( ( pindex = fopen( indexFileName , "r" ) ) == NULL )
  {
    printf("\nUser defined index file does not exist ... \n");
   
    exit( 3 );
  }

  if( ( pinpcrd = fopen( inpCRDFileName , "r" ) ) == NULL )
  {
    printf("\nUser defined crd file does not exist ... \n");
   
    exit( 3 );
  }
  
  
  if( velDump == YES )
  {
    if( ( pinpvel = fopen( inpVELFileName , "r" ) ) == NULL )
    {
      printf("\nUser defined vel file does not exist ... \n");
   
      exit( 3 );
    }
  
  }
  
  
  // -------------------------------> Reading GRO File ... <---------------------------------- //
  
  // =====> pre-Loading .gro file to get NAtom info ...  
  
  char grotitlestring[MAXLINE];
  
  int iline = 3 ;
  
  int iload = 0 ;
  
  int blank_signal , groinfo ;

  char buffer[ MAXCHARINLINE ] ;
  
  char cache[ MAXCHARINLINE ] ;

  char tmp_char ;
  
  int natomgroline , natomgrotitle ;
  

  

  
  
  if( ( pEqGRO = fopen( eqGROFileName , "r" ) ) == NULL )
  {
    printf("\nUser defined .gro file %s does not exist ... \n" , eqGROFileName );
   
    exit( 3 );
  }
  else
  {
    rewind( pEqGRO );
    
    //printf("\nCurrent character is %c ... \n" , fgetc( pEqGRO ) );
    
    fskip( pEqGRO , 1 );

    fscanf( pEqGRO , "%d" , &natomgrotitle );
    
    fskip( pEqGRO , 1 );
    
    printf("\n Second line of .gro file says it is describing %d atoms ... \n\n" , natomgrotitle );
    
    
    while( ( groinfo = freadline( buffer , MAXCHARINLINE , pEqGRO , ';' ) ) != 0 )
    {
      itmp = inLineWC( buffer ) ;
      
      break ;
    }
    
    if( itmp == 6 )
    {
      exVelocity = NO ;
      
      printf("\nI see there is no velocity information in .gro file ...\n\n") ;
      
      if( velDump == YES )
      {
        printf("\nEh... There is NO velocity information in your .gro file ... I don't think I can output for you ...\n") ;
        
        velDump = NO ;
      }
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
    
  }
 
  printf("\nNow let's pre-load the .gro file ans see how many atoms it is describing ... \n");

  rewind( pEqGRO ) ;
    
  fskip( pEqGRO , 2 ) ;

  while( ( groinfo = freadline( buffer , MAXCHARINLINE , pEqGRO , ';' ) ) != 0 )
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
      
      if( ( tmp_char = getfirst( buffer ) ) == ';' )
      {
        //printf("\nThis is a comment line ... So nothing will be loaded ...\n");
        
        //fskip( pinputfile , 1 );
        
        continue ;
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

  natomgroline = iload - 1 ; // Because the last line is boxvector line .
  
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
  
  rewind( pEqGRO );

  ncart = 3 * natom ;
 
  if( exs == 18 )
  {
    natomsoluteSelect = natom ;
  }
  else if( exs == 19 && natomsoluteSelect > natom ) // Will be dead ... 
  {
    printf("\nThere are only %d atoms in this system ... you cannot select more than that ... \n" , natom );
    
    if( natomgrotitle > natomgroline && natomsoluteSelect <= natomgrotitle )
    {
      printf("\nAlthough ... the second line of your initial .gro file did indicate there were supposed to be %d atoms in system ... So go back and make sure what you are trying to do ... \n" , natomgrotitle );
    }
    else if( natomgrotitle < natomgroline && natomsoluteSelect <= natomgroline )
    {
      printf("\nAlthough ... your initial .gro file did describe %d atoms in system ... So go back and make sure what you are trying to do ... \n" , natomgroline );
    }
    
    exit( 78 );
  }
  else if( exs == 19 && natomsoluteSelect <= natom )
  {
    printf("\nYou have selected %d atoms for rotation ... There are %d atoms en toto in this system ... \n" , natomsoluteSelect , natom );
  
  }
  else
  {
    printf("\nSomething is wrong with the atom selection process ... NAtom = %d , NAtomSelect = %d ... \n" , natom , natomsoluteSelect );
    
    exit( 78 );
  }
  


  // =====> Summarizing natom situations
  
  printf("\n%d atoms are selected as solute ...\n" , natomsoluteSelect ) ;




  // =====> Actually Reading the GRO File ...
  
  GRO EqAtomList[ natom ] ; 
  
  rewind( pEqGRO );
    
  fskip( pEqGRO , 1 );

  fskip( pEqGRO , 1 );
    
  printf("\nNow let's read the actual .gro file  ... \n");
  
  iload = 0 ; iline = 0 ;
  
  while( ( groinfo = freadline( buffer , MAXCHARINLINE , pEqGRO , ';' ) ) != 0 )
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

        //printf( "%s\t" , atomlist[ iatom ].resname );

        //sscanf( pEqGRO , "%s" , EqAtomList[ iload ].atomname );
        
        strpickword( buffer , 2 , cache ) ;
        
        strcpy( EqAtomList[ iload ].atomname , cache ) ;

        //printf( "%s" , atomlist[ iatom ].atomname );

        //sscanf( pEqGRO , "%d" , &EqAtomList[ iload ].atomnumber );
        
        strpickword( buffer , 3 , cache ) ;
        
        EqAtomList[ iload ].atomnumber = atoi( cache ) ;

        //printf( "\nWorking on No. %d atom ...\n" , EqAtomList[ iatom ].atomnumber );

        //sscanf( pEqGRO , "%lf" , &EqAtomList[ iload ].cx ); //printf("\n Cx is %lf ...\t" , EqAtomList[ iatom ].cx);
        
        strpickword( buffer , 4 , cache ) ; EqAtomList[ iload ].cx = atof( cache ) ;
        
        //sscanf( pEqGRO , "%lf" , &EqAtomList[ iload ].cy ); //printf("\n Cy is %lf ...\t" , EqAtomList[ iatom ].cy);
        
        strpickword( buffer , 5 , cache ) ; EqAtomList[ iload ].cy = atof( cache ) ;
        
        //sscanf( pEqGRO , "%lf" , &EqAtomList[ iload ].cz ); //printf("\n Cz is %lf ...\n\n" , EqAtomList[ iatom ].cz);
        
        strpickword( buffer , 6 , cache ) ; EqAtomList[ iload ].cz = atof( cache ) ;
        
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
    
    if( iload == natom ) break ;
    
    
    iline ++ ;

  }
  
  if( natomgrotitle < natomgroline )
  {
    fskip( pEqGRO , natomgroline - natomgrotitle ) ;
  }
  
  double boxvector[ 3 ]; 
  
  fscanf( pEqGRO , "%lf" , boxvector + 0 );

  fscanf( pEqGRO , "%lf" , boxvector + 1 );

  fscanf( pEqGRO , "%lf" , boxvector + 2 );



  // -------------------------------> Reading Solute ITP File For Mass Information ... <---------------------------------- //
  
  
  // =====> Pre-Loading ... 
  
  int natomSoluteDataBase ;
  
  printf("\n---> Now let's pre-load the .itp file for SOLUTE and see how many atoms it is describing ... <---\n");
  
  rewind( psoluteITP ) ;
  
  iline = 1 ;
  
  iload = 0 ;
  
  fsearch( psoluteITP , "atoms" ) ;
  
  fskip( psoluteITP , 1 );

  int itpinfo ;
  
  while( ( itpinfo = freadline( buffer , MAXCHARINLINE , psoluteITP , ';' ) ) != 0 )
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
      
      if( ( tmp_char = getfirst( buffer ) ) == ';' )
      {
        printf("\nThis is a comment line ... So nothing will be loaded ...\n");
        
        //fskip( pinputfile , 1 );
      }
      else if( ( tmp_char = getfirst( buffer ) ) == '[' )
      {
        printf("\nSTARTING OF NEW DIRECTIVE ... END READING ... \n\n");
        
        break ;
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

  natomSoluteDataBase = iload ;
  
  printf("\nThere are %d atoms described in itp file ...\n" , natomSoluteDataBase );

  if( natomSoluteDataBase != natomsoluteSelect )
  {
    printf("\nWhile you specified you wanted %d atoms as the solute, in your itp file, there are only %d atoms being described ...\n" , natomsoluteSelect , natomSoluteDataBase );
    
    printf("\nPlease check your itp file and make sure the # of atoms match ...\n");
    
    printf("\nWe will continue anyway but once we find any required info does not exist in your itp file , we will terminate this process immediately ...\n\n") ;
  }



  // =====> Loading ... 
  
  printf("\n---> Now let's actually load the .itp file for SOLUTE ... <---\n");
  
  ITP soluteAtomDataBase[ natomSoluteDataBase ] , * psoluteAtomDataBase ;
  
  psoluteAtomDataBase = soluteAtomDataBase ;
  
  iline = 1 ;
  
  iload = 0 ;
  
  rewind( psoluteITP ) ;
  
  fsearch( psoluteITP , "atoms" ) ;
  
  fskip( psoluteITP , 1 );
  
  while( ( itpinfo = freadline( buffer , MAXCHARINLINE , psoluteITP , ';' ) ) != 0 )
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
      else if( ( tmp_char = getfirst( buffer ) ) == '[' )
      {
        //printf("\nSTARTING OF NEW DIRECTIVE ... END READING ... \n\n");
        
        break ;
      }
      else
      {
        //printf("\nLine reads : %s ...\n" , buffer );
        
        iload ++ ;
        
        strpickword( buffer , 1 , cache );
        
        soluteAtomDataBase[ iload -1 ].nr = atoi( cache ) ;
        
        strpickword( buffer , 2 , cache );
        
        strcpy( soluteAtomDataBase[ iload -1 ].type , cache );
        
        strpickword( buffer , 3 , cache );
        
        soluteAtomDataBase[ iload -1 ].resnr = atoi( cache ) ;
        
        strpickword( buffer , 4 , cache );
        
        strcpy( soluteAtomDataBase[ iload -1 ].residue , cache );
        
        strpickword( buffer , 5 , cache );
        
        strcpy( soluteAtomDataBase[ iload -1 ].atomname , cache );
        
        strpickword( buffer , 6 , cache );
        
        soluteAtomDataBase[ iload -1 ].cgnr = atoi( cache ) ;
        
        strpickword( buffer , 7 , cache );
        
        soluteAtomDataBase[ iload -1 ].charge = atof( cache ) ;
        
        strpickword( buffer , 8 , cache );
        
        soluteAtomDataBase[ iload -1 ].mass = atof( cache ) ;
        
        //printf("\nCharge of this atom is %lf ...\n" , soluteAtomDataBase[ iload -1 ].charge ) ;
 
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
  
  
  // =====> Pre-Loading ... Solvent
    
  int natomSolventDataBase ;
  
  printf("\n\n\n---> Now let's pre-load the .itp file for SOLVENT and see how many atoms it is describing ... <---\n");
  
  iline = 1 ;
  
  iload = 0 ;
  
  rewind( psolventITP ) ;
  
  fsearch( psolventITP , "atoms" ) ;
  
  fskip( psolventITP , 1 );
  
  while( ( itpinfo = freadline( buffer , MAXCHARINLINE , psolventITP , ';' ) ) != 0 )
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
      
      if( ( tmp_char = getfirst( buffer ) ) == ';' )
      {
        printf("\nThis is a comment line ... So nothing will be loaded ...\n");
        
        //fskip( pinputfile , 1 );
      }
      else if( ( tmp_char = getfirst( buffer ) ) == '[' )
      {
        printf("\nSTARTING OF NEW DIRECTIVE ... END READING ... \n\n");
        
        break ;
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

  natomSolventDataBase = iload ;
  
  printf("\nThere are %d atoms described in solvent itp file ...\n" , natomSolventDataBase );


  // =====> Loading ... Solvent
  
  printf("\n---> Now let's actually load the .itp file for SOLVENT ... <---\n\n\n");
  
  ITP solventAtomDataBase[ natomSolventDataBase ] , * psolventAtomDataBase;
  
  psolventAtomDataBase = solventAtomDataBase ;
  
  iline = 1 ;
  
  iload = 0 ;
  
  rewind( psolventITP ) ;
  
  fsearch( psolventITP , "atoms" ) ;
  
  fskip( psolventITP , 1 );
  
  while( ( itpinfo = freadline( buffer , MAXCHARINLINE , psolventITP , ';' ) ) != 0 )
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
      else if( ( tmp_char = getfirst( buffer ) ) == '[' )
      {
        //printf("\nSTARTING OF NEW DIRECTIVE ... END READING ... \n\n");
        
        break ;
      }
      else
      {
        //printf("\nLine reads : %s ...\n" , buffer );
        
        iload ++ ;
        
        strpickword( buffer , 1 , cache );
        
        solventAtomDataBase[ iload -1 ].nr = atoi( cache ) ;
        
        strpickword( buffer , 2 , cache );
        
        strcpy( solventAtomDataBase[ iload -1 ].type , cache );
        
        strpickword( buffer , 3 , cache );
        
        solventAtomDataBase[ iload -1 ].resnr = atoi( cache ) ;
        
        strpickword( buffer , 4 , cache );
        
        strcpy( solventAtomDataBase[ iload -1 ].residue , cache );
        
        strpickword( buffer , 5 , cache );
        
        strcpy( solventAtomDataBase[ iload -1 ].atomname , cache );
        
        strpickword( buffer , 6 , cache );
        
        solventAtomDataBase[ iload -1 ].cgnr = atoi( cache ) ;
        
        strpickword( buffer , 7 , cache );
        
        solventAtomDataBase[ iload -1 ].charge = atof( cache ) ;
        
        strpickword( buffer , 8 , cache );
        
        solventAtomDataBase[ iload -1 ].mass = atof( cache ) ;
        
        //printf("\nCharge of this atom is %lf ...\n" , soluteAtomDataBase[ iload -1 ].charge ) ;
 
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

  
  // -------------------------------> Reading both ITP File For vdW Information ... <---------------------------------- //
  
  
  // =====> Pre-Loading ... 
  
  rewind( psoluteITP ) ;
  
  int natomSoluteTypes ; 
  
  printf("\n---> Now let's pre-load the .itp file for SOLUTE again and see how many atoms types it is describing ... <---\n");
  
  iline = 1 ;
  
  iload = 0 ;
  
  rewind( psoluteITP ) ;
  
  fsearch( psoluteITP , "atomtypes" ) ;
  
  fskip( psoluteITP , 1 );
  
  while( ( itpinfo = freadline( buffer , MAXCHARINLINE , psoluteITP , ';' ) ) != 0 )
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
      
      if( ( tmp_char = getfirst( buffer ) ) == ';' )
      {
        printf("\nThis is a comment line ... So nothing will be loaded ...\n");
        
        //fskip( pinputfile , 1 );
      }
      else if( ( tmp_char = getfirst( buffer ) ) == '[' )
      {
        printf("\nSTARTING OF NEW DIRECTIVE ... END READING ... \n\n");
        
        break ;
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

  natomSoluteTypes = iload ;
  
  printf("\nThere are %d atomtypes described in itp file ...\n" , natomSoluteTypes );



  // =====> Loading ... 
  
  printf("\n---> Now let's actually load the .itp file for SOLUTE ... <---\n\n\n");
  
  TOP soluteAtomTypes[ natomSoluteTypes ];
  
  TOP * psoluteAtomTypes ; psoluteAtomTypes = soluteAtomTypes ;
  
  iline = 1 ;
  
  iload = 0 ;
  
  rewind( psoluteITP ) ;
  
  fsearch( psoluteITP , "atomtypes" ) ;
  
  fskip( psoluteITP , 1 );
  
  while( ( itpinfo = freadline( buffer , MAXCHARINLINE , psoluteITP , ';' ) ) != 0 )
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
      else if( ( tmp_char = getfirst( buffer ) ) == '[' )
      {
        //printf("\nSTARTING OF NEW DIRECTIVE ... END READING ... \n\n");
        
        break ;
      }
      else
      {
        //printf("\nLine reads : %s ...\n" , buffer );
        
        iload ++ ;
        
        strpickword( buffer , 1 , cache );
        
        strcpy( soluteAtomTypes[ iload -1 ].typename , cache );
        
        strpickword( buffer , 2 , cache );
        
        strcpy( soluteAtomTypes[ iload -1 ].bondtype , cache );
        
        strpickword( buffer , 3 , cache );
        
        soluteAtomTypes[ iload -1 ].mass = atof( cache ) ;
        
        strpickword( buffer , 4 , cache );
        
        soluteAtomTypes[ iload -1 ].charge = atof( cache ) ;
        
        strpickword( buffer , 5 , cache );
        
        strcpy( soluteAtomTypes[ iload -1 ].ptype , cache );
        
        strpickword( buffer , 6 , cache );
        
        soluteAtomTypes[ iload -1 ].sigma = atof( cache ) ;
        
        strpickword( buffer , 7 , cache );
        
        soluteAtomTypes[ iload -1 ].epsilon = atof( cache ) ;
                
        //printf("\nCharge of this atom is %lf ...\n" , soluteAtomDataBase[ iload -1 ].charge ) ;
 
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
  
  
  /*
  // =====> Pre-Loading ... Solvent
  
  int natomSolventTypes ; 
  
  printf("\n---> Now let's pre-load the .itp file for SOLVENT again and see how many atoms types it is describing ... <---\n");
  
  iline = 1 ;
  
  iload = 0 ;
  
  rewind( psolventITP ) ;
  
  fsearch( psolventITP , "atoms" ) ;
  
  fskip( psolventITP , 1 );
  
  while( ( itpinfo = freadline( buffer , MAXCHARINLINE , psolventITP , ';' ) ) != 0 )
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
      
      if( ( tmp_char = getfirst( buffer ) ) == ';' )
      {
        printf("\nThis is a comment line ... So nothing will be loaded ...\n");
        
        //fskip( pinputfile , 1 );
      }
      else if( ( tmp_char = getfirst( buffer ) ) == '[' )
      {
        printf("\nSTARTING OF NEW DIRECTIVE ... END READING ... \n\n");
        
        break ;
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

  natomSolventTypes = iload ;
  
  printf("\nThere are %d atomtypes described in itp file ...\n" , natomSolventTypes );



  // =====> Loading ... 
  
  printf("\n---> Now let's actually load the .itp file for SOLVENT ... <---\n\n\n");
  
  TOP solventAtomTypes[ natomSolventTypes ];
  
  TOP * psolventAtomTypes ; psolventAtomTypes = solventAtomTypes ;
  
  iline = 1 ;
  
  iload = 0 ;
  
  rewind( psolventITP ) ;
  
  fsearch( psolventITP , "atoms" ) ;
  
  fskip( psolventITP , 1 );
  
  while( ( itpinfo = freadline( buffer , MAXCHARINLINE , psolventITP , ';' ) ) != 0 )
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
      else if( ( tmp_char = getfirst( buffer ) ) == '[' )
      {
        //printf("\nSTARTING OF NEW DIRECTIVE ... END READING ... \n\n");
        
        break ;
      }
      else
      {
        //printf("\nLine reads : %s ...\n" , buffer );
        
        iload ++ ;
        
        strpickword( buffer , 1 , cache );
        
        strcpy( solventAtomTypes[ iload -1 ].typename , cache );
        
        strpickword( buffer , 2 , cache );
        
        strcpy( solventAtomTypes[ iload -1 ].bondtype , cache );
        
        strpickword( buffer , 3 , cache );
        
        solventAtomTypes[ iload -1 ].mass = atof( cache ) ;
        
        strpickword( buffer , 4 , cache );
        
        solventAtomTypes[ iload -1 ].charge = atof( cache ) ;
        
        strpickword( buffer , 5 , cache );
        
        strcpy( solventAtomTypes[ iload -1 ].ptype , cache );
        
        strpickword( buffer , 6 , cache );
        
        solventAtomTypes[ iload -1 ].sigma = atof( cache ) ;
        
        strpickword( buffer , 7 , cache );
        
        solventAtomTypes[ iload -1 ].epsilon = atof( cache ) ;
                
        //printf("\nCharge of this atom is %lf ...\n" , soluteAtomDataBase[ iload -1 ].charge ) ;
 
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
  
  
  
  
  // -------------------------------> Read index file and acquire info of P and Q ... <---------------------------------- //
  
  int next_atom_index , current_atom_index ;
  
  int natomPGroup , natomQGroup ;
  
  int ndxinfo , fsearchFlag ;
  
  // =====> P Group , Pre-Loading 
  
  rewind( pindex ) ;
  
  if( ( fsearchFlag = fsearch( pindex , pGroupName ) ) == 0 ) 
  {
    printf("\nERROR : User specified group [ %s ] NOT FOUND in ndx file ...\n\n" , pGroupName ) ;
    
    exit( 39 ) ;
  }
  
  fskip( pindex , 1 ) ;
  
  fscanf( pindex , "%s" , tmpString );
  
  //printf("\nFor this group , the 1st grabbed tmpstring is %s ... \n" , tmpstring );
    
  //while( strcmp( tmpstring , "[") != 0 )
   
  for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
  {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
    if( strcmp( tmpString , "-" ) == 0 )
    {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
      fscanf( pindex , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
        
      itmp = itmp + ( next_atom_index - current_atom_index );
        
      current_atom_index = next_atom_index ;
        
      if( ( ndxinfo = fscanf( pindex , "%s" , tmpString ) ) == EOF )
      {
         break ;  
      }

    }
    else
    {
      itmp ++ ;
        
      current_atom_index = atoi( tmpString );
        
      //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
      fscanf( pindex , "%s" , tmpString );
    }

    

  }

  natomPGroup = itmp ;
  
  printf("\nBased on index file , P group %s has %d atoms ...\n" , pGroupName , natomPGroup );


  int * pGroupAtoms = calloc( natomPGroup , sizeof( int ) ) ;
  
  rewind( pindex ) ;
  
  fsearch( pindex , pGroupName ) ;
  
  fskip( pindex , 1 ) ;
  
  fscanf( pindex , "%s" , tmpString );

  for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
  {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
    if( strcmp( tmpString , "-" ) == 0 )
    {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
      fscanf( pindex , "%d" , &next_atom_index );
        
      //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
      
      for( iatom = current_atom_index + 1 ; iatom <= next_atom_index ; iatom ++ )
      {
        *( pGroupAtoms + itmp + iatom - current_atom_index - 1 ) = iatom ;

      }
        
      itmp = itmp + ( next_atom_index - current_atom_index );
        
      current_atom_index = next_atom_index ;
        
      if( ( ndxinfo = fscanf( pindex , "%s" , tmpString ) ) == EOF )
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
        
      if( ( ndxinfo = fscanf( pindex , "%s" , tmpString ) ) == EOF )
      {
         break ;  
      }

    }

    

  }


  printf("\nFinish loading information for group [ %s ] ...\n" , pGroupName );

  // =====> Q Group , Pre-Loading 
  
  rewind( pindex ) ;
  
  
  if( ( fsearchFlag = fsearch( pindex , qGroupName ) ) == 0 ) 
  {
    printf("\nERROR : User specified group [ %s ] NOT FOUND in ndx file ...\n\n" , qGroupName ) ;
    
    exit( 39 ) ;
  }
  
  
  fskip( pindex , 1 ) ;
  
  fscanf( pindex , "%s" , tmpString );
  
  //printf("\nFor this group , the 1st grabbed tmpstring is %s ... \n" , tmpstring );
    
  //while( strcmp( tmpstring , "[") != 0 )
   
  for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
  {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
    if( strcmp( tmpString , "-" ) == 0 )
    {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
      fscanf( pindex , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
        
      itmp = itmp + ( next_atom_index - current_atom_index );
        
      current_atom_index = next_atom_index ;
        
      if( ( ndxinfo = fscanf( pindex , "%s" , tmpString ) ) == EOF )
      {
         break ;  
      }

    }
    else
    {
      itmp ++ ;
        
      current_atom_index = atoi( tmpString );
        
      //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
      if( ( ndxinfo = fscanf( pindex , "%s" , tmpString ) ) == EOF )
      {
         break ;  
      }

    }

    

  }

  natomQGroup = itmp ;
  
  printf("\nBased on index file , Q group %s has %d atoms ...\n" , qGroupName , natomQGroup );


  int * qGroupAtoms = calloc( natomQGroup , sizeof( int ) ) ;
  
  rewind( pindex ) ;
  
  fsearch( pindex , qGroupName ) ;
  
  fskip( pindex , 1 ) ;
  
  fscanf( pindex , "%s" , tmpString );

  for( itmp = 0 ; strcmp( tmpString , "[" ) != 0 ;   )
  {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
    if( strcmp( tmpString , "-" ) == 0 )
    {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
      fscanf( pindex , "%d" , &next_atom_index );
        
      //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
      
      for( iatom = current_atom_index + 1 ; iatom <= next_atom_index ; iatom ++ )
      {
        *( qGroupAtoms + itmp + iatom - current_atom_index - 1 ) = iatom ;

      }
        
      itmp = itmp + ( next_atom_index - current_atom_index );
        
      current_atom_index = next_atom_index ;
        
      fscanf( pindex , "%s" , tmpString ) ;
    }
    else
    {   
      current_atom_index = atoi( tmpString );
      
      *( qGroupAtoms + itmp ) = current_atom_index ;
      
      itmp ++ ;
        
      //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
      if( ( ndxinfo = fscanf( pindex , "%s" , tmpString ) ) == EOF )
      {
         break ;  
      }

    }

    

  }

  printf("\nFinish loading information for group [ %s ] ...\n" , qGroupName );


  // -------------------------------> Calculate COM for P and Q in Eq... <---------------------------------- //

//double comEqP[ 3 ] , comEqQ[ 3 ] , coordinateN[ 3 ] , comNEqP[ 3 ] , comNEqQ[ 3 ] ;
//  double * comEqP , * comEqQ , * coordinateN , * comNEqP , * comNEqQ ;
  
//  double * vectorEqPQ , * vectorNEqPQ;

  comEqP = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comEqP ) ;
  
  comEqQ = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comEqQ ) ;

  comVelocityEqP = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comVelocityEqP ) ;

  comNEqP = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comNEqP ) ;
  
  comNEqQ = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comNEqQ ) ;
  
  coordinateEqN = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , coordinateEqN ) ;
  
  coordinateNEqN = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , coordinateNEqN ) ;
  
  vectorEqPQ = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , vectorEqPQ ) ;
  
  vectorNEqPQ = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , vectorNEqPQ ) ;

  vectorEqPN = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , vectorEqPN ) ;
  
  vectorNEqPN = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , vectorNEqPN ) ;




  double molMassP , molMassQ ;
  
  molMassP = 0.000 ; molMassQ = 0.000 ;
  
  for( iatom = 0 ; iatom < natomPGroup ; iatom ++ )
  {
    if( * ( pGroupAtoms + iatom ) > natomsoluteSelect )
    {
      printf("\nWe need the information about atom # %d , but you itp file does not contain that info. Mission Aborted ...\n" , * ( pGroupAtoms + iatom ) ) ;
      
      exit( 683 ) ;
    
    }
  
    molMassP = molMassP + soluteAtomDataBase[ * ( pGroupAtoms + iatom ) - 1 ].mass ;
  
    *( comEqP + 0 ) = *( comEqP + 0 ) + soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cx ;
    
    *( comEqP + 1 ) = *( comEqP + 1 ) + soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cy ;
  
    *( comEqP + 2 ) = *( comEqP + 2 ) + soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cz ;
    
    if( exVelocity == YES )
    {
      *( comVelocityEqP + 0 ) = *( comVelocityEqP + 0 ) + soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].vx ;
      
      *( comVelocityEqP + 1 ) = *( comVelocityEqP + 1 ) + soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].vy ;
    
      *( comVelocityEqP + 2 ) = *( comVelocityEqP + 2 ) + soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].vz ;    
    }
  
  }
  
  *( comEqP + 0 ) = *( comEqP + 0 ) / molMassP ;
  
  *( comEqP + 1 ) = *( comEqP + 1 ) / molMassP ;
  
  *( comEqP + 2 ) = *( comEqP + 2 ) / molMassP ;
  
  if( exVelocity == YES )
  {
    *( comVelocityEqP + 0 ) = *( comVelocityEqP + 0 ) / molMassP ;
    
    *( comVelocityEqP + 1 ) = *( comVelocityEqP + 1 ) / molMassP ;
    
    *( comVelocityEqP + 2 ) = *( comVelocityEqP + 2 ) / molMassP ;
  }

  
  printf("\nCOM coordinates of P in EqMD is % 10.8E , % 10.8E , % 10.8E ...\n\n" , *( comEqP + 0 ) , *( comEqP + 1 ) , *( comEqP + 2 ) ) ;

  if( exVelocity == YES )
  {
    printf("\nCOM velocities of P in EqMD is % 10.8E , % 10.8E , % 10.8E ...\n\n" , *( comVelocityEqP + 0 ) , *( comVelocityEqP + 1 ) , *( comVelocityEqP + 2 ) ) ;
  }



  for( iatom = 0 ; iatom < natomQGroup ; iatom ++ )
  {
    if( * ( qGroupAtoms + iatom ) > natomsoluteSelect )
    {
      printf("\nWe need the information about atom # %d , but you itp file does not contain that info. Mission Aborted ...\n" , * ( qGroupAtoms + iatom ) ) ;
      
      exit( 683 ) ;
    
    }

    molMassQ = molMassQ + soluteAtomDataBase[ * ( qGroupAtoms + iatom ) - 1 ].mass ;
  
    *( comEqQ + 0 ) = *( comEqQ + 0 ) + soluteAtomDataBase[ *( qGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cx ;
    
    *( comEqQ + 1 ) = *( comEqQ + 1 ) + soluteAtomDataBase[ *( qGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cy ;
  
    *( comEqQ + 2 ) = *( comEqQ + 2 ) + soluteAtomDataBase[ *( qGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cz ;
  
  }
  
  *( comEqQ + 0 ) = *( comEqQ + 0 ) / molMassQ ;
  
  *( comEqQ + 1 ) = *( comEqQ + 1 ) / molMassQ ;
  
  *( comEqQ + 2 ) = *( comEqQ + 2 ) / molMassQ ;

  printf("\nCOM coordinates of Q in EqMD is % 10.8E , % 10.8E , % 10.8E ...\n\n" , *( comEqQ + 0 ) , *( comEqQ + 1 ) , *( comEqQ + 2 ) ) ;

  
  *( vectorEqPQ + 0 ) = *( comEqQ + 0 ) - *( comEqP + 0 ) ;
  
  *( vectorEqPQ + 1 ) = *( comEqQ + 1 ) - *( comEqP + 1 ) ;
  
  *( vectorEqPQ + 2 ) = *( comEqQ + 2 ) - *( comEqP + 2 ) ;



  
  *( coordinateEqN + 0 ) = EqAtomList[ referenceAtomNumber - 1 ].cx ;
  
  *( coordinateEqN + 1 ) = EqAtomList[ referenceAtomNumber - 1 ].cy ;
  
  *( coordinateEqN + 2 ) = EqAtomList[ referenceAtomNumber - 1 ].cz ;
  
  
  *( vectorEqPN + 0 ) = ( *( coordinateEqN + 0 ) ) - ( *( comEqP + 0 ) ) ;
  
  *( vectorEqPN + 1 ) = ( *( coordinateEqN + 1 ) ) - ( *( comEqP + 1 ) ) ;
  
  *( vectorEqPN + 2 ) = ( *( coordinateEqN + 2 ) ) - ( *( comEqP + 2 ) ) ;
  

  

  double * crdFitEqSolute = calloc( 3 * natomsoluteSelect , sizeof( double ) ) ;
  
  double * velFitEqSolute = calloc( 3 * natomsoluteSelect , sizeof( double ) ) ;
  
  dzeros( natomsoluteSelect , 3 , crdFitEqSolute ) ;
    
  dzeros( natomsoluteSelect , 3 , velFitEqSolute ) ;
  
  
  for( iatom = 0 ; iatom < natomsoluteSelect ; iatom ++ )
  {
    *( crdFitEqSolute + 3 * iatom + 0 ) = EqAtomList[ iatom ].cx - *( comEqP + 0 ) ;
    
    *( crdFitEqSolute + 3 * iatom + 1 ) = EqAtomList[ iatom ].cy - *( comEqP + 1 ) ;
    
    *( crdFitEqSolute + 3 * iatom + 2 ) = EqAtomList[ iatom ].cz - *( comEqP + 2 ) ;
    
    if( exVelocity == YES )
    {
      *( velFitEqSolute + 3 * iatom + 0 ) = EqAtomList[ iatom ].vx - *( comVelocityEqP + 0 ) ;
    
      *( velFitEqSolute + 3 * iatom + 1 ) = EqAtomList[ iatom ].vy - *( comVelocityEqP + 1 ) ;
    
      *( velFitEqSolute + 3 * iatom + 2 ) = EqAtomList[ iatom ].vz - *( comVelocityEqP + 2 ) ;
    
    }

  }



  // -------------------------------> Loading crd and Calculate COM for NEq ... <---------------------------------- //

  int length_crd_file = flength( pinpcrd ) ;
  
  int length_vel_file ;
  
  if( velDump == YES )
  {
    length_vel_file = flength( pinpvel ) ;
  }
  else
  {
    length_vel_file = 3 * natomsoluteSelect ;
  }
  
  if( length_crd_file != 3 * natomsoluteSelect || length_vel_file != 3 * natomsoluteSelect )
  {
    printf("\nPlease check your crd and vel file . The length is not the # of Cart of solute as user-defined ..\n");
  
    exit( 1094 );
  }

  double * crdNEq , * velNEq ;
  
  crdNEq = ( double * ) calloc( 3 * natomsoluteSelect , sizeof( double ) ) ;
  
  if( velDump == YES )
  {
    velNEq = ( double * ) calloc( 3 * natomsoluteSelect , sizeof( double ) ) ;
  }
  
  rewind( pinpcrd ) ; 
  
  fload( pinpcrd , crdNEq ) ;
  
  if( velDump == YES )
  {
    rewind( pinpvel ) ;
    
    fload( pinpvel , velNEq ) ;
  }
  
  
  
  
  
  comVelocityNEqP = calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comVelocityNEqP ) ;
  
  for( iatom = 0 ; iatom < natomPGroup ; iatom ++ )
  {  
    if( * ( pGroupAtoms + iatom ) > natomsoluteSelect )
    {
      printf("\nWe need the information about atom # %d , but you itp file does not contain that info. Mission Aborted ...\n" , * ( pGroupAtoms + iatom ) ) ;
      
      exit( 683 ) ;
    
    }
    
    *( comNEqP + 0 ) = *( comNEqP + 0 ) + soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass * ( *( crdNEq + 3 * ( *( pGroupAtoms + iatom ) - 1 ) + 0 ) ) ;
    
    *( comNEqP + 1 ) = *( comNEqP + 1 ) + soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass * ( *( crdNEq + 3 * ( *( pGroupAtoms + iatom ) - 1 ) + 1 ) ) ;
  
    *( comNEqP + 2 ) = *( comNEqP + 2 ) + soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass * ( *( crdNEq + 3 * ( *( pGroupAtoms + iatom ) - 1 ) + 2 ) ) ;
    
    if( velDump == YES )
    {
      *( comVelocityNEqP + 0 ) = *( comVelocityNEqP + 0 ) + soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass * ( *( velNEq + 3 * iatom + 0 ) ) ;
      
      *( comVelocityNEqP + 1 ) = *( comVelocityNEqP + 1 ) + soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass * ( *( velNEq + 3 * iatom + 1 ) ) ;
    
      *( comVelocityNEqP + 2 ) = *( comVelocityNEqP + 2 ) + soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass * ( *( velNEq + 3 * iatom + 2 ) ) ;
    }
  
  }
  
  *( comNEqP + 0 ) = *( comNEqP + 0 ) / molMassP ;
  
  *( comNEqP + 1 ) = *( comNEqP + 1 ) / molMassP ;
  
  *( comNEqP + 2 ) = *( comNEqP + 2 ) / molMassP ;
  
  if( velDump == YES )
  {
    *( comVelocityNEqP + 0 ) = *( comVelocityNEqP + 0 ) / molMassP ;
    
    *( comVelocityNEqP + 1 ) = *( comVelocityNEqP + 1 ) / molMassP ;
    
    *( comVelocityNEqP + 2 ) = *( comVelocityNEqP + 2 ) / molMassP ;
  }



  for( iatom = 0 ; iatom < natomQGroup ; iatom ++ )
  {
    if( * ( qGroupAtoms + iatom ) > natomsoluteSelect )
    {
      printf("\nWe need the information about atom # %d , but you itp file does not contain that info. Mission Aborted ...\n" , * ( qGroupAtoms + iatom ) ) ;
      
      exit( 683 ) ;
    
    }
    
    *( comNEqQ + 0 ) = *( comNEqQ + 0 ) + soluteAtomDataBase[ *( qGroupAtoms + iatom ) - 1 ].mass * ( *( crdNEq + 3 * ( *( qGroupAtoms + iatom ) - 1 ) + 0 ) ) ;
    
    *( comNEqQ + 1 ) = *( comNEqQ + 1 ) + soluteAtomDataBase[ *( qGroupAtoms + iatom ) - 1 ].mass * ( *( crdNEq + 3 * ( *( qGroupAtoms + iatom ) - 1 ) + 1 ) ) ;
  
    *( comNEqQ + 2 ) = *( comNEqQ + 2 ) + soluteAtomDataBase[ *( qGroupAtoms + iatom ) - 1 ].mass * ( *( crdNEq + 3 * ( *( qGroupAtoms + iatom ) - 1 ) + 2 ) ) ;
  
  }
  
  *( comNEqQ + 0 ) = *( comNEqQ + 0 ) / molMassQ ;
  
  *( comNEqQ + 1 ) = *( comNEqQ + 1 ) / molMassQ ;
  
  *( comNEqQ + 2 ) = *( comNEqQ + 2 ) / molMassQ ;


  
  *( coordinateNEqN + 0 ) = *( crdNEq + 3 * ( referenceAtomNumber - 1 ) + 0 ) ;  
  
  *( coordinateNEqN + 1 ) = *( crdNEq + 3 * ( referenceAtomNumber - 1 ) + 1 ) ;  
  
  *( coordinateNEqN + 2 ) = *( crdNEq + 3 * ( referenceAtomNumber - 1 ) + 2 ) ;  
  
  
  *( vectorNEqPN + 0 ) = ( *( coordinateNEqN + 0 ) ) - ( *( comNEqP + 0 ) ) ;
  
  *( vectorNEqPN + 1 ) = ( *( coordinateNEqN + 1 ) ) - ( *( comNEqP + 1 ) ) ;
  
  *( vectorNEqPN + 2 ) = ( *( coordinateNEqN + 2 ) ) - ( *( comNEqP + 2 ) ) ;
  
  printf("\nIn NEq , PN vector is % 10.8E\t% 10.8E\t% 10.8E" , *( vectorNEqPN + 0 )  , *( vectorNEqPN + 1 )  , *( vectorNEqPN + 2 )  ) ;


  double * crdFitNEqSolute = calloc( 3 * natomsoluteSelect , sizeof( double ) ) ;

  double * velFitNEqSolute = calloc( 3 * natomsoluteSelect , sizeof( double ) ) ;
  
  dzeros( natomsoluteSelect , 3 , crdFitNEqSolute ) ;
    
  dzeros( natomsoluteSelect , 3 , velFitNEqSolute ) ;
  
  
  
  for( iatom = 0 ; iatom < natomsoluteSelect ; iatom ++ )
  {
    *( crdFitNEqSolute + 3 * iatom + 0 ) = *( crdNEq + 3 * iatom + 0 ) - *( comNEqP + 0 ) ;
    
    *( crdFitNEqSolute + 3 * iatom + 1 ) = *( crdNEq + 3 * iatom + 1 ) - *( comNEqP + 1 ) ;
    
    *( crdFitNEqSolute + 3 * iatom + 2 ) = *( crdNEq + 3 * iatom + 2 ) - *( comNEqP + 2 ) ;
    
    if( velDump == YES )
    {
      *( velFitNEqSolute + 3 * iatom + 0 ) = *( velNEq + 3 * iatom + 0 ) ; //- *( comVelocityNEqP + 0 ) ;
    
      *( velFitNEqSolute + 3 * iatom + 1 ) = *( velNEq + 3 * iatom + 1 ) ; //- *( comVelocityNEqP + 1 ) ;
    
      *( velFitNEqSolute + 3 * iatom + 2 ) = *( velNEq + 3 * iatom + 2 ) ; //- *( comVelocityNEqP + 2 ) ;
    
    }
  }


  *( vectorNEqPQ + 0 ) = *( comNEqQ + 0 ) - *( comNEqP + 0 ) ;
  
  *( vectorNEqPQ + 1 ) = *( comNEqQ + 1 ) - *( comNEqP + 1 ) ;
  
  *( vectorNEqPQ + 2 ) = *( comNEqQ + 2 ) - *( comNEqP + 2 ) ;






  
/*
  *( coordinateEqN + 0 ) = EqAtomList[ referenceAtomNumber - 1 ].cx ;
  
  *( coordinateEqN + 1 ) = EqAtomList[ referenceAtomNumber - 1 ].cy ;
  
  *( coordinateEqN + 2 ) = EqAtomList[ referenceAtomNumber - 1 ].cz ;
*/








  // -------------------------------> Pausing ... debugging outputs ... <---------------------------------- //

  /*

  debug = fopen( "pEqCoordinates.deb" , "wb+") ;
  
  for( iatom = 0 ; iatom < natomPGroup ; iatom ++ )
  {
    fprintf( debug , "%s\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n" , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].atomname , soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cx , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cy , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cz ) ;

  }
  
  fprintf( debug , "COM\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n\n" , molMassP , *( comEqP + 0 ) , *( comEqP + 1 ) , *( comEqP + 2 ) ) ;
  
  fclose( debug ) ;




  

  debug = fopen( "qEqCoordinates.deb" , "wb+") ;
  
  for( iatom = 0 ; iatom < natomQGroup ; iatom ++ )
  {
    fprintf( debug , "%s\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n" , EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].atomname , soluteAtomDataBase[ *( qGroupAtoms + iatom ) - 1 ].mass , EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cx , EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cy , EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cz ) ;

  }
  
  fprintf( debug , "COM\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n" , molMassQ , *( comEqQ + 0 ) , *( comEqQ + 1 ) , *( comEqQ + 2 ) ) ;
  
  fclose( debug ) ;
    
  
  
  
  
  debug = fopen( "pNEqCoordinates.deb" , "wb+") ;
  
  for( iatom = 0 ; iatom < natomPGroup ; iatom ++ )
  {
    fprintf( debug , "%s\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n" , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].atomname , soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass , *( crdNEq + 3 * ( *( pGroupAtoms + iatom ) - 1 ) + 0 ) , *( crdNEq + 3 * ( *( pGroupAtoms + iatom ) - 1 ) + 1 ) , *( crdNEq + 3 * ( *( pGroupAtoms + iatom ) - 1 ) + 2 ) ) ;

  }
  
  fprintf( debug , "COM\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n" , molMassP , *( comNEqP + 0 ) , *( comNEqP + 1 ) , *( comNEqP + 2 ) ) ;
  
  fclose( debug ) ;





  debug = fopen( "qNEqCoordinates.deb" , "wb+") ;

  for( iatom = 0 ; iatom < natomQGroup ; iatom ++ )
  {
    fprintf( debug , "%s\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n" , EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].atomname , soluteAtomDataBase[ *( qGroupAtoms + iatom ) - 1 ].mass , *( crdNEq + 3 * ( *( qGroupAtoms + iatom ) - 1 ) + 0 ) , *( crdNEq + 3 * ( *( qGroupAtoms + iatom ) - 1 ) + 1 ) , *( crdNEq + 3 * ( *( qGroupAtoms + iatom ) - 1 ) + 2 ) ) ;

  }
  
  fprintf( debug , "COM\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n" , molMassQ , *( comNEqQ + 0 ) , *( comNEqQ + 1 ) , *( comNEqQ + 2 ) ) ;
  
  fclose( debug ) ;
  

  
  if( velDump == YES )
  {
    
    
    debug = fopen( "pEqVelocities.deb" , "wb+") ;
    
    for( iatom = 0 ; iatom < natomPGroup ; iatom ++ )
    {
      fprintf( debug , "%s\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n" , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].atomname , soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].vx , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].vy , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].vz ) ;

    }
  
    fprintf( debug , "COM\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n\n" , molMassP , *( comVelocityEqP + 0 ) , *( comVelocityEqP + 1 ) , *( comVelocityEqP + 2 ) ) ;
  
    fclose( debug ) ;
  
  
    debug = fopen( "pNEqVelocities.deb" , "wb+") ;
    
    for( iatom = 0 ; iatom < natomPGroup ; iatom ++ )
    {
      fprintf( debug , "%s\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n" , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].atomname , soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass , *( velNEq + 3 * ( *( pGroupAtoms + iatom ) - 1 ) + 0 ) , *( velNEq + 3 * ( *( pGroupAtoms + iatom ) - 1 ) + 1 ) , *( velNEq + 3 * ( *( pGroupAtoms + iatom ) - 1 ) + 2 ) ) ;

    }
  
    fprintf( debug , "COM\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n\n" , molMassP , *( comVelocityNEqP + 0 ) , *( comVelocityNEqP + 1 ) , *( comVelocityNEqP + 2 ) ) ;
  
    fclose( debug ) ;
  
  
    debug = fopen( "velFitNEqSolute.deb" , "wb+" ) ;
    
    for( iatom = 0 ; iatom < natomPGroup ; iatom ++ )
    {
      fprintf( debug , "%s\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n" , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].atomname , soluteAtomDataBase[ *( pGroupAtoms + iatom ) - 1 ].mass , *( velFitNEqSolute + 3 * ( *( pGroupAtoms + iatom ) - 1 ) + 0 ) , *( velFitNEqSolute + 3 * ( *( pGroupAtoms + iatom ) - 1 ) + 1 ) , *( velFitNEqSolute + 3 * ( *( pGroupAtoms + iatom ) - 1 ) + 2 ) ) ;

    }
    
    fclose( debug ) ;
  
  
  }


  */



  // -------------------------------> Building Eq.Rotation Matrices ... <---------------------------------- //

  double a , b , c , rab , rabc ;
  
  double * Ms1 , * Ms2 , * Ms3 , * tmpM , * M1 ;
  
  Ms1 = ( double * ) calloc( 9 , sizeof( double ) ) ; dzeros( 3 , 3 , Ms1 ) ;
  
  Ms2 = ( double * ) calloc( 9 , sizeof( double ) ) ; dzeros( 3 , 3 , Ms2 ) ;
  
  Ms3 = ( double * ) calloc( 9 , sizeof( double ) ) ; dzeros( 3 , 3 , Ms3 ) ;
  
  tmpM = ( double * ) calloc( 9 , sizeof( double ) ) ; dzeros( 3 , 3 , tmpM ) ;

  M1 = ( double * ) calloc( 9 , sizeof( double ) ) ; dzeros( 3 , 3 , M1 ) ;


  a = *( vectorEqPQ + 0 ) ; b = *( vectorEqPQ + 1 ) ; c = *( vectorEqPQ + 2 ) ;
  
  rab = sqrt( a * a + b * b ) ;
  
  rabc = sqrt( a * a + b * b + c * c ) ;
  
  printf("\nOkay ... some parameters ... :\na = % 10.8f , b = % 10.8f , c = % 10.8f \n" , a , b , c ) ;
  
  printf("\nrab = % 10.8f , rabc = % 10.8f \n" , rab , rabc ) ;
  
  *( Ms1 + 0 ) = a / rab          ;   *( Ms1 + 1 ) = b / rab ;   *( Ms1 + 2 ) = 0.000 ;
  
  *( Ms1 + 3 ) = -1.000 * b / rab ;   *( Ms1 + 4 ) = a / rab ;   *( Ms1 + 5 ) = 0.000 ;
  
  *( Ms1 + 6 ) = 0.000            ;   *( Ms1 + 7 ) = 0.000   ;   *( Ms1 + 8 ) = 1.000 ;
  
  
  *( Ms2 + 0 ) = c / rabc         ;   *( Ms2 + 1 ) = 0.000   ;   *( Ms2 + 2 ) = -1.000 * rab / rabc ;
  
  *( Ms2 + 3 ) = 0.000            ;   *( Ms2 + 4 ) = 1.000   ;   *( Ms2 + 5 ) = 0.000               ;
  
  *( Ms2 + 6 ) = rab / rabc       ;   *( Ms2 + 7 ) = 0.000   ;   *( Ms2 + 8 ) = c / rabc            ;
  





  dgemm_( "T" , "T" , &ithree , &ithree , &ithree , &done , Ms2 , &ithree , Ms1 , &ithree , &dzero , tmpM , &ithree ) ;
  
  dtranspose( 3 , tmpM , tmpM ) ;
  


  //*( vectorEqPN + 0 ) 
  
  double * vectorFitEqPN = calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , vectorFitEqPN ) ;
  
  double aNFit , bNFit , cNFit , rNabFit ;
  
  *( vectorFitEqPN + 0 ) = ( *( tmpM + 0 ) ) * ( *( vectorEqPN + 0 ) ) + ( *( tmpM + 1 ) ) * ( *( vectorEqPN + 1 ) ) + ( *( tmpM + 2 ) ) * ( *( vectorEqPN + 2 ) ) ; 

  *( vectorFitEqPN + 1 ) = ( *( tmpM + 3 ) ) * ( *( vectorEqPN + 0 ) ) + ( *( tmpM + 4 ) ) * ( *( vectorEqPN + 1 ) ) + ( *( tmpM + 5 ) ) * ( *( vectorEqPN + 2 ) ) ; 

  *( vectorFitEqPN + 2 ) = ( *( tmpM + 6 ) ) * ( *( vectorEqPN + 0 ) ) + ( *( tmpM + 7 ) ) * ( *( vectorEqPN + 1 ) ) + ( *( tmpM + 8 ) ) * ( *( vectorEqPN + 2 ) ) ; 


  aNFit = *( vectorFitEqPN + 0 ) ; 
  
  bNFit = *( vectorFitEqPN + 1 ) ; 

  cNFit = *( vectorFitEqPN + 2 ) ; 

  rNabFit = sqrt( aNFit * aNFit + bNFit * bNFit ) ;
  
  
  
  *( Ms3 + 0 ) = aNFit / rNabFit          ;   *( Ms3 + 1 ) = bNFit / rNabFit ;   *( Ms3 + 2 ) = 0.000 ;
  
  *( Ms3 + 3 ) = -1.000 * bNFit / rNabFit ;   *( Ms3 + 4 ) = aNFit / rNabFit ;   *( Ms3 + 5 ) = 0.000 ;
  
  *( Ms3 + 6 ) = 0.000                    ;   *( Ms3 + 7 ) = 0.000           ;   *( Ms3 + 8 ) = 1.000 ;
  
  
  dgemm_( "T" , "T" , &ithree , &ithree , &ithree , &done , Ms3 , &ithree , tmpM , &ithree , &dzero , M1 , &ithree ) ;
  
  dtranspose( 3 , M1 , M1 ) ;
  
  

  /*
  debug = fopen("EqMs.deb" , "wb+") ;
  
  doutput( debug , 3 , 3 , Ms1 ) ;
  
  fprintf( debug , "\n\n\n\n\n" ) ;
  
  doutput( debug , 3 , 3 , Ms2 ) ;
  
  fprintf( debug , "\n\n\n\n\n" ) ;
  
  doutput( debug , 3 , 3 , Ms3 ) ;
  
  fprintf( debug , "\n\n\n\n\n" ) ;
  
  doutput( debug , 3 , 3 , M1 ) ;
  
  fclose( debug ) ;
  
  */













  // -------------------------------> Building NEq.Rotation Matrices ... <---------------------------------- //

  //double a , b , c , rab , rabc ;
  
  //double * Ms1 , * Ms2 , * Ms3 , * tmpM , * M1 ;
  
  //Ms1 = ( double * ) calloc( 9 , sizeof( double ) ) ; 
  
  double * M2 = ( double * ) calloc( 9 , sizeof( double ) ) ; 
  
  dzeros( 3 , 3 , M2 ) ;
  
  dzeros( 3 , 3 , Ms1 ) ;
  
  //Ms2 = ( double * ) calloc( 9 , sizeof( double ) ) ; 
  
  dzeros( 3 , 3 , Ms2 ) ;
  
  //Ms3 = ( double * ) calloc( 9 , sizeof( double ) ) ; 
  
  dzeros( 3 , 3 , Ms3 ) ;
  
  //tmpM = ( double * ) calloc( 9 , sizeof( double ) ) ; 
  
  dzeros( 3 , 3 , tmpM ) ;

  //M1 = ( double * ) calloc( 9 , sizeof( double ) ) ; 
  



  a = *( vectorNEqPQ + 0 ) ; b = *( vectorNEqPQ + 1 ) ; c = *( vectorNEqPQ + 2 ) ;
  
  rab = sqrt( a * a + b * b ) ;
  
  rabc = sqrt( a * a + b * b + c * c ) ;
  
  printf("\nOkay ... some parameters ... :\na = % 10.8f , b = % 10.8f , c = % 10.8f \n" , a , b , c ) ;
  
  printf("\nrab = % 10.8f , rabc = % 10.8f \n" , rab , rabc ) ;
  
  *( Ms1 + 0 ) = a / rab          ;   *( Ms1 + 1 ) = b / rab ;   *( Ms1 + 2 ) = 0.000 ;
  
  *( Ms1 + 3 ) = -1.000 * b / rab ;   *( Ms1 + 4 ) = a / rab ;   *( Ms1 + 5 ) = 0.000 ;
  
  *( Ms1 + 6 ) = 0.000            ;   *( Ms1 + 7 ) = 0.000   ;   *( Ms1 + 8 ) = 1.000 ;
  
  
  *( Ms2 + 0 ) = c / rabc         ;   *( Ms2 + 1 ) = 0.000   ;   *( Ms2 + 2 ) = -1.000 * rab / rabc ;
  
  *( Ms2 + 3 ) = 0.000            ;   *( Ms2 + 4 ) = 1.000   ;   *( Ms2 + 5 ) = 0.000               ;
  
  *( Ms2 + 6 ) = rab / rabc       ;   *( Ms2 + 7 ) = 0.000   ;   *( Ms2 + 8 ) = c / rabc            ;
  





  dgemm_( "T" , "T" , &ithree , &ithree , &ithree , &done , Ms2 , &ithree , Ms1 , &ithree , &dzero , tmpM , &ithree ) ;
  
  dtranspose( 3 , tmpM , tmpM ) ;
  
  
/*
  debug = fopen("NEqtmpM.deb" , "wb+") ;
  
  doutput( debug , 3 , 3 , tmpM ) ;
  
  fclose( debug ) ;
*/  
  

  //*( vectorEqPN + 0 ) 
  
  double * vectorFitNEqPN = calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , vectorFitNEqPN ) ;
  

  
  *( vectorFitNEqPN + 0 ) = ( *( tmpM + 0 ) ) * ( *( vectorNEqPN + 0 ) ) + ( *( tmpM + 1 ) ) * ( *( vectorNEqPN + 1 ) ) + ( *( tmpM + 2 ) ) * ( *( vectorNEqPN + 2 ) ) ; 

  *( vectorFitNEqPN + 1 ) = ( *( tmpM + 3 ) ) * ( *( vectorNEqPN + 0 ) ) + ( *( tmpM + 4 ) ) * ( *( vectorNEqPN + 1 ) ) + ( *( tmpM + 5 ) ) * ( *( vectorNEqPN + 2 ) ) ; 

  *( vectorFitNEqPN + 2 ) = ( *( tmpM + 6 ) ) * ( *( vectorNEqPN + 0 ) ) + ( *( tmpM + 7 ) ) * ( *( vectorNEqPN + 1 ) ) + ( *( tmpM + 8 ) ) * ( *( vectorNEqPN + 2 ) ) ; 


  aNFit = *( vectorFitNEqPN + 0 ) ; 
  
  bNFit = *( vectorFitNEqPN + 1 ) ; 

  cNFit = *( vectorFitNEqPN + 2 ) ; 

  rNabFit = sqrt( aNFit * aNFit + bNFit * bNFit ) ;
  
  printf("\nOkay ... some parameters ... :\naNFit = % 10.8f , bNFit = % 10.8f , cNFit = % 10.8f \n" , aNFit , bNFit , cNFit ) ;
  
  printf("\nrab = % 10.8f \n" , rNabFit  ) ;
  
  
  *( Ms3 + 0 ) = aNFit / rNabFit          ;   *( Ms3 + 1 ) = bNFit / rNabFit ;   *( Ms3 + 2 ) = 0.000 ;
  
  *( Ms3 + 3 ) = -1.000 * bNFit / rNabFit ;   *( Ms3 + 4 ) = aNFit / rNabFit ;   *( Ms3 + 5 ) = 0.000 ;
  
  *( Ms3 + 6 ) = 0.000                    ;   *( Ms3 + 7 ) = 0.000           ;   *( Ms3 + 8 ) = 1.000 ;
  
  
  dgemm_( "T" , "T" , &ithree , &ithree , &ithree , &done , Ms3 , &ithree , tmpM , &ithree , &dzero , M2 , &ithree ) ;
  
  dtranspose( 3 , M2 , M2 ) ;
  
  



  double * M = calloc( 9 , sizeof( double) ) ; dzeros( 3 , 3 , M ) ;
  
  dgemm_( "N" , "T" , &ithree , &ithree , &ithree , &done , M1 , &ithree , M2 , &ithree , &dzero , M , &ithree ) ;
  
  dtranspose( 3 , M , M ) ;




//   debug = fopen( "allFitVelocities.deb" , "wb+" ) ;
  
  for( iatom = 0 ; iatom < natomsoluteSelect ; iatom ++ )
  {
    *( dtmpArray + 0 ) = ( *( M + 0 ) ) * ( *( crdFitNEqSolute + 3 * iatom + 0 ) ) + ( *( M + 1 ) ) * ( *( crdFitNEqSolute + 3 * iatom + 1 ) ) + ( *( M + 2 ) ) * ( *( crdFitNEqSolute + 3 * iatom + 2 ) ) ; 

    *( dtmpArray + 1 ) = ( *( M + 3 ) ) * ( *( crdFitNEqSolute + 3 * iatom + 0 ) ) + ( *( M + 4 ) ) * ( *( crdFitNEqSolute + 3 * iatom + 1 ) ) + ( *( M + 5 ) ) * ( *( crdFitNEqSolute + 3 * iatom + 2 ) ) ; 

    *( dtmpArray + 2 ) = ( *( M + 6 ) ) * ( *( crdFitNEqSolute + 3 * iatom + 0 ) ) + ( *( M + 7 ) ) * ( *( crdFitNEqSolute + 3 * iatom + 1 ) ) + ( *( M + 8 ) ) * ( *( crdFitNEqSolute + 3 * iatom + 2 ) ) ; 

    EqAtomList[ iatom ].cx = *( dtmpArray + 0 ) + *( comEqP + 0 ) ; // Here EqAtomList is actually going to take all the transformed coordinates for solute and substitute into the corresponding position ;
    
    EqAtomList[ iatom ].cy = *( dtmpArray + 1 ) + *( comEqP + 1 ) ;
  
    EqAtomList[ iatom ].cz = *( dtmpArray + 2 ) + *( comEqP + 2 ) ;
    
    
    if( exVelocity == YES )
    {
      *( dtmpArray + 0 ) = ( *( M + 0 ) ) * ( *( velFitNEqSolute + 3 * iatom + 0 ) ) + ( *( M + 1 ) ) * ( *( velFitNEqSolute + 3 * iatom + 1 ) ) + ( *( M + 2 ) ) * ( *( velFitNEqSolute + 3 * iatom + 2 ) ) ; 
  
      *( dtmpArray + 1 ) = ( *( M + 3 ) ) * ( *( velFitNEqSolute + 3 * iatom + 0 ) ) + ( *( M + 4 ) ) * ( *( velFitNEqSolute + 3 * iatom + 1 ) ) + ( *( M + 5 ) ) * ( *( velFitNEqSolute + 3 * iatom + 2 ) ) ; 
  
      *( dtmpArray + 2 ) = ( *( M + 6 ) ) * ( *( velFitNEqSolute + 3 * iatom + 0 ) ) + ( *( M + 7 ) ) * ( *( velFitNEqSolute + 3 * iatom + 1 ) ) + ( *( M + 8 ) ) * ( *( velFitNEqSolute + 3 * iatom + 2 ) ) ; 
  
      EqAtomList[ iatom ].vx = *( dtmpArray + 0 ) ; //+ *( comVelocityEqP + 0 ) ;
      
      EqAtomList[ iatom ].vy = *( dtmpArray + 1 ) ; //+ *( comVelocityEqP + 1 ) ;
    
      EqAtomList[ iatom ].vz = *( dtmpArray + 2 ) ; //+ *( comVelocityEqP + 2 ) ;
      
      //fprintf( debug , )
    
    }
  

  }


  /*
  debug = fopen("NEqMs.deb" , "wb+") ;
  
  doutput( debug , 3 , 3 , Ms1 ) ;
  
  fprintf( debug , "\n\n\n\n\n" ) ;
  
  doutput( debug , 3 , 3 , Ms2 ) ;
  
  fprintf( debug , "\n\n\n\n\n" ) ;
  
  doutput( debug , 3 , 3 , Ms3 ) ;
  
  fprintf( debug , "\n\n\n\n\n" ) ;
  
  doutput( debug , 3 , 3 , M2 ) ;
  
  fclose( debug ) ;
  
  
  
  debug = fopen("M.deb" , "wb+") ;
  
  doutput( debug , 3 , 3 , M ) ;
  
  fclose( debug ) ;

  */

  
  
  // -------------------------------> Telling the Minimum Distance between solute and solvent molecules ... <---------------------------------- //

  int natomSolvent = natom - natomsoluteSelect ;

  double * distances = calloc( natomSolvent * natomsoluteSelect , sizeof( double ) );
  
  int pos , posSolute , posSolvent ;
  
  char solventHotAtomName[ 3 ] , soluteHotAtomName[ 3 ] ;
  
  int success = NO ; 
  
  double distThresholdScalingFactor = 1.00 ;
  
  double minimumDistance ; // which needs to be .lt. distThreshold
  
  
  if( natomsoluteSelect == natom )
  {
    success = YES ;
  }
  else
  {
    for( iatom = 0 ; iatom < natomsoluteSelect ; iatom ++ )
    {
      dtmp = 0.000 ;
      
      for( itmp = natomsoluteSelect ; itmp < natom ; itmp ++ )
      {
        dtmp = 0.000 ;
        
        dtmp = dtmp + ( EqAtomList[ iatom ].cx - EqAtomList[ itmp ].cx ) * ( EqAtomList[ iatom ].cx - EqAtomList[ itmp ].cx ) ;
      
        dtmp = dtmp + ( EqAtomList[ iatom ].cy - EqAtomList[ itmp ].cy ) * ( EqAtomList[ iatom ].cy - EqAtomList[ itmp ].cy ) ;
        
        dtmp = dtmp + ( EqAtomList[ iatom ].cz - EqAtomList[ itmp ].cz ) * ( EqAtomList[ iatom ].cz - EqAtomList[ itmp ].cz ) ;
        
        //printf("\n# %d atom , also # %d atom in solvent , coordinate is [ % 10.6f\t% 10.6f\t% 10.6f\t] \n\n" , itmp , itmp - natomsoluteSelect , EqAtomList[ itmp ].cx , EqAtomList[ itmp ].cy , EqAtomList[ itmp ].cz ) ;
        
        *( distances + iatom * natomSolvent + itmp - natomsoluteSelect ) = sqrt( dtmp ) ;
      
      }
    

    }


    /*
    debug = fopen("distances.deb" , "wb+");
    
    doutput( debug , natomsoluteSelect , natomSolvent , distances ) ;
    
    fclose( debug ) ;
    */


    
    if( distThreshold > 0.000 )
    {
      minimumDistance = dmin( natomSolvent * natomsoluteSelect , distances ) ;
    
      pos = dminID( natomSolvent * natomsoluteSelect , distances ) + 1 ; // Human Label
    
      itmp = pos % natomSolvent ;
    
      if( itmp == 0 )
      {
        itmp = itmp + natomSolvent ;
      }
    
      posSolvent = itmp + natomsoluteSelect - 1 ; // C-Label
    
      posSolute = ( pos - itmp ) / natomSolvent ; 
    
      tellsymbol( EqAtomList[ posSolute ].atomname , soluteHotAtomName ) ;
    
      tellsymbol( EqAtomList[ posSolvent ].atomname , solventHotAtomName ) ;    
    
      printf("\nAfter rotation , minimum distance between solute and solvent atoms is % 10.6f Angstrom between # %d solute atom [ %s ] and # %d solvent atom [ %s ] ... \n\n" , minimumDistance * 10 , posSolute + 1 , soluteHotAtomName , posSolvent + 1 , solventHotAtomName ) ;
    
    
    
  
    

    debug = fopen("distances.deb" , "wb+");
    
    doutput( debug , natomsoluteSelect , natomSolvent , distances ) ;
    
    fclose( debug ) ;


    
    
      if( strcmp( soluteHotAtomName , "H" ) == 0 && strcmp( solventHotAtomName , "H" ) == 0 )
      {
        distThresholdScalingFactor = 0.60;
      }
      else //if( strcmp( soluteHotAtomName , "H" ) != 0 || strcmp( solventHotAtomName , "H" ) != 0 )
      {
        distThresholdScalingFactor = 1.00 ;
      }
    
    
      distThreshold = distThreshold * distThresholdScalingFactor ;
    

    
    
      if( minimumDistance >= distThreshold )
      {
        success = YES ;
      }
      else
      {
        success = NO ;
      }
  
  
    }
    else
    {
      printf("\nNote : Universal minimum distance test was skipped ...\n") ;
      
      success = YES ;
    }

  
  
  }
  
  
  
  
  
  // -------------------------------> Telling the vdW contact distances ... <---------------------------------- //
  
  double * vdwdistances = calloc( natomSolvent * natomsoluteSelect , sizeof( double ) );
  
  dzeros( natomsoluteSelect , natomSolvent , vdwdistances ) ;
  
  int vdwSuccess = YES ;
  
  double rvdw1 , rvdw2 , rvdwContact ;
  
  char name1[ 20 ] , name2[ 20 ] , type1[ 20 ] , type2[ 20 ] ;
  
  
  /*
    double tellvdwradius( TOP * p , int howmanytypes , char * whattype ) ;
  
    int tellatomtype( ITP * p , int howmanytypes , char * whatatom , char * at ) ; 
  */
  
  if( vdwCheck == NO )
  {
    printf("\nAs mentioned before , van der Waals radius will not be checked ... \n") ;
    
    vdwSuccess = YES ;
  }
  else if( natomsoluteSelect == natom )
  {
    vdwSuccess = YES ;
  }
  else
  {
    for( iatom = 0 ; iatom < natomsoluteSelect ; iatom ++ )
    {
      //printf("\nSEARCHING for # %d atom in solute &2705& \n" , iatom + 1 ) ; 
      if( ( info = tellatomtype( psoluteAtomDataBase , natomSoluteDataBase , EqAtomList[ iatom ].atomname , type1 ) ) != 0 )
      {
        tellsymbol( EqAtomList[ iatom ].atomname , name1 ) ; // name1 and name2 are actually symbols ... 
        
        //rvdw1 = tellvdwradius( psoluteAtomTypes , natomSoluteTypes , type1 ) ;
        
        rvdw1 = tellvdwradius( name1 ) ;
        
        //printf("\n# %d atom in solute , atomtype is %s , vdw radius is % 10.6f ...\n" , iatom + 1 , type1 , rvdw1 ) ;
      }
      else
      {
        exit( 365 ) ;
      }
      
      for( itmp = natomsoluteSelect ; itmp < natom ; itmp ++ )
      {
        dtmp = *( distances + iatom * natomSolvent + itmp - natomsoluteSelect ) ;
        
        if( ( info = tellatomtype( psolventAtomDataBase , natomSolventDataBase , EqAtomList[ itmp ].atomname , type2 ) ) != 0 )
        {
          tellsymbol( EqAtomList[ itmp ].atomname , name2 ) ; // name1 and name2 are actually symbols ... 
          
          //rvdw2 = tellvdwradius( psoluteAtomTypes , natomSoluteTypes , type2 ) ; // Attention! Here since all the itp atomtypes are in solute itp file ... so we have to use psoluteAtomTypes and natomSoluteTypes ;
          
          rvdw2 = tellvdwradius( name2 ) ;
          
          //printf("\n# %d atom in solvent , atomtype is %s , vdw radius is % 10.6f ...\n" , itmp + 1 , type2 , rvdw2 ) ;
        }
        else
        {
          exit( 365 ) ;
        }
        
        rvdwContact = ( rvdw1 + rvdw2 ) / 2 * vdwFactor ;
        
        *( vdwdistances + iatom * natomSolvent + itmp - natomsoluteSelect ) = rvdwContact ;
        
        if( rvdwContact > dtmp )
        {
          vdwSuccess = NO ;
          
          tellsymbol( EqAtomList[ iatom ].atomname , soluteHotAtomName ) ;
    
          tellsymbol( EqAtomList[ itmp ].atomname , solventHotAtomName ) ;    
    
          printf("\nAfter rotation , distance between # %d solute atom [ %s / %s ] and # %d solvent atom [ %s / %s ] is % 10.6f Angstrom, which is SMALLER than corresponding scaled van der Waals contact radius % 10.6f Angstrom with vdw-scaling-factor % 10.6f ... \n\n" , iatom + 1 , soluteHotAtomName , type1 , itmp + 1 , solventHotAtomName , type2 , dtmp * 10 , rvdwContact * 10 , vdwFactor ) ;
          
          
          //break ;
          
        }
      
      }
    
    
      /*
      if( vdwSuccess == NO )
      {
        break ;
      }
      */
    
    }  
  
  }
  
  

  /* 
  debug = fopen( "vdwContact.deb" , "wb+" ) ;
  
  doutput( debug , natomsoluteSelect , natomSolvent , vdwdistances ) ;
  
  fclose( debug ) ;
  */

  
  // -------------------------------> Output-ing in .gro format with option of printing velocities ... <---------------------------------- //


  //if( success == YES )
  //{
  FILE * poutGROFile = fopen( outGROFileName , "wb+" );

  fprintf( poutGROFile , "Rotated Molecule Geom. Generated by G_NEQALIGN_D : \n"  );

  fprintf( poutGROFile , " %d  \n" , natom );
  
  
  if( velDump == YES )
  {
    for( iatom = 0 ; iatom < natom ; iatom ++ )
    {
      fprintf( poutGROFile , "%5d%-5s%5s%5d%18.13f%18.13f%18.13f%18.14f%18.14f%18.14f\n" , EqAtomList[ iatom ].resnumber , EqAtomList[ iatom ].resname , EqAtomList[ iatom ].atomname , EqAtomList[ iatom ].atomnumber , EqAtomList[ iatom ].cx , EqAtomList[ iatom ].cy , EqAtomList[ iatom ].cz , EqAtomList[ iatom ].vx , EqAtomList[ iatom ].vy , EqAtomList[ iatom ].vz );
    }
  
  }
  else
  {
    for( iatom = 0 ; iatom < natom ; iatom ++ )
    {
      fprintf( poutGROFile , "%5d%-5s%5s%5d%18.13f%18.13f%18.13f\n" , EqAtomList[ iatom ].resnumber , EqAtomList[ iatom ].resname , EqAtomList[ iatom ].atomname , EqAtomList[ iatom ].atomnumber , EqAtomList[ iatom ].cx , EqAtomList[ iatom ].cy , EqAtomList[ iatom ].cz );
    }
  
  }
  
  fprintf( poutGROFile , "%10.6f   %10.6f   %10.6f  \n" , *( boxvector + 0 ) , *( boxvector + 1 ) , *( boxvector + 2 ) );

  /* Note : Because here we are creating this gro file for generating dat file, so the actual boxvector does not matter.
  
     In fact, Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0 and the actual format in .gro file is following : 
     
     " box vectors (free format, space separated reals), values: v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y) "
     
     So , here we output box vector values just to finish this gro format and the values don't matter .                 */


  //}







  if( vdwCheck == NO && distThreshold < 0.00 )
  {
    printf("\n=====> Code : GREY --- Neither van der Waals radius nor universal radius check was performed  ... <=====\n");
    
    exit( 115 ) ;
  } 
  else if( success == YES && vdwSuccess == YES )
  {
    return( 0 ) ;
  }
  else if( success != YES && vdwSuccess == YES )
  {
    printf("\n=====> Code : BLUE --- Fitted configuration does not meet the minimum distance requirement ... <=====\n");
    
    exit( 117 ) ;
  }
  else if( success == YES && vdwSuccess != YES )
  {
    printf("\n=====> Code : RED --- Fitted configuration does not meet the van der Waals distance requirement ... <=====\n");
    
    exit( 119 ) ;
  }
  else if( success != YES && vdwSuccess != YES )
  {
    printf("\n=====> Code : PURPLE --- Fitted configuration does not meet ANY distance requirement ... <=====\n");
    
    exit( 121 ) ;
  }
  else
  {
    printf("\n=====> UNKNOWN ERROR <=====\n") ;
    
    exit( 2811 ) ;
  }



}




//------------------------------------------------------------//



void tellsymbol( char * atomMDname , char * atomSymbol)
{
  // -------> Initiating by buffering atom name from MD to a local char-array ...
  
  int namelength = strlen( atomMDname );
  
  char * namebuffer = calloc( ( namelength + 1 ) , sizeof( char ) ) ;
  
  strcpy( namebuffer , atomMDname );
  
  //int atomlabel ;
  
  
  
  // -------------> Deciding the corresponding atom symbol for current atom  ... 
  
  char firstLetter = * namebuffer ;
  
  char secondLetter = *( namebuffer + 1 );
  
  if( secondLetter <= '9' && secondLetter >= '0' )
  {
    *( atomSymbol + 0 ) = firstLetter ;
    
    *( atomSymbol + 1 ) = '\0';

  } 
  else
  {
    *( atomSymbol + 0 ) = firstLetter ;
    
    *( atomSymbol + 1 ) = secondLetter ;  
    
    *( atomSymbol + 2 ) = '\0' ;
    
   }
  
  
  
  // -------------> Finishing up ...
  
  //return( atomlabel );
  
  


}




//------------------------------------------------------------//


double tellvdwradius( char * whatsymbol )
{
  double vdwradius ;
  
  if( strcmp( whatsymbol , "H" ) == 0 )
  {
    vdwradius = 0.12 ;
  }
  else if( strcmp( whatsymbol , "C" ) == 0 )
  {
    vdwradius = 0.17 ;
  }
  else if( strcmp( whatsymbol , "N" ) == 0 )
  {
    vdwradius = 0.155 ;
  }
  else if( strcmp( whatsymbol , "O" ) == 0 )
  {
    vdwradius = 0.152 ;
  }
  else if( strcmp( whatsymbol , "F" ) == 0 )
  {
    vdwradius = 0.147 ;
  }
  else if( strcmp( whatsymbol , "P" ) == 0 )
  {
    vdwradius = 0.18 ;
  }
  else if( strcmp( whatsymbol , "S" ) == 0 )
  {
    vdwradius = 0.18 ;
  }
  else if( strcmp( whatsymbol , "Cl" ) == 0 )
  {
    vdwradius = 0.175 ;
  }
  else
  {
    vdwradius = 0.839723487 ;
    
    printf("\nUn-recognized atom symbol [ %s ] ...\n" , whatsymbol ) ;
  }
  

  return( vdwradius ) ;


}


//------------------------------------------------------------//


int tellatomtype( ITP * p , int howmanytypes , char * whatatom , char * at )
{
  int itype ;
  
  int flag = 0 ;
  
  for( itype = 0 ; itype < howmanytypes ; itype ++ )
  { 
    //printf("\nThis one , atom name in itp file is [ %s ] ... \n" , p -> atomname ) ;
    
    if( strcmp( whatatom , p -> atomname ) == 0 )
    {
      strcpy( at , p -> type ) ;
      
      flag = 1 ;
      
      //printf("\nFound it ... in itp file ...\n") ;
      
      break ;
    }
    
    p ++ ;
  
  }

  if( flag == 0 )
  { 
    strcpy( at , "!" );
    
    printf("\nNo matching atom name \" %s \"from .itp file ...\n" , whatatom ) ;
  }

  return( flag ) ;


}


//------------------------------------------------------------//









