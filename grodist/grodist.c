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
  double atommass ;

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



int main( int argc, char * argv[] )
{
  // -----------------------------------------> Variables ... <--------------------------------------------- //
  
  char ** pcmd ; pcmd = argv ;
  
  FILE * psoluteITP , * pindex , * pEqGRO , * pinpcrd , * pinpvel ;
  
  FILE * debug ;
  
  char indexFileName[ 100 ] , eqGROFileName[ 100 ] ;
  
  char outGROFileName[ 100 ] ;
  
  //char outCRDFileName[ 100 ] , outVELFileName[ 100 ] ;
  
  char pGroupName[ 50 ] , qGroupName[ 50 ] ;
  
  int natom , natomsoluteGRO ;
  
  int ncart ;
  
  int iatom , icart ;
  
  double * comEqP , * comEqQ ;
  
  double * comVelocityEqP , * comVelocityEqQ ;
  
  double * vectorEqPQ , * vectorVelocityEqPQ ; //, * vectorNEqPQ , * vectorNEqPN ;
  
  
  
  double done = 1.0000 ; double dzero = 0.0000 ; 
  
  int ithree = 3 ; int ione = 1 ; int izero = 0 ;



  double dtmp ; 
  
  double dtmpArray[ 50 ] ;
  
  int itmp ;
  
  char ctmp ;
  
  char tmpString[ 100 ] ;

  
  int icmd ;
  
  // -----------------------------------------> Functions ... <--------------------------------------------- //
  
  double tellmass( char * atomMDname ) ;
  
  void tellsymbol( char * atomMDname , char * atomSymbol ) ;

  // ----------------------------------> Recording Command-Line Arguments ... <---------------------------------- //
  
  time_t current_time;

  time( &current_time );

  char now[ 300 ] ;

  strcpy( now , ctime( &current_time ) );

  int lennow = strlen( now ) ;

  *( now + lennow - 1 ) = ' ';

    
    
    printf("\n***************************************************************************\n");
      printf("* G_GRODIST_D : Calculating COM Distances between Groups from .gro file.  *\n");
      printf("*                                                                         *\n");
      printf("*  ");
  for( icmd = 0 ; icmd < argc ; icmd ++ )
  {
    printf("%s " , *( pcmd + icmd ) );
  }
  printf("\n");
      printf("*                                                                        *\n");
      printf("*                                                                        *\n");
      printf("* Current Time : %s                               *\n" , now );
      printf("*                                                                        *\n");
      printf("*                                                                        *\n");
      printf("**************************************************************************\n");

 
  
  // ------------------------------------------> Defaults ... <------------------------------------------ //
  
  strcpy( pGroupName , "Donor" ) ;
  
  strcpy( qGroupName , "Acceptor" ) ;
  
  strcpy( indexFileName , "system.index" );
  
  //strcpy( inpCRDFileName , "solute.crd" );
  
  //strcpy( inpVELFileName , "solute.vel" );
  
  
  
  
  
  // -------------------------------> Parsing command-line arguments ... <---------------------------------- //
  
  int exo = 22 ;
  
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
	      case 'c' : strcpy( eqGROFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Input .gro File name : %s ...\n" , eqGROFileName ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 //exgro = 19 ;
	                 
	                 break ; 


	      case 'o' : strcpy( outGROFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Output .gro File name : %s ...\n" , outGROFileName ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 exo = 23 ;
	                 
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
	                 

	      case 'h' : printf("\nUsage:  %s [ -n 'input index file name (in GROMACS .ndx format)' ] [ -c 'input EqMD gro file' ] [ -o output file name ] [ -p P Group Name ][ -q Q Group Name ]\n\n" , * argv ); 
	                 
	                 //printf("\n===> NOTE : 1) For -v flag , you can choose from [ YES / Y / yes / y / 1 ] OR [ NO / N / no / n / 0 ] ... \n\n");
	                 	      
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
  
  
  
  int lenEqGROFileName = strlen( eqGROFileName ) ;
  
  if( exo == 22 )
  {
    strncpy( outGROFileName , eqGROFileName , lenEqGROFileName - 4 ) ;
    
    *( outGROFileName + lenEqGROFileName - 4 ) = '\0' ;
    
    strcat( outGROFileName , ".dist" ) ;
    
    printf("\nBy default , output GRO file name will be %s ...\n\n" , outGROFileName ) ;
  
  }
  
  
  
  
  
  // -------------------------------> Verify File Access ... <---------------------------------- //
  
/*
FILE * psoluteITP , * pindex , * pEqGRO , * pinpcrd , * pinpvel ;
char soluteITPName[ 100 ] , indexFileName[ 100 ] , eqGROFileName[ 100 ] , inpCRDFileName[ 100 ] , inpVELFileName[ 100 ] ;
*/

  if( ( pindex = fopen( indexFileName , "r" ) ) == NULL )
  {
    printf("\nUser defined index file does not exist ... \n");
   
    exit( 3 );
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
  
  int exVelocity = 0 ;
  

  
  
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
        
        EqAtomList[ iload ].atommass = tellmass( EqAtomList[ iload ].atomname );

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


  /*
  debug = fopen( "allMass.deb" , "wb+" ) ;
  
  for( iatom = 0 ; iatom < natom ; iatom ++ )
  {
    fprintf( debug , "%5d%5s\t%5s\t% 10.6f\n\n" , EqAtomList[ iatom ].resnumber , EqAtomList[ iatom ].resname , EqAtomList[ iatom ].atomname , EqAtomList[ iatom ].atommass ) ;
  
  }
  
  fclose( debug ) ;
  */

  
  // -------------------------------> Read index file and acquire info of P and Q ... <---------------------------------- //
  
  int next_atom_index , current_atom_index ;
  
  int natomPGroup , natomQGroup ;
  
  int ndxinfo ;
  
  // =====> P Group , Pre-Loading 
  
  rewind( pindex ) ;
  
  fsearch( pindex , pGroupName ) ;
  
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


  // =====> Q Group , Pre-Loading 
  
  rewind( pindex ) ;
  
  fsearch( pindex , qGroupName ) ;
  
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
        
      if( ( ndxinfo = fscanf( pindex , "%s" , tmpString ) ) == EOF )
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
        
      if( ( ndxinfo = fscanf( pindex , "%s" , tmpString ) ) == EOF )
      {
         break ;  
      }

    }

    

  }




  // -------------------------------> Calculate COM for P and Q in Eq... <---------------------------------- //

  //double comEqP[ 3 ] , comEqQ[ 3 ] , coordinateN[ 3 ] , comNEqP[ 3 ] , comNEqQ[ 3 ] ;
  
  //  double * comEqP , * comEqQ , * coordinateN , * comNEqP , * comNEqQ ;
  
  //  double * vectorEqPQ , * vectorNEqPQ;
  
  double distanceBetweenPandQ = 0.000 ;

  comEqP = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comEqP ) ;
  
  comEqQ = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comEqQ ) ;
  
  comVelocityEqP = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comVelocityEqP ) ;
  
  comVelocityEqQ = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comVelocityEqQ ) ;
  
  //comNEqP = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comNEqP ) ;
  
  //comNEqQ = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comNEqQ ) ;
  
  //coordinateEqN = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , coordinateEqN ) ;
  
  //coordinateNEqN = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , coordinateNEqN ) ;
  
  vectorEqPQ = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , vectorEqPQ ) ;
  
  //vectorNEqPQ = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , vectorNEqPQ ) ;

  vectorVelocityEqPQ = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , vectorVelocityEqPQ ) ;
  
  //vectorNEqPN = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , vectorNEqPN ) ;




  double molMassP , molMassQ ;
  
  // ==> P
  
  molMassP = 0.000 ; molMassQ = 0.000 ;
  
  for( iatom = 0 ; iatom < natomPGroup ; iatom ++ )
  {
    molMassP = molMassP + EqAtomList[ * ( pGroupAtoms + iatom ) - 1 ].atommass ;
  
    *( comEqP + 0 ) = *( comEqP + 0 ) + EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].atommass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cx ;
    
    *( comEqP + 1 ) = *( comEqP + 1 ) + EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].atommass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cy ;
  
    *( comEqP + 2 ) = *( comEqP + 2 ) + EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].atommass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cz ;
    
    if( exVelocity == YES )
    {
      *( comVelocityEqP + 0 ) = *( comVelocityEqP + 0 ) + EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].atommass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].vx ;
      
      *( comVelocityEqP + 1 ) = *( comVelocityEqP + 1 ) + EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].atommass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].vy ;
    
      *( comVelocityEqP + 2 ) = *( comVelocityEqP + 2 ) + EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].atommass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].vz ;    
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

  
  printf("\nMass of P fragment is % 10.8E ...\n" , molMassP ) ;
  
  printf("\nCOM coordinates of P in EqMD is % 10.8E , % 10.8E , % 10.8E ...\n\n" , *( comEqP + 0 ) , *( comEqP + 1 ) , *( comEqP + 2 ) ) ;

  if( exVelocity == YES )
  {
    printf("\nCOM velocities of P in EqMD is % 10.8E , % 10.8E , % 10.8E ...\n\n" , *( comVelocityEqP + 0 ) , *( comVelocityEqP + 1 ) , *( comVelocityEqP + 2 ) ) ;
  }

  // ==> Q

  for( iatom = 0 ; iatom < natomQGroup ; iatom ++ )
  {
    molMassQ = molMassQ + EqAtomList[ * ( qGroupAtoms + iatom ) - 1 ].atommass ;
  
    *( comEqQ + 0 ) = *( comEqQ + 0 ) + EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].atommass * EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cx ;
    
    *( comEqQ + 1 ) = *( comEqQ + 1 ) + EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].atommass * EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cy ;
  
    *( comEqQ + 2 ) = *( comEqQ + 2 ) + EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].atommass * EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cz ;
  
    if( exVelocity == YES )
    {
      *( comVelocityEqQ + 0 ) = *( comVelocityEqQ + 0 ) + EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].atommass * EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].vx ;
      
      *( comVelocityEqQ + 1 ) = *( comVelocityEqQ + 1 ) + EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].atommass * EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].vy ;
    
      *( comVelocityEqQ + 2 ) = *( comVelocityEqQ + 2 ) + EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].atommass * EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].vz ;    
    }
  
  }
  
  *( comEqQ + 0 ) = *( comEqQ + 0 ) / molMassQ ;
  
  *( comEqQ + 1 ) = *( comEqQ + 1 ) / molMassQ ;
  
  *( comEqQ + 2 ) = *( comEqQ + 2 ) / molMassQ ;

  if( exVelocity == YES )
  {
    *( comVelocityEqQ + 0 ) = *( comVelocityEqQ + 0 ) / molMassQ ;
    
    *( comVelocityEqQ + 1 ) = *( comVelocityEqQ + 1 ) / molMassQ ;
    
    *( comVelocityEqQ + 2 ) = *( comVelocityEqQ + 2 ) / molMassQ ;
  }


  printf("\nMass of Q fragment is % 10.8E ...\n" , molMassQ ) ;
  
  printf("\nCOM coordinates of Q in EqMD is % 10.8E , % 10.8E , % 10.8E ...\n\n" , *( comEqQ + 0 ) , *( comEqQ + 1 ) , *( comEqQ + 2 ) ) ;

  if( exVelocity == YES )
  {
    printf("\nCOM velocities of Q in EqMD is % 10.8E , % 10.8E , % 10.8E ...\n\n" , *( comVelocityEqQ + 0 ) , *( comVelocityEqQ + 1 ) , *( comVelocityEqQ + 2 ) ) ;
  }

  
  // ==> Between P and Q 
  
  *( vectorEqPQ + 0 ) = *( comEqQ + 0 ) - *( comEqP + 0 ) ;
  
  *( vectorEqPQ + 1 ) = *( comEqQ + 1 ) - *( comEqP + 1 ) ;
  
  *( vectorEqPQ + 2 ) = *( comEqQ + 2 ) - *( comEqP + 2 ) ;
  
  if( exVelocity == YES )
  {
    *( vectorVelocityEqPQ + 0 ) = *( comVelocityEqP + 0 ) - *( comVelocityEqQ + 0 ) ;
    
    *( vectorVelocityEqPQ + 1 ) = *( comVelocityEqP + 1 ) - *( comVelocityEqQ + 1 ) ;
    
    *( vectorVelocityEqPQ + 2 ) = *( comVelocityEqP + 2 ) - *( comVelocityEqQ + 2 ) ;
  
  
  }


  distanceBetweenPandQ = distanceBetweenPandQ + ( *( vectorEqPQ + 0 ) ) * ( *( vectorEqPQ + 0 ) ) ;
  
  distanceBetweenPandQ = distanceBetweenPandQ + ( *( vectorEqPQ + 1 ) ) * ( *( vectorEqPQ + 1 ) ) ;
  
  distanceBetweenPandQ = distanceBetweenPandQ + ( *( vectorEqPQ + 2 ) ) * ( *( vectorEqPQ + 2 ) ) ;
  
  distanceBetweenPandQ = sqrt( distanceBetweenPandQ) ; 

  /* 
  debug = fopen( "donor.deb" , "wb+" ) ;
  
  for( iatom = 0 ; iatom < natomPGroup ; iatom ++ )
  {
    fprintf( debug , "%d\t%5s\t% 10.6f     % 12.8E\t% 12.8E\t% 12.8E\n\n" , *( pGroupAtoms + iatom ) , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].atomname , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].atommass , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cx , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cy , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cz ) ;
  
  }
  
  fclose( debug ) ;


  debug = fopen( "acceptor.deb" , "wb+" ) ;
  
  for( iatom = 0 ; iatom < natomQGroup ; iatom ++ )
  {
    fprintf( debug , "%d\t%5s\t% 10.6f     % 12.8E\t% 12.8E\t% 12.8E\n\n" , *( qGroupAtoms + iatom ) , EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].atomname , EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].atommass , EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cx , EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cy , EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cz ) ;
  
  }
  
  fclose( debug ) ;
  */


  // -------------------------------> Output-ing in .gro format with option of printing velocities ... <---------------------------------- //


  FILE * poutGROFile = fopen( outGROFileName , "wb+" );

  if( exVelocity == YES )
  {
    fprintf( poutGROFile , "% 12.8E\t% 12.8E\t% 12.8E\t% 12.8E\n\n" , distanceBetweenPandQ , *( vectorVelocityEqPQ + 0 ) , *( vectorVelocityEqPQ + 1 ) , *( vectorVelocityEqPQ + 2 ) ) ;
  }
  else
  {
    fprintf( poutGROFile , "% 12.8E\n\n\n" , distanceBetweenPandQ ) ;
  }

  fclose( poutGROFile ) ;

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  return( 0 ) ;






}
























double tellmass( char * atomMDname )
{
  // -------> Initiating by buffering atom name from MD to a local char-array ...
  
  int namelength = strlen( atomMDname );
  
  char * namebuffer = malloc( ( namelength + 1 ) * sizeof( char ) ) ;
  
  strcpy( namebuffer , atomMDname );
  
  double atomMass;
  
  
  
  // -------------> Deciding the corresponding atom symbol for current atom  ... 
  
  char firstLetter = * namebuffer ;
  
  char secondLetter = *( namebuffer + 1 );
  
  if( secondLetter <= '9' && secondLetter >= '0' )
  {
    switch ( firstLetter )
    {
      case 'H' : atomMass = 1.000 ; break ;
      
      case 'h' : atomMass = 1.000 ; break ;
      
      case 'B' : atomMass = 10.81 ; break ; 
      
      case 'b' : atomMass = 10.81 ; break ; 
      
      case 'C' : atomMass = 12.00 ; break ; 
      
      case 'c' : atomMass = 12.00 ; break ; 
      
      case 'N' : atomMass = 14.00 ; break ; 
      
      case 'n' : atomMass = 14.00 ; break ; 
      
      case 'O' : atomMass = 16.00 ; break ; 
      
      case 'o' : atomMass = 16.00 ; break ;
      
      case 'F' : atomMass = 19.00 ; break ;
      
      case 'f' : atomMass = 19.00 ; break ;
      
      case 'P' : atomMass = 30.97 ; break ;
      
      case 'p' : atomMass = 30.97 ; break ;
    
      case 'S' : atomMass = 32.06 ; break ;
      
      case 's' : atomMass = 32.06 ; break ;
      
      case 'K' : atomMass = 39.10 ; break ;
      
      case 'k' : atomMass = 39.10 ; break ;
      
      default  : printf("\nUnkown Atom detected : %c... Mission aborting ...\n\n" , firstLetter );
       
                 exit(1);

    }
  } 
  else
  {
      switch ( firstLetter )
      {
        case 'S' : atomMass = 28.09 ; break ; // ---> Si <--- //
        
        case 's' : atomMass = 28.09 ; break ; // ---> Si <--- //
        
        case 'C' : atomMass = 35.50 ; break ; // ---> Cl <--- //
        
        case 'c' : atomMass = 35.50 ; break ; // ---> Cl <--- //
        
        case 'M' : atomMass = 24.31 ; break ; // ---> Mg <--- //
        
        case 'm' : atomMass = 24.31 ; break ; // ---> Mg <--- //
        
        case 'Z' : atomMass = 65.38 ; break ; // ---> Zn <--- //
        
        case 'z' : atomMass = 65.38 ; break ; // ---> Zn <--- //
        
        default  : printf("\nUnkown Atom detected : %c%c... Mission aborting ...\n\n" , firstLetter , secondLetter );
        
                   exit(1);

      }
   
   }
  
  
  
  // -------------> Finishing up ...
  
  return( atomMass );
  

}




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



