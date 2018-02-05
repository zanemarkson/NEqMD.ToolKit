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
  
  char soluteITPName[ 500 ] , indexFileName[ 500 ] , eqGROFileName[ 500 ] ;
  
  char outGROFileName[ 500 ] ;
  
  //char outCRDFileName[ 100 ] , outVELFileName[ 100 ] ;
  
  char pGroupName[ 150 ] , qGroupName[ 150 ] ;
  
  int referenceAtomNumber , velDump ;
  
  int natom , natomsoluteGRO , natomsoluteSelect , natomITP ;
  
  int ncart ;
  
  int iatom , icart ;
  
  double * comEqP , * comEqQ , * coordinateEqN ; //, * coordinateNEqN , * comNEqP , * comNEqQ ;
  
  double * comVelocityEqP ;
  
  double * vectorEqPQ , * vectorEqPN ; //, * vectorNEqPQ , * vectorNEqPN ;
  
  
  
  double done = 1.0000 ; double dzero = 0.0000 ; 
  
  int ithree = 3 ; int ione = 1 ; int izero = 0 ;



  double dtmp ; 
  
  double dtmpArray[ 150 ] ;
  
  int itmp ;
  
  char ctmp ;
  
  char tmpString[ 500 ] ;

  
  int icmd ;

  // ----------------------------------> Recording Command-Line Arguments ... <---------------------------------- //
  
  time_t current_time;

  time( &current_time );

  char now[ 300 ] ;

  strcpy( now , ctime( &current_time ) );

  int lennow = strlen( now ) ;

  *( now + lennow - 1 ) = ' ';

    
    
    printf("\n**********************************************************************\n");
      printf("* G_ROTATEMOL_D : Align mol. to reference orientation to avoid clash. *\n");
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
  
  referenceAtomNumber = 1 ;
  
  velDump = NO ;
  
  strcpy( pGroupName , "Donor" ) ;
  
  strcpy( qGroupName , "Acceptor" ) ;
  
  strcpy( soluteITPName , "solute.itp" );
  
  strcpy( indexFileName , "system.index" );
  
  //strcpy( inpCRDFileName , "solute.crd" );
  
  //strcpy( inpVELFileName , "solute.vel" );
  
  
  
  
  
  // -------------------------------> Parsing command-line arguments ... <---------------------------------- //
  
  int exs = 18 ; int exo = 22 ;
  
  
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
	                 
	      case 'i' : strcpy( soluteITPName , *( ++ pcmd ) ) ; 
			 
			         printf("\nCommand-line argument indicates : Input solute itp File name : %s ...\n" , soluteITPName ); 
	      
	                 icmd = icmd + 2 ; 
	                 
	                 break ;
  
          case 'n' : strcpy( indexFileName , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Input index File name : %s ...\n" , indexFileName ); 
	         
	                 icmd = icmd + 2 ;
	         
	                 break ;
          /*
	      case 'r' : strcpy( inpCRDFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Input crd file name : %s ...\n" , inpCRDFileName ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ; 
	      */
	      case 'v' : strcpy( tmpString , *( ++ pcmd ) );
	      
	                 if( strcmp( tmpString , "yes" ) == 0 || strcmp( tmpString , "YES" ) == 0 || strcmp( tmpString , "Y" ) == 0 || strcmp( tmpString , "y" ) == 0 || strcmp( tmpString , "1" ) == 0 )
	                 {  
	                   velDump = YES ;
	                   
	                   printf("\nCommand-line argument indicates : When velocities are available , they WILL be printed out ... \n") ;
	                 }  
	                 else if( strcmp( tmpString , "no" ) == 0 || strcmp( tmpString , "NO" ) == 0 || strcmp( tmpString , "N" ) == 0 || strcmp( tmpString , "n" ) == 0 || strcmp( tmpString , "0" ) == 0 )
	                 {
	                   velDump = NO ; 
	                   
	                   printf("\nCommand-line argument indicates : When velocities are available , they WILL NOT be printed out ... \n") ;
	                 }  
	                 else
	                 {
	                   printf("\nInvalid choice of velocity dumping option : \" %s \" ... You can choose from [ YES / Y / yes / y / 1 ] OR [ NO / N / no / n / 0 ] ...\n\n" , tmpString ) ;
	                   
	                   exit( 249 ) ;
	                 }
	                 
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
	                 
	      case 'N' : referenceAtomNumber = atoi( *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : # %d atom will be used as the L-shape reference ...\n" , referenceAtomNumber );

                     icmd = icmd + 2 ;
                     
                     break;
          
	      case 's' : strcpy( tmpString , *( ++ pcmd ) ) ;
	      
	                 if( strcmp( tmpString , "all" ) == 0 || strcmp( tmpString , "All" ) == 0 || strcmp( tmpString , "ALL" ) == 0 )
	                 {
	                   exs = 20 ;
	                   
	                   printf("\nCommand-line argument indicates : All available atoms will be treated as solute ...\n" );
	                 
	                 }
	                 else
	                 {
	                   natomsoluteSelect = atoi( tmpString );
	      
	                   printf("\nCommand-line argument indicates : The first %d atom will be treated as solute ...\n" , natomsoluteSelect );
	                   
	                   exs = 19 ;
	                   
	                 }  

                     icmd = icmd + 2 ;
                     
                     break;
         

	      case 'h' : printf("\nUsage:  %s [ -i 'input solute itp file name' ] [ -n 'input index file name ( GROMACS .ndx compatible )' ] [ -c 'input EqMD gro file' ][ -p P Group Name ][ -q Q Group Name ][ -N atom number of L-Shape reference atom ] [-s # of atoms to be rotated ] [ -v Whether to print out velocities information ]\n\n" , * argv ); 
	                 
	                 printf("\n===> NOTE : 1) For -v flag , you can choose from [ YES / Y / yes / y / 1 ] OR [ NO / N / no / n / 0 ] ... \n\n");
	                 
	                 printf("\n            2) For -s flag , you can \"All\" or \"all\" or \"ALL\" to indicate all atoms should be rotated ... \n\n");
	                 	      
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
    
    strcat( outGROFileName , ".fit.gro" ) ;
    
    printf("\nBy default , output GRO file name will be %s ...\n\n" , outGROFileName ) ;
  
  }
  
  
  
  
  
  // -------------------------------> Verify File Access ... <---------------------------------- //
  
/*
FILE * psoluteITP , * pindex , * pEqGRO , * pinpcrd , * pinpvel ;
char soluteITPName[ 100 ] , indexFileName[ 100 ] , eqGROFileName[ 100 ] , inpCRDFileName[ 100 ] , inpVELFileName[ 100 ] ;
*/

  if( ( psoluteITP = fopen( soluteITPName , "r" ) ) == NULL )
  {
    printf("\nUser defined .itp file does not exist ... \n");
   
    exit( 3 );
  }
  
  if( ( pindex = fopen( indexFileName , "r" ) ) == NULL )
  {
    printf("\nUser defined index file does not exist ... \n");
   
    exit( 3 );
  }


  /*
  if( ( pinpcrd = fopen( inpCRDFileName , "r" ) ) == NULL )
  {
    printf("\nUser defined crd file does not exist ... \n");
   
    exit( 3 );
  }
  
  if( ( pinpvel = fopen( inpVELFileName , "r" ) ) == NULL )
  {
    printf("\nUser defined vel file does not exist ... \n");
   
    exit( 3 );
  }
  */
  
  
  

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
  else if( exs == 20 )
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



  // -------------------------------> Reading ITP File ... <---------------------------------- //
  
  
  // =====> Pre-Loading ... 
  
  printf("\nNow let's pre-load the .itp file ans see how many atoms it is describing ... \n");
  
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

  natomITP = iload ;
  
  printf("\nThere are %d atoms described in itp file ...\n" , natomITP );

  if( natomITP != natomsoluteSelect )
  {
    printf("\nWhile you specified you wanted %d atoms as the solute, in your itp file, there are only %d atoms being described ...\n" , natomsoluteSelect , natomITP );
    
    printf("\nPlease check your itp file and make sure the # of atoms match ...\n");
    
    printf("\nWe will continue anyway but once we find any required info does not exist in your itp file , we will terminate this process immediately ...\n\n") ;
  }


  // =====> Loading ... 
  
  printf("\nNow let's actually load the .itp file  ... \n");
  
  ITP atomdatabase[ natomsoluteSelect ];
  
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
        
        atomdatabase[ iload -1 ].nr = atoi( cache ) ;
        
        strpickword( buffer , 2 , cache );
        
        strcpy( atomdatabase[ iload -1 ].type , cache );
        
        strpickword( buffer , 3 , cache );
        
        atomdatabase[ iload -1 ].resnr = atoi( cache ) ;
        
        strpickword( buffer , 4 , cache );
        
        strcpy( atomdatabase[ iload -1 ].residue , cache );
        
        strpickword( buffer , 5 , cache );
        
        strcpy( atomdatabase[ iload -1 ].atomname , cache );
        
        strpickword( buffer , 6 , cache );
        
        atomdatabase[ iload -1 ].cgnr = atoi( cache ) ;
        
        strpickword( buffer , 7 , cache );
        
        atomdatabase[ iload -1 ].charge = atof( cache ) ;
        
        strpickword( buffer , 8 , cache );
        
        atomdatabase[ iload -1 ].mass = atof( cache ) ;
        
        printf("\nCharge of this atom is %lf ...\n" , atomdatabase[ iload -1 ].charge ) ;
 
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

  comEqP = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comEqP ) ;
  
  comEqQ = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comEqQ ) ;
  
  comVelocityEqP = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comVelocityEqP ) ;
  
  //comNEqP = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comNEqP ) ;
  
  //comNEqQ = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , comNEqQ ) ;
  
  coordinateEqN = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , coordinateEqN ) ;
  
  //coordinateNEqN = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , coordinateNEqN ) ;
  
  vectorEqPQ = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , vectorEqPQ ) ;
  
  //vectorNEqPQ = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , vectorNEqPQ ) ;

  vectorEqPN = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , vectorEqPN ) ;
  
  //vectorNEqPN = ( double * ) calloc( 3 , sizeof( double ) ) ; dzeros( 3 , 1 , vectorNEqPN ) ;




  double molMassP , molMassQ ;
  
  molMassP = 0.000 ; molMassQ = 0.000 ;
  
  for( iatom = 0 ; iatom < natomPGroup ; iatom ++ )
  {
    if( * ( pGroupAtoms + iatom ) > natomsoluteSelect )
    {
      printf("\nWe need the information about atom # %d , but you itp file does not contain that info. Mission Aborted ...\n" , * ( pGroupAtoms + iatom ) ) ;
      
      exit( 683 ) ;
    
    }
  
    molMassP = molMassP + atomdatabase[ * ( pGroupAtoms + iatom ) - 1 ].mass ;
  
    *( comEqP + 0 ) = *( comEqP + 0 ) + atomdatabase[ *( pGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cx ;
    
    *( comEqP + 1 ) = *( comEqP + 1 ) + atomdatabase[ *( pGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cy ;
  
    *( comEqP + 2 ) = *( comEqP + 2 ) + atomdatabase[ *( pGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cz ;
    
    if( exVelocity == YES )
    {
      *( comVelocityEqP + 0 ) = *( comVelocityEqP + 0 ) + atomdatabase[ *( pGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].vx ;
      
      *( comVelocityEqP + 1 ) = *( comVelocityEqP + 1 ) + atomdatabase[ *( pGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].vy ;
    
      *( comVelocityEqP + 2 ) = *( comVelocityEqP + 2 ) + atomdatabase[ *( pGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].vz ;    
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

    molMassQ = molMassQ + atomdatabase[ * ( qGroupAtoms + iatom ) - 1 ].mass ;
  
    *( comEqQ + 0 ) = *( comEqQ + 0 ) + atomdatabase[ *( qGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cx ;
    
    *( comEqQ + 1 ) = *( comEqQ + 1 ) + atomdatabase[ *( qGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cy ;
  
    *( comEqQ + 2 ) = *( comEqQ + 2 ) + atomdatabase[ *( qGroupAtoms + iatom ) - 1 ].mass * EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cz ;
  
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
  

  printf("\nCOM vector P-> N is % 10.8E , % 10.8E , % 10.8E ...\n\n" , *( vectorEqPN + 0 ) , *( vectorEqPN + 1 ) , *( vectorEqPN + 2 ) ) ;
  

  double * crdFitEqSolute = calloc( 3 * natomsoluteSelect , sizeof( double ) ) ;
  
  double * velFitEqSolute = calloc( 3 * natomsoluteSelect , sizeof( double ) ) ;
  
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



  // -------------------------------> Pausing ... debugging outputs ... <---------------------------------- //
  /*

  debug = fopen( "pEqCoordinates.deb" , "wb+") ;
  
  for( iatom = 0 ; iatom < natomPGroup ; iatom ++ )
  {
    fprintf( debug , "%s\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n" , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].atomname , atomdatabase[ *( pGroupAtoms + iatom ) - 1 ].mass , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cx , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cy , EqAtomList[ *( pGroupAtoms + iatom ) - 1 ].cz ) ;

  }
  
  fprintf( debug , "COM\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n\n" , molMassP , *( comEqP + 0 ) , *( comEqP + 1 ) , *( comEqP + 2 ) ) ;
  
  fclose( debug ) ;






  debug = fopen( "qEqCoordinates.deb" , "wb+") ;
  
  for( iatom = 0 ; iatom < natomQGroup ; iatom ++ )
  {
    fprintf( debug , "%s\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n" , EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].atomname , atomdatabase[ *( qGroupAtoms + iatom ) - 1 ].mass , EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cx , EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cy , EqAtomList[ *( qGroupAtoms + iatom ) - 1 ].cz ) ;

  }
  
  fprintf( debug , "COM\t% 10.8E\t% 10.8E\t% 10.8E\t%10.8E\n" , molMassQ , *( comEqQ + 0 ) , *( comEqQ + 1 ) , *( comEqQ + 2 ) ) ;
  
  fclose( debug ) ;
  
  */

  // -------------------------------> Building Rotation Matrices ... <---------------------------------- //

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
  

  printf("\nMs1 Reads : \n% 16.8f\t% 16.8f\t% 16.8f\n\n% 16.8f\t% 16.8f\t% 16.8f\n\n% 16.8f\t% 16.8f\t% 16.8f\n\n" , *( Ms1 + 0 ) , *( Ms1 + 1 ) , *( Ms1 + 2 ) , *( Ms1 + 3 ) , *( Ms1 + 4 ) , *( Ms1 + 5 ) , *( Ms1 + 6 ) , *( Ms1 + 7 ) , *( Ms1 + 8 ) ) ;
  
  printf("\nMs2 Reads : \n% 16.8f\t% 16.8f\t% 16.8f\n\n% 16.8f\t% 16.8f\t% 16.8f\n\n% 16.8f\t% 16.8f\t% 16.8f\n\n" , *( Ms2 + 0 ) , *( Ms2 + 1 ) , *( Ms2 + 2 ) , *( Ms2 + 3 ) , *( Ms2 + 4 ) , *( Ms2 + 5 ) , *( Ms2 + 6 ) , *( Ms2 + 7 ) , *( Ms2 + 8 ) ) ;



  dgemm_( "T" , "T" , &ithree , &ithree , &ithree , &done , Ms2 , &ithree , Ms1 , &ithree , &dzero , tmpM , &ithree ) ;
  
  dtranspose( 3 , tmpM , tmpM ) ;
  
  printf("\nMs2*Ms1 Reads : \n% 16.8f\t% 16.8f\t% 16.8f\n\n% 16.8f\t% 16.8f\t% 16.8f\n\n% 16.8f\t% 16.8f\t% 16.8f\n\n" , *( tmpM + 0 ) , *( tmpM + 1 ) , *( tmpM + 2 ) , *( tmpM + 3 ) , *( tmpM + 4 ) , *( tmpM + 5 ) , *( tmpM + 6 ) , *( tmpM + 7 ) , *( tmpM + 8 ) ) ;

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
  
  
  printf("\nAt 3rd step , parameters are : \nan = % 12.8f , bn = % 12.8f , cn = % 12.8f \nrabn = % 10.8f \n" , aNFit , bNFit , cNFit , rNabFit ) ;
  
  
  *( Ms3 + 0 ) = aNFit / rNabFit          ;   *( Ms3 + 1 ) = bNFit / rNabFit ;   *( Ms3 + 2 ) = 0.000 ;
  
  *( Ms3 + 3 ) = -1.000 * bNFit / rNabFit ;   *( Ms3 + 4 ) = aNFit / rNabFit ;   *( Ms3 + 5 ) = 0.000 ;
  
  *( Ms3 + 6 ) = 0.000                    ;   *( Ms3 + 7 ) = 0.000           ;   *( Ms3 + 8 ) = 1.000 ;
  
  
  printf("\nMs3 Reads : \n% 16.8f\t% 16.8f\t% 16.8f\n\n% 16.8f\t% 16.8f\t% 16.8f\n\n% 16.8f\t% 16.8f\t% 16.8f\n\n" , *( Ms3 + 0 ) , *( Ms3 + 1 ) , *( Ms3 + 2 ) , *( Ms3 + 3 ) , *( Ms3 + 4 ) , *( Ms3 + 5 ) , *( Ms3 + 6 ) , *( Ms3 + 7 ) , *( Ms3 + 8 ) ) ;
  
  
  dgemm_( "T" , "T" , &ithree , &ithree , &ithree , &done , Ms3 , &ithree , tmpM , &ithree , &dzero , M1 , &ithree ) ;
  
  dtranspose( 3 , M1 , M1 ) ;
  
  
  for( iatom = 0 ; iatom < natomsoluteSelect ; iatom ++ )
  {
    *( dtmpArray + 0 ) = ( *( M1 + 0 ) ) * ( *( crdFitEqSolute + 3 * iatom + 0 ) ) + ( *( M1 + 1 ) ) * ( *( crdFitEqSolute + 3 * iatom + 1 ) ) + ( *( M1 + 2 ) ) * ( *( crdFitEqSolute + 3 * iatom + 2 ) ) ; 

    *( dtmpArray + 1 ) = ( *( M1 + 3 ) ) * ( *( crdFitEqSolute + 3 * iatom + 0 ) ) + ( *( M1 + 4 ) ) * ( *( crdFitEqSolute + 3 * iatom + 1 ) ) + ( *( M1 + 5 ) ) * ( *( crdFitEqSolute + 3 * iatom + 2 ) ) ; 

    *( dtmpArray + 2 ) = ( *( M1 + 6 ) ) * ( *( crdFitEqSolute + 3 * iatom + 0 ) ) + ( *( M1 + 7 ) ) * ( *( crdFitEqSolute + 3 * iatom + 1 ) ) + ( *( M1 + 8 ) ) * ( *( crdFitEqSolute + 3 * iatom + 2 ) ) ; 

    EqAtomList[ iatom ].cx = *( dtmpArray + 0 ) ;
    
    EqAtomList[ iatom ].cy = *( dtmpArray + 1 ) ;
  
    EqAtomList[ iatom ].cz = *( dtmpArray + 2 ) ;
    
    
    if( exVelocity == YES )
    {
      *( dtmpArray + 0 ) = ( *( M1 + 0 ) ) * ( *( velFitEqSolute + 3 * iatom + 0 ) ) + ( *( M1 + 1 ) ) * ( *( velFitEqSolute + 3 * iatom + 2 ) ) + ( *( M1 + 2 ) ) * ( *( velFitEqSolute + 3 * iatom + 2 ) ) ; 
  
      *( dtmpArray + 1 ) = ( *( M1 + 3 ) ) * ( *( velFitEqSolute + 3 * iatom + 0 ) ) + ( *( M1 + 4 ) ) * ( *( velFitEqSolute + 3 * iatom + 5 ) ) + ( *( M1 + 2 ) ) * ( *( velFitEqSolute + 3 * iatom + 2 ) ) ; 
  
      *( dtmpArray + 2 ) = ( *( M1 + 6 ) ) * ( *( velFitEqSolute + 3 * iatom + 0 ) ) + ( *( M1 + 7 ) ) * ( *( velFitEqSolute + 3 * iatom + 8 ) ) + ( *( M1 + 2 ) ) * ( *( velFitEqSolute + 3 * iatom + 2 ) ) ; 
  
      EqAtomList[ iatom ].vx = *( dtmpArray + 0 ) ;
      
      EqAtomList[ iatom ].vy = *( dtmpArray + 1 ) ;
    
      EqAtomList[ iatom ].vz = *( dtmpArray + 2 ) ;
    
    
    }
  
  
  
  
  }



  // -------------------------------> Output-ing in .gro format with option of printing velocities ... <---------------------------------- //


  FILE * poutGROFile = fopen( outGROFileName , "wb+" );

  fprintf( poutGROFile , "Rotated Molecule Geom. Generated by G_ROTATEMOL_D : \n"  );

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









  return( 0 ) ;




}



