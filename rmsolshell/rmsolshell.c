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
  double atommass ;
  

} GRO ;


typedef struct datform
{
  double cx ;
  double cy ;
  double cz ;
  //char atomsymbol[3];
  int atomlabel ;
  double atommass ;
  //char which_layer ;
  double atomcharge ;

} DAT ;

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
  
  // ==============> Declaration : Variables ... <========= //
  

  char ** pcmd ;

  pcmd = argv ; //pcmd ++ ; 

  int icmd = 1 ;
  
  FILE * pgroinput , * poutGROFile , * pinpHessian ;
  
  char inpgroname[ 100 ] , outGROFileName[ 100 ] ; 
  
  int natom , ncart , natomselect ; 
  
  int iatom , icart ;
  
  double distThreshold ;
  
  int vdwCheck = NO , unidistCheck = NO ; 
  
  double vdwFactor = 1.000 ;
  
  
  int debuggingMode = NO ;


  int iline , iload ;
  
  int irow , icol ;
  
  int blank_signal , groinfo ;

  char buffer[ MAXCHARINLINE ] ;
  
  char cache[ MAXCHARINLINE ] ;

  



  double dtmp ; 
  
  double dtmpArray[ 100 ] ;
  
  int itmp , info ; 
  
  char tmpString[ 150 ] , tmp_char ;


  FILE * debug ;




  // ==============> Declaration : In-Place Utility Functions ... <========= //


  void tellsymbol( char * atomMDname , char * atomSymbol) ;
  
  int tellatomtype( ITP * p , int howmanytypes , char * whatatom , char * at ) ;
  
  double tellvdwradius( char * whatsymbol ) ;
  
  double tellmass( char * atomMDname ) ;
  

  // ----------------------------------> Recording Command-Line Arguments ... <---------------------------------- //
  
  time_t current_time;

  time( &current_time );

  char now[ 300 ] ;

  strcpy( now , ctime( &current_time ) );

  int lennow = strlen( now ) ;

  *( now + lennow - 1 ) = ' ';

    
    
    printf("\n**********************************************************************\n");
      printf("* G_RMSOLSHELL_D : Remove solvent shell for faster solvation.        *\n");
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

 
  
  // ==============> Handling the file names and Charge & Multiplicity ... <========= //
  
  
  // -------> Parsing the Command Line Arguments ... 
  
  pcmd = argv ; icmd = 1 ;
  
  //int exn = 10 ; int exr  = 16 ; int exR = 18 ; int exH  = 22 ; int exM = 28 ; int exL = 30 ; 

  int exf = 1 ; int exo = 22 ; 
  
  int exs = 18 ; int ext = 88 ;
  
  int exH = 70 ; int exvdw = 60 ;

  char * flag ;
  
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
	      
	      case 'f' : strcpy( inpgroname , *( ++ pcmd ) ) ; 
			 
                     printf("\nCommand-line argument indicates : Input File name : %s ...\n" , inpgroname ); 
	      
	                 exf = 7 ; 
	                 
	                 icmd = icmd + 2 ; 
	                 
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
	      
	      case 't' : strcpy( tmpString , *( ++ pcmd ) ) ;
	      
	                 if( strcmp( tmpString , "YES" ) == 0 || strcmp( tmpString , "yes" ) == 0 || strcmp( tmpString , "Yes" ) == 0 )
	                 {
	                   printf("\nCommand-line argument indicates : Universal distance check WILL be performed with default threshold = 1.20 Angstrom\n" );
	                   
	                   distThreshold = 1.200 ;
	                   
	                   unidistCheck = YES ;
	                   	                   
	                 }
                     else if( strcmp( tmpString , "NO" ) == 0 || strcmp( tmpString , "no" ) == 0 || strcmp( tmpString , "No" ) == 0 )
                     {
	                   printf("\nCommand-line argument indicates : Universal distance check will NOT be performed" );
	                   
	                   distThreshold = -0.4000 ;
	                   
	                   unidistCheck = NO ;

                     }
			         else if( *( tmpString ) <= '9' && *( tmpString ) >= '0' )
			         {
			           distThreshold = atof( tmpString ) ;
			 
	                   printf("\nCommand-line argument indicates : Universal distance check WILL be performed with user-defined threshold = % 10.6f\n" , distThreshold );
			 
			           unidistCheck = YES ;
			         
			         }
			         else
			         {
			           printf("\nInvalid choice of universal distance checking ... Aborting ...\n") ;

			           exit( 285 ) ;
			 
			         }

                     ext = 89 ;
                     
                     distThreshold = distThreshold / 10.00 ; // Change to nm unit since all the calculations here are in nm unit ... I know, it IS weird ...
                     
                     icmd = icmd + 2 ;
                     
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

		           exit( 286 ) ;
			 
			         }

                     exvdw = 61 ;
                     
                     icmd = icmd + 2 ;
                     
                     break;
                     
                     

          case 'g' : strcpy( tmpString , *( ++ pcmd ) );

	                 if( strcmp( tmpString , "YES" ) == 0 || strcmp( tmpString , "yes" ) == 0 || strcmp( tmpString , "Yes" ) == 0 )
	                 {
	                   printf("\nCommand-line argument indicates : Debugging mode ON ... Extra info will be printed out ... \n" );
	                   
	                   debuggingMode = YES ;
	                   
	                 }
                     else if( strcmp( tmpString , "NO" ) == 0 || strcmp( tmpString , "no" ) == 0 || strcmp( tmpString , "No" ) == 0 )
                     {
	                   printf("\nCommand-line argument indicates : Production mode ... minimum info will be printed ...\n" );
	                   
	                   debuggingMode = NO ;

                     }
                     else
                     {
			           printf("\nInvalid choice of debugging mode switch option ... Aborting ...\n") ;
			           
			           exit( 286 ) ;

                     }

                     icmd = icmd + 2 ;
                     
                     break;



	      case 'o' : strcpy( outGROFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Output GRO file name : %s ...\n" , outGROFileName ); 
	                 
	                 exo = 23 ;
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ; 

		                 
	      case 'h' : printf("\nUsage:  % 15s [ -f 'input gro file name' ] [ -s # of atoms chosen as the solute ]" , * argv) ;
	                 printf("\n                        [ -w whether to perform van der Waals contacting check ( YES / Yes / yes ) or ( NO / No / no ) or floating point number to indicate scaling factor ]");
	                 printf("\n                        [ -t threshold for universal distance check ] [ -o Output GRO File Name ]\n");
	                 //printf("\n           [ -c 'input EqMD gro file of solvent' ][ -p P Group Name ][ -q Q Group Name ][ -N atom number of L-Shape reference atom ]");
	                 //printf("\n           [ -w whether to perform van der Waals contacting check ( YES / Yes / yes ) or ( NO / No / no ) or floating point number to indicate scaling factor ]");
	                 //printf("\n           [ -s # of atoms in solute molecule ] [ -t distance threshold to accept one alignment , unit = Angstrom ]\n\n" ); 
	                 
	                 printf("\nNote : 1) For \"-w\" option, YES/Yes/yes will cause the vdw-check to perform with default scaling 1.00 while NO/No/no will shut down the vdw-check.\n");
	                 printf("\n          Specifying a floating number will also cause the vdw-check to perform but the floating number will be the user-defined scaling factor (vdwFactor).\n");

	                 printf("\n       2) For \"-t\" option, YES/Yes/yes will cause the universal-distance-check to perform with default threshold 1.20Å while NO/No/no will shut down the unidist-check.\n");
	                 printf("\n          Specifying a floating number will also cause the unidist-check to perform but the floating number will be the user-defined distance threshold (vdwFactor).\n");

	                 printf("\n       3) For \"-s\" option, [ -s all ] or [ -s All ] indicates all atoms chosen as solute;\n");
	                 printf("\n          Default for -s is all atoms when nresidue = 1  or natom in 1st residue when nresidue != 1 \n\n\n");

	                 printf("\n       4) For \"-g\" option, YES/Yes/yes will turn on the debugging mode which allows additional information to be printed out; Default is NO ;\n\n\n");

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

    
  
  // --------------> Summarizing and File Access ...
  
  printf("\n\nInput GMX .gro file name is : %s ...\n" , inpgroname );
  
  // --------> Setting DEFAULT value for important variables ... 


  if( exf == 1 )
  {
    strcpy( inpgroname , "system.gro" );
    
    printf("\nNo input .gro file provided, default \"system.gro\" in play ...\n") ;
  
  }

  
  int lenEqGROFileName = strlen( inpgroname ) ;
  
  if( exo == 22 )
  {
    strncpy( outGROFileName , inpgroname , lenEqGROFileName - 4 ) ;
    
    *( outGROFileName + lenEqGROFileName - 4 ) = '\0' ;
    
    strcat( outGROFileName , ".fit.gro" ) ;
    
    printf("\nBy default , output GRO file name will be %s ...\n\n" , outGROFileName ) ;
  
  }
  

  if( exvdw == 60 )
  {
    vdwCheck = YES ;
    
    vdwFactor = 1.000 ;
  
  }

  if( ext == 88 )
  {
    distThreshold = -0.400 ; // Unit = nm ;
    
    unidistCheck = NO ;  
  
  }
  
  if( unidistCheck + vdwCheck == NO )
  {
    printf("\nFatal Error : There must be one type of criterion to define the shell ... \n") ;
    
    exit( 40 ) ;
  
  }
  else if( unidistCheck + vdwCheck == 2 * YES )
  {
    printf("\nWARNING : universal distance check and van der Waal contact distance check cannot be invoked at the same time ...\n") ;
    
    printf("\nBy default, only van der Waals contact check will be kept ...\n\n") ;
    
    distThreshold = -0.400 ;
    
    unidistCheck = NO ;
  
  }
  
  
  
  
  

  
  // ==============> Reading information from GMX .gro file ... <========= //
  
  char grotitlestring[MAXLINE];
  
  iline = 3 ;
  
  iload = 0 ;
  
  int natomgroline , natomgrotitle ;
  
  int exVelocity = NO ;
  

  
  
  if( ( pgroinput = fopen( inpgroname , "r" ) ) == NULL )
  {
    printf("\nUser defined .gro file %s does not exist ... \n" , inpgroname );
   
    exit( 3 );
  }
  else
  {
    rewind( pgroinput );
    
    //printf("\nCurrent character is %c ... \n" , fgetc( pEqGRO ) );
    
    fskip( pgroinput , 1 );

    fscanf( pgroinput , "%d" , &natomgrotitle );
    
    fskip( pgroinput , 1 );
    
    printf("\n Second line of .gro file says it is describing %d atoms ... \n\n" , natomgrotitle );
    
    
    while( ( groinfo = freadline( buffer , MAXCHARINLINE , pgroinput , ';' ) ) != 0 )
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

  /*
  rewind( pgroinput ) ;
    
  fskip( pgroinput , 2 ) ;
  
  while( ( groinfo = freadline( buffer , MAXCHARINLINE , pgroinput , ';' ) ) != 0 )
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
  */
  
  
  natomgroline = preLoadGRO( pgroinput ) ;
  
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
 




  // =====> Actually Reading the GRO File ...
  
  GRO EqAtomList[ natom ] ; 
  
  rewind( pgroinput );
    
  fskip( pgroinput , 1 );

  fskip( pgroinput , 1 );
    
  printf("\nNow let's read the actual .gro file  ... \n");
  
  iload = 0 ; iline = 0 ;
  
  while( ( groinfo = freadline( buffer , MAXCHARINLINE , pgroinput , ';' ) ) != 0 )
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
    fskip( pgroinput , natomgroline - natomgrotitle ) ;
  }
  
  double boxvector[ 3 ]; 
  
  fscanf( pgroinput , "%lf" , boxvector + 0 );

  fscanf( pgroinput , "%lf" , boxvector + 1 );

  fscanf( pgroinput , "%lf" , boxvector + 2 );
  
  
  


  // Debugging output ...
  
  /*
  debug = fopen("allVelocities.deb" , "wb+" ) ;
  
  for( iatom = 0 ; iatom < natom ; iatom ++ )
  {
    fprintf( debug , "\n% 10.6f\t% 10.6f\t% 10.6f\n" , EqAtomList[ iatom ].vx , EqAtomList[ iatom ].vy , EqAtomList[ iatom ].vz );
  
  }
  
  fclose( debug ) ;
  */
  




  // ---------------> Figuring out how many solvent molecules (residues) are in the system ...
  
  int nresidue = EqAtomList[ natom - 1 ].resnumber ;
  
  // ---------------> Figuring out how many atoms are in each  residue ...
  
  int * n_of_atom_in_residues = calloc( nresidue , sizeof( int ) ) ;
  
  izeros( nresidue , 1 , n_of_atom_in_residues );
  
  //double * molecularmass = calloc( nresidue , sizeof( double ) ) ;
  
  //dzeros( nresidue , 1 , molecularmass );
  
  int iresidue = 1 ;
  
  for( iatom = 0 ; iatom < natom ; iatom ++ )
  {
     if( EqAtomList[ iatom ].resnumber == ( iresidue + 1 ) ) 
     {
       //printf("\n#%d residue has %d atoms ...\n" , iresidue , *( n_of_atom_in_residues + iresidue - 1 ) );
       
       iresidue ++ ;
       
       //printf("\n\nStarting ... # %d residue (molecule) ... \n\n" , iresidue );
       
     }
     
     else if( EqAtomList[ iatom ].resnumber == iresidue )
     {
       //printf("\nStill in this residue ... \n");
       
       //continue ;
     }
     
     else 
     {
       printf("\nSomething is wrong with the EqAtomList  ... Mission Aborting ...\n\n");
       
       exit(1);
     }
     
     
     *( n_of_atom_in_residues + iresidue - 1 ) = *( n_of_atom_in_residues + iresidue - 1 ) + 1 ;
     
     ///*( molecularmass + iresidue - 1 ) = *( molecularmass + iresidue - 1 ) + atomcast[ iatom ].atommass ;
  
  } 
    
  
  
  
  //---> Dealing with the "Default senario " of Chosen Atoms ... 
  
  int nsoluteResidue = 1 ;
  
  int natomUncoveredBySolute = 0 ;
  
  itmp = 0 ; 
 
  if( exs == 18 && nresidue == 1 )
  {
    natomselect = natom ;
  }
  else if( exs == 18 && nresidue > 1 )
  {
    natomselect = *( n_of_atom_in_residues + 0 ) ;
    
    nsoluteResidue = 1 ;
    
    natomUncoveredBySolute = 0 ;
    
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
  
    for( iresidue = 0 ; iresidue < nresidue ; iresidue ++ )
    {
      itmp = itmp + ( *( n_of_atom_in_residues + iresidue ) ) ;
      
      if( natomselect <= itmp )
      {
        nsoluteResidue = iresidue + 1 ;
        
        natomUncoveredBySolute = itmp - natomselect ;
        
        break ;
      
      }
    
    }
    
    printf("\n[ %d ] residue(s) is/are involved in solute ... with [ %d ] left-over atoms...\n\n" , nsoluteResidue , natomUncoveredBySolute ) ;
    
    if( natomUncoveredBySolute != 0 && nsoluteResidue == 1 )
    {
      printf("\n=> WARNING : ONLY PART OF THE FIRST RESIDUE IS SELECTED AS SOLUTE ... <=\n") ;
    }
    else if( natomUncoveredBySolute != 0 && nsoluteResidue > 1 )
    {
      printf("\n=> WARNING : SELECTED SOLUTE CONTAINS PARTIAL RESIDUES ... <=\n");
    
    }
    else if( natomUncoveredBySolute == 0 )
    {
      printf("\nSo ... it looks that 1 or multiple WHOLE residue(s) is/are considered as solute ... \n\n") ;
    }
    else
    {
      printf("\nSomething is wrong with the selection of solute ... @826\n\n");
      
      exit( 78 ) ;
    
    }
  
  

  }
  else if( exs == 20 )
  {
    natomselect = natom ;
    
    nsoluteResidue = nresidue ;
    
    natomUncoveredBySolute = 0 ;
  }
  else
  {
    printf("\nSomething is wrong with the atom selection process ... NAtom = %d , NAtomSelect = %d ... \n" , natom , natomselect );
    
    exit( 78 );
  }
  



  
  // -------------------------------> Telling the Distance between solute and solvent atoms ... <---------------------------------- //

  int natomSolvent = natom - natomselect ;

  double * distances = calloc( natomSolvent * natomselect , sizeof( double ) );
  
  int pos , posSolute , posSolvent ;
  
  char solventHotAtomName[ 3 ] , soluteHotAtomName[ 3 ] ;
  
  int success = NO ; 
  
  double distThresholdScalingFactor = 1.00 ;
  
  double ratomDistance ; 
  
  
  
  
  for( iatom = 0 ; iatom < natomselect ; iatom ++ ) // Calculating all solute-solvent distances ... 
  {
    dtmp = 0.000 ;
    
    for( itmp = natomselect ; itmp < natom ; itmp ++ )
    {
      dtmp = 0.000 ;
      
      dtmp = dtmp + ( EqAtomList[ iatom ].cx - EqAtomList[ itmp ].cx ) * ( EqAtomList[ iatom ].cx - EqAtomList[ itmp ].cx ) ;
    
      dtmp = dtmp + ( EqAtomList[ iatom ].cy - EqAtomList[ itmp ].cy ) * ( EqAtomList[ iatom ].cy - EqAtomList[ itmp ].cy ) ;
      
      dtmp = dtmp + ( EqAtomList[ iatom ].cz - EqAtomList[ itmp ].cz ) * ( EqAtomList[ iatom ].cz - EqAtomList[ itmp ].cz ) ;
      
      //printf("\n# %d atom , also # %d atom in solvent , coordinate is [ % 10.6f\t% 10.6f\t% 10.6f\t] \n\n" , itmp , itmp - natomselect , EqAtomList[ itmp ].cx , EqAtomList[ itmp ].cy , EqAtomList[ itmp ].cz ) ;
      
      *( distances + iatom * natomSolvent + itmp - natomselect ) = sqrt( dtmp ) ;
    
    }
    

  }
  
  
  /*
  debug = fopen( "distances.deb" , "wb+" ) ;
  
  doutput( debug , natomselect , natomSolvent , distances ) ;
  
  fclose( debug ) ;
  
  */
  
  // -------------------------------> Telling the vdW contact distances ... <---------------------------------- //
  
  
  double * vdwdistances = calloc( natomSolvent * natomselect , sizeof( double ) );
  
  dzeros( natomselect , natomSolvent , vdwdistances ) ;
  
  double * vdwradii_solute = calloc( natomselect , sizeof( double ) ) ;
  
  dzeros( natomselect , 1 , vdwradii_solute ) ;
  
  double * vdwradii_solvent = calloc( natomSolvent , sizeof( double ) ) ;
  
  dzeros( natomSolvent , 1 , vdwradii_solvent ) ;
  
  
  double rvdw1 , rvdw2 , rvdwContact ;
  
  char name1[ 20 ] , name2[ 20 ] , type1[ 20 ] , type2[ 20 ] ;

  
  
  if( vdwCheck == YES && unidistCheck == NO )
  {
    
    for( iatom = 0 ; iatom < natomselect ; iatom ++ )
    {
      tellsymbol( EqAtomList[ iatom ].atomname , name1 ) ; // name1 and name2 are actually symbols ... 
      
      rvdw1 = tellvdwradius( name1 ) ;
    
      *( vdwradii_solute + iatom ) = rvdw1 ;
    }
    
    for( itmp = natomselect ; itmp < natom ; itmp ++ )
    {
        tellsymbol( EqAtomList[ itmp ].atomname , name2 ) ; // name1 and name2 are actually symbols ... 
        
        rvdw2 = tellvdwradius( name2 ) ;
        
        *( vdwradii_solvent + itmp - natomselect ) = rvdw2 ;
    
    }
    
    
    
    
    
    
    for( iatom = 0 ; iatom < natomselect ; iatom ++ )
    {
      //printf("\nSEARCHING for # %d atom in solute &2705& \n" , iatom + 1 ) ;
      
      //tellsymbol( EqAtomList[ iatom ].atomname , name1 ) ; // name1 and name2 are actually symbols ... 
      
      //rvdw1 = tellvdwradius( name1 ) ;
      
      rvdw1 = *( vdwradii_solute + iatom ) ;
      
      /*
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
      */
      
      for( itmp = natomselect ; itmp < natom ; itmp ++ )
      {
        //dtmp = *( distances + iatom * natomSolvent + itmp - natomselect ) ;
        
        //tellsymbol( EqAtomList[ itmp ].atomname , name2 ) ; // name1 and name2 are actually symbols ... 
        
        //rvdw2 = tellvdwradius( name2 ) ;
        
        rvdw2 = *( vdwradii_solvent + itmp - natomselect ) ;
        
        /*
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
        */
        
        rvdwContact = ( rvdw1 + rvdw2 ) / 2 * vdwFactor ;
        
        *( vdwdistances + iatom * natomSolvent + itmp - natomselect ) = rvdwContact ;
        
      
      }

    
    }  
  
  }
  else if( vdwCheck == NO && unidistCheck == YES ) 
  {
    for( iatom = 0 ; iatom < natomselect ; iatom ++ )
    {
      for( itmp = natomselect ; itmp < natom ; itmp ++ )
      {
        *( vdwdistances + iatom * natomSolvent + itmp - natomselect ) = distThreshold ;
      }
  
    }
  
  }
  



  if( debuggingMode == YES )
  {
    debug = fopen( "vdwradii_solute.deb" , "wb+") ;
  
    doutput( debug , natomselect , 1 , vdwradii_solute ) ;
  
    fclose( debug ) ;
  
  
    debug = fopen( "vdwradii_solvent.deb" , "wb+") ;
  
    doutput( debug , natomSolvent , 1 , vdwradii_solvent ) ;
  
    fclose( debug ) ;
  
  
    debug = fopen( "vdwdistances.deb" , "wb+" ) ;
  
    doutput( debug , natomselect , natomSolvent , vdwdistances ) ;
  
    fclose( debug ) ;
  
  }


  
  // -------------------------------> Defining Solvent Shell ... <---------------------------------- //
  
  int * stayInShell_residue = calloc( nresidue , sizeof( int ) ) ;
  
  for( iresidue = 0 ; iresidue < nresidue ; iresidue ++ ) *( stayInShell_residue + iresidue ) = NO ;
  
  int * stayInShell_atom = calloc( natom , sizeof( int ) ) ;
  
  for( iatom = 0 ; iatom < natom ; iatom ++ )  *( stayInShell_atom + iatom ) = NO ;
  
  /*
  int * stayInShell_cart = calloc( ncart , sizeof( int ) ) ;
  
  for( icart = 0 ; icart < ncart ; icart ++ )  *( stayInShell_cart + icart ) = NO ;
  */
  
  int iatom_acc = 0 ; int icart_acc = 0 ; int iresidue_acc = 0 ;
  
  int currentResidue_atomstart , currentResidue_atomend ;
  
  int currentResidue_cartstart , currentResidue_cartend ;
  
  int nresidueInShell , natomInShell , ncartInShell ;
  
    
  
  //---> Firstly , do it for each atom ... All the "uncovered atoms" will be automatically included into the solute ...   
  
  for( iatom = 0 ; iatom < natomselect + natomUncoveredBySolute ; iatom ++ )
  {
    *( stayInShell_atom + iatom ) = NO ;
  }
  
  /*
  for( icart = 0 ; icart < 3 * ( natomselect + natomUncoveredBySolute ) ; icart ++ )
  {
    *( stayInShell_cart + icart ) = NO ;
  }
  */
  
  
  for( itmp = natomselect + natomUncoveredBySolute ; itmp < natom ; itmp ++ )
  {
    for( iatom = 0 ; iatom < natomselect ; iatom ++ )
    {
      rvdwContact = *( vdwdistances + iatom * natomSolvent + itmp - natomselect ) ;
      
      ratomDistance = *( distances + iatom * natomSolvent + itmp - natomselect ) ;
      
      if( rvdwContact >= ratomDistance )
      {
        *( stayInShell_atom + itmp ) = YES ;
        
        break ;
      }

    }
  
  }
  
  
    
  
  
  //---> Secondly , reflect these results for each residue ... 
  
  for( iresidue = 0 ; iresidue < nsoluteResidue ; iresidue ++ )
  {
    *( stayInShell_residue + iresidue ) = NO ; 
    
    // This is for the solute itself. Of course solute will not be in the to-be-removed
    // shell.
  }
  
  iatom_acc = natomselect + natomUncoveredBySolute ; // C-Label ... 
  
  for( iresidue = nsoluteResidue ; iresidue < nresidue ; iresidue ++ )
  {
    currentResidue_atomstart = iatom_acc ; // C-Label again ... 
    
    iatom_acc = iatom_acc + ( *( n_of_atom_in_residues + iresidue ) ) ;
    
    currentResidue_atomend = iatom_acc ; // C-Label ... 
    
    for( iatom = currentResidue_atomstart ; iatom < currentResidue_atomend ; iatom ++ )
    {
      if( *( stayInShell_atom + iatom ) == YES )
      {
        *( stayInShell_residue + iresidue ) = YES ;
        
        break ;
       
      }
    
    }
  
  
  
  }
  
  
  
  
  //---> Third step : clear out the stayInShell_atom & re-assign ( solute atoms not touched )
  
  for( iatom = natomselect + natomUncoveredBySolute ; iatom < natom ; iatom ++ ) 
  {
    *( stayInShell_atom + iatom ) = NO ;
  }
  
  
  
  iatom_acc = natomselect + natomUncoveredBySolute ; // C-Label ... 
  
  icart_acc = 3 * iatom_acc ; // C-Label ... 
  
  for( iresidue = nsoluteResidue ; iresidue < nresidue ; iresidue ++ )
  {
    currentResidue_atomstart = iatom_acc ; // C-Label again ... 
    
    currentResidue_cartstart = 3 * currentResidue_atomstart ; // C-Label again ...
    
    iatom_acc = iatom_acc + ( *( n_of_atom_in_residues + iresidue ) ) ;
    
    currentResidue_atomend = iatom_acc ; // C-Label ... 
    
    currentResidue_cartend = 3 * currentResidue_atomend ; // C-Label again ...
    
    if( *( stayInShell_residue + iresidue ) == YES )
    {
      for( iatom = currentResidue_atomstart ; iatom < currentResidue_atomend ; iatom ++ )
      {
        *( stayInShell_atom + iatom ) = YES ;
      }
      
      /*
      for( icart = currentResidue_cartstart ; icart < currentResidue_cartend ; icart ++ )
      {
        *( stayInShell_cart + icart ) = YES ;
      }
      */

    }
  
  
  
  }
  
  
  
  //---> Summarizing ... 
  
  nresidueInShell = isum( stayInShell_residue , nresidue , 1 ) ;
  
  natomInShell = isum( stayInShell_atom , natom , 1 ) ;
  
  //ncartInShell = isum( stayInShell_cart , ncart , 1 ) ;
  
  printf("\nIn total,  there are [ %d ] residues including solute in-shell , [ %d ] atoms in-shell  ...\n\n" , nresidueInShell , natomInShell ) ;
  
  printf("\n<SHELL>\n\n  %d  \n\n</SHELL>\n\n" , nresidueInShell ) ; 
  
  //---> debugging output
  
  char component ;
  
  if( debuggingMode == YES )
  {
    debug = fopen("stayInShell.residue" , "wb+") ;
  
    iatom_acc = 0 ;
  
    for( iresidue = 0 ; iresidue < nresidue ; iresidue ++ )
    {
      fprintf( debug , "\n[ %5d ] [ %s ] [ %d ]\n" , iresidue + 1 , EqAtomList[ iatom_acc ].resname , *( stayInShell_residue + iresidue ) ) ;
    
      iatom_acc = iatom_acc + ( *( n_of_atom_in_residues + iresidue ) ) ;
  
    }
  
    fclose( debug ) ;
  
  
    debug = fopen("stayInShell.atom" , "wb+") ;
  
    for( iatom = 0 ; iatom < natom ; iatom ++ )
    {
      fprintf( debug , "\n[ %5d ] [ %s ] [ %d ]" , iatom + 1 , EqAtomList[ iatom ].atomname , *( stayInShell_atom + iatom ) ) ;
  
    }
  
    fclose( debug ) ;
  

    /*(
    debug = fopen("stayInShell.cart" , "wb+") ;
  
    for( icart = 0 ; icart < ncart ; icart ++ )
    {
      itmp = icart % 3 ; 
    
      iatom = ( icart - itmp ) / 3 ;
    
      if( itmp == 0 )
      {
        component = 'X' ;
      }
      else if( itmp == 1 )
      {
        component = 'Y' ;
      }
      else
      {
        component = 'Z' ;
      }
    
      fprintf( debug , "\n[ %5d ] [ %s ].[ %c ] [ %d ]" , icart + 1 , EqAtomList[ iatom ].atomname , component , *( stayInShell_cart + icart ) ) ;
  
    }
  
    fclose( debug ) ;
    */
  
  }
  
  
  

  
  // -------------------------------> Output the GRO File ... <---------------------------------- //

  poutGROFile = fopen( outGROFileName , "wb+") ;
  
  iatom_acc = 1 ;
  
  iresidue_acc = 1 ;
  
  int iatom_acc_fromBegining = 0 ;
  
  int iatom_real ;


  fprintf( poutGROFile , "Shell-Removed Geom. Generated by G_RMSOLSHELL_D : \n"  );

  fprintf( poutGROFile , " %d  \n" , natom - natomInShell );
  


  if( exVelocity == NO )
  {
    for( iresidue = 0 ; iresidue < nresidue ; iresidue ++ )
    {
      printf("\n===> Checking %d residue ... <===\n" , iresidue ) ;
      
      if( *( stayInShell_residue + iresidue ) == NO )
      {
        printf("\n=> And it is not in shell ... SO KEEP IT <=\n") ;
      
        for( iatom = 0 ; iatom < *( n_of_atom_in_residues + iresidue ) ; iatom ++ )
        {
          iatom_real = iatom_acc_fromBegining + iatom ; // C-Label ...
        
          //printf("\nOutputing No.[ %d ] atom in whole system ...\n" , iatom_real ) ;
        
          fprintf( poutGROFile , "%5d%-5s%5s%5d%18.13f%18.13f%18.13f\n" , iresidue_acc , EqAtomList[ iatom_real ].resname , EqAtomList[ iatom_real ].atomname , iatom_acc , EqAtomList[ iatom_real ].cx , EqAtomList[ iatom_real ].cy , EqAtomList[ iatom_real ].cz );
      
          iatom_acc ++ ;
        
        }
            
        iresidue_acc ++ ;
      
      }
  
      iatom_acc_fromBegining = iatom_acc_fromBegining + ( *( n_of_atom_in_residues + iresidue ) ) ;
  
    }
  
  }
  else if( exVelocity == YES )
  {
    for( iresidue = 0 ; iresidue < nresidue ; iresidue ++ )
    {
      printf("\n===> Checking %d residue ... <===\n" , iresidue ) ;
      
      if( *( stayInShell_residue + iresidue ) == NO )
      {
        printf("\n=> And it is not in shell ... SO KEEP IT <=\n") ;
      
        for( iatom = 0 ; iatom < *( n_of_atom_in_residues + iresidue ) ; iatom ++ )
        {
          iatom_real = iatom_acc_fromBegining + iatom ; // C-Label ...
        
          //printf("\nOutputing No.[ %d ] atom in whole system ...\n" , iatom_real ) ;
        
          fprintf( poutGROFile , "%5d%-5s%5s%5d%18.13f%18.13f%18.13f%18.14f%18.14f%18.14f\n" , iresidue_acc , EqAtomList[ iatom_real ].resname , EqAtomList[ iatom_real ].atomname , iatom_acc , EqAtomList[ iatom_real ].cx , EqAtomList[ iatom_real ].cy , EqAtomList[ iatom_real ].cz , EqAtomList[ iatom_real ].vx , EqAtomList[ iatom_real ].vy , EqAtomList[ iatom_real ].vz );
      
          iatom_acc ++ ;
        
        }
            
        iresidue_acc ++ ;
      
      }
  
      iatom_acc_fromBegining = iatom_acc_fromBegining + ( *( n_of_atom_in_residues + iresidue ) ) ;
  
    }
  
  }
  
  fprintf( poutGROFile , "%10.6f   %10.6f   %10.6f  \n" , *( boxvector + 0 ) , *( boxvector + 1 ) , *( boxvector + 2 ) );
  


  // -----------------> The End ... Really???
  
  return(0);
  
  
  
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
    
      case 'S' : atomMass = 32.00 ; break ;
      
      case 's' : atomMass = 32.00 ; break ;
      
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



  
  
  










