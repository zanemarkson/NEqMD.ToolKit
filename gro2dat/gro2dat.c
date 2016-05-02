#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
  
  FILE * pgroinput , * pitpinput , * pCNDOoutput ;

  char inpgroname[ MAXCHARINLINE ] , inpitpname[ MAXCHARINLINE ] , outCNDOname[ MAXCHARINLINE ] ;

  int multiplicity , charge ;

  int natom , ncart , natomselect , natomDump ; 
  
  int iatom ;


  double dtmp ; 
  
  double dtmpArray[ MAXCHARINLINE ] ;
  
  int itmp ; 
  
  char tmpString[ MAXCHARINLINE ] ;

  // --------> Declaring utility functions ...

  double tellmass( char * atomMDname ) ;
  
  int tellatom( char * atomMDname ) ;

  void dzeros( int dimrow, int dimcol, double *p ) ;
  
  void izeros( int dimrow, int dimcol, int *p ) ;

  // --------> Some default values  ... 

  multiplicity = 1 ;
  
  charge = 0 ;
  

  
  // ==============> Handling the file names and Charge & Multiplicity ... <========= //
  
  
  // -------> Parsing the Command Line Arguments ... 

  char ** pcmd ;

  pcmd = argv ; //pcmd ++ ; 

  int icmd = 1 ;
  
  
  //int exn = 10 ; int exr  = 16 ; int exR = 18 ; int exH  = 22 ; int exM = 28 ; int exL = 30 ; 

  int exf = 1 ; int exo = 3 ; int exx = 36 ; 
  
  int exs = 18 ; int exD = 28 ;

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

	      case 'o' : strcpy( outCNDOname , *( ++ pcmd ) ) ;
	      
	                 printf("\nCommand-line argument indicates : Output File name : %s ...\n" , outCNDOname ); 
	                 
	                 exo = 9 ; 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ; 
	                 
	      case 'x' : strcpy( inpitpname , *( ++ pcmd ) ) ;
	      
	                 if( strcmp( inpitpname , "dummy") == 0 || strcmp( inpitpname , "none") == 0 || strcmp( inpitpname , "None") == 0 )
	                 {
	                   exx = 99 ;
	                   
	                   printf("\nLooks like we are not using point charge here ... \n") ;
	                 }
	                 
	                 else
	                 {
	                   printf("\nCommand-line argument indicates : Input .itp File name : %s ...\n" , inpitpname ); 
	                 
	                   exx = 37 ;
	                   
	                 }  
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;
         
	      case 's' : strcpy( tmpString , *( ++ pcmd ) ) ;
	      
	                 if( strcmp( tmpString , "all" ) == 0 || strcmp( tmpString , "All" ) == 0 || strcmp( tmpString , "ALL" ) == 0 )
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


	      case 'D' : strcpy( tmpString , *( ++ pcmd ) ) ;
	      
	                 if( strcmp( tmpString , "all" ) == 0 || strcmp( tmpString , "All" ) == 0 || strcmp( tmpString , "ALL" ) == 0 )
	                 {
	                   exD = 30 ;
	                   
	                   printf("\nCommand-line argument indicates : All atoms will be chosen to be dumped in output file ...\n" );
	                 }
	                 else if( strcmp( tmpString , "solute" ) == 0 || strcmp( tmpString , "Solute" ) == 0 || strcmp( tmpString , "SOLUTE" ) == 0 )
	                 {
	                   exD = 31 ;
	                   
	                   printf("\nCommand-line argument indicates : Only selected solute atoms will be chosen to be dumped in output file ...\n" );
	                 }
	                 else
	                 {
	                   printf("\nReceived information : %s ...\n" , tmpString ) ;
	                   
	                   natomDump = atoi( tmpString ); 
	                  
	                   exD = 29 ;
	                   
	                   printf("\nCommand-line argument indicates : First %d atoms will be dumped into output file ...\n" , natomDump );
	                 }
	                 
	                  
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ;
          
          
          /*
	      case 'r' : radiusM = atof( *( ++ pcmd ) ); 
	                 
	                 exr = 17 ; 
	      
	                 printf("\nCommand-line argument indicates : User defined radius for middle layer : %12.8f ...\n" , radiusM ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ;

	      case 'R' : radiusL = atof( *( ++ pcmd ) ); 
	      
	                 exR = 19 ;
	      
	                 printf("\nCommand-line argument indicates : User defined radius for lower layer : %12.8f ...\n" , radiusL ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ;
	
          
	      case 'H' : strcpy( highest_method , *( ++ pcmd ) ); 
	      
	                 exH = 23 ;
	      
	                 printf("\nCommand-line argument indicates : User defined method for Highest layer : %s ...\n" , highest_method ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ;
	
	
          
	      case 'M' : strcpy( middle_method , *( ++ pcmd ) ); 
	      
	                 exM = 29 ;
	      
	                 printf("\nCommand-line argument indicates : User defined method for Middle layer : %s ...\n" , middle_method ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ;
	

	      case 'L' : strcpy( lower_method , *( ++ pcmd ) ); 
	      
	                 exL = 31 ;
	      
	                 printf("\nCommand-line argument indicates : User defined method for Lower layer : %s ...\n" , lower_method ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ;
	       */
	
	
	      case 'h' : printf("\nUsage:  %s [ -f 'input gro file name' ] [(optional) -o 'output CNDO input file name' ] [ -x input GMX .itp file ] [ -s # of atoms chosen as the solute ] [ -D # of atoms chosen to be Dumped in output file ]\n\n" , * argv ); 
	                 
	                 printf("\nNOTE : 1) [ -s all ] or [ -s All ] indicates all atoms chosen as solute;\n\n       2) Default for -s is all atoms when nresidue = 1  or natom in 1st residue when nresidue != 1 \n");
	                 
	                 printf("\n       3) User can use [ -x none / None / dummy ] to indicate all selected atoms ( may be less than total number of atoms ) will be treated as solute \n\n");
	                 
	                 //printf("\m       3) ") ;
	                 
	                 //printf("\nUsage:  %s [ -t G09 calculation type : 1=ONIOM ; 2=Point Charge ] [ -f 'input gro file name' ] [(optional) -o 'output g09 file name' ] [ -n # of layers (integer) ] [ (optional) -r radius of middle layer (real) ] [-R radius of lower layer (real) ] [ -H method for Highest layer (string) ] [ (optional) -M method for Middle layer (string) ] [ -L method for Lower layer (string) ] [ -x input GMX .itp file ]\n\n" , * argv ); 
	      
	                 //exh = 9 ;
	                 
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
  
  // --------> Setting DEFAULT value for important variables ... 

    
  
  char defoutname[ MAXCHARINLINE ] , definpname[ MAXCHARINLINE ] , defitpname[ MAXCHARINLINE ];

  int inputnamelength , outputnamelength , inputitpnamelength;

  printf("\nID = %d ...\n" , exo * exf * exx );

  switch( exo * exf * exx )  // File Names ... int exf = 1 / 7 ; int exo = 3 / 9; int exx = 36 / 37 - 99 ; 
  {
     case 108  : strcpy( inpgroname , "sys.gro" ); 
     
                 strcpy( inpitpname , "sys.itp" ); 
                 
                 strcpy( outCNDOname , "sys.dat" ); 
                 
                 break ;
     
     case 1*3*99  : strcpy( inpgroname , "sys.gro" ); 
     
                 //strcpy( inpitpname , "sys.itp" ); 
                 
                 strcpy( outCNDOname , "sys.dat" ); 
                 
                 break ;


     case 756  : inputnamelength = strlen( inpgroname ) ; 

                 strncpy( defoutname, inpgroname , inputnamelength - 4 ) ;
             
                 *( defoutname + inputnamelength - 4 ) = '\0' ;
             
                 strcat( defoutname, ".dat") ;
              
                 strcpy( outCNDOname , defoutname ) ;
                 
                 strncpy( defitpname , inpgroname , inputnamelength - 4 );
                 
                 *( defitpname + inputnamelength - 4 ) = '\0' ;
                 
                 strcat( defitpname , ".itp" );
                 
                 strcpy( inpitpname , defitpname );
              
              
                 break ;
 
     case 3*7*99  : inputnamelength = strlen( inpgroname ) ; 

                 strncpy( defoutname, inpgroname , inputnamelength - 4 ) ;
             
                 *( defoutname + inputnamelength - 4 ) = '\0' ;
             
                 strcat( defoutname, ".dat") ;
              
                 strcpy( outCNDOname , defoutname ) ;
                 
                 break ;
              
    
     case 7*9*99 : printf("\n\nHoorayyyyyyyy ... Both input and output name are specified !!!\n\n");
     
                 printf("\nAnd ... We don't need itp file ! \n");
               
                  
                 break ;
             
     case 324  : outputnamelength = strlen( outCNDOname ) ;   
     
                 strncpy( definpname, outCNDOname , outputnamelength - 4 ) ;
             
                 *( definpname + outputnamelength - 4 ) = '\0' ;
             
                 strcat( definpname, ".gro") ;
             
                 strcpy( inpgroname , definpname );
               
               
                 strncpy( defitpname , outCNDOname , outputnamelength - 4 );
  
                 *( defitpname + outputnamelength - 4 ) = '\0' ;
  
                 strcat( defitpname , ".itp" );
  
                 strcpy( inpitpname , defitpname );
 
                 break ;


     case 1*9*99  : outputnamelength = strlen( outCNDOname ) ;   
     
                 strncpy( definpname, outCNDOname , outputnamelength - 4 ) ;
             
                 *( definpname + outputnamelength - 4 ) = '\0' ;
             
                 strcat( definpname, ".gro") ;
             
                 strcpy( inpgroname , definpname );
               
 
                 break ;


     case 111  : strcpy( inpgroname , "sys.gro" ); 
     
                 strcpy( outCNDOname , "sys.inp" ); 
                 
                 break ;




     case 777  : inputnamelength = strlen( inpgroname ) ; 
      
                 strncpy( defoutname, inpgroname , inputnamelength - 4 ) ;
             
                 *( defoutname + inputnamelength - 4 ) = '\0' ;
             
                 strcat( defoutname, ".dat") ;
              
                 strcpy( outCNDOname , defoutname ) ;
                 
                 break ;

                 
     case 2331 : printf("\n\nHoorayyyyyyyy ... All file names are specified !!!\n\n");
     
                 break ;
                 
     case 333  : outputnamelength = strlen( outCNDOname ) ;   
     
                 strncpy( definpname, outCNDOname , outputnamelength - 4 ) ;
             
                 *( definpname + outputnamelength - 4 ) = '\0' ;
             
                 strcat( definpname, ".gro") ;
             
                 strcpy( inpgroname , definpname );
                 
                 break ;
  

  
  }
  


  
  // --------------> Summarizing and File Access ...
  
  printf("\n\nInput GMX .gro file name is : %s ...\n" , inpgroname );
  
  printf("\nOutput CNDO dat file name is %s ...\n" , outCNDOname );
    
  
  

  if( ( pitpinput = fopen( inpitpname , "r" ) ) == NULL && exx != 99 )
  {
    printf("\nUser defined .itp file %s does not exist ... \n" , inpitpname );
   
    exit( 3 );
  }
  
  
  // ==============> Reading information from GMX .gro file ... <========= //
  
  char grotitlestring[MAXLINE];
  
  int iline = 3 ;
  
  int iload = 0 ;
  
  int blank_signal , groinfo ;

  char buffer[ MAXCHARINLINE ] ;
  
  char cache[ MAXCHARINLINE ] ;

  char tmp_char ;
  
  int natomgroline , natomgrotitle ;
  
  int exVelocity = 0 ;
  

  
  
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
  
  GRO atomlist[ natom ] ; 
  
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
              
        sscanf( buffer , "%5d%5s" , &atomlist[ iload ].resnumber , atomlist[ iload ].resname );

        //printf( "%s\t" , atomlist[ iatom ].resname );

        //sscanf( pEqGRO , "%s" , EqAtomList[ iload ].atomname );
        
        strpickword( buffer , 2 , cache ) ;
        
        strcpy( atomlist[ iload ].atomname , cache ) ;

        //printf( "%s" , atomlist[ iatom ].atomname );

        //sscanf( pEqGRO , "%d" , &EqAtomList[ iload ].atomnumber );
        
        strpickword( buffer , 3 , cache ) ;
        
        atomlist[ iload ].atomnumber = atoi( cache ) ;

        //printf( "\nWorking on No. %d atom ...\n" , EqAtomList[ iatom ].atomnumber );

        //sscanf( pEqGRO , "%lf" , &EqAtomList[ iload ].cx ); //printf("\n Cx is %lf ...\t" , EqAtomList[ iatom ].cx);
        
        strpickword( buffer , 4 , cache ) ; atomlist[ iload ].cx = atof( cache ) ;
        
        //sscanf( pEqGRO , "%lf" , &EqAtomList[ iload ].cy ); //printf("\n Cy is %lf ...\t" , EqAtomList[ iatom ].cy);
        
        strpickword( buffer , 5 , cache ) ; atomlist[ iload ].cy = atof( cache ) ;
        
        //sscanf( pEqGRO , "%lf" , &EqAtomList[ iload ].cz ); //printf("\n Cz is %lf ...\n\n" , EqAtomList[ iatom ].cz);
        
        strpickword( buffer , 6 , cache ) ; atomlist[ iload ].cz = atof( cache ) ;
        
        if( exVelocity == YES )
        {
          strpickword( buffer , 7 , cache ) ; atomlist[ iload ].vx = atof( cache ) ;
          
          strpickword( buffer , 8 , cache ) ; atomlist[ iload ].vy = atof( cache ) ;
          
          strpickword( buffer , 9 , cache ) ; atomlist[ iload ].vz = atof( cache ) ;
        
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
  
  
  

  // ---------------> Figuring out how many solvent molecules (residues) are in the system ...
  
  int nresidue = atomlist[ natom - 1 ].resnumber ;
  
  // ---------------> Figuring out how many atoms are in each  residue ...
  
  int * n_of_atom_in_residues = calloc( nresidue , sizeof( int ) ) ;
  
  izeros( nresidue , 1 , n_of_atom_in_residues );
  
  //double * molecularmass = calloc( nresidue , sizeof( double ) ) ;
  
  //dzeros( nresidue , 1 , molecularmass );
  
  int iresidue = 1 ;
  
  for( iatom = 0 ; iatom < natom ; iatom ++ )
  {
     if( atomlist[ iatom ].resnumber == ( iresidue + 1 ) ) 
     {
       //printf("\n#%d residue has %d atoms ...\n" , iresidue , *( n_of_atom_in_residues + iresidue - 1 ) );
       
       iresidue ++ ;
       
       //printf("\n\nStarting ... # %d residue (molecule) ... \n\n" , iresidue );
       
     }
     
     else if( atomlist[ iatom ].resnumber == iresidue )
     {
       //printf("\nStill in this residue ... \n");
       
       //continue ;
     }
     
     else 
     {
       printf("\nSomething is wrong with the atomlist or atomcast ... Mission Aborting ...\n\n");
       
       exit(1);
     }
     
     
     *( n_of_atom_in_residues + iresidue - 1 ) = *( n_of_atom_in_residues + iresidue - 1 ) + 1 ;
     
     ///*( molecularmass + iresidue - 1 ) = *( molecularmass + iresidue - 1 ) + atomcast[ iatom ].atommass ;
  
  } 
    
  
  
  
  //Just debugging ...
  
  for( iresidue = 1 ; ( iresidue - 1 ) < nresidue ; iresidue ++ )
  {                                                   
  
    printf("\nThere are %d atoms in No. %d residue ... \n" , *( n_of_atom_in_residues + iresidue - 1 ) , iresidue );
           
    //printf("\nThe molecular mass of No. %d residue is %lf ...\n" , iresidue , *( molecularmass + iresidue - 1 ) );
  
  }
 
  
  
  
  //---> Dealing with the "Default senario " of Chosen Atoms ... 
 
  if( exs == 18 && nresidue == 1 )
  {
    natomselect = natom ;
  }
  else if( exs == 18 && nresidue > 1 )
  {
    natomselect = *( n_of_atom_in_residues + 0 ) ;
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
    printf("\nYou have selected %d atoms for rotation ... There are %d atoms en toto in this system ... \n" , natomselect , natom );
  
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
  


  if( exD == 28 )
  {
    natomDump = natom ;
  }
  else if( exD == 30 )
  {
    natomDump = natom ;
  }
  else if( exD == 31 )
  {
    natomDump = natomselect ;
    
    printf("\nBased on command-line input , only [ %d ] solute atoms will be dumped ...\n\n" , natomselect ) ;
  }
  else if( exD == 29 )
  {
    if( natomDump <= natom && natomDump > natomselect )
    {
      printf("\nBased on command-line input , besides solute , [ %d ] atoms from solvent will be dumped ...\n\n" , natomDump - natomselect ) ;
    }
    else if( natomDump == natomselect )
    {
      printf("\nBased on command-line input , only [ %d ] solute atoms will be dumped ...\n\n" , natomselect ) ;
    }
    else if( natomDump < natomselect && natomDump >= 0 )
    {
      printf("\nWARNING : You chose to dump ONLY PART OF SOLUTE atoms [ %d ] into your cndo .dat file ...\n\n" , natomDump ) ;
    }
    else if( natomDump < 0 )
    {
      natomDump = natomselect ;
      
      printf("\nBased on command-line input , only [ %d ] solute atoms will be dumped ...\n\n" , natomselect ) ;
    }
    else if( natomDump > natom )
    {
      printf("\nWARNING : There are in total only [ %d ] atoms in the system , so you cannot request more than [ %d ] atoms dumped ...\n\n" , natom , natom ) ;
      
      printf("\nNow we re-set the natomDump to be natom ...\n\n" ) ;
      
      natomDump = natom ; 
    }
    
  
  
  }
  





  // -------------------------------> Reading ITP File ... <---------------------------------- //

  int itpinfo , natomITP ;

  if( exx != 99 )
  {
    // =====> Pre-Loading ... 
    
    printf("\nNow let's pre-load the .itp file and see how many atoms it is describing ... \n");
    
    iline = 1 ;
    
    iload = 0 ;
    
    fsearch( pitpinput , "atoms" ) ;
    
    fskip( pitpinput , 1 );

    
    
    while( ( itpinfo = freadline( buffer , MAXCHARINLINE , pitpinput , ';' ) ) != 0 )
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

    
    
    /*
    if( natomITP != natomselect )
    {
      printf("\nWhile you specified you wanted %d atoms as the solute, in your itp file, there are only %d atoms being described ...\n" , natomsoluteSelect , natomITP );
      
      printf("\nPlease check your itp file and make sure the # of atoms match ...\n");
      
      printf("\nWe will continue anyway but once we find any required info does not exist in your itp file , we will terminate this process immediately ...\n\n") ;
    }
    */

  }
  else
  {
    natomITP = natom ;  
    
    natomselect = natom ;
  }

  ITP atomdatabase[ natomITP ];
  
  if( exx != 99 )
  {
    // =====> Loading ... 
    
    printf("\nNow let's actually load the .itp file  ... \n");
    
    iline = 1 ;
    
    iload = 0 ;
    
    //ITP atomdatabase[ natomITP ];
    
    rewind( pitpinput ) ;
    
    fsearch( pitpinput , "atoms" ) ;
    
    fskip( pitpinput , 1 );
    
    while( ( itpinfo = freadline( buffer , MAXCHARINLINE , pitpinput , ';' ) ) != 0 )
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
  
  }

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 // ==============> Organizing Info and Defining Layers ... <========= //
 
 // ---------------> Allocating mem for dat array ...
 
  DAT atomcast[ natom ];

  for( iatom = 0 ; iatom < natom ; iatom ++ )
  {
    atomcast[ iatom ].cx = atomlist[ iatom ].cx * 10.00 ;
    
    atomcast[ iatom ].cy = atomlist[ iatom ].cy * 10.00 ;
    
    atomcast[ iatom ].cz = atomlist[ iatom ].cz * 10.00 ;
  
    atomcast[ iatom ].atomlabel = tellatom( atomlist[ iatom ].atomname );
    
    //atomcast[ iatom ].atommass = tellmass( atomlist[ iatom ].atomname );
    
    atomcast[ iatom ].atomcharge = 0.00 ;
    
    //printf("\nNo. %d atom ; X = %lf ; Y = %lf ; Z = %lf ; Atom = %d Mass = %lf ...\n" , iatom+1 , atomcast[ iatom ].cx , atomcast[ iatom ].cy , atomcast[ iatom ].cz , atomcast[ iatom ].atomlabel , atomcast[ iatom ].atommass );
  
  
  }

  

  // ==============> Loading information into atomcast ... <========= //

  int itype = 0 ;

  /*
  
  FILE * debug ;
  
  debug = fopen( "atomnames.deb" , "wb+" );
  
  for( iatom = 0 ; iatom < natom ; iatom ++ )
  {
    fprintf( debug , "\n%d\t%s\n" , iatom + 1 , atomlist[ iatom ].atomname );
  
  }

  fclose( debug );
  
  debug = fopen( "atomtypes.deb" , "wb+");
  
  for( itype = 0 ; itype < ntype ; itype ++ )
  {
    fprintf( debug , "%8d   %s   %lf\t\n" , atomdatabase[ itype ].nr , atomdatabase[ itype ].atomname , atomdatabase[ itype ].charge );
  }
  
  fclose( debug );

  */

  
  int info_res , info_atomname ;
  
  int ntype ;
  
  if( exx != 99 )
  {
    ntype = natomITP ;
    
    for( iatom = natomselect  ; iatom < natom ; iatom ++ )
    {
      for( itype = 0 ; itype < ntype ; itype ++ )
      {
        if( ( info_res = strcmp( atomlist[ iatom ].resname , atomdatabase[ itype ].residue ) ) == 0 && ( info_atomname = strcmp( atomlist[ iatom ].atomname , atomdatabase[ itype ].atomname ) ) == 0 )
        {
          atomcast[ iatom ].atomcharge = atomdatabase[ itype ].charge ;
          
          //printf("\nNo.%d atom , charge is %lf ...\n" , iatom + 1 , atomcast[ iatom ].atomcharge );
          
          break ;
        
        }

      }
    
    }

  }




  
  // ---------------> Outputing ...
 
  pCNDOoutput = fopen( outCNDOname , "wb+" );
  
   
  fprintf( pCNDOoutput , "%s\n" , inpgroname );

  fprintf( pCNDOoutput , "%s\n" , "HAMILT= INDO" );
  
  fprintf( pCNDOoutput , "%s\n" , "STOP= CI" );
  
  fprintf( pCNDOoutput , "%s\n" , "ROTINV= YES" );
  
  //fprintf( pCNDOoutput , "%s\n" , "SHIFT= 20" );
  
  fprintf( pCNDOoutput , "%s\n" , "BETA= INDO/S" );
  
  fprintf( pCNDOoutput , "%s\n" , "POINTGRP= C1" );
  
  fprintf( pCNDOoutput , "%s\n" , "EX_FROM= 60" );
  
  fprintf( pCNDOoutput , "%s\n" , "MAX_CI= 60" );
  
  fprintf( pCNDOoutput , "%s\n" , "CI_DUMP= 60" );
  
  //fprintf( pCNDOoutput , "%s\n" , "DUMP= MAX" );
  
  fprintf( pCNDOoutput , "%s%d\n" , "CHARGE= " , charge );
  
  fprintf( pCNDOoutput , "%s%d\n" , "MULT_CI= " , multiplicity );
  
  fprintf( pCNDOoutput , "%s\n" , "MAX_ITS= 300" );
  
  fprintf( pCNDOoutput , "%s\n\n" , "RESTART= MO" );
  
  

  
  if( natomDump >= natomselect )
  {
    for( iatom = 0 ; iatom < natomselect ; iatom ++ )
    {
      fprintf( pCNDOoutput , "%7.3f   %7.3f   %7.3f   %5d\n" , atomcast[ iatom ].cx , atomcast[ iatom ].cy , atomcast[ iatom ].cz , atomcast[ iatom ].atomlabel );
  
    }  
  }
  else
  {
    for( iatom = 0 ; iatom < natomDump ; iatom ++ )
    {
      fprintf( pCNDOoutput , "%7.3f   %7.3f   %7.3f   %5d\n" , atomcast[ iatom ].cx , atomcast[ iatom ].cy , atomcast[ iatom ].cz , atomcast[ iatom ].atomlabel );
  
    }  
  }

  
  //fprintf( pCNDOoutput , "\n\n\n" );


  if( exx != 99 )
  {
    if( natomDump > natomselect && natomDump <= natom )
    {
      for( iatom = natomselect ; iatom < natomDump ; iatom ++ )
      {
          //printf("\nPrinting out No.%d atom which is in residue : %s into CNDO input file ... \n\n" , iatom + 1 , atomlist[ iatom ].resname );
      
          fprintf( pCNDOoutput , "%7.3f   %7.3f   %7.3f   %5d     %10.6f\n" , atomcast[ iatom ].cx , atomcast[ iatom ].cy , atomcast[ iatom ].cz , -1*atomcast[ iatom ].atomlabel , atomcast[ iatom ].atomcharge );
    
      }
    }
    
  } 


  fprintf( pCNDOoutput , "\n\n\n\n\n" ) ;
   /*
  */
  
  
  // -----------------> The End ... Really???
  
  return(0);
  
  
  
  }
  
  
  
  
  
  
  
  
  
  
  











int tellatom( char * atomMDname )
{
  // -------> Initiating by buffering atom name from MD to a local char-array ...
  
  int namelength = strlen( atomMDname );
  
  char * namebuffer = malloc( ( namelength + 1 ) * sizeof( char ) ) ;
  
  strcpy( namebuffer , atomMDname );
  
  int atomlabel ;
  
  
  
  // -------------> Deciding the corresponding atom symbol for current atom  ... 
  
  char firstLetter = * namebuffer ;
  
  char secondLetter = *( namebuffer + 1 );
  
  if( secondLetter <= '9' && secondLetter >= '0' )
  {
    switch ( firstLetter )
    {
      case 'H' : atomlabel = 1 ; break ;
      
      case 'h' : atomlabel = 1 ; break ;
      
      case 'B' : atomlabel = 5 ; break ; 
      
      case 'b' : atomlabel = 5 ; break ; 
      
      case 'C' : atomlabel = 6 ; break ; 
      
      case 'c' : atomlabel = 6 ; break ; 
      
      case 'N' : atomlabel = 7 ; break ; 
      
      case 'n' : atomlabel = 7 ; break ; 
      
      case 'O' : atomlabel = 8 ; break ; 
      
      case 'o' : atomlabel = 8 ; break ;
      
      case 'F' : atomlabel = 9 ; break ;
      
      case 'f' : atomlabel = 9 ; break ;
      
      case 'P' : atomlabel = 15 ; break ;
      
      case 'p' : atomlabel = 15 ; break ;
    
      case 'S' : atomlabel = 16 ; break ;
      
      case 's' : atomlabel = 16 ; break ;
      
      case 'K' : atomlabel = 19 ; break ;
      
      case 'k' : atomlabel = 19 ; break ;
      
      default  : printf("\nUnkown Atom detected : %c... Mission aborting ...\n\n" , firstLetter );
       
                 exit(1);

    }
  } 
  else
  {
      switch ( firstLetter )
      {
        case 'S' : atomlabel = 14 ; break ; // ---> Si <--- //
        
        case 's' : atomlabel = 14 ; break ; // ---> Si <--- //
        
        case 'C' : atomlabel = 17 ; break ; // ---> Cl <--- //
        
        case 'c' : atomlabel = 17 ; break ; // ---> Cl <--- //
        
        case 'M' : atomlabel = 12 ; break ; // ---> Mg <--- //
        
        case 'm' : atomlabel = 12 ; break ; // ---> Mg <--- //
        
        case 'Z' : atomlabel = 30 ; break ; // ---> Zn <--- //
        
        case 'z' : atomlabel = 30 ; break ; // ---> Zn <--- //
        
        default  : printf("\nUnkown Atom detected : %c%c... Mission aborting ...\n\n" , firstLetter , secondLetter );
        
                   exit(1);

      }
   
   }
  
  
  
  // -------------> Finishing up ...
  
  return( atomlabel );
  
  


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


void dzeros(int dimrow, int dimcol, double *p)
{
  int i,j;
  
  for (i=0;i<dimrow;i++)
  {
    for (j=0;j<dimcol;j++)
    {
      *(p+dimcol*i+j)=0.0000;
    }
  }
  
}

void izeros(int dimrow, int dimcol, int *p)
{
  int i,j;
  
  for (i=0;i<dimrow;i++)
  {
    for (j=0;j<dimcol;j++)
    {
      *(p+dimcol*i+j)=0;
    }
  }
  
}



