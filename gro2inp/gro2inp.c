#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
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

typedef struct g09form
{
  double cx ;
  double cy ;
  double cz ;
  //char atomsymbol[3];
  int atomlabel ;

} G09 ;




int main( int argc, char * argv[] )
{
  
  // =====> Some functions ...
  
  int tellatom( char * atomMDname ) ;
  
  
  // =====> Some Constants ...

  double done = 1.0000;

  double dzero = 0.0000;

  time_t current_time;

  time( &current_time );

  char now[ 300 ] ;

  strcpy( now , ctime( &current_time ) );

  int lennow = strlen( now ) ;

  *( now + lennow - 1 ) = ' ';

  // =====> Declaring variables ... 

  int max_cmd_args = 40 ;

  if( argc > max_cmd_args )
  {
    printf("\nToo many command-line arguments ... \n");

    exit( 6 );
  }
   
  int icmd ;

  char ** pcmd ; pcmd = argv ;

  printf("\n****************************************************************\n");
    printf("*  G_GRO2INP_D : GMX .gro file to Gaussian .inp file converter *\n");
    printf("*                                                              *\n");
  printf("*  ");
  for( icmd = 0 ; icmd < argc ; icmd ++ )
  {
    printf("%s " , *( pcmd + icmd ) );
  }
  printf("\n");
    printf("*                                                              *\n");
    printf("*                                                              *\n");
    printf("* Current Time : %s                     *\n" , now );
    printf("*                                                              *\n");
    printf("*                                                              *\n");
    printf("****************************************************************\n");




  FILE * pgro ;
  
  FILE * pg09 ;
  
  /////FILE * pfreq , * pdxdr , * pmass , * phess ;
  
  FILE * debug ;
  
  char inpgroname[ 100 ] , outg09name[ 100 ] ;
  
  char defoutname[ 100 ];
  
  char method[ 200 ];
  
  int exgro , exg09 , exselect , exmethod ;

  int natom , ncart , nmode , ntraj ;
  
  int natomselect , ncartselect , nmodeselect ;

  double * cart_q ;


  int itmp , iatom , icart , imode , itraj ; 
  
  double dtmp ;
  
  char ctmp ;
  
  char tmpString[ MAXCHARINLINE ] ;

  // =====> Parsing command line input arguments ...

  // -----> Defaults ...
  
  strcpy( inpgroname , "sys.gro" ) ;
  
  strcpy( method , "#P ZIndo( Singlets , NStates = 20 ) NoSymm pop=full test" ) ;

  
  // -----> Command-Line arguments ...
  
  exgro = 18 ; exg09 = 20 ; exselect = 26 ; exmethod = 28 ;

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

	      case 'c' : strcpy( inpgroname , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Input .gro File name : %s ...\n" , inpgroname ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 exgro = 19 ;
	                 
	                 break ; 
	                 
	      case 'o' : strcpy( outg09name , *( ++ pcmd ) ) ; 
			 
			         printf("\nCommand-line argument indicates : Input dxdr File name : %s ...\n" , outg09name ); 
	      
	                 icmd = icmd + 2 ; 
	                 
	                 exg09 = 21 ;
	                 
	                 break ;


	      case 's' : strcpy( tmpString , *( ++ pcmd ) ) ;
	      
	                 if( strcmp( tmpString , "all" ) == 0 || strcmp( tmpString , "ALL" ) == 0 || strcmp( tmpString , "All" ) == 0 )
	                 {
	                   printf("\nCommand-line argument indicates : ALL ATOMS are selected for generating G09 input ... \n" ); 
	                   
	                   exselect = 28 ;
	                 
	                 }
	                 else
	                 {
	                   natomselect = atoi( tmpString ) ;
	      
	                   printf("\nCommand-line argument indicates : The first %d atoms are selected for generating G09 input ... \n" , natomselect ); 
	                   
	                   exselect = 27 ;
	                 
	                 }
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;
	                 
	      case 'L' : strcpy( method , *( ++ pcmd ) ) ;
	      
	                 printf("\nCommand-line argument indicates : QC method is ===> %s <===\n" , method );
	                 
	                 icmd = icmd + 2 ;
	                 
	                 exmethod = 29 ;
	                 
	                 break ;

	                 
	      case 'h' : printf("\nUsage:  %s [ -c 'input .gro file name' ] [ -o 'Output G09 .inp file name'] [ -s # of atom selected ] [ -L 'Method used in G09 QC calculation']\n\n" , * argv ); 
	                 
	                 printf("\n        NOTE : 1) when -s option is specified as [ all / ALL / All ] , then all atoms will be selected for generating G09 input. This is also the default. \n\n") ;
	                 
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
  

  // =====> Explain the File Name Situations ...
  
  if( exgro == 18 )
  {
    printf("\nNo Name of input .gro file specified, go with %s ...\n" , inpgroname );
  }
  
  int inputnamelength = strlen( inpgroname ) ;
  
  if( exg09 == 20 )
  {
    strncpy( defoutname, inpgroname , inputnamelength - 4 ) ;
  
    *( defoutname + inputnamelength - 4 ) = '\0' ;
  
    strcat( defoutname, ".inp") ;

    strcpy( outg09name , defoutname );
    
    printf("\nNo Name of output G09 .inp file specified, go with %s ...\n" , outg09name );

  }







  // =====> pre-Loading .gro file to get NAtom info ...  
  
  char grotitlestring[MAXLINE];
  
  int iline = 3 ;
  
  int iload = 0 ;
  
  int blank_signal , groinfo ;

  char buffer[ MAXCHARINLINE ] ;
  
  char cache[ MAXCHARINLINE ] ;

  char tmp_char ;
  
  int natomgroline , natomgrotitle ;
  
  int exVelocity = NO ;
  
  
  if( ( pgro = fopen( inpgroname , "r" ) ) == NULL )
  {
    printf("\nUser defined .gro file does not exist ... \n");
   
    exit( 3 );
  }
  else
  {

    rewind( pgro );
    
    printf("\nCurrent character is %c ... \n" , fgetc( pgro ) );
    
    fskip( pgro , 1 );

    fscanf( pgro , "%d" , &natomgrotitle );
    
    fskip( pgro , 1 );
    
    printf("\n Second line of .gro file says it is describing %d atoms ... \n\n" , natomgrotitle );


    while( ( groinfo = freadline( buffer , MAXCHARINLINE , pgro , ';' ) ) != 0 )
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

  rewind( pgro ) ;
    
  fskip( pgro , 2 ) ;


  while( ( groinfo = freadline( buffer , MAXCHARINLINE , pgro , ';' ) ) != 0 )
  { 
    printf("\n//--------------> WORKING ON NO. %d LINE ... <-------------//\n" , iline );
    
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

  natomgroline = iload - 1 ;
  
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
  
  rewind( pgro );



  ncart = 3 * natom ;
 
  switch ( natom )
  {
    case 1 : nmode = 0 ; break;

    case 2 : nmode = 1 ; break;

    default : nmode = ncart - 0 ; break;
 
  }
  
  if( exselect == 26 )
  {
    natomselect = natom ;
  }
  else if( exselect == 28 )
  {
    natomselect = natom ;
  }
  else if( exselect == 27 && natomselect > natom )
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
  else if( exselect == 27 && natomselect <= natom )
  {
    printf("\nYou have selected %d atoms for NEIC ... There are %d atoms en toto in this system ... \n" , natomselect , natom );
  
  }
  else
  {
    printf("\nSomething is wrong with the atom selection process ... NAtom = %d , NAtomSelect = %d ... \n" , natom , natomselect );
    
    exit( 78 );
  }
  
  
  ncartselect = 3 * natomselect ;
  
  switch ( natomselect )
  {
    case 1 : nmodeselect = 0 ; break;

    case 2 : nmodeselect = 1 ; break;

    default : nmodeselect = ncartselect - 6 ; break;
 
  }
  

 
  // =====> Allocating memories ... 

  cart_q = calloc( ncartselect , sizeof(double) ); 
  
  dzeros( ncartselect , 1 , cart_q );



  // =====> Loading Configuration for selected atoms , i.e. Loading .gro file ... Yikes ... 

  GRO atomlist[ natomselect ] ; 
  
  rewind( pgro );
  
  fskip( pgro , 2 );
  
  iload = 0 ; iline = 0 ;
  
  while( ( groinfo = freadline( buffer , MAXCHARINLINE , pgro , ';' ) ) != 0 )
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
    
    if( iload == natomselect ) break ;
    
    
    iline ++ ;

  }
  
  
  rewind( pgro );
  
  fskip( pgro , natomgroline + 1 );

  double boxvector[3];

  fscanf( pgro , "%lf" , boxvector + 0 );

  fscanf( pgro , "%lf" , boxvector + 1 );

  fscanf( pgro , "%lf" , boxvector + 2 );

  rewind( pgro );
  
  
  debug = fopen( "geometry.only.nm" , "wb+");
  
  for( iatom = 0 ; iatom < natomselect ; iatom ++ )
  {
    *( cart_q + iatom * 3 + 0 ) = atomlist[ iatom ].cx ; // Convert into Bohr ;
    
    *( cart_q + iatom * 3 + 1 ) = atomlist[ iatom ].cy ;
    
    *( cart_q + iatom * 3 + 2 ) = atomlist[ iatom ].cz ;
    
    fprintf( debug , "% 16.12f\t% 16.12f\t% 16.12f\n\n" , *( cart_q + iatom * 3 + 0 ) , *( cart_q + iatom * 3 + 1 ) , *( cart_q + iatom * 3 + 2 ) );
  }
  
  fclose(debug);
  
  
  
  
  G09 atomcast[ natomselect ];

  for( iatom = 0 ; iatom < natomselect ; iatom ++ )
  {
    atomcast[ iatom ].cx = atomlist[ iatom ].cx * 10.00 ;
    
    atomcast[ iatom ].cy = atomlist[ iatom ].cy * 10.00 ;
    
    atomcast[ iatom ].cz = atomlist[ iatom ].cz * 10.00 ;
  
    atomcast[ iatom ].atomlabel = tellatom( atomlist[ iatom ].atomname );
    
    //printf("\nNo. %d atom ; X = %lf ; Y = %lf ; Z = %lf ; Atom = %d Mass = %lf ...\n" , iatom+1 , atomcast[ iatom ].cx , atomcast[ iatom ].cy , atomcast[ iatom ].cz , atomcast[ iatom ].atomlabel , atomcast[ iatom ].atommass );
  
  
  }


  
  
  // =====> Deal with the Gaussian input stream ...
  
  char chkname[ 100 ] , rwfname[ 100 ] ;
  
  int outg09namelength = strlen( outg09name );
  
  
  strncpy( chkname , outg09name , outg09namelength - 4 ) ;
  
  *( chkname + outg09namelength - 4 ) = '\0' ;
  
  strcat( chkname, ".chk") ;
  
  
  strncpy( rwfname , outg09name , outg09namelength - 4 ) ;
  
  *( rwfname + outg09namelength - 4 ) = '\0' ;
  
  strcat( rwfname, ".rwf") ;
  
  
  if( exmethod == 28 )
  {
    printf("\nNo QC method specified in command-line , go with ===> %s <===\n" , method );
  }
  
  
  
  
  
  pg09 = fopen( outg09name , "wb+" );
  
  fprintf( pg09 , "%%chk=%s\n%%rwf=%s\n%%mem=4GB\n%%nprocshared=2\n\n" , chkname , rwfname );
  
  fprintf( pg09 , "%s\n\n" , method );
  
  fprintf( pg09 , "%s corresponding QC\n\n" , inpgroname );
  
  fprintf( pg09 , "0\t1\t\n" );
  
  for( iatom = 0 ; iatom < natomselect ; iatom ++ )
  {
    fprintf( pg09 , "%d\t% 16.12f\t% 16.12f\t% 16.12f\n" , atomcast[ iatom ].atomlabel , atomcast[ iatom ].cx , atomcast[ iatom ].cy , atomcast[ iatom ].cz );
  
  }

  
  
  
  int nwords , signal ;
  
  nwords = inLineWC( method ) ;

  signal = 0 ;


  for( itmp = 0 ; itmp < nwords ; itmp ++ )
  {
    strpickword( method , itmp + 1 , buffer ) ;

    if( strcmp( buffer , "dftba") == 0 || strcmp( buffer , "DFTBA") == 0 )
    {
      signal = 1 ;

      break ;
    }
  
  }
  
  
  if( signal == 1 )
  {
    fprintf( pg09 , "\n@GAUSS_EXEDIR:dftba.prm\n\n" );
  }

  fprintf( pg09 , "\n\n\n" ) ;


  fclose( pg09 );
  
  
  







/*-o-o-o-o-o-o-o ... Em...What else? ... o-o-o-o-o-o-o-o-*/



return( 0 );

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



