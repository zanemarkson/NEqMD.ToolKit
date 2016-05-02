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



typedef struct g09input
{   
  char atomsymbol[3] ;
  double cx ;
  double cy ;
  double cz ;
  int dump ;
  
} INP ;



int main( int argc , char * argv[ ] )
{


  char ** pcmd = argv ; 
  
  int icmd ;
  
  char originalINPFileName[ 100 ] , indexFileName[ 100 ] , outputINPName[ 100 ] ;

  int len_originalINPFileName ;
  
  FILE * poriginalINPFile , * pIndexFile , * poutINPFile ; 

  char pGroupName[ 50 ] , qGroupName[ 50 ] , rGroupName[ 50 ] , xGroupName[ 50 ] ;
  
  int lmo ; // lmo is just for convenience = nbasis^2
  
  int HOMO , LUMO ;
  
  int nbasis , natom , nelectron , ncistate ;
  
  int ibasis , iatom , iao ;
  
  int iload , iline ;
  
  double done = 1.0000 ; double dzero = 0.0000 ; 
  
  int ithree = 3 ; int ione = 1 ; int izero = 0 ;
  


  char buffer[ MAXCHARINLINE ] ;
  
  int lenBuff ;
  
  char cache[ MAXCHARINLINE ] ;

  
  int itmp , itmp2 ; double dtmp ; char ctmp , tmp_char ; 
  
  char stmp[ 150 ] , stmp2[ 150 ] , tmpString[ 150 ];
  
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
      printf("* G_REORGXYZ_D : Re-Organize Coordinate With Respect To Index(ZM)    *\n");
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
  
  strcpy( originalINPFileName , "orginial.inp" ) ;
  
  strcpy( indexFileName , "system.index" ) ;
  
  strcpy( outputINPName , "reorganized.inp" ) ;
  
  
  strcpy( pGroupName , "Donor-H" ) ;
  
  strcpy( qGroupName , "Acceptor-H" ) ;
  
  strcpy( rGroupName , "Bridge-H" ) ;
  
  strcpy( xGroupName , "IRModes" ) ;
  
  
  // =====> Parsing cmd-line arguments ...
  
  int exX = YES , exR = YES ;
  
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
	    	      
	      case 'c' : strcpy( originalINPFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Input Coordinate File name : %s ...\n" , originalINPFileName ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 //exgro = 19 ;
	                 
	                 break ; 

	     
	      case 'o' : strcpy( outputINPName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Output File name : %s ...\n" , outputINPName ); 
	                 
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
	      

	      case 'h' : printf("\nUsage:  %s [ -c input xyz file name ][ -o output file name ] [ -n index file name ]\n     [ -p Group P ] [ -q Group Q ] [ -r Group R ] [ -x Group X ]\n" , * argv ); 
	                 
                     printf("\n                Note : 1) when -m is specified as \"None\" or \"none\" or \"NONE\", no reference will be checked and orbital numbers provided in CI file will be directly used.\n");
	                 
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
  
  if( ( poriginalINPFile = fopen( originalINPFileName , "r" ) ) == NULL )
  {
    printf("\nUser specified original coordinate file %s NOT FOUND ...\n" , originalINPFileName );
    
    exit( 63 ) ;
  
  }
  
  if( ( pIndexFile = fopen( indexFileName , "r" ) ) == NULL )
  {
    printf("\nUser specified index File %s NOT FOUND ...\n" , indexFileName );
    
    exit( 63 ) ;
  
  }
  
  
  // -------------------------------> Loading Coordinates from Input inp File ... <---------------------------------- //
  
  int len_xyzInput , natomInput ;


  while( ( info = freadline( buffer , MAXCHARINLINE , poriginalINPFile , '!' ) ) != 0 )
  { 
    blank_signal = stellblank( buffer ) ;
    
    if( blank_signal != 0 )
    {
      tmp_char = getfirst( buffer ) ;
      
      printf("\nLine reads : %s \n\n" , buffer ) ;
      
      if( tmp_char == '#' )
      {        
        break ;
      }
    
    }

  }
  
  
  while( ( info = freadline( buffer , MAXCHARINLINE , poriginalINPFile , '!' ) ) != 0 )
  { 
    blank_signal = stellblank( buffer ) ;
    
    if( blank_signal == 0 )  break ;
    
    
  }
  
  fskip( poriginalINPFile , 3 ) ;
  
  
  
  iload = 0 ; iline = 0 ;
  
  while( ( info = freadline( buffer , MAXCHARINLINE , poriginalINPFile , '!' ) ) != 0 )
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
      
      if( ( tmp_char = getfirst( buffer ) ) == '!' )
      {
        //printf("\nThis is a comment line ... So nothing will be loaded ...\n");
        
        //fskip( pinputfile , 1 );
        
        continue ;
      }
      else
      {
        iload ++ ;
        
      }
      
      //printf("\n%s\n" , buffer );
    }
    else
    {
      printf("\nSomething is wrong with the reading file part ...\n");
      
      exit(1);
    }
    
  }
  
  natomInput = iload ;

  
  
  
  
  
  INP originalINP[ natomInput ] ;
  
  rewind( poriginalINPFile ) ;
    
  while( ( info = freadline( buffer , MAXCHARINLINE , poriginalINPFile , '!' ) ) != 0 )
  { 
    blank_signal = stellblank( buffer ) ;
    
    printf("\nLine reads : %s \n\n" , buffer ) ;
    
    if( blank_signal != 0 )
    {
      tmp_char = getfirst( buffer ) ;
      
      if( tmp_char == '#' )
      {
        break ;    
      }
    
    }

  }
  
  
  while( ( info = freadline( buffer , MAXCHARINLINE , poriginalINPFile , '!' ) ) != 0 )
  { 
    blank_signal = stellblank( buffer ) ;
    
    if( blank_signal == 0 )  break ;
    
    
  }
  
  fskip( poriginalINPFile , 3 ) ;
  
  iload = 0 ; iline = 0 ;
  
  while( ( info = freadline( buffer , MAXCHARINLINE , poriginalINPFile , '!' ) ) != 0 )
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
      
      if( ( tmp_char = getfirst( buffer ) ) == '!' )
      {
        //printf("\nThis is a comment line ... So nothing will be loaded ...\n");
        
        //fskip( pinputfile , 1 );
        
        continue ;
      }
      else
      {
        strpickword( buffer , 1 , cache ) ;
        
        strcpy( originalINP[ iload ].atomsymbol , cache ) ;
        
        strpickword( buffer , 2 , cache ) ;
        
        originalINP[ iload ].cx = atof( cache ) ;
        
        strpickword( buffer , 3 , cache ) ;
        
        originalINP[ iload ].cy = atof( cache ) ;
        
        strpickword( buffer , 4 , cache ) ;
        
        originalINP[ iload ].cz = atof( cache ) ;
        
        

        
        iload ++ ;
        
      }
      
      //printf("\n%s\n" , buffer );
    }
    else
    {
      printf("\nSomething is wrong with the reading file part ...\n");
      
      exit(1);
    }
    
  
    if( iload == natomInput ) break ;
  
    iline ++ ;
  
  }
  
  
  

  // -------------------------------> Read index file and acquire info of P, Q, R and X ... <---------------------------------- //
  
  
  int next_atom_index , current_atom_index ;
  
  int natomPGroup , natomQGroup , natomRGroup , natomXGroup ;
  
  // =====> P Group , Pre-Loading 
  
  rewind( pIndexFile ) ;
  
  info = fsearch( pIndexFile , pGroupName ) ;
  
  if( info == 0 )
  {
    printf("\nYour index file does not contain entry with name [ %s ] \n\n" , pGroupName ) ;
    
    exit( 61 ) ;
  }
  
  fskip( pIndexFile , 1 ) ;
  
  fscanf( pIndexFile , "%s" , tmpString );
  
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
        
      fscanf( pIndexFile , "%s" , tmpString ) ;
    }
    else
    {
      itmp ++ ;
        
      current_atom_index = atoi( tmpString );
        
      //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
      fscanf( pIndexFile , "%s" , tmpString );
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
  
  fscanf( pIndexFile , "%s" , tmpString );
  
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
        
      fscanf( pIndexFile , "%s" , tmpString ) ;
    }
    else
    {
      itmp ++ ;
        
      current_atom_index = atoi( tmpString );
        
      //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
      fscanf( pIndexFile , "%s" , tmpString );
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
  
    fscanf( pIndexFile , "%s" , tmpString );
  
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
        
        fscanf( pIndexFile , "%s" , tmpString ) ;
      }
      else
      {
        itmp ++ ;
        
        current_atom_index = atoi( tmpString );
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%s" , tmpString );
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
  
    fscanf( pIndexFile , "%s" , tmpString );
  
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
        
        fscanf( pIndexFile , "%s" , tmpString ) ;
      }
      else
      {
        itmp ++ ;
        
        current_atom_index = atoi( tmpString );
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%s" , tmpString );
      }

    

    }

    natomXGroup = itmp ;
  
    printf("\nBased on index file , X group %s has %d atoms ...\n" , xGroupName , natomXGroup );
  
  }
  else
  {
    natomXGroup = 0 ;
  }

  
  
  int natomOutput = natomPGroup + natomQGroup + natomRGroup + natomXGroup ;
  
  if( natomOutput > natomInput )
  {
    printf("\nWe only have %d atoms from input ... After adding all available groups together ... we have a system with %d atoms ... Please check you index file ...\n" , natomInput , natomOutput );
    
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
        
      fscanf( pIndexFile , "%s" , tmpString ) ;
    }
    else
    {   
      current_atom_index = atoi( tmpString );
      
      *( pGroupAtoms + itmp ) = current_atom_index ;
      
      itmp ++ ;
        
      //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
      fscanf( pIndexFile , "%s" , tmpString );
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
        
      fscanf( pIndexFile , "%s" , tmpString ) ;
    }
    else
    {   
      current_atom_index = atoi( tmpString );
      
      *( qGroupAtoms + itmp ) = current_atom_index ;
      
      itmp ++ ;
        
      //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
      fscanf( pIndexFile , "%s" , tmpString );
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
        
        fscanf( pIndexFile , "%s" , tmpString ) ;
      }
      else
      {     
        current_atom_index = atoi( tmpString );
      
        *( rGroupAtoms + itmp ) = current_atom_index ;
      
        itmp ++ ;
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%s" , tmpString );
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
        
        fscanf( pIndexFile , "%s" , tmpString ) ;
      }
      else
      {   
        current_atom_index = atoi( tmpString );
      
        *( xGroupAtoms + itmp ) = current_atom_index ;
      
        itmp ++ ;
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        fscanf( pIndexFile , "%s" , tmpString );
      }

    }
  
  }

  
  
  
  // ===> Allocate Memory For Storing Re-Organized Coordinates ...
  
  
  INP reorgINP[ natomOutput ] ;
  

  
  // ===> Put It IN => P 
  
  int id ;
  
  for( iatom = 0 ; iatom < natomPGroup ; iatom ++ )
  {
    id = iatom ; 
    
    strcpy( reorgINP[ id ].atomsymbol , originalINP[ *( pGroupAtoms + iatom ) - 1 ].atomsymbol ) ;
    
    reorgINP[ id ].cx = originalINP[ *( pGroupAtoms + iatom ) - 1 ].cx ;
    
    reorgINP[ id ].cy = originalINP[ *( pGroupAtoms + iatom ) - 1 ].cy ;
    
    reorgINP[ id ].cz = originalINP[ *( pGroupAtoms + iatom ) - 1 ].cz ;
  
  }

  
  // ===> Put It IN => Q
  
  
  for( iatom = 0 ; iatom < natomQGroup ; iatom ++ )
  {
    id = iatom + natomPGroup ; 
    
    strcpy( reorgINP[ id ].atomsymbol , originalINP[ *( qGroupAtoms + iatom ) - 1 ].atomsymbol ) ;
    
    reorgINP[ id ].cx = originalINP[ *( qGroupAtoms + iatom ) - 1 ].cx ;
    
    reorgINP[ id ].cy = originalINP[ *( qGroupAtoms + iatom ) - 1 ].cy ;
    
    reorgINP[ id ].cz = originalINP[ *( qGroupAtoms + iatom ) - 1 ].cz ;

  }
  
  // ===> Put It IN => R 
  
  
  if( exR == YES )
  {
    for( iatom = 0 ; iatom < natomRGroup ; iatom ++ )
    {
      id = iatom + natomPGroup + natomQGroup ; 
    
      strcpy( reorgINP[ id ].atomsymbol , originalINP[ *( rGroupAtoms + iatom ) - 1 ].atomsymbol ) ;
    
      reorgINP[ id ].cx = originalINP[ *( rGroupAtoms + iatom ) - 1 ].cx ;
    
      reorgINP[ id ].cy = originalINP[ *( rGroupAtoms + iatom ) - 1 ].cy ;
  
      reorgINP[ id ].cz = originalINP[ *( rGroupAtoms + iatom ) - 1 ].cz ;
  
    }
  
  }
  
  // ===> Put It IN => X 
  
  
  if( exX == YES )
  {
    for( iatom = 0 ; iatom < natomXGroup ; iatom ++ )
    {
      id = iatom + natomPGroup + natomQGroup + natomRGroup ; 
    
      strcpy( reorgINP[ id ].atomsymbol , originalINP[ *( xGroupAtoms + iatom ) - 1 ].atomsymbol ) ;
    
      reorgINP[ id ].cx = originalINP[ *( xGroupAtoms + iatom ) - 1 ].cx ;
    
      reorgINP[ id ].cy = originalINP[ *( xGroupAtoms + iatom ) - 1 ].cy ;
  
      reorgINP[ id ].cz = originalINP[ *( xGroupAtoms + iatom ) - 1 ].cz ;
  
    }
  
  
  }
  
  
  // ===> Output
  
  
  poutINPFile = fopen( outputINPName , "wb+" ) ;
  
  
  rewind( poriginalINPFile ) ;
    
  while( ( info = freadline( buffer , MAXCHARINLINE , poriginalINPFile , '!' ) ) != 0 )
  { 
    blank_signal = stellblank( buffer ) ;
    
    fprintf( poutINPFile , "%s\n" , buffer ) ;
    
    if( blank_signal != 0 )
    {
      tmp_char = getfirst( buffer ) ;
      
      if( tmp_char == '#' )
      {
        break ;    
      }
    
    }

  }
  
  
  while( ( info = freadline( buffer , MAXCHARINLINE , poriginalINPFile , '!' ) ) != 0 )
  { 
    blank_signal = stellblank( buffer ) ;
    
    fprintf( poutINPFile , "%s\n" , buffer ) ;
    
    if( blank_signal == 0 )  break ;
    
    
  }
  
  iline = 0 ;
  
  while( ( info = freadline( buffer , MAXCHARINLINE , poriginalINPFile , '!' ) ) != 0 )
  { 
    fprintf( poutINPFile , "%s\n" , buffer ) ;
    
    if( iline == 2 )  break ;
    
    iline ++ ;
    
    
  }
 
 
 
  
  
  for( iatom = 0 ; iatom < natomOutput ; iatom ++ )
  {
    fprintf( poutINPFile , "%s   % 10.6f   % 10.6f   % 10.6f\n" , reorgINP[ iatom ].atomsymbol , reorgINP[ iatom ].cx , reorgINP[ iatom ].cy , reorgINP[ iatom ].cz ) ;
  
  }
  
  fprintf( poutINPFile , "\n\n\n" ) ;
  
  fclose( poutINPFile ) ;
  
  
  
  
  
  
  

  return( 0 ) ;






}



