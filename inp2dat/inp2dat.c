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
  int atomnumber ;
  int dump ;
  
} INP ;



int main( int argc , char * argv[ ] )
{


  char ** pcmd = argv ; 
  
  int icmd ;
  
  char inpFileName[ 100 ] ,  datFileName[ 100 ] ;

  int len_inpFileName , len_datFileName ;
  
  FILE * pinpFile , * pdatFile ; 

  int HOMO , LUMO ;
  
  int nbasis , natom ;
  
  int charge , multiplicity ;
  
  int natomSelect ;
  
  int ibasis , iatom , iao ;
  
  int iload , iline ;
  
  double done = 1.0000 ; double dzero = 0.0000 ; 
  
  int ithree = 3 ; int ione = 1 ; int izero = 0 ;
  
  int debuggingMode = NO ;
  

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
      printf("* G_INP2DAT_D : Change file format from G09 inp to CNDO .dat file    *\n");
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

 
  
  // =====> Those existences ...
  
  int exinp  = 16 ;
  
  int exdat = 18 ;
  
  int exselect = 20 ;
  
  
  // =====> Default Parameters ...
  
  charge = 0 ; 
  
  multiplicity = 1 ;
  
  
  
  // =====> Parsing cmd-line arguments ...
  
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
	    	      
	      case 'f' : strcpy( inpFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Input Coordinate File name : %s ...\n" , inpFileName ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 exinp = 17 ;
	                 
	                 break ; 

	     
	      case 'o' : strcpy( datFileName , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Output File name : %s ...\n" , datFileName ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 exdat = 19 ;
	                 
	                 break ; 

          
	      case 's' : strcpy( tmpString , *( ++ pcmd ) ) ;
	      
	                 if( strcmp( tmpString , "all" ) == 0 || strcmp( tmpString , "All" ) == 0 )
	                 {
	                   exselect = 22 ;
	                   
	                   printf("\nCommand-line argument indicates : All atoms will be chosen as solute ...\n" );
	                 }
	                 else
	                 {
	                   printf("\nReceived information : %s ...\n" , tmpString ) ;
	                   
	                   natomSelect = atoi( tmpString ); 
	                  
	                   exselect = 21 ;
	                   
	                   printf("\nCommand-line argument indicates : First %d atoms will be chosen as solute ...\n" , natomSelect );
	                 
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
          
	      case 'h' : printf("\nUsage:  %s [ -f input G09 .inp/.gjf/.com file name ][ -o CNDO format .dat file name ] [ -s # of atom selected as solute ] [ -g YES/yes/Yes/Y/y or NO/no/No/N/n to invoke or not invoke debugging mode ]\n" , * argv ); 
	                 
	                 printf("\nNOTE : 1) [ -s all ] or [ -s All ] indicates all atoms chosen as solute;\n\n\n");
	                 
	                 icmd = icmd + 1 ; 
	                 
	                 exit( 0 ) ;
	                 
	      

	      default : printf("\n\nInvalid option ' %s ' ... Please refer to the usage by typing ' %s -h '\n\n" , flag , * argv ); 
	      
	                icmd = argc ; 
	                
	                exit( 1 );

      
      
      
      }
    
    }
    else
    {
        printf("\n\nInvalid option ' %s ' ... Please refer to the usage by typing ' %s -h '\n\n" , flag , * argv );

	    exit(1);
      
      
    }
    
 
  } 
  

  
  // =====> Taking care of file names ... <===== //
  
  int fileNameIndicator = exinp * exdat ;
  
  switch( fileNameIndicator )
  {
    case 16 * 18 : strcpy( inpFileName , "structure.inp" ) ;
    
                   strcpy( datFileName , "structure.dat" ) ;
    
                   printf("\nNeither G09 .inp file name nor CNDO .dat file name provided ... Trying default names [ %s ] & [ %s ] ...\n\n" , inpFileName , datFileName ) ;
  
                   break ;
                   
    case 16 * 19 : len_datFileName = strlen( datFileName ) ;
                   
                   strncpy( inpFileName , datFileName , len_datFileName - 4 ) ;
                   
                   *( inpFileName + len_datFileName - 4 ) = '\0' ;
                   
                   printf("\nOutput G09 .inp file name NOT provided ... The default name based on input CNDO .dat file name is [ %s ]...\n" , inpFileName ) ;
                   
                   break ;
                   
    case 17 * 18 : len_inpFileName = strlen( inpFileName ) ;
    
                   strncpy( datFileName , inpFileName , len_inpFileName - 4 ) ;
                   
                   *( datFileName + len_inpFileName - 4 ) = '\0' ;
                   
                   printf("\nOutput CNDO .dat file name NOT provided ... The default name based on input G09 .inp file name is [ %s ]...\n" , datFileName ) ;
                   
                   break ;
                   
    case 17 * 19 : printf("\nGood! Both G09 .inp and CNDO .dat file names are provided ...[ %s ] & [ %s ] \n\n" , inpFileName , datFileName ) ;
    
                   break ;
  
  
  }
  
  
  
  
  
  // =====> Verify File Access <===== //
  
  if( ( pinpFile = fopen( inpFileName , "r" ) ) == NULL )
  {
    printf("\nUser specified original coordinate file %s NOT FOUND ...\n" , inpFileName );
    
    exit( 63 ) ;
  
  }
  
  
  // =====> Pre-Loading Coordinates from Input inp File ... <===== //
  
  int len_xyzInput , natomInput ;


  while( ( info = freadline( buffer , MAXCHARINLINE , pinpFile , '!' ) ) != 0 )
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

  } // Searching for the 1st line of #Route section ...
  
  
  while( ( info = freadline( buffer , MAXCHARINLINE , pinpFile , '!' ) ) != 0 )
  { 
    blank_signal = stellblank( buffer ) ;
    
    if( blank_signal == 0 )  break ;

  }  // Searching for the end of #Route section ...
  
  
  while( ( info = freadline( buffer , MAXCHARINLINE , pinpFile , '!' ) ) != 0 )
  { 
    blank_signal = stellblank( buffer ) ;
    
    if( blank_signal == 0 )  break ;

  }  // Searching for the end of Comment card section ...
  
  
  /*
  if( ( info = freadline( buffer , MAXCHARINLINE , pinpFile , '!' ) ) != 0 )
  {
    strpickword( buffer , 1 , cache ) ;
    
    charge = atoi( cache ) ;
    
    strpickword( buffer , 2 , cache ) ;
    
    multiplicity = atoi( cache ) ;
  
  }
  */


  fskip( pinpFile , 1 ) ;
  
  
  
  iload = 0 ; iline = 0 ;
  
  while( ( info = freadline( buffer , MAXCHARINLINE , pinpFile , '!' ) ) != 0 )
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
      
      exit( 77 );
    }
    
  }
  
  natomInput = iload ;

  printf("\nYour input G09 .inp file provides geometry for %d atoms ...\n" , natomInput ) ;
  
  
  
  if( exselect == 20 )
  {
    printf("\nBy default , all %d atoms are selected as solute ... \n\n" , natomInput ) ;
    
    natomSelect = natomInput ;
  
  }
  else if( exselect == 22 )
  {
    printf("\nPer user's request , all %d atoms are selected as solute ... \n\n" , natomInput ) ;
    
    natomSelect = natomInput ;
  }
  else if( exselect == 21 )
  {
    if( natomSelect > natomInput )
    {
      printf("\nOops ... you requested %d atoms as solute but we only have %d atoms in G09 .inp file ... By default , all atoms are selected as solute ... \n\n" , natomSelect , natomInput ) ;
      
      natomSelect = natomInput ;

    }
    else if( natomSelect <= natomInput && natomSelect > 0 )
    {
      printf("\nPer user's request , the first %d atoms are selected as solute ... \n\n" , natomSelect ) ;
    }
    else if( natomSelect <= 0 )
    {
      printf("\n# of solute atom must be a positive integer!!!\n\n") ;
      
      exit( 52 ) ;
    }
  
  
  }
  
  
  
  // =====> Actually Loading Coordinates from Input inp File ... <===== //
  
  INP originalINP[ natomInput ] ;
  
  rewind( pinpFile ) ;
    
  while( ( info = freadline( buffer , MAXCHARINLINE , pinpFile , '!' ) ) != 0 )
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

  } // Searching for the 1st line of #Route section ...
  
  
  while( ( info = freadline( buffer , MAXCHARINLINE , pinpFile , '!' ) ) != 0 )
  { 
    blank_signal = stellblank( buffer ) ;
    
    if( blank_signal == 0 )  break ;

  }  // Searching for the end of #Route section ...
  
  
  while( ( info = freadline( buffer , MAXCHARINLINE , pinpFile , '!' ) ) != 0 )
  { 
    blank_signal = stellblank( buffer ) ;
    
    if( blank_signal == 0 )  break ;

  }  // Searching for the end of Comment card section ...
  
  
  
  if( ( info = freadline( buffer , MAXCHARINLINE , pinpFile , '!' ) ) != 0 )
  {
    strpickword( buffer , 1 , cache ) ;
    
    charge = atoi( cache ) ;
    
    strpickword( buffer , 2 , cache ) ;
    
    multiplicity = atoi( cache ) ;
    
    printf("\nSystem under study has : charge = [ %d ] , multiplicity = [ %d ] ... \n\n" , charge , multiplicity ) ;
  
  }


  //fskip( pinpFile , 1 ) ;
  
  
  iload = 0 ; iline = 0 ;
  
  while( ( info = freadline( buffer , MAXCHARINLINE , pinpFile , '!' ) ) != 0 )
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
        
        printf("\nWe have an atom [ %s ] ...\n" , cache ) ;
        
        strpickword( buffer , 2 , cache ) ;
        
        originalINP[ iload ].cx = atof( cache ) ;
        
        strpickword( buffer , 3 , cache ) ;
        
        originalINP[ iload ].cy = atof( cache ) ;
        
        strpickword( buffer , 4 , cache ) ;
        
        originalINP[ iload ].cz = atof( cache ) ;
        
        
        originalINP[ iload ].atomnumber = 0 ;
        

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
  
  if( debuggingMode == YES )
  {
    debug = fopen( "all.coordinates.deb" , "wb+" ) ;
  
    for( iatom = 0 ; iatom < natomSelect ; iatom ++ )
    {
      fprintf( debug , "%s	%7.3f   %7.3f   %7.3f   %5d\n" , originalINP[ iatom ].atomsymbol , originalINP[ iatom ].cx , originalINP[ iatom ].cy , originalINP[ iatom ].cz , originalINP[ iatom ].atomnumber );
    }

  
    fclose( debug ) ;
  }

  // =====> Figuring out the atom number of all atoms ... <===== //
  
  for( iatom = 0 ; iatom < natomInput ; iatom ++ )
  {
    originalINP[ iatom ].atomnumber = tellG09atom( originalINP[ iatom ].atomsymbol ) ;
  }
  
  
  
  // =====> Output ... <===== //
  

  pdatFile = fopen( datFileName , "wb+" ) ;

  fprintf( pdatFile , "%s\n" , inpFileName );
  
  fprintf( pdatFile , "%s\n" , "HAMILT= INDO" );
  
  fprintf( pdatFile , "%s\n" , "STOP= SCF" );
  
  fprintf( pdatFile , "%s\n" , "ROTINV= YES" );
  
  //fprintf( pdatFile , "%s\n" , "SHIFT= 20" );
  
  fprintf( pdatFile , "%s\n" , "BETA= INDO/S" );
  
  fprintf( pdatFile , "%s\n" , "POINTGRP= C1" );
  
  fprintf( pdatFile , "%s\n" , "EX_FROM= 60" );
  
  fprintf( pdatFile , "%s\n" , "MAX_CI= 60" );
  
  fprintf( pdatFile , "%s\n" , "CI_DUMP= 60" );
  
  fprintf( pdatFile , "%s%d\n" , "CHARGE= " , charge );
  
  fprintf( pdatFile , "%s%d\n" , "MULT_CI= " , multiplicity );
  
  fprintf( pdatFile , "%s\n" , "MAX_ITS= 300" );
  
  fprintf( pdatFile , "%s\n\n" , "RESTART= AUTO" );

  
  
  for( iatom = 0 ; iatom < natomSelect ; iatom ++ )
  {
    fprintf( pdatFile , "%7.3f   %7.3f   %7.3f   %5d\n" , originalINP[ iatom ].cx , originalINP[ iatom ].cy , originalINP[ iatom ].cz , originalINP[ iatom ].atomnumber );
  }

  
  fprintf( pdatFile , "\n\n\n" ) ;

  fclose( pdatFile ) ;


}





int tellG09atom( char * atomname )
{
  // -------> Initiating by buffering atom name from MD to a local char-array ...
  
  int namelength = strlen( atomname );
  
  char * namebuffer = malloc( ( namelength + 1 ) * sizeof( char ) ) ;
  
  strcpy( namebuffer , atomname );
  
  int atomlabel ;
  
  
  
  // -------------> Deciding the corresponding atom symbol for current atom  ... 
  
  char firstLetter = * namebuffer ;
  
  char secondLetter ; 
  
  
  if( namelength == 1 )
  {
    if( firstLetter > '9' || firstLetter < '0' )
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
       
                   exit( 79 );

      }
    } 
    else
    {
      atomlabel = atoi( namebuffer ) ;
    }
  
  }
  else if( namelength == 2 )
  {
    secondLetter = *( namebuffer + 1 ) ;
    
    if( firstLetter > '9' || firstLetter < '0' )
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
        
                   exit( 79 );

      }
    
    }
    else
    {
      atomlabel = atoi( namebuffer ) ;
    }
  
  }
  else
  {
    printf("\nUnkown Atom detected : [ %s ] ... Mission aborting ...\n\n" , namebuffer ) ;
    
    exit( 79 ) ;
  
  }
  
  
  // -------------> Finishing up ...
  
  return( atomlabel );
  
  


}


