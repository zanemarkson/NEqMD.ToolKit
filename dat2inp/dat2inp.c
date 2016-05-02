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


int main( int argc, char * argv[] )
{
  
  FILE * pDATinput , * pGJFoutput ;

  char inpDATname[ 100 ] , outGJFname[ 100 ] ;
  
  char method[ 200 ] ;

  int multiplicity , charge ;

  int natom , ncart , natomselect , natomPrint ; 
  
  int iatom ;


  double dtmp ; 
  
  double dtmpArray[ 100 ] ;

  char buffer[ MAXCHARINLINE ] ;
  
  char cache[ MAXCHARINLINE ] ;
  
  char tmpString[ 150 ] ;
  
  char tmp_char ;

  int itmp , itmp2 ; 
  
  

  // --------> Declaring utility functions ...


  // --------> Some default values  ... 

  multiplicity = 1 ;
  
  charge = 0 ;
  
  strcpy( method , "#P ZIndo( Singlets , NStates = 15 ) NoSymm Pop=Full Test" ) ;
  

    // ==============> Handling the file names and Charge & Multiplicity ... <========= //
  
  
  // -------> Parsing the Command Line Arguments ... 

  char ** pcmd ;

  pcmd = argv ; //pcmd ++ ; 

  int icmd = 1 ;
  
  
  //int exn = 10 ; int exr  = 16 ; int exR = 18 ; int exH  = 22 ; int exM = 28 ; int exL = 30 ; 

  int exf = 1 ; int exo = 3 ;
  
  int exs = 18 ; int exmethod = 28 ;

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
	      
	      case 'f' : strcpy( inpDATname , *( ++ pcmd ) ) ; 
			 
                     printf("\nCommand-line argument indicates : Input File name : %s ...\n" , inpDATname ); 
	      
	                 exf = 7 ; 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ;

	      case 'o' : strcpy( outGJFname , *( ++ pcmd ) ) ;
	      
	                 printf("\nCommand-line argument indicates : Output File name : %s ...\n" , outGJFname ); 
	                 
	                 exo = 9 ; 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ; 
	                 
         
	      case 'L' : strcpy( method , *( ++ pcmd ) ) ;
	      
	                 printf("\nCommand-line argument indicates : G09 Method is : %s ...\n" , method ); 
	                 
	                 exmethod = 29 ;
	                 
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
	
	
	      case 'h' : printf("\nUsage:  %s [ -f 'input dat file name' ] [(optional) -o 'generated g09 input file name' ] [ -s # of atoms chosen as the solute ] [ -L Method for G09 Calculations ]\n\n" , * argv ); 
	                 
	                 printf("\nNOTE : 1) [ -s all ] or [ -s All ] indicates all atoms chosen as solute;\n\n       2) Default for -s is all atoms when nresidue = 1  or natom in 1st residue when nresidue != 1 \n");
	                 
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

  
  char defGJFname[100] , defDATname[100] ;

  int inputnamelength , outputnamelength ;
  
  switch( exf * exo )  // File Names ... int exf = 1 / 7 ; int exo = 3 / 9; 
  {
    case 1 * 3 : strcpy( inpDATname , "sys.dat" ); 
    
                 strcpy( outGJFname , "sys.inp" ) ;
                 
                 break ;
                 
                 
    case 1 * 9 : outputnamelength = strlen( outGJFname ) ;   
     
                 strncpy( defDATname, outGJFname , outputnamelength - 4 ) ;
             
                 *( defDATname + outputnamelength - 4 ) = '\0' ;
             
                 strcat( defDATname, ".dat") ;
             
                 strcpy( inpDATname , defDATname );
                 
                 break ;
                 
                 
    case 7 * 3 : inputnamelength = strlen( inpDATname ) ; 
      
                 strncpy( defGJFname, inpDATname , inputnamelength - 4 ) ;
             
                 *( defGJFname + inputnamelength - 4 ) = '\0' ;
             
                 strcat( defGJFname, ".inp") ;
              
                 strcpy( outGJFname , defGJFname ) ;
                 
                 break ;
                 
    case 7 * 9 : printf("\n\nHoorayyyyyyyy ... All file names are specified !!!\n\n");
     
                 break ;
                 

  
  }


  // --------------> Summarizing and File Access ...
  
  printf("\n\nInput CNDO .dat file name is : %s ...\n" , inpDATname );
  
  printf("\nOutput G09 inp file name is %s ...\n" , outGJFname );
    
  
  

  if( ( pDATinput = fopen( inpDATname , "r" ) ) == NULL )
  {
    printf("\nUser defined .dat file %s does not exist ... \n" , inpDATname );
   
    exit( 3 );
  }
  
  // pGJFoutput = fopen( outGJFname , "wb+" ) ;

  // -------------->  Reading information from CNDO .dat file ...

  int iline = 0 , nRouteLines ;
  
  int iload = 0 , nload ;
  
  int blank_signal , info ;
  
  int isoluteAtom = 0 , nsoluteAtom = 0 ;
  
  int ipointCharge = 0 , npointCharge  = 0 ;
  
  
  double * molecularSpecification ;


  if( fsearch( pDATinput , "CHARGE=" ) == 1 )
  {
     fscanf( pDATinput , "%s" , cache ) ;
     
     charge = atoi( cache ) ;
  
  }

  rewind( pDATinput ) ;
  
  if( fsearch( pDATinput , "MULT_CI=" ) == 1 )
  {
     fscanf( pDATinput , "%s" , cache ) ;
     
     multiplicity = atoi( cache ) ;
  
  }

  printf("\nCharge is %d , Multiplicity is %d ...\n" , charge , multiplicity ) ;
  
  rewind( pDATinput ) ;
  
  
  while( ( info = freadline( buffer , MAXCHARINLINE , pDATinput , '$' ) ) != 0 )
  {
    blank_signal = stellblank( buffer ) ;
    
    printf("\nLine reads : %s ...\n" , buffer );
    
    if( blank_signal == 0 )
    {
      break ;
    } 
    else
    {
      iline ++ ;
    }
     
  }
  
  nRouteLines = iline ;
  
  
  
  while( ( info = freadline( buffer , MAXCHARINLINE , pDATinput , '$' ) ) != 0 )
  {
    blank_signal = stellblank( buffer ) ;
    
    if( blank_signal == 0 )
    {
      //printf("\nNo.%d line is a blank line ... Moving on ...\n" , iline ) ;
      
      continue ;
    }
    else if( blank_signal == 1 )
    {  
      //printf("\nNo.%d line is NOT a blank line ... loading ...\n" , iline );
      
      if( ( tmp_char = getfirst( buffer ) ) == '$' )
      {
        //printf("\nThis is a comment line ... So nothing will be loaded ...\n");
        
        //fskip( pinputfile , 1 );
        
        continue ;
      }
      else
      {
        printf("\nLine reads : %s ...\n" , buffer );
        
        itmp = inLineWC( buffer ) ;
        
        if( itmp == 4 )
        {
          isoluteAtom ++ ;
        
        }
        else if( itmp == 5 )
        {
          ipointCharge ++ ;
        }
        else
        {
          continue ;
        }
        
      }
      
      //printf("\n%s\n" , buffer );
    }
    else
    {
      printf("\nSomething is wrong with the reading file part ...\n");
      
      exit(1);
    }
    
    iload ++ ;
  }
  
  nload = iload ; nsoluteAtom = isoluteAtom ; npointCharge = ipointCharge ;
  
  printf("\n En TOTO %d coordinates and other specifications loaded ... %d are atoms and %d are point charges ... \n" , nload , isoluteAtom , ipointCharge ) ;
  
  
  
  
  DAT atomlist[ nsoluteAtom + npointCharge ] ;
  
  rewind( pDATinput ) ;
  
  fskip( pDATinput , nRouteLines ) ;
  
  iload = 0 ; isoluteAtom = 0 ; ipointCharge = 0 ;
  
  
  while( ( info = freadline( buffer , MAXCHARINLINE , pDATinput , '$' ) ) != 0 )
  {
    blank_signal = stellblank( buffer ) ;
    
    if( blank_signal == 0 )
    {
      //printf("\nNo.%d line is a blank line ... Moving on ...\n" , iline ) ;
      
      continue ;
    }
    else if( blank_signal == 1 )
    {  
      //printf("\nNo.%d line is NOT a blank line ... loading ...\n" , iline );
      
      if( ( tmp_char = getfirst( buffer ) ) == '$' )
      {
        //printf("\nThis is a comment line ... So nothing will be loaded ...\n");
        
        //fskip( pinputfile , 1 );
        
        continue ;
      }
      else
      {
        //printf("\nLine reads : %s ...\n" , buffer );
        
        itmp = inLineWC( buffer ) ;
        
        if( itmp == 4 )
        {
          strpickword( buffer , 1 , cache ) ;
          
          atomlist[ isoluteAtom ].cx  = atof( cache ) ;
          
          strpickword( buffer , 2 , cache ) ;
          
          atomlist[ isoluteAtom ].cy  = atof( cache ) ;
          
          strpickword( buffer , 3 , cache ) ;
          
          atomlist[ isoluteAtom ].cz  = atof( cache ) ;
          
          strpickword( buffer , 4 , cache ) ;
          
          atomlist[ isoluteAtom ].atomlabel  = atoi( cache ) ;
          
          printf("\n No. %d atom ...\n" , isoluteAtom ) ;
          
          isoluteAtom ++ ;
        
        }
        else if( itmp == 5 )
        {
          strpickword( buffer , 1 , cache ) ;
          
          atomlist[ nsoluteAtom + ipointCharge ].cx  = atof( cache ) ;
          
          strpickword( buffer , 2 , cache ) ;
          
          atomlist[ nsoluteAtom + ipointCharge ].cy  = atof( cache ) ;
          
          strpickword( buffer , 3 , cache ) ;
          
          atomlist[ nsoluteAtom + ipointCharge ].cz  = atof( cache ) ;
          
          strpickword( buffer , 4 , cache ) ;
          
          atomlist[ nsoluteAtom + ipointCharge ].atomlabel  = -1 * atoi( cache ) ;
          
          strpickword( buffer , 5 , cache ) ;
          
          atomlist[ nsoluteAtom + ipointCharge ].atomcharge = atof( cache ) ;
          
          printf("\n No. %d atom ...\n" , nsoluteAtom + ipointCharge ) ;
                    
          ipointCharge ++ ;
        }
        else
        {
          continue ;
        }
        
      }
      //printf("\n%s\n" , buffer );
    }
    else
    {
      printf("\nSomething is wrong with the reading file part ...\n");
      
      exit(1);
    }
    
    iload ++ ;
  }
  
  
  
  
  
  if( iload == nload && ipointCharge == npointCharge && isoluteAtom == nsoluteAtom )
  {
    printf("\n En TOTO %d coordinates and other specifications loaded ... %d are atoms and %d are point charges ... \n" , nload , nsoluteAtom , npointCharge ) ;
  }
  else
  {
    printf("\nEh...oh... Something is wrong during the loading ... \n");
    
    printf("\n En TOTO %d coordinates and other specifications loaded ... %d are atoms and %d are point charges ... \n" , nload , nsoluteAtom , npointCharge ) ;
    
    exit( 515 ) ;
  }
  
  
  natom = nsoluteAtom + npointCharge ;
  
  switch ( exs )
  {
    case 19 : natomPrint = natomselect ; break ;
    
    case 18 : natomPrint = nsoluteAtom ; break ;
    
    case 20 : natomPrint = natom ; break ;
    
    default : printf("\nSomething is wrong during the printing-out ...\n") ; exit( 538 ) ;
  
  }  
  
  
  
  // -------------->  Outputing ... 
  
  char chkname[ 100 ] , rwfname[ 100 ] ;
  
  outputnamelength = strlen( outGJFname ) ;   
  
  
  strncpy( chkname , outGJFname , outputnamelength - 4 ) ;
  
  *( chkname + outputnamelength - 4 ) = '\0' ;
  
  strcat( chkname, ".chk") ;
  
  
  strncpy( rwfname , outGJFname , outputnamelength - 4 ) ;
  
  *( rwfname + outputnamelength - 4 ) = '\0' ;
  
  strcat( rwfname, ".rwf") ;
  
  
  if( exmethod == 28 )
  {
    printf("\nNo QC method specified in command-line , go with ===> %s <===\n" , method );
  }
  

  
  pGJFoutput = fopen( outGJFname , "wb+" );
  
  fprintf( pGJFoutput , "%%chk=%s\n%%rwf=%s\n%%mem=4GB\n%%nprocshared=2\n\n" , chkname , rwfname );
  
  fprintf( pGJFoutput , "%s\n\n" , method );
  
  fprintf( pGJFoutput , "%s corresponding QC\n\n" , inpDATname );
  
  fprintf( pGJFoutput , "%d\t%d\t\n" , charge , multiplicity );
  
  for( iatom = 0 ; iatom < natomPrint ; iatom ++ )
  {
    fprintf( pGJFoutput , "%d\t% 16.12f\t% 16.12f\t% 16.12f\n" , atomlist[ iatom ].atomlabel , atomlist[ iatom ].cx , atomlist[ iatom ].cy , atomlist[ iatom ].cz );
  
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
    fprintf( pGJFoutput , "\n@GAUSS_EXEDIR:dftba.prm\n\n" );
  }

  fprintf( pGJFoutput , "\n" ) ;


  fclose( pGJFoutput );
  









}







