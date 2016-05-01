#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "futil.h"
#include "mutil.h"
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



int main( int argc , char * argv[ ] )
{
  FILE * pdxdr , * pfreq , * patomlist , * pfakefreq ;
  
  FILE * debug ;
  
  char dxdrname [ 50 ] , atomlistname [ 50 ] , freqname [ 50 ], outfreqname [ 50 ] ;
  
  char ** pcmd ; pcmd = argv ;
  
  int icmd , itmp ;
  
  int natom , ncart , nmode , nfreqprovide;
  
  int natomSelect , ncartSelect , nmodeSelect ;
  
  int icart , iatom , imode , ifreq ;
  
  int exfreq , exatomlist , exselect ;
  
  
  char buffer[ MAXCHARINLINE ] ;

  int lenBuff ;

  char cache[ MAXCHARINLINE ] ;

  char stmp[ 150 ] , stmp2[ 150 ] , tmpString[ 150 ];
  
  

  time_t current_time;

  time( &current_time );

  char now[ 300 ] ;

  strcpy( now , ctime( &current_time ) );

  int lennow = strlen( now ) ;

  *( now + lennow - 1 ) = ' ';

   
  // ========> Recording Command-Line Arguments ...
  
    printf("\n**********************************************************************\n");
      printf("* G_FAKEFREQ_D : Organize (DX/DR) into G09 Format for visualization. *\n");
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

 

  // =====> Setting up default file names ... 
  
  
  strcpy( dxdrname , "dxdr.deb" );
  
  strcpy( freqname , "vibfrequency.deb" );
  
  strcpy( outfreqname , "fakeg09.freq" );
  
  strcpy( atomlistname , "atom.list" );
  
  // =====> Parsing command line input arguments ...
  
  exselect = 20 ;
  
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
	      
	      case 'f' : strcpy( dxdrname , *( ++ pcmd ) ) ; 
			 
			         printf("\nCommand-line argument indicates : Input dxdr File name : %s ...\n" , dxdrname ); 
	      
	                 icmd = icmd + 2 ; 
	                 
	                 break ;

	      case 'o' : strcpy( outfreqname , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Output File name : %s ...\n" , outfreqname ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ; 
	                 
	      case 'w' : strcpy( freqname , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Input frequency File name : %s ...\n" , freqname ); 
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;

	      case 'l' : strcpy( atomlistname , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Input atom list File name : %s ...\n" , atomlistname ); 
	                 
	                 icmd = icmd + 2 ;
	                 
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
          

          
	      case 'h' : printf("\nUsage:  %s [ -f 'input dxdr file name' ] [(optional) -o 'output NM assignment file name' ] [ -l atom list file (listing the atomic number of all atoms in system)] [ -s # of atom selected to output ] [ -w 'input vibrational frequency file' ]\n\n" , * argv ); 
	                 
	                 printf("\n\n==> NOTE : In order to be fool-proof , if there is mis-match on the NAtom between dxdr file and any other file , this code will set NAtom to be the number from dxdr automatically ...\n\n");
	                 
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
  

  // =====> Open the files ...


  if( ( pdxdr = fopen( dxdrname , "r" ) ) == NULL )
  {
    printf("\nUser-defined dxdr file %s does not exist ...\n\n" , dxdrname );
    
    exit( 1 );
  
  }

  double * freq ;

  if( ( pfreq = fopen( freqname , "r" ) ) == NULL )
  {
    printf("\nUser-defined frequency file %s does not exist ...\n\n" , freqname );
    
    printf("\nBut it's OK , we will procede without frequency information ... \n");
    
    exfreq = 0 ;
  
  }
  else
  {
    nfreqprovide = flength( pfreq );
    
    rewind( pfreq ) ;
    
    freq = ( double * ) calloc( nfreqprovide , sizeof( double ) ) ;
    
    fload( pfreq , freq );
    
    rewind( pfreq );
    
    exfreq = 1 ;

  }


  int * atomlist , natomlistprovide ;
  

  if( ( patomlist = fopen( atomlistname , "r" ) ) == NULL )
  {
    printf("\nUser-defined atom list file %s does not exist ...\n\n" , atomlistname );
    
    printf("\nMeaning ... All of your atoms will be carbon ...\n");
    
    exatomlist = 0 ;
  
  }
  else
  {
    natomlistprovide = flength( patomlist );
    
    rewind( patomlist );
    
    atomlist = ( int * ) calloc( natomlistprovide , sizeof( int ) ) ;
    
    int_fload( patomlist , atomlist );
    
    exatomlist = 1 ;
  }




  
  // =====> Loading the eigenvectors ( dxdr file ) ... 
  
  printf("\n--------------------------> Begin reading dxdr file ...  <----------------------------" );
  
  int len_dxdr = flength( pdxdr ) ;
  
  printf("\nThere are %d numbers in eigenvectors file ... \n" , len_dxdr );
  
  ncart = sqrt( len_dxdr );
  
  natom = ncart / 3 ;
  
  double * dxdr = calloc( len_dxdr , sizeof( double ) ) ;
  
  rewind( pdxdr );
  
  fload( pdxdr , dxdr ) ;
  
  rewind( pdxdr ) ;  
  
  printf("\nDone with loading dxdr file ... \n");
  
  switch( natom )
  {
    case 1 : printf("\nSeriously? Only ONE atom? ...\n") ; exit( 9 ) ; break ;
    
    case 2 : nmode = 6 ; break ;
    
    default : nmode = ncart - 0 ; break ;
  
  }

  if( exfreq == 0 )
  {
    freq = ( double * ) calloc( nmode , sizeof( double ) ) ;
    
    dzeros( nmode , 1 , freq ) ;
  }
  else
  {
    if( nfreqprovide > nmode )
    {
      printf("\nOkay ... The number of vibrational frequency you provided is more than the normal mode we have from eigenvectors ...\n");
      
      printf("\nSo I will take only the first #_of_eigenvectors frequencies from the file you provided ...\n");
    }
    else if( nfreqprovide < nmode )
    {
      printf("\nOkay ... The number of vibrational frequency you provided is less than the normal mode we have from eigenvectors ...\n");
      
      printf("\nThe rest will be padded with zeros ... Sorry, I really don't want to do the diagonalization ...\n");
      
      free( freq ) ;
      
      freq = ( double * ) calloc( nmode , sizeof( double ) ) ;
      
      rewind( pfreq );
      
      fload( pfreq , freq ) ;
      
      dzeros( nmode - nfreqprovide , 1 , freq + nfreqprovide ) ;

    }
  
  
  }
  
  
  
  if( exatomlist == 0 )
  {
    for( iatom = 0 ; iatom < natom ; iatom ++ ) 
    {
      *( atomlist + iatom ) = 6 ;
    }
    
  }
  else
  {
    if( natomlistprovide > natom )
    {
      printf("\nOkay ... The number of atom you provided is more than the NAtom we have from eigenvectors ...\n");
      
      printf("\nSo I will take only the first #_of_atoms from the file you provided ...\n");
    
    }
    else if( natomlistprovide < natom )
    {
      printf("\nOkay ... The number of atom you provided is less than the NAtom we have from eigenvectors ...\n");
      
      printf("\nSo all the rest will be automatically set as Carbon ... \n");
      
      free( atomlist );
      
      atomlist = ( int * ) calloc( natom , sizeof( int ) ) ;
      
      rewind( patomlist );
      
      int_fload( patomlist , atomlist );
      
      for( iatom = natom - natomlistprovide ; iatom < natom ; iatom ++ )
      {
        *( atomlist + iatom ) = 6 ;
      }
      
      
    
    }
  
  }
  
  if( exselect == 20 )
  {
    printf("\nBy default , all %d atoms are selected as solute ... \n\n" , natom ) ;

    natomSelect = natom ;

  }
  else if( exselect == 22 )
  {
    printf("\nPer user's request , all %d atoms are selected as solute ... \n\n" , natom ) ;

    natomSelect = natom ;
  }
  else if( exselect == 21 )
  {
    if( natomSelect > natom )
    {
      printf("\nOops ... you requested %d atoms as solute but we only have %d atoms in G09 .inp file ... By default , all atoms are selected as solute ... \n\n" , natomSelect , natom ) ;

      natomSelect = natom ;

    }
    else if( natomSelect <= natom && natomSelect > 0 )
    {
      printf("\nPer user's request , the first %d atoms are selected as solute ... \n\n" , natomSelect ) ;
    }
    else if( natomSelect <= 0 )
    {
      printf("\n# of solute atom must be a positive integer!!!\n\n") ;

      exit( 52 ) ;
    }


  }
  
  
  ncartSelect = natomSelect * 3 ;
  
  switch( natomSelect )
  {
    case 1 : printf("\nSeriously? Only ONE atom selected ? ...\n") ; exit( 9 ) ; break ;
    
    case 2 : nmodeSelect = 6 ; break ;
    
    default : nmodeSelect = ncartSelect - 0 ; break ;
  
  }

  

  // =====> Writing the output file ( MATLAB loadable pure number array file ) ... 
  
 printf("\n--------------------------> Begin organizing dxdr file into G09 format...  <----------------------------\n\n\n" );

 pfakefreq = fopen( outfreqname , "wb+" );
 
 int ibatch = 0 ;  
 
 int nbatch = nmodeSelect / 3 ;
 
 
 
 if( natomSelect == 1 )
 {
   printf("\nReally? Really?? Really???\n\n");
   
   exit( 1 ) ;
 }
 else if( natomSelect == 2 )
 {
   fprintf( pfakefreq , "                  1\n                    A\nFrequencies -- %10.4f\n" , *( freq + 5 ) ) ;
   
   fprintf( pfakefreq , "Red. masses -- %10.4f\nFrc consts  -- %10.4f\nIR Inten    -- %10.4f\n Atom  AN      X      Y      Z     \n" , 1.00 , 1.00 , 1.00 );
 
   for( iatom = 0 ; iatom < natomSelect ; iatom ++ )
   {
     fprintf( pfakefreq , "%5d%4d  % 8.4f  % 8.4f  %8.4f\n" , iatom + 1 , *( atomlist + iatom ) , *( dxdr + 3 * iatom + 0 ) , *( dxdr + 3 * iatom + 1 ) , *( dxdr + 3 * iatom + 2 ) );
   
   }
 
 }
 else
 {
   for( ibatch = 2 ; ibatch < nbatch ; ibatch ++ )
   {
     fprintf( pfakefreq , "                      %5d                             %5d                             %5d\n" , 3 * ( ibatch - 2 ) + 1 , 3 * ( ibatch - 2 ) + 2 , 3 * ( ibatch - 2 ) + 3 );
     
     fprintf( pfakefreq , "                      %5c                             %5c                             %5c\n" , 'A' , 'A', 'A' );
     
     fprintf( pfakefreq , "%-12s --  %10.4f                        %10.4f                        %10.4f\n" , " Frequencies" , *( freq + 3 * ibatch + 0 ) , *( freq + 3 * ibatch + 1 ) , *( freq + 3 * ibatch + 2 ) );
     
     fprintf( pfakefreq , "%-12s --  %10.4f                        %10.4f                        %10.4f\n" , " Red. masses" , 1.00 , 1.00 , 1.00 );
     
     fprintf( pfakefreq , "%-12s --  %10.4f                        %10.4f                        %10.4f\n" , " Frc consts" , 10.00 , 10.00 , 10.00 );
     
     fprintf( pfakefreq , "%-12s --  %10.4f                        %10.4f                        %10.4f\n" , " IR Inten" , 10.00 , 10.00 , 10.00 );
     
     fprintf( pfakefreq , "   Atom  AN     X         Y         Z             X         Y         Z             X         Y         Z\n");
     
     for( iatom = 0 ; iatom < natomSelect ; iatom ++ )
     {
       fprintf( pfakefreq , "%5d%4d" , iatom + 1 , *( atomlist + iatom ) ) ;
     
       fprintf( pfakefreq , "   % 8.4f  % 8.4f  %8.4f" , *( dxdr + nmode * ( 3 * iatom + 0 ) + ( 3 * ( ibatch + 0 ) + 0 ) ) , *( dxdr + nmode * ( 3 * iatom + 1 ) + ( 3 * ( ibatch + 0 ) + 0 ) ) , *( dxdr + nmode * ( 3 * iatom + 2 ) + ( 3 * ( ibatch + 0 ) + 0 ) )  ) ;
       
       fprintf( pfakefreq , "      % 8.4f  % 8.4f  %8.4f" , *( dxdr + nmode * ( 3 * iatom + 0 ) + ( 3 * ( ibatch + 0 ) + 1 ) ) , *( dxdr + nmode * ( 3 * iatom + 1 ) + ( 3 * ( ibatch + 0 ) + 1 ) ) , *( dxdr + nmode * ( 3 * iatom + 2 ) + ( 3 * ( ibatch + 0 ) + 1 ) )  ) ;    

       fprintf( pfakefreq , "      % 8.4f  % 8.4f  %8.4f\n" , *( dxdr + nmode * ( 3 * iatom + 0 ) + ( 3 * ( ibatch + 0 ) + 2 ) ) , *( dxdr + nmode * ( 3 * iatom + 1 ) + ( 3 * ( ibatch + 0 ) + 2 ) ) , *( dxdr + nmode * ( 3 * iatom + 2 ) + ( 3 * ( ibatch + 0 ) + 2 ) )  ) ;    
       
     }
   
   }

   
 }

  printf("\n==============> ALL DONE !!! ^-^ !!! <===============\n\n\n\n\n\n");


  return( 0 );

} // The End ... 


