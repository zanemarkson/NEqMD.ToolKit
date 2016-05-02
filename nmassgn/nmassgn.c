#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "futil.h"
#include "mutil.h"
#include "gutil.h"

#ifndef A2BOHR
#define A2BOHR 0.529177249
#endif



int main( int argc , char * argv[ ] )
{
  FILE * pndx , * pdxdr , * pfreq , * patomlist , * passgn ;
  
  FILE * debug ;
  
  char dxdrname [ 50 ] , ndxname [ 50 ] , atomlistname [ 50 ] , freqname [ 50 ], outassgnname [ 50 ] ;
  
  char ** pcmd ; pcmd = argv ;
  
  int icmd , itmp ;
  
  int ncart , natom , nmode , nfreqprovide;
  
  int icart , iatom , imode , ifreq ;
  
  int exfreq , exatomlist ;
  
  int fakeOrNot = 0 ;

  time_t current_time;

  time( &current_time );

  char now[ 300 ] ;

  strcpy( now , ctime( &current_time ) );

  int lennow = strlen( now ) ;

  *( now + lennow - 1 ) = ' ';

   
  // ========> Recording Command-Line Arguments ...
  
    printf("\n**********************************************************************\n");
      printf("* G_NMASSGN_D : Listing the localization of all vibrational modes.   *\n");
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
  
  strcpy( ndxname , "system.index" );
  
  strcpy( freqname , "vibfrequency.deb" );
  
  strcpy( outassgnname , "system.assgn" );
  
  strcpy( atomlistname , "atom.list" );
  
  // =====> Parsing command line input arguments ...
  
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

	      case 'o' : strcpy( outassgnname , *( ++ pcmd ) ); 
	      
	                 printf("\nCommand-line argument indicates : Output File name : %s ...\n" , outassgnname ); 
	                 
	                 icmd = icmd + 2 ; 
	                 
	                 break ; 
	                 
	      case 'w' : strcpy( freqname , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Input frequency File name : %s ...\n" , freqname ); 
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;

	      case 'n' : strcpy( ndxname , *( ++ pcmd ) );
	      
	                 printf("\nCommand-line argument indicates : Input index File name : %s ...\n" , ndxname ); 
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;

	      case 'l' : strcpy( atomlistname , *( ++ pcmd ) );
	      
	                 if( ( strcmp( atomlistname , "none" ) ) == 0 )
	                 {
	                   printf("\nCommand-line argument indicates : No AtomList needed ... Will not produce fake frequency file ...\n" );
	                   
	                   fakeOrNot = 0 ; 
	                 }
	                 else
	                 {
	                   printf("\nCommand-line argument indicates : Input atom list File name : %s ...\n" , atomlistname ); 
	                   
	                   fakeOrNot = 1 ;
	                 }
	                 
	                 icmd = icmd + 2 ;
	                 
	                 break ;
          
	      case 'h' : printf("\nUsage:  %s [ -f 'input dxdr file name' ] [ -n 'input index file name (in GROMACS .ndx format)' ] [(optional) -o 'output NM assignment file name' ] [ -l atom list file (listing the atomic number of all atoms in system ; \"none\" if no frequency file desired)][ -w 'input vibrational frequency file ; best case length = nmode' ]\n\n" , * argv ); 
	                 
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
  
  if( ( pndx = fopen( ndxname , "r" ) ) == NULL )
  {
    printf("\nUser-defined index file %s does not exist ...\n\n" , ndxname );
    
    exit( 1 );
  
  }

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
  
  if( fakeOrNot == 1 )
  {
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
  
  }



  // =====> Read the index file (*****.ndx) ... 
  
  int countlrb , countrrb ;
  
  int ngroup ;
  
  int ilrb , irrb , igroup ;
  
  char c ;
  
  char ** groupnames ;
  
  char tmpstring[ 100 ] ;
  
  char tmpchar ;
  
  // -----> Finding out how many groups are defined in index file ...
  
  ilrb = 0 ; irrb = 0 ; 
  
  //fsearch( pndx , "begin" );
  
  //fscanf( pndx , "%s" , tmpstring );
  
  //if( strcmp( tmpstring , "]" ) != 0 )
  //{
  //  printf("\nSomething is wrong with the ndx file format at the beginning ... Mission Aborting ... \n") ;
  //  
  //  exit(3);
  //
  //}
  
  //fskip( pndx , 1 ) ;
  
  rewind( pndx ) ;
  
  while( ( c = fgetc( pndx ) ) != EOF )
  {
    if( c == '[' )
    {
      ilrb ++ ;
    }
    else if( c == ']')
    {
      irrb ++ ;
    }
    else
    {
      continue ;
    }
  
  }
  
  if( ilrb == irrb )
  {
    ngroup = ilrb ;
    
    printf("\nThere are %d defined groups in the index file ...\n" , ngroup );
    
  }
  else  
  {
    printf("\nSomething is wrong in the index file ... There are unmatched brackets or broken entries ...\n");
    
    exit( 2 );
    
  }
  
  if( ngroup == 0 )
  {
    printf("\nThis index file is an empty file ... Mission Aborting ...\n");
    
    exit(4);
  }
  

  // -----> Finding out how many atoms are in each group 
  
  int * natom_in_each_group = calloc( ngroup , sizeof( int ) ) ;
  
  int current_atom_index , next_atom_index ;
  
  int info ;
  
  rewind( pndx ) ;
  
  //fsearch( pndx , "begin" );
  
  //fskip( pndx , 1 ) ;
  
  for( igroup = 0 ; igroup < ngroup ; igroup ++ )
  {
    printf("\n--------------------------> GROUP #%d <----------------------------" , igroup + 1 );
    
    //itmp = 0 ;
    
    fsearch( pndx , "]" ) ;
    
    //fskip( pndx , 1 ) ;
    
    //printf("\nHere we are ... : %c " , fgetc( pndx ) );
    
    fscanf( pndx , "%s" , tmpstring );
    
    //printf("\nFor this group , the 1st grabbed tmpstring is %s ... \n" , tmpstring );
    
    //while( strcmp( tmpstring , "[") != 0 )
   
    for( itmp = 0 ; strcmp( tmpstring , "[" ) != 0 && info != EOF ;   )
    {
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
      if( strcmp( tmpstring , "-" ) == 0 )
      {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
        fscanf( pndx , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
        
        itmp = itmp + ( next_atom_index - current_atom_index );
        
        *( natom_in_each_group + igroup ) = itmp ;
        
        current_atom_index = next_atom_index ;
        
        info = fscanf( pndx , "%s" , tmpstring ) ;
      }
      else
      {
        itmp ++ ;
        
        *( natom_in_each_group + igroup ) = itmp ;
        
        current_atom_index = atoi( tmpstring );
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        info = fscanf( pndx , "%s" , tmpstring );
      }
      
      
    
    }
    
    if( info == EOF )
    {
      printf("\nHit the bottom of file ...\n\n") ;
      
      break ;
    }
  
    printf("\nDone with #%d group ... \n" , igroup + 1 );
    

  }

  printf("\n--------------------------> Done checking up groups <----------------------------");


  /* Debugging output for this part ... */
  
  printf("\n== Group Information Summary ==\n");
  
  for( igroup = 0 ; igroup < ngroup ; igroup ++ )
  {
    printf("\nNo.%d Group Has %d atoms ...\n" , igroup + 1 , *( natom_in_each_group + igroup ) );
  
  }

  // -----> Recording the name of each group ...
  
  groupnames = ( char ** ) calloc( ngroup , sizeof( char * ) ) ;
  
  int tmp_groupname_length ;
  
  rewind( pndx ) ;
  
  //fsearch( pndx , "begin" );
  
  //fskip( pndx , 1 );
  
  for( igroup = 0 ; igroup < ngroup ; igroup ++ )
  {
    fsearch( pndx , "[" ) ;
    
    fscanf( pndx , "%s" , tmpstring );
    
    tmp_groupname_length = strlen( tmpstring ) ;
    
    *( groupnames + igroup ) = ( char * ) calloc( tmp_groupname_length + 3 , sizeof( char ) ) ;
    
    strcpy( *( groupnames + igroup ) , tmpstring ) ;
    
    printf("\nName of %d group is %s ...\n" , igroup + 1 , *( groupnames + igroup ) );
  
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
    if( fakeOrNot == 1 )
    {
      for( iatom = 0 ; iatom < natom ; iatom ++ ) 
      {
        *( atomlist + iatom ) = 6 ;
      }
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
  
  
  
  
  // =====> Loading the eigenvectors ( dxdr file ) ... 
  
  
  
  // -----> Based on group information, prepare the recording structures (array) ... 
  
  printf("\n=================================> _ <=====================================\n");
  
  printf("\nStarting from here, we will actually read index file and give the assignment ... \n");
  
  double * assgninfo = calloc( ngroup * nmode , sizeof( double ) ) ;
  
  rewind( pndx ) ;
  
  //fsearch( pndx , "begin" );
  
  //fskip( pndx , 1 ) ;
  
  double tmp_xcomponent , tmp_ycomponent , tmp_zcomponent , tmp_magnitude ;
  
  //printf("\nDebugging ... Debugging ... % 12.8E\t% 12.8E\t% 12.8E\n\n" , *( dxdr + 3 ) , *( dxdr + 4 ) , *( dxdr + 5 ) );


  info = 9 ;
  
  for( igroup = 0 ; igroup < ngroup ; igroup ++ )
  { 
    printf("\n--------------------------> GROUP #%d : %s <----------------------------" , igroup + 1 , *( groupnames + igroup ) );
    
    fsearch( pndx , "]" ) ;
    
    fscanf( pndx , "%s" , tmpstring );
    
    //printf("\nFor this group , the 1st grabbed tmpstring is %s ... \n" , tmpstring );
    
    for( itmp = 0 ; strcmp( tmpstring , "[" ) != 0 && info != EOF ;   )
    {   
      //printf("\nGrabbed tmpstring is %s ... \n" , tmpstring ) ;
      
      if( strcmp( tmpstring , "-" ) == 0 )
      {
        //printf("\nThe '-' situation happened ... before '-' we are at No.%d atom ... \n" , current_atom_index );
        
        fscanf( pndx , "%d" , &next_atom_index );
        
        //printf("\nAnd after '-' we are at No.%d atom ... \n" , next_atom_index );
        
        itmp = itmp + ( next_atom_index - current_atom_index );
        
        for( iatom = current_atom_index ; iatom < next_atom_index ; iatom ++ )
        {  
           //tmp_xcomponent = 0.0000 ; tmp_ycomponent = 0.0000 ; tmp_zcomponent = 0.0000 ; 
      
           //tmp_magnitude = tmp_xcomponent * tmp_xcomponent + tmp_ycomponent * tmp_ycomponent + tmp_zcomponent * tmp_zcomponent ;
           
           printf("\nWorking on No.%d atom ... in No.%d group ...\n" , iatom + 1 , igroup + 1 );
           
           for( imode = 0 ; imode < nmode ; imode ++ )
           {
             tmp_xcomponent = *( dxdr + ( 3 * ( iatom - 0 ) + 0 ) * nmode + imode ) ;
             
             //printf("\nNO.%d mode has % 12.8E on #%d atom X ... \n" , imode + 1 , tmp_xcomponent , iatom + 1 );
             
             tmp_ycomponent = *( dxdr + ( 3 * ( iatom - 0 ) + 1 ) * nmode + imode ) ;
             
             //printf("\nNO.%d mode has % 12.8E on #%d atom Y ... \n" , imode + 1 , tmp_ycomponent , iatom + 1 );
           
             tmp_zcomponent = *( dxdr + ( 3 * ( iatom - 0 ) + 2 ) * nmode + imode ) ;
             
             //printf("\nNO.%d mode has % 12.8E on #%d atom Z ... \n" , imode + 1 , tmp_zcomponent , iatom + 1 );
             
             tmp_magnitude = tmp_xcomponent * tmp_xcomponent + tmp_ycomponent * tmp_ycomponent + tmp_zcomponent * tmp_zcomponent ;
             
             //printf("\nNO.%d mode has % 12.8E on #%d atom ... \n" , imode + 1 , tmp_magnitude , iatom + 1 ) ;
             
             *( assgninfo + igroup * nmode + imode ) = *( assgninfo + igroup * nmode + imode ) + tmp_magnitude ;
             
             printf("\nUp to No.%d atom , No.%d mode has % 12.8E on #%d group ...\n" , iatom + 1 , imode + 1 , *( assgninfo + igroup * nmode + imode ) , igroup + 1 );
           
           }

        }
        
        current_atom_index = next_atom_index ;
        
        info = fscanf( pndx , "%s" , tmpstring ) ;
      }
      else
      {
        itmp ++ ;
        
        current_atom_index = atoi( tmpstring );
        
        //printf("\nNormal situation ... current_atom_index is %d ... \n" , current_atom_index );
        
        printf("\nWorking on No.%d atom ... in No.%d group ...\n" , current_atom_index  , igroup + 1 );
        
        for( imode = 0 ; imode < nmode ; imode ++ )
        {
          tmp_xcomponent = *( dxdr + ( 3 * ( current_atom_index - 1 ) + 0 ) * nmode + imode ) ;
          
          //printf("\nNO.%d mode has % 12.8E on #%d atom X ... \n" , imode + 1 , tmp_xcomponent , current_atom_index );
          
          tmp_ycomponent = *( dxdr + ( 3 * ( current_atom_index - 1 ) + 1 ) * nmode + imode ) ;
          
          //printf("\nNO.%d mode has % 12.8E on #%d atom Y ... \n" , imode + 1 , tmp_ycomponent , current_atom_index );
          
          tmp_zcomponent = *( dxdr + ( 3 * ( current_atom_index - 1 ) + 2 ) * nmode + imode ) ;
          
          //printf("\nNO.%d mode has % 12.8E on #%d atom Z ... \n" , imode + 1 , tmp_zcomponent , current_atom_index );
          
          tmp_magnitude = tmp_xcomponent * tmp_xcomponent + tmp_ycomponent * tmp_ycomponent + tmp_zcomponent * tmp_zcomponent ;
          
          //printf("\nNO.%d mode has % 12.8E on #%d atom ... \n" , imode + 1 , tmp_magnitude , current_atom_index ) ;
          
          *( assgninfo + igroup * nmode + imode ) = *( assgninfo + igroup * nmode + imode ) + tmp_magnitude ;
          
          printf("\nUp to No.%d atom , No.%d mode has % 12.8E on #%d group ...\n" , current_atom_index , imode + 1 , *( assgninfo + igroup * nmode + imode ) , igroup + 1 );
        
        }
        
        info = fscanf( pndx , "%s" , tmpstring );
        
      }
      
    }
    
    if( info == EOF )
    {
      printf("\nHit the bottom of file ...\n\n") ;
      
      break ;
    }
  
    printf("\nDone Loading Information of [ % 5d ] group ... \n" , igroup + 1 );
    
  
  }


  /* Debugging output ... */
  
  debug = fopen( "assgnfull.info" , "wb+" );
  
  fprintf( debug , "Mode#\t" ) ;
  
  for( igroup = 0 ; igroup < ngroup ; igroup ++ )
  {
    fprintf( debug , "% 12s\t" , *( groupnames + igroup ) ) ;
  
  }
  
  fprintf( debug , "vib-Frequencies" );
  
  fprintf( debug , "\n\n" );
  
  for( imode = 0 ; imode < nmode ; imode ++ )
  {
    fprintf( debug , "%d\t" , imode + 1 );
    
    for( igroup = 0 ; igroup < ngroup ; igroup ++ )
    {
      fprintf( debug , "% 12.8E\t" , *( assgninfo + igroup * nmode + imode ) );
    
    }
  
    fprintf( debug , "% 12.8E\n\n" , *( freq + imode ) ) ;
  
  }




  // =====> Writing the output file ( MATLAB loadable pure number array file ) ... 
  
  
  passgn = fopen( outassgnname , "wb+" ) ;
  
  fprintf( passgn , "\n\n" );
  
  for( imode = 0 ; imode < nmode ; imode ++ )
  {
    fprintf( passgn , "%d\t" , imode + 1 );
    
    for( igroup = 0 ; igroup < ngroup ; igroup ++ )
    {
      fprintf( passgn , "% 12.8E\t" , *( assgninfo + igroup * nmode + imode ) );
    
    }
  
    fprintf( passgn , "% 12.8E\n\n" , *( freq + imode ) ) ;
  
  }


  // =====> Writing the output file ( MATLAB loadable pure number array file ) ... 

 FILE * pfakefreq ;
 
 int ibatch = 0 ;
 
 int nbatch = nmode / 3 ;
 
 
 

 if( fakeOrNot == 1 )
 {
   pfakefreq = fopen( "fakeg09.freq" , "wb+" );
   
   if( natom == 1 )
   {
     printf("\nReally? Really?? Really???\n\n");
   
     exit( 1 ) ;
   }
   else if( natom == 2 )
   {
     fprintf( pfakefreq , "                  1\n                    A\nFrequencies -- %10.4f\n" , *( freq + 5 ) ) ;
   
     fprintf( pfakefreq , "Red. masses -- %10.4f\nFrc consts  -- %10.4f\nIR Inten    -- %10.4f\n Atom  AN      X      Y      Z     \n" , 1.00 , 1.00 , 1.00 );
 
     for( iatom = 0 ; iatom < natom ; iatom ++ )
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
     
       for( iatom = 0 ; iatom < natom ; iatom ++ )
       {
         fprintf( pfakefreq , "%5d%4d" , iatom + 1 , *( atomlist + iatom ) ) ;
     
         fprintf( pfakefreq , "   % 8.4f  % 8.4f  %8.4f" , *( dxdr + nmode * ( 3 * iatom + 0 ) + ( 3 * ( ibatch + 0 ) + 0 ) ) , *( dxdr + nmode * ( 3 * iatom + 1 ) + ( 3 * ( ibatch + 0 ) + 0 ) ) , *( dxdr + nmode * ( 3 * iatom + 2 ) + ( 3 * ( ibatch + 0 ) + 0 ) )  ) ;
       
         fprintf( pfakefreq , "      % 8.4f  % 8.4f  %8.4f" , *( dxdr + nmode * ( 3 * iatom + 0 ) + ( 3 * ( ibatch + 0 ) + 1 ) ) , *( dxdr + nmode * ( 3 * iatom + 1 ) + ( 3 * ( ibatch + 0 ) + 1 ) ) , *( dxdr + nmode * ( 3 * iatom + 2 ) + ( 3 * ( ibatch + 0 ) + 1 ) )  ) ;    

         fprintf( pfakefreq , "      % 8.4f  % 8.4f  %8.4f\n" , *( dxdr + nmode * ( 3 * iatom + 0 ) + ( 3 * ( ibatch + 0 ) + 2 ) ) , *( dxdr + nmode * ( 3 * iatom + 1 ) + ( 3 * ( ibatch + 0 ) + 2 ) ) , *( dxdr + nmode * ( 3 * iatom + 2 ) + ( 3 * ( ibatch + 0 ) + 2 ) )  ) ;    
       
       }
   
     }
   

   }

 }









  return( 0 );

} // The End ... 


