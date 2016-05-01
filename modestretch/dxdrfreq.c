#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gutil.h"
#include "mutil.h"
#include "futil.h"



#ifndef A2BOHR
#define A2BOHR 0.529177249
#endif






void dxdrfreq( char massname[] , char hessname[] , char hessunit[] , int natom , int natomselect , double * dxdr , double * freq )
{

FILE * pmass , * phess ; //, *pmo2ao;

int ncart , nmode ;

int natomprovide , ncartselect , ncartprovide , nmodeselect ;

int len_hess_cart, len_tri_hess_cart;

int j,k,m;

int itmp; double dtmp;

char * keyword;


double * hess_cart , * hess_cart_select ;

double * mass_provide , * mass; 

double * vib_dxdr , * vib_freq ;

/*-o-o-o-o-o-o-o-o For debugging output o-o-o-o-o-o-o-o-o-*/

FILE * debug;





/*-o-o-o-o-o-o-o Read input command-line arguments and input files while deciding some parameters ...o-o-o-o-o-o-o-o-o-o-*/

// ====> Mass File ... <==== //

if( ( pmass = fopen( massname ,"r" ) ) == NULL )
{
  printf("\nUser defined mass info file name is %s ...\n\n" , massname );
  
  printf("\nFile containing mass information does not exist...\n");

  exit(1);
}
else
{
  printf("\nUser defined mass info file name is %s ...\n\n" , massname );

}

itmp = flength( pmass ) ;

rewind( pmass );


// =====> Input NAtom and NAtom_Select : Two integer numbers <===== //
//

printf("\nSo there are %d atoms in total in this system ...\n" , natom );

printf("\n%d atoms are selected ...\n" , natomselect );


if( itmp != natom && itmp != natomselect )
{
  printf("\nSomething is wrong with the mass file ... There are %d atomic mass in the file while the total NAtom in this system is %d and %d atoms are selected ... \n" , itmp , natom , natomselect );

  exit(1);
}
else if( itmp == natomselect )
{
  printf("\nOKay ... I see you provided the mass info for just the selected atoms ... \n");
}
else if( itmp == natom )
{
  printf("\nOKay ... I see you provided the mass for all the atom in this system ... \n"); 
}
else if( itmp > natomselect && itmp < natom )
{
  printf("\nOKay ... you provided %d numbers for mass. The total NAtom in system is %d while %d atoms are selected ...\n" , itmp , natom , natomselect );

  printf("\nSo ... I will assume the first %d numbers in your mass file correspond to the selected atoms ...\n" , natomselect );
}
else
{
  printf("\nSomething is wrong with the mass file ... There are %d atomic mass in the file while the total NAtom in this system is %d and %d atoms are selected ... \n" , itmp , natom , natomselect );

  exit(1);

}


ncart = 3 * natom ; 

nmode = ncart - 0 ;

ncartselect = 3 * natomselect ;

nmodeselect = ncartselect - 0 ;

natomprovide = itmp ;


// =====> Hessian File :

double hessconvert ;

if( strcmp( hessunit , "au" ) == 0 )
{
  printf("\nThe unit this Hessian file in is ATOMIC UNIT : Hartree/(Bohr^2) ... \n");
  
  hessconvert = 1.0000 ;
}
else if( strcmp( hessunit , "gmx" ) == 0 )
{
  printf("\nThe unit this Hessian file in is GMX-Unit : kJ/mol/(nm^2) ... \n");
  
  hessconvert = 9.37582E5 ;
}
else
{
  printf("\nWhat did you say about the unit of you frequency again ??? \n");

  exit( 107 );
} 


if( ( phess = fopen( hessname , "r" ) ) == NULL )
{
  printf("\nUser defined Hessian file name is %s ...\n\n" , hessname );
  
  printf("\nCartesian Hessian file does not exist...\n");

  exit(1);
}
else
{
  printf("\nUser defined Hessian file name is %s ...\n\n" , hessname );

}

int gmxhesslength ;

fskip( phess , 1 );

fscanf( phess , "%d" , &gmxhesslength ) ;


fskip( phess , 1 );

len_hess_cart = flength( phess );

if( len_hess_cart != gmxhesslength * gmxhesslength )
{
  printf("\nSomething is wrong with reading the hessian file ... Hessian file says itself the dimension is %d * %d but the actuall length of file is %d ... \n" , gmxhesslength , gmxhesslength , len_hess_cart );

  printf("\nI am gonna try to process the whole thing anyway ... But mark this WARNING here please ... \n");

}
else
{
  printf("\nAccording to user-provided Hessian file, the dimension of provided Hessian is %d * %d ... \n" , gmxhesslength , gmxhesslength );
}


if( len_hess_cart == ncart * ncart )
{
  printf("\nOKay ... I see you provided the Hessian for all %d atom in system ... \n" , natom );
}
else if( len_hess_cart == ncartselect * ncartselect )
{
  printf("\nOKay ... I see you only provided the Hessian for the %d selected atoms ... \n" , natomselect );
}
else if( len_hess_cart < ncart * ncart && len_hess_cart > ncartselect * ncartselect )
{
  printf("\nThere are %d numbers in the Hessian file which is less than the square of total %d Cartesian coordinates in system but more than that of %d Cartesian coordinates of selected atoms ... \n" , len_hess_cart , ncart , ncartselect );

  printf("\nSo I am taking the first %d number ( %d * %d ) as the Hessian for selected atoms ...\n" , ncartselect * ncartselect , ncartselect , ncartselect );

}
else
{
  printf("\nSomething is wrong with the Hessian file ... There are %d numbers in the Hessian file ... \n" , len_hess_cart );

  exit(1);

}

rewind( phess );

/*-o-o-o-o-o-o-o Allocating space ... o-o-o-o-o-o-o-o-*/

ncartprovide = sqrt( len_hess_cart ); printf("\nAlthough ... the Hessian file you provided is a %d * %d square matrix ...\n" , ncartprovide , ncartprovide );

hess_cart = calloc( len_hess_cart , sizeof(double));  dzeros( len_hess_cart , 1, hess_cart );

hess_cart_select = calloc( ncartselect * ncartselect , sizeof(double) ) ; dzeros( ncartselect * ncartselect , 1 , hess_cart_select );

//tri_hess_cart = calloc( ncart*(ncart+1)/2 , sizeof(double)); dzeros(ncart*(ncart+1)/2, 1, tri_hess_cart);

vib_freq = calloc( nmodeselect , sizeof(double));  dzeros( nmodeselect , 1 , vib_freq );

vib_dxdr = calloc( nmodeselect * nmodeselect , sizeof(double));  dzeros( nmodeselect * nmodeselect , 1 , vib_dxdr );

mass_provide = calloc( natomprovide , sizeof( double ) );

dzeros( natomprovide , 1 , mass_provide );

mass = calloc( natomselect , sizeof( double ) );

dzeros( natomselect , 1 , mass );



/*-o-o-o-o-o-o-o Importing data ... o-o-o-o-o-o-o-o-*/


/* Importing Hessian Matrix (in Cartesian coordinate) ... */

//keyword="Full";

//fsearch( phess , keyword ) ;

rewind( phess );

fskip( phess , 2 );

fload( phess , hess_cart );

for( j = 0 ; j < len_hess_cart ; j ++ )
{
  *( hess_cart + j ) = *( hess_cart + j )  ; // For G09 Hessian ( already in atomic units , Hartree/Bohr^2 )
  //*( hess_cart + j ) = *( hess_cart + j ) / 9.37582E5 ; // For GROMACS Hessian ( kj/mol/nm^2 )
  *( hess_cart + j ) = *( hess_cart + j )  / hessconvert ; // For All ... Controlled by hessconvert ;
}


/* Selecting the desired part of Hessian matrix */

for( j = 0 ; j < ncartselect ; j ++ )
{
  for( k = 0 ; k < ncartselect ; k ++ )
  {
    *( hess_cart_select + j * ncartselect + k ) = *( hess_cart + j * ncartprovide + k );
  }

}


/* Change triangular Hessian(Cart) into symmetric matrix  */

//dtri2sym( 'L', ncart, tri_hess_cart, hess_cart);



debug = fopen("hess_cart.deb", "wb+");

doutput(debug, ncart, ncart, hess_cart);

fclose(debug);


debug = fopen("hess_cart_select.deb", "wb+");

doutput( debug , ncartselect , ncartselect , hess_cart_select );

fclose(debug);





/* Importing the mass of each atom */

rewind( pmass );


fload( pmass , mass_provide );

debug = fopen( "mass_provide.deb", "wb+" );

doutput( debug , natomselect , 1 , mass_provide );

fclose(debug);



for( j = 0 ; j < natomselect ; j ++ ) *( mass + j ) = *( mass_provide + j );

debug = fopen( "mass.deb", "wb+" );

doutput( debug , natomselect , 1 , mass );

fclose(debug);



for( j = 0 ; j < natomselect ; j ++ )   *( mass + j ) = ( *( mass + j ) ) * AMU2AU;

debug = fopen("mass_au.deb", "wb+");

doutput( debug , natomselect , 1 , mass );

fclose(debug);


/* Performing mass-weighting for the Force constant matrix  */





masswt( natomselect , mass , hess_cart_select , hess_cart_select );



debug = fopen("hess_masswt_select.deb", "wb+");

doutput( debug , ncartselect , ncartselect , hess_cart_select );

fclose(debug);





/* Performing matrix-diagonalization for Hess_Cart*/

dsyev_f2c( ncartselect , hess_cart_select , vib_dxdr, vib_freq  );

dtranspose( ncartselect , vib_dxdr , vib_dxdr );

mark(2);


if( natomselect == 1 )
{
  printf("\nSeriously? Only one atom ? NO WAYYYYYY ... \n");

  exit( 90 );
}
else if( natomselect == 2  )
{
  *( vib_freq + 5 ) = sqrt( *( vib_freq + 5 ) )  ;
}

else
{
  for( j = 0 ; j < nmodeselect ; j ++ )
  {
    if( *( vib_freq + j ) >= 0.000 )
	    
	    *( vib_freq + j ) = sqrt( *( vib_freq + j ) ) ;
    else
	    *( vib_freq + j ) = -1.000 * sqrt( -1.0000 * ( *( vib_freq + j ) ) ) ; // See, we are in ATOMIC UNITS ... ;
  
  }

}



// =====> Output ... to calling routine and file ... 

for( j = 0 ; j < nmodeselect * nmodeselect ; j ++ )  *( dxdr + j ) = *( vib_dxdr + j );

for( j = 0 ; j < nmodeselect ; j ++ )  *( freq + j ) = *( vib_freq + j );







debug = fopen("dxdr.deb", "wb+");

doutput(debug, ncartselect, ncartselect, vib_dxdr);

fclose(debug);


debug = fopen("vibfrequency.deb", "wb+");

doutput(debug, nmodeselect, 1, vib_freq);

fclose(debug);




/* So ... There is an additional 1/sqrt(mass) thing to do ... AHHHHHH */


/*
massdivide( 'R', natom, mass, dxdr, dxdr);




debug = fopen("dxdr_massdvd.deb", "wb+");

doutput(debug, ncart, ncart, dxdr);

fclose(debug);

*/




/* Do the DXDR*supervector(Hessian_Cart) thing...  */






/* Output the desired dfdr matrix   */  



/* I assume it's done... AHHHHHHHH*/


printf("\n\nI assume it's Allllllll Done   ...\n\n");



















/*

free(nbomo);

free(hess_cart);

free(tri_hess_cart);

free(vib_freq);

free(dxdr);

free(dfdx_array);

free(tri_dfdx_array);

free(poltri_dfdx_array);

free(dfdr_array);

free(dfdr_nbo);

free(mass);

*/









//return(1);

}






