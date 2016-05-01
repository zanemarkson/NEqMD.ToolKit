#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "mutil.h"
#include "futil.h"


//-----------------  Macros  -------------------//

#ifndef PI
#define PI 3.1415926535
#endif








//------------------------------------------------------------//


/* * * * * * * * * * * * * * * * * * * * * * */
/* May. 9th, 2013, Zheng Ma                  */
/* This routine is specificly written for    */
/* the calculation of :                      */
/*       lambda = dF / dR                    */
/* to convert the derivative matrix of Fock  */
/* with respect to cartisian coordinate to   */
/* that with respect to normal mode, i.e.    */
/*   dF / dX ---->  dF / dR                  */ 
/*                                           */
/* * * * * * * * * * * * * * * * * * * * * * */





void dx2dr( double * u, int ulda, double * h, int hlda, double * hnorm)
{
  //int ic, ir, k;

  int iblock, inorm, icart, k;


  double * temp = malloc(hlda*hlda*sizeof(double));
  
  //double * collect = malloc(ulda*hlda*hlda*sizeof(double));

  


  for(inorm=0;inorm<ulda;inorm++)   // icart=iblock;
  {
    //icart=iblock;
    
    for(icart =0;icart<ulda;icart++)
    {
       
      printf("\nicart=%d,  inorm=%d,  num = % 12.10E  ", icart, inorm, *(u+icart*ulda+inorm));

      printf("\nBlock view, hlda = %d: \n", hlda);

      dmatprint(hlda, hlda, h+icart*hlda*hlda );


    //  dnummat( *(u+iblock*ulda+inorm), h, hlda*ulda, hlda, inorm*hlda, 0, hlda, hlda, temp );
      dnummat( *(u+icart*ulda+inorm), h+icart*hlda*hlda, hlda*ulda, hlda, hlda, hlda, temp );
      
      
      
      printf("\nOnly temp ... \n");
      
      dmatprint( hlda, hlda, temp);
//void dmatplus( int rdim, int cdim, double * left, double * right, int leftlda, int rightlda, double * res);
      dmatplus( hlda, hlda, hnorm+inorm*hlda*hlda, temp, hlda, hlda, hnorm+inorm*hlda*hlda);
      
      printf("\nCollective ... \n");
      
      dmatprint(hlda*ulda, hlda, hnorm);
      
      for(k=0; k<hlda*hlda; k++) *(temp+k)=0.00000;

    }
   
    


  }


}





//------------------------------------------------------------//

void massdivide( char direction, int na, double * mass, double * undivided, double * divided)
{
  int j,k;

  int ncart = 3*na;


  double * massexp = malloc( ncart* sizeof(double));
  
  double * buffer = malloc( ncart*ncart* sizeof(double));


  for(j=0;j<na;j++)
  {
    *(massexp+3*j+0)=*(mass+j);
    *(massexp+3*j+1)=*(mass+j);
    *(massexp+3*j+2)=*(mass+j);
  }


 if( (direction == 'C') || (direction == 'c') )
 {
   for(k=0;k<ncart;k++)
   {
     for(j=0;j<ncart;j++)
     {
       *(buffer+j*ncart+k)= (    (*(undivided+j*ncart+k))/sqrt(*(massexp+k))    ); 
     
     }
   
   }
 }  
 else if ( (direction == 'R') || (direction =='r') )
 {
   for(j=0;j<ncart;j++)                                                     
   {                                                                        
     for(k=0;k<ncart;k++)                                                   
     {                                                                      
       *(buffer+j*ncart+k)= (    (*(undivided+j*ncart+k)) / sqrt(*(massexp+j))    ); 
                                                                         
     }                                                                      
                                                                         
    }                                                                    
   
 }
 else
 {
   printf("\nPlease specify which direction the mass division should be made.\n");
   exit(1);
 }


  for(j=0; j<ncart*ncart; j++) 
          *(divided+j)=*(buffer+j);

 
 
 
 
 
 
 





}





//------------------------------------------------------------//


void mo2nbo(int natom, int nbasis, double * u, double * pre, double * post)
{
  int j,k;
  
  int imode;
  
  int bsd = nbasis;
  
  double done = 1.00000;  double dzero = 0.0000;
  
  
  
  
  
  
  double * basetmp = malloc( bsd*bsd*sizeof(double) ); 
  
  dzeros(bsd, bsd, basetmp);
  
  double * inter1 = malloc( bsd*bsd*sizeof(double) );
  
  dzeros(bsd, bsd, inter1);
  
  double * inter2 = malloc( bsd*bsd*sizeof(double) );
  
  dzeros(bsd, bsd, inter2);
  

  
  
  for( imode = 0; imode < 3*natom; imode++)
  {
    
  //dmatassgn(int dimrow, int dimcol,  double  *porigin, double  *ptarget)
    dmatassgn( bsd, bsd, pre+bsd*bsd*imode, basetmp);
    
//SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
  
    dgemm_("T","T", &bsd, &bsd, &bsd, &done, u, &bsd,
           basetmp, &bsd, &dzero, inter1, &bsd);
    
    dgemm_("N","N", &bsd, &bsd, &bsd, &done, inter1, &bsd,
           u, &bsd, &dzero, inter2, &bsd);
          
    dtranspose(bsd, inter2, inter2);

    dmatassgn( bsd, bsd, inter2, post+bsd*bsd*imode);
    
    dzeros(bsd, bsd, basetmp);
    
    dzeros(bsd, bsd, inter1);
    
    dzeros(bsd, bsd, inter2);
  

  }













free(basetmp);  

free(inter1);

free(inter2);


}




//------------------------------------------------------------//



void masswt( int na, double * mass, double * unwtmat, double * wtmat)
{

  int j,k;

  int ncart=3*na;

  double * massexp = malloc( ncart* sizeof(double));

  double * buffer = malloc( ncart*ncart* sizeof(double));


  for(j=0;j<na;j++)
  {
    *(massexp+3*j+0)=*(mass+j);
    *(massexp+3*j+1)=*(mass+j);
    *(massexp+3*j+2)=*(mass+j);
  }


  for(j=0;j<ncart;j++)
	  for(k=0;k<ncart;k++)
		  *(buffer+j*ncart+k)=*(unwtmat+j*ncart+k)/sqrt(*(massexp+j))/sqrt(*(massexp+k));


  for(j=0;j<ncart*ncart;j++)
		  *(wtmat+j)=*(buffer+j);



}





//------------------------------------------------------------//


void pick( FILE * p_of_info_file, char kw[], int nskip, double * database )
{
  
  
  
  int j,length;
  
  int idskip, idstep;
  
  char * ctmp=malloc(99*sizeof(char));
  
  char c;
  
  double * buffer;
  
  
  
  

  //FILE * tmp1=fopen("import_info.tmp","wb+");
  
  FILE * tmp2=fopen("modified_info.tmp","wb+");
  
  fsearch( p_of_info_file, kw);
  
  //fscanf(p_of_info_file, "%s", ctmp); printf("\n%s\n\n", ctmp);
  
  idskip=0; 
  
  c=fgetc(p_of_info_file);
 
  //printf("\nidskip=%d,   idstep=%d\n, c=  %c \n", idskip, idstep, c); 
  
  fskip(p_of_info_file,  nskip);

  //fcopy(p_of_info_file, tmp1);
  
  fsub(p_of_info_file, 'D', 'E', tmp2);
  
  rewind(tmp2);
    
  length=flength(tmp2);

  printf("\n\nThere are %d numbers in this file \n\n", length );
  

  rewind(tmp2);
  
  buffer=malloc(length*sizeof(double));
  
  fload(tmp2, buffer);
  
  
  
  for(j=0;j<length;j++)
  {  
     *(database+j) = *(buffer+j) ;
      
     //printf("\nThe No. %d number is %16.13f\n\n", j, *(buffer+j));
  
  }
  
  
  
  fdel( "modified_info.tmp"  );

  
  

}


//------------------------------------------------------------//



/* * * * * * * * * * * * * * * * * * * * * */
/* May. 24th, 2013, Zheng Ma               */
/* Added function: Read NBasis and return  */
/* as the variable "ires.                  */
/*                                         */
/* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * */
/* May. 24th, 2013, Zheng Ma               */
/* Changed Feature: Now any integer return */
/* value will be placed in *ires, double   */
/* precision value in *dres.               */
/* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * */
/* Jun. 4th, 2013, Zheng Ma                */
/* Added feature for grabbing reduced mass */
/* from Gaussian output file (Freq job).   */
/* Being tested currently ...              */
/* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * */
/* Jun. 6th, 2013, Zheng Ma                */
/* Added feature for grabbing reduced mass */
/* from Gaussian output file (Freq job).   */
/* Testing Done... Ready to work           */
/* * * * * * * * * * * * * * * * * * * * * */



void readgauss( char * flag, FILE * gaussoutput, int * ires, double * dres )
{
  int j,k;
  
  int number_of_atom, number_of_mode;

  int iatom, imode;
  
  char c;
  
  double tmp;
  
  
  
  
  rewind(gaussoutput);
  
  fsearch( gaussoutput, "NAtoms=" );
  
  /* 
  for(j=0;j<6;j++) 
  {
    c=fgetc(gaussoutput);

    printf("\nNext character is %c ... \n",c);
  }
  */



  fscanf(gaussoutput, "%lf", &tmp);

  //printf("\n picked number is %13.10f ...\n\n", tmp);

  rewind(gaussoutput);

  number_of_atom =(int)tmp;
  
  switch (number_of_atom)
  {
  case 1 : number_of_mode = 0; break;

  case 2 : number_of_mode = 1; break;

  default : number_of_mode = 3*number_of_atom - 6; break;
  }


  
  //printf("\n picked number is %d ...\n\n", *number_of_atom);



  if ( strcmp( flag, "mass" ) == 0 )
  {
    //printf("\nTASK : Grabbing the mass list \n");
	  
    fsearch(gaussoutput, "AtmWgt=");

    //c=fgetc(gaussoutput);
    //
    
    iatom=0;

    for( iatom=0; iatom < number_of_atom; iatom++ )
    {
      c=fgetc( gaussoutput );
     
      if( c != '\n')
      {
        fscanf( gaussoutput, "%lf", &tmp );
   
        //printf("\nNo. %d atom, mass is %13.10f ...\n", iatom, tmp);

        *(dres+iatom)=tmp;

      }
      else
      {
        fsearch(gaussoutput, "AtmWgt=");

        c=fgetc(gaussoutput);

        fscanf( gaussoutput, "%lf", &tmp );

	    //printf("\nNo. %d atom, mass is %13.10f ...\n", iatom, tmp);

        *(dres+iatom) = tmp;

      }
	     
	     
    } 
  
  
  }
  
  
  
  else if ( strcmp( flag, "natom" ) == 0 )
  {

    *ires = number_of_atom;


  }
  
  
  
  else if ( strcmp( flag, "nbasis" ) == 0 )
  {

   rewind(gaussoutput);
   
   fsearch( gaussoutput, "NBasis=" );
   
   fscanf(gaussoutput, "%lf", &tmp);
   
   *ires = (int)tmp;

  }
  
  
  
  else if ( strcmp( flag, "reducedmass" ) == 0 )
  {
    //printf("\nTASK : Grabbing the reduced-mass list \n");
	  
    char * word1=malloc(99*sizeof(char));

    char * word2=malloc(99*sizeof(char));
    
    fsearch(gaussoutput, "Red.");
    
    fscanf(gaussoutput, "%s", word1); //printf("\n --- %s ---\n", word1);
    
    fscanf(gaussoutput, "%s", word2); //printf("\n --- %s ---\n", word2);


    while( ( strcmp( word1, "masses" ) != 0 ) ||  ( strcmp( word2, "--" ) != 0 )   )
    {
      fsearch(gaussoutput, "Red.");
    
      fscanf(gaussoutput, "%s", word1); //printf("\n --- %s ---\n", word1);
    
      fscanf(gaussoutput, "%s", word2); //printf("\n --- %s ---\n", word2);

    }
    

    for( imode=0; imode < number_of_mode; imode++ )
    {
      c=fgetc( gaussoutput );
      
      if( c != '\n')
      {
        fscanf( gaussoutput, "%lf", &tmp );

        *(dres+imode)=tmp;
        
        //printf("\nNo. %d normal mode, reduced-mass is % 10.8E ...\n", imode, *(dres+imode));
        

      }
      else
      {
        fsearch(gaussoutput, "Red.");
    
        fscanf(gaussoutput, "%s", word1); //printf("\n --- %s ---\n", word1);
    
        fscanf(gaussoutput, "%s", word2); //printf("\n --- %s ---\n", word2);


        while( ( strcmp( word1, "masses" ) != 0 ) ||  ( strcmp( word2, "--" ) != 0 )   )
        {
          fsearch(gaussoutput, "Red.");
    
          fscanf(gaussoutput, "%s", word1); //printf("\n --- %s ---\n", word1);
    
          fscanf(gaussoutput, "%s", word2); //printf("\n --- %s ---\n", word2);

        }
    


        c=fgetc(gaussoutput);

        fscanf( gaussoutput, "%lf", &tmp );

        *(dres+imode) = tmp;
        
        //printf("\nNo. %d normal mode, reduced-mass is % 10.8E ...\n", imode, *(dres+imode));

      }
	     
	     
    } 
  
  
  }
  
  
  
  else if ( strcmp( flag, "optgeom" ) == 0 )
  {
    rewind(gaussoutput);
    
    fsearch( gaussoutput, "Stationary"  );
    
    fsearch( gaussoutput, "Input"  );
    
    fskip( gaussoutput, 5 );
    
    double * cartcoord = malloc( 3*number_of_atom * sizeof(double) );
    
    
    for(iatom=0; iatom<number_of_atom; iatom++)
    {
      fscanf( gaussoutput, "%lf", &tmp); //printf("\nThis number is %10.8f\n", tmp);
      fscanf( gaussoutput, "%lf", &tmp); //printf("\nThis number is %10.8f\n", tmp);
      fscanf( gaussoutput, "%lf", &tmp); //printf("\nThis number is %10.8f\n", tmp);
      
      fscanf( gaussoutput, "%lf", cartcoord+3*iatom+0 ); 
      //printf("\nX_%d %10.8f\n", iatom, *(cartcoord+3*iatom+1) );
      
      fscanf( gaussoutput, "%lf", cartcoord+3*iatom+1 );
      //printf("\nY_%d %10.8f\n", iatom, *(cartcoord+3*iatom+2) );
      
      fscanf( gaussoutput, "%lf", cartcoord+3*iatom+2 );
      //printf("\nY_%d %10.8f\n", iatom, *(cartcoord+3*iatom+3) );
    
    }
    
    for(j=0; j<3*number_of_atom; j++)
    {
      //printf("\nNo. %d picked number : % 10.8f\n", j, *(cartcoord+j) );
      
      *(dres+j) = *(cartcoord+j);
    
    }
  
  
  
  }

  else
  {
    printf("\nPlease specify the flag of the info you need ...\n\n");
  }

  








}








//------------------------------------------------------------//

double boltzmannfactor( double freq, int quanta, double T )
{
  double weight;
  
  double energy = freq * ( quanta + 0.0000 );
  
  double kb = 3.1668114E-6 ; // IN Atomic Unit , i.e. Hartree ...
  
  weight = exp( -1.000 * energy / (kb * T) ) * ( 1.00 - exp( -1.000 * freq / (kb * T) ) );
  
  return( weight );

}


//------------------------------------------------------------//

double boltzmannAvgEnergy( double freq , double T )
{
  double kb = 3.1668114E-6; // IN Atomic Unit , i.e. Hartree ...
  
  double fb_cutoff = 1E-12 ;
  
  double e_cutoff = 1E-12 ;
  
  double gap_cutoff = 1E-12 ;

  int istate = 0 ;

  double avgEnergy = 0.0000 ;
  
  double accfb = 0.0000 ;
  
  double weight , currentStateEnergy ;
 

  if( T == 0.00 )
  {
    avgEnergy = 0.5000 * freq;
  } 
  
  else
  {
    // -----> when istate  = 0 
  
    weight = exp( -1.000 * freq * ( istate + 0.0000 ) / ( kb * T ) ) * ( 1.00 - exp( -1.000 * freq / ( kb * T ) ) ) ;
  
    currentStateEnergy = freq * ( istate + 0.5000 ) * weight ;
  
    //printf("\nCalculating Boltzmann average energy for mode with frequecy % 12.8E ... \n" , freq );
  
    while( weight >= fb_cutoff  &&  currentStateEnergy >= e_cutoff && ( 1 - accfb ) >= gap_cutoff )
    {
      //printf("\nIteration #%d , weight = % 12.8E , currentStateEnergy = % 12.8E ...\n" , istate , weight , currentStateEnergy );
    
      avgEnergy = avgEnergy +  currentStateEnergy ;
    
      accfb = accfb + weight ;
    
      istate ++ ;
    
      weight = exp( -1.000 * freq * ( istate + 0.0000 ) / ( kb * T ) ) * ( 1.00 - exp( -1.000 * freq / ( kb * T ) ) ) ;
    
      currentStateEnergy = freq * ( istate + 0.5000 ) * weight ;
  
    }
  }
  
  
  return( avgEnergy );



}




















//------------------------------------------------------------//


void effquanta( int nmode, int irmode, double temperature, 
                double * freqlist, double cutoff, double * quanta)
{

  int j,istep;
  
  double weight;
  


  for(j=0;j<nmode;j++)
  {
    if( j == irmode-1)
    {
      printf("\nNo. %d mode, this is the ir-excited mode, n = 1\n\n", j);
      
      istep = 1;
      
      *(quanta + j) = sqrt(3.0000); 
    }
    
    else
    {
      istep = 0;
      
      while( (weight = boltzmannfactor( *(freqlist+j), istep, temperature ) )>= cutoff)
      {
        printf("\n No.%d mode, No.%d level, weight = % 10.8E\n\n", j, istep, weight);
        
        *(quanta+j) = (*(quanta+j) ) + sqrt(2.000 * istep + 1.000) * weight;
        
        //*(quanta+j) = (*(quanta+j) ) + sqrt(2.000 * istep + 1.000) * weight;
        
        istep++;
      
      }
   
      //*(quanta + j) = (*(quanta + j)) * (*(freqlist + j));

      istep--;
    
    }
  
    printf("\nNo. %d normal mode, energy level truncated at n = %d\n", j, istep);
  
  }


}


//------------------------------------------------------------//

int dropball( double rnd , double freq , double temperature )
{
  //double boltzmannfactor( double freq, int quanta, double T );
  
  int istate, location;
  
  double pnow, pnext;
  
  double tmp = rnd;
  
  istate = 0; 
  
  if( temperature == 0.000 )
  {
    istate = 0 ;
    
    printf("\nIt is freezing out here ...\n");
  }
  else
  {
     while( ( pnow = boltzmannfactor( freq, istate, temperature ) ) <= tmp )
     {
        //printf("\n\n\t\tBoltzmann factor on level %d is %lf , tmp is %lf ...\n", istate, pnow, tmp);
        //printf("\n\t\tPassing No. %d vib-level, still searching ...\n\n", istate);
    
        tmp = tmp - pnow ;
    
        istate ++ ;
    
     }
  

  }



  //printf("\n\n\t\tBoltzmann factor on level %d is %lf , tmp is %lf ...\n\n", istate, pnow, tmp);
  
  
  if ( istate == 0)
  
    location = 0 ;
  
  else
    
    location = istate ;
  
  
  
  
  return( location );
  
  






}





//------------------------------------------------------------//

void thetasequencegen( int nmode, int ntraj, int * seed, double * theta)
{
int j,k,m;

double random_tmp;

printf("\nAttention! There are en toto %d trajectories ...\n\n",ntraj);

for(j=0;j<nmode;j++)
{
  printf("\n----------------------No. %d Mode --------------------\n\n", j); 

  //seed[j] = time(NULL);

  srand(seed[j]);
    
  for(k=0;k<ntraj;k++)
  {
    random_tmp = (double) rand() / (double) RAND_MAX;
  
    printf("\nNo. %d sequence, No. %d random number, value is % 10.8f ---, corresponding phase angle sinus is % 10.8f ---\n\n", k, j, random_tmp, sin(random_tmp*2*PI) );
  
    *(theta + j*ntraj+k) = random_tmp * 2.000 * PI;
  
  
  
  }

  printf("\n\n");

}


}













//------------------------------------------------------------//

void thetagen( int nmode, int * seed, double * theta)
{
int j,k,m;

double random_tmp;



for(j=0;j<nmode;j++)
{
  printf("\n----------------------No. %d Mode --------------------\n\n", j); 

  //seed[j] = time(NULL);

  srand(seed[j]);
    
  random_tmp = (double) rand() / (double) RAND_MAX;
  
    printf("\nNo. %d sequence, No. %d random number, value is % 10.8f ---, corresponding phase angle sinus is % 10.8f ---\n\n", j, k, random_tmp, sin(random_tmp*2*PI) );
  
  *(theta + j) = random_tmp * 2.000 * PI;
  

  //printf("\n\n");

}


}




//------------------------------------------------------------//

void reducedmasscalc( int natom, double * mass, double * hesscart, char mwcornot, 
                      double * reducedmass)
{
  int ncart, nmode;
  
  int j,k;
  
  double * buffer, * fmwc, * dxdr, * freq;
  
  double tmp;
  
  
  
  ncart=3*natom;

  switch (natom)
  {
    case 1 : nmode = 0; break;

    case 2 : nmode = 1; break;

    default : nmode = ncart - 6; break;
  }

  buffer = malloc( ncart * sizeof(double) );
  /* Here buffer is ncart long, which means the first 5 or 6 number should not be
     be passed into the reducedmass*/
  
  fmwc = malloc( ncart*ncart * sizeof(double)  );
  
  dxdr = malloc( ncart*ncart * sizeof(double)  );
  
  freq = malloc( ncart * sizeof(double) );
  


  if( ( mwcornot == 'Y') || ( mwcornot == 'y') )
  {
    for(j=0; j<ncart*ncart; j++)
    {
      *(fmwc + j) = *(hesscart + j);
    }

  }
  
  else if( ( mwcornot == 'N') || ( mwcornot == 'n') )
  {
    masswt( natom, mass, hesscart, fmwc);
  
  }
  
  else
  {
    printf("\nIn reduced-mass calculation ...\n");
    printf("\nPlease specify whether provided Hessian is mass-weighted or not !!!\n\n");
    
    exit(1);
    
  }


  dsyev_f2c( ncart, fmwc, dxdr, freq  );
  
  printf("\nDiag in red. calc done ...\n");

  dtranspose(ncart, dxdr, dxdr);
  
  
  
  FILE * debug = fopen("dxdr_in_red.deb", "wb+");
  
  doutput(debug, ncart, ncart, dxdr);

  fclose(debug);



  massdivide( 'R', natom, mass, dxdr, dxdr);
  
  
  for(j=0; j<ncart; j++)
  {
    tmp = 0.00000;
    
    for(k=0; k<ncart; k++)
    {
      tmp = tmp + (*(dxdr + k*ncart+j)) * (*(dxdr + k*ncart+j));

    }
  
    *(buffer + j) = 1.000/tmp;
  
  }


  for(j=ncart-nmode; j<ncart; j++)
  {
    *(reducedmass + j-(ncart-nmode) ) = *(buffer + j);
  }





}



//------------------------------------------------------------//

void checkDistribution( int ntraj , int nbin , double * sample , 
                        double *xx , int *yy )
{
  double whomax , whomin ;
  
  double dbin , idnfull;

  double * xAxis;
  
  int idn , itraj , * collective ;
  
  int j,k;
  
  
  
  whomax = dmax( ntraj , sample ); 
  
  whomin = dmin( ntraj , sample );
  
  dbin = ( whomax - whomin ) / nbin ;
  
  printf("\nIN the SAMPLE, max value is % 10.8E , min value is % 10.8E , interval is % 10.8E ...\n" , whomax , whomin , dbin );
  
  
  collective = malloc( nbin * sizeof(int) ); izeros( nbin , 1 , collective );
  
  xAxis = malloc( nbin * sizeof(double) ); dzeros( nbin , 1 , xAxis );
  
     
  for( j = 0 ; j < nbin; j ++ )
  {
     *( xAxis + j ) = whomin + ( j + 0.500 ) * dbin ;
     
     //printf("\nxAxis( %d ) = % 10.8E ... \n", j, whomin + ( j + 0.500 ) * dbin );

  }
  
  
  
  for( j = 0 ; j < ntraj ; j ++ )
  {
    idnfull = ( *( sample + j ) - whomin ) / dbin ; 
    
    idn = (int) idnfull ; 
    
    *( collective + idn ) = 1 + ( *( collective + idn ) );
    
    //printf("\nInside SAMPLE , # %d trajectory , idnfull is % 10.8E , idn is %d ...\n", j , idnfull , idn );
  
  }
  

  for( j = 0 ; j < nbin; j ++ )
  {
     *( xx + j ) = *( xAxis + j );
     
     *( yy + j ) = *( collective + j );
  
  }









}










//------------------------------------------------------------//

int tellIR_range( double position , double front , double tail )
{
  int result ;
  
  if( ( front < 0.00 ) && ( tail < 0.00 ) )
  {
    result = 0 ;
  }
  else
  {
    if ( ( position <= tail ) && ( position >= front ) )
    {
      result = 1 ;
    }
    else
    {
      result = 0 ;
    }
  }
  
  
  return( result );

}



//------------------------------------------------------------//

int tellIR_list( int place , int * list , int length_of_list )
{
  int itmp ; 
  
  int acc = 1 ;
  
  int result ;
  
  for ( itmp = 0 ; itmp < length_of_list ; itmp ++ )
  {
    acc = acc * ( place - *( list + itmp ) );
  }

  if( acc == 0 )
  {
    result = 1 ;
  }
  else
  {
    result = 0 ;
  }
  
  return( result );





}


//------------------------------------------------------------//

int tellActive_list( int place , int * list , int length_of_list )
{
  int itmp ; 
  
  int acc = 1 ;
  
  int result ;
  
  for ( itmp = 0 ; itmp < length_of_list ; itmp ++ )
  {
    acc = acc * ( place - *( list + itmp ) );
  }

  if( acc == 0 )
  {
    result = 1 ;
  }
  else
  {
    result = -1 ;
  }
  
  return( result );


}










//------------------------------------------------------------//


void loadCNDOdipole( char direction , FILE * pcndoLogFile , double * buffer )
{
  // ---> Fock in AO
  
  int irow , icol , iblock , idMO ;
  
  int rowid , colid ;
  
  int nbasis , lmo ;

  int moPerBlock = 15 ;
  
  int moLeftOver ;
  
  int nblock ;
  
  int itmp ;
  
  char cache[ 200 ] ;
  
  char stmp[ 200 ] ;
  
     
  
  double * dipole ;
  
  // ---> Learn NBasis ...
  
  rewind( pcndoLogFile ) ;
  
  fsearch( pcndoLogFile , "Basis" ) ;
  
  fscanf( pcndoLogFile , "%s" , stmp );
  
  if( strcmp( stmp , "fns" ) == 0 )
  {
    fscanf( pcndoLogFile , "%s" , stmp ) ; 
   
    fscanf( pcndoLogFile , "%d" , &nbasis );
  }
  else
  {
    fsearch( pcndoLogFile , "SHELL" );
    
    fsearch( pcndoLogFile , "TOTAL" ) ;
    
    fscanf( pcndoLogFile , "%d" , &nbasis ) ;
  
  }
  
  printf("\nDimension of the basis set used in this calculation is %d ! \n" , nbasis ) ;
  
  lmo = nbasis * nbasis ;
  
  moLeftOver = nbasis % moPerBlock ;
  
  nblock = ( nbasis - moLeftOver ) / moPerBlock ;
  
  
  
  if( moLeftOver != 0 )
  {
    nblock = nblock + 1 ;
  }      
  else
  {
    moLeftOver = moPerBlock ;
  }
  
  printf("\nThere are %d blocks and %d orbitals leftover ... \n" , nblock , moLeftOver ) ;
  

  rewind( pcndoLogFile ) ;
  
  dipole = ( double * ) calloc( lmo , sizeof( double ) ) ; dzeros( lmo , 1 , dipole ) ;
  



  rewind( pcndoLogFile ) ;
  
  fsearch( pcndoLogFile , "DIPOLE" ) ;
  
  if( direction == 'X' || direction == 'x' )
  {
    fsearch( pcndoLogFile , "Dipole" ) ;
    
    printf("\nReading X-Direction Dipole Integral Matrix ...\n\n");
  }
  else if( direction == 'Y' || direction == 'y' )
  {
    fsearch( pcndoLogFile , "Dipole" ) ;
    
    fsearch( pcndoLogFile , "Dipole" ) ;
    
    printf("\nReading Y-Direction Dipole Integral Matrix ...\n\n");

  }
  else if( direction == 'Z' || direction == 'z' )
  {
    fsearch( pcndoLogFile , "Dipole" ) ;
    
    fsearch( pcndoLogFile , "Dipole" ) ;
    
    fsearch( pcndoLogFile , "Dipole" ) ;
    
    printf("\nReading Z-Direction Dipole Integral Matrix ...\n\n");
  }
  else
  {
    printf("\nPlease only say one of these choices : X ( or x ) , Y ( or y ) OR Z ( or z )\n");
    
    exit( 152 ) ; 
  }
  

  fskip( pcndoLogFile , 3 ) ;
  
  for( iblock = 0 ; iblock < nblock - 1 ; iblock ++ )
  {
    printf("\n# %d Block ... \n" , iblock ) ;
    
    /* First 15 Lines */
    
    for( irow = 0 ; irow < moPerBlock ; irow ++ )
    {
      rowid = iblock * moPerBlock + irow ;
      
      // printf("\nReading # %d line of Dipole ...\n" , rowid + 1 ) ;
      
      fscanf( pcndoLogFile , "%d" , &itmp ) ;
      
      if( rowid != itmp - 1 )
      {
        printf("\nSomething is wrong with reading the 15-triangle in %d ( Human label) block ...\n" , iblock + 1 ) ;
        
        printf("\nLog file says it is the No. %d line ... you are actually @ No. %d line ...\n" , itmp , rowid + 1 ) ;
        
        exit( 382 ) ;
      }
      
      fscanf( pcndoLogFile , "%s" , stmp ) ;  fscanf( pcndoLogFile , "%s" , stmp ) ;  fscanf( pcndoLogFile , "%s" , stmp ) ;
      
      for( icol = 0 ; icol < irow + 1 ; icol ++ )
      {
        colid = iblock * moPerBlock + icol ;
        
        fscanf( pcndoLogFile , "%s" , cache ) ;
        
        *( dipole + nbasis * rowid + colid ) = atof( cache ) ;
        
      }
      
      // printf("\nThe last number of this line is %lf ... " , atof( cache ) ) ;
    

    }
    
    /* The Rest */
    
    for( irow = rowid + 1 ; irow < nbasis ; irow ++ )
    {
      fscanf( pcndoLogFile , "%d" , &itmp ) ;
      
      // printf("\nReading # %d line of Dipole ...\n" , irow + 1 ) ;
      
      if( irow != itmp - 1 )
      {
        printf("\nSomething is wrong with reading the big-rectangle in %d ( Human label) block ...\n" , iblock + 1 ) ;
        
        exit( 382 ) ;
      }
      
      fscanf( pcndoLogFile , "%s" , stmp ) ;  fscanf( pcndoLogFile , "%s" , stmp ) ;  fscanf( pcndoLogFile , "%s" , stmp ) ;
      

      for( icol = iblock * moPerBlock ; icol < ( iblock + 1 ) * moPerBlock ; icol ++ )
      {
        fscanf( pcndoLogFile , "%s" , cache ) ;
        
        *( dipole + nbasis * irow + icol ) = atof( cache ) ;
      }
    
      // printf("\nThe last number of this line is %lf ... " , atof( cache ) ) ;
    
    }
  
    
    fskip( pcndoLogFile , 4 ) ;
    
  }
  
  /* Now only a moLeftOver-by-moLeftOver triangle left */
  
  //fskip( pcndoLogFile , 4 ) ;
  
  printf("\nTHE FINAL BLOCK\n" ) ;
  
  for( irow = 0 ; irow < moLeftOver ; irow ++ )
  {
    rowid = ( nblock - 1 ) * moPerBlock + irow ;
    
    fscanf( pcndoLogFile , "%d" , &itmp ) ;
    
    if( rowid != itmp - 1 )
    {
      printf("\nSomething is wrong with reading the 15-triangle in the final block ...\n" ) ;
      
      printf("\nLog file says it is the No. %d line ... you are actually @ No. %d line ...\n" , itmp , rowid + 1 ) ;
      
      exit( 382 ) ;
    }
    
    fscanf( pcndoLogFile , "%s" , stmp ) ;  fscanf( pcndoLogFile , "%s" , stmp ) ;  fscanf( pcndoLogFile , "%s" , stmp ) ;
  
    for( icol = 0 ; icol < irow + 1 ; icol ++ )
    {
      colid = ( nblock - 1 ) * moPerBlock + icol ;
      
      fscanf( pcndoLogFile , "%s" , cache ) ;
      
      *( dipole + nbasis * rowid + colid ) = atof( cache ) ;
      
    }
  
    // printf("\nThe last number of this line is %lf ... " , atof( cache ) ) ;
  }
  
  //void dtri2sym( char uplo, int lda, double * tri, double * sym)
  
  //dtri2sym( 'L' , nbasis , fao , fao ) ;
  
  for( irow = 0 ; irow < nbasis ; irow ++ )
  {
    for( icol = irow ; icol < nbasis ; icol ++ )
    {
      *( dipole + irow * nbasis + icol ) = *( dipole + icol * nbasis + irow ) ;
    
    }
  
  
  }
  
  for( itmp = 0 ; itmp < nbasis * nbasis ; itmp ++ )
  {
    *( buffer + itmp ) = *( dipole + itmp ) ;
  }
  
  
  
  
  
  

}

//------------------------------------------------------------//


























