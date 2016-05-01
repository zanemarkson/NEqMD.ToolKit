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

#ifndef KBAU
#define KBAU 3.1668107400994983e-06
#endif

/*
double hcPsi( int n , double xp , double w , double A )
{
  double hc = A * hermite( n , xp ) * exp( -0.50 * xp * xp ) ;
  //printf( "\nhermite(%d , %lf) = %lf\n" , n , xp , hermite( n , xp ) ) ;
  return( hc ) ;

}

double gaussianDistribution( double sigma , double mu , double x )
{
  double g = 1.00 / sigma / sqrt( 2.00 * PI ) * exp( -1.00 * ( x - mu ) * ( x - mu ) / 2.00 / ( sigma * sigma ) ) ;
  
  return g ;

}
*/



double wignerPosition( int n , double p , double w , double A )
{
  double scaledOmega = w ;
  double leftCutOff = -1.00 * sqrt( ( 2.0 * n + 1.00 ) / scaledOmega );
  int nstepOneSide = 1200 ;
  int nsteps = 2 * nstepOneSide ;
  double stepsize = -1.00 * leftCutOff / nstepOneSide ;
  //printf("\nleftCutOff = %lf , stepsize = %lf ...\n\n" , leftCutOff , stepsize ) ;
  
  int istep = 0 ;
  double accumulatedP = 1E-16 ;
  double wavefunc ;
  
  for( istep = 0 ; accumulatedP < p ; istep ++ )
  {
    //printf("\nBefore stepping % 10.6E forward , accumulatedP = % 10.8E ...\n" , stepsize , accumulatedP ) ;
    wavefunc = hcPsi( n , ( leftCutOff + stepsize * istep ) * sqrt(w) , w , A ) ;
    //wavefunc = A * hermite( n , ( leftCutOff + stepsize * istep ) * sqrt(w) ) * exp( -0.50 * ( leftCutOff + stepsize * istep ) * sqrt(w) * ( leftCutOff + stepsize * istep ) * sqrt(w) ) ;
    //printf("\nAdded % 12.8E * % 12.8E * % 12.8E = % 12.8E ... \n\n" , stepsize , wavefunc , wavefunc , stepsize * wavefunc * wavefunc ) ;
    accumulatedP = accumulatedP + stepsize * wavefunc * wavefunc ;
    //printf("\nUp to istep = %d , accumulatedP = % 10.8E , difference towards target is % 16.12E... \n" , istep , accumulatedP , p - accumulatedP ) ;
  }
  
  return ( leftCutOff + istep * stepsize ) ;
 
}

double wignerVelocity( int n , double p , double w , double A )
{
  double scaledOmega = w ;
  double leftCutOff = -1.00 * sqrt( ( 2.0 * n + 1.0 ) * scaledOmega ) ;
  int nstepOneSide = 400 ;
  int nsteps = 2 * nstepOneSide ;
  double stepsize = -1.00 * leftCutOff / nstepOneSide ;
  
  int istep = 0 ;
  double accumulatedP = 1E-16  ;
  double wavefunc = 0.00 ;
  
  for( istep = 0 ; accumulatedP < p ; istep ++ )
  {
    wavefunc = hcPsi( n , ( leftCutOff + stepsize * istep ) / sqrt(w) , w , A ) ;
    accumulatedP = accumulatedP + stepsize * wavefunc * wavefunc ;
    //printf("\nUp to istep = %d , accumulatedP = % 10.8E , difference towards half is % 16.12E... \n" , istep , accumulatedP , 0.50 - accumulatedP ) ;
  }
  
  return ( leftCutOff + istep * stepsize ) ;

}






double canonicalPosition( int n , double p , double w , double T , double A )
{
  double scaledOmega = w ;
  double ctp = sqrt( ( 2.0 * n + 1.00 ) / scaledOmega );
  double leftCutOff = -1.00 * ctp;
  //double leftCutOff = -60.00 ;
  
  int nstepOneSide = 1200 ;
  int nsteps = 2 * nstepOneSide ;
  double stepsize = -1.00 * leftCutOff / nstepOneSide ;
  
  double beta = 1.00 / ( KBAU * T ) ;
  double u = tanh( scaledOmega * beta / 2.00 ) ;
  double mu = 0.00 ;
  double sigma = sqrt( ( 1.00 / scaledOmega ) / ( 2.00 * u ) ) ;
  
  int istep = 0 ;
  double accumulatedP = 1E-16 ;
  double wavefunc = 0.00 ;
  
  for( istep = 0 ; accumulatedP < p ; istep ++ )
  {
    wavefunc = gaussianDistribution( sigma , mu , ( leftCutOff + stepsize * istep ) ) ;
    //printf("\nAdded % 12.8E ... \n\n" , hc ) ;
    accumulatedP = accumulatedP + stepsize * wavefunc * A ;
    //printf("\nUp to istep = %d , accumulatedP = % 10.8E , difference towards half is % 16.12E... \n" , istep , accumulatedP , 0.50 - accumulatedP ) ;
  }
  
  return( leftCutOff + ( istep - 1 ) * stepsize ) ;
 
}

double canonicalVelocity( int n , double p , double w , double T , double A )
{
  double scaledOmega = w ;
  double ctp = sqrt( ( 2.0 * n + 1.0 ) * scaledOmega ) ; 
  double leftCutOff = -1.00 * ctp ;
  //double leftCutOff = -1.00 ;

  int nstepOneSide = 400 ;
  int nsteps = 2 * nstepOneSide ;
  double stepsize = -1.00 * leftCutOff / nstepOneSide ;
  
  double beta = 1.00 / ( KBAU * T ) ;
  double u = tanh( scaledOmega * beta / 2.00 ) ;
  double mu = 0.00 ;
  double sigma = sqrt( ( 1.00 * scaledOmega ) / ( 2.00 * u ) ) ;
  
  int istep = 0 ;
  double accumulatedP = 1E-16  ;
  double wavefunc = 0.00 ;
  
  for( istep = 0 ; accumulatedP < p ; istep ++ )
  {
    wavefunc = gaussianDistribution( sigma , mu , ( leftCutOff + stepsize * istep ) ) ;
    
    accumulatedP = accumulatedP + stepsize * wavefunc * A ;
    //printf("\nUp to istep = %d , accumulatedP = % 10.8E , difference towards half is % 16.12E... \n" , istep , accumulatedP , 0.50 - accumulatedP ) ;
  }
  
  return( leftCutOff + ( istep - 1 ) * stepsize ) ;

}




double xiGenWigner( int n , double ZPEscalingFactor , double rd , double omega , double A , double * angle )
{
  double xi ;
  double angle_local ;
  double p = rd / (double) IM ;
  double scaledOmega = 0.00 ;
  
  if( n <= 2 )
  {
    scaledOmega = ( n + 0.50 * ZPEscalingFactor ) / ( n + 0.50 ) * omega ;
    printf("\nNormalization factor is % 12.8E ...\n" , A ) ;
    xi = wignerPosition( n , p , scaledOmega , sqrt( A ) ) ;
    printf("\nTarget Probability is %lf , found position is %lf ...\n" , p , xi ) ;
    * angle = -1.00 * PI ;
  }
  else
  {
    angle_local = p * 2.000 * PI ;
    * angle = angle_local ;
    xi = sqrt( 2.000 * n + 1.000 * ZPEscalingFactor ) / sqrt( omega ) * cos(angle_local) ;
  }

  return( xi ); 

}


double xidotGenWigner( int n , double ZPEscalingFactor , double rd , double omega , double A , double * angle )
{
  double xidot ;
  double angle_local ;
  double p = rd / (double) IM ;
  double scaledOmega = 0.00 ;
  
  if( n <= 2 )
  {
    scaledOmega = ( n + 0.50 * ZPEscalingFactor ) / ( n + 0.50 ) * omega ;
    xidot = wignerVelocity( n , p , scaledOmega , sqrt( A ) ) ;
    * angle = -1.00 * PI ;
  }
  else
  {
    angle_local = p * 2.000 * PI ;
    * angle = angle_local ;
    xidot = sqrt( 2.000 * n + 1.000 * ZPEscalingFactor ) * sqrt( omega ) * sin(angle_local);
  }
  
  return( xidot ) ;

}





double xiGenCanonical( int n , double ZPEscalingFactor , double rd , double omega , double A , double T , double * angle )
{
  double xi ;
  double angle_local ;
  double p = rd / (double) IM ;
  double scaledOmega = 0.00 ;
  
  if( n <= 2 )
  {
    scaledOmega = ( n + 0.50 * ZPEscalingFactor ) / ( n + 0.50 ) * omega ;
    printf("\nNormalization factor is % 12.8E ...\n" , A ) ;
    xi = canonicalPosition( n , p , scaledOmega , T , A ) ;
    printf("\nTarget Probability is %lf , found position is %lf ...\n" , p , xi ) ;
    * angle = -1.00 * PI ;
  }
  else
  {
    angle_local = p * 2.000 * PI ;
    * angle = angle_local ;
    xi = sqrt( 2.000 * n + 1.000 * ZPEscalingFactor ) / sqrt( omega ) * cos(angle_local) ;
  }

  return( xi ); 

}


double xidotGenCanonical( int n , double ZPEscalingFactor , double rd , double omega , double A , double T , double * angle )
{
  double xidot ;
  double angle_local ;
  double p = rd / (double) IM ;
  double scaledOmega = 0.00 ;
  
  if( n <= 2 )
  {
    scaledOmega = ( n + 0.50 * ZPEscalingFactor ) / ( n + 0.50 ) * omega ;
    xidot = canonicalVelocity( n , p , scaledOmega , T , A ) ;
    * angle = -1.00 * PI ;
  }
  else
  {
    angle_local = p * 2.000 * PI ;
    * angle = angle_local ;
    xidot = sqrt( 2.000 * n + 1.000 * ZPEscalingFactor ) * sqrt( omega ) * sin(angle_local);
  }
  
  return( xidot ) ;

}










void icgen( int natomselect , int * irstatus_list , long * seeda , long * seedn , double fixedPhaseAngle ,
            double nquanta , double irfront , double irtail , int asinitial , int asfinal , 
            double T , double ZPEscalingFactor , double * mass , double * vibfreq, 
            double * cart_q0, double * U , double * cart_q,  double * cart_p , 
            double * wignerPositionNormArray , double * wignerVelocityNormArray , 
            double * canonicalPositionNormArray , double * canonicalVelocityNormArray ,
            int debuggingMode )
{

  //-o-o-o-o-o-o-o some constants... o-o-o-o-o-o-o-o-//
  
  double done = 1.0000; int ione = 1;
  
  double dzero = 0.0000; int izero = 0;
  
  
  
  //-o-o-o-o-o-o-o Variables... o-o-o-o-o-o-o-o-//
  
  
  
  int ncartselect , nmodeselect ;
  
  int imode , itraj , irstatus ;

  int iatom , icart ;

  int itmp; double dtmp;
  
  int atomID ;
  
  
  int * irmodelist ;
  
  long seedacurr , seedncurr ;
  
  double level;
  
  double * levelassgn ;
  
  double * mass3n;
  
  double angle, locator ;
  
  double xicurr , xidotcurr;
  
  double * xi , * xidot;
  
  double * kineticEnergy , * potentialEnergy ;
  
  
  
  
  
  //-o-o-o-o-o-o-o-o For debugging output o-o-o-o-o-o-o-o-o-//
  
  register FILE * debug;
  
  
  
  
  
  //-o-o-o-o-o-o-o Read Gaussian output and load basic parameters... o-o-o-o-o-o-o-o-//

  //readgauss( "natom", plog, &natom, &dtmp );
  
  
  

  ncartselect = 3 * natomselect ;

  switch ( natomselect )
  {
    case 1 : nmodeselect = 0 ; break;
  
    case 2 : nmodeselect = 1 ; break;
  
    default : nmodeselect = ncartselect - 6 ; break;
  
  
  }
  

  
  //-o-o-o-o-o-o-o-o Allocating space and set up original seeds... o-o-o-o-o-o-o-o-o-o-o-//
  //locator = (double) ranzm( &seedncurr ) / (double) IM ;
  
  
  
  
  xi = calloc( nmodeselect * 1 , sizeof(double) ); dzeros( nmodeselect , 1 , xi );
  
  xidot = calloc( nmodeselect * 1 , sizeof(double) ); dzeros( nmodeselect , 1 , xidot );
  

  
  if( debuggingMode == YES )
  {
    debug = fopen( "wignerPositionNormArray.deb" , "wb+" );

    doutput( debug , nmodeselect , 3 , wignerPositionNormArray );

    fclose( debug );
  
  
  
    debug = fopen( "wignerVelocityNormArray.deb" , "wb+" );

    doutput( debug , nmodeselect , 3 , wignerVelocityNormArray );

    fclose( debug );


  }
  
    
      
  
  /*
  if( debuggingMode == YES )
  {
    debug = fopen( "canonicalPositionNormArray.deb" , "wb+" );

    doutput( debug , nmodeselect , 3 , canonicalPositionNormArray );

    fclose( debug );
  
  
  
    debug = fopen( "canonicalVelocityNormArray.deb" , "wb+" );

    doutput( debug , nmodeselect , 3 , canonicalVelocityNormArray );

    fclose( debug );


  }
  */
  
  //-o-o-o-o-o-o-o Generating normal mode coordinate and momentum ... o-o-o-o-o-o-o-o-//
  
  double * energy ;
  
  double * angleList ;
  
  double * arcdq , * arcq , * arcp ;
  
  double * xiveccurr , * xidotveccurr ;
  
  double * dqcurr , * pcurr ; 
  
  
  energy = calloc( nmodeselect * 1 , sizeof( double ) ); dzeros( nmodeselect , 1 , energy );
  
  angleList = calloc( nmodeselect * 1 , sizeof( double ) ); dzeros( nmodeselect , 1 , angleList );
  
  xiveccurr = calloc( nmodeselect , sizeof( double ) ); 
  
  xidotveccurr = calloc( nmodeselect , sizeof( double ) );
  
  kineticEnergy = calloc( nmodeselect * 1 , sizeof( double ) );
  
  potentialEnergy = calloc( nmodeselect * 1 , sizeof( double ) );
  
  levelassgn =  calloc( nmodeselect * 1 , sizeof( double ) );
  
  for( imode = 0 ; imode < nmodeselect ; imode ++ ) *(levelassgn + imode ) = -1 ;
  
  dqcurr = calloc( ncartselect , sizeof( double ) ); 
  
  pcurr = calloc( ncartselect , sizeof( double ) ); 
  
  arcdq = calloc( ncartselect * 1 , sizeof( double ) ); dzeros( ncartselect , 1 , arcdq );
  
  arcq = calloc( ncartselect * 1 , sizeof( double ) ); dzeros( ncartselect , 1 , arcq );
  
  arcp = calloc( ncartselect * 1 , sizeof( double ) ); dzeros( ncartselect , 1 , arcp );
  
  

  
  dzeros( nmodeselect , 1 , xiveccurr );
    
  dzeros( nmodeselect , 1 , xidotveccurr );
    
  dzeros( ncartselect , 1 , dqcurr );
    
  dzeros( ncartselect , 1 , pcurr );
  
  
  
  
  
    
  // Only for debugging purposes ... printing out the execution process ... //
  
  //printf("\n-o-o-o-o-o-o-o-o Working on Traj# %d -o-o-o-o-o-o-o\n" , itraj );
  
  //printf("\nAnything wrong here?\n");

  //fprintf( debug, "\n-o-o-o-o-o-o-o-o Working on Traj# %d -o-o-o-o-o-o-o\n", itraj );

 
  // Only for debugging purposes ... printing out the execution process ... 
  
  
  
  double rdcurr = 0.00 ;
  double scaledOmega = 0.00 ;
  //for ( imode = asfinal - 1 ; imode > asinitial - 2 ; imode -- ) // asinitial and asfinal are in Human Labelling ...
  for ( imode = asinitial - 1 ; imode < asfinal ; imode ++ ) // asinitial and asfinal are in Human Labelling ...
  //for( imode = 0 ; imode < nmodeselect ; imode ++ )
  {
    //fprintf( debug, "\n\t\t\tWorking on Mode# %d\t......\n", imode + 1 );
    
    printf( "\n\t\t\tWorking on Mode# %d\t......\n", imode + 1 );
    
    xicurr = 0.0000;  xidotcurr = 0.0000;   
    
    irstatus = *( irstatus_list + imode ) ;
    
    int ilevel ;
  
    //printf("\n-o-o-o-o-o-o-o-o-o-o Working on NORMAL MODE %d -o-o-o-o-o-o-o-o-\n", imode );
  
    if ( irstatus == 1 )
    {
      printf("\n\nBOOM... This is THE IR-EXCITED MODE ... \n\n");
    
      //for ( itraj = 0 ; itraj < ntraj ; itraj ++ )
      //{
        // --------  Initializing this trajectory -------- //
      
        //printf("\nNORMAL MODE # %d ... \n" , imode );
      
        seedacurr = *( seeda + imode ); 
        
        //printf("\nNow seed.theta.current is %ld ... ", seedacurr );
        
        seedncurr = *( seedn + imode ); 
        
        //printf("Now seed.n.current is %ld\n\n", seedncurr);

      
        
        // --------  Generating energy level location -------- //
        
        locator = (double) ranzm( &seedncurr ) / (double) IM ;
    
        //printf("\nNow locator.current[%d] is ... it does not matter ...", imode );
        
        level = nquanta ; 
        
        *( levelassgn + imode ) = level ;
     
        //*( levelassgn + imode * ntraj + itraj ) = level ; 

        printf("\nNow level.current is %lf\n\n", level);
        
        // --------  Generating phase angle -------- //
        
        rdcurr = (double) ranzm( &seedacurr ) ;
        //angle = rdcurr / (double) IM * 2.000 * PI ;
        //angle = -1.00 * PI ;
        

        //printf("\nNow theta.current[%d] is %lf ... ", imode , angle );

        //printf("After update seed.theta[%d] is %ld\n", imode, seedacurr);
      
            
        // --------  Calculating and Accumulating Xi and Xidot ------ //
        // !!!!!!!!!!!!  Mar. 31st, 2015, modified for fixedPhaseAngle !!!!!!!!!! //
        
        if( fixedPhaseAngle == 9.00 )
        {
          //xicurr = sqrt( 2.000 * level + 1.000 * ZPEscalingFactor ) / sqrt( *( vibfreq + imode ) ) * cos(angle) ;
        
          //xidotcurr = sqrt( 2.000 * level + 1.000 * ZPEscalingFactor ) * sqrt( *( vibfreq + imode ) ) * sin(angle);
          ilevel = (int)level ;
          printf("\nNormalization factor is % 10.8E and % 10.8E ...\n\n" , *( wignerPositionNormArray + 3 * imode + ilevel ) , *( wignerVelocityNormArray + 3 * imode + ilevel ) ) ;
          xicurr = xiGenWigner( level , ZPEscalingFactor , rdcurr , *( vibfreq + imode ) , *( wignerPositionNormArray + 3 * imode + ilevel ) , &angle ) ;
          xidotcurr = xidotGenWigner( level , ZPEscalingFactor , rdcurr , *( vibfreq + imode ) , *( wignerVelocityNormArray + 3 * imode + ilevel ) , &angle ) ;
          scaledOmega = ( level + 0.5 * ZPEscalingFactor ) / ( level + 0.5 ) * ( *( vibfreq + imode ) ) ;
          *( kineticEnergy + imode  ) = 0.50000 * xidotcurr * xidotcurr ;
          *( potentialEnergy + imode  ) = 0.50000 * xicurr * xicurr * ( scaledOmega * scaledOmega ) ;

        
        }
        else
        {
          angle = fixedPhaseAngle * PI ; 
          
          xicurr = level ; //sqrt( 2.000 * level + 1.000 * ZPEscalingFactor ) / sqrt( *( vibfreq + imode ) ) * cos( fixedPhaseAngle * PI ) ;
         
          xidotcurr = 0.00 ; //sqrt( 2.000 * level + 1.000 * ZPEscalingFactor ) * sqrt( *( vibfreq + imode ) ) * sin( fixedPhaseAngle * PI );
          
          *( kineticEnergy + imode  ) = 0.50000 * xidotcurr * xidotcurr ;
        
          *( potentialEnergy + imode  ) = 0.50000 * xicurr * xicurr * ( *( vibfreq + imode ) ) * ( *( vibfreq + imode ) ) ;
        
          
        }
        
        *( energy + imode  ) = ( *( vibfreq + imode ) ) * ( level + 0.5000 * ZPEscalingFactor );
        
        *( xiveccurr + imode ) =  xicurr ;
    
        *( xidotveccurr + imode ) =  xidotcurr ;
        
        *( angleList + imode ) = angle ;
     
        
    
    
        // --------  Updating the seeds --------//
       
        *( seeda + imode ) = seedacurr ; 
      
        //printf("\nAfter update seed.theta[%d] is %d\n", imode, seedacurr);
    
        *( seedn + imode ) = seedncurr;
    
        //printf("\nAfter update seed.n[%d] is %d\n", imode, seedncurr);
        
      
        // --------  Calculating the energy --------//
            
        printf("\nCurrent energy is %lf with vibration frequency % 10.8E ...\n\n", 
            ( *( vibfreq + imode ) ) * ( level + 0.50 * ZPEscalingFactor ) , *( vibfreq + imode ) );
      
        // --------  And ... what else ? --------//
  
      
      

    //}
    }
    else if( irstatus == -1 )
    {
      printf("\n\nThis is FROZEN MODE ... \n\n");
      
       // --------  Initializing this trajectory --------//
        
        //printf("\nNORMAL MODE # %d ... \n" , imode );
    
        seedacurr = *( seeda + imode ); 
        //printf("\nNow seed.theta.current is %ld ... ", seedacurr);
       
        seedncurr = *( seedn + imode ); 
        //printf("Now seed.n.current is %ld\n\n", seedncurr);
    
        
        rdcurr = (double) ranzm( &seedacurr ) ;
        
        
        // --------  Generating energy level location -------- //
        
        locator = (double) ranzm( &seedncurr ) / (double) IM ;
    
        //printf("\nNow locator.current[%d] is %lf ... ", imode , locator);
    
        //printf("After update seed.n[%d] is %ld ... ", imode, seedncurr);
    
        /*#1*///fprintf( debug , "\nGood till here...\n" );
    
        level = -1.00 ;
    
        *( levelassgn + imode ) = level ;
	
	    printf("Now level.current[%d] is %lf\n\n", imode , level);
        
        
        /*#2*///fprintf( debug, "\nGood till here...\n" );
        
        // --------  Generating phase angle -------- //
        
        angle = (double) ranzm( &seedacurr ) / (double) IM * 2.000 * PI ;
    
        //printf("\nNow theta.current[%d] is %lf ... ", imode , angle);
    
        //printf("After update seed.theta[%d] is %ld\n\n", imode, seedacurr);
    
        /* --------  Calculating and Accumulating Xi and Xidot --------*/
        
        xicurr = 0.000000 ;
    
        xidotcurr = 0.000000 ;
    
    
        *( xiveccurr + imode ) =  xicurr ;
    
        *( xidotveccurr + imode ) =  xidotcurr ;
        
    
        /*#3*///fprintf(debug, "\nGood till here...\n");
    
        /* --------  Updating the seeds --------*/
       
        *( seeda + imode ) = seedacurr ; 
    
        //printf("\nAfter update seed.theta[%d] is %d\n", imode, seedacurr);
    
        *( seedn + imode ) = seedncurr;
    
        //printf("\nAfter update seed.n[%d] is %d\n", imode, seedncurr);
        
       
        /* --------  Calculating the energy --------*/
  
        
        *( energy + imode  ) = 0.000000 ;
        
        *( angleList + imode ) = 0.000000 ;
  
        *( kineticEnergy + imode  ) = 0.00000  ;
        
        *( potentialEnergy + imode  ) = 0.000000 ;
        
        /*#4*///fprintf(debug, "\nGood till here...\n");
        
        
        //printf("\nCurrent energy is %lf with vibration frequency % 10.8E ...\n\n", 
        //      ( *( vibfreq + imode ) ) * ( level + 0.5000 ) , *(vibfreq + imode));
        
       /* --------  And ... what else ? --------*/

    }
    else
    {
      printf("\n\nThis is NOT IR-EXCITED MODE ... \n\n");
      
       // --------  Initializing this trajectory --------//
        
        //printf("\nNORMAL MODE # %d ... \n" , imode );
    
        seedacurr = *( seeda + imode ); 
        //printf("\nNow seed.theta.current is %ld ... ", seedacurr);
       
        seedncurr = *( seedn + imode ); 
        //printf("Now seed.n.current is %ld\n\n", seedncurr);
    
        
        // --------  Generating energy level location -------- //
        
        locator = (double) ranzm( &seedncurr ) / (double) IM ;
    
        //printf("\nNow locator.current[%d] is %lf ... ", imode , locator);
    
        //printf("After update seed.n[%d] is %ld ... ", imode, seedncurr);
    
        /*#1*///fprintf( debug , "\nGood till here...\n" );
    
        level = dropball( locator , *( vibfreq + imode ) , T );
  
        //==============> 02/03/2016 For nquanta > 2 , reset it to be 2 <============= //
	
		level = ( level > 2.00 ) ? 2.00 : level ;
        
	//==============> <===========================================> <============= //
    
        *( levelassgn + imode ) = level ;
	
        printf("Now level.current[%d] is %lf\n\n", imode , level);
        
        
        /*#2*///fprintf( debug, "\nGood till here...\n" );
        
        
        // --------  Generating phase angle -------- //
        
        rdcurr = (double) ranzm( &seedacurr ) ;
        
        //angle = (double) ranzm( &seedacurr ) / (double) IM * 2.000 * PI ;
    
        //printf("\nNow theta.current[%d] is %lf ... ", imode , angle);
    
        //printf("After update seed.theta[%d] is %ld\n\n", imode, seedacurr);
    
        
        
        /* --------  Calculating and Accumulating Xi and Xidot --------*/
        // !!!!!!!!!!!!  Mar. 31st, 2015, modified for fixedPhaseAngle !!!!!!!!!! //

        
        if( fixedPhaseAngle == 9.00 )
        {
          //xicurr = sqrt( 2.000 * level + 1.000 * ZPEscalingFactor ) / sqrt( *( vibfreq + imode ) ) * cos(angle) ;
        
          //xidotcurr = sqrt( 2.000 * level + 1.000 * ZPEscalingFactor ) * sqrt( *( vibfreq + imode ) ) * sin(angle);
          ilevel = (int)level ;
          xicurr = xiGenWigner( level , ZPEscalingFactor , rdcurr , *( vibfreq + imode ) , *( wignerPositionNormArray + 3 * imode + ilevel ) , &angle ) ;
          xidotcurr = xidotGenWigner( level , ZPEscalingFactor , rdcurr , *( vibfreq + imode ) , *( wignerVelocityNormArray + 3 * imode + ilevel ) , &angle ) ;
          scaledOmega = ( level + 0.5 * ZPEscalingFactor ) / ( level + 0.5 ) * ( *( vibfreq + imode ) ) ;
          *( kineticEnergy + imode  ) = 0.50000 * xidotcurr * xidotcurr ;
          *( potentialEnergy + imode  ) = 0.50000 * xicurr * xicurr * ( scaledOmega * scaledOmega ) ;

        
        }
        else
        {
          angle = fixedPhaseAngle * PI ; 
          
          xicurr = level ; //sqrt( 2.000 * level + 1.000 * ZPEscalingFactor ) / sqrt( *( vibfreq + imode ) ) * cos( fixedPhaseAngle * PI ) ;
         
          xidotcurr = 0.00 ; //sqrt( 2.000 * level + 1.000 * ZPEscalingFactor ) * sqrt( *( vibfreq + imode ) ) * sin( fixedPhaseAngle * PI );
          
          *( kineticEnergy + imode  ) = 0.50000 * xidotcurr * xidotcurr ;
        
          *( potentialEnergy + imode  ) = 0.50000 * xicurr * xicurr * ( *( vibfreq + imode ) ) * ( *( vibfreq + imode ) ) ;
        
          
        }
        
        *( energy + imode  ) = ( *( vibfreq + imode ) ) * ( level + 0.5000 * ZPEscalingFactor );
        
        *( xiveccurr + imode ) =  xicurr ;
    
        *( xidotveccurr + imode ) =  xidotcurr ;
        
        *( angleList + imode ) = angle ;
     
        
    
        /*#3*///fprintf(debug, "\nGood till here...\n");
    
        /* --------  Updating the seeds --------*/
       
        *( seeda + imode ) = seedacurr ;   
    
        //printf("\nAfter update seed.theta[%d] is %d\n", imode, seedacurr);
    
        *( seedn + imode ) = seedncurr;
    
        //printf("\nAfter update seed.n[%d] is %d\n", imode, seedncurr);
        
       
        /* --------  Calculating the energy --------*/
  
        
        *( energy + imode  ) = ( *( vibfreq + imode ) ) * ( level + 0.5000 * ZPEscalingFactor );
        
        *( angleList + imode ) = angle ;
  
        /*#4*///fprintf(debug, "\nGood till here...\n");
        
        
        //printf("\nCurrent energy is %lf with vibration frequency % 10.8E ...\n\n", 
        //      ( *( vibfreq + imode ) ) * ( level + 0.5000 ) , *(vibfreq + imode));
        
       /* --------  And ... what else ? --------*/


      }

 
  //}

  //printf("\nAccumulative.Energy.%d after %d trajectory is % 10.8E ... \n\n", imode , itraj+1 , *(energy + imode) );
  


  
  } //=====> Loop with imode ends ...

  

  for( imode = 0 ; imode < nmodeselect ; imode ++ )
  {
    *( xi + imode ) = *( xiveccurr + imode ) ;

    *( xidot + imode ) = *( xidotveccurr + imode );
  
  }





  printf("\nDone dealing with normal modes for this trajectory ... Starting transformation into Cartesian ... \n");

  if( ( nmodeselect != 1 ) && ( natomselect != 1 ) )
  {
    dgemv_( "T" , &nmodeselect , &ncartselect , &done , U , &nmodeselect , xiveccurr , &ione , &dzero , dqcurr , &ione );

    dgemv_( "T" , &nmodeselect , &ncartselect , &done , U , &nmodeselect , xidotveccurr , &ione , &dzero , pcurr , &ione );
  }
  
  else
  {
    for( icart = 0 ; icart < ncartselect ; icart ++)
    {
      *( dqcurr + icart ) = *( U + icart ) * ( * xiveccurr );
      
      *( pcurr + icart ) = *( U + icart ) * ( * xidotveccurr );
    
    }
  
  }
  
  printf("\nDone with transformation ... \n");
  
 
  
  for( icart = 0 ; icart < ncartselect ; icart ++)
  { 
    atomID = ( icart - icart % 3 ) / 3 ;
    
    *( arcp + icart ) = ( *( pcurr + icart ) ) / sqrt( *( mass + atomID) );
    //*( arcp + icart * ntraj + itraj ) = ( *( pcurr + icart ) ) ;
    
    *( arcdq + icart ) = ( *( dqcurr + icart ) ) / sqrt( *( mass + atomID) );
    //*( arcdq + icart * ntraj + itraj ) = ( *( dqcurr + icart ) ) ;
    
    *( arcq + icart ) = *( cart_q0 + icart ) + ( *( dqcurr + icart ) ) / sqrt( *( mass + atomID) );
    //*( arcq + icart * ntraj + itraj ) = *( cart_q0 + icart ) + ( *( dqcurr + icart ) ) ;

  }

  
  
  
  







  for( itmp = 0 ; itmp < ncartselect  ; itmp ++ )
  {
    *( cart_q + itmp ) = *( arcq + itmp );

    *( cart_p + itmp ) = *( arcp + itmp );
  }


  //-o-o-o-o-o-o-o ... Em...What else? ... o-o-o-o-o-o-o-o-//
    
  
  if( debuggingMode == YES )
  {
    debug = fopen( "xi.deb" , "a+" );

    doutput( debug , 1 , nmodeselect , xi );

    fclose( debug );



    debug = fopen("xidot.deb" , "a+");
    
    doutput( debug , 1 , nmodeselect , xidot );
    
    fclose( debug );



    debug = fopen("arcq.deb", "a+");
    
    doutput( debug, 1 , ncartselect , cart_q);
    
    fclose( debug );



    debug = fopen("arcp.deb", "a+");
    
    doutput( debug , 1 , ncartselect , cart_p);
    
    fclose(debug);



    debug = fopen("energy.deb" , "a+");
    
    doutput(debug , 1 , nmodeselect , energy);
    
    fclose( debug );


    debug = fopen("kineticEnergy.deb", "a+");
    
    doutput(debug, 1 , nmodeselect , kineticEnergy );
    
    fclose(debug);


    debug = fopen("potentialEnergy.deb", "a+");
    
    doutput(debug, 1 , nmodeselect , potentialEnergy );
    
    fclose(debug);
    

    debug = fopen("angleList.deb" , "a+");
    
    doutput( debug , 1 , nmodeselect , angleList );
    
    fclose( debug );
    
    

    debug = fopen("levelassgn.deb", "a+");
      
    doutput(debug, 1 , nmodeselect , levelassgn );
    
    fclose(debug);

  }


//-o-o-o-o-o-o-o ... All right... Free-Up the memory blocks ... o-o-o-o-o-o-o-o-//



  free( energy ) ;

  free( xiveccurr ) ;

  free( xidotveccurr ) ;

  free( kineticEnergy ) ;

  free( potentialEnergy ) ;

  free( levelassgn ) ;

  free( dqcurr ) ;

  free( pcurr ) ;

  free( arcdq ) ;

  free( arcq ) ;

  free( arcp ) ;








}






