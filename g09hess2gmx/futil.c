#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


//-----------------  Macros  -------------------//

#ifndef PI
#define PI 3.1415926535
#endif

#ifndef YES
#define YES 1
#endif

#ifndef NO
#define NO 0 
#endif






//------------------------------------------------------------//












// -------------------------> Deck * fcopy <--------------------------- //

void fcopy(FILE * ctrlc, FILE * ctrlv)
{
  int c;
  
  while ( (c=fgetc(ctrlc)) != EOF  )
  
       fputc(c, ctrlv);



}


// -------------------------> Deck * fdel <--------------------------- //

void fdel( char filename[] )
{

  int rm_status=remove(  filename );
  
  if( rm_status == 0 )
  {
    printf("%s file deleted successfully.\n",filename);
  }  
  else
  {
    printf("Unable to delete the file\n");
    perror("Error");
  }



}




// -------------------------> Deck * flength <--------------------------- //

/* * * * * * * * * * * * * * * * * * * * * */
/* Feb. 21st, 2013, Zheng Ma               */
/* For a file containing randomly arranged */
/* numbers, this routine reveals how many  */
/* numbers are contained in specified file */
/* * * * * * * * * * * * * * * * * * * * * */



int flength( FILE * pfile)
{
  
  //FILE * recover; recover = pfile;
  
  int length=0;
 
  int c; 
  
  double temp;

  
  do
  {
    fscanf(pfile, "%lf", &temp);
    
    c=fgetc(pfile);
    
    if( c != EOF )
    {
      length++;
      
      //printf("\nNo. %d number ----- %16.13f\n", length, temp);
    
    }
  
  }
  while( c != EOF);
  
  //pfile=recover;

  return(length);  
  
  



}



// -------------------------> Deck * fskip <--------------------------- //

void fskip( FILE * pfile, int nskip )
{
  int idskip = 0 ;
  
  char c;
  
  while( idskip < nskip )
  { 
    c = fgetc( pfile );
    //printf("\nidskip=%d, c=  %c \n", idskip , c ); 
    //idstep++;
    
    if( c == '\n' )
    {
      idskip ++ ;
    
    }
    else if( c == EOF )
    {
      printf("\nNothing to skip ... Hit the bottom of this file already ...\n\n");
      
      exit(1);
    
    }
    else
    {
    
    }
  
  }
    
    


}



// -------------------------> Deck * fload <--------------------------- //

/* * * * * * * * * * * * * * * * * * * * * */
/* Feb. 21st, 2013, Zheng Ma               */
/* For a file containing randomly arranged */
/* numbers, this routine copies all the    */
/* numbers into a user specified array     */
/* * * * * * * * * * * * * * * * * * * * * */



void fload( FILE * pfile, double * buffer)
{
  //FILE * recover; recover = pfile;
  
  int length=0;
 
  int c; 
  
  double temp;

  
  do
  {
    fscanf(pfile, "%lf", &temp);
    
    c=fgetc(pfile);
    
    if( c != EOF )
    {
      *(buffer+length)=temp;
      
      temp=0.0000;
      
      length++;
      
      //printf("\nIt is No. %d number in this file ...\n", length);
    
    }
    

  
  }
  while( c != EOF);
  
  //pfile=recover;

  
  



}




// -------------------------> Deck * fsearch <--------------------------- //

void fsearch( FILE * pfile, char kw[]  )
{ 
       char * word_sc=malloc(99*sizeof(char));
  
       char c;
 
       int status=0;


       while(feof(pfile)!=1)
       {
            fscanf(pfile, "%s", word_sc);
    
            if(strcmp(word_sc, kw)==0)
	    {
                status=1; 
            
	//	printf("\nCurrent word : %s ; keyword : %s\n", word_sc, kw);
            
	    	break; 
            }
            
        //    printf("\nCurrent word : %s ; keyword : %s\n", word_sc, kw);
            
       }
  
       if( status == 0 ) 
       {
       
       printf("\nSpecified phrase : %s  was not found in specified file\n\n", kw);

       //exit(1);
       
       }

       
       /*
       c='s';
  
       while(c != EOF)
       { 
           c=fgetc(pfile);
    
           if(c =='\n') break;
       }
       */
  

}


// -------------------------> Deck * fsub <--------------------------- //

void fsub(FILE * ctrlc, char ctrlf, char ctrls, FILE * ctrlv)
{
  int c;
  
  while ( (c=fgetc(ctrlc)) != EOF  )
  {
       if( c == ctrlf)
          
          fputc(ctrls, ctrlv);
       
       else
          
          fputc(c, ctrlv);
       
  }
  
  




}


// -------------------------> Deck * stellblankline <--------------------------- //


int stellblank( char * line )
{
  int signal = 0 ;
  
  int icharacter ;
  
  char c ;
  
  char * cs ; cs = line ;
  
  int length = strlen( cs ) ;
  
  if( length == 0 )
  {
    signal = 0 ;
    
    printf("\nThis line is a blank line because there is no character in it ...\n\n");
  
  }
  else
  {
    for( icharacter = 0 ; icharacter < length ; icharacter ++ )
    {
      c = * ( cs + icharacter );
      
      if( c != ' ' && c != '\t' && c != '\n' )
      {
        printf("\nCharacter %c found ... This line is not a blank line ...\n\n" , c );
        
        signal = 1 ;
        
        break ;
      }
    
    }
  
    if( signal == 0 )
      printf("\nThis line is a blank line ...\n\n");
    else
      printf("\nThis line is a NON-blank line ... \n\n");
  
  }
  
  return( signal ) ;


}



// -------------------------> Deck * freadline <--------------------------- //

int freadline( char * s , int n , register FILE * iop , char symbol )
{
  register int c , cc ; 
  
  register char * cs ;
  
  int info ;
  
  int ichar = -1 ;
  
  cs = s ;
  
  while( -- n > 0 && ( c = getc( iop ) ) != EOF )//&& c != '\t' && c != ' ' )
  {
    //printf("\nCurrent character is %c ...\n" , c );
    
    ichar ++ ;
    
    if( c == '\n' )
    {
      break ;
    }
    else if( c == symbol )
    {
      if( ichar == 0 )
        * cs = c ;
      else
        * cs = ' ' ;
      
      cs ++ ;
      
      cc = c ;
      
      while( cc != '\n' && cc != EOF )
      {
        cc = getc( iop );
        
        //printf("\nSkipping character %c ...\n" , cc );
      }
      
      c = cc ;
      
      break ;
    
    }
    //else if( c == ' ' || c == '\t')
    //{}
    else
    {
      * cs = c ;
      
      cs ++ ;
      
    }
    
    //if ( ( * cs ++ = c ) == '\n' )  break ;
  }
  
  * cs = '\0' ;
  
  info = ( ( c == EOF && cs == s ) ? 0 : 1 );
  
  return( info );

}


// -------------------------> Deck * getfirst <--------------------------- //

char getfirst( char * line )
{
  int signal = 0 ;
  
  int icharacter ;
  
  char c ;
  
  char * cs ; cs = line ;
  
  int length = strlen( cs ) ;
  
  if( length == 0 )
  {
    signal = 0 ;
    
    printf("\nThis line is a blank line ...\n\n");
    
    c = '\n' ;
  
  }
  else
  {
    for( icharacter = 0 ; icharacter < length ; icharacter ++ )
    {
      c = * ( cs + icharacter );
      
      if( c != ' ' && c != '\t' && c != '\n' )
      {
        //printf("\nCharacter %c found ... This line is not a blank line ...\n\n" , c );
        
        signal = 1 ;
        
        break ;
      }
    
    }
  
    if( signal == 0 )
    {
      printf("\nThis line is a blank line ...\n\n");
    
      c = '\n' ;
    }
    else
      printf("\nThe first character is %c ... \n\n" , c );
  
  }
  
  return( c ) ;



}


// -------------------------> Deck * strwordcount <--------------------------- //

int strwordcount( char * string )
{
  char * cs ; cs = string ;
  
  int length = strlen( cs ) ;
  
  //printf("\nThis line has %d characters ...\n" , length );
  
  int c , nw , inword , ichar ;
  
  inword = NO ;
  
  nw = 0 ;
  
  
  while( ( length -- > 0 ) && ( c = * cs ++ ) != '\0' )
  {
    if( c == ' ' || c == '\n' || c == '\t' )
    
      inword = NO ;
    
    else if( inword == NO )
    {
      inword = YES ;
      
      ++ nw ;
    }

  }

  return( nw ) ;


}


// -------------------------> Deck * strpickword <--------------------------- //

void strpickword( char * string , int position , char * cache )
{
  char * cs ; cs = string ;
    
  int length = strlen( cs ) ;
    
  char * tmp_cache = calloc( ( length + 1 ) , sizeof( char ) ) ;
  
  //printf("\nThis line has %d characters ...\n" , length );
  
  int c , nw , inword , ichar , iword ;
  
  inword = NO ;
  
  iword = 0 ;
  
  ichar = 0 ;
  
  
  while( ( length -- > 0 ) && ( c = * cs ++ ) != '\0' )
  {
    if( c == ' ' || c == '\n' || c == '\t' )
    
      inword = NO ;
    
    else if( inword == NO )
    {
      inword = YES ;
      
      iword = iword + 1 ;
      
      if( iword == position )
      {
        * ( tmp_cache + ichar ) = c ;
        
        //printf("%c\t" , * ( tmp_cache + ichar ) );
        
        ichar ++ ;
      }
      
      
    }
    else
    {
      if( iword == position )
      {
        * ( tmp_cache + ichar ) = c ;
        
        //printf("%c\t" , * ( tmp_cache + ichar ) );
        
        ichar ++ ;
      }
    }

  }

  if( iword < position )
  {
    printf("\nThere are only %d words in this string ... So there is no #%d word ...\n" , iword , position );
  }
  else
  {
    * ( tmp_cache + ichar ) = '\0' ;
    
    printf("\nGrabbed word is %s\n" , tmp_cache );
    
    strcpy( cache , tmp_cache ) ;
  }






}






