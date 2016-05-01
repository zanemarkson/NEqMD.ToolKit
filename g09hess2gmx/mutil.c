#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>


//------------------------------------------------------------//

void dmatplus( int rdim, int cdim, double * left, double * right, int leftlda, int rightlda, double * res)
{
  int j,k;  

  double * temp = malloc(rdim*cdim*sizeof(double));
    
  //for(j=0;j<rdim*cdim;j++)   *(temp+j)=(*(left+j)) + (*(right+j));
  
  for(j=0;j<rdim;j++)
  {
    for(k=0;k<cdim;k++)
    {
      *(temp+j*cdim+k) = (*(left+j*leftlda+k)) + (*(right+j*rightlda+k));
      
      //printf();
    }
  }
   
  //for(k=0;k<rdim*cdim;k++) *(res+k)=*(temp+k);

  for(j=0;j<rdim;j++)
  {
    for(k=0;k<cdim;k++)
    {
      *(res+j*cdim+k) = *(temp+j*cdim+k);
    }
  }


  

  free(temp);




}

//------------------------------------------------------------//

/* * * * * * * * * * * * * * * * * * * * * * */
/* May. 8th, 2013, Zheng Ma                  */
/* This routine performs :                   */
/*           B = c * A                       */
/* where A and B are matrices and c is a     */
/* number. c is provided as parameter of     */
/* function. Here, the multiplication could  */
/* be applied on a part of the matrix A      */
/* instead of the whole matrix. The part     */
/* of A which the multiplication will apply  */
/* on is specified by rstart, cstart, rdim   */
/* and cdim. The parameter "rlda" and        */
/* "clda" are specifying the dimension of    */
/* input matrix.                             */
/* * * * * * * * * * * * * * * * * * * * * * */



void dnummat( double coeff, double * inp, int rlda, int clda, int rdim, int cdim, double * out )
{
  int j,k;

  //double * small = malloc(rdim*cdim*sizeof(double));

  //double * big = malloc(rlda*clda*sizeof(double));


  //for(j=0;j<rlda*clda;j++)  *(out+j)=*(inp+j);

  printf("\nThe multiplier coefficient is % 10.8E ...\n", coeff);
  
  for(j=0; j<rdim; j++)
  {
    for(k=0; k<cdim; k++)
    {
      //*(out+j*rlda+k)= (*(inp+j*rlda+k))*coeff;
      
      //printf("\nj=%d, k=%d, inpnumber=%16.13f\n", j, k, *(inp+j*clda+k)  );

      *(out+j*cdim+k)=(*(inp+j*clda+k))*coeff;

      
      //printf("\nj=%d, k=%d, outnumber=%16.13f\n", j, k, *(out+j*cdim+k)   );


    }



  }




}


//------------------------------------------------------------//

void double_allocate( int size, double * head)
{
  int j;

  head=malloc(size*sizeof(double));

  for( j=0; j<size; j++)  
  {
    *(head+j)=0.000000;
  }



}



//------------------------------------------------------------//


/* * * * * * * * * * * * * * * * * * * * * */
/* May. 2nd, 2013, Zheng Ma                */
/* This is a routine which serves as       */
/* an interface for LAPACK routine DSYEV   */
/*                                         */
/* DSYEV is used to do the matrix          */
/* diagnoalization, i.e. calculate the eig-*/
/* values and eigen-vectors for a real     */
/* symmetric matrix.                       */
/* * * * * * * * * * * * * * * * * * * * * */



void dsyev_f2c( int dimension, double * sym, double * eigvectors,  double * eigvalues  )
{
  int n, lda, lwork, info;

  double * dummy, * work, *w;

  char * jobz="V";

  char * uplo="U";

  int j,k,l;

  
  
  n=dimension; 
  
  lda=n;

  lwork=3*n-1;

  dummy=malloc(n*n*sizeof(double));

  work=malloc(lwork*sizeof(double));

  w=malloc(n*sizeof(double));

  for(j=0;j<n*n;j++)  
  {
    //printf("\nAllocating No. %d element\n",j+1);

    *(dummy+j)=*(sym+j);

  }


  for(j=0;j<n;j++)  *(work+j)=0.0000;
   
  for(j=0;j<n;j++)  *(w+j)=0.0000;



  dsyev_(jobz, uplo, &n, dummy, &lda, w, work, &lwork, &info);

  //for(j=0;j<n;j++) printf("\nNo. %d eigenvalue of this matrix is :  %16.13f\n", j+1, *(w+j));



  

  if(info != 0)
  {
   printf("\nMatrix Diagnolization failed... The value of info is %d...\n", info);
   exit(1);
  }
  else
  {
   for(j=0;j<n;j++) {*(eigvalues+j)=*(w+j);}

   for(j=0;j<n*n;j++) {*(eigvectors+j)=*(dummy+j);}


  }
 
  
  //free(dummy);
  //free(work);
  //free(w);
  

  printf("\nMatrix Diagonalization is Done ...\n\n");



}




//------------------------------------------------------------//


/* * * * * * * * * * * * * * * * * * * * * */
/* May. 6th, 2013, Zheng Ma                */
/* This routine is written to transform    */
/* a real triangular matrix from (U)pper   */
/* or (L)ower to a full symmetric matrix.  */
/*                                         */
/* Input parameter "uplo" specifies the    */
/* type of the triangular matrix           */
/*                                         */
/* * * * * * * * * * * * * * * * * * * * * */


void dtri2sym( char uplo, int lda, double * tri, double * sym)
{

  int j,k;


  if( (uplo == 'L') || ( uplo == 'l') )
  {
   //printf("\nThe requirement is %c, and the lda is %d ...\n\n", uplo, lda );
   
   for(k=0;k<lda;k++)
   {
     for(j=k;j<lda;j++)
     {
       *(sym+j*lda+k)=*(tri+j*(j+1)/2+k);

       *(sym+k*lda+j)=*(sym+j*lda+k);

     }

   }



  }
  else if ( (uplo == 'U') || ( uplo == 'u') )
  {
   //printf("\nThe requirement is %c, and the lda is %d ...\n\n", uplo, lda );

   for(k=0;k<lda;k++)
   {
     for(j=0;j<=k;j++)
     {
       *(sym+j*lda+k)=*(tri+(2*lda-j+1)*j/2+(k-j));

       *(sym+k*lda+j)=*(sym+j*lda+k);

     }


   }


  }
  

  else
  {
     printf("\nPlease specify the type of triangular matrix: (U)pper or (L)ower ...\n");
  } 





}





//------------------------------------------------------------//




/* * * * * * * * * * * * * * * * * * * * * */
/* May. 6th, 2013, Zheng Ma                */
/* This routine is written to transform    */
/* a cmplx triangular matrix from (U)pper  */
/* or (L)ower to a full symmetric matrix.  */
/*                                         */
/* Input parameter "uplo" specifies the    */
/* type of the triangular matrix           */
/*                                         */
/* * * * * * * * * * * * * * * * * * * * * */


void ztri2sym( char uplo, int lda, double _Complex * tri, double _Complex * sym)
{

  int j,k;


  if( (uplo == 'L') || ( uplo == 'l') )
  {
   //printf("\nThe requirement is %c, and the lda is %d ...\n\n", uplo, lda );
   
   for(k=0;k<lda;k++)
   {
     for(j=k;j<lda;j++)
     {
       *(sym+j*lda+k)=*(tri+j*(j+1)/2+k);

       *(sym+k*lda+j)=*(sym+j*lda+k);

     }

   }



  }
  else if ( (uplo == 'U') || ( uplo == 'u') )
  {
   //printf("\nThe requirement is %c, and the lda is %d ...\n\n", uplo, lda );

   for(k=0;k<lda;k++)
   {
     for(j=0;j<=k;j++)
     {
       *(sym+j*lda+k)=*(tri+(2*lda-j+1)*j/2+(k-j));

       *(sym+k*lda+j)=*(sym+j*lda+k);

     }


   }


  }
  

  else
  {
     printf("\nPlease specify the type of triangular matrix: (U)pper or (L)ower ...\n");
  } 




}

//------------------------------------------------------------//

void dmatprint(int dimrow, int dimcol, double *p)
{

   int i,j;
   
   for (i=0;i<dimrow;i++)
   {
     //printf("\n\n");
     
     for (j=0;j<dimcol;j++)
     {
        printf("% 20.8E\t", *(p+i*dimcol+j)  );
     }
     
     printf("\n");
     
   }
   
   printf("\n\n");
}







//------------------------------------------------------------//
//------------------------------------------------------------//
//------------------------------------------------------------//
//------------------------------------------------------------//
//------------------------------------------------------------//
//------------------------------------------------------------//
//------------------------------------------------------------//


void dkron( double *M1, int dimrow1, int dimcol1, 
            double *M2, int dimrow2, int dimcol2, 
            double *R )
{
  int i,j,k,l;
  
  for (i=0;i<dimrow1;i++)
  {
    for (j=0;j<dimcol1;j++)
    {
      for (k=0;k<dimrow2;k++)
      {
        for (l=0;l<dimcol2;l++)
        {
          *(R+(i*dimrow2+k)*(dimcol1*dimcol2)+(j*dimcol2+l))=
                                           (*(M1+i*dimcol1+j)) * (*(M2+k*dimcol2+l));
        }
      
      }
    
    }
  
  }
  


}

//------------------------------------------------------------//

void dmatmul( int dimrowleft, int dimcolright, int dimcommon, 
              double *left, double *right, double *result )

{
  int i,j,k;
  
  for (k = 0; k < dimcolright; k++)
  {
    for (j = 0; j < dimcommon; j++)
    {
      for (i = 0; i < dimrowleft; i++ )
      {
        *(result+i*dimcolright+k)=(*(result+i*dimcolright+k)) + 
                                  (*(left+i*dimcommon+j)) * (*(right+j*dimcolright+k));
      
      }
      
    }
  
  }




}

//------------------------------------------------------------//

void dmatassgn(int dimrow, int dimcol,  double  *porigin, double  *ptarget)
{

  int i,j;
  
  
  
  for (i=0;i<dimrow;i++)
  {
    for (j=0;j<dimcol;j++)
    {
       *(ptarget+i*dimcol+j)=*(porigin+i*dimcol+j);
    }
    
  }


 }

//------------------------------------------------------------//

void dzeros(int dimrow, int dimcol, double *p)
{
  int i,j;
  
  for (i=0;i<dimrow;i++)
  {
    for (j=0;j<dimcol;j++)
    {
      *(p+dimcol*i+j)=0.0000;
    }
  }
  
}

//------------------------------------------------------------//

void dones(int dimrow, int dimcol, double *p)
{
  int i,j;
  
  for (i=0;i<dimrow;i++)
  {
    for (j=0;j<dimcol;j++)
    {
      *(p+dimcol*i+j)=1.0000;
    }
  }
  
}

//------------------------------------------------------------//

int double_factorial(int x)
{

  int df,j;
  

  if (x<-1)
  {
     exit(1);
  }
  
  else if ((x==-1)||(x==0))
  {
     df=1;
  }
  
  else
  {
    df=1;
    
    for (j=x;j>0;j=j-2)
    {
      df=df*j;
    }
    
   }
   
   return df;
   
}

//------------------------------------------------------------//
double dtrace(int dimofmat, double *p)
{

int i;

double result=0.000;

for(i=0;i<dimofmat;i++)
{
  result+=*(p+i*dimofmat+i);

}

return result;



}

//------------------------------------------------------------//
void dtranspose( int dimofmat, double *porigin, double *ptarget )
{

  int i,j;
  
  double * temp=malloc(dimofmat*dimofmat*sizeof(double));
  
  for (i=0;i<dimofmat;i++)
  {
    for (j=0;j<dimofmat;j++)
    {
      *(temp+dimofmat*i+j)=*(porigin+dimofmat*j+i);
    }
    
  }
  
  for (i=0;i<dimofmat;i++)
  {
    for (j=0;j<dimofmat;j++)
    {
      *(ptarget+dimofmat*i+j)=*(temp+dimofmat*i+j);
    }
    
  }
  
  
}


//------------------------------------------------------------//
void dvec2mat(int dimrow, int dimcol, double *pvec, double *pmat)
{
  int i,j,k;
  
  /*for(i=0;i<dimrow;i++)
  {
    for(j=0;j<dimcol;j++)
    {
      *(pmat+i*dimcol+j)=*(pvec+j*dimrow+i);
    }
  }*/
  
  int N=dimrow*dimcol;
  
  i=0;
  
  j=0;
  
  for(k=0;k<N;k++)
  {
    *(pmat+i*dimcol+j)=*(pvec+k);
    
    i++;
    
    if(i==dimrow)
    {
      j++;
      
      i=0;
      
    }
      
  }



}



//------------------------------------------------------------//
int factorial(int x)
{
  int f,j;

  if (x==-1)
    
    f=1;
    
  else if (x==0)
  
    f=1;
    
  else
  {
    f=1;

    for (j=x;j>0;j--)
    {
       f=f*j;
    }

  }
  
  return f;
  
}


//------------------------------------------------------------//
void zkron( double _Complex *M1, int dimrow1, int dimcol1, 
            double _Complex *M2, int dimrow2, int dimcol2, 
            double _Complex *R )
{
  int i,j,k,l;
  
  for (i=0;i<dimrow1;i++)
  {
    for (j=0;j<dimcol1;j++)
    {
      for (k=0;k<dimrow2;k++)
      {
        for (l=0;l<dimcol2;l++)
        {
          *(R+(i*dimrow2+k)*(dimcol1*dimcol2)+(j*dimcol2+l))=
                                                     (*(M1+i*dimcol1+j)) * (*(M2+k*dimcol2+l));
        }
      
      }
    
    }
  
  }
  


}

//------------------------------------------------------------//
void zmatmul( int dimrowleft, int dimcolright, int dimcommon, 
              double _Complex *left, double _Complex *right, double _Complex *result )

{
  int i,j,k;
  
  double _Complex temp;
  
  for (k = 0; k < dimcolright; k++)
  {
    for (j = 0; j < dimcommon; j++)
    {
      temp=(*(right+j*dimcolright+k));
      
      for (i = 0; i < dimrowleft; i++ )
      {
        *(result+i*dimcolright+k)=(*(result+i*dimcolright+k)) + 
                                  (*(left+i*dimcommon+j)) * temp;
      
      }
      
    }
  
  }

}
//------------------------------------------------------------//
void zmatassgn(int dimrow, int dimcol,  double _Complex *porigin, double _Complex *ptarget)
{

  int i,j;
  
  
  
  for (i=0;i<dimrow;i++)
  {
    for (j=0;j<dimcol;j++)
    {
       *(ptarget+i*dimcol+j)=*(porigin+i*dimcol+j);
    }
    
  }


 }


//------------------------------------------------------------//
void zmatprint(int dimrow, int dimcol, double _Complex *p)
{

   int i,j;
   
   for (i=0;i<dimrow;i++)
   {
     //printf("\n");
     
     for (j=0;j<dimcol;j++)
     {
        printf("%16.13f\t%16.13f\t\t", creal(*(p+i*dimcol+j)), cimag(*(p+i*dimcol+j))  );
     }
     
     printf("\n");
     
   }
   
   printf("\n\n");
}

//------------------------------------------------------------//
void zones(int dimrow, int dimcol, double _Complex *p)
{
  int i,j;
  
  for (i=0;i<dimrow;i++)
  {
    for (j=0;j<dimcol;j++)
    {
      *(p+dimcol*i+j)=1.0000+0.0000*I;
    }
  }
  
}

//------------------------------------------------------------//
double _Complex ztrace(int dimofmat, double _Complex *p)
{

int i;

double _Complex result=0.0000+0.00000*I;

for(i=0;i<dimofmat;i++)
{
  result+=*(p+i*dimofmat+i);

}

return result;



}

//------------------------------------------------------------//
void ztranspose( int dimofmat, double _Complex *porigin, double _Complex *ptarget )
{

  int i,j;
  
  double _Complex * temp=malloc(dimofmat*dimofmat*sizeof(double _Complex));
  
  for (i=0;i<dimofmat;i++)
  {
    for (j=0;j<dimofmat;j++)
    {
      *(temp+j*dimofmat+i)=*(porigin+dimofmat*i+j);
    }
    
  }
  
  for (i=0;i<dimofmat;i++)
  {
    for (j=0;j<dimofmat;j++)
    {
      *(ptarget+dimofmat*i+j)=*(temp+i*dimofmat+j);
    }
    
  }
  
  free(temp);
  
  
}

//------------------------------------------------------------//
void zvec2mat(int dimrow, int dimcol, double _Complex *pvec, double _Complex *pmat)
{
  int i,j,k;
  
  /*for(i=0;i<dimrow;i++)
  {
    for(j=0;j<dimcol;j++)
    {
      *(pmat+i*dimcol+j)=*(pvec+j*dimrow+i);
    }
  }*/
  
  int N=dimrow*dimcol;
  
  i=0;
  
  j=0;
  
  for(k=0;k<N;k++)
  {
    *(pmat+i*dimcol+j)=*(pvec+k);
    
    i++;
    
    if(i==dimrow)
    {
      j++;
      
      i=0;
      
    }
      
  }



}

//------------------------------------------------------------//
void zzeros(int dimrow, int dimcol, double _Complex *p)
{
  int i,j;
  
  for (i=0;i<dimrow;i++)
  {
    for (j=0;j<dimcol;j++)
    {
      *(p+dimcol*i+j)=0.0000+0.0000*I;
    }
  }
  
}

//------------------------------------------------------------//
void tril(int nofline, int dimofmat, double *porigin, double *ptarget)
{
  int i,j;
  
  if (nofline<0)
  {
    for (i=(-nofline);i<dimofmat;i++)
    {
      for (j=0;j<(i+nofline);j++)
      {
        *(ptarget+dimofmat*i+j)=*(porigin+dimofmat*i+j);
      }
      
    }
  }
  
  else 
  {
    for (i=0;i<(dimofmat-nofline);i++)
    {
      for (j=0;j<=(i+nofline);j++)
      {
        *(ptarget+dimofmat*i+j)=*(porigin+dimofmat*i+j);
      }
    }
    
    for (i=(dimofmat-nofline);i<dimofmat;i++)
    {
      for (j=0;j<dimofmat;j++)
      {
        *(ptarget+dimofmat*i+j)=*(porigin+dimofmat*i+j);
      }
    }
    
  }
  
  
  
}

//------------------------------------------------------------//
void triu(int nofline, int dimofmat, double *porigin, double *ptarget)
{

  int i,j;
  
  if (nofline<0)
  {
    for (i=0;i<(-nofline);i++)
    {
      for (j=0;j<dimofmat;j++)
      {
        *(ptarget+dimofmat*i+j)=*(porigin+dimofmat*i+j);
      }
    }
 
    for (i=(-nofline);i<dimofmat;i++)
    {
      for (j=(i+nofline);j<dimofmat;j++)
      {
        *(ptarget+dimofmat*i+j)=*(porigin+dimofmat*i+j);
      }
      
    }
    
   }
   
   
   else // if (nofline>=0)
   {
     for (i=0;i<(dimofmat-nofline);i++)
     {
       for (j=(i+nofline);j<dimofmat;j++)
       {
         *(ptarget+dimofmat*i+j)=*(porigin+dimofmat*i+j);
       }
     }
     
    }
    
    
}

//------------------------------------------------------------//

//------------------------------------------------------------//
void dLiouvillian(int dim, double *hilbert, double *L)
{
  //printf("\n\nGood till here\n\n");
            
  //void dtranspose( int dimofmat, double *porigin, double *ptarget );
  
  //double temp3[dim*dim][dim*dim];

  double temp1[dim*dim][dim*dim], temp2[dim*dim][dim*dim];
  
  
  
  double eye[dim][dim];
  
  int i,j;
  
  for (i=0;i<dim;i++)
  {
    for (j=0;j<dim;j++)
    {
      if (i==j)
        
        eye[i][j]=1.00;
        
      else
      
        eye[i][j]=0.00;
   
    }
    
  }
  
  dkron( &eye[0][0], dim, dim, hilbert, dim, dim, &temp1[0][0]);
  
  dkron( hilbert, dim, dim, &eye[0][0], dim, dim, &temp2[0][0]);
  
  dtranspose( dim*dim, &temp2[0][0], &temp2[0][0]);
  
  for (i=0;i<dim*dim;i++)
  {
    for (j=0;j<dim*dim;j++)
    {
      *(L+i*dim*dim+j)=temp1[i][j]-temp2[i][j];
    }
  
  }
  
  
}


//------------------------------------------------------------//
void zLiouvillian(int dim, double *hilbert, double _Complex *L)
{
  //printf("\n\nGood till here\n\n");
  
  //void dtranspose( int dimofmat, double *porigin, double *ptarget );
  
  //double temp3[dim*dim][dim*dim];

  //double temp1[dim*dim][dim*dim], temp2[dim*dim][dim*dim];
  
  int dim2=dim*dim;
  
  double * temp1=malloc(dim2*dim2*sizeof(double));
  
  double * temp2=malloc(dim2*dim2*sizeof(double));
  
  //printf("\nGood till here!!!\n\n");
  
  double eye[dim][dim];
  
  int i,j;
  
  for (i=0;i<dim;i++)
  {
    for (j=0;j<dim;j++)
    {
      if (i==j)
        
        eye[i][j]=1.00;
        
      else
      
        eye[i][j]=0.00;
   
    }
    
  }
 
  dkron( &eye[0][0], dim, dim, hilbert, dim, dim, temp1);
  
  dkron( hilbert, dim, dim, &eye[0][0], dim, dim, temp2);
  
  dtranspose( dim2, temp2, temp2);
  
  for (i=0;i<dim*dim;i++)
  {
    for (j=0;j<dim*dim;j++)
    {
      *(L+i*dim2+j)=(*(temp1+i*dim2+j))-(*(temp2+i*dim2+j))+0.000*I;
    }
  
  }
  
  
}

//------------------------------------------------------------//
void imatrixprint(int dimrow, int dimcol, int *p)
{

   int i,j;
   
   for (i=0;i<dimrow;i++)
   {
     printf("\n\n");
     
     for (j=0;j<dimcol;j++)
     {
        printf("%d\t", *(p+i*dimcol+j)  );
     }
     
     printf("\n");
     
   }
   
   printf("\n\n");
}


//------------------------------------------------------------//
void doutput(FILE *pout, int dimrow, int dimcol, double *M)
{
  int i,j;
  
  for(i=0;i<dimrow;i++)
  {
    for(j=0;j<dimcol;j++)
    {
      //printf("\nThis number is : % 10.8E\n", (*(M+i*dimcol+j))  );
      
      fprintf(pout, "% 10.8E   ", (*(M+i*dimcol+j))  ); 
  
    }
  
  fprintf(pout, "\n\n");


  }



}

//------------------------------------------------------------//

void ioutput(FILE *pout, int dimrow, int dimcol, int *M)
{
  int i,j;
  
  for(i=0;i<dimrow;i++)
  {
    for(j=0;j<dimcol;j++)
    {
      //printf("\nThis number is : % 10.8E\n", (*(M+i*dimcol+j))  );
      
      fprintf(pout, "% d   ", (*(M+i*dimcol+j))  ); 
  
    }
  
  fprintf(pout, "\n\n");


  }



}

//------------------------------------------------------------//

void loutput(FILE *pout, int dimrow, int dimcol, long *M)
{
  int i,j;
  
  for(i=0;i<dimrow;i++)
  {
    for(j=0;j<dimcol;j++)
    {
      //printf("\nThis number is : % 10.8E\n", (*(M+i*dimcol+j))  );
      
      fprintf(pout, "% ld   ", (*(M+i*dimcol+j))  ); 
  
    }
  
  fprintf(pout, "\n\n");


  }



}

//------------------------------------------------------------//


void zoutput(FILE *pout, int dimrow, int dimcol, double _Complex *M)
{
  int i,j;
  
  for(i=0;i<dimrow;i++)
  {
    for(j=0;j<dimcol;j++)
    {
      fprintf(pout, "%lf+I*%lf\t", creal(*(M+i*dimcol+j)), cimag(*(M+i*dimcol+j))  ); 
  
    }
  
  fprintf(pout, "\n\n");


  }



}

//-----------------------------------------------------------------------------------------------//


void dzmatmul( int dimrowleft, int dimcolright, int dimcommon, 
               double  *left, double _Complex *right, double _Complex *result )

{
  int i,j,k;
  
  double _Complex temp;
  
  for (k = 0; k < dimcolright; k++)
  {
    for (j = 0; j < dimcommon; j++)
    {
      temp=(*(right+j*dimcolright+k));
      
      for (i = 0; i < dimrowleft; i++ )
      {
        *(result+i*dimcolright+k)=(*(result+i*dimcolright+k)) + 
                                  (*(left+i*dimcommon+j)) * temp;
      
      }
      
    }
  
  }

}




//-----------------------------------------------------------------------------------------------//




void zdmatmul( int dimrowleft, int dimcolright, int dimcommon, 
               double _Complex *left, double *right, double _Complex *result )

{
  int i,j,k;
  
  double _Complex temp;
  
  for (k = 0; k < dimcolright; k++)
  {
    for (j = 0; j < dimcommon; j++)
    {
      temp=(*(right+j*dimcolright+k));
      
      for (i = 0; i < dimrowleft; i++ )
      {
        *(result+i*dimcolright+k)=(*(result+i*dimcolright+k)) + 
                                  (*(left+i*dimcommon+j)) * temp;
      
      }
      
    }
  
  }

}

//-----------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------//
void mark(int choice)
{
  if(choice==2)
  {
    printf("\nGood till here!!!\n\n");
  }
  
  else if (choice==1)
  {
    printf("\nTerminating...\n\n");
    
    exit(1);
    
  }


}

//-----------------------------------------------------------------------------------------------//

void izeros(int dimrow, int dimcol, int *p)
{
  int i,j;
  
  for (i=0;i<dimrow;i++)
  {
    for (j=0;j<dimcol;j++)
    {
      *(p+dimcol*i+j)=0;
    }
  }
  
}

//------------------------------------------------------------//

void iones(int dimrow, int dimcol, int *p)
{
  int i,j;
  
  for (i=0;i<dimrow;i++)
  {
    for (j=0;j<dimcol;j++)
    {
      *(p+dimcol*i+j)=1;
    }
  }
  
}

//------------------------------------------------------------//

double dmax( int length , double * sample )
{
  int j;
  
  double tmp = *sample ;
  
  for(j=1; j<length; j++)
  {
    if( tmp <= (*(sample + j)) )
     
      tmp = *(sample + j);  
  
  }
  
  return(tmp);


}


//------------------------------------------------------------//




double dmin( int length , double * sample )
{
  int j;
  
  double tmp = *sample ;
  
  for(j=1; j<length; j++)
  {
    if( tmp >= (*(sample + j)) )
     
      tmp = *(sample + j);  
  
  }
  
  return(tmp);


}



//------------------------------------------------------------//

int imax( int length , int * sample )
{
  int j;
  
  int tmp = *sample ;
  
  for(j=1; j<length; j++)
  {
    if( tmp <= (*(sample + j)) )
     
      tmp = *(sample + j);  
  
  }
  
  return(tmp);


}


//------------------------------------------------------------//




int imin( int length , int * sample )
{
  int j;
  
  int tmp = *sample ;
  
  for(j=1; j<length; j++)
  {
    if( tmp >= (*(sample + j)) )
     
      tmp = *(sample + j);  
  
  }
  
  return(tmp);


}



//------------------------------------------------------------//

#define IA 16807 

#define IM 2147483647 

#define AM (1.0/IM) 

#define IQ 127773 

#define IR 2836 

#define MASK 123459876




long ranzm( long * idum )
{
  long k;
  
  long ans;
  
  *idum ^= MASK; 
  
  k=(*idum)/IQ; 
  
  *idum=IA*(*idum-k*IQ)-IR*k; 
  
  if (*idum < 0) 
  {
    *idum += IM; 
  }
  
  
  ans=1*(*idum);

  *idum ^= MASK;

  return ans; 


}







//------------------------------------------------------------//






