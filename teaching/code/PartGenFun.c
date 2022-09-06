/*
------------------------------------------------------------------------------------------
By Alexander M Kasprzyk, January 2003
http://www.math.unb.ca/~kasprzyk/
------------------------------------------------------------------------------------------
This program lists the coefficients of x^i in the expansion of
	(1 + x)(1 + x^2)...(1 + x^k)
for k = 1,..., kval, where 'kval' is defined by the user.
------------------------------------------------------------------------------------------
*/

/* standard libraries */
#include <stdio.h>
#include <stdlib.h>

/* useful constants */
#define	kRuleOff		"---------------------------------------------\n\n"
#define	kTrue		1
#define	kFalse		0

/* function prototypes */
int						main				( void );
static unsigned long	doAppInit			( void );
static char				doAllocateMemory	( unsigned long, double **, double ** );
static FILE *			doCreateFile		( unsigned long );
static void				doCalculate			( unsigned long, double *, double *, FILE * );

/* main -	the program entry/exit point */
int main( void )
{
	unsigned long	kval;
	
	/* input the range over which to calcuate the coefficients */
	kval = doAppInit();
	if( kval > 1 )
	{
		double		*oldArray, *newArray;
		
		/* allocate the memory */
		if( doAllocateMemory( kval, &oldArray, &newArray ) )
		{
			FILE		*file;
			
			/* create the LaTeX output file (if requested) */
			file = doCreateFile( kval );
			
			/* calculate the coefficients */
			doCalculate( kval, oldArray, newArray, file );
			
			/* close the LaTeX output file */
			if( file )	fclose( file );
			
			/* free the memory */
			free( (void *)oldArray );
			free( (void *)newArray );
		}
	}
}	

/* doAppInit -	call to initialize the application */
static unsigned long doAppInit( void )
{
	int	kval;
	
	/* display the introductary text */
	printf( "%sProgrammed by Alexander M Kasprzyk, Jan 2003.\n\n", kRuleOff );
	printf( "\thttp://www.math.unb.ca/~kasprzyk/\n\n%s", kRuleOff );
	printf( "This program lists the coefficients of x^i in the expansion of\n\t(1 + x)(1 + x^2)...(1 + x^k)" );
	
	/* input the range over which to calculate the results */
	printf( "\n\n%sList coefficients for 1 <= k <= ", kRuleOff );
	scanf( "%d", &kval );
	
	return( kval );
}

/* doAllocateMemory -	call to assign the memory needed for the search */
static char doAllocateMemory( unsigned long kval, double **oldArray, double **newArray )
{
	unsigned long	deln = kval * (kval + 1) / 2 + 1;
	
	/* allocate the memory */
	if( !(*oldArray = (double *)malloc( sizeof( double ) * deln )) )
	{
		printf( "Not enough memory!\n" );
		return( kFalse );
	}
	if( !(*newArray = (double *)malloc( sizeof( double ) * deln )) )
	{
		free( (void *)(*oldArray) );
		printf( "Not enough memory!\n" );
		return( kFalse );
	}
	
	return( kTrue );
}

/* doCreateFile -	call to create the LaTeX output file (if required) */
static FILE *doCreateFile( unsigned long kval )
{
	char			latex[5];
	FILE			*file = kFalse;
			
	/* ask whether to save the results to a file or not */
	printf("\nDo you wish to output these results to a LaTeX file? (y/n) " );
	scanf( "%s", latex );
	
	/* if requested, create the file */
	if( (latex[0] == 'y') || (latex[0] == 'Y') )
	{
		char		name[50];
		
		sprintf( name, "Coeff_%d.tex", kval );
		if( file = fopen( name, "w" ) )
		{
			fprintf( file, "The coefficients of $x^i$ in $\\prod_{i=1}^k(1+x^i)$ for $k=1,\\ldots,%d$.\n", kval );
			fprintf( file, "\\begin{longtable}{|r|p{6in}|}\n\\hline$k$&Coefficients\\\\\n\\hline\\endhead\n" );
		}
		else
			printf( "\nUnable to create the file '%s'!\n", name );
	}
	printf( "\n%s", kRuleOff );
	
	/* return the file (if any) */
	return( file );
}

/* doCalculate -	call to calculate the coefficients */
static void doCalculate( unsigned long kval, double *oldArray, double *newArray, FILE *file )
{
	unsigned long	deln = kval * (kval + 1) / 2 + 1, count;
	
	/* set the case of n=1 by hand */
	for( count = 0; count < deln; count++ )
		*(newArray + count) = 0;
	*newArray = *(newArray + 1) = *oldArray = *(oldArray + 1) = 1;
	deln = 2;
	printf( "n =\t1\n1\t1\t" );
	if( file )	fprintf( file, "1&1,1\\\\\n\\hline\n" );
	
	/* now calculate the remaining cases inductively */
	for( count = 2; count <= kval; count++ )
	{
		unsigned long		counter;
		
		printf( "\n\nn =\t%d\n", count );
		if( file )	fprintf( file, "%d&", count );
		deln += count;
		
		printf( "%.0f\t", *newArray );
		if( file )	fprintf( file, "%.0f", *newArray );
		
		for( counter = 1; counter < count; counter++ )
		{
			printf( "%.0f\t", *(newArray + counter) );
			if( file )	fprintf( file, ", %.0f", *(newArray + counter) );
		}
			
		for( counter = count; counter < deln; counter++ )
		{
			*(oldArray + counter) = *(newArray + counter);
			*(newArray + counter) += *(oldArray + counter - count);
			printf( "%.0f\t", *(newArray + counter) );
			if( file )	fprintf( file, ", %.0f", *(newArray + counter) );
		}
		*(oldArray + count) = *(newArray + count);
		if( file )	fprintf( file, "\\\\\n\\hline\n" );
	}
	
	/* finish off */
	printf( "\n\n%sFinished.\n", kRuleOff );
	if( file )	fprintf( file, "\\end{longtable}\n", kval );
}
