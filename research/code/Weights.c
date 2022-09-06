/*
--------------------------------------------------------------------------------------------------
This program calculates all the possible weights of Gorenstein Fano weighted projective space with
at worst terminal, or canonical, singularities. The user specifies the dimension and the type of
singularities permitted.
--------------------------------------------------------------------------------------------------
First programmed by Alexander M Kasprzyk, Sept. 2005.
Modified July, 2007.
    http://www.math.unb.ca/~kasprzyk/
--------------------------------------------------------------------------------------------------
*/

/* standard libraries */
#include <stdio.h>
#include <stdlib.h>

/* application constants */
#define	kTrue			1		/* useful truth values */
#define	kFalse			0
#define	kRuleOff		"-------------------------------------------------\n"

/* structure definitions */
typedef struct
{
	long		a, b;			/* a fractional number */
} fnum;

typedef struct
{
	FILE		*file;			/* the output file (if any) */
	char		terminal;		/* are we restricting ourselves to terminal singularities or not? */
	long		n;				/* the dimension we're working in */
	long		*k;				/* the array of k-values */
} WeightRec, *WeightPtr;

/* function prototypes */
int					main							( void );
static FILE	 *		doAppInit						( long *, char * );
static void			doAddWeightToList				( WeightPtr, long );
static FILE *		doCreateOutputFile				( long, char );
static void			doCloseOutputFile				( FILE * );
static WeightPtr	doNewWeight						( long, char );
static void			doDisposeWeight					( WeightPtr );
static void			doCalculatek					( WeightPtr, long );
static long			doTightenUpper					( WeightPtr, long, long );
static void			doFinishCalculationCanonical	( WeightPtr );
static void			doFinishCalculationTerminal		( WeightPtr );
static void			doCalculate3Lambdas				( WeightPtr, long, long );
static char			doCheckSumValid					( WeightPtr, long );
static void			doLowestForm					( fnum * );
static long			doFindHCF						( long, long );
static long			doFindLCM						( long, long );
static long			doNextPrime						( long );

/* main -	the program entry/exit point */
int main( void )
{
	WeightPtr		w = kFalse;
	long			n;
	char			terminal;
	FILE			*file;
	
	/* initialize the application */
	file = doAppInit( &n, &terminal );
	
	/* allocate the memory */
	if( !(w = doNewWeight( n, terminal )) )
		return( kFalse );
	w->file = file;
	
	/* start the calculation */
	doCalculatek( w, w->n );
	
	/* free the memory */
	doDisposeWeight( w );
	
	/* finished */
	printf( "\nFinished.\n" );
	
	return( kTrue );
}

/* doAppInit -	call to initialize the application */
static FILE *doAppInit( long *n, char *terminal )
{	
	int		temp;
	char		str[60];
	FILE		*file = kFalse;
	
	/* display some introductorary verbiage */
	printf( "%s\nThis program calculates all the possible weights\nof Gorenstein Fano weighted projective space with\n", kRuleOff );
	printf( "at worst terminal, or canonical, singularities.\nThe user specifies the dimension and the type of\nsingularities permitted.\n\n" );
	printf( "%s\nProgrammed by Alexander M Kasprzyk, Sept. 2005.\nLast modified July, 2007.\n", kRuleOff );
	printf( "\thttp://www.math.unb.ca/~kasprzyk/\n\n%s\n", kRuleOff );
	
	/* obtain the values from the user */
	/* first the dimension */
	printf( "Dimension = " );
	scanf( "%d", &temp );
	*n = (long)temp;
	
	/* now whether to limit to terminal singularities or not */
	printf( "Terminal singularities only? (y/n) " );
	scanf( "%s", str );
	if( str[0] == 'y' )	*terminal = kTrue;
	else				*terminal = kFalse;
	
	/* shall we save the output to a LaTeX file or not? */
	printf( "Output results to LaTeX file? (y/n) " );
	scanf( "%s", str );
	if( str[0] == 'y' )	file = doCreateOutputFile( *n, *terminal );
	
	/* leave a blank line */
	printf( "\n" );
	
	return( file );
}

/* doAddWeightToList -call to add the given weight to the list we're maintaining */
static void doAddWeightToList( WeightPtr w, long h )
{
	long		i;
	
	/* output the \lambda_i to the screen */
	for( i = 0; i <= w->n; i++ )		printf( "%d, ", h / w->k[i] );
	printf( "h=%d\n", h );
	
	/* output the \lambda_i to the file */
	if( w->file )
	{
		for( i = 0; i <= w->n; i++ )	fprintf( w->file, "$%d$&", h / w->k[i] );
		fprintf( w->file, "$%d$\\\\\n", h );
	}
}

/* doCreateOutputFile -	call to create a LaTeX output file, and write the headers */
static FILE *doCreateOutputFile( long n, char terminal )
{
	char		name[25];
	FILE		*file;
	
	/* create the file name */
	if( terminal )	sprintf( name, "dim_%d_term.tex", n );
	else			sprintf( name, "dim_%d_canon.tex", n );
	
	/* create the file */
	if( !(file = fopen( name, "w" )) )
		printf( "Error! Unable to create file '%s'.\n", name );
	else
	{	/* write the file headers */
		long		i;
		
		if( terminal )	fprintf( file, "$n=%d$; terminal singularities.\n\n\\begin{tabular}{|", n );
		else			fprintf( file, "$n=%d$; canonical singularities.\n\n\\begin{tabular}{|", n );
		for( i = 0; i <= n; i++ )	fprintf( file, "c" );
		fprintf( file, "|c|}\n\\hline\n" );
		for( i = 0; i <= n; i++ )	fprintf( file, "$\\lambda_%d$&", i );
		fprintf( file, "$h$\\\\\n\\hline\n" );
	}
	
	return( file );
}

/* doCloseOutputFile -	call to close the output LaTeX file */
static void doCloseOutputFile( FILE *file )
{
	/* write the footers */
	fprintf( file, "\\hline\n\\end{tabular}\n" );
	
	/* finally, close the file */
	fclose( file );
}

/* doNewWeight -	call to allocate the memory for a new weight of dimension n */
static WeightPtr doNewWeight( long n, char terminal )
{
	WeightPtr	w = kFalse;
	long		i;
	
	/* check the dimension is valid */
	if( n < 2 )
	{
		printf( "The dimension must be greater than 1! (You have entered %d.)\n", n );
		return( kFalse );
	}
	else if( (n < 3) && terminal )
	{
		printf( "In the terminal case, the dimension must be greater than 2! (You have entered 2.)\n" );
		return( kFalse );
	}
	
	/* allocate the memory for the weight */
	if( !(w = (WeightPtr)malloc( sizeof( WeightRec ) )) )
	{
		printf( "Not enough memory to create the weight structure!\n" );
		return( kFalse );
	}
	
	/* allocate the memory for the k-value array */
	if( !(w->k = (long *)malloc( sizeof( long ) * (n + 1) )) )
	{
		free( (void *)w );
		printf( "Not enough memory to create a k array!\n" );
		return( kFalse );
	}
	
	/* set the dimension and singularity type, and set the file pointer to 0 */
	w->terminal = terminal;
	w->n = n;
	w->file = kFalse;
	
	/* if we're in the terminal case, add P^n by hand, since we use sharp bounds which assume that \lambda_n > 1. */
	if( terminal )
	{
		for( i = 0; i <= n; i++ )	w->k[i] = n + 1;
		doAddWeightToList( w, n + 1 );
	}
	
	/* finally, zero the array; we're ready to begin */
	for( i = 0; i <= n; i++ )	w->k[i] = 0;
	
	return( w );
}
	
/* doDisposeWeight -	call to free the memory of a weight */
static void doDisposeWeight( WeightPtr w )
{
	/* check that we actually have a weight before we try to dispose of it */
	if( w )
	{
		/* dispose of the k-array */
		if( w->k )		free( (void *)w->k );
		
		/* close any open files */
		if( w->file )	doCloseOutputFile( w->file );
		
		/* dispose of the weight */
		free( (void *)w );
	}
}

/* doCalculatek -	calculates recursively the possible k-values */
static void doCalculatek( WeightPtr w, long i )
{
	long		lower, upper, j;
	fnum		sum = {0,1}, temp;
	
	/* first check that we're not already calculated enough */
	if( w->terminal && (i <= 2) )
	{	/* terminal bound */
		doFinishCalculationTerminal( w );
		return;
	}
	else if( !w->terminal && !i )
	{	/* canonical bound */
		doFinishCalculationCanonical( w );
		return;
	}
	
	/* otherwise, step through the possible values of k_i */
	/* first establish the lower bounds */
	lower = w->n-i+2;
	if( w->terminal )	lower++;
	if( (i < w->n) && (lower < w->k[i+1]) )	lower = w->k[i+1];
	
	/* now the upper bounds */
	for( j = i + 1; j <= w->n; j++ )
	{
		sum.a = sum.a * w->k[j] + sum.b;
		sum.b = sum.b * w->k[j];
		doLowestForm( &sum );
	}
	temp.a = (i + 1) * sum.b;
	temp.b = sum.b - sum.a;
	doLowestForm( &temp );
	upper = temp.a / temp.b;
			
	/* is this range sensible? */
	if( upper < lower )	return;
	
	/* now check whether the upper bound can be bettered by direct calculation */
	if( i < w-> n )
		upper = doTightenUpper( w, i, upper );
	else if( w->terminal )	/* \lambda_n <= n-1 in the terminal case */
		upper = w->n - 1;
	
	/* is the range still sensible? */
	if( upper < lower )	return;
	
	/* iterate on this range */
	for( j = lower; j <= upper; j++ )
	{
		w->k[i] = j;
		doCalculatek( w, i - 1 );
	}
}

/* doTightenUpper -	call to check whether the upper bound can be tightened by direct calculation of the sum */
static long doTightenUpper( WeightPtr w, long i, long upper )
{
	long		kappa, sum = 1;
	
	/* work through the \kappa */
	for( kappa = 2; kappa <= upper - 1; kappa++ )
	{
		long		j;
		
		/* calculate the sum as it currently stands */
		sum++;
		for( j = i + 1; j <= w->n; j++ )
			if( !(kappa % w->k[j]) )	sum--;
		
		/* check that the values are valid */
		if( w->terminal )
		{	/* the terminal case */
			if( (kappa < upper - 2) && (sum > w->n - 1) )			upper = kappa + 2;
			if( (kappa < upper - 3) && (sum > w->n - 2) )			upper = kappa + 3;
			if( (sum < 2) && (kappa >= 2) && (kappa < upper - 2) )	upper = kappa + 2;
		}
		else
		{	/* the canonical case */
			if( (kappa < upper - 2) && ((sum > w->n) || !sum) )		upper = kappa + 2;
		}
	}
	
	/* return the bound we found */
	return( upper );	
}

/* doFinishCalculationCanonical -	call to check the weight is valid and output the k-values; canonical case */
static void doFinishCalculationCanonical( WeightPtr w )
{
	long		h = 1, i, total = 0, lam;
	
	/* calculate the value of h */
	for( i = 1; i <= w->n; i++ )	h = doFindLCM( w->k[i], h );
	
	/* calculate the sum \lambda_1 + ... + \lambda_n, in order to find \lambda_0 */
	for( i = 1; i <= w->n; i++ )	total += h / w->k[i];
	lam = h - total;
	
	/* check whether the lambda_0 value is sensible */
	if( lam < 1 )			return;
	if( h % lam )			return;
	if( lam > h / w->k[1] )	return;
	
	/* all's appears well; we have our weights */
	w->k[0] = h / lam;
	if( doCheckSumValid( w, h ) )
		doAddWeightToList( w, h );
}

/* doFinishCalculationTerminal -	call to check the weight is valid and output the k-values; terminal case */
static void doFinishCalculationTerminal( WeightPtr w )
{
	long		h = 1, i, total = 0, lam;
	
	/* calculate the value of h */
	for( i = 3; i <= w->n; i++ )	h = doFindLCM( w->k[i], h );
	
	/* calculate the sum \lambda_3 + ... + \lambda_n, in order to find \lambda_0, \lambda_1, and \lambda_2 */
	for( i = 3; i <= w->n; i++ )	total += h / w->k[i];
	lam = h - total;
	
	/* check whether the value of \lambda_0+\lambda_1+\lambda_2 is sensible */
	if( lam < 3 )	return;
	
	/* check through the possible \lambda_i, i = 0,1,2, outputting any possibilities we find */
	doCalculate3Lambdas( w, h, lam );
}

/* doCalculate3Lambdas -	call to check through the possible values of \Lambda_0,\lambda_1,\lambda_2 in the terminal case */
static void doCalculate3Lambdas( WeightPtr w, long h, long lam )
{
	long		lambda2, max2;
	
	/* establish upper bounds on \lambda_2 */
	max2 = lam - 2;
	if( max2 > h / w->k[3] )	max2 = h / w->k[3];
	
	/* now check each possible value of \lambda_2 in turn */
	for( lambda2 = max2; lambda2 >= lam / 3; lambda2-- )
		if( !(h % lambda2) )	/* we require \lambda_2 | h */
		{
			long		lambda1, max1;
			
			/* estabilsh upper bounds on \lambda_1 */
			max1 = lam - lambda2 - 1;
			if( max1 > lambda2 )		max1 = lambda2;
			
			/* now check each possible value of \lambda_1 in turn */
			for( lambda1 = max1; lambda1 >= (lam - lambda2) / 2; lambda1-- )
				if( !(h % lambda1) )	/* we require \lambda_1 | h */
				{
					long		lambda0 = lam - lambda2 - lambda1;
					
					/* check that the value of \lambda_0 is valid */
					if( (lambda0 <= lambda1) && (lambda0 > 0) && !(h % lambda0) )
					{
						/* calculate k_0, k_1, and k_2 */
						w->k[0] = h / lambda0;
						w->k[1] = h / lambda1;
						w->k[2] = h / lambda2;
						
						/* check that the k-values are valid */
						if( (w->k[2] > w->n) && (w->k[1] > w->n + 1) && (w->k[0] > w->n + 1) && doCheckSumValid( w, h ) )							
							doAddWeightToList( w, h );	/* we've found a possible weight */
					}
				}
		}
}

/* doCheckSumValid -	call to check that s(\kappa) and \Sigma(\kappa) are behaving as we expect */
static char doCheckSumValid( WeightPtr w, long h )
{
	long		kappa, sum = 1;
	
	/* step through all \kappa\in {1, ... ,h-1} calculating the s value and sum */
	for( kappa = 2; kappa < h; kappa++ )
	{
		long		s = 0, i;
		
		/* calculate s(\kappa) */
		for( i = 0; i <= w->n; i++ )
			if( !(kappa % w->k[i]) )	s++;
		
		/* calculate \Sigma(\kappa) */
		sum = sum + 1 - s;
		
		/* check that the values are valid */
		if( w->terminal )
		{	/* the terminal case */
			if( s > w->n - 3 )									return( kFalse );
			if( (kappa == 2) && (s > 0) )						return( kFalse );
			if( (kappa <= h - 3) && (sum > w->n - 2) )			return( kFalse );
			if( (kappa == h - 2) && (sum > w->n - 1) )			return( kFalse );
			if( (sum < 2) && (kappa >= 2) && (kappa <= h - 2) )	return( kFalse );
		}
		else
		{	/* the canonical case */
			if( s > w->n - 1 )									return( kFalse );
			if( (kappa <= h - 2) && ((sum > w->n) || !sum) )	return( kFalse );
		}
		
	}
	
	return( kTrue );
}

/* doLowestForm -	call to put the fraction in its lowest form */
static void doLowestForm( fnum *fr )
{
	long			test = 2, neg = 1;
	
	/* get the -'ve tidy */
	if( fr->b < 0 )			{ fr->b *= -1; neg *= -1; }
	if( fr->a < 0 )			{ fr->a *= -1; neg *= -1; }
	
	/* clear out the obvious possibilities */
	if( !fr->a )			{ fr->b = 1; return; }
	if( !(fr->a % fr->b) )	{ fr->a /= (neg * fr->b); fr->b = 1; return; }
		
	/* start clearing the cofactors */
	while( (test <= fr->a) && (test <= fr->b) )
	{
		if( !(fr->a % test) && !(fr->b % test) )	{ fr->a /= test; fr->b /= test; }
		else	test = doNextPrime( test );
	}
	fr->a *= neg;
}

/* doFindHCF -	call to calculate the HCF of two integers */
static long doFindHCF( long a, long b )
{
	long	found = 1, test = 2;
	
	/* get the -'ve tidy */
	if( a < 0 )			a*= -1;
	if( b < 0 )			b *= -1;
	
	/* do some obvious checks */
	if( !a )			return( b );
	if( !b )			return( a );
	if( !(a % b) )		return( b );
	if( !(b % a) )		return( a );
	
	/* start checking the primes */
	while( (test < a) && (test < b) )
	{
		if( !(a % test) && !(b % test) )		{ a /= test; b /= test; found *= test; }
		else	test = doNextPrime( test );
	}
	
	return( found );
}

/* doFindLCM -	call to calculate the LCM of two integers */
static long doFindLCM( long a, long b )
{
	return( (a * b) / doFindHCF( a, b ) );
}

/* doNextPrime -	call to get the next prime (or else the next odd numer) */
/*	Prime number list from "The Prime Pages", http://primes.utm.edu/	*/
static long doNextPrime( long p )
{
	switch( p )
	{
		case 2:		p = 3;		break;	case 3:		p = 5;		break;	case 5:		p = 7;		break;
		case 7:		p = 11;		break;	case 11:	p = 13;		break;	case 13:	p = 17;		break;
		case 17:	p = 19;		break;	case 19:	p = 23;		break;	case 23:	p = 29;		break;
		case 29:	p = 31;		break;	case 31:	p = 37;		break;	case 37:	p = 41;		break;
		case 41:	p = 43;		break;	case 43:	p = 47;		break;	case 47:	p = 53;		break;
		case 53:	p = 59;		break;	case 59:	p = 61;		break;	case 61:	p = 67;		break;
		case 67:	p = 71;		break;	case 71:	p = 73;		break;	case 73:	p = 79;		break;
		case 79:	p = 83;		break;	case 83:	p = 89;		break;	case 89:	p = 97;		break;
		case 97:	p = 101;	break;	case 101:	p = 103;	break;	case 103:	p = 107;	break;
		case 107:	p = 109;	break;	case 109:	p = 113;	break;	case 113:	p = 127;	break;
		case 127:	p = 131;	break;	case 131:	p = 137;	break;	case 137:	p = 139;	break;
		case 139:	p = 149;	break;	case 149:	p = 151;	break;	case 151:	p = 157;	break;
		case 157:	p = 163;	break;	case 163:	p = 167;	break;	case 167:	p = 173;	break;
		case 173:	p = 179;	break;	case 179:	p = 181;	break;	case 181:	p = 191;	break;
		case 191:	p = 193;	break;	case 193:	p = 197;	break;	case 197:	p = 199;	break;
		case 199:	p = 211;	break;	case 211:	p = 223;	break;	case 223:	p = 227;	break;
		case 227:	p = 229;	break;	case 229:	p = 233;	break;	case 233:	p = 239;	break;
		case 239:	p = 241;	break;	case 241:	p = 251;	break;	case 251:	p = 257;	break;
		case 257:	p = 263;	break;	case 263:	p = 269;	break;	case 269:	p = 271;	break;
		case 271:	p = 277;	break;	case 277:	p = 281;	break;	case 281:	p = 283;	break;
		case 283:	p = 293;	break;	case 293:	p = 307;	break;	case 307:	p = 311;	break;
		case 311:	p = 313;	break;	case 313:	p = 317;	break;	case 317:	p = 331;	break;
		case 331:	p = 337;	break;	case 337:	p = 347;	break;	case 347:	p = 349;	break;
		case 349:	p = 353;	break;	case 353:	p = 359;	break;	case 359:	p = 367;	break;
		case 367:	p = 373;	break;	case 373:	p = 379;	break;	case 379:	p = 383;	break;
		case 383:	p = 389;	break;	case 389:	p = 397;	break;	case 397:	p = 401;	break;
		case 401:	p = 409;	break;	case 409:	p = 419;	break;	case 419:	p = 421;	break;
		case 421:	p = 431;	break;	case 431:	p = 433;	break;	case 433:	p = 439;	break;
		case 439:	p = 443;	break;	case 443:	p = 449;	break;	case 449:	p = 457;	break;
		case 457:	p = 461;	break;	case 461:	p = 463;	break;	case 463:	p = 467;	break;
		case 467:	p = 479;	break;	case 479:	p = 487;	break;	case 487:	p = 491;	break;
		case 491:	p = 499;	break;	case 499:	p = 503;	break;	case 503:	p = 509;	break;
		case 509:	p = 521;	break;	case 521:	p = 523;	break;	case 523:	p = 541;	break;
		case 541:	p = 547;	break;	case 547:	p = 557;	break;	case 557:	p = 563;	break;
		case 563:	p = 569;	break;	case 569:	p =	571;	break;	case 571:	p = 577;	break;
		case 577:	p = 587;	break;	case 587:	p = 593;	break;	case 593:	p = 599;	break;
		case 599:	p = 601;	break;	case 601:	p = 607;	break;	case 607:	p = 613;	break;
		case 613:	p = 617;	break;	case 617:	p = 619;	break;	case 619:	p = 631;	break;
		case 631:	p = 641;	break;	case 641:	p = 643;	break;	case 643:	p = 647;	break;
		case 647:	p = 653;	break;	case 653:	p = 659;	break;	case 659:	p = 661;	break;
		case 661:	p = 673;	break;	case 673:	p = 677;	break;	case 677:	p = 683;	break;
		case 683:	p = 691;	break;	case 691:	p = 701;	break;	case 701:	p = 709;	break;
		case 709:	p = 719;	break;	case 719:	p = 727;	break;	case 727:	p = 733;	break;
		case 733:	p = 739;	break;	case 739:	p = 743;	break;	case 743:	p = 751;	break;
		case 751:	p = 757;	break;	case 757:	p = 761;	break;	case 761:	p = 769;	break;
		case 769:	p = 773;	break;	case 773:	p = 787;	break;	case 787:	p = 797;	break;
		case 797:	p = 809;	break;	case 809:	p = 811;	break;	case 811:	p = 821;	break;
		case 821:	p = 823;	break;	case 823:	p = 827;	break;	case 827:	p = 829;	break;
		case 829:	p = 839;	break;	case 839:	p = 853;	break;	case 853:	p = 857;	break;
		case 857:	p = 859;	break;	case 859:	p = 863;	break;	case 863:	p = 877;	break;
		default:	p += 2;		break;
	}
	
	return( p );
}
