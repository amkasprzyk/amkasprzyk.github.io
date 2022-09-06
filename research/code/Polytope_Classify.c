/*
----------------------------------------------------------------------------------------------------------
Programmed by Alexander M Kasprzyk, May 2003.
http://www.math.unb.ca/~kasprzyk/
----------------------------------------------------------------------------------------------------------
For an explanation of the mathematics, see:
Toric Fano threefolds with terminal singularities, Tohoku Mathematical Journal, 58 (2006), no. 1, 101-121.
----------------------------------------------------------------------------------------------------------
*/

/* standard libraries */
#include <stdio.h>
#include <stdlib.h>

/* constants */
#define	kTrue					1		/* handy truth values */
#define	kFalse					0

#define	kNoError				0		/* no error */
#define	kMemError				1		/* not enough memory */

#define	kNumMin					13		/* the number of minimal polytopes */
#define	kRuleOff				"---------------------------------------------\n\n"

/* macro functions */
#define	MPointsEqual( pt1, pt2 )		(((pt1).x == (pt2).x) && ((pt1).y == (pt2).y) && ((pt1).z == (pt2).z))
#define	MPointsNotEqual( pt1, pt2 )		(((pt1).x != (pt2).x) || ((pt1).y != (pt2).y) || ((pt1).z != (pt2).z))
#define	MNormal( pt1, pt2, nor )		(nor).x = (pt1).y * (pt2).z - (pt1).z * (pt2).y; (nor).y = (pt1).z * (pt2).x - (pt1).x * (pt2).z; (nor).z = (pt1).x * (pt2).y - (pt1).y * (pt2).x
#define	MSubtract( pt1, pt2, res )		(res).x = (pt1).x - (pt2).x; (res).y = (pt1).y - (pt2).y; (res).z = (pt1).z - (pt2).z
#define	MDot( pt1, pt2 )				((pt1).x * (pt2).x + (pt1).y * (pt2).y + (pt1).z * (pt2).z)
#define	MSetPoint( pt1, l1, l2, l3 )	(pt1).x = l1; (pt1).y = l2; (pt1).z = l3
#define	MSPt( a, b, c, d, e )			MSetPoint( p[a]->vertices[b], c, d, e )
#define	M3DTo2D( wd, ht, pt )			wd = 30.0 * (double)(pt).x + 12.0 * (double)(pt).z; ht = 30.0 * (double)(pt).y -18.0 * (double)(pt).z

/* data structures */
typedef struct
{
	short		x, y, z;				/* a 3-dimentional point */
} Point3DRec, *Point3DPtr;

typedef struct
{
	short		top, front, right,		/* a box */
				bottom, back, left;
} BoundsRec, *BoundsPtr;

typedef struct	PolyListRec	PolyListRec, *PolyListPtr, **PolyListHandle;

typedef struct
{
	short		numVertices,			/* the number of vertices */
				numChildren,			/* the number of children */
				numParents,				/* the number of parents */
				id;						/* the polytope ID (assigned at the end) */
	char		simplicial;				/* is the polytope simplicial? */
	Point3DPtr	vertices;				/* the list of vertices */
	PolyListPtr	children,				/* the list of child polytopes */
				parents;				/* the list of parent polytopes */
} PolytopeRec, *PolytopePtr, **PolytopeHandle;

struct PolyListRec
{
	PolytopePtr	p;						/* the polytope */
	PolyListPtr	next;					/* the next polytope in the list */
};

/* global variables */
PolyListPtr		gPolyList;				/* the list of found polytopes */

/* function prototypes */
int					main						( void );
static void			doAppInit					( void );
static char			doClassifyPolytopes			( void );
static void			doDisposePolytopeList		( void );
static char			doCreateMinimalPolytopes	( PolytopeHandle );
static char			doAddPolytopeToList			( PolytopePtr );
static PolytopePtr	doNewPolytope				( short );
static PolytopePtr	doNewPolytopeChild			( PolytopePtr );
static void			doDisposePolytope			( PolytopePtr );
static char			doIsSimplicial				( PolytopePtr );
static char			doAreCoplanar				( Point3DPtr, Point3DPtr, Point3DPtr, Point3DPtr );
static char			doIsFace					( PolytopePtr, Point3DPtr, Point3DPtr, Point3DPtr );
static char			doIsEdge					( PolytopePtr, Point3DPtr, Point3DPtr );
static char			doIsFreeTetrahedron			( Point3DPtr, Point3DPtr, Point3DPtr );
static void			doAddPointToBoundingBox		( BoundsPtr, Point3DPtr );
static char			doIsInternal				( Point3DPtr, Point3DPtr, Point3DPtr, Point3DPtr );
static char			doIsChildFano				( PolytopePtr, Point3DPtr );
static char			doEnlargePolytope			( PolytopePtr );
static char			doAddOverVertex				( PolytopePtr );
static char			doAddOverEdge				( PolytopePtr );
static char			doAddOverFace				( PolytopePtr );
static char			doFindNewVertex				( PolytopePtr, Point3DPtr, Point3DPtr, Point3DPtr );
static char			doCheckBarrycentric			( PolytopePtr, short, short, short, short, Point3DPtr, Point3DPtr, Point3DPtr );
static char			doCheckBarrycentricPerm		( PolytopePtr, short, short, short, short, Point3DPtr, Point3DPtr, Point3DPtr );
static char			doUpdateList				( PolytopePtr, PolyListHandle );
static void			doAddChildToList			( PolytopePtr, PolytopePtr );
static PolytopePtr	doIsNewPolytope				( PolytopePtr, char * );
static char			doArePolytopesSimilar		( PolytopePtr, PolytopePtr );
static char			doRotatePolytope			( PolytopePtr, PolytopePtr, PolytopePtr, short, short, short );
static char			doArePolytopesSame			( PolytopePtr, PolytopePtr );
static void			doAssignIDs					( void );
static void			doSaveResults				( void );
static void			doWriteVertices				( PolytopePtr, FILE * );
static void			doWriteList					( PolyListPtr, FILE * );

/* main -	the program entry/exit point */
int main( void )
{	
	/* initialize the application */
	doAppInit();
		
	/* generate the polytope list */
	if( doClassifyPolytopes() == kNoError )
	{
		/* save the results */
		doSaveResults();
		
		/* finally dispose of the polytope list */
		doDisposePolytopeList();
	}
	else
		printf( "Calculation aborted!!!\n" );
		
	/* finished */
	printf( "\n%sFinished.\n", kRuleOff );
}


/* doAppInit -	call to initialize the application */
static void doAppInit( void )
{	
	printf( "%sProgrammed by Alexander M Kasprzyk, May 2003.\n", kRuleOff );
	printf( "\thttp://www.math.unb.ca/~kasprzyk/\n\n%s", kRuleOff );
	printf( "Classification can take up to 30 minutes.\n\n%s", kRuleOff );
	
	gPolyList = kFalse;
}

/* doClassifyPolytopes -	call to classify the polytopes */
static char doClassifyPolytopes( void )
{
	PolytopePtr	p[kNumMin];
	char		err;
	short		i;
	
	/* create the seeds */
	if( err = doCreateMinimalPolytopes( p ) )
		return( err );
	
	/* grow from each seed in turn */
	for( i = 0; i < kNumMin; i++ )
	{
		printf( "Growing Minimal Polytope %d of %d...\n", i + 1, kNumMin );
		if( err = doEnlargePolytope( p[i] ) )
			return( err );
	}
	
	/* assign the polytope ID's */
	doAssignIDs();
	
	return( kNoError );
}

/* doDisposePolytopeList -	call to dispose of the polytope list */
static void doDisposePolytopeList( void )
{
	while( gPolyList )
	{
		PolyListPtr	temp = gPolyList->next;
		
		doDisposePolytope( gPolyList->p );
		free( (void *)gPolyList );
		gPolyList = temp;
	}
}

/* doCreateMinimalPolytopes -	call to create the minimal polytopes */
static char doCreateMinimalPolytopes( PolytopeHandle p )
{
	short 		i;
	char		err;
	
	/* allocate the memory for the minimal polytopes */
	if( !(p[0] = doNewPolytope( 4 )) )	return( kMemError );
	if( !(p[1] = doNewPolytope( 4 )) )	return( kMemError );
	if( !(p[2] = doNewPolytope( 4 )) )	return( kMemError );
	if( !(p[3] = doNewPolytope( 4 )) )	return( kMemError );
	if( !(p[4] = doNewPolytope( 4 )) )	return( kMemError );
	if( !(p[5] = doNewPolytope( 4 )) )	return( kMemError );
	if( !(p[6] = doNewPolytope( 4 )) )	return( kMemError );
	if( !(p[7] = doNewPolytope( 4 )) )	return( kMemError );
	
	if( !(p[8] = doNewPolytope( 5 )) )	return( kMemError );
	if( !(p[9] = doNewPolytope( 5 )) )	return( kMemError );
	if( !(p[10] = doNewPolytope( 5 )) )	return( kMemError );
	
	if( !(p[11] = doNewPolytope( 6 )) )	return( kMemError );
	if( !(p[12] = doNewPolytope( 6 )) )	return( kMemError );
	
	/* set the vertices of the polytopes */
	MSPt( 0, 0, 1, 0, 0 );	MSPt( 0, 1, 0, 1, 0 );	MSPt( 0, 2, 0, 0, 1 );	MSPt( 0, 3, -1, -1, -1 );
	MSPt( 1, 0, 1, 0, 0 );	MSPt( 1, 1, 0, 1, 0 );	MSPt( 1, 2, 1, -3, 5 );	MSPt( 1, 3, -2, 2, -5 );
	MSPt( 2, 0, 1, 0, 0 );	MSPt( 2, 1, 0, 1, 0 );	MSPt( 2, 2, 1, 1, 2 );	MSPt( 2, 3, -1, -1, -1 );
	MSPt( 3, 0, 1, 0, 0 );	MSPt( 3, 1, 0, 1, 0 );	MSPt( 3, 2, 1, -2, 3 );	MSPt( 3, 3, -1, 1, -2 );
	MSPt( 4, 0, 1, 0, 0 );	MSPt( 4, 1, 0, 1, 0 );	MSPt( 4, 2, -2, 1, 5 );	MSPt( 4, 3, 1, -1, -3 );
	MSPt( 5, 0, 1, 0, 0 );	MSPt( 5, 1, 0, 1, 0 );	MSPt( 5, 2, 1, -2, 5 );	MSPt( 5, 3, -1, 1, -4 );
	MSPt( 6, 0, 1, 0, 0 );	MSPt( 6, 1, 0, 1, 0 );	MSPt( 6, 2, 1, -2, 7 );	MSPt( 6, 3, -1, 1, -5 );
	MSPt( 7, 0, 1, 0, 0 );	MSPt( 7, 1, 0, 1, 0 );	MSPt( 7, 2, -2, 2, 7 );	MSPt( 7, 3, 1, -2, -5 );
	
	MSPt( 8, 0, 1, 0, 0 );	MSPt( 8, 1, 0, 1, 0 );	MSPt( 8, 2, 0, 0, 1 );	MSPt( 8, 3, -1, -1, 0 );	MSPt( 8, 4, 0, 0, -1 );
	MSPt( 9, 0, 1, 0, 0 );	MSPt( 9, 1, 0, 1, 0 );	MSPt( 9, 2, 1, 2, 3 );	MSPt( 9, 3, -1, -1, 0 );	MSPt( 9, 4, -1, -2, -3 );
	MSPt( 10, 0, 1, 0, 0 );	MSPt( 10, 1, 0, 1, 0 );	MSPt( 10, 2, 1, 1, 1 );	MSPt( 10, 3, -1, -1, 0 );	MSPt( 10, 4, 0, 0, -1 );
	
	MSPt( 11, 0, 1, 0, 0 );	MSPt( 11, 1, 0, 1, 0 );	MSPt( 11, 2, 0, 0, 1 );	MSPt( 11, 3, -1, 0, 0 );	MSPt( 11, 4, 0, -1, 0 );	MSPt( 11, 5, 0, 0, -1 );
	MSPt( 12, 0, 1, 0, 0 );	MSPt( 12, 1, 0, 1, 0 );	MSPt( 12, 2, 1, 1, 2 );	MSPt( 12, 3, -1, 0, 0 );	MSPt( 12, 4, 0, -1, 0 );	MSPt( 12, 5, -1, -1, -2 );

	
	/* add the minimal polytopes to the list */
	for( i = 0; i < kNumMin; i++ )
		if( err = doAddPolytopeToList( p[i] ) )	return( err );
	
	/* return success */
	return( kNoError );
}

/* doAddPolytopeToList -	call to add the polytope to the found list */
static char doAddPolytopeToList( PolytopePtr p )
{
	/* find the end of the list and add on the polytope */
	if( gPolyList )
	{
		PolyListPtr	foundList = gPolyList;
		
		while( foundList->next )	foundList = foundList->next;
		
		if( !(foundList->next = (PolyListPtr)malloc( sizeof( PolyListRec ) )) )
			return( kMemError );
		p->id = foundList->p->id + 1;
		foundList->next->p = p;
		foundList->next->next = kFalse;
	}
	else
	{
		if( !(gPolyList = (PolyListPtr)malloc( sizeof( PolyListRec ) )) )
			return( kMemError );
		p->id = 1;
		gPolyList->p = p;
		gPolyList->next = kFalse;
	}
	
	/* check whether the polytope is simplicial or not */
	p->simplicial = doIsSimplicial( p );
	
	return( kNoError );
}	

/* doNewPolytope -	call to create a new polytope with the given number of vertices */
static PolytopePtr doNewPolytope( short numVertices )
{
	PolytopePtr	p;
	
	/* allocate the memory for the polytope */
	if( !(p = (PolytopePtr)malloc( sizeof( PolytopeRec ) )) )
	{
		printf( "\nNot enough memory to create new polytope!!!\n\n" );
		return( kFalse );
	}
	
	/* set the polytope's details */
	p->numVertices = numVertices;
	p->numChildren = 0;
	p->numParents = 0;
	p->children = kFalse;
	p->parents = kFalse;

	/* allocate the memory for the vertices */
	if( !(p->vertices = (Point3DPtr)malloc( sizeof( Point3DRec ) * numVertices )) )
	{
		free( (void *)p );
		printf( "\nNot enough memory to create new polytope!!!\n\n" );
		return( kFalse );
	}
		
	/* return the polytope */
	return( p );
}

/* doNewPolytopeChild -	call to create a new child polytope */
static PolytopePtr doNewPolytopeChild( PolytopePtr p )
{
	PolytopePtr	q;
	short		i;
	
	/* allocate the memory for the polytope */
	if( !(q = doNewPolytope( p->numVertices + 1 )) )
		return( kFalse );
	
	/* set the known vertices */
	for( i = 0; i < p->numVertices; i++ )
		q->vertices[i] = p->vertices[i];
	
	/* return the child */
	return( q );
}

/* doDisposePolytope -	call to dispose of a polytope */
static void doDisposePolytope( PolytopePtr p )
{
	if( p )
	{
		PolyListPtr		temp;
		
		/* dispose of the vertices */
		if( p->vertices )	free( (void *)(p->vertices) );
		
		/* dispose of the children list */
		temp = p->children;
		while( temp )
		{
			PolyListPtr	next = temp->next;
			
			free( (void *)temp );
			temp = next;
		}
		
		/* dispose of the parent list */
		temp = p->parents;
		while( temp )
		{
			PolyListPtr	next = temp->next;
			
			free( (void *)temp );
			temp = next;
		}
		
		/* dispose of the polytope */
		free( (void *)p );
	}
}

/* doIsSimplicial -	call to check that all the faces are composed of triangles */
static char doIsSimplicial( PolytopePtr p )
{
	short		i, j, k, l;
	
	/* scan through the vertices in groups of four, looking for coplanar ones */
	for( i = 0; i < p->numVertices; i++ )
		for( j = i + 1; j < p->numVertices; j++ )
			for( k = j + 1; k < p->numVertices; k++ )
				for( l = k + 1; l < p->numVertices; l++ )
					if( doAreCoplanar( p->vertices + i, p->vertices + j, p->vertices + k, p->vertices + l ) )
						if( doIsFace( p, p->vertices + i, p->vertices + j, p->vertices + k ) )
							return( kFalse );
						
	return( kTrue );
}

/* doAreCoplanar -	checks whether the four points are coplanar */
static char doAreCoplanar( Point3DPtr a, Point3DPtr b, Point3DPtr c, Point3DPtr d )
{
	Point3DRec	bma, cma, dma, nabc;
	
	/* calculate the result of shifting the origin to a */
	MSubtract( *b, *a, bma );
	MSubtract( *c, *a, cma );
	MSubtract( *d, *a, dma );
	
	/* calculate the normal */
	MNormal( bma, cma, nabc );
	
	/* check whether d is in the plane */
	if( MDot( dma, nabc ) == 0 )
		return( kTrue );
	
	return( kFalse );
}

/* doIsFace -	checks whether the three points define a face of the polytope */
static char doIsFace( PolytopePtr p, Point3DPtr a, Point3DPtr b, Point3DPtr c )
{
	Point3DRec	o = {0,0,0}, bma, cma, nabc;
	short		i;
	
	/* calculate the result of shifting the origin to a */
	MSubtract( *b, *a, bma );
	MSubtract( *c, *a, cma );
	
	/* calculate the normal */
	MNormal( bma, cma, nabc );
	
	/* check whether it is a face */
	for( i = 0; i < p->numVertices; i++ )
		if( MPointsNotEqual( *a, p->vertices[i] ) && MPointsNotEqual( *b, p->vertices[i] ) && MPointsNotEqual( *c, p->vertices[i] ) )
			if( !doIsInternal( p->vertices + i, &nabc, &o, a ) )
				return( kFalse );
	
	return( kTrue );
}

/* doIsEdge -	checks whether the two points define an edge of the polytope */
static char doIsEdge( PolytopePtr p, Point3DPtr a, Point3DPtr b )
{
	short	i, j;
	
	/* first we wish to find a two more points giving faces which are not coplanar */
	for( i = 0; i < p->numVertices; i++ )
		if( MPointsNotEqual( *a, p->vertices[i] ) && MPointsNotEqual( *b, p->vertices[i] ) && doIsFace( p, p->vertices + i, a, b ) )
			for( j = i + 1; j < p->numVertices; j++ )
				if( MPointsNotEqual( *a, p->vertices[j] ) && MPointsNotEqual( *b, p->vertices[j] ) && doIsFace( p, p->vertices + j, a, b ) && !doAreCoplanar( p->vertices + i, p->vertices + j, a, b ) )
					return( kTrue );
	
	return( kFalse );
}

/* doIsFreeTetrahedron -	call to test whether the tetrahedron {a,b,c,0} is lattice-point free */
static char doIsFreeTetrahedron( Point3DPtr a, Point3DPtr b, Point3DPtr c )
{
	BoundsRec		bbox = {0,0,0,0,0,0};
	Point3DRec		o = {0,0,0}, nabc, noab, noac, nobc, bma, cma, count;
	
	/* set the bounding box of the tetrahedron */
	doAddPointToBoundingBox( &bbox, a );
	doAddPointToBoundingBox( &bbox, b );
	doAddPointToBoundingBox( &bbox, c );
	
	/* calculate the result of shifting the origin to a */
	MSubtract( *b, *a, bma );
	MSubtract( *c, *a, cma );
	
	/* calculate the normal data for the faces */
	MNormal( bma, cma, nabc );
	MNormal( *a, *b, noab );
	MNormal( *a, *c, noac );
	MNormal( *b, *c, nobc );
	
	/* start checking the points in the bounding box, checking that the point isn't a vertex point or the origin */
	for( count.x = bbox.back; count.x <= bbox.front; count.x++ )
		for( count.y = bbox.left; count.y <= bbox.right; count.y++ )
			for( count.z = bbox.bottom; count.z <= bbox.top; count.z++ )
				if( MPointsNotEqual( count, *a ) && MPointsNotEqual( count, *b ) && MPointsNotEqual( count, *c ) && ((count.x != 0) || (count.y != 0) || (count.z != 0)) )
					if( doIsInternal( &count, &nabc, &o, a ) && doIsInternal( &count, &noab, c, &o ) && doIsInternal( &count, &noac, b, &o ) && doIsInternal( &count, &nobc, a, &o ) )
						return( kFalse );
	
	/* return true */
	return( kTrue );
}

/* doAddPointToBoundingBox -	call to add the point to the bounding box */
static void doAddPointToBoundingBox( BoundsPtr bbox, Point3DPtr a )
{
	if( a->x < bbox->back )			bbox->back = a->x;
	else if( a->x > bbox->front )	bbox->front = a->x;
	if( a->y < bbox->left )			bbox->left = a->y;
	else if( a->y > bbox->right )	bbox->right = a->y;
	if( a->z < bbox->bottom )		bbox->bottom = a->z;
	else if( a->z > bbox->top )		bbox->top = a->z;
}

/* doIsInternal -	call to test whether the given point is on the inside of the face or not (d = internal point, a = point on face, n = normal to face) */
static char doIsInternal( Point3DPtr x, Point3DPtr n, Point3DPtr d, Point3DPtr a )
{
	short		parDot = MDot( *d, *n ) - MDot( *a, *n ), norDot = MDot( *x, *n ) - MDot( *a, *n );
	
	if( !parDot )				return( kFalse );
	if( !norDot )				return( kTrue );
	if( parDot * norDot < 0 )	return( kFalse );
	
	return( kTrue );
}

/* doIsChildFano -	call to test whether the child polytope is Fano, if so we recurse on the child */
static char doIsChildFano( PolytopePtr p, Point3DPtr newVertex )
{
	char			wasNew;
	short			i, j;
	PolytopePtr		q, child;
	
	/* check the new vertex isn't actually an old vertex or the origin */
	if( !newVertex->x && !newVertex->y && !newVertex->z )
		return( kNoError );
	for( i = 0; i < p->numVertices; i++ )
		if( MPointsEqual( *newVertex, p->vertices[i] ) )
			return( kNoError );
		
	/* scan through all the possible tetrahedra checking for non-zero, non-vertex lattice points */
	for( i = 0; i < p->numVertices; i++ )
		for( j = i + 1; j < p->numVertices; j++ )
			if( !doIsFreeTetrahedron( p->vertices + i, p->vertices + j, newVertex ) )
				return( kNoError );
	
	/* create the memory for the child polytope and add in the new vertex */
	if( !(q = doNewPolytopeChild( p )) )
		return( kMemError );
	q->vertices[p->numVertices] = *newVertex;
	
	/* check the new polytope against the polytope list and save the results */
	child = doIsNewPolytope( q, &wasNew );
	doAddChildToList( p, child );
		
	/* finish up by inducting if required */
	if( !wasNew )	doDisposePolytope( q );
	else			return( doEnlargePolytope( q ) );
	
	return( kNoError );
}

/* doEnlargePolytope -	call to try and enlarge the given Fano polytope */
static char doEnlargePolytope( PolytopePtr p )
{	
	char			err;
	
	/* try to extend the polytope by adding in a new vertex */
	if( (err = doAddOverVertex( p )) == kNoError )
		if( (err = doAddOverEdge( p )) == kNoError )
			err = doAddOverFace( p );
	
	/* return any errors */
	return( err );
}

/* doAddOverVertex -	call to try to add a vertex over a previous vertex */
static char doAddOverVertex( PolytopePtr p )
{
	short		i;
	
	/* the new possible point is determined using the first possibility of condition (2) */
	for( i = 0; i < p->numVertices; i++ )
	{
		char			err;
		Point3DRec		newPoint;
		
		/* calculate the new vertex */
		newPoint.x = -p->vertices[i].x;
		newPoint.y = -p->vertices[i].y;
		newPoint.z = -p->vertices[i].z;
		
		/* recurse on the new polytope */
		if( err = doIsChildFano( p, &newPoint ) )
			return( err );
	}
	
	return( kNoError );
}

/* doAddOverEdge -	call to try to add a vertex over an edge created by two vertices */
static char doAddOverEdge( PolytopePtr p )
{
	short		i, j;
	
	/* the new possible point is determined using the second possibility of condition (2) */
	for( i = 0; i < p->numVertices; i++ )
		for( j = i + 1; j < p->numVertices; j++ )
		{
			char			err;
			Point3DRec		newPoint;
			
			/* calculate the new vertex */
			newPoint.x = -p->vertices[i].x - p->vertices[j].x;
			newPoint.y = -p->vertices[i].y - p->vertices[j].y;
			newPoint.z = -p->vertices[i].z - p->vertices[j].z;
			
			/* recurse on the new polytope */
			if( err = doIsChildFano( p, &newPoint ) )
				return( err );
		}
	
	return( kNoError );
}

/* doAddOverFace -	call to try to add a vertex over a face created by three vertices */
static char doAddOverFace( PolytopePtr p )
{
	short		i, j, k;
	
	/* the new possible point is determined using the third possibility of condition (2) */
	for( i = 0; i < p->numVertices; i++ )
		for( j = i + 1; j < p->numVertices; j++ )
			for( k = j + 1; k < p->numVertices; k++ )
			{
				char			err;
				
				/* calculate the new vertex and recurse on the new polytope */
				if( (err = doFindNewVertex( p, p->vertices + i, p->vertices + j, p->vertices + k )) != kNoError )
					return( err );
			}
	
	return( kNoError );
}

/* doFindNewVertex -	call to try and satisfy the barycentric conditions */
static char doFindNewVertex( PolytopePtr p, Point3DPtr a, Point3DPtr b, Point3DPtr c )
{
	char			err;

	/* check through all possible barycentric coordinates */
	/* (1,1,1,1) */
	if( (err = doCheckBarrycentric( p, 1, 1, 1, 1, a, b, c )) != kNoError )		return( err );
	/* (1,1,1,2) */
	if( (err = doCheckBarrycentric( p, 1, 1, 1, 2, a, b, c )) != kNoError )		return( err );
	if( (err = doCheckBarrycentric( p, 1, 1, 2, 1, a, b, c )) != kNoError )		return( err );
	if( (err = doCheckBarrycentric( p, 1, 2, 1, 1, a, b, c )) != kNoError )		return( err );
	if( (err = doCheckBarrycentric( p, 2, 1, 1, 1, a, b, c )) != kNoError )		return( err );
	/* (1,1,2,3) */
	if( (err = doCheckBarrycentric( p, 1, 1, 2, 3, a, b, c )) != kNoError )		return( err );
	if( (err = doCheckBarrycentric( p, 1, 2, 1, 3, a, b, c )) != kNoError )		return( err );
	if( (err = doCheckBarrycentric( p, 2, 1, 1, 3, a, b, c )) != kNoError )		return( err );
	if( (err = doCheckBarrycentric( p, 1, 1, 3, 2, a, b, c )) != kNoError )		return( err );
	if( (err = doCheckBarrycentric( p, 1, 2, 3, 1, a, b, c )) != kNoError )		return( err );
	if( (err = doCheckBarrycentric( p, 2, 1, 3, 1, a, b, c )) != kNoError )		return( err );
	if( (err = doCheckBarrycentric( p, 1, 3, 1, 2, a, b, c )) != kNoError )		return( err );
	if( (err = doCheckBarrycentric( p, 1, 3, 2, 1, a, b, c )) != kNoError )		return( err );
	if( (err = doCheckBarrycentric( p, 2, 3, 1, 1, a, b, c )) != kNoError )		return( err );
	if( (err = doCheckBarrycentric( p, 3, 1, 1, 2, a, b, c )) != kNoError )		return( err );
	if( (err = doCheckBarrycentric( p, 3, 1, 2, 1, a, b, c )) != kNoError )		return( err );
	if( (err = doCheckBarrycentric( p, 3, 2, 1, 1, a, b, c )) != kNoError )		return( err );
	/* (1,2,3,5) */
	if( (err = doCheckBarrycentricPerm( p, 1, 2, 3, 5, a, b, c )) != kNoError )	return( err );
	/* (1,3,4,5) */
	if( (err = doCheckBarrycentricPerm( p, 1, 3, 4, 5, a, b, c )) != kNoError )	return( err );
	/* (2,3,5,7) */
	if( (err = doCheckBarrycentricPerm( p, 2, 3, 5, 7, a, b, c )) != kNoError )	return( err );
	/* (3,4,5,7) */
	if( (err = doCheckBarrycentricPerm( p, 3, 4, 5, 7, a, b, c )) != kNoError )	return( err );
	
	return( kNoError );
}

/* doCheckBarrycentric -	call to check with the given barrycentric coordinates */
static char doCheckBarrycentric( PolytopePtr p, short l1, short l2, short l3, short l4, Point3DPtr a, Point3DPtr b, Point3DPtr c )
{	
	Point3DRec		r;
								
	r.x = -l1 * a->x - l2 * b->x - l3 * c->x;
	r.y = -l1 * a->y - l2 * b->y - l3 * c->y;
	r.z = -l1 * a->z - l2 * b->z - l3 * c->z;
	
	/* check that the new vertex is in Z^3, and if so recurse */
	if( !(r.x % l4) && !(r.y % l4) && !(r.z % l4) )
	{
		char				err;
		
		r.x /= l4; r.y /= l4; r.z /= l4;
		
		/* recurse with the new vertex */
		if( err = doIsChildFano( p, &r ) )
			return( err );
	}

	return( kNoError );
}

/* doCheckBarrycentricPerm -	call to check all possible permutations with the given barrycentric coordinates */
static char doCheckBarrycentricPerm( PolytopePtr p, short l1, short l2, short l3, short l4, Point3DPtr a, Point3DPtr b, Point3DPtr c )
{	
	short		bar[4], c1, c2, c3, c4;
	
	/* build a quick table */
	bar[0] = l1; bar[1] = l2; bar[2] = l3; bar[3] = l4;
	
	/* work through all possible permutations of barycentric coordinates */
	for( c1 = 0; c1 < 4; c1++ )
		for( c2 = 0; c2 < 4; c2++ )
			if( c1 != c2 )
				for( c3 = 0; c3 < 4; c3++ )
					if( (c1 != c3) && (c2 != c3) )
						for( c4 = 0; c4 < 4; c4++ )
							if( (c1 != c4) && (c2 != c4) && (c3 != c4) )
							{
								char		err;
								
								if( (err = doCheckBarrycentric( p, bar[c1], bar[c2], bar[c3], bar[c4], a, b, c )) != kNoError )
									return( err );
							}
	return( kNoError );
}

/* doUpdateList -	call to add and update a polytope and the list */
static char doUpdateList( PolytopePtr p, PolyListHandle list )
{
	PolyListPtr		temp = *list, entry;
	
	/* check the polytope isn't already in the list */
	while( temp )
	{
		if( temp->p->id == p->id )		return( kFalse );
		temp = temp->next;
	}
	
	/* allocate the memory */
	if( !(entry = (PolyListPtr)malloc( sizeof( PolyListRec ) )) )
		return( kFalse );
	
	/* assign the child to the list */
	entry->p = p;
	entry->next = kFalse;
	if( *list )
	{
		temp = *list;
		while( temp->next )	temp = temp->next;
		temp->next = entry;
	}
	else
		*list = entry;
	
	return( kTrue );
}

/* doAddChildToList -	call to add the child to the list of children (also adds the parent to the child's list) */
static void doAddChildToList( PolytopePtr p, PolytopePtr child )
{
	/* first we add the child to the parent's list */
	if( doUpdateList( child, &(p->children) ) )
		p->numChildren++;
	
	/* now we add the parent to the child's list */
	if( doUpdateList( p, &(child->parents) ) )
		child->numParents++;
}

/* doIsNewPolytope -	call to check the polytope against the found list and, if necessary, add it to the list */
static PolytopePtr doIsNewPolytope( PolytopePtr p, char *wasNew )
{
	PolyListPtr	foundList = gPolyList;
	
	/* scan through the list looking for polytopes with the same number of vertices and compairing them to p */
	*wasNew = kFalse;
	while( foundList )
	{
		if( foundList->p->numVertices == p->numVertices )
			if( doArePolytopesSimilar( p, foundList->p ) )
				return( foundList->p );
		foundList = foundList->next;
	}
	
	/* the polytope must be new; add it to the list */
	*wasNew = kTrue;
	doAddPolytopeToList( p );
	
	return( p );
}

/* doArePolytopesSimilar -	call to check whether the two given polytopes are the same up to GL(3,Z) */
static char doArePolytopesSimilar( PolytopePtr p, PolytopePtr q )
{
	PolytopePtr	c;
	short		i, j, k;
	
	/* create a copy of p to apply transformations to */
	if( !(c = doNewPolytope( p->numVertices )) )
		return( kFalse );
	
	/* try finding a rotation to switch between the two polytopes */
	for( i = 0; i < p->numVertices; i++ )
		for( j = 0; j < p->numVertices; j++ )
			if( i != j )
				for( k = 0; k < p->numVertices; k++ )
					if( (i != k) && (j != k) )
						if( doRotatePolytope( c, p, q, i, j, k ) )
							if( doArePolytopesSame( c, q ) )
							{
								doDisposePolytope( c );
								return( kTrue );
							}
	
	doDisposePolytope( c );
	
	return( kFalse );
}

/* doRotatePolytope -	call to rotate the polytope with respect to the given data */
static char doRotatePolytope( PolytopePtr c, PolytopePtr p, PolytopePtr q, short i, short j, short k )
{
	short	t[3][3], g = q->vertices[k].x - q->vertices[i].x * p->vertices[2].x - q->vertices[j].x * p->vertices[2].y;
						
	/* construct the candidate rotation */
	if( !(g % p->vertices[2].z) )
	{
		t[2][0] = g / p->vertices[2].z;
		g = q->vertices[k].y - q->vertices[i].y * p->vertices[2].x - q->vertices[j].y * p->vertices[2].y;
		if( !(g % p->vertices[2].z) )
		{
			t[2][1] = g / p->vertices[2].z;
			g = q->vertices[k].z - q->vertices[i].z * p->vertices[2].x - q->vertices[j].z * p->vertices[2].y;
			if( !(g % p->vertices[2].z) )
			{
				short	det;
				
				t[2][2] = g / p->vertices[2].z;
				t[0][0] = q->vertices[i].x;
				t[0][1] = q->vertices[i].y;
				t[0][2] = q->vertices[i].z;
				t[1][0] = q->vertices[j].x;
				t[1][1] = q->vertices[j].y;
				t[1][2] = q->vertices[j].z;
				
				det = t[0][0] * (t[1][1] * t[2][2] - t[1][2] * t[2][1]) + t[0][1] * (t[1][2] * t[2][0] - t[1][0] * t[2][2]) + t[0][2] * (t[1][0] * t[2][1] - t[1][1] * t[2][0]);
				
				if( (det == 1) || (det == -1) )
				{
					/* adjust c accordingly */
					for( g = 0; g < p->numVertices; g++ )
					{
						c->vertices[g].x = p->vertices[g].x * t[0][0] + p->vertices[g].y * t[1][0] + p->vertices[g].z * t[2][0];
						c->vertices[g].y = p->vertices[g].x * t[0][1] + p->vertices[g].y * t[1][1] + p->vertices[g].z * t[2][1];
						c->vertices[g].z = p->vertices[g].x * t[0][2] + p->vertices[g].y * t[1][2] + p->vertices[g].z * t[2][2];
					}
				}
				
				return( kTrue );
			}
		}
	}
	
	return( kFalse );
}	

/* doArePolytopesSame -	call to check whether the two given polytopes are the same */
static char doArePolytopesSame( PolytopePtr p, PolytopePtr q )
{
	short	i, j;
	
	/* check whether the vertices are the same (up to ordering) */
	for( i = 0; i < p->numVertices; i++ )
	{
		char		found = 0;
		
		for( j = 0; j < p->numVertices; j++ )
			if( MPointsEqual( p->vertices[i], q->vertices[j] ) )
				found = 1;
		if( !found )
			return( kFalse );
	}
	
	return( kTrue );
}

/* doAssignIDs -	call to assign the polytope ID numbers */
static void doAssignIDs( void )
{
	short			curID = 1, curVertex = 4, numPolys;
	
	/* step through the list, first in vertex order */
	do
	{
		PolyListPtr		temp = gPolyList;
		short			maxNumChildren = 0, count;
		
		/* find out how many polytopes there are with the current number of vertices */
		numPolys = 0;
		while( temp )
		{
			if( temp->p->numVertices == curVertex )
			{
				numPolys++;
				if( temp->p->numChildren > maxNumChildren )
					maxNumChildren = temp->p->numChildren;
			}
			temp = temp->next;
		}
		
		/* if there were some polytopes, step throug assigning IDs based on the number of children */
		for( count = 0; count <= maxNumChildren; count++ )
		{
			temp = gPolyList;
			while( temp )
			{
				if( (temp->p->numVertices == curVertex) && (temp->p->numChildren == count) )
				{
					temp->p->id = curID;
					curID++;
				}
				temp = temp->next;
			}
		}
		
		/* move on */
		curVertex++;
	} while( numPolys );
}

/* doSaveResults -	call to output the raw data (as a text file) */
static void doSaveResults( void )
{
	short		i, maxID = 0;
	FILE		*dataFile;
	PolyListPtr	temp = gPolyList;
	
	/* find the maximum ID */
	while( temp )
	{
		if( temp->p->id > maxID )		maxID = temp->p->id;
		temp = temp->next;
	}
	
	/* create the new file */
	if( !(dataFile = fopen( "Polytope_Data.txt", "w" )) )
	{
		printf( "Unable to create the polytope data file!!!\n" );
		return;
	}
	
	/* write the data header */
	fprintf( dataFile, "Polytope ID\tNum Vertices\tNum Parent\tNum Children\tSimplicial\tMinimal\tMaximal\nVertex List\nParent List (if any)\nChild List (if any)\n---\n" );
	
	/* output the list data */
	for( i = 1; i <= maxID; i++ )
	{		
		char			found = kFalse;
		
		temp = gPolyList;
		while( temp && !found )
		{
			if( temp->p->id == i )
			{
				PolyListPtr	child = temp->p->children, parent = temp->p->parents;
				
				/* set the found flag to true */
				found = kTrue;
				
				/* fill in the data for the polytope */
				fprintf( dataFile, "%d\t%d\t%d\t%d\t", i, temp->p->numVertices, temp->p->numParents, temp->p->numChildren );
				if( temp->p->simplicial )
				{
					if( !temp->p->numParents )			fprintf( dataFile, "1\t1\t0\n" );
					else if( !temp->p->numChildren )	fprintf( dataFile, "1\t0\t1\n" );
					else								fprintf( dataFile, "1\t0\t0\n" );
				}
				else if( !temp->p->numParents )			fprintf( dataFile, "0\t1\t0\n" );
				else if( !temp->p->numChildren )		fprintf( dataFile, "0\t0\t1\n" );
				else									fprintf( dataFile, "0\t0\t0\n" );
				
				/* write the vertices */
				doWriteVertices( temp->p, dataFile );
				
				/* write the parent list */
				doWriteList( parent, dataFile );
				
				/* write the children list */
				doWriteList( child, dataFile );
				
				/* rule off */
				fprintf( dataFile, "---\n" );
			}
			temp = temp->next;
		}
	}
	
	/* close the file */
	fclose( dataFile );
}

/* doWriteVertices -	call to write the vertices to the given file */
static void doWriteVertices( PolytopePtr p, FILE *dataFile )
{
	short	j;
	
	/* the x entries */
	for( j = 0; j < p->numVertices; j++ )
	{
		fprintf( dataFile, "%d", p->vertices[j].x );
		if( j != p->numVertices - 1 )		fprintf( dataFile, "\t" );
	}
	fprintf( dataFile, "\n" );
	
	/* the y entries */
	for( j = 0; j < p->numVertices; j++ )
	{
		fprintf( dataFile, "%d", p->vertices[j].y );
		if( j != p->numVertices - 1 )		fprintf( dataFile, "\t" );
	}
	fprintf( dataFile, "\n" );
	
	/* the z entries */
	for( j = 0; j < p->numVertices; j++ )
	{
		fprintf( dataFile, "%d", p->vertices[j].z );
		if( j != p->numVertices - 1 )		fprintf( dataFile, "\t" );
	}
	fprintf( dataFile, "\n" );
}

/* doWriteList -	call to write the given list to the given file */
static void doWriteList( PolyListPtr list, FILE *dataFile )
{
	/* write the list */
	if( list )
	{
		short		maxID = 0, count;
		PolyListPtr	temp = list;
		
		/* find the maximum ID */
		while( temp )
		{
			if( temp->p->id > maxID )	maxID = temp->p->id;
			temp = temp->next;
		}
		
		/* now output the list in numerical order */
		for( count = 1; count < maxID; count++ )
		{
			char		found = kFalse;
			
			temp = list;
			while( temp && !found )
			{
				if( temp->p->id == count )
				{
					fprintf( dataFile, "%d\t", count );
					found = kTrue;
				}
				temp = temp->next;
			}
		}
		fprintf( dataFile, "%d\n", maxID );
	}
}
