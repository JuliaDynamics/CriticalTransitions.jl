// Ordered Line Integral Method - Midpoint (OLIM-M) for computing the quasi-potential
// Copyright: Maria Cameron and Daisy Dahiya, June 22, 2017

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define PI 3.141592653589793
#define mabs(a) ((a) >= 0 ? (a) : -(a))
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 ))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define INFTY 1.0e+6
#define TOL 1.0e-12
#define NX 1024
#define NY 1024
#define K 22
#define NCURVEMAX 399


struct mycurve {
	double x;
	double y;
	struct mycurve *next;
	struct mycurve *prev;
};


struct myvector {
  double x;  
  double y;
};


struct mymatrix {
  double a11;
  double a12;
  double a21;
  double a22;
};

struct mysol {
  double g;
  char c;
};  

struct sol_info {
	char type; // solution type:	1 = 1ptupdate, 2 = 2ptupdate, 0 = initialization, 'n' = never reached
	int ind0;
	int ind1;
	double s;
};


int main(void);
struct myvector myfield(struct myvector x); /* B */
double angle(double x,double y);
double length(double x,double y);
void param(void);
void initial_curve(void);
void olim(void);
struct mysol triangle_update(int ind,int ind0,int ind1);
double one_pt_update(int ind,int ind0);			   
void addtree(int ind); /* adds a node to the binary tree
                                 of the "considered" points */ 
void updatetree(int ind); /* updates the binary tree */
void deltree(void); /* deletes the root of the binary tree */
struct mymatrix matrix_inverse(struct mymatrix matr);
struct mymatrix matrix_product(struct mymatrix a,struct mymatrix b);
struct mymatrix transpose(struct mymatrix matr);
struct myvector matr_vec(struct mymatrix matr,struct myvector vec);
double dot_product(struct myvector a,struct myvector b); 
struct myvector vec_lin_comb(struct myvector v1,struct myvector v2,double a,double b);
struct myvector vec_difference(struct myvector v1,struct myvector v2);
struct myvector vec_sum(struct myvector v1,struct myvector v2);
struct myvector getpoint(int ind);
double length_vec(struct myvector x);
double exact_solution(struct myvector x);
double init(struct myvector x);
struct mysol hybrid_nonlin_solver(double a,double b,double u0,double u1,
			struct myvector x0,struct myvector x1,struct myvector b0,struct myvector b1,struct myvector x);
double myfun(double s,double u0,double u1,struct myvector x0,struct myvector x1,
						struct myvector b0,struct myvector b1,
						struct myvector x1mx0,struct myvector b0mb1,struct myvector x);
double geometric_action_line(struct myvector x0, struct myvector x1);
double init(struct myvector x);
/***************************************/

const int nx1=NX-1, ny1=NY-1, nxy=NX*NY; 
int NCURVE;
const int KK = K*K;
int count=0; /* # of considered points */
double h,hx,hy;
int ms[NX*NY]; /* 0 = 'Unknown', 1 = 'Considered', 2 = "in Accepted Front", 3 = "Accepted" */
double g[NX*NY]; /* function to be computed */
int pos[NX*NY]; /* pos(index of mesh pt) = position in binary tree */
int tree[NX*NY]; /* tree(position in the tree) = index of mesh pt */
// double UPS, HUPS; /* unisotropy ratio */
const int neii[8]={1, NX+1, NX, NX-1, -1, -NX-1, -NX, -NX+1 }; /* neighbor's indices */
double Uexact[NX*NY];
struct sol_info solinfo[NX*NY];
int NXY = NX*NY, NMESH;
struct myvector x_ipoint;
struct mycurve *icurve;

// variables for the potential 
double XMIN,XMAX,YMIN,YMAX; // define the computational domain

// choose the vector field b
const char chfield='l'; 
// Vector fields with asymptotically stable equilibrium
// chfield = 'l': b = [-2, -alin; 2*alin, -1][ x; y]; (in the Matlab notations),
// analytic solution  U = 2x^2 + y^2
const double alin = 10.0;
// chfield = 'q' ->  Maier-Stein model; 
// chfield = 'r' -> FitzHugh-Nagumo model with a = 1.2

// Vector fields with stable limit cycle
// chfield = 'b' ->  Brusselator; 
// CHANGE const char *ficurve_name = "circle.txt"; to
// const char *ficurve_name = "icurve_brusselator.txt"; 
// chfield = 'c' -> analytic solution U = (1/2)(r^2-1)^2, l = (y,-x)

// FILES
const char *ficurve_name = "circle.txt"; // input file with points of the initial curve
const char *f_qpot_name = "Qpot.txt"; // output file with the quasipotential
const char *f_solinfo_name = "stype.txt"; // output file with the solution type: init, 1 or 2 point update


/**************************************/
struct myvector myfield(struct myvector x) {
  struct myvector v;
  double aux;

  switch( chfield ) {
	  case 'b':
		  v.x = 1.0 + x.x*x.x*x.y - 4.0*x.x;
		  v.y = 3*x.x - x.x*x.x*x.y;
		  break;
	  case 'c':
	  	  aux = 1.0 - x.x*x.x - x.y*x.y;
		  v.x = x.y + x.x*aux;
		  v.y = -x.x + x.y*aux;
		  break;
	  case 'q':
		  v.x = x.x - x.x*x.x*x.x - 10.0*x.x*x.y*x.y;
		  v.y = -(1.0 + x.x*x.x)*x.y;
		  break;
	  case 'r':
		  v.x = x.x - x.x*x.x*x.x/3.0-x.y;
		  v.y = x.x + 1.2;
		  break;
	  case 'l':
		  v.x = -2*x.x - alin*x.y;
		  v.y = 2*alin*x.x - x.y;
		  break;
	  default:
    printf("chfield = %c, please correct\n",chfield);
    exit(1);
    break;
 }   
  return v;
}

/*************************************/

double exact_solution(struct myvector x) {
  
  if( chfield == 'l' ) return 2.0*x.x*x.x + x.y*x.y;
  else if( chfield == 'c' ) return 0.5*(x.x*x.x+x.y*x.y-1.0)*(x.x*x.x+x.y*x.y-1.0);
  else return 0.0;

}


/*************************************/

void param() {
  int i,j,ind;  

  switch( chfield ) {
	  case 'q':
		  x_ipoint.x=-1.0; x_ipoint.y=0.0; 
		  XMIN = -2.0; XMAX = 2.0;
		  YMIN = -2.0; YMAX = 2.0;
		  break;
	  case 'r':
		  x_ipoint.x=-1.2; x_ipoint.y=-1.2*(1.0-1.2*1.2/3.0); 
		  XMIN = -2.5; XMAX = 2.5;
		  YMIN = -2.5; YMAX = 2.5;
		  break;
	  case 'l':
		  x_ipoint.x = 0.0; x_ipoint.y = 0.0; 
		  XMIN = -1.0; XMAX = 1.0;
		  YMIN = -1.0; YMAX = 1.0;
		  break;
	  case 'c':	  		  		  
		  XMIN = -2.0; XMAX = 2.0;
		  YMIN = -2.0; YMAX = 2.0;
		  break;
	  case 'b':	  		  		  
		  XMIN = 0.0; XMAX = 7.0;
		  YMIN = 0.0; YMAX = 7.0;
		  break;
	  default:
		  printf("chfield = %c, please correct\n",chfield);
		  exit(1);
		  break;
  }		
  
  printf("in param()\n");
  hx = (XMAX - XMIN)/(NX-1);
  hy = (YMAX - YMIN)/(NY-1);
  h = sqrt(hx*hx + hy*hy);
  for( i=0; i<NX; i++ ) {
    for( j=0; j<NY; j++ ) {
      ind = i + NX*j;
	  ms[ind] = 0;
	  g[ind] = INFTY;
	  solinfo[ind].type = 'n';
    }
  }
  // compute the exact solution if it is available
  if( chfield == 'l' || chfield == 'c' ) {
	  for( i=0; i<NX; i++ ) {
		for( j=0; j<NY; j++ ) {
		  ind = i + NX*j;
		  Uexact[ind] = exact_solution(getpoint(ind));
		}
	  }
  }
}

/************************************/

void ipoint() {
  int i,j,ind,ind0,m,n;
  double gtemp;
  const int isur[4]={0, 1, NX+1, NX};
  struct myvector x,l,b;

  i = floor((x_ipoint.x - XMIN)/hx);
  j = floor((x_ipoint.y - YMIN)/hy);
  ind0 = i+j*NX;
  for( m=0; m<4; m++ ) {
    ind = ind0+isur[m];
	x = getpoint(ind);
    gtemp = init(x);
	g[ind] = gtemp;
	ms[ind] = 1;
	addtree(ind);
	solinfo[ind].type = 0;
  }
}


/***********************************/

void initial_curve() {
  int i,j,ind,i0,j0,i1,j1,m,n;
  double xc,yc,xc0,xc1,yc0,yc1;
  FILE *fcurve;
  int Nac = 0, iac[NX*100];
  
   printf("in initial_curve()\n");
 
  fcurve = fopen(ficurve_name,"r");	
  if( fcurve == NULL ) {
  	printf("Please provide a data file for the initial curve\n");
  	exit(1);
  }
  icurve = (struct mycurve *)malloc(NCURVEMAX*sizeof(struct mycurve));
  n = 0;
  while( !feof(fcurve) ) {
  	fscanf(fcurve,"%le\t%le\n",&xc,&yc);
  	(icurve + n) -> x = xc;
  	(icurve + n) -> y = yc;
    n++;
  }
  if( n > NCURVEMAX ) {
  	printf("Error: j = %i while NCURVEMAX = %i\n",n,NCURVEMAX);
  	exit(1);
  }
  fclose(fcurve);
  NCURVE = n;
  for( n = 0; n < NCURVE; n++ ) {
  	(icurve + n) -> next = (n == NCURVE-1) ? icurve : (icurve + n + 1);
  	(icurve + n) -> prev = (n == 0) ? (icurve + NCURVE - 1) : (icurve + n - 1);
  }
  printf("start\n");
  for( n = 0; n < NCURVE; n++ ) {
  	  xc0 = (icurve + n) -> x;
  	  yc0 = (icurve + n) -> y; 	
  	  xc1 = ((icurve + n) -> next) -> x;
  	  yc1 = ((icurve + n) -> next) -> y;
	  i0 = min(floor((xc0 - XMIN)/hx),floor((xc1 - XMIN)/hx));
	  j0 = min(floor((yc0 - YMIN)/hy),floor((yc1 - YMIN)/hy));
	  i1 = max(ceil((xc0 - XMIN)/hx),ceil((xc1 - XMIN)/hx));
	  j1 = max(ceil((yc0 - YMIN)/hy),ceil((yc1 - YMIN)/hy));
	  for( i = i0; i <= i1; i++ ) {
	  	for( j = j0; j <= j1; j++ ) {
	  		ind = i + j*NX;
	  		if( ms[ind] == 0 ) {	  			
	  			g[ind] = init(getpoint(ind)); 
	  			ms[ind] = 1;
	  			solinfo[ind].type = 0;
 	  			addtree(ind);
	  			iac[Nac] = ind;
	  			Nac++;
	  		}
	  	}
	 }
  }
}		
//************************************/

double init(struct myvector x) {
	double a, b, c, d, A, B, C, aux1, aux2, aux;
	int n,nmin,step,imin;
	double dmin,dtemp,lb0,temp;
	struct mycurve *cpoint;
	struct myvector x0,x1,v0,v1,v,bvec0,bvec1,bvec;
	
	switch( chfield ) {
		case 'l':
			x = vec_difference(x,x_ipoint);
			a= -2.0; b = -alin; c = 2.0*alin; d = -1.0;
			aux1 = c - b;
			aux2 = a + d;
			aux = aux1*aux1 + aux2*aux2;
			aux1 *= aux2/aux;
			aux2 *= aux2/aux; 
			A = -(a*aux2 + c*aux1);
			B = -(b*aux2 + d*aux1);
			C = -(d*aux2 - b*aux1);
			return A*x.x*x.x + 2.0*B*x.x*x.y + C*x.y*x.y;
			break;
		case 'c': 
			// find the closest point on the curve to x
			step = max(1,NCURVE/12);
			nmin = 0; dmin = INFTY;
			for( n = 0; n < NCURVE; n+=step ) {
				dtemp = length(((icurve + n) -> x) - x.x,((icurve + n) -> y) - x.y);
				if( dtemp < dmin ) {
					dmin = dtemp;
					nmin = n;
				}
			}
			cpoint = icurve + nmin;
			imin = nmin;
			step = 2*step;
			for( n = 1; n < step; n++ ) { 
				cpoint = cpoint -> next; 
				dtemp = length((cpoint -> x) - x.x,(cpoint -> y) - x.y);
				if( dtemp < dmin ) {
					dmin = dtemp;
					imin = (nmin + n)%NCURVE;
				}
			}	 
			cpoint = icurve + nmin;
			for( n = 1; n < step; n++ ) { 
				cpoint = cpoint -> prev; 
				dtemp = length((cpoint -> x) - x.x,(cpoint -> y) - x.y);
				if( dtemp < dmin ) {
					dmin = dtemp;
					imin = (nmin + NCURVE - n)%NCURVE;
				}
			}
			x0.x = (icurve + imin) -> x;
			x0.y = (icurve + imin) -> y;	 
			x1.x = ((icurve + imin) -> next) -> x;
			x1.y = ((icurve + imin) -> next) -> y;
			v0 = vec_difference(x0,x);
			v1 = vec_difference(x1,x);
			if( dot_product(v0,vec_difference(v0,v1)) < 0.0 ) {
				x1.x = ((icurve + imin) -> prev) -> x;
				x1.y = ((icurve + imin) -> prev) -> y;
				v1 = vec_difference(x1,x);
			}
			v = vec_difference(x1,x0);
			c = length_vec(v);
			dmin = fabs(v0.x*v1.y - v0.y*v1.x)/c;
			a = sqrt(v0.x*v0.x + v0.y*v0.y - dmin*dmin);
			aux = a/c;
			v.x *= aux;
			v.y *= aux;
			x0 = vec_sum(x0,v);
			x1 = vec_lin_comb(x,x0,0.5,0.5);
			bvec0 = myfield(x0);
			bvec1 = myfield(x1);
			bvec = myfield(x);
			lb0 = length_vec(bvec0);
			aux = dot_product(bvec0,bvec1)/(lb0*lb0);
			bvec0.x *= aux;
			bvec0.y *= aux;
			aux = dot_product(bvec0,bvec)/(lb0*lb0);
			bvec.x *= aux;
			bvec.y *= aux;
			temp = dmin*(4.0*length_vec(vec_difference(bvec1,bvec0)) + length_vec(vec_difference(bvec,bvec0)))/3.0;
			return temp;
			break;
		default:
			return 0.0;
			break;
	}	
}			 




/**********************************************/
/*** ordered line integral method ***/

void olim(void) {
  int i,j,k,m,ind,ind0,ii,jj,ind1,indupdate,i0,j0,imin,n;
  int mycount = 0; /* counter for accepted points */
  double gmin,gtemp,gold;
  int Naf, AFneib[8], Nc, NCneib[8]; /* accepted front neighbors and new considered of the newly accepted point */
  struct mysol sol; 
  struct myvector vec,b,b0,b1,v0,v1,vnew; 
  double s, s1;
  char pr = 'n'; // a logic variable: if pr == 'y', solution will be printed out in the terminal window
  
  printf("in olim()\n");

  while( count > 0 ) {
    ind = tree[1];
    vnew = getpoint(ind);
    j = ind/NX;
    i = ind%NX;
    /* x and y of the newly accepted point */
    ms[ind] = 2;
    deltree();
	mycount++;
    if( i<=2 || i>=nx1-2 || j<=2 || j>= ny1-2 || g[ind] >= INFTY-1) {
	  printf("%d\t(%i\t%i) is accepted, g=%.4e\n",mycount,i,j,g[ind]);
	  break; /* quit if we reach the boundary of the computational domain */
	}
	
	if( pr == 'y' ) {
		printf("%i: count = %i, %i is accepted (%.3e,%.3e), %.3e, g=%.4e, Uexact = %.4e, stype = %i",
				mycount,count,ind,vnew.x,vnew.y,length_vec(vnew),g[ind],Uexact[ind],solinfo[ind].type);
		if( solinfo[ind].type == 1 ) printf(", ind0 = %i\n", solinfo[ind].ind0);
		else if( solinfo[ind].type == 2 ) printf(", ind0 = %i, ind1 = %i, s = %.4e\n", solinfo[ind].ind0,solinfo[ind].ind1,solinfo[ind].s);
		else printf("\n");
 	}

    /* Inspect the neighbors of the new Accepted point */
    Naf = 0;
    Nc = 0;
    for( k=0; k<8; k++ ) {
		  ind1 = ind + neii[k]; // neighbor of the new Accepted point
// 		  printf("ind1 = %i, ms = %i\n",ind1,ms[ind1]);
		  // update AcceptedFront
		  if( ms[ind1] == 2 ) {
				m = 0;
				for( n=0; n < 8; n++ ) {
				   ind0 = ind1 + neii[n];
				   if( ms[ind0] < 2 ) m++;
				}   
				if( m == 0 ) { /* ind1 has no considered neighbors */
				  ms[ind1] = 3;
				}
				else {
				  AFneib[Naf] = ind1;
				  Naf++;
				}  
		  }
		  else if( ms[ind1] == 0 ) { // the neighbor ind1 will be a new Considered point
		    	vec = getpoint(ind);
				NCneib[Nc] = ind1;
				Nc++;
		  }  
	}
//	printf("New Accepted: (%i,%i)\n",i,j);
	/* update considered points */
	for( i0 = max(i-K,0); i0 <= min(i+K,nx1); i0++)  {
		ii = i - i0;
		jj = ceil(sqrt(labs(KK - ii*ii)));
		for( j0 = max(j-jj,0); j0 <= min(j+jj,ny1); j0++ ) {
			indupdate = i0 + NX*j0;
						
			if( ms[indupdate] == 1 ) {
				gold = g[indupdate];
				gtemp  = one_pt_update(indupdate,ind);
				vec = getpoint(indupdate);
				if( gtemp < g[indupdate] ) {
					g[indupdate] = gtemp;
					solinfo[indupdate].type = 1;
					solinfo[indupdate].ind0 = ind;
				}
				b0 = myfield(vec_lin_comb(vnew,vec,0.5,0.5)); 
				for( m = 0; m < Naf; m++ ) {
				  ind1 = AFneib[m];
				  v1 = getpoint(ind1);
				  b1 = myfield(vec_lin_comb(v1,vec,0.5,0.5)); 
				  sol = hybrid_nonlin_solver(0.0,1.0,g[ind],g[ind1],vnew,v1,b0,b1,vec);
				  if( sol.c == 'y' ) {
					  s = sol.g;
					  s1 = 1.0 - s;	
					  gtemp = s*g[ind] + s1*g[ind1] + geometric_action_line(vec_lin_comb(vnew,v1,s,s1),vec);	
					  if( gtemp < g[indupdate] ) {
					  	g[indupdate] = gtemp;
						solinfo[indupdate].type = 2;
						solinfo[indupdate].ind0 = ind;
						solinfo[indupdate].ind1 = ind1;
						solinfo[indupdate].s = s;
					  }	
				  }	
				}
				if( gold > g[indupdate] ) {
				  updatetree(indupdate);
				}   
			} 
		}
	}		

		
     /* Shift Unknown neighbors of the new Accepted point to Considered and compute values at them */ 			  
	 for( m = 0; m < Nc; m++ ) { /* for all points that switched from unknown to considered */
		   indupdate = NCneib[m];
		   i = indupdate%NX;
		   j = indupdate/NX;
		   vec = getpoint(indupdate);
// 		   printf("New Considered: i = %i, j = %i, %.3e\n",i,j,length_vec(vec));
		   gmin = INFTY;
		   imin = ind;
		   for( i0 = max(0,i-K); i0 <= min(nx1,i+K); i0++ ) {
		   		ii = i - i0;
				jj = ceil(sqrt(labs(KK - ii*ii)));
		   	   for( j0 = max(j-jj,0); j0 <= min(ny1,j+jj); j0++ ) {
					 ind0 = i0+NX*j0;
					 /* compute tentative values using poins of the accepted front or close accepted poonts */
					 if( ms[ind0] == 2 ) {//|| (ms[ind0]==3 && labs(i-i0)<1.5 && labs(j-j0)<1.5) ) {
						 gtemp = one_pt_update(indupdate,ind0);
						 if( gtemp < gmin ) {
							 gmin = gtemp;
							 imin = ind0;
						 }
					 }
				}
		   }
		   ind0 = imin;	 
		   g[indupdate] = gmin;
		   solinfo[indupdate].type = 1;
		   solinfo[indupdate].ind0 = ind0;
		   v0 = getpoint(ind0);

		   
		   for( k=0; k<8; k++ ) {
				 ind1 = ind0 + neii[k];
				 if( ms[ind1] == 2 ) {
				     v1 = getpoint(ind1);
				     b0 = myfield(vec_lin_comb(v0,vec,0.5,0.5)); 
				     b1 = myfield(vec_lin_comb(v1,vec,0.5,0.5));				     
				     sol = hybrid_nonlin_solver(0.0,1.0,g[ind0],g[ind1],v0,v1,b0,b1,vec);
					  if( sol.c == 'y' ) {
						  s = sol.g;
						  s1 = 1.0 - s;	
						  gtemp = s*g[ind0] + s1*g[ind1] + geometric_action_line(vec_lin_comb(v0,v1,s,s1),vec);
						  if( gtemp < g[indupdate] ) {
						  	g[indupdate] = gtemp;
							solinfo[indupdate].type = 2;
							solinfo[indupdate].ind0 = ind0;
							solinfo[indupdate].ind1 = ind1;
							solinfo[indupdate].s = s;
						  }	
					  }				  
				  } /* end if( ms[ind1] == 2 ) */
		  }	      
		  addtree(indupdate);
		  ms[indupdate] = 1;

	} /* end for( m = 0; m < Nc; m++ ) */

  } /* end while ( count > 0 ) */
}

  
/*********************************************/

double one_pt_update(int ind,int ind0) {
  struct myvector x,x0,l,b,bb;
  double gtemp,ab,len;
  
  x = getpoint(ind);
  x0 = getpoint(ind0);
  gtemp = g[ind0] + geometric_action_line(x0,x);
  
//  fprintf(fup,"%i\t%.6e\t%i\t%i\t%i\t%i\t%.4e\n",ind,gtemp,1,ind0,-1,-1,length(l.x,l.y));

  return gtemp;
}
/*-------------*/

double geometric_action_line(struct myvector x0, struct myvector x1) {
  struct myvector l,b;
  double s,ab,len;
  
  l = vec_difference(x1,x0);
  b = myfield(vec_lin_comb(x0,x1,0.5,0.5));

  return length_vec(b)*length_vec(l) - dot_product(b,l);
}


/*-------------*/  
  
struct myvector getpoint(int ind) {
  struct myvector x0,x1,l;
  int i,ind0,ind1;
  double s;
  
  if( ind < NXY ) {
	  l.x = hx*(ind%NX) + XMIN;
	  l.y = hy*(ind/NX) + YMIN;
//  	  printf("ind = %i, hx = %.4e, hy = %.4e, i = %i, j = %i, (%.4e,%.4e)\n",ind,hx,hy,ind%NX,ind/NX,l.x,l.y);
  }
  return l;
}

			
/***** N o n l i n e a r   1D    s o l v e r *****/


struct mysol hybrid_nonlin_solver(double a,double b,double u0,double u1,
			struct myvector x0,struct myvector x1,struct myvector b0,struct myvector b1,struct myvector x) {
	double c, fa, fb, fc, d, fd, dd, df, dm, ds, t;
	struct myvector x1mx0,b0mb1;
	struct mysol sol;
	double NONLIN_TOL = 1.0e-6;
	int iter = 0, itermax = 100;

	x1mx0 = vec_difference(x1,x0);
	b0mb1 = vec_difference(b0,b1);

	c = a;
	fa = myfun(a,u0,u1,x0,x1,b0,b1,x1mx0,b0mb1,x); 
	fb = myfun(b,u0,u1,x0,x1,b0,b1,x1mx0,b0mb1,x);  
	fc = fa;
	if( (fa > 0 && fb > 0 ) || (fa < 0 && fb < 0) ) {
//	 root is not bracketed 
		sol.c = 'n';
		sol.g = INFTY;
		return sol;
	}
	while( iter < itermax ) {
		if( fabs(fc) < fabs(fb) ) { 
			t = c; c = b; b = t;
			t = fc; fc = fb; fb = t;
			a = c; fa = fc;
		}		
		if( fabs(b - c) < NONLIN_TOL ) break;		
		dm = 0.5*(c - b);
		df = fa - fb;

		if( fabs(df) < NONLIN_TOL ) ds = dm;
		else ds = -fb*(a - b)/df;
		if( (ds > 0 && dm < 0) || (ds < 0 && dm > 0) || (fabs(ds) > fabs(dm)) ) dd = dm;
		else dd = ds;
		
		if( fabs(dd) < NONLIN_TOL ) dd = 0.5*sgn(dm)*NONLIN_TOL;

		d = b + dd;
		fd = myfun(d,u0,u1,x0,x1,b0,b1,x1mx0,b0mb1,x); 
		if( fabs(fd) < NONLIN_TOL ) {
			b = d;
			break;
		}
		a = b; b = d; fa = fb; fb = fd;
		if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
			c = a; fc = fa;
		}
		iter++;
	}
	sol.c = 'y';
	sol.g = b;
	
	return sol;
}


/*--------------*/

double myfun(double s,double u0,double u1,struct myvector x0,struct myvector x1,
						struct myvector b0,struct myvector b1,
						struct myvector x1mx0,struct myvector b0mb1,struct myvector x) {
	double s1 = 1.0 - s,ls,lbs;
	struct myvector xs,xmxs,bs;
	
	xs = vec_lin_comb(x0,x1,s,s1);
	xmxs = vec_difference(x,xs);
	bs = vec_lin_comb(b0,b1,s,s1);
	ls = length_vec(xmxs);
	lbs = length_vec(bs);
	
	if( lbs > TOL ) {	
	return u0 - u1 + dot_product(xmxs,x1mx0)*lbs/ls + dot_product(bs,b0mb1)*ls/lbs 
			- dot_product(x1mx0,bs) - dot_product(b0mb1,xmxs);
	}
	else {
	return u0 - u1 + dot_product(xmxs,x1mx0)*lbs/ls 
			- dot_product(x1mx0,bs) - dot_product(b0mb1,xmxs);
	}
	
}
			



/**********************************************/
/*** linear algebra ***/

double angle(double x,double y) {
  double ang;
  if( y >= 0.0 ) ang = acos(x/sqrt(x*x + y*y));
  else ang = 2.0*PI - acos(x/sqrt(x*x + y*y));
  return ang;
} 

double length(double x,double y) {
  return sqrt(x*x+y*y);
}

double length_vec(struct myvector x) {
  return sqrt(x.x*x.x+x.y*x.y);
}

struct mymatrix matrix_inverse(struct mymatrix matr) {
  struct mymatrix mi;
  double rdet;

  rdet=1.0/(matr.a11*matr.a22-matr.a12*matr.a21);
  mi.a11=matr.a22*rdet;
  mi.a12=-matr.a12*rdet;
  mi.a21=-matr.a21*rdet;
  mi.a22=matr.a11*rdet;
  return mi;
}

struct mymatrix matrix_product(struct mymatrix a,struct mymatrix b) {
  struct mymatrix c;

  c.a11=a.a11*b.a11+a.a12*b.a21;
  c.a12=a.a11*b.a12+a.a12*b.a22;
  c.a21=a.a21*b.a11+a.a22*b.a21;
  c.a22=a.a21*b.a12+a.a22*b.a22;
  return c;
}

struct mymatrix transpose(struct mymatrix matr) {
  struct mymatrix mt;

  mt.a11=matr.a11;
  mt.a22=matr.a22;
  mt.a12=matr.a21;
  mt.a21=matr.a12;
  return mt;
}

struct myvector matr_vec(struct mymatrix matr,struct myvector vec) {
  struct myvector c;

  c.x=matr.a11*vec.x+matr.a12*vec.y;
  c.y=matr.a21*vec.x+matr.a22*vec.y;
  return c;
}

double dot_product(struct myvector a,struct myvector b) {
   return a.x*b.x+a.y*b.y;
}

struct myvector ga_plus_b(double lam,struct myvector a,struct myvector b) {
  struct myvector c;

  c.x=lam*a.x+b.x;
  c.y=lam*a.y+b.y;
  return c;
}

double solve_quadratic(double a,double b,double c) {
  double discr,sol;

  discr=b*b-4.0*a*c;
  if( discr < 0.0 ) sol=INFTY;
  else sol=0.5*(-b+sqrt(discr))/a;
  return sol;
}

struct myvector vec_lin_comb(struct myvector v1,struct myvector v2,double a,double b) {
	struct myvector v;
	
	v.x = a*v1.x + b*v2.x;
	v.y = a*v1.y + b*v2.y;
	
	return v;
}
	
struct myvector vec_difference(struct myvector v1,struct myvector v2) {
	struct myvector v;
	
	v.x = v1.x - v2.x;
	v.y = v1.y - v2.y;
	
	return v;
}
	
struct myvector vec_sum(struct myvector v1,struct myvector v2) {
	struct myvector v;
	
	v.x = v1.x + v2.x;
	v.y = v1.y + v2.y;
	
	return v;
}

/****************************************/
/**************************************************************/
/************ FUNCTIONS RELATED TO THE BINARY TREE ***************/

void addtree(int ind) {
  int loc, ptemp;
  int indp, indc;
  char ch;

  count++;
  tree[count]=ind;
  pos[ind]=count;
  if( count > 1 ) {
    loc=count;
    indc=tree[loc];
    indp=tree[loc/2];
    ch=( g[indc] < g[indp] ) ? 'y' : 'n';
    while( ch == 'y' ) {
      ptemp=pos[indc];
      pos[indc]=pos[indp];
      tree[loc/2]=indc;
      pos[indp]=ptemp;
      tree[loc]=indp;
      loc=loc/2;
      if( loc > 1 ) {
        indc=tree[loc];
        indp=tree[loc/2];
        ch=( g[indc] < g[indp] ) ? 'y' : 'n';
      }
      else ch='n';
    }
  }
}

/*------------------------------------------------------------------*/

void updatetree(int ind) {
  int loc, lcc;
  double g0,g1,g2;

  g0=g[ind];
  loc=pos[ind];
  while( loc > 1 && g0 < g[tree[loc/2]] ) {
    tree[loc]=tree[loc/2];
    pos[tree[loc]]=loc;
    loc=loc/2;
    tree[loc]=ind;
    pos[tree[loc]]=loc;
  }  
  g1=g[tree[loc*2]];
  g2=g[tree[loc*2+1]];
  lcc=count;
  while( (loc*2 <= count && g0 > g[tree[loc*2]]) || (loc*2+1 <= count && g0 > g[tree[loc*2+1]]) )  {
    lcc=( loc*2+1 <=count && g[tree[loc*2+1]] < g[tree[loc*2]] ) ? loc*2+1 : loc*2;
    tree[loc]=tree[lcc];
    pos[tree[loc]]=loc;
    loc=lcc;
    tree[loc]=ind; 
    pos[tree[loc]]=loc;
  }
}

/*---------------------------------------------------------------------*/


/* deletes root of the binary tree */
void deltree() {
  int loc, ptemp, ind, lcc, ic, ic1, ic2, mind;
  char chd, ch='n';;

  mind=tree[1];
  pos[tree[1]]=0;
  tree[1]=tree[count];
  pos[tree[1]]=1;
  count--;
  loc=1;
  ind=tree[1];
  lcc=2*loc;
  if( lcc < count )  {
    ic1=tree[lcc];
    ic2=tree[lcc+1];
    if( (g[ind]) > (g[ic1]) || (g[ind]) > (g[ic2]) ) {
      if( (g[ic1]) <= (g[ic2]) )  {
        chd='l';
	    ic=ic1;
      }
      else {
        chd='r';
	    ic=ic2;
	    lcc++;
      }
    }
    else chd='n';
  }
  else if( lcc == count ) {
    ic=tree[lcc];
    if( (g[ind]) > (g[ic]) ) {chd='l'; if(ch=='y') printf("left\n");}
    else chd='n';
  }
  else chd='n';
  while( chd != 'n' ) {    
    ptemp=pos[ind];
    pos[ind]=pos[ic];
    tree[loc]=ic;
    pos[ic]=ptemp;
    tree[lcc]=ind;
    loc=lcc;
    lcc=2*loc;
    if( lcc < count )  {
      ic1=tree[lcc];
      ic2=tree[lcc+1];
      if( (g[ind]) > (g[ic1]) || (g[ind]) > (g[ic2]) ) {
        if( (g[ic1]) <= (g[ic2]) )  {
          chd='l';
	      ic=ic1;
        }
        else {
          chd='r';
	      ic=ic2;
	      lcc++;
        }
      }
      else chd='n';
    }
    else if( lcc == count ) {
      ic=tree[lcc];
      if(ch=='y') printf("child: loc(%i)=%i, t1=%.12e\n",ic1,lcc,g[ic1]);
      if( (g[ind]) > (g[ic]) ) { chd='l';if(ch=='y') printf("left\n");}
      else chd='n';
    }
    else chd='n';
  } /* end while( chd != 'n' ) */
}


/********************************************************/		    
/*** main ***/

 int main() {
  int i,j,ind,k,si; 
  clock_t CPUbegin;
  double cpu,dd,errmax = 0.0,erms = 0.0;
  FILE *fg, *fs;
  
  param();
  switch( chfield ) {
  	case 'l': case 'q': case 'r':
  		ipoint();
  		break;
  	case 'c': case 'b':
  		initial_curve();
  		break;
  	default:
  		printf("Enter a correct value of the char variable chfield\n");
  		exit(1);
  		break;
  	}			
  CPUbegin=clock();
  olim();
  cpu=(clock()-CPUbegin)/((double)CLOCKS_PER_SEC);
  printf("cputime of olim() = %g\n",cpu);
  /* compute the intencity of the reactive current */
  fg = fopen(f_qpot_name,"w");
  fs = fopen(f_solinfo_name,"w");
  ind=0;
  k = 0;
  for( j=0; j<NY; j++ ) {
    for( i=0; i<NX; i++ ) {
    	if( ms[ind] < 2 ) {
    		g[ind] = INFTY;
    		si = -1;
    	}	
    	else {
    		si = solinfo[ind].type;
			if( chfield == 'l' || chfield == 'c') {
    			dd = fabs(g[ind] - Uexact[ind]);
    			errmax = max(errmax,dd);		
    			erms += dd*dd;
    			k++;
    		}
    	}
        fprintf(fg,"%.12e\t",g[ind]);
        fprintf(fs,"%i\t",si);
        ind++;
    }
    fprintf(fg,"\n");
    fprintf(fs,"\n");
  }
  fclose(fg);
  fclose(fs);
  printf("NX = %i, NY = %i\n",NX,NY);
  if( chfield == 'l' || chfield == 'c') {
  	printf("errmax = %.4e, erms = %.4e\n",errmax,sqrt(erms/k));
  }

  return 0;
}  
