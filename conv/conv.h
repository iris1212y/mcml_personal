/***********************************************************
 *	Program Name: conv.
 *
 *	A program used to process the data from mcml -
 *	A Monte Carlo simulation of photon distribution in
 *	multi-layered turbid media in ANSI Standard C.
 ****
 *	Creation Date:	11/1991.
 *	Current Date:	6/1992.
 *
 *	Lihong Wang, Ph. D.
 *	Steven L. Jacques, Ph. D.
 *	Laser Biology Research Laboratory - 17
 *	M.D. Anderson Cancer Center
 *	University of Texas
 *	1515 Holcombe Blvd.
 *	Houston, TX 77030
 *	USA
 *
 ****
 *	General Naming Conventions:
 *	Preprocessor names: all capital letters,
 *		e.g. #define PREPROCESSORS
 *	Globals: first letter of each word is capital, no
 *		underscores,
 *		e.g. short GlobalVar;
 *	Dummy variables:  first letter of each word is capital,
 *		and words are connected by underscores,
 *		e.g. void NiceFunction(char Dummy_Var);
 *	Local variables:  all lower cases, words are connected
 *		by underscores,
 *		e.g. short local_var;
 *	Function names or data types:  same as Globals.
 *
 ****
 *	Dimension of length: cm.
 *
 ****/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>

#define PI 3.1415926
#define STRLEN 256		/* String length. */

#define Boolean char

#define SIGN(x) ((x)>=0 ? 1:-1)
#define MAX(x, y) ((x)>(y) ? (x):(y))
#define MIN(x, y) ((x)<(y) ? (x):(y))

#define NULLCONVSTRUCT {0, 0, 0, 0, 0}
#define NULLOUTSTRUCT {0.0,\
		NULL,NULL,NULL,0.0,\
		NULL,NULL,NULL,0.0,\
		NULL,NULL,NULL,0.0,\
		NULL,NULL,NULL,NULL,NULL,\
		0, NULLCONVSTRUCT}

/****************** Stuctures *****************************/

/****
 *	Structure used to describe the geometry and optical
 *	properties of a layer.
 *	z0 and z1 are the z coordinates for the upper boundary
 *	and lower boundary respectively.
 ****/
typedef struct {
  float      z0, z1;		/* z coordinates of a layer. [cm] */
  float      n;		/* refractive index of a layer. */
  float      mua;		/* absorption coefficient. [1/cm] */
  float      mus;		/* scattering coefficient. [1/cm] */
  float      g;		/* anisotropy. */
}           LayerStruct;

/****
 *	Parameters to describe a photon beam.
 *  Pencil: infinitely narrow beam. This is default for the
 *			beam from the mcml output.
 *	Flat:	Flat beam with radius R.
 *	Gaussian:	Gaussian with 1/e2 radius R.
 *	Others: general beam described by points with interpolation.
 ****/
typedef struct {
  enum {
    pencil, flat, gaussian, others
  }           type;		/* beam type. */
  float      P;		/* total power. [J] */
  float      R;		/* radius. [cm]     */
}           BeamStruct;

/****
 *	Input parameters for each independent run.
 *
 *	z and r are for the cylindrical coordinate system. [cm]
 *	a is for the angle alpha between the photon exiting
 *	direction and the surface normal. [radian]
 *
 *	The grid line separations in z, r, and alpha
 *	directions are dz, dr, and da respectively.  The numbers
 *	of grid lines in z, r, and alpha directions are
 *	nz, nr, and na respectively.
 *
 *	The member layerspecs will point to an array of
 *	structures which store parameters of each layer.
 *	This array has (number_layers + 2) elements. One
 *	element is for a layer.
 *	The layers 0 and (num_layers + 1) are for top ambient
 *	medium and the bottom ambient medium respectively.
 *
 *	For convolution, the grid line separations in z, and alpha
 *	directions are still dz, and da respectively.  The numbers
 *	of grid lines in z, and alpha directions are still
 *	nz, and na respectively. However, the grid line separation
 *	and the number of grid lines in r direction are drc and
 *	nrc respectively.
 ****/
typedef struct {
  char        in_fname[STRLEN];	/* name of mcml output file. */
  char        in_fformat;	/* format of mcml output file . */
  /* 'A' for ASCII, */
  /* 'B' for binary. */
  long        num_photons;	/* to be traced. */
  float      Wth;		/* play roulette if photon */
  /* weight < Wth. */

  float      dz;		/* z grid separation.[cm] */
  float      dr;		/* r grid separation.[cm] */
  float      da;		/* alpha grid separation. */
  /* [radian] */
  short       nz;		/* array range 0..nz-1. */
  short       nr;		/* array range 0..nr-1. */
  short       na;		/* array range 0..na-1. */

  short       num_layers;	/* number of layers. */
  LayerStruct *layerspecs;	/* layer parameters. */

  BeamStruct  beam;		/* incident beam of finite size. */
  float      drc;		/* convolution r grid separation.[cm] */
  short       nrc;		/* convolution array range 0..nrc-1. */
  float       eps;		/* relative error in convolution. */
}           InputStruct;

/****
 *	Structure to keep track of what quantities are convolved.
 *	They are initialized to 0, and are set to 1 when the
 *	represented quantity is convolved.
 *
 *	When the beam profile is changed, they are all reset to 0.
 ****/
typedef struct {
  Boolean     Rd_ra;
  Boolean     Rd_r;
  Boolean     A_rz;
  Boolean     Tt_ra;
  Boolean     Tt_r;
}           ConvStruct;

/****
 *	Structures for scored physical quantities
 *	from mcml and to be convolved for photon
 *	beams of finite size.  Therefore, "Out"
 *	here means the output of both mcml and conv.
 *
 *	The member allocated is used to keep the status
 *	of the arrays.  It is set to 1 if all the arrays
 *	are allocated and assigned values.  It is set to
 *	0 otherwise.
 *
 *	z and r represent z and r coordinates of the
 *	cylindrical coordinate system. [cm]
 *	a is the angle alpha between the photon exiting
 *	direction and the normal to the surfaces. [radian]
 *	See comments of the InputStruct.
 *	See manual for the physcial quantities.
 ****/
typedef struct {
  float      Rsp;		/* specular reflectance. [-] */
  float    **Rd_ra;		/* 2D distribution of diffuse */
  /* reflectance. [1/(cm2 sr)] */
  float     *Rd_r;		/* 1D radial distribution of diffuse */
  /* reflectance. [1/cm2] */
  float     *Rd_a;		/* 1D angular distribution of diffuse */
  /* reflectance. [1/sr] */
  float      Rd;		/* total diffuse reflectance. [-] */

  float    **A_rz;		/* 2D probability density in turbid */
  /* media over r & z. [1/cm3] */
  float     *A_z;		/* 1D probability density over z. */
  /* [1/cm] */
  float     *A_l;		/* each layer's absorption */
  /* probability. [-] */
  float      A;		/* total absorption probability. [-] */

  float    **Tt_ra;		/* 2D distribution of total */
  /* transmittance. [1/(cm2 sr)] */
  float     *Tt_r;		/* 1D radial distribution of */
  /* transmittance. [1/cm2] */
  float     *Tt_a;		/* 1D angular distribution of */
  /* transmittance. [1/sr] */
  float      Tt;		/* total transmittance. [-] */

  float    **Rd_rac;		/* convolved data. [J/(cm2 sr)] */
  float     *Rd_rc;		/* [J/cm2]  	 */
  float    **A_rzc;		/* [J/cm3]  	 */
  float    **Tt_rac;		/* [J/(cm2 sr)] */
  float     *Tt_rc;		/* [J/cm2]  	 */

  char        allocated;	/* set to 1 when arrays are allocated. */
  ConvStruct  conved;
}           OutStruct;

/***********************************************************
 *	Routine prototypes for dynamic memory allocation and
 *	release of arrays and matrices.
 *	Modified from Numerical Recipes in C.
 ****/
float     *AllocVector(short, short);
float    **AllocMatrix(short, short, short, short);
void        FreeVector(float *, short, short);
void        FreeMatrix(float **, short, short, short, short);
void        nrerror(char *);

/***********************************************************
 *	Other prototypes.
 ****/
void 
IsoPlot(float **Z,		/* the 2D array Z[i][j]. */
	long int IXmax,
	long int IYmax,		/* the 0<=i<=IXmax, 0<=j<=IYmax. */
	float Dx,
	float Dy);		/* the gridline separations. */
short       GetShort(short, short);
float       GetFloat(float, float);
