/**************************************************************************/
/*                                                                        */
/*  ws.x (ws = Wigner-Seitz): Generates points of a given lattice that    */
/*                            are inside the WS cell of a given superl.   */
/*  Version 0.1                                                           */
/*                                                                        */
/*  Copyright (C) 1998,2001 PrhoPhi                                       */
/*                                                                        */
/*  This program is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU General Public License as published by  */
/*  the Free Software Foundation; either version 2 of the License, or	  */
/*  (at your option) any later version.                                   */
/*									  */
/*  This program is distributed in the hope that it will be useful,	  */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of	  */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  */
/*  GNU General Public License for more details.			  */
/*                                                                        */
/*  You should have received a copy of the GNU General Public License	  */
/*  along with this program; if not, write to the Free Software		  */
/*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.		  */
/*                                                                        */
/**************************************************************************/


/* Known bug: the algorithm for generating the cluster points could fail to work
   properly in the some special cases. Be careful with the results! */
/* (for cubic systems it works fine) */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* square funtcion */
#define sqr(x)   ( (x) * (x) )

/* square distance between two points */
#define dist(x, y)  (sqr((x)[0]-(y)[0])+sqr((x)[1]-(y)[1])+sqr((x)[2]-(y)[2]))

/* dot product */
#define dot(x, y) ( (x)[0]*(y)[0] + (x)[1]*(y)[1] + (x)[2]*(y)[2] )


#define ARRAY_SIZE 100         /* Inital size of arrays */
#define ARRAY_INC 10         

#define NR_BASIS_VEC 4         /* initial number of basis vectors */

#define SORT_LENGTH 1          /* sorted input */
#define ELIMINATE_DEGEN 1      /* print degenerate points only once */
#define PRINT_BORDER 1         /* print border points of WS cell */
#define PRINT_LENGTH 1         /* print distance from origin for cl. points */
#define PRINT_CLUSTER 1        /* print cluster points */
#define COMPACT 1              /* cluster as compact as possible */

const double epsilon=1.0e-6;   /* allowed inaccurancy */

typedef double vector[3];      
typedef struct
{
  int size;
  int nrel;
  vector *array;
} varray;


#define NR_CELLTYPES 7
const char cellnames[NR_CELLTYPES][80]={ "userdef", "transformation", "cubic",
					 "fcc", "bcc", 
					 "hexagonal", "rotated_hex" };


const vector lattice_vectors[NR_CELLTYPES-2][3]={     /* built in lattices */
  /* cubic */
   {{  1.000000000000, 0.0000000000000, 0.0000000000000 },
    {  0.000000000000, 1.0000000000000, 0.0000000000000 },
    {  0.000000000000, 0.0000000000000, 1.0000000000000 }},
  /* fcc */
   {{  0.500000000000,  0.000000000000,  0.500000000000 },
    {  0.500000000000,  0.500000000000,  0.000000000000 },
    {  0.000000000000,  0.500000000000,  0.500000000000 }},
  /* bcc */
  {{  0.500000000000, -0.500000000000,  0.500000000000 },
   {  0.500000000000,  0.500000000000, -0.500000000000 },
   { -0.500000000000,  0.500000000000,  0.500000000000 }},
  /* hexagonal */
  {{  1.000000000000,  0.000000000000,  0.000000000000 },
   {  0.500000000000,  0.866025403785,  0.000000000000 },
   {  0.000000000000,  0.000000000000,  1.000000000000 }},
   /* rotated_hex */
   {{ 0.866025403785, -0.500000000000,  0.000000000000 },
    { 0.866025403785,  0.500000000000,  0.000000000000 },
    { 0.000000000000,  0.000000000000,  1.000000000000 }}};

#define MAX_LINE_LENGTH 81                   /* max. length of an input line */



const vector zero_vec={0.0, 0.0, 0.0};                 /* zero vector */



/* ws_border: calculates the middle points of the flats of the WS cell */

varray ws_border(const varray , double *);



/* ws_cluster: calculates the points of a given lattice located inside */
/*    of a polygon defined through the middlepoints of the flats    */

varray ws_cluster(const varray, const varray, const varray, const varray, 
		  const vector, const double);



/* delcomment: remove space characters from the begin of the line & */
/*   everything after the command mark (#)                          */

int delcomment(char *);



/* getlatticetype: find out from input line what kind of lattice is meant */

int getlatticetype(const char *);



/* scalelattice: scale lattice with the cell parameters */

int scalelattice(varray, const int, const double, const double);



/* invert_matrix: invert the 3x3 matrix */

int invert_matrix(const varray, varray);



/* error: gives an error message & exits */

int error(const char * );




int main( int argc, char *argv[])
{



  
  
  vector shift;                   /* origin of the lattice */


  varray superl;                  /* super lattice */
  varray superl_rec;              /* inverse matrix of superl */
  varray lattice, basis;             /* lattice and basis vectors */

  varray borderpoints;            /* middle points of the flats of the super
				     lattices WS cell */
  
  varray cluster;                 /* lattice points fitting in the WS cell of
				     the super lattice (=cluster points) */
  
  double min_radius;              /* radius of the biggest sphere fitting in
				     the WS cell of the superlattice */

  double *distances;              /* vector containing the distance of the
				     cluster points from the origin */

  int *indexes;                   /* indexes of array cluster, when items are
				     sorted by their distance from origin */
  
  int super_lattice_type, lattice_type;    /* type of the super lattice
					      & the lattice */

  double cellparam1, cellparam2;     /* cell parameters: edge lengths of the
					Bravais cell */

  int nr_basis;                     /* the number of the basis vectors */



  FILE *input_file;                 /* Pointer to the input file */

  char line[82];                    /* contains the last read line */
  int line_number;                  /* number of the line in line[81] */


  int read_vec;                     /* indicates how many vectors had been
				       read in (for details see code) */

  double vector_x, vector_y, vector_z; /* stores temporaly vector components */

  int coeff1, coeff2, coeff3;       /* border pt in lattice vectors */

  int i1, i2, i3;
  double d1;




  /* Read in the input data from the file given as the fist (and only) arg */


  
  /*  if (argc<2) error("no input file specified!");*/
  if (argc<2) input_file=stdin;
  else {

    input_file= fopen(argv[1], "rt");
    if (!input_file)
      {
	sprintf(line, "can't open input file %s !", argv[1]);
	error(line);
      }
  }

  superl.array=(vector *) malloc(3*sizeof(vector));
  superl.size=3;
  superl.nrel=3;

  superl_rec.array=(vector *) malloc(3*sizeof(vector));
  superl_rec.size=3;
  superl_rec.nrel=3;

  lattice.array=(vector *) malloc(3*sizeof(vector));
  lattice.size=3;
  lattice.nrel=3;

  basis.array=(vector *) malloc(4*sizeof(vector));
  basis.size=NR_BASIS_VEC;
  basis.nrel=0;

  
  if ((!superl.array) || (!superl_rec.array) || (!lattice.array) 
      || (!basis.array))
    error("can't create inital arrays! (Perhaps more memory required.)");



  /* Read in the information from the input file (line for line) */


  line_number=1;
  while (!(feof(input_file)))
    {
      if (fgets(line, MAX_LINE_LENGTH, input_file))
	{
	  delcomment(line);
	  if (!strlen(line)) continue;
	  switch (line_number)
	    {


	      /* supper lattice type */

	    case 1:
	      super_lattice_type=getlatticetype(line);
	      if (super_lattice_type==-1)
		error("wrong cell typpe!");

	      if (super_lattice_type>1) 
		{
		  for (i1=0; i1<3; i1++)
		    for (i2=0; i2<3; i2++)
		      superl.array[i1][i2]=
			lattice_vectors[super_lattice_type-2][i1][i2];
		  line_number++;
		}
	      else 
		{
		  read_vec=0;
		}
	      break;



	      /* super lattice vectors or transformation vectors (only
	         if super lattice type = userdef or transformation) */

	    case 2:
	      if (sscanf(line, "%lf %lf %lf", &vector_x, &vector_y, &vector_z)
		  !=3)
		error("wrong super lattice vectors!");

	      superl.array[read_vec][0]=vector_x;
	      superl.array[read_vec][1]=vector_y;
	      superl.array[read_vec][2]=vector_z;
	      read_vec++;
	      if (read_vec<3) line_number--;
	      else line_number++;
	      break;



	      /* super lattice parameters */

	    case 3:
	      if (sscanf(line, "%lf %lf", &cellparam1, &cellparam2)!=2)
		error("wrong cell parameters!");
	      scalelattice(superl, super_lattice_type, cellparam1, cellparam2);
	      break;



	      /* (normal) lattice type */	     

	    case 4:
	      lattice_type=getlatticetype(line);
	      if (lattice_type==-1 || lattice_type==1)
		error("wrong cell type!");

	      if (lattice_type>1)
		{
		  for (i1=0; i1<3; i1++)
		    for (i2=0; i2<3; i2++)
		      lattice.array[i1][i2]=
			lattice_vectors[lattice_type-2][i1][i2];
		  line_number++;
		}
	      else
		{
		  read_vec=0;
		}
	      break;



	      /*  lattice vectors (only if lattice type = userdef) */

	    case 5:
	      if (sscanf(line, "%lf %lf %lf", &vector_x, &vector_y, &vector_z)
		  !=3)
		error("wrong lattice vectors!");
	      lattice.array[read_vec][0]=vector_x;
	      lattice.array[read_vec][1]=vector_y;
	      lattice.array[read_vec][2]=vector_z;
	      read_vec++;
	      if (read_vec<3) line_number--;
	      break;



	      /* nr. of basis vectors */

	    case 6:
	      if (!sscanf(line, "%d\n", &nr_basis))
		error("wrong number of basis vectors!");
	      if (nr_basis>basis.size)
		{
		  free(basis.array);
		  basis.size=nr_basis;
		  basis.array=(vector *) malloc(nr_basis*sizeof(vector));
		  if (!basis.array)
		    error("can't create array for basis vectors!");
		}
	      read_vec=0;
  	      break;



	      /* basis vectors (at least one) */

	    case 7:
	      if (sscanf(line, "%lf %lf %lf", &vector_x, &vector_y, &vector_z)
		  !=3)
		error("wrong basis vectors!");
	      basis.array[read_vec][0]=vector_x;
	      basis.array[read_vec][1]=vector_y;
	      basis.array[read_vec][2]=vector_z;
	      basis.nrel++;
	      read_vec++;
	      if (read_vec<nr_basis) line_number--;
	      else
		if (lattice_type<2) line_number++;
	      break;


	      
	      /* lattice parameters for the (normal) lattice */

	    case 8:
	      if (sscanf(line, "%lf %lf", &cellparam1, &cellparam2)!=2)
		error("wrong cell parameters (you have to specicify 2 numbers even in\n case of cubic, fcc, bcc lattices)");

	      scalelattice(lattice, lattice_type, cellparam1, cellparam2);
	      scalelattice(basis, lattice_type, cellparam1, cellparam2);

	      break;



	      /* shift vector */

	    case 9:
	      if (sscanf(line, "%lf %lf %lf", &vector_x, &vector_y, &vector_z)
		  !=3)
		error("wrong shift vectors!");

	      shift[0]=vector_x;
	      shift[1]=vector_y;
	      shift[2]=vector_z;
	      break;
	    }
	  line_number++;
	}
    }

  if (line_number!=10) error("input file has not the correct length!");



  /* Make the transformation lattice -> superlattice if necessary */

  if (super_lattice_type==1)
    {
      for (i1=0; i1<3; i1++)
	for (i2=0; i2<3; i2++)
	  superl_rec.array[i1][i2]=superl.array[i1][0]*lattice.array[0][i2]
	    +superl.array[i1][1]*lattice.array[1][i2]
	    +superl.array[i1][2]*lattice.array[2][i2];

      for (i1=0; i1<3; i1++)
	for (i2=0; i2<3; i2++)
	  superl.array[i1][i2]=superl_rec.array[i1][i2];
    }

  if (invert_matrix(superl, superl_rec))
    error("superlattice vectors are not independent!");



  /* Print out input data */

  fprintf(stdout,"super lattice: \n");
  for (i1=0; i1<3; i1++)
    {
      for (i2=0; i2<3; i2++)
	fprintf(stdout,"%14.8lf", superl.array[i1][i2]);
      fprintf(stdout,"\n");
    }


  fprintf(stdout, "reciprocal super lattice (in 2Pi units): \n");
  for (i1=0; i1<3; i1++)
    {
      for (i2=0; i2<3; i2++)
	fprintf(stdout, "%14.8lf", superl_rec.array[i1][i2]);
      fprintf(stdout, "\n");
    }

  fprintf(stdout, "lattice: \n");
  for (i1=0; i1<3; i1++)
    {
      for (i2=0; i2<3; i2++)
	fprintf(stdout, "%14.8lf", lattice.array[i1][i2]);
      fprintf(stdout, "\n");
    }

  fprintf(stdout, "basis: \n");
  for (i1=0; i1<basis.nrel; i1++)
    {
      for (i2=0; i2<3; i2++)
	fprintf(stdout, "%14.8lf", basis.array[i1][i2]);
      fprintf(stdout, "\n");
    }



  
 /* calculating and (if demanded) printing the middle points of neighbour
    WS cells */
  
  if (PRINT_BORDER || PRINT_CLUSTER)
    {
      borderpoints=ws_border(superl, &min_radius);
      if (!borderpoints.nrel)
	error("can't create border points!");
    }


  if (PRINT_BORDER)
    {
      fprintf(stdout, "Neighbour WS cell's middle points:\n");
      for (i1=0; i1<borderpoints.nrel; i1++)
	{
	  vector_x=2*borderpoints.array[i1][0];
	  vector_y=2*borderpoints.array[i1][1];
	  vector_z=2*borderpoints.array[i1][2];
	  coeff1=floor((superl_rec.array[0][0]*vector_x
	    +superl_rec.array[1][0]*vector_y
	    +superl_rec.array[2][0]*vector_z)+0.5);
	  coeff2=floor((superl_rec.array[0][1]*vector_x
	    +superl_rec.array[1][1]*vector_y
	    +superl_rec.array[2][1]*vector_z)+0.5);
	  coeff3=floor((superl_rec.array[0][2]*vector_x
	    +superl_rec.array[1][2]*vector_y
	    +superl_rec.array[2][2]*vector_z)+0.5);

	  fprintf(stdout, 
		  "%4d:  %14.8lf %14.8lf %14.8lf       ( %3d,%3d,%3d )\n",
		  i1+1,vector_x, vector_y, vector_z,
		  coeff1, coeff2, coeff3);
	}
      fprintf(stdout, "Minimal radius: %10.7lf\n", sqrt(min_radius));
    }




  /* if demanded calculating and printing the cluster points */

  if (PRINT_CLUSTER)
    {
  
      cluster=ws_cluster(borderpoints, superl_rec, lattice, basis, shift,
			 min_radius);
      if (!cluster.nrel)
	error("can't create cluster points!");


      if (SORT_LENGTH || PRINT_LENGTH)
	{
	  distances=(double *) malloc(cluster.nrel*sizeof(double));
	  if (!distances)
	      error("can't create distance or index vector!");


          for (i1=0; i1<cluster.nrel; i1++)
            distances[i1]=sqrt(dist(cluster.array[i1], zero_vec));
	}

      indexes=(int *) malloc(cluster.nrel*sizeof(int));
      if (!indexes)
	error("can't create index array!");
      for (i1=0; i1<cluster.nrel; i1++)
        indexes[i1]=i1;



      /*  Sorting the cluster point indexes by the distance of the cluster */
      /*  points to the origin  (bubble sort algorithm)                    */
      
      if (SORT_LENGTH)
	{
	  for (i1=0; i1<cluster.nrel; i1++)
	    {
	      for (i2=i1+1; i2<cluster.nrel; i2++)
		{
		  if (distances[i2]<distances[i1])
		    {
		      d1=distances[i1];
		      distances[i1]=distances[i2];
		      distances[i2]=d1;
		      i3=indexes[i1];
		      indexes[i1]=indexes[i2];
		      indexes[i2]=i3;
		    }
		}
	    }

	}


      /* Printing the cluster points */
      
      fprintf(stdout, "Cluster points:\n");
      for (i1=0; i1<cluster.nrel; i1++)
	{
	  i2=indexes[i1];
	  if (PRINT_LENGTH)
	    {
	      fprintf(stdout, "%4d:  %14.8lf %14.8lf %14.8lf %19.8lf\n",
		      i1+1, cluster.array[i2][0], cluster.array[i2][1],
		     cluster.array[i2][2], distances[i1]);
	    }
	  else
	    {
	      fprintf(stdout, "%3d:   %14.8lf   %14.8lf   %14.8lf\n", i1+1,
		     cluster.array[i2][0],
		     cluster.array[i2][1],
		     cluster.array[i2][2]);
	    }
	}
    }


  /* Everybody is born free, so make everybody free! */
  
  free(cluster.array);
  free(superl.array);
  free(superl_rec.array);
  free(lattice.array);
  free(borderpoints.array);
  free(basis.array);

  return 0;
}





/* ws_border: calculates the middle points of the flats of the WS cell */
/*     input: lattice_vec: lattice vectors                             */
/*    return: ws_border                                                */
/*            min_radius: radius of the biggest sphare fitting in the  */
/*                        WS cell                                      */
       


varray ws_border(const varray lattice_vec, double *min_radius)
{

   varray borderpt;                /* middle points of the WS cells flats */
   
   int coord_max;                  /* */
   
   int coord1, coord2, coord3;     /* coordinates of the current point in
				      lattice vectors */

   int nr_new_pt;                  /* number of new points found in the WS cell
				      when coord_max had a given worth */
   
   int is_inside;                  /* flags if the current point is inside of
				      the WS cell */

   double dist_orig;               /* distance of current point from origin */

   vector *actual_bpt;            /* cartesian coordinates of current point */
   vector *old_bpt;               /* cartesian coords of an other already
				     calculated points in the WS cell */

   int i1, i2;

   /* initializing borderpt */
   
   borderpt.nrel=0;
   borderpt.size=ARRAY_SIZE;
   borderpt.array=(vector *) malloc(borderpt.size*sizeof(vector));
   if (!borderpt.array)
     return borderpt;

   *min_radius=3e30;
   borderpt.array[0][0]=1e15;
   borderpt.array[0][1]=1e15;
   borderpt.array[0][2]=1e15;
   borderpt.nrel=1;

   
   coord_max=1;

   /* Calculating border points */
   
  do
    {
      nr_new_pt=0;
      for (coord1=-coord_max; coord1<coord_max+1; coord1++)
	{
          for (coord2=-coord_max; coord2<coord_max+1; coord2++)
	    {
	      for (coord3=-coord_max; coord3<coord_max+1; coord3++)
		{
		  /* calculate only points that haven't calculated yet */
		  
		  if ((coord1!=-coord_max && coord1!=coord_max)
		      && (coord2!=-coord_max && coord2!=coord_max)
		      && (coord3!=-coord_max && coord3!=coord_max))
		    continue;

		  actual_bpt=borderpt.array+borderpt.nrel;

		  (*actual_bpt)[0]=0.5*(coord1*lattice_vec.array[0][0]
					+coord2*lattice_vec.array[1][0]
					+coord3*lattice_vec.array[2][0]);
		  (*actual_bpt)[1]=0.5*(coord1*lattice_vec.array[0][1]
					+coord2*lattice_vec.array[1][1]
					+coord3*lattice_vec.array[2][1]);
		  (*actual_bpt)[2]=0.5*(coord1*lattice_vec.array[0][2]
					+coord2*lattice_vec.array[1][2]
					+coord3*lattice_vec.array[2][2]);

		  is_inside=1;


		  for (i1=0, old_bpt=borderpt.array; i1<borderpt.nrel;
		       i1++, old_bpt++)
		    {
		      /* point (a) is in the innerside of the poligon
		         defined through the points in borderpt, if
			 <a|b> <= <b|b> for every b in borderpt */
		     
		      if (dot(*actual_bpt,*old_bpt)+epsilon>
                         dot(*old_bpt, *old_bpt))
			{
			  is_inside=0;
			  break;
			}
		    }

		  if (is_inside)
		    {
                      dist_orig=dot(*actual_bpt, *actual_bpt);
		      if (dist_orig<*min_radius)
			*min_radius=dist_orig;

		      
		      /* if an old point in borderpt is outside from the plane
		         going through actual_bpt and perpendicular to its
			 space vector => not in WS cell => remove it */

		      for (i1=0, old_bpt=borderpt.array; i1<borderpt.nrel;
			   i1++, old_bpt++)
			{
			  if (dot(*actual_bpt, *old_bpt)+epsilon
			      >dist_orig)
			    {
                              for (i2=i1; i2<borderpt.nrel; i2++)
                                {
                                  borderpt.array[i2][0]=
				    borderpt.array[i2+1][0];
                                  borderpt.array[i2][1]=
				    borderpt.array[i2+1][1];
                                  borderpt.array[i2][2]=
				    borderpt.array[i2+1][2];
                                }
                              borderpt.nrel--;
                              i1--;
                              old_bpt--;
			    }
			}
                      nr_new_pt+=1;
		      borderpt.nrel+=1;


		      /* if borderpt is too small, let's enlarge it */
		      
                      if (borderpt.nrel==borderpt.size)
                        {
                          borderpt.size+=ARRAY_INC;
                          borderpt.array=(vector *) realloc(borderpt.array,
                                                 borderpt.size*sizeof(vector));
			  if (!borderpt.array)
			    {
			      return borderpt;
			      exit(-1);
			    }
                        }
		      
		    }
		}
	    }
	}
    coord_max+=1;
    } while (nr_new_pt);

  return borderpt;
}




/* ws_cluster: calculates the points of a lattice (a) located inside of  */
/*    the WS cell of lattice (b)                                         */
/*                                                                       */
/*    input: border_pt : middle points of the WS cells of lattice (b)    */
/*           super_rec : reciprocal vectors of lattice (b) in 2*Pi units */
/*           lattice   : lattice vectors of lattice (a)                  */
/*           basis     : basis vectors of lattice (a)                    */
/*           shift     : coordinates of the origin of lattice (a)        */
/*           min_radius: radius of the bigges sphere fitting in          */
/*	                   the WS cell of lattice (a)                    */
/*   return: ws_cluster                                                  */


varray ws_cluster(const varray border_pt, const varray super_rec, 
		  const varray lattice, const varray basis, 
		  const vector shift, const double min_radius)
{

  int coord_max;
  int coord1, coord2, coord3;          /* coords of the current lattice point
					  in lattice vectors */

  int bvnr;                          /* current basis vector number */

  int is_inside;                     /* flags if the current point is in the
					WS cell */
  varray cluster;                    /* points in the WS cells (cluster pts) */
  
  vector *actual_pt;                 /* current border point */
  vector *actual_cpt;                /* current cluster point */
  vector *old_cpt, *old_cpt2;        /* an already calculated cluster point */
  
  double dx, dy, dz;                 /* x,y,z component of the difference
				        between two cluster points */

  double coeff1, coeff2, coeff3;    /* coefficients when difference vector
				       expressed as sum of super gr. vectors */

  double frac1, frac2, frac3;       /* fractional part of coeff(i) */

  int nr_new_pt;                    /* number of new cluster points found
				       in the last cycle */

  double dist_old, dist_new;       /* sum of distance of an old/new cluster 
				      member from every other cluster member */

  int i1, i2;
  // int debug;
  double tmp;


  /* initializing array */

  cluster.size=ARRAY_SIZE;
  cluster.nrel=0;
  cluster.array=(vector *) malloc(cluster.size*sizeof(vector));
  if (!cluster.array)
    return cluster;


  coord_max=0;



  /* Calculating cluster points */
  

  do
    {
      nr_new_pt=0;
      for (coord1=coord_max; coord1>-coord_max-1; coord1--)
	{
	  for (coord2=coord_max; coord2>-coord_max-1; coord2--)
	    {
	      for (coord3=coord_max; coord3>-coord_max-1; coord3--)
		{

		  if ((coord1!=-coord_max && coord1!=coord_max)
		       && (coord2!=-coord_max && coord2!=coord_max)
		       && (coord3!=-coord_max && coord3!=coord_max))
		    continue;

                  for (bvnr=0; bvnr<basis.nrel; bvnr++)
		    {

		      actual_cpt=cluster.array+cluster.nrel;
		      (*actual_cpt)[0]=
			coord1*lattice.array[0][0]
                        +coord2*lattice.array[1][0]
			+coord3*lattice.array[2][0]
                        +shift[0]+basis.array[bvnr][0];
		      (*actual_cpt)[1]=
			coord1*lattice.array[0][1]
                        +coord2*lattice.array[1][1]
			+coord3*lattice.array[2][1]
                        +shift[1]+basis.array[bvnr][1];
		      (*actual_cpt)[2]=
			coord1*lattice.array[0][2]
                        +coord2*lattice.array[1][2]
			+coord3*lattice.array[2][2]
                        +shift[2]+basis.array[bvnr][2];


		      if (dot(*actual_cpt, *actual_cpt)+epsilon>min_radius)
			{

			  is_inside=1;

			  /* actual_pt (a) is in the WS cell if for every
			     border point (b): <a|b> <= <b|b> */

			  for (i1=0, actual_pt=border_pt.array;
                               i1<border_pt.nrel; i1++, actual_pt++)
			    {
			      
			      if (dot(*actual_pt, *actual_cpt)-epsilon>
				  dot(*actual_pt, *actual_pt))
				{
				  is_inside=0;
				  break;
				}
			    }

			  if (is_inside)
			    {
			      /* cluster point is degenerate if the difference
				 vector to any other cluster point is a sum of
				 integer number * super lattice vectors
				 => forget it */
			      
                              if (ELIMINATE_DEGEN)
				{
			       
				  for (i1=0, old_cpt=cluster.array;
				       i1<cluster.nrel; i1++, old_cpt++)
				    {
				      dx=(*old_cpt)[0]-(*actual_cpt)[0];
				      dy=(*old_cpt)[1]-(*actual_cpt)[1];
				      dz=(*old_cpt)[2]-(*actual_cpt)[2];

				      coeff1=super_rec.array[0][0]*dx
					+super_rec.array[1][0]*dy
					+super_rec.array[2][0]*dz;
				      coeff2=super_rec.array[0][1]*dx
					+super_rec.array[1][1]*dy
					+super_rec.array[2][1]*dz;
				      coeff3=super_rec.array[0][2]*dx
					+super_rec.array[1][2]*dy
					+super_rec.array[2][2]*dz;

				      frac1=fabs(modf(coeff1, &dx));
				      if (frac1>0.5) frac1=1.0-frac1;
				      frac2=fabs(modf(coeff2, &dy));
				      if (frac2>0.5) frac2=1.0-frac2;
				      frac3=fabs(modf(coeff3, &dz));
				      if (frac3>0.5) frac3=1.0-frac3;
				      
				      if (frac1<epsilon && frac2<epsilon
					  && frac3<epsilon)
					{
					  if (COMPACT)
					    {
					      dist_old=1e15;
					      dist_new=1e15;
					      for (i2=0, old_cpt2=cluster.array;
						   i2<cluster.nrel; 
						   i2++, old_cpt2++)
						{
						  if (i2!=i1) {
						    tmp =
						      sqr((*old_cpt2)[0]
							  -(*old_cpt)[0])
						      +sqr(*(old_cpt2)[1]
							   -(*old_cpt)[1])
						      +sqr((*old_cpt2)[2]
							   -(*old_cpt)[2]);
						    if (dist_old > tmp) 
						      dist_old = tmp;
						    tmp = 
						      sqr((*old_cpt2)[0]
							  -(*actual_cpt)[0])
						      +sqr((*old_cpt2)[1]
							   -(*actual_cpt)[1])
						      +sqr((*old_cpt2)[2]
							   -(*actual_cpt)[2]);
    						    if (dist_new > tmp)
						      dist_new = tmp;
						  }
						}
					      if (dist_old>dist_new) 
						{
						  (*old_cpt)[0]=
						    (*actual_cpt)[0];
						  (*old_cpt)[1]=
						    (*actual_cpt)[1];
						  (*old_cpt)[2]=
						    (*actual_cpt)[2];
						}   
					    }
					  cluster.nrel-=1;
					  nr_new_pt-=1;
					  
					  break;
					}
				      			
				    }
				}
                              cluster.nrel+=1;
                              nr_new_pt+=1;
			    }
			}
		      else
			{
			  cluster.nrel+=1;
			  nr_new_pt+=1;
			}

		      /* if array cluster is to small => enlarge it */
		      
		      if (cluster.nrel==cluster.size)
                        {
                          cluster.size+=ARRAY_INC;
                          cluster.array=(vector *) realloc(cluster.array,
							   sizeof(vector)*
							   cluster.size);
			  if (!cluster.array)
			    {
			      return cluster;
			      exit(-1);
			    }
                        }
		    }

		}
	    }
	}
      coord_max+=1;
    } while (coord_max<30);/* while (nr_new_pt);*/

  return cluster;

}


    
int getlatticetype(const char *input_line)
{

  int celltype;

  for (celltype=0; celltype<NR_CELLTYPES; celltype++)
    {
      if (strstr(input_line, cellnames[celltype]))
	break;
    }
  if (celltype==NR_CELLTYPES)
    return(-1);
  else
    return(celltype);

}



/* delcomment: remove space characters from the begin of the line &     */
/*             everything after the command mark (#)                    */
/*     input:  input_line the line last read in                         */
/*    return:  delcomment: 0                                            */
/*             input_line: the modified input_line                      */

int delcomment(char *input_line)
{
  char *char_pos1;
  // char *char_pos2;
  int i1;

  char_pos1=strchr(input_line, '#');
  if (char_pos1)
    *char_pos1='\0';

  for (i1=0, char_pos1=input_line; i1<strlen(input_line); i1++, char_pos1++)
    {
      if (!isspace(*char_pos1)) break;
    }
  strcpy(input_line, char_pos1) ;

  return 0;
}



/* scalelattice: scale lattice with the cell parameters                    */
/*     input: lattice:       lattice to be scaled                          */
/*            lattice_type:  type of the lattice (0 userdef, 1 cubic ...)  */
/*            cellparam1: multiplicating factor for x and y coordinates    */
/*            cellparam2: multiplicating factor for z cooordinate          */
/*    return: scalelattice: 0                                              */
/*            lattice:      scaled lattice                                 */

int scalelattice(varray lattice, const int lattice_type,
		 const double cellparam1, const double cellparam2)
{
  int i1, i2;

  switch (lattice_type)
    {
    case 0: break;
    case 1: break;
    case 2:
    case 3:
    case 4:
      for (i1=0; i1<lattice.nrel; i1++)
	for (i2=0; i2<3; i2++)
	  lattice.array[i1][i2]*=cellparam1;
      break;
    case 5:
    case 6:
      for (i1=0; i1<lattice.nrel; i1++)
	{
	  lattice.array[i1][0]*=cellparam1;
	  lattice.array[i1][1]*=cellparam1;
	  lattice.array[i1][2]*=cellparam2;
	}
      break;
    }

  return 0;
}



/* invert_matrix: invert the 3x3 matrix                            */
/*         input: superl:     3x3 matrix to be inversed            */
/*        return: invert_matrix: 1 if determinant=0.0, 0 otherwise */
/*                superl_rec:    inverse matrix                    */

int invert_matrix(const varray superl, varray superl_rec)
{
  int i1, i2;
  int index1p, index1m, index2p, index2m;
  double determinant;


  for (i1=0; i1<3; i1++)
    {
      index1m=(i1-1+3) % 3;
      index1p=(i1+1) % 3;
      for (i2=0; i2<3; i2++)
	{
	  index2p=(i2+1) % 3;
	  index2m=(i2-1+3) % 3;
	  superl_rec.array[i2][i1]=
	    superl.array[index1p][index2p]*superl.array[index1m][index2m]
	    -superl.array[index1m][index2p]*superl.array[index1p][index2m];
	}
    }

  determinant=superl.array[0][0]*superl_rec.array[0][0]
              +superl.array[0][1]*superl_rec.array[1][0]
	      +superl.array[0][2]*superl_rec.array[2][0];

  if (fabs(determinant)<epsilon) return 1;

  for (i1=0; i1<3; i1++)
    for (i2=0; i2<3; i2++)
      superl_rec.array[i1][i2]/=determinant;

  return 0;

}



/*  error: gives an error message & exits                */
/*  input: errormsg: the string containing the error msg */
/* return: error: 0                                      */

int error(const char *errormsg)
{
  fprintf(stderr, "ws.x::error:%s\n", errormsg);
  exit(-1);
  
  return 0;
}
