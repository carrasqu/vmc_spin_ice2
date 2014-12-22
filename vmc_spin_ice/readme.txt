First, two definitions 
nh is the total number of sites nh=pow(L,3)*16; 

ntetra is the total number of tetrahedra  ntetra=pow(L,3)*2*4;

These are some of the conventions I used:
// even tetrahedron: up (+1) means in; down (-1) means out.
// odd tetrahedron: up (+1) means out; down (-1) means in.


I defined these tables:

int ivic[nh][6];  // each site its 6 neighbours

ivic[i][j] where i=0, nh-1 are each of the sites in the lattice organized exactly as you organized them in your code. j are the nearest neighbors of site j. For each site i there are 6 of them. maybe this is useful for the calculation of the energy.  


int tetra[ntetra][4]; // each tetrahedron and their 4 sites

tetra[i][j]: for a given tetrahedron i=0,ntetra-1, j=0,3 specifies each of the four sites belonging to it.  The tetrahedra i are organized in a similar order the sites. See the picture. 

int connect[nh][2]; // each site connects two tetrahedra

connect[i][j] tells us which two tetrahedra are connected by site i=0,nh-1. 

For instance,  connect will tell me that site 3 connects tetrahedra 0 and tetrahedra 1 and so on.

All the tables are for periodic boundary conditions based on the cubic cell, in the same way that your code is written (I hope!) 

All subroutines I wrote are in the file lat.cpp because  I wanted to have some checks without messing up your files. I have a main there, but I think we should move this functions to your routines.h file once we know everything is find. Last thing, I have used a crappy random number generator for the loop update but we should be using your MersenneTwister algorithm. 
