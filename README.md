# Smush

Smush is a program for computeing the 
[persistent homology](http://en.wikipedia.org/wiki/Persistent_homology) 
of a point cloud.  It attempts to gain a performance boost by
collapsing simplices in order to reduce the dimensions of the 
vector spaces involved.  It was written by Sean Bowman and Itamar 
Gal at the University of Texas, Austin between 2010-2012.

The ideas in smush are roughly based on the paper [The Tidy Set: A Minimal Simplicial Set for
Computing Homology of Clique Complexes](http://www.cs.dartmouth.edu/~afra/papers/manuscript/tidy.pdf) by A. Zomorodian.  We extend his ideas to the computation of barcodes for filtered complexes.
We extend his techniques to computing the persistent homology of 
filtred complexes.  There are other differences, too.

## Building and so forth

First of all, this software is highly expermental, might not work 
for you, is unnecessarily verbose, probably contains horrible bugs,
etc., etc.

Smush uses DJB's redo system for building.  There is a
[nice Python implementation](https://github.com/apenwarr/redo) by 
Avery Pennarun.  (Really worth checking out, redo is super cool.)
You'll need to edit the .do files to point to the directories where
your Boost (and possibly other) libraries live.

Smush also uses Boost and the C++ standard library extensively,
including c++11 features, so you'll need a moderately recent
compiler to get things to work.

## More info

To run the tests, you can do 

     redo test

To build the main program (main_collapse), do

     redo main_collapse

The program takes as arguments the maximum dimension of homology
to compute (k) and a list of scale values to compute the homology.
The inpute format is a file containing a textual list of points in 
R^n, one point per line.  So, to compute the 0th and 1st persistent 
homology of the points in test.dat at the 
scales 0.1, 0.2, and 0.3, you'd run

     ./main_collapse -v 0.1 0.2 0.3 < test.dat

Included in the directory is a file `sphere-300.dat` of points 
sampled uniformly randomly from a 2--sphere in R^3.

## More warnings

This program can easily take a long time and lots of memory 
for even fairly simple point clouds.  The scale parameters are
especially finicky.  Again, this program is not for the 
inexperienced.

Good luck, and enjoy!

