### The polymatroid XM is linear

This part deals with the case when the left polymatroid on *xyab* is linear.
It happens if and only if all six Ingleton inequalities are non-negative.

When both *xyab* and *abuv* are linear, they are adhesive. (This is true in
general when M has two points, and both XM and MY are linear.) Thus exactly
one of the six Ingleton inequalities for the right polmatroid *abuv* is 
negative. This case is
indicated by the `$ingleton` variable with value from 1 to 6.

All vertices of the adhesive cone satisfy *x*⟂*aby*, *y*⟂*abx*, *u*⟂*abv*,
*v*⟂*abu* as it is proved in the paper.  Also, one of *a*⟂*xyb* and
*a*⟂*buv*, and one of *b*⟂*axy* and *b*⟂*auv* also holds. This latter
possibility makes four cases indicated by `$case`.

Then we proceed as follows. For each possible value of `$ingleton` and
`$case` all inequalities are generated coming from the Shannon inequalities
and the claimed facets. Next this set is truncated to be independent.  The
resulting inequalities are printed out, and then the program
[POLCO](https://csb.ethz.ch/tools/software/polco.html) is
executed which generates all extremal vertices of this part. The extremal vertices are
pulled back to (XM,MY) and checked if they are indeed adhesive.

Intermediate LP problems are solved using
[GLPK](https://www.gnu.org/software/glpk/) with the C program [gspx.c](../utils/gspx.c)
as wrapper. That program accepts the constraint matrix *M* and the left
hand side, and returns whether there is a feasible solution when all
variables are non-negative. 

The perl program [handle.pl](handle.pl) goes over all cases and checks them.
It aborts if any of the returned vertices is not an internal point.



