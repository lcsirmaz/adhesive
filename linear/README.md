### The polymatroid XM is linear

This part deals with the case when the left polymatroid on *xyab* is linear.
It happens of and only if all six Ingleton inequalities are non-negative.

When both *xyab* and *abuv* are linear, they are adhesive. (This is true in
general when M has two points, and both XM and MY are linear.) Thus exactly
one of the six Ingleton inequalities for *abuv* is negative. This case is
indicated by the `$ingleton` variable with value from 1 to 6.

All vertices of the adhesive cone satisfy x|aby, y|abx, u|abv, v|abu as it
is proved in the paper. Also, one of a|xyb and a|buv, and one of b|axy and b|auv
also holds. This latter possiblity makes four cases indicated by `$case`.

Then we proceed as follows.  For each possible value of `$ingleton` and
`$case` all inequalities are generated coming from the Shannon inequalities
and the claimed facets.  Next this set is truncated to be independent.  The
resulting inequalities are printed out, and then the program `POLCO` is
executed to generate all extremal vertices of this part.  These vertices are
pulled back to (XM,MY) and checked if they are indeed adhesive.


