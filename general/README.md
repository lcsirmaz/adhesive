### The polymatroids XM and MY are not linear

This part deals with the case when none of the polymatroids *xyab* and
*abuv* are linear, that is, exactly one of the six Ingleton inequalities is
negative.

As in the linear case, we also have that vertices of the adhesive cone
satisfy *x*⟂*aby*, *y*⟂*abx*, *u*⟂*abv*, *v*⟂*abu* as wee as one of
*a*⟂*xyb* and *a*⟂*buv*, and one of *b*⟂*axy* and *b*⟂*auv*. Now the *xyab*
and *abuv* parts are symmetric, as well as swapping *a* and *b*, this we
have to consider two cases only: either *a*⟂*xyb* and *b*⟂*axy*, or *a*⟂*xyb*
and *b*⟂*auv* indicated by the PERL variable `$case`. The 36 possibilities
for which Ingleton inequality is violated is indicated by `$ingelton1` and
`$ingleton2`.

From here the procedure is similar to that of the
[linear](../linear/README.md) case. For each combination all facets of the
truncated cone is generated and then reduced to be minimal one. From this
list [POLCO](https://csb.ethz.ch/tools/software/polco.html) is used to
generate all extremal vertices (rays), and the the utility program
[vertex](../utils/vertex.pl) is used to check if they are indeed internal
points.


