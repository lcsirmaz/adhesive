### Utility perl programs

* [mkalll](mkall.pl) generates all symmetric versions of the inequalities 
* [gspx](gspx.c) C frontent to the [glpk](https://www.gnu.org/software/glpk/) LP solver
* [vertex](vertex.pm) Perl module checking whether (XM) and (MY) are adhesive
* [polco.jar](polco.jar) Copy of [POLCO](https://csb.ethz.ch/tools/software/polco.html) 

Use [mkall](mkall.pl) to generate all bounding facets from
[conditions](../conditions.txt) containing one representative from the 
symmetric variants.

[gspx](gspx.c) is a C helper program which reads the constrain matrix and
right hand side of an LP problem, and returns whether is has a feasible
solution when all variables are non-negative.

The [vertex](vertex.pm) module checks whether two polymatroids specified on
(xyab) and (abuv) are adhesive, namely whether there is a solution for the
missing ranks satisfying the required polymatroid inequalities.


