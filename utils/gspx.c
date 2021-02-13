/* simplex program to do the hard work
*
* usage: [-n | -x ] <datafile> <outputfile>
*  -n    noshuffle (don't shuffle rows and columns)
*  -x    feasible solution only -- not an error if no solution exists
*
*  each variable is >=0, 
*  datafile:
*     <columns> <rows> <V>
*   matrix: first row, second row, ... last row (rowwise) (double numbers)
*   right hand side (rows many):
*     <num> ( if =), g<num> (if >=) l<num> (if <=)
*  
*   If '-x' argument is given, that's all (asking whenther a feasible solution exists)
*   otherwise the goal: cols many coeffs of the goal
*
*  result file:
*    not created if there was an error (res!=0), or no solution exists (escept for -x)
*    first line: "round=%d, V=<value>"   (XXXX if -x and no feasible solution exists)
*    second line: base: <idx>:<value>, for non-zero values only.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static inline double mfloor(double x){
    return x<0.0 ? -1.0 + ((long int)(-x)): (double)((long int)(x)); 
}
static inline int mrandom(int v)
{  return v<=1? 0 : (random()&0xffffff)%v;
}

/*-------------------------------------------------------------------------*/
#include "glpk.h"
static glp_prob *P=NULL;	/* glpk simplex structure */
static glp_smcp parm;		/* parameters */
/*-----------------------------------------------------------------------*/
#include <malloc.h>

int cols,rows;
double V;
double *spx_M;
int *spx_rowidx, *spx_colidx, *glp_rowidx;
int shuffle=1;
int feasible=0;

#define M(r,c)	spx_M[(r)+((c)-1)*rows-1] /* as indices go from 1..max */

void read_file(char *fname)
{FILE *f; int i,j,t; double v;
    f=fopen(fname,"r");
    if(f==NULL){
        fprintf(stderr,"Cannot open file %s for coefficients\n",fname);
        exit(1);
    }
    if(3!=fscanf(f,"%d%d%lg",&cols,&rows,&V)){
        fprintf(stderr,"Wrong data file %s (first 3 values: cols/rows/V)\n",fname);
        exit(1);
    }
    spx_colidx=malloc((cols+1)*sizeof(int));
    spx_rowidx=malloc((rows+1)*sizeof(int));
    spx_M=malloc(cols*rows*sizeof(double));
    glp_rowidx=malloc((rows+1)*sizeof(int));
    if(!spx_M || !spx_colidx || !spx_rowidx || !glp_rowidx){
        fprintf(stderr,"Out of memory...\n"); exit(1);
    }
    /* put a random permutation ... */
    for(j=0;j<=cols;j++) spx_colidx[j]=j;
    for(i=0;i<=rows;i++) spx_rowidx[i]=i;
    for(i=0;i<=rows;i++) glp_rowidx[i]=i;
    if(shuffle){ /* idx[0]=0 */
        for(i=1;i<cols;i++){
            j=i+mrandom(1+cols-i);
            t=spx_colidx[i]; spx_colidx[i]=spx_colidx[j];
            spx_colidx[j]=t;
        }
        for(j=1;j<rows;j++){
            i=j+mrandom(1+rows-j);
            t=spx_rowidx[j]; spx_rowidx[j]=spx_rowidx[i];
            spx_rowidx[i]=t;
        }
    }
    /* initialize glpk structure */
    P = glp_create_prob();
    glp_add_cols(P,cols); glp_add_rows(P,rows);
    /* each column is >=0 */
    for(j=1;j<=cols;j++) glp_set_col_bnds(P,j,GLP_LO,0.0,0.0);
    /* each row is == */
    for(i=1;i<=rows;i++)for(j=1;j<=cols;j++){
        if(1!=fscanf(f,"%lg",&M(spx_rowidx[i],spx_colidx[j]))){
            fprintf(stderr,"Cannot read M[%d][%d]\n",i,j); exit(1);
        }
    }
    /* add the matrix to the glp problem */
    for(j=1;j<=cols;j++){
        glp_set_mat_col(P,j,rows,glp_rowidx,&M(1,j)-1);
    }
    /* each row might start with "g" or "l", which changes the row stuatus */
    for(i=1;i<=rows;i++){ int dir=GLP_FX; // GLP_LO, GLP_UP
        if(1!=fscanf(f,"%lg",&v)){int c;
            while((c=fgetc(f))==' '|| c=='\n'||c=='\r'); // skip spaces
            if(c=='g'||c=='G'){ dir=GLP_LO; }
            else if(c=='l'||c=='L'){ dir=GLP_UP; }
            else { fprintf(stderr,"Cannot read c[%d] (char=%c)\n",i,c); exit(1); }
            if(1!=fscanf(f,"%lg",&v)){
               fprintf(stderr,"Cannot read c[%d] (number)\n",i); exit(1);
            }
        }
        glp_set_row_bnds(P,spx_rowidx[i],dir,v,v);
    }
    /* set the objective */
    if(feasible){
       for(j=1;j<=cols;j++){
          glp_set_obj_coef(P,j,0.0);
       }
    } else {
        for(j=1;j<=cols;j++){
            if(1!=fscanf(f,"%lg",&v)){
                fprintf(stderr,"Cannot read b[%d]\n",j); exit(1);
            }
            glp_set_obj_coef(P,spx_colidx[j],-v);
        }
    }
    free(glp_rowidx); free(spx_M);
    free(spx_rowidx); /* maybe keep these for result */
    fclose(f);
}
/*-------------------------------------------------------------------------*/
static char *glp_status_msg(int stat)
{static char *statmsg[] = {
"solution is undefined",	// GLP_UNDEF
"solution is feasible",		// GLP_FEAS
"solution is infeasible",	// GLP_INFEAS
"no feasible solution exists",	// GLP_NOFEAS
"solution is optimal",		// GLP_OPT
"solution is unbounded",	// GLP_UNBND
};
    if(1<=stat && stat<=6) return statmsg[stat-1];
    return "unknown";
}

static char *glp_return_msg(int retval)
{static char *retmsg[] = {
"invalid basis",			// GLP_EBADB	 *
"singular matrix",			// GLP_ESING	 *
"ill-conditioned matrix",		// GLP_ECOND	 *
"invalid bounds",			// GLP_EBOUND
"solver failed",			// GLP_EFAIL	 *
"objective lower limit reached",	// GLP_EOBJLL
"objective upper limit reached",	// GLP_EOBJUL
"iteration limit exceeded",		// GLP_EITLIM	 *
"time limit exceeded",			// GLP_ETMLIM	 *
"no primal feasible solution",		// GLP_ENOPFS
"no dual feasible solution",		// GLP_ENODFS
"root LP optimum not provided",		// GLP_EROOT
"search terminated by application",	// GLP_ESTOP
"relative mip gap tolerance reached",	// GLP_EMIPGAP
"no primal/dual feasible solution",	// GLP_ENOFEAS
"no convergence",			// GLP_ENOCVG
"numerical instability",		// GLP_EINSTAB
"invalid data",				// GLP_EDATA
"result out of range",			// GLP_ERANGE
};
    if(1<=retval && retval<=0x13) return retmsg[retval-1];
    return "unknown";
}

/*-------------------------------------------------------------------------*/

FILE *resfile=NULL;


int frac(double d, int*n,int*m)
/* computes d=n/m */
{int i; double dd,ddi; int sign;
   sign=1; if(d<0){d=-d; sign=-1;}
   for(i=1;i<200;i++){ dd=i*d; ddi=dd-mfloor(dd+0.5);
       if(ddi<1e-7 && ddi>-1e-7){
            *m=i; *n=(int)(mfloor(dd+0.5)); if(sign<0) *n= -(*n);
            return 1;
       }
   }
   return 0;
}

void showfrac(double c)
{int n,m;
    if(frac(c,&n,&m)){
        if(m==1){ fprintf(resfile,"%d",n); }
        else { fprintf(resfile, "%d/%d",n,m); }
    } else {
        fprintf(resfile,"%18.14lf",c);
    }
}

void print_base(void)
{int i; double v;
    for(i=1;i<=cols;i++){
        v=glp_get_col_prim(P,spx_colidx[i]);
        if(v>1e-9||v<-1e-9){
           fprintf(resfile,"%d:",i-1); showfrac(v); fprintf(resfile,",");
        }
    }
    fprintf(resfile,"\n");
}

/*===============================================================*/
#include <time.h>

int main(int argc, char *argv[])
{int ret=0;
    shuffle=1; feasible=0;
again:    
    if(argc>1 && argv[1][0]=='-'){
        if(argv[1][1]=='x'){
           feasible=1;
        } else if(argv[1][1]=='n'){
           shuffle=0;
        } else {
           fprintf(stderr,"unknown first argument %s\n",argv[1]); return 1;
        }
        argc--; argv++;
        goto again;
    }
    if(argc<2){
        fprintf(stderr, "usage: %s [-noshuffle | -x] <data> <output>\n",argv[0]);
        return 1;
    }
    if(shuffle) srand(time(NULL));
    read_file(argv[1]);
    // set glp parameters
    glp_init_smcp(&parm);
    parm.meth = GLP_PRIMAL;	// DUAL,DUALP
    parm.msg_lev = GLP_MSG_ERR; // OFF,ON,ALL,ERR
    parm.pricing = GLP_PT_PSE;  // PSE,STD
    parm.r_test = GLP_RT_STD;   // HAR,STD
//    parm.r_test = GLP_RT_HAR;	// HAR,STD
    parm.it_lim = 8000000;	// iteration limit
    parm.tm_lim = 1000000;	// time limit in milliseconds
    parm.out_frq =500;		// output frquency in interation numbers
//    parm.presolve = GLP_ON;	// try presolve
    // make initial bases
    glp_set_obj_dir(P,GLP_MIN); // minimize
    glp_term_out(GLP_OFF);	// don't print out messages
    glp_adv_basis(P,0);
    ret=glp_simplex(P,&parm);
    if(ret==0 && glp_get_status(P)==GLP_OPT){ // well done
        resfile=fopen(argv[2],"w"); if(!resfile){
            fprintf(stderr,"Cannot open resfile %s\n",argv[2]); exit(1);
        }
        fprintf(resfile,"round=%d, V=",glp_get_it_cnt(P));
          showfrac(V+glp_get_obj_val(P)); fprintf(resfile,"\n");
        print_base();
        fclose(resfile);
        return 0;
    }
    if(feasible!=0 && ret==0 && glp_get_status(P)==GLP_NOFEAS){
        resfile=fopen(argv[2],"w"); if(!resfile){
            fprintf(stderr,"Cannot open resfile %s\n",argv[2]); exit(1);
        }
        fprintf(resfile,"round=%d, V=XXXX\n",glp_get_it_cnt(P));
        fclose(resfile);
        return 0;
    }
    // some error 
    if(ret){
        fprintf(stderr,"The glp solver says: %s (%d)\n",glp_return_msg(ret),ret);
    } else {
        int status=glp_get_status(P);
        fprintf(stderr,"The glp solver says: %s (%d)\n",glp_status_msg(status),status);
    }
    return 0;
}

