/* PCG64 and binomial routines derived from NumPy 2.2.3 (BSD-3-Clause).
   Initial state is SeedSequence(20260829) -> PCG64, not saved random draws. */
#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define RAND_INT_TYPE int64_t
#define MIN(a,b) ((a)<(b)?(a):(b))
typedef struct {int has_binomial; double psave; int64_t nsave; double r,q,fm; int64_t m; double p1,xm,xl,xr,c,laml,lamr,p2,p3,p4;} binomial_t;
typedef struct {__uint128_t state,inc; int has32; uint32_t cached32; binomial_t binomial;} rng_state;
typedef struct {void* state; double (*next_double)(void*);} bitgen_t;
static uint64_t raw64(rng_state* r) {
 r->state=r->state*(((__uint128_t)2549297995355413924ULL<<64)|4865540595714422341ULL)+r->inc;
 uint64_t v=(uint64_t)(r->state>>64)^(uint64_t)r->state;unsigned rot=(unsigned)(r->state>>122);
 return (v>>rot)|(v<<((-rot)&63));
}
static uint32_t raw32(rng_state* r) {if(r->has32){r->has32=0;return r->cached32;}uint64_t v=raw64(r);r->has32=1;r->cached32=(uint32_t)(v>>32);return (uint32_t)v;}
static double uniform_double(void* s){return (raw64((rng_state*)s)>>11)*0x1.0p-53;}
static double next_double(bitgen_t* b){return b->next_double(b->state);}
RAND_INT_TYPE random_binomial_btpe(bitgen_t *bitgen_state, RAND_INT_TYPE n,
                                   double p, binomial_t *binomial) {
  double r, q, fm, p1, xm, xl, xr, c, laml, lamr, p2, p3, p4;
  double a, u, v, s, F, rho, t, A, nrq, x1, x2, f1, f2, z, z2, w, w2, x;
  RAND_INT_TYPE m, y, k, i;

  if (!(binomial->has_binomial) || (binomial->nsave != n) ||
      (binomial->psave != p)) {
    /* initialize */
    binomial->nsave = n;
    binomial->psave = p;
    binomial->has_binomial = 1;
    binomial->r = r = MIN(p, 1.0 - p);
    binomial->q = q = 1.0 - r;
    binomial->fm = fm = n * r + r;
    binomial->m = m = (RAND_INT_TYPE)floor(binomial->fm);
    binomial->p1 = p1 = floor(2.195 * sqrt(n * r * q) - 4.6 * q) + 0.5;
    binomial->xm = xm = m + 0.5;
    binomial->xl = xl = xm - p1;
    binomial->xr = xr = xm + p1;
    binomial->c = c = 0.134 + 20.5 / (15.3 + m);
    a = (fm - xl) / (fm - xl * r);
    binomial->laml = laml = a * (1.0 + a / 2.0);
    a = (xr - fm) / (xr * q);
    binomial->lamr = lamr = a * (1.0 + a / 2.0);
    binomial->p2 = p2 = p1 * (1.0 + 2.0 * c);
    binomial->p3 = p3 = p2 + c / laml;
    binomial->p4 = p4 = p3 + c / lamr;
  } else {
    r = binomial->r;
    q = binomial->q;
    fm = binomial->fm;
    m = binomial->m;
    p1 = binomial->p1;
    xm = binomial->xm;
    xl = binomial->xl;
    xr = binomial->xr;
    c = binomial->c;
    laml = binomial->laml;
    lamr = binomial->lamr;
    p2 = binomial->p2;
    p3 = binomial->p3;
    p4 = binomial->p4;
  }

/* sigh ... */
Step10:
  nrq = n * r * q;
  u = next_double(bitgen_state) * p4;
  v = next_double(bitgen_state);
  if (u > p1)
    goto Step20;
  y = (RAND_INT_TYPE)floor(xm - p1 * v + u);
  goto Step60;

Step20:
  if (u > p2)
    goto Step30;
  x = xl + (u - p1) / c;
  v = v * c + 1.0 - fabs(m - x + 0.5) / p1;
  if (v > 1.0)
    goto Step10;
  y = (RAND_INT_TYPE)floor(x);
  goto Step50;

Step30:
  if (u > p3)
    goto Step40;
  y = (RAND_INT_TYPE)floor(xl + log(v) / laml);
  /* Reject if v==0.0 since previous cast is undefined */
  if ((y < 0) || (v == 0.0))
    goto Step10;
  v = v * (u - p2) * laml;
  goto Step50;

Step40:
  y = (RAND_INT_TYPE)floor(xr - log(v) / lamr);
  /* Reject if v==0.0 since previous cast is undefined */
  if ((y > n) || (v == 0.0))
    goto Step10;
  v = v * (u - p3) * lamr;

Step50:
  k = llabs(y - m);
  if ((k > 20) && (k < ((nrq) / 2.0 - 1)))
    goto Step52;

  s = r / q;
  a = s * (n + 1);
  F = 1.0;
  if (m < y) {
    for (i = m + 1; i <= y; i++) {
      F *= (a / i - s);
    }
  } else if (m > y) {
    for (i = y + 1; i <= m; i++) {
      F /= (a / i - s);
    }
  }
  if (v > F)
    goto Step10;
  goto Step60;

Step52:
  rho =
      (k / (nrq)) * ((k * (k / 3.0 + 0.625) + 0.16666666666666666) / nrq + 0.5);
  t = -k * k / (2 * nrq);
  /* log(0.0) ok here */
  A = log(v);
  if (A < (t - rho))
    goto Step60;
  if (A > (t + rho))
    goto Step10;

  x1 = y + 1;
  f1 = m + 1;
  z = n + 1 - m;
  w = n - y + 1;
  x2 = x1 * x1;
  f2 = f1 * f1;
  z2 = z * z;
  w2 = w * w;
  if (A > (xm * log(f1 / x1) + (n - m + 0.5) * log(z / w) +
           (y - m) * log(w * r / (x1 * q)) +
           (13680. - (462. - (132. - (99. - 140. / f2) / f2) / f2) / f2) / f1 /
               166320. +
           (13680. - (462. - (132. - (99. - 140. / z2) / z2) / z2) / z2) / z /
               166320. +
           (13680. - (462. - (132. - (99. - 140. / x2) / x2) / x2) / x2) / x1 /
               166320. +
           (13680. - (462. - (132. - (99. - 140. / w2) / w2) / w2) / w2) / w /
               166320.)) {
    goto Step10;
  }

Step60:
  if (p > 0.5) {
    y = n - y;
  }

  return y;
}

RAND_INT_TYPE random_binomial_inversion(bitgen_t *bitgen_state, RAND_INT_TYPE n,
                                        double p, binomial_t *binomial) {
  double q, qn, np, px, U;
  RAND_INT_TYPE X, bound;

  if (!(binomial->has_binomial) || (binomial->nsave != n) ||
      (binomial->psave != p)) {
    binomial->nsave = n;
    binomial->psave = p;
    binomial->has_binomial = 1;
    binomial->q = q = 1.0 - p;
    binomial->r = qn = exp(n * log(q));
    binomial->c = np = n * p;
    binomial->m = bound = (RAND_INT_TYPE)MIN(n, np + 10.0 * sqrt(np * q + 1));
  } else {
    q = binomial->q;
    qn = binomial->r;
    np = binomial->c;
    bound = binomial->m;
  }
  X = 0;
  px = qn;
  U = next_double(bitgen_state);
  while (U > px) {
    X++;
    if (X > bound) {
      X = 0;
      px = qn;
      U = next_double(bitgen_state);
    } else {
      U -= px;
      px = ((n - X + 1) * p * px) / (X * q);
    }
  }
  return X;
}

int64_t random_binomial(bitgen_t *bitgen_state, double p, int64_t n,
                        binomial_t *binomial) {
  double q;

  if ((n == 0LL) || (p == 0.0f))
    return 0;

  if (p <= 0.5) {
    if (p * n <= 30.0) {
      return random_binomial_inversion(bitgen_state, n, p, binomial);
    } else {
      return random_binomial_btpe(bitgen_state, n, p, binomial);
    }
  } else {
    q = 1.0 - p;
    if (q * n <= 30.0) {
      return n - random_binomial_inversion(bitgen_state, n, q, binomial);
    } else {
      return n - random_binomial_btpe(bitgen_state, n, q, binomial);
    }
  }
}


static void finalize_rng(SEXP p){free(R_ExternalPtrAddr(p));R_ClearExternalPtr(p);}
SEXP ampath_rng_new(SEXP seed){if(asInteger(seed)!=20260829)error("Supported manuscript seed is 20260829");rng_state* r=calloc(1,sizeof(rng_state));if(!r)error("Allocation failed");
r->state=((__uint128_t)12290389383239710926ULL<<64)|2779326069082736607ULL;
r->inc=((__uint128_t)4148745151490402713ULL<<64)|13528396886862849435ULL;

SEXP p=PROTECT(R_MakeExternalPtr(r,R_NilValue,R_NilValue));R_RegisterCFinalizerEx(p,finalize_rng,TRUE);UNPROTECT(1);return p;}
static uint32_t bounded(rng_state* r,uint32_t n){if(n==1)return 0;uint64_t m=(uint64_t)raw32(r)*n;uint32_t lo=(uint32_t)m;if(lo<n){uint32_t threshold=(-n)%n;while(lo<threshold){m=(uint64_t)raw32(r)*n;lo=(uint32_t)m;}}return m>>32;}
SEXP ampath_rng_counts(SEXP pointer,SEXP size,SEXP reps,SEXP method){rng_state* r=R_ExternalPtrAddr(pointer);if(!r)error("Invalid RNG pointer");int n=asInteger(size),B=asInteger(reps),multi=asLogical(method);if(n<1||B<1)error("Positive sizes required");SEXP ans=PROTECT(allocMatrix(INTSXP,n,B));memset(INTEGER(ans),0,sizeof(int)*(size_t)n*B);bitgen_t b={r,uniform_double};
for(int k=0;k<B;k++){if(k%32==0)R_CheckUserInterrupt();int* out=INTEGER(ans)+(size_t)k*n;
if(!multi){for(int j=0;j<n;j++)out[bounded(r,n)]++;}
else {double remaining=1.0,p=1.0/n;int64_t dn=n;for(int j=0;j<n-1;j++){out[j]=(int)random_binomial(&b,p/remaining,dn,&r->binomial);dn-=out[j];if(dn<=0)break;remaining-=p;}if(dn>0)out[n-1]=(int)dn;}}
UNPROTECT(1);return ans;}
SEXP ampath_rng_integers(SEXP pointer,SEXP size,SEXP count){rng_state* r=R_ExternalPtrAddr(pointer);int n=asInteger(size),N=asInteger(count);SEXP ans=PROTECT(allocVector(INTSXP,N));for(int j=0;j<N;j++)INTEGER(ans)[j]=bounded(r,n);UNPROTECT(1);return ans;}
/* CPython 3.12 set(defaultdict) iteration, with incremental table growth. */
static void set_insert(uint64_t key,int index,uint64_t* table,int* indices,size_t mask){uint64_t hash=key%2305843009213693951ULL;size_t i=hash&mask;uint64_t perturb=hash;for(;;){if(!indices[i]){table[i]=key;indices[i]=index;return;}if(i+9<=mask){for(int k=1;k<=9;k++){if(!indices[i+k]){table[i+k]=key;indices[i+k]=index;return;}}}perturb>>=5;i=(i*5+1+perturb)&mask;}}
SEXP ampath_python_set_order(SEXP keys){int n=LENGTH(keys);size_t size=8;uint64_t* table=calloc(size,sizeof(uint64_t));int* indices=calloc(size,sizeof(int));if(!table||!indices)error("Allocation failed");
for(int j=0;j<n;j++){set_insert((uint64_t)REAL(keys)[j],j+1,table,indices,size-1);if((size_t)(j+1)*5>=(size-1)*3){size_t newsize=8,minused=(j+1)*((j+1)>50000?2:4);while(newsize<=minused)newsize<<=1;uint64_t* nt=calloc(newsize,sizeof(uint64_t));int* ni=calloc(newsize,sizeof(int));if(!nt||!ni)error("Allocation failed");for(size_t i=0;i<size;i++)if(indices[i])set_insert(table[i],indices[i],nt,ni,newsize-1);free(table);free(indices);table=nt;indices=ni;size=newsize;}}
SEXP ans=PROTECT(allocVector(INTSXP,n));int k=0;for(size_t i=0;i<size;i++)if(indices[i])INTEGER(ans)[k++]=indices[i];free(table);free(indices);UNPROTECT(1);return ans;}
