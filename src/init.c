#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(d2norm)(int*, double*, int*, double*);
extern void F77_NAME(es1e)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(es1v)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(eseee)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(eseei)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(eseev)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(eseii)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(eseve)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(esevi)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(esevv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(esvee)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(esvei)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(esvev)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(esvii)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(esvve)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(esvvi)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(esvvv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hc1e)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hc1v)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hceee)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hceii)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hcvii)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hcvvv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mcltrw)(void *, void *, void *, void *, void *);
extern void F77_NAME(me1e)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(me1ep)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(me1v)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(me1vp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(meeee)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(meeeep)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(meeei)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(meeeip)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(meeev)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(meeevp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(meeii)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(meeiip)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(meeve)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(meevi)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(meevip)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(meevv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mevee)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mevei)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(meveip)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mevev)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mevevp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mevii)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(meviip)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mevve)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mevvi)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mevvip)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mevvv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mevvvp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mnxiip)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mnxxip)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mnxxxp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ms1e)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ms1ep)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ms1v)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ms1vp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mseee)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mseeep)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mseei)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mseeip)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mseev)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mseevp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mseii)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mseiip)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mseve)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msevi)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msevip)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msevv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msvee)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msvei)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msveip)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msvev)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msvevp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msvii)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msviip)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msvve)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msvvi)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msvvip)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msvvv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msvvvp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mvn1d)(void *, void *, void *, void *, void *);
extern void F77_NAME(mvn1p)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mvnxii)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mvnxxi)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mvnxxx)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(shapeo)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(uncholf)(void *, void *, void *, void *, void *);//
extern void F77_NAME(covwf)(double*, double*, int*, int*, int*, double*, double*, double*);
extern void F77_NAME(crossprodf)(void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"d2norm",     (DL_FUNC) &F77_NAME(d2norm),     4},
    {"es1e",       (DL_FUNC) &F77_NAME(es1e),       9},
    {"es1v",       (DL_FUNC) &F77_NAME(es1v),       9},
    {"eseee",      (DL_FUNC) &F77_NAME(eseee),     12},
    {"eseei",      (DL_FUNC) &F77_NAME(eseei),     11},
    {"eseev",      (DL_FUNC) &F77_NAME(eseev),     14},
    {"eseii",      (DL_FUNC) &F77_NAME(eseii),     10},
    {"eseve",      (DL_FUNC) &F77_NAME(eseve),     14},
    {"esevi",      (DL_FUNC) &F77_NAME(esevi),     11},
    {"esevv",      (DL_FUNC) &F77_NAME(esevv),     14},
    {"esvee",      (DL_FUNC) &F77_NAME(esvee),     14},
    {"esvei",      (DL_FUNC) &F77_NAME(esvei),     11},
    {"esvev",      (DL_FUNC) &F77_NAME(esvev),     14},
    {"esvii",      (DL_FUNC) &F77_NAME(esvii),     10},
    {"esvve",      (DL_FUNC) &F77_NAME(esvve),     14},
    {"esvvi",      (DL_FUNC) &F77_NAME(esvvi),     11},
    {"esvvv",      (DL_FUNC) &F77_NAME(esvvv),     12},
    {"hc1e",       (DL_FUNC) &F77_NAME(hc1e),       7},
    {"hc1v",       (DL_FUNC) &F77_NAME(hc1v),       8},
    {"hceee",      (DL_FUNC) &F77_NAME(hceee),     12},
    {"hceii",      (DL_FUNC) &F77_NAME(hceii),      9},
    {"hcvii",      (DL_FUNC) &F77_NAME(hcvii),     10},
    {"hcvvv",      (DL_FUNC) &F77_NAME(hcvvv),     14},
    {"mcltrw",     (DL_FUNC) &F77_NAME(mcltrw),     5},
    {"me1e",       (DL_FUNC) &F77_NAME(me1e),      12},
    {"me1ep",      (DL_FUNC) &F77_NAME(me1ep),     16},
    {"me1v",       (DL_FUNC) &F77_NAME(me1v),      12},
    {"me1vp",      (DL_FUNC) &F77_NAME(me1vp),     16},
    {"meeee",      (DL_FUNC) &F77_NAME(meeee),     14},
    {"meeeep",     (DL_FUNC) &F77_NAME(meeeep),    18},
    {"meeei",      (DL_FUNC) &F77_NAME(meeei),     14},
    {"meeeip",     (DL_FUNC) &F77_NAME(meeeip),    18},
    {"meeev",      (DL_FUNC) &F77_NAME(meeev),     18},
    {"meeevp",     (DL_FUNC) &F77_NAME(meeevp),    22},
    {"meeii",      (DL_FUNC) &F77_NAME(meeii),     13},
    {"meeiip",     (DL_FUNC) &F77_NAME(meeiip),    17},
    {"meeve",      (DL_FUNC) &F77_NAME(meeve),     26},
    {"meevi",      (DL_FUNC) &F77_NAME(meevi),     14},
    {"meevip",     (DL_FUNC) &F77_NAME(meevip),    18},
    {"meevv",      (DL_FUNC) &F77_NAME(meevv),     22},
    {"mevee",      (DL_FUNC) &F77_NAME(mevee),     26},
    {"mevei",      (DL_FUNC) &F77_NAME(mevei),     17},
    {"meveip",     (DL_FUNC) &F77_NAME(meveip),    21},
    {"mevev",      (DL_FUNC) &F77_NAME(mevev),     18},
    {"mevevp",     (DL_FUNC) &F77_NAME(mevevp),    22},
    {"mevii",      (DL_FUNC) &F77_NAME(mevii),     13},
    {"meviip",     (DL_FUNC) &F77_NAME(meviip),    17},
    {"mevve",      (DL_FUNC) &F77_NAME(mevve),     26},
    {"mevvi",      (DL_FUNC) &F77_NAME(mevvi),     14},
    {"mevvip",     (DL_FUNC) &F77_NAME(mevvip),    18},
    {"mevvv",      (DL_FUNC) &F77_NAME(mevvv),     15},
    {"mevvvp",     (DL_FUNC) &F77_NAME(mevvvp),    19},
    {"mnxiip",     (DL_FUNC) &F77_NAME(mnxiip),    10},
    {"mnxxip",     (DL_FUNC) &F77_NAME(mnxxip),    11},
    {"mnxxxp",     (DL_FUNC) &F77_NAME(mnxxxp),    11},
    {"ms1e",       (DL_FUNC) &F77_NAME(ms1e),       7},
    {"ms1ep",      (DL_FUNC) &F77_NAME(ms1ep),     11},
    {"ms1v",       (DL_FUNC) &F77_NAME(ms1v),       7},
    {"ms1vp",      (DL_FUNC) &F77_NAME(ms1vp),     11},
    {"mseee",      (DL_FUNC) &F77_NAME(mseee),      9},
    {"mseeep",     (DL_FUNC) &F77_NAME(mseeep),    13},
    {"mseei",      (DL_FUNC) &F77_NAME(mseei),      9},
    {"mseeip",     (DL_FUNC) &F77_NAME(mseeip),    13},
    {"mseev",      (DL_FUNC) &F77_NAME(mseev),     12},
    {"mseevp",     (DL_FUNC) &F77_NAME(mseevp),    16},
    {"mseii",      (DL_FUNC) &F77_NAME(mseii),      8},
    {"mseiip",     (DL_FUNC) &F77_NAME(mseiip),    12},
    {"mseve",      (DL_FUNC) &F77_NAME(mseve),     18},
    {"msevi",      (DL_FUNC) &F77_NAME(msevi),      9},
    {"msevip",     (DL_FUNC) &F77_NAME(msevip),    13},
    {"msevv",      (DL_FUNC) &F77_NAME(msevv),     14},
    {"msvee",      (DL_FUNC) &F77_NAME(msvee),     17},
    {"msvei",      (DL_FUNC) &F77_NAME(msvei),     14},
    {"msveip",     (DL_FUNC) &F77_NAME(msveip),    18},
    {"msvev",      (DL_FUNC) &F77_NAME(msvev),     14},
    {"msvevp",     (DL_FUNC) &F77_NAME(msvevp),    18},
    {"msvii",      (DL_FUNC) &F77_NAME(msvii),      8},
    {"msviip",     (DL_FUNC) &F77_NAME(msviip),    12},
    {"msvve",      (DL_FUNC) &F77_NAME(msvve),     18},
    {"msvvi",      (DL_FUNC) &F77_NAME(msvvi),      9},
    {"msvvip",     (DL_FUNC) &F77_NAME(msvvip),    13},
    {"msvvv",      (DL_FUNC) &F77_NAME(msvvv),     10},
    {"msvvvp",     (DL_FUNC) &F77_NAME(msvvvp),    14},
    {"mvn1d",      (DL_FUNC) &F77_NAME(mvn1d),      5},
    {"mvn1p",      (DL_FUNC) &F77_NAME(mvn1p),      9},
    {"mvnxii",     (DL_FUNC) &F77_NAME(mvnxii),     6},
    {"mvnxxi",     (DL_FUNC) &F77_NAME(mvnxxi),     7},
    {"mvnxxx",     (DL_FUNC) &F77_NAME(mvnxxx),     6},
    {"shapeo",     (DL_FUNC) &F77_NAME(shapeo),     7},
    {"uncholf",    (DL_FUNC) &F77_NAME(uncholf),    5},
    {"covwf",      (DL_FUNC) &F77_NAME(covwf),      8},
    {"crossprodf", (DL_FUNC) &F77_NAME(crossprodf), 6}, //
    {NULL, NULL, 0}
};

void R_init_mclust(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
