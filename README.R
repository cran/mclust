The file blapack.f contains a number of functions from the LAPACK
and BLAS libraries, necessary for mclust.f. However, on SUN Solaris
systems and maybe some others, f77 will give problems 
(function i_dnnt is not defined; R CMD check will give an error).

A workaround is to use f2c, which worked for my system...
Issue the following commands in the src directory:

cat blapack.f | f2c > blapack.c
mv blapack.f blapack.f.bak
mv Makefile.SUN Makefile

Edit the makefile to indicate the location of libf2c.a 
(here: /usr/local/lib) and do a R INSTALL mclust at the top of the
mclust tree.

If you compiled R using BLAS, you might not even need blapack.f
