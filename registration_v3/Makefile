DEBUG= -g
#OPT= -O4
INCLUDES= -I. 
CFLAGS= `pkg-config --cflags opencv` $(INCLUDES) $(OPT) $(DEBUG)
CPPFLAGS=$(CFLAGS)
LIBS=`pkg-config --libs opencv` -lm -ljpeg -lpthread
CPP=g++
CC=gcc

OBJS=orderedpointset.o unorderedpointset.o symmetry.o opt.o point3d.o rangeimage.o nose.o filter.o holes.o register.o main.o roi.o coordlist.o spline.o

#programs=register abs2wrl addglasses sfidiff transformabs abs2range maskglasses
programs=register abs2wrl absscale sfidiff sfi2pgm

#all: register abs2wrl addglasses sfidiff transformabs abs2range maskglasses
#all: register sficurv pnm2sfi sfidiff
all: register abs2wrl absscale sfidiff sfi2pgm

.c.o: 
	$(CC) $(CFLAGS) -c $<

.cc.o: 
	$(CPP) $(CPPFLAGS) -c $<

register: $(OBJS)
	g++ -o $@ $(OBJS) $(LIBS)
	
sfidiff: sfidiff.o rangeimage.o point3d.o unorderedpointset.o
	g++ -o $@ sfidiff.o rangeimage.o point3d.o unorderedpointset.o $(LIBS)
        
absscale: absscale.o orderedpointset.o
	g++ -o $@ absscale.o orderedpointset.o $(LIBS)
        
abs2wrl: abs2wrl.o orderedpointset.o
	g++ -o $@ abs2wrl.o orderedpointset.o $(LIBS)
	
sfi2pgm: sfi2pgm.o 
	g++ -o $@ sfi2pgm.o $(LIBS)
	
abs2range: abs2range.o orderedpointset.o unorderedpointset.o point3d.o
	g++ -o $@ abs2range.o orderedpointset.o unorderedpointset.o point3d.o $(LIBS)
	
maskglasses: maskglasses.o rangeimage.o point3d.o
	g++ -o $@ maskglasses.o rangeimage.o point3d.o $(LIBS)
	
addglasses: addglasses.o rangeimage.o point3d.o
	g++ -o $@ addglasses.o rangeimage.o point3d.o $(LIBS)
	
sfidiff: sfidiff.o rangeimage.o point3d.o unorderedpointset.o
	g++ -o $@ sfidiff.o rangeimage.o point3d.o unorderedpointset.o $(LIBS)
	
pnm2sfi: pnm2sfi.o rangeimage.o point3d.o unorderedpointset.o
	g++ -o $@ pnm2sfi.o rangeimage.o point3d.o unorderedpointset.o $(LIBS)
	
sficurv: sficurv.o rangeimage.o point3d.o unorderedpointset.o
	g++ -o $@ sficurv.o rangeimage.o point3d.o unorderedpointset.o $(LIBS)
	
transformabs: unorderedpointset.o transformabs.o rangeimage.o orderedpointset.o point3d.o
	g++ -o $@ unorderedpointset.o transformabs.o orderedpointset.o rangeimage.o point3d.o $(LIBS)

install: all
	cp $(programs) ../bin

clean:
	rm -f *.o

depend:
	makedepend $(INCLUDES) *.cc *.h

# DO NOT DELETE

abs2wrl.o: /usr/include/stdio.h /usr/include/features.h
abs2wrl.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
abs2wrl.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
abs2wrl.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
abs2wrl.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
abs2wrl.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
abs2wrl.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
abs2wrl.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
abs2wrl.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
abs2wrl.o: /usr/include/bits/byteswap-16.h /usr/include/sys/types.h
abs2wrl.o: /usr/include/time.h /usr/include/sys/select.h
abs2wrl.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
abs2wrl.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
abs2wrl.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
abs2wrl.o: /usr/include/bits/stdlib-float.h /usr/include/unistd.h
abs2wrl.o: /usr/include/bits/posix_opt.h /usr/include/bits/environments.h
abs2wrl.o: /usr/include/bits/confname.h /usr/include/getopt.h
abs2wrl.o: /usr/include/string.h /usr/include/xlocale.h /usr/include/math.h
abs2wrl.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
abs2wrl.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
abs2wrl.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
abs2wrl.o: /usr/include/bits/mathcalls.h orderedpointset.h point3d.h
abs2wrl.o: plane3d.h
coordlist.o: /usr/include/stdio.h /usr/include/features.h
coordlist.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
coordlist.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
coordlist.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
coordlist.o: /usr/include/libio.h /usr/include/_G_config.h
coordlist.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
coordlist.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
coordlist.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
coordlist.o: /usr/include/endian.h /usr/include/bits/endian.h
coordlist.o: /usr/include/bits/byteswap.h /usr/include/bits/byteswap-16.h
coordlist.o: /usr/include/sys/types.h /usr/include/time.h
coordlist.o: /usr/include/sys/select.h /usr/include/bits/select.h
coordlist.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
coordlist.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
coordlist.o: /usr/include/alloca.h /usr/include/bits/stdlib-float.h
coordlist.o: /usr/include/string.h /usr/include/xlocale.h /usr/include/math.h
coordlist.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
coordlist.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
coordlist.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
coordlist.o: /usr/include/bits/mathcalls.h compat.h /usr/include/sys/times.h
coordlist.o: coordlist.h
filter.o: /usr/include/stdio.h /usr/include/features.h
filter.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
filter.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
filter.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
filter.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
filter.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
filter.o: /usr/include/math.h /usr/include/bits/huge_val.h
filter.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
filter.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
filter.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h compat.h
filter.o: /usr/include/time.h /usr/include/sys/times.h roi.h
filter.o: unorderedpointset.h point3d.h rangeimage.h /usr/include/stdlib.h
filter.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
filter.o: /usr/include/endian.h /usr/include/bits/endian.h
filter.o: /usr/include/bits/byteswap.h /usr/include/bits/byteswap-16.h
filter.o: /usr/include/sys/types.h /usr/include/sys/select.h
filter.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
filter.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
filter.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
filter.o: /usr/include/bits/stdlib-float.h /usr/include/string.h
filter.o: /usr/include/xlocale.h filter.h
holes.o: /usr/include/stdio.h /usr/include/features.h
holes.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
holes.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
holes.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
holes.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
holes.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
holes.o: holes.h rangeimage.h /usr/include/stdlib.h
holes.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
holes.o: /usr/include/endian.h /usr/include/bits/endian.h
holes.o: /usr/include/bits/byteswap.h /usr/include/bits/byteswap-16.h
holes.o: /usr/include/sys/types.h /usr/include/time.h
holes.o: /usr/include/sys/select.h /usr/include/bits/select.h
holes.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
holes.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
holes.o: /usr/include/alloca.h /usr/include/bits/stdlib-float.h
holes.o: /usr/include/string.h /usr/include/xlocale.h compat.h
holes.o: /usr/include/sys/times.h unorderedpointset.h point3d.h
holes.o: /usr/include/math.h /usr/include/bits/huge_val.h
holes.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
holes.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
holes.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h nose.h
holes.o: opt.h
main.o: /usr/include/stdio.h /usr/include/features.h
main.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
main.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
main.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
main.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
main.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
main.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
main.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
main.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
main.o: /usr/include/bits/byteswap-16.h /usr/include/sys/types.h
main.o: /usr/include/time.h /usr/include/sys/select.h
main.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
main.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
main.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
main.o: /usr/include/bits/stdlib-float.h /usr/include/string.h
main.o: /usr/include/xlocale.h /usr/include/math.h
main.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
main.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
main.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
main.o: /usr/include/bits/mathcalls.h point3d.h plane3d.h orderedpointset.h
main.o: unorderedpointset.h rangeimage.h compat.h /usr/include/sys/times.h
main.o: nose.h opt.h filter.h holes.h symmetry.h register.h roi.h
nose.o: compat.h /usr/include/time.h /usr/include/bits/types.h
nose.o: /usr/include/features.h /usr/include/stdc-predef.h
nose.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
nose.o: /usr/include/gnu/stubs.h /usr/include/bits/typesizes.h
nose.o: /usr/include/sys/times.h /usr/include/math.h
nose.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
nose.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
nose.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
nose.o: /usr/include/bits/mathcalls.h nose.h /usr/include/stdio.h
nose.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
nose.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h opt.h
nose.o: point3d.h unorderedpointset.h rangeimage.h /usr/include/stdlib.h
nose.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
nose.o: /usr/include/endian.h /usr/include/bits/endian.h
nose.o: /usr/include/bits/byteswap.h /usr/include/bits/byteswap-16.h
nose.o: /usr/include/sys/types.h /usr/include/sys/select.h
nose.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
nose.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
nose.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
nose.o: /usr/include/bits/stdlib-float.h /usr/include/string.h
nose.o: /usr/include/xlocale.h roi.h spline.h coordlist.h
opt.o: /usr/include/stdlib.h /usr/include/features.h
opt.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
opt.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
opt.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
opt.o: /usr/include/endian.h /usr/include/bits/endian.h
opt.o: /usr/include/bits/byteswap.h /usr/include/bits/types.h
opt.o: /usr/include/bits/typesizes.h /usr/include/bits/byteswap-16.h
opt.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h
opt.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
opt.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
opt.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
opt.o: /usr/include/bits/stdlib-float.h /usr/include/stdio.h
opt.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
opt.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
opt.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
opt.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
opt.o: /usr/include/bits/posix2_lim.h /usr/include/math.h
opt.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
opt.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
opt.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
opt.o: /usr/include/bits/mathcalls.h opt.h
orderedpointset.o: /usr/include/stdio.h /usr/include/features.h
orderedpointset.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
orderedpointset.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
orderedpointset.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
orderedpointset.o: /usr/include/libio.h /usr/include/_G_config.h
orderedpointset.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
orderedpointset.o: /usr/include/bits/sys_errlist.h /usr/include/unistd.h
orderedpointset.o: /usr/include/bits/posix_opt.h
orderedpointset.o: /usr/include/bits/environments.h
orderedpointset.o: /usr/include/bits/confname.h /usr/include/getopt.h
orderedpointset.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
orderedpointset.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
orderedpointset.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
orderedpointset.o: /usr/include/bits/byteswap-16.h /usr/include/sys/types.h
orderedpointset.o: /usr/include/time.h /usr/include/sys/select.h
orderedpointset.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
orderedpointset.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
orderedpointset.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
orderedpointset.o: /usr/include/bits/stdlib-float.h /usr/include/string.h
orderedpointset.o: /usr/include/xlocale.h /usr/include/ctype.h
orderedpointset.o: /usr/include/math.h /usr/include/bits/huge_val.h
orderedpointset.o: /usr/include/bits/huge_valf.h
orderedpointset.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
orderedpointset.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
orderedpointset.o: /usr/include/bits/mathcalls.h compat.h
orderedpointset.o: /usr/include/sys/times.h orderedpointset.h point3d.h
orderedpointset.o: plane3d.h
pnm2sfi.o: /usr/include/stdio.h /usr/include/features.h
pnm2sfi.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
pnm2sfi.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
pnm2sfi.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
pnm2sfi.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
pnm2sfi.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
pnm2sfi.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
pnm2sfi.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
pnm2sfi.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
pnm2sfi.o: /usr/include/bits/byteswap-16.h /usr/include/sys/types.h
pnm2sfi.o: /usr/include/time.h /usr/include/sys/select.h
pnm2sfi.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
pnm2sfi.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
pnm2sfi.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
pnm2sfi.o: /usr/include/bits/stdlib-float.h /usr/include/unistd.h
pnm2sfi.o: /usr/include/bits/posix_opt.h /usr/include/bits/environments.h
pnm2sfi.o: /usr/include/bits/confname.h /usr/include/getopt.h rangeimage.h
pnm2sfi.o: /usr/include/string.h /usr/include/xlocale.h compat.h
pnm2sfi.o: /usr/include/sys/times.h unorderedpointset.h point3d.h
pnm2sfi.o: /usr/include/math.h /usr/include/bits/huge_val.h
pnm2sfi.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
pnm2sfi.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
pnm2sfi.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
point3d.o: /usr/include/stdio.h /usr/include/features.h
point3d.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
point3d.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
point3d.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
point3d.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
point3d.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
point3d.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
point3d.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
point3d.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
point3d.o: /usr/include/bits/byteswap-16.h /usr/include/sys/types.h
point3d.o: /usr/include/time.h /usr/include/sys/select.h
point3d.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
point3d.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
point3d.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
point3d.o: /usr/include/bits/stdlib-float.h /usr/include/string.h
point3d.o: /usr/include/xlocale.h /usr/include/math.h
point3d.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
point3d.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
point3d.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
point3d.o: /usr/include/bits/mathcalls.h point3d.h
rangeimage.o: /usr/include/stdio.h /usr/include/features.h
rangeimage.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
rangeimage.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
rangeimage.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
rangeimage.o: /usr/include/libio.h /usr/include/_G_config.h
rangeimage.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
rangeimage.o: /usr/include/bits/sys_errlist.h /usr/include/jpeglib.h
#rangeimage.o: /usr/include/jconfig.h /usr/include/jmorecfg.h
rangeimage.o: /usr/include/x86_64-linux-gnu/jconfig.h /usr/include/jmorecfg.h
rangeimage.o: /usr/include/setjmp.h /usr/include/bits/setjmp.h
rangeimage.o: /usr/include/bits/sigset.h /usr/include/X11/Xarch.h
rangeimage.o: rangeimage.h /usr/include/stdlib.h
rangeimage.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
rangeimage.o: /usr/include/endian.h /usr/include/bits/endian.h
rangeimage.o: /usr/include/bits/byteswap.h /usr/include/bits/byteswap-16.h
rangeimage.o: /usr/include/sys/types.h /usr/include/time.h
rangeimage.o: /usr/include/sys/select.h /usr/include/bits/select.h
rangeimage.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
rangeimage.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
rangeimage.o: /usr/include/bits/stdlib-float.h /usr/include/string.h
rangeimage.o: /usr/include/xlocale.h compat.h /usr/include/sys/times.h
rangeimage.o: unorderedpointset.h point3d.h /usr/include/math.h
rangeimage.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
rangeimage.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
rangeimage.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
rangeimage.o: /usr/include/bits/mathcalls.h
register.o: /usr/include/time.h /usr/include/bits/types.h
register.o: /usr/include/features.h /usr/include/stdc-predef.h
register.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
register.o: /usr/include/gnu/stubs.h /usr/include/bits/typesizes.h
register.o: /usr/include/stdio.h /usr/include/libio.h
register.o: /usr/include/_G_config.h /usr/include/wchar.h
register.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
register.o: opt.h /usr/include/math.h /usr/include/bits/huge_val.h
register.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
register.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
register.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
register.o: rangeimage.h /usr/include/stdlib.h /usr/include/bits/waitflags.h
register.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
register.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
register.o: /usr/include/bits/byteswap-16.h /usr/include/sys/types.h
register.o: /usr/include/sys/select.h /usr/include/bits/select.h
register.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
register.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
register.o: /usr/include/alloca.h /usr/include/bits/stdlib-float.h
register.o: /usr/include/string.h /usr/include/xlocale.h compat.h
register.o: /usr/include/sys/times.h unorderedpointset.h point3d.h register.h
register.o: plane3d.h filter.h nose.h symmetry.h holes.h roi.h
roi.o: unorderedpointset.h point3d.h /usr/include/math.h
roi.o: /usr/include/features.h /usr/include/stdc-predef.h
roi.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
roi.o: /usr/include/gnu/stubs.h /usr/include/bits/huge_val.h
roi.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
roi.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
roi.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h rangeimage.h
roi.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
roi.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
roi.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
roi.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
roi.o: /usr/include/bits/byteswap-16.h /usr/include/sys/types.h
roi.o: /usr/include/time.h /usr/include/sys/select.h
roi.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
roi.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
roi.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
roi.o: /usr/include/bits/stdlib-float.h /usr/include/stdio.h
roi.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
roi.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
roi.o: /usr/include/string.h /usr/include/xlocale.h compat.h
roi.o: /usr/include/sys/times.h roi.h
sficurv.o: /usr/include/stdio.h /usr/include/features.h
sficurv.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
sficurv.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
sficurv.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
sficurv.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
sficurv.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
sficurv.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
sficurv.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
sficurv.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
sficurv.o: /usr/include/bits/byteswap-16.h /usr/include/sys/types.h
sficurv.o: /usr/include/time.h /usr/include/sys/select.h
sficurv.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
sficurv.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
sficurv.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
sficurv.o: /usr/include/bits/stdlib-float.h /usr/include/unistd.h
sficurv.o: /usr/include/bits/posix_opt.h /usr/include/bits/environments.h
sficurv.o: /usr/include/bits/confname.h /usr/include/getopt.h rangeimage.h
sficurv.o: /usr/include/string.h /usr/include/xlocale.h compat.h
sficurv.o: /usr/include/sys/times.h unorderedpointset.h point3d.h
sficurv.o: /usr/include/math.h /usr/include/bits/huge_val.h
sficurv.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
sficurv.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
sficurv.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
sfidiff.o: /usr/include/stdio.h /usr/include/features.h
sfidiff.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
sfidiff.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
sfidiff.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
sfidiff.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
sfidiff.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
sfidiff.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
sfidiff.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
sfidiff.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
sfidiff.o: /usr/include/bits/byteswap-16.h /usr/include/sys/types.h
sfidiff.o: /usr/include/time.h /usr/include/sys/select.h
sfidiff.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
sfidiff.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
sfidiff.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
sfidiff.o: /usr/include/bits/stdlib-float.h /usr/include/unistd.h
sfidiff.o: /usr/include/bits/posix_opt.h /usr/include/bits/environments.h
sfidiff.o: /usr/include/bits/confname.h /usr/include/getopt.h rangeimage.h
sfidiff.o: /usr/include/string.h /usr/include/xlocale.h compat.h
sfidiff.o: /usr/include/sys/times.h unorderedpointset.h point3d.h
sfidiff.o: /usr/include/math.h /usr/include/bits/huge_val.h
sfidiff.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
sfidiff.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
sfidiff.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
spline.o: /usr/include/stdio.h /usr/include/features.h
spline.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
spline.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
spline.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
spline.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
spline.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
spline.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
spline.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
spline.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
spline.o: /usr/include/bits/byteswap-16.h /usr/include/sys/types.h
spline.o: /usr/include/time.h /usr/include/sys/select.h
spline.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
spline.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
spline.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
spline.o: /usr/include/bits/stdlib-float.h /usr/include/string.h
spline.o: /usr/include/xlocale.h /usr/include/math.h
spline.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
spline.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
spline.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
spline.o: /usr/include/bits/mathcalls.h compat.h /usr/include/sys/times.h
spline.o: spline.h coordlist.h
symmetry.o: /usr/include/stdio.h /usr/include/features.h
symmetry.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
symmetry.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
symmetry.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
symmetry.o: /usr/include/libio.h /usr/include/_G_config.h
symmetry.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
symmetry.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
symmetry.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
symmetry.o: /usr/include/endian.h /usr/include/bits/endian.h
symmetry.o: /usr/include/bits/byteswap.h /usr/include/bits/byteswap-16.h
symmetry.o: /usr/include/sys/types.h /usr/include/time.h
symmetry.o: /usr/include/sys/select.h /usr/include/bits/select.h
symmetry.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
symmetry.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
symmetry.o: /usr/include/alloca.h /usr/include/bits/stdlib-float.h
symmetry.o: /usr/include/pthread.h /usr/include/sched.h
symmetry.o: /usr/include/bits/sched.h /usr/include/bits/setjmp.h
symmetry.o: unorderedpointset.h point3d.h /usr/include/math.h
symmetry.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
symmetry.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
symmetry.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
symmetry.o: /usr/include/bits/mathcalls.h rangeimage.h /usr/include/string.h
symmetry.o: /usr/include/xlocale.h compat.h /usr/include/sys/times.h nose.h
symmetry.o: opt.h symmetry.h
transformabs.o: /usr/include/stdio.h /usr/include/features.h
transformabs.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
transformabs.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
transformabs.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
transformabs.o: /usr/include/libio.h /usr/include/_G_config.h
transformabs.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
transformabs.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
transformabs.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
transformabs.o: /usr/include/endian.h /usr/include/bits/endian.h
transformabs.o: /usr/include/bits/byteswap.h /usr/include/bits/byteswap-16.h
transformabs.o: /usr/include/sys/types.h /usr/include/time.h
transformabs.o: /usr/include/sys/select.h /usr/include/bits/select.h
transformabs.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
transformabs.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
transformabs.o: /usr/include/alloca.h /usr/include/bits/stdlib-float.h
transformabs.o: /usr/include/unistd.h /usr/include/bits/posix_opt.h
transformabs.o: /usr/include/bits/environments.h /usr/include/bits/confname.h
transformabs.o: /usr/include/getopt.h point3d.h /usr/include/math.h
transformabs.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
transformabs.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
transformabs.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
transformabs.o: /usr/include/bits/mathcalls.h orderedpointset.h plane3d.h
transformabs.o: rangeimage.h /usr/include/string.h /usr/include/xlocale.h
transformabs.o: compat.h /usr/include/sys/times.h unorderedpointset.h
tst.o: /usr/include/stdio.h /usr/include/features.h
tst.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
tst.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
tst.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
tst.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
tst.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h spline.h
tst.o: coordlist.h
unorderedpointset.o: /usr/include/stdio.h /usr/include/features.h
unorderedpointset.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
unorderedpointset.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
unorderedpointset.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
unorderedpointset.o: /usr/include/libio.h /usr/include/_G_config.h
unorderedpointset.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
unorderedpointset.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
unorderedpointset.o: /usr/include/bits/waitflags.h
unorderedpointset.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
unorderedpointset.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
unorderedpointset.o: /usr/include/bits/byteswap-16.h /usr/include/sys/types.h
unorderedpointset.o: /usr/include/time.h /usr/include/sys/select.h
unorderedpointset.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
unorderedpointset.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
unorderedpointset.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
unorderedpointset.o: /usr/include/bits/stdlib-float.h /usr/include/string.h
unorderedpointset.o: /usr/include/xlocale.h /usr/include/math.h
unorderedpointset.o: /usr/include/bits/huge_val.h
unorderedpointset.o: /usr/include/bits/huge_valf.h
unorderedpointset.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
unorderedpointset.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
unorderedpointset.o: /usr/include/bits/mathcalls.h compat.h
unorderedpointset.o: /usr/include/sys/times.h unorderedpointset.h point3d.h
compat.o: /usr/include/time.h /usr/include/bits/types.h
compat.o: /usr/include/features.h /usr/include/stdc-predef.h
compat.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
compat.o: /usr/include/gnu/stubs.h /usr/include/bits/typesizes.h
compat.o: /usr/include/sys/times.h
filter.o: rangeimage.h /usr/include/stdlib.h /usr/include/features.h
filter.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
filter.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
filter.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
filter.o: /usr/include/endian.h /usr/include/bits/endian.h
filter.o: /usr/include/bits/byteswap.h /usr/include/bits/types.h
filter.o: /usr/include/bits/typesizes.h /usr/include/bits/byteswap-16.h
filter.o: /usr/include/sys/types.h /usr/include/time.h
filter.o: /usr/include/sys/select.h /usr/include/bits/select.h
filter.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
filter.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
filter.o: /usr/include/alloca.h /usr/include/bits/stdlib-float.h
filter.o: /usr/include/stdio.h /usr/include/libio.h /usr/include/_G_config.h
filter.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
filter.o: /usr/include/bits/sys_errlist.h /usr/include/string.h
filter.o: /usr/include/xlocale.h compat.h /usr/include/sys/times.h
filter.o: unorderedpointset.h point3d.h /usr/include/math.h
filter.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
filter.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
filter.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
filter.o: /usr/include/bits/mathcalls.h
holes.o: rangeimage.h /usr/include/stdlib.h /usr/include/features.h
holes.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
holes.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
holes.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
holes.o: /usr/include/endian.h /usr/include/bits/endian.h
holes.o: /usr/include/bits/byteswap.h /usr/include/bits/types.h
holes.o: /usr/include/bits/typesizes.h /usr/include/bits/byteswap-16.h
holes.o: /usr/include/sys/types.h /usr/include/time.h
holes.o: /usr/include/sys/select.h /usr/include/bits/select.h
holes.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
holes.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
holes.o: /usr/include/alloca.h /usr/include/bits/stdlib-float.h
holes.o: /usr/include/stdio.h /usr/include/libio.h /usr/include/_G_config.h
holes.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
holes.o: /usr/include/bits/sys_errlist.h /usr/include/string.h
holes.o: /usr/include/xlocale.h compat.h /usr/include/sys/times.h
holes.o: unorderedpointset.h point3d.h /usr/include/math.h
holes.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
holes.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
holes.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
holes.o: /usr/include/bits/mathcalls.h nose.h opt.h
nose.o: /usr/include/stdio.h /usr/include/features.h
nose.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
nose.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
nose.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
nose.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
nose.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h opt.h
nose.o: /usr/include/math.h /usr/include/bits/huge_val.h
nose.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
nose.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
nose.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h point3d.h
nose.o: unorderedpointset.h rangeimage.h /usr/include/stdlib.h
nose.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
nose.o: /usr/include/endian.h /usr/include/bits/endian.h
nose.o: /usr/include/bits/byteswap.h /usr/include/bits/byteswap-16.h
nose.o: /usr/include/sys/types.h /usr/include/time.h
nose.o: /usr/include/sys/select.h /usr/include/bits/select.h
nose.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
nose.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
nose.o: /usr/include/alloca.h /usr/include/bits/stdlib-float.h
nose.o: /usr/include/string.h /usr/include/xlocale.h compat.h
nose.o: /usr/include/sys/times.h
opt.o: /usr/include/math.h /usr/include/features.h /usr/include/stdc-predef.h
opt.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
opt.o: /usr/include/gnu/stubs.h /usr/include/bits/huge_val.h
opt.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
opt.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
opt.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
orderedpointset.o: point3d.h /usr/include/math.h /usr/include/features.h
orderedpointset.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
orderedpointset.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
orderedpointset.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
orderedpointset.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
orderedpointset.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
orderedpointset.o: /usr/include/bits/mathcalls.h plane3d.h
point3d.o: /usr/include/math.h /usr/include/features.h
point3d.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
point3d.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
point3d.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
point3d.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
point3d.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
point3d.o: /usr/include/bits/mathcalls.h
rangeimage.o: /usr/include/stdlib.h /usr/include/features.h
rangeimage.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
rangeimage.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
rangeimage.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
rangeimage.o: /usr/include/endian.h /usr/include/bits/endian.h
rangeimage.o: /usr/include/bits/byteswap.h /usr/include/bits/types.h
rangeimage.o: /usr/include/bits/typesizes.h /usr/include/bits/byteswap-16.h
rangeimage.o: /usr/include/sys/types.h /usr/include/time.h
rangeimage.o: /usr/include/sys/select.h /usr/include/bits/select.h
rangeimage.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
rangeimage.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
rangeimage.o: /usr/include/alloca.h /usr/include/bits/stdlib-float.h
rangeimage.o: /usr/include/stdio.h /usr/include/libio.h
rangeimage.o: /usr/include/_G_config.h /usr/include/wchar.h
rangeimage.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
rangeimage.o: /usr/include/string.h /usr/include/xlocale.h compat.h
rangeimage.o: /usr/include/sys/times.h unorderedpointset.h point3d.h
rangeimage.o: /usr/include/math.h /usr/include/bits/huge_val.h
rangeimage.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
rangeimage.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
rangeimage.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
register.o: unorderedpointset.h point3d.h /usr/include/math.h
register.o: /usr/include/features.h /usr/include/stdc-predef.h
register.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
register.o: /usr/include/gnu/stubs.h /usr/include/bits/huge_val.h
register.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
register.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
register.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
register.o: rangeimage.h /usr/include/stdlib.h /usr/include/bits/waitflags.h
register.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
register.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
register.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
register.o: /usr/include/bits/byteswap-16.h /usr/include/sys/types.h
register.o: /usr/include/time.h /usr/include/sys/select.h
register.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
register.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
register.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
register.o: /usr/include/bits/stdlib-float.h /usr/include/stdio.h
register.o: /usr/include/libio.h /usr/include/_G_config.h
register.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
register.o: /usr/include/bits/sys_errlist.h /usr/include/string.h
register.o: /usr/include/xlocale.h compat.h /usr/include/sys/times.h
roi.o: unorderedpointset.h point3d.h /usr/include/math.h
roi.o: /usr/include/features.h /usr/include/stdc-predef.h
roi.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
roi.o: /usr/include/gnu/stubs.h /usr/include/bits/huge_val.h
roi.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
roi.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
roi.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h rangeimage.h
roi.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
roi.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
roi.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
roi.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
roi.o: /usr/include/bits/byteswap-16.h /usr/include/sys/types.h
roi.o: /usr/include/time.h /usr/include/sys/select.h
roi.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
roi.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
roi.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
roi.o: /usr/include/bits/stdlib-float.h /usr/include/stdio.h
roi.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
roi.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
roi.o: /usr/include/string.h /usr/include/xlocale.h compat.h
roi.o: /usr/include/sys/times.h
spline.o: coordlist.h
symmetry.o: /usr/include/stdio.h /usr/include/features.h
symmetry.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
symmetry.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
symmetry.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
symmetry.o: /usr/include/libio.h /usr/include/_G_config.h
symmetry.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
symmetry.o: /usr/include/bits/sys_errlist.h opt.h /usr/include/math.h
symmetry.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
symmetry.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
symmetry.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
symmetry.o: /usr/include/bits/mathcalls.h
unorderedpointset.o: point3d.h /usr/include/math.h /usr/include/features.h
unorderedpointset.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
unorderedpointset.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
unorderedpointset.o: /usr/include/bits/huge_val.h
unorderedpointset.o: /usr/include/bits/huge_valf.h
unorderedpointset.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
unorderedpointset.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
unorderedpointset.o: /usr/include/bits/mathcalls.h
