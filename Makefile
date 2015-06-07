PROG =	test_mlfv

SRCS =	mlfv_mod.f90 test_mlfv.f90 

OBJS =	mlfv_mod.o test_mlfv.o

OBJS_LIB = mlfv_mod.o

LIBS =  

F90 = gfortran
F90FLAGS = -O3 -g -Wall 
LDFLAGS =

#
# ... Cancel .mod.o rule
#
%.o : %.mod
# 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS) $(INCS)

clean:
	rm -f $(PROG) $(OBJS) *.mod *.a

install:
	mv -f $(PROG) /usr/local/bin
	
libmittag.a: $(OBJS)
	ar rc $@ $^ && ranlib $@

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

