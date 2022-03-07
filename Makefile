ASAN = -static-libasan -fsanitize=address
CFLAGS = -Wall -Wextra -O3 -flto=auto -mtune=native $(ASAN)
LDFLAGS= -lgsl -lm -lblas -lgslcblas 
CC=gcc

MKLROOT=/opt/intel/lib
MKLLIB= -Wl,--start-group ${MKLROOT}/libmkl_intel_lp64.a ${MKLROOT}/libmkl_sequential.a ${MKLROOT}/libmkl_core.a -Wl,--end-group 
 
maxdet_completion.o: maxdet_completion.c
	$(CC) $(CFLAGS) -c $< 

read_and_complete: read_and_complete.c maxdet_completion.o
	$(CC) $(CFLAGS)  $^ $(LDFLAGS)  -o $@ 

read_and_complete2: read_and_complete.o maxdet_completion.o
	$(CC) $(CFLAGS)  $^ -lgsl $(MKLLIB) -lm -o $@ 

clean:
	rm maxdet_completion.o read_and_complete read_and_complete2
