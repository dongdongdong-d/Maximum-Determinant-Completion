CFLAGS= -Wall -Wextra -g3
LDFLAGS= -lgsl -lm -lblas
CC=gcc


maxdet_completion.o: maxdet_completion.c
	$(CC) $(CFLAGS)  $(LDFLAGS)  -c $< 
	
read_and_complete: read_and_complete.c maxdet_completion.o
	$(CC) $(CFLAGS)  $^ $(LDFLAGS)  -o $@ 

clean:
	rm maxdet_completion.o read_and_complete
