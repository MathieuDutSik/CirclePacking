CFILE = $(wildcard *.C)
OFILE = $(wildcard *.o)
PROG = $(patsubst %.C,%.prog, $(CFILE))
OBJS = $(patsubst %.C,%.o, $(CFILE))

%.o: %.C ; g++ -O2 -Wall $< -c -o $@

%.prog: %.o ; g++ $< -o $@



all: $(PROG) $(OBJS)

clean:
	rm -f $(PROG) $(OBJS)
