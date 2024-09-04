CCSRCFILES	= CirclePackSphere.c

PROGRAM		= CirclePackSphere

OBJECTS		= $(PROGRAM:%=%.o) \
		  $(CCSRCFILES:.C=.o) 

WARNING         = -Wall
#WARNING         = 

LINKOBJECTS	= $(WARNING) $(CCSRCFILES:.c=.o)

#
# flags for compiler and linker
# =============================
#

#CC		= gcc -g
CC		= gcc -g

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

OPTLEVEL	= -O3

CCFLAGS 	= $(WARNING) $(OPTLEVEL)

LIBS		= -lm

#
# dependencies and rules
# ======================
#

%.o: %.c
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CCFLAGS) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CC) $(LINKOBJECTS) $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f *.c~
	rm -f *~
	rm -f $(PROGRAM)
#
# more dependencies
# =================
#

CirclePackSphere.o:
