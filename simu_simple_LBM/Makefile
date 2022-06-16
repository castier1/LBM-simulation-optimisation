#Commands
CC=gcc
MPICC=mpicc
RM=rm -f
MAKEDEPEND=makedepend

#flags
CFLAGS=-Wall -g -fopenmp
LDFLAGS=-lm
OFLAGS=-O3 -funroll-loops -march=x86-64
VARS=#-DTRACE

#Files
LBM_SOURCES=main.c lbm_phys.c lbm_init.c lbm_struct.c lbm_comm.c lbm_config.c
LBM_HEADERS=$(wildcards:*.h)
LBM_OBJECTS=$(LBM_SOURCES:.c=.o)

TARGET=lbm display

all: $(TARGET)

%.o: %.c
	$(MPICC) $(CFLAGS) $(VARS) $(OFLAGS) -c -o $@ $<

lbm: $(LBM_OBJECTS)
	$(MPICC) $(CFLAGS) $(OFLAGS) -o $@ $^ $(LDFLAGS)

display: display.c
	$(CC) $(CFLAGS) -o $@ display.c

clean:
	$(RM) $(LBM_OBJECTS)
	$(RM) $(TARGET)

depend:
	$(MAKEDEPEND) -Y. $(LBM_SOURCES) display.c

.PHONY: clean all depend