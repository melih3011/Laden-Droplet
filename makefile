CC = gcc-mp-13
CFLAGS = -O2 -fopenmp
CLIBS = -lm 

OBJS = sp_one.o white_bear.o

OBJS2 = sp_one.o white_bear_q3.o

OBJS3 = sp_one.o rosenfeld.o

OBJS4 = sp_one.o lutsko.o

OBJS5 = sp_one.o lutsko_wb.o

OBJS6 = sp_one.o lutsko_wb_mark_II.o

all: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(CLIBS) -o spone_wb

wbq3: $(OBJS2)
	$(CC) $(CFLAGS) $(OBJS2) $(CLIBS) -o spone_wbq3

rf: $(OBJS3)
	$(CC) $(CFLAGS) $(OBJS3) $(CLIBS) -o spone_rf

lk: $(OBJS4)
	$(CC) $(CFLAGS) $(OBJS4) $(CLIBS) -o spone_lk

lkwb: $(OBJS5)
	$(CC) $(CFLAGS) $(OBJS5) $(CLIBS) -o spone_lkwb

lkwb2: $(OBJS6)
	$(CC) $(CFLAGS) $(OBJS6) $(CLIBS) -o spone_lkwb2

clean : 
	rm -f *.o core 
