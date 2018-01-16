TARGET	= PartA
OBJS	= PartA.o
TARGET_ext = .prog

FC  	= f95
FFLAGS 	= -free -ftrap=%none
FC 		= gfortran
FFLAGS 	= -free
all: $(TARGET).prog

$(TARGET).prog: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(TARGET)$(TARGET_ext)


clean:
	@/bin/rm -fv $(TARGET)$(TARGET_ext)
	@/bin/rm -fv $(OBJS)


%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# dependencies
Time_loop.o : PartA.f90

