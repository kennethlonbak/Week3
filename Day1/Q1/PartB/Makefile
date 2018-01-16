TARGET	= PartB
OBJS	= $(TARGET).o
TARGET_ext = .prog

FC  	= f95
#FFLAGS 	= -free -ftrap=%none
FC 		= gfortran
FFLAGS 	= -free -O3
all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(TARGET)$(TARGET_ext)


clean:
	@/bin/rm -fv $(TARGET)$(TARGET_ext)
	@/bin/rm -fv $(OBJS)


%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# dependencies
$(TARGET).o : $(TARGET).f90


