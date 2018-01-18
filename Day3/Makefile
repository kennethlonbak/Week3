TARGET	= GS_SOLVER
OBJS	= m_$(TARGET).o $(TARGET).o
TARGET_EXT = .prog

FC  	= mpif90
FFLAGS = -free -g -O0

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(TARGET)$(TARGET_EXT)

clean:
	@/bin/rm -fv $(TARGET)$(TARGET_EXT)
	@/bin/rm -fv $(OBJS)
	@/bin/rm -fv m_$(TARGET).mod


%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# dependencies
GS_SOLVER.o : GS_SOLVER.f90
m_GS_SOLVER.o : m_GS_SOLVER.f90

