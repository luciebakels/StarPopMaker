CXX = mpic++
RM = rm -f
CPPFLAGS = -g -std=c++11
MPICHLIB = -lmpich

SRCS = run.cpp radiation.cpp population.cpp allvars.cpp read.cpp write.cpp star.cpp gas.cpp general.cpp profiles.cpp
OBJS = $(subst .cpp,.o,$(SRCS))

all: mpipop

mpipop: $(OBJS)
	$(CXX) $(CPPFLAGS) $(MPICHLIB) -o mpipop $(OBJS)

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) $(MPICHLIB) -MM $^>>./.depend;

clean:
	$(RM) mpipop $(OBJS)

distclean: clean
	$(RM) *~ .depend

include .depend
