SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio128/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio128/concert
CPOPTDIR      = /opt/ibm/ILOG/CPLEX_Studio128/cpoptimizer

CCC = g++ -O0
CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -Wno-ignored-attributes

CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPOPTLIBDIR   = $(CPOPTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) -L$(CPOPTLIBDIR)
CCLNFLAGS = -lcp -lconcert -lilocplex -lcplex -lm -lpthread -ldl -Wall -std=c++17

CPLEXINCDIR   = $(CPLEXDIR)/include
CONCERTINCDIR = $(CONCERTDIR)/include
CPOPTINCDIR   = $(CPOPTDIR)/include

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(CPOPTINCDIR)

all:
	$(CCC) main.cc problem.cc $(CCFLAGS) $(CCLNDIRS) -o ipb $(CCLNFLAGS)

clean:
	rm *.o ipb
