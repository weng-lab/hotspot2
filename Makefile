BINDIR = bin
SRCDIR = src
CXX = g++
CXXFLAGS = -O3

default:
	$(CXX) $(CXXFLAGS) $(SRCDIR)/hotspot2.cpp -o $(BINDIR)/hotspot2
	$(CXX) $(CXXFLAGS) $(SRCDIR)/tallyCountsInSmallWindows.cpp -o $(BINDIR)/tallyCountsInSmallWindows

clean:
	rm $(BINDIR)/*
