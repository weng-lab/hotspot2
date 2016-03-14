BINDIR = bin
SRCDIR = src
CXX = g++
CXXFLAGS = -O3 -pedantic -Wall -ansi -static

default:
	mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) $(SRCDIR)/hotspot2.cpp -o $(BINDIR)/hotspot2
	$(CXX) $(CXXFLAGS) $(SRCDIR)/tallyCountsInSmallWindows.cpp -o $(BINDIR)/tallyCountsInSmallWindows

clean:
	rm -f $(BINDIR)/hotspot2
	rm -f $(BINDIR)/tallyCountsInSmallWindows
