BINDIR = bin
SRCDIR = src
CXX = g++
CXXFLAGS = -O2 -march=native -mtune=native -funroll-loops -DNDEBUG -static

TARGETS = hotspot2_part1 hotspot2_part2
EXE = $(addprefix $(BINDIR)/,$(TARGETS))

default: $(EXE)

$(BINDIR)/% : $(SRCDIR)/%.cpp
	mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(EXE)
