// To compile this code into an executable,
// simply enter the command
//
// $ g++ -O3 tallyCountsInSmallWindows.cpp -o tallyCountsInSmallWindows
//
// or substitute any desired name for the executable for the last argument.
// The argument -O3 (capital "oh") generates optimized code;
// it can be omitted if desired.
// Any C++ compiler can be used in place of g++.
//
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

struct Site {
  string chr;
  long beg; // end - 1
  long end;
  string id;
  int numCuts;
};

string _g_chr;
int _g_start, _g_end, _g_count;
void outputSite(const string chr, const int end, const int cutcount)
{
  //cout << chr << "\t" << end-1 << "\t" << end << "\ti\t" << cutcount << "\n";
  if (chr != _g_chr || end != _g_end + 1 || cutcount != _g_count)
    {
      if (_g_chr != "")
        {
          cout << _g_chr << "\t" << _g_start << "\t" << _g_end << "\ti\t" << _g_count << "\n";
        }
      _g_chr = chr;
      _g_start = end - 1;
      _g_count = cutcount;
    }
  _g_end = end;
}

bool reportNonoverlappingTallies(const char* exeName, const map<string, int>& chromSizes, const int& halfWindowSize, const bool& reportSomethingForEachWin);
bool reportNonoverlappingTallies(const char* exeName, const map<string, int>& chromSizes, const int& halfWindowSize, const bool& reportSomethingForEachWin)
{
  const int BUFSIZE(100);
  char buf[BUFSIZE], *p;
  map<string, int>::const_iterator it = chromSizes.begin();
  Site prevSite, curSite;
  const int windowSize(1 + 2 * halfWindowSize);
  int posL(1), posC(1 + halfWindowSize), posR(windowSize);
  int sum(0);
  int linenum(0), fieldnum;

  prevSite.chr = string("xxxNONExxx");

  while (cin.getline(buf, BUFSIZE))
    {
      linenum++;
      fieldnum = 1;
      if (!(p = strtok(buf, "\t")))
        {
        MissingField:
          cerr << "Error:  Missing required field " << fieldnum
               << " on line " << linenum << " of the input file." << endl
               << endl;
          return false;
        }
      curSite.chr = string(p);
      if (!chromSizes.empty() && curSite.chr != it->first)
        {
          if (chromSizes.end() == chromSizes.find(curSite.chr))
            {
              cerr << "Error:  " << exeName << " encountered unrecognized chromosome \""
                   << curSite.chr << "\" on line " << linenum << '.' << endl;
              return false;
            }
        }
      fieldnum++;
      if (!(p = strtok(NULL, "\t")))
        goto MissingField;
      curSite.beg = atol(p);
      fieldnum++;
      if (!(p = strtok(NULL, "\t")))
        goto MissingField;
      curSite.end = atol(p);
      fieldnum++;
      if (!(p = strtok(NULL, "\t")))
        goto MissingField;
      curSite.id = string(p);
      fieldnum++;
      if (!(p = strtok(NULL, "\t")))
        goto MissingField;
      curSite.numCuts = atoi(p);
      // Ignore any fields beyond the fifth.

      // Safety check
      if (curSite.end != curSite.beg + 1)
        {
          cerr << "Error:  " << exeName << " expects each line of input to represent a single site (col3 = col2 + 1).\n"
               << "The entry in line " << linenum << " was:\n"
               << curSite.chr << '\t' << curSite.beg << '\t' << curSite.end << '\t' << curSite.id
               << '\t' << p << endl;
          return false;
        }

      if (curSite.chr != prevSite.chr || curSite.end > posR)
        {
          // The incoming site lies beyond the bounds of the current window,
          // so write the sum for the current window.
          if (linenum > 1)
            {
              if (chromSizes.empty())
                outputSite(prevSite.chr, posC, sum);
              else
                {
                  if (posC <= it->second)
                    outputSite(prevSite.chr, posC, sum);
                  else // Report the last position in the chromosome as the "center of the window."
                    outputSite(prevSite.chr, it->second, sum);
                }
            }
          else if (!chromSizes.empty())
            it = chromSizes.find(curSite.chr); // guaranteed to succeed if we reach here
          if (reportSomethingForEachWin)
            {
              // Slide the current window rightward, by one full window width at a time (non-overlapping),
              // until the window contains the incoming site.
              // Report sum = 0 until we reach the incoming site.
              if (curSite.chr != prevSite.chr)
                {
                  if (linenum > 1)
                    {
                      if (!chromSizes.empty())
                        {
                          // Write zeroes until the end of the chromosome is reached.
                          posC += windowSize;
                          while (posC <= it->second)
                            {
                              outputSite(prevSite.chr, posC, 0);
                              posC += windowSize;
                            }
                          it = chromSizes.find(curSite.chr); // guaranteed to succeed if we reach here
                        }
                      posL = 1;
                      posC = 1 + halfWindowSize;
                      posR = windowSize;
                    }
                  if (posR < curSite.end)
                    outputSite(curSite.chr, posC, 0);
                }
              // (If there's no data for a whole chromosome, don't bother writing a chromosome's worth of zeroes.)
              while (true)
                {
                  posL += windowSize;
                  posC += windowSize;
                  posR += windowSize;
                  if (posR < curSite.end)
                    outputSite(curSite.chr, posC, 0);
                  else
                    break;
                }
            }
          else
            {
              // Place the left edge of the current window at the location of the incoming site (i.e., leap ahead to it).
              posL = curSite.end;
              posC = posL + halfWindowSize;
              posR = posC + halfWindowSize;
              if (!chromSizes.empty() && curSite.chr != prevSite.chr)
                it = chromSizes.find(curSite.chr); // guaranteed to succeed if we reach here
            }
          sum = 0;
        }

      sum += curSite.numCuts;
      prevSite = curSite;
    }

  // Report the tally corresponding to the end of the file.
  if (chromSizes.empty())
    outputSite(prevSite.chr, posC, sum);
  else
    {
      if (posC <= it->second)
        {
          outputSite(prevSite.chr, posC, sum);
          if (reportSomethingForEachWin)
            {
              // Write zeroes until the end of the chromosome is reached.
              posC += windowSize;
              while (posC <= it->second)
                {
                  outputSite(prevSite.chr, posC, 0);
                  posC += windowSize;
                }
            }
        }
      else
        outputSite(prevSite.chr, it->second, sum);
    }

  return true;
}

struct SiteData {
  int endPos;
  int numCuts;
  bool isMissingData;
};

class OverlordOfOverlapping {
public:
  void initialize(const int& halfWindowSize, const bool& everyBp, const map<string, int>& chromSizes, const char* pExeName);
  bool parseAndProcess(void);

private:
  void processStoredSites(void);
  vector<SiteData> m_storedSites;
  string m_chrom;
  vector<SiteData>::size_type m_idxInsertHere;
  vector<SiteData>::size_type m_halfWindowSize;
  bool m_reportSomethingForEveryBp;
  map<string, int> m_chromSizes;
  int m_sizeOfCurChrom;
  const char* m_exeName;
};

void OverlordOfOverlapping::initialize(const int& halfWindowSize, const bool& everyBp, const map<string, int>& chromSizes, const char* pExeName)
{
  const int VECSIZE(10000000);
  m_halfWindowSize = halfWindowSize;
  m_reportSomethingForEveryBp = everyBp;
  m_chromSizes = chromSizes;
  m_storedSites.resize(VECSIZE);
  m_idxInsertHere = 0;
  m_chrom = string("xxxNONExxx");
  m_sizeOfCurChrom = 0;
  m_exeName = pExeName;
}

void OverlordOfOverlapping::processStoredSites(void)
{
  vector<SiteData>::size_type idxL(0), idxC(m_halfWindowSize), idxR;
  const int windowSize(2 * m_halfWindowSize + 1);
  const bool endOfChrom = m_idxInsertHere < m_storedSites.size() ? true : false;
  int sum(0);

  // First, compute the sum for the leftmost window.
  for (idxR = idxL; idxR < m_idxInsertHere && idxR < idxL + windowSize; idxR++)
    sum += m_storedSites[idxR].numCuts;
  if (m_idxInsertHere != 0)
    {
      if (!m_storedSites[idxC].isMissingData || m_reportSomethingForEveryBp) // guaranteed idxC < idxInsertHere
        outputSite(m_chrom, m_storedSites[idxC].endPos, sum);
    }

  while (idxR < m_idxInsertHere)
    {
      sum -= m_storedSites[idxL++].numCuts;
      sum += m_storedSites[idxR++].numCuts;
      if (!m_storedSites[++idxC].isMissingData || m_reportSomethingForEveryBp)
        outputSite(m_chrom, m_storedSites[idxC].endPos, sum);
    }

  if (endOfChrom)
    {
      if (m_reportSomethingForEveryBp && m_idxInsertHere != 0 && m_sizeOfCurChrom != 0)
        {
          int curPosC = m_storedSites[idxC].endPos, curPosR = curPosC + m_halfWindowSize;
          while (curPosR < m_sizeOfCurChrom)
            {
              if (idxL < m_idxInsertHere)
                sum -= m_storedSites[idxL++].numCuts;
              curPosR++;
              curPosC++;
              outputSite(m_chrom, curPosC, sum);
            }
        }
      m_idxInsertHere = 0;
    }
  else
    {
      // Move the unprocessed data at the end of the vector to the beginning of the vector,
      // for it to be processed during the next call to this function.
      int i = 0;
      idxL++;
      while (idxL < idxR)
        m_storedSites[i++] = m_storedSites[idxL++];
      m_idxInsertHere = i;
    }
}

bool OverlordOfOverlapping::parseAndProcess(void)
{
  const int BUFSIZE(100);
  char buf[BUFSIZE], *p;
  map<string, int>::const_iterator it = m_chromSizes.begin();
  Site curSite;
  int curPos(1);
  int linenum(0), fieldnum;

  while (cin.getline(buf, BUFSIZE))
    {
      linenum++;
      fieldnum = 1;
      if (!(p = strtok(buf, "\t")))
        {
        MissingField:
          cerr << "Error:  " << m_exeName << " encountered a missing required field " << fieldnum
               << " on line " << linenum << " of the input file." << endl
               << endl;
          return false;
        }
      curSite.chr = string(p);
      if (!m_chromSizes.empty() && curSite.chr != it->first)
        {
          if (m_chromSizes.end() == m_chromSizes.find(curSite.chr))
            {
              cerr << "Error:  " << m_exeName << " encountered unrecognized chromosome \""
                   << curSite.chr << "\" on line " << linenum << '.' << endl;
              return false;
            }
        }
      if (curSite.chr != m_chrom) // by design, this also gets triggered when processing the first line of input
        {
          processStoredSites(); // resets m_idxInsertHere to 0
          m_chrom = curSite.chr;
          if (!m_chromSizes.empty())
            {
              it = m_chromSizes.find(curSite.chr); // guaranteed to succeed
              m_sizeOfCurChrom = it->second;
            }
          curPos = 1;
        }
      fieldnum++;
      if (!(p = strtok(NULL, "\t")))
        goto MissingField;
      curSite.beg = atol(p);
      fieldnum++;
      if (!(p = strtok(NULL, "\t")))
        goto MissingField;
      curSite.end = atol(p);
      fieldnum++;
      if (!(p = strtok(NULL, "\t")))
        goto MissingField;
      curSite.id = string(p);
      fieldnum++;
      if (!(p = strtok(NULL, "\t")))
        goto MissingField;
      curSite.numCuts = atoi(p);
      // Ignore any fields beyond the fifth.

      // Safety check
      if (curSite.end != curSite.beg + 1)
        {
          cerr << "Error:  " << m_exeName << " expects each line of input to represent a single site (col3 = col2 + 1).\n"
               << "The entry in line " << linenum << " was:\n"
               << curSite.chr << '\t' << curSite.beg << '\t' << curSite.end << '\t' << curSite.id
               << '\t' << p << endl;
          return false;
        }

      while (curPos < curSite.end)
        {
          // Store a count of 0 for each bp of missing data.
          m_storedSites[m_idxInsertHere].endPos = curPos++;
          m_storedSites[m_idxInsertHere].numCuts = 0;
          m_storedSites[m_idxInsertHere++].isMissingData = true;
          if (m_storedSites.size() == m_idxInsertHere)
            processStoredSites(); // resets m_idxInsertHere to 0
        }
      m_storedSites[m_idxInsertHere].endPos = curPos++;
      m_storedSites[m_idxInsertHere].numCuts = curSite.numCuts;
      m_storedSites[m_idxInsertHere++].isMissingData = false;
      if (m_storedSites.size() == m_idxInsertHere)
        processStoredSites(); // resets m_idxInsertHere to 0
    }

  processStoredSites();

  return true;
}

bool loadChromSizes(ifstream& infile, map<string, int>& chromSizes, const char* exeName);
bool loadChromSizes(ifstream& infile, map<string, int>& chromSizes, const char* exeName)
{
  const int BUFSIZE(100);
  char buf[BUFSIZE], *p;
  string chrom;
  int linenum(0), fieldnum;

  while (infile.getline(buf, BUFSIZE))
    {
      linenum++;
      fieldnum = 1;
      if (!(p = strtok(buf, "\t")))
        {
        MissingData:
          cerr << "Error:  " << exeName << " failed to find required field " << fieldnum
               << " on line " << linenum << " of the file of chromosome names and sizes." << endl;
          return false;
        }
      chrom = string(p);
      fieldnum++;
      if (!(p = strtok(NULL, "\t")))
        goto MissingData;
      chromSizes[chrom] = atoi(p); // Assume 2-column input:  chromosome name, length.
      fieldnum++;
      if (p = strtok(NULL, "\t"))
	{
	  // Assume BED-format input:  chromosome name, start coordinate, end coordinate; ignore any additional columns.
	  int chrLength(atoi(p));
	  if ((0 != chromSizes[chrom] && 1 != chromSizes[chrom]) || chrLength <= chromSizes[chrom])
	    cerr << "Error:  " << exeName << " encountered unexpected or invalid BED input on line "
		 << linenum << " of the file of chromosome names and sizes.\n"
		 << "Expected 2-column input (name and size) or BED format (either name, 1, size or name, 0, size)." << endl;
	  return false;
	  chromSizes[chrom] = chrLength;
	}
    }

  return true;
}

int main(int argc, const char* argv[])
{
  if (argc < 3 || argc > 5)
    {
    Usage:
      cerr << "Usage:  " << argv[0] << " halfWindowSize \"overlapping\"/\"nonoverlapping\" [\"reportEachUnit\"] [fileOfChromSizes]\n"
           << "where\n"
           << "\t* halfWindowSize is the width in bp on each side of each position within which to tally counts (e.g. 100 for 201bp window)\n"
           << "\t* the word \"overlapping\" requests overlapping windows sliding 1bp at a time, \"nonoverlapping\" requests non-overlapping windows\n"
           << "\t* (optional) the string \"reportEachUnit\" reports a tally for every non-overlapping window or for every central bp\n"
           << "\t* (optional) a file of chromosome sizes can be provided for error-checking and guidance at the ends of chromosomes"
           << endl
           << endl;
      return -1;
    }

  ios_base::sync_with_stdio(false); // calling this static method in this way turns off checks, speeds up I/O

  int halfWinSize(atoi(argv[1]));
  bool overlappingWindowsSliding1bp;
  string thisArgString;
  const char* p = argv[2];
  while (*p)
    thisArgString += tolower(*p++);
  if (string("overlapping") == thisArgString)
    overlappingWindowsSliding1bp = true;
  else
    {
      if (string("nonoverlapping") == thisArgString)
        overlappingWindowsSliding1bp = false;
      else
        {
          cerr << "Unrecognized 2nd argument \"" << argv[2] << "\" to " << argv[0]
               << "; must be \"overlapping\" or (with no hyphen) \"nonoverlapping\"." << endl;
          goto Usage;
        }
    }
  ifstream ifsChromSizes;
  bool reportEachUnit(false);
  if (argc > 3)
    {
      thisArgString = string("");
      p = argv[3];
      while (*p)
        thisArgString += tolower(*p++);
      if (string("reporteachunit") == thisArgString)
        reportEachUnit = true;
      else
        {
          ifsChromSizes.open(argv[3]);
          if (!ifsChromSizes)
            {
              cerr << "Unrecognized 3rd argument \"" << argv[3] << "\" to " << argv[0]
                   << "; must be \"reportEachUnit\" or the path to a file of chromosome names and sizes." << endl;
              goto Usage;
            }
        }
      if (5 == argc)
        {
          if (reportEachUnit)
            {
              ifsChromSizes.open(argv[4]);
              if (!ifsChromSizes)
                {
                  cerr << argv[0] << " was unable to open file \"" << argv[4] << "\" for reading." << endl;
                  goto Usage;
                }
            }
          else
            {
              thisArgString = string("");
              p = argv[4];
              while (*p)
                thisArgString += tolower(*p++);
              if (string("reporteachunit") == thisArgString)
                reportEachUnit = true;
              else
                {
                  cerr << "Unrecognized 4th argument \"" << argv[4] << "\" to " << argv[0] << '.' << endl;
                  goto Usage;
                }
            }
        }
    }

  map<string, int> chrSizes;
  if (ifsChromSizes)
    {
      if (!loadChromSizes(ifsChromSizes, chrSizes, argv[0]))
        return -1;
      ifsChromSizes.close();
      ifsChromSizes.clear();
    }

  if (overlappingWindowsSliding1bp)
    {
      OverlordOfOverlapping o;
      o.initialize(halfWinSize, reportEachUnit, chrSizes, argv[0]);
      if (!o.parseAndProcess())
        return -1;
    }
  else
    {
      if (!reportNonoverlappingTallies(argv[0], chrSizes, halfWinSize, reportEachUnit)) // chrSizes can be empty
        return -1;
    }

  // Force flushing
  outputSite("", 0, 0);
  return 0;
}
