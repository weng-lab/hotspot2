// To compile this code into an executable,
// simply enter the command
//
// $ g++ -O3 hotspot2_part2.cpp -o hotspot2_part2
//
// or substitute any desired name for the executable for the last argument.
// The argument -O3 (capital "oh") generates optimized code;
// it can be omitted if desired.
// Any C++ compiler can be used in place of g++.
//
#include "hotspot2_version.h" // for versioning
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime> // for seeding the random number generator
#include <deque>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <limits> // for epsilon()
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility> // for pair
#include <vector>

const long double CHANGE_OF_SCALE(10000.);

struct SiteRangeData {
  int negLog10P_scaled;
  int begPos;
  int width;
  int chromID;
  long double FDR;
};

bool NegLog10Pscaled_GT(const SiteRangeData& a, const SiteRangeData& b)
{
  return a.negLog10P_scaled > b.negLog10P_scaled;
}

bool OutputOrder_LT(const SiteRangeData& a, const SiteRangeData& b)
{
  if (a.chromID == b.chromID)
    return a.begPos < b.begPos;
  return a.chromID < b.chromID;
}

bool buildFDRmapping(std::ifstream& ifs, std::vector<std::pair<int, long double> >& p_to_q)
{
  const int BUFSIZE(100);
  char buf[BUFSIZE], *p;
  long numPvalues(0);
  std::map<int, int> tempMap;
  std::map<int, int>::iterator it;
  std::pair<int, int> negLog10PscaledAndNumOccs;
  int linenum(0), fieldnum;

  while (ifs.getline(buf,BUFSIZE))
    {
      linenum++;
      fieldnum = 1;
      if (!(p = strtok(buf, "\t")) || !*p)
	{
	MissingField:
	  std::cerr << "Error:  Failed to find field " << fieldnum
	       << " on line " << linenum << " of the file of P-values."
	       << std::endl << std::endl;
	  return false;
	}
      negLog10PscaledAndNumOccs.first = atoi(p);
      fieldnum++;
      if (!(p = strtok(NULL, "\t")))
	goto MissingField;
      negLog10PscaledAndNumOccs.second = atoi(p);
      
      it = tempMap.lower_bound(negLog10PscaledAndNumOccs.first); // >=
      if (tempMap.end() == it)
	{
	  // no lower bound found in the map; every map element is smaller, or the map is empty
	  if (!tempMap.empty())
	    tempMap.insert(--it, negLog10PscaledAndNumOccs);
	  else
	    tempMap.insert(negLog10PscaledAndNumOccs);
	}
      else
	{
	  if (it->first == negLog10PscaledAndNumOccs.first)
	    it->second += negLog10PscaledAndNumOccs.second;
	  else
	    {
	      if (it != tempMap.begin())
		it--;
	      tempMap.insert(it, negLog10PscaledAndNumOccs);
	    }
	}
      numPvalues += negLog10PscaledAndNumOccs.second;
    }

  if (0 == numPvalues)
    {
      std::cerr << "Error:  Received an empty file of P-values." << std::endl << std::endl;
      return false;
    }

  p_to_q.resize(tempMap.size()); // assume one entry per map entry; final size will almost surely be a bit smaller
  
  it = tempMap.end();
  it--;
  int i = 0;
  long double prevFDR(-1.), FDR, N(static_cast<long double>(numPvalues)), numThisExtremeOrMoreExtreme(0);
  while (it != tempMap.begin())
    {
      p_to_q[i].first = it->first;
      numThisExtremeOrMoreExtreme += static_cast<long double>(it->second);
      FDR = static_cast<long double>(std::pow(10., -it->first/CHANGE_OF_SCALE)) * N / numThisExtremeOrMoreExtreme;
      if (FDR < prevFDR)
	FDR = prevFDR;
      if (FDR > 0.999)
	{
	  p_to_q[i++].second = 1.;
	  p_to_q.resize(i);
	  break;
	}
      p_to_q[i].second = FDR;
      prevFDR = FDR;
      i++;
      it--;
    }

  // It's extremely unlikely, but possible that only the smallest observation has B-H FDR = 1.
  if (i > 0 && p_to_q[i-1].second < 1.)
    {
      p_to_q[i].first = tempMap.begin()->first;
      p_to_q[i].second = 1.;
    }

  return true;
}

bool buildIntToChromNameMap(std::ifstream& infile, std::map<int, std::string*>& mapOut)
{
  const int BUFSIZE(100);
  char buf[BUFSIZE], *p;
  int ID;
  std::map<int, std::string*>::const_iterator it;
  int linenum(0), fieldnum;

  while (infile.getline(buf, BUFSIZE))
    {
      linenum++;
      fieldnum = 1;
      if (!(p = strtok(buf, "\t")) || !*p)
	{
	MissingField:
	  std::cerr << "Error:  Failed to find field " << fieldnum
	       << " on line " << linenum
	       << " of the number-to-chromosomeName mapping file."
	       << std::endl << std::endl;
	  return false;
	}
      ID = atoi(p);
      it = mapOut.find(ID);
      if (it != mapOut.end())
	{
	  std::cerr << "Error:  Multiple entries found in the number-to-chromosomeName mapping file\n"
	       << "with value " << ID << " in the first column; each mapping must be unique."
	       << std::endl << std::endl;
	  return false;
	}
      fieldnum++;
      if (!(p = strtok(NULL, "\t")) || !*p)
	goto MissingField;
      mapOut[ID] = new std::string(p);
    }
  
  return true;
}

bool parseAndProcessInput(const std::map<int, std::string*>& intToChromNameMap,
			  const std::vector<std::pair<int, long double> >& PvalToFDRmapping,
			  const int& N, const long double& FDRthreshold, const bool& writePvals)
{
  const int BUFSIZE(1000);
  char buf[BUFSIZE], *p;
  std::vector<SiteRangeData> vec;
  std::vector<SiteRangeData>::iterator it;
  int linenum(0), fieldnum;

  vec.resize(N);

  while (std::cin.getline(buf, BUFSIZE))
    {
      linenum++;
      if (linenum > N)
	{
	  std::cerr << "Error:  Expected to find exactly " << N
	       << " lines of data in the file of location and P-value data, but at least one additional line was found."
	       << std::endl << std::endl;
	  return false;
	}
      if (1 == linenum)
	it = vec.begin();
      else
	it++;
      fieldnum = 1;
      if (!(p = strtok(buf, "\t")) || !*p)
	{
	MissingField:
	  std::cerr << "Error:  Failed to find required field " << fieldnum
	       << " on line " << linenum << " of the file of location and P-value data."
	       << std::endl << std::endl;
	  return false;
	}
      it->chromID = atoi(p);
      fieldnum++;
      if (!(p = strtok(NULL, "\t")) || !*p)
	goto MissingField;
      it->begPos = atoi(p);
      fieldnum++;
      if (!(p = strtok(NULL, "\t")) || !*p)
	goto MissingField;
      it->width = atoi(p);
      fieldnum++;
      if (!(p = strtok(NULL, "\t")) || !*p)
	goto MissingField;
      it->negLog10P_scaled = atoi(p);
    }

  if (linenum < N)
    {
      std::cerr << "Error:  Expected to find exactly " << N
	   << " lines of data in the file of location and P-value data,"
	   << " but only " << linenum << " lines were found."
	   << std::endl << std::endl;
      return false;
    }

  sort(vec.begin(), vec.end(), NegLog10Pscaled_GT);

  if (vec[0].negLog10P_scaled != PvalToFDRmapping[0].first)
    {
      std::cerr << "Coding error:  Line " << __LINE__
	   << ", expected scaled -log10(P) values in the initial elements\n"
	   << "of the genomic data and P-to-FDR vectors to be equal, but they are not\n"
	   << '(' << vec[0].negLog10P_scaled << " for the former and "
	   << PvalToFDRmapping[0].first << " for the latter.)"
	   << std::endl << std::endl;
      return false;
    }

  it = vec.begin();
  for (std::vector<std::pair<int, long double> >::const_iterator itPtoQ = PvalToFDRmapping.begin();
       itPtoQ != PvalToFDRmapping.end(); itPtoQ++)
    {
      while (it != vec.end() && it->negLog10P_scaled == itPtoQ->first)
	{
	  it->FDR = itPtoQ->second;
	  it++;
	}
    }
  while (it != vec.end())
    {
      it->FDR = 1.;
      it++;
    }
  
  sort(vec.begin(), vec.end(), OutputOrder_LT);  
  int prevChromID(-1);
  std::string *prevChromString(NULL);
  
  for (it = vec.begin(); it != vec.end(); it++)
    {
      // Even though the map is tiny, perform as few map lookups as possible.
      if (it->chromID != prevChromID)
	{
	  prevChromString = intToChromNameMap.at(it->chromID);
	  prevChromID = it->chromID;
	}
      if (it->FDR <= FDRthreshold)
	{
	  std::cout << *prevChromString << '\t'
	       << it->begPos << '\t'
	       << it->begPos + it->width << "\ti\t"
	       << it->FDR;
	  if (writePvals)
	    std::cout << '\t' << std::pow(10., -it->negLog10P_scaled/CHANGE_OF_SCALE);
	  std::cout << '\n';
	}
    }
      
  return true;
}


int main(int argc, char* argv[])
{
  // Option defaults
  long double fdr_threshold = 1.00;
  int write_pvals = 0;
  int print_help = 0;
  int print_version = 0;
  std::string infilename = "";
  std::string infilePvals = "";
  std::string infileChromNames = "";
  std::string outfilename = "";
  int numEntries = 0;

  // Long-opt definitions
  static struct option long_options[] = {
    { "fdr_threshold", required_argument, 0, 'f' },
    { "num_entries", required_argument, 0, 'n' },
    { "write_pvals", no_argument, &write_pvals, 1 },
    { "input", required_argument, 0, 'i' },
    { "infileChromNames", required_argument, 0, 'c' },
    { "infilePvalData", required_argument, 0, 'p' },
    { "output", required_argument, 0, 'o' },
    { "help", no_argument, &print_help, 1 },
    { "version", no_argument, &print_version, 1 },
    { 0, 0, 0, 0 }
  };

  // Parse options
  char c;
  std::stringstream ss; // Used for parsing long doubles (allows scientific notation)
  while ((c = getopt_long(argc, argv, "f:n:c:p:i:o:hvV", long_options, NULL)) != -1)
    {
      switch (c)
        {
        case 'f':
          ss << optarg;
          ss >> fdr_threshold;
          break;
        case 'n':
          numEntries = atoi(optarg);
          break;
        case 'i':
          infilename = optarg;
          break;
        case 'c':
          infileChromNames = optarg;
          break;
        case 'p':
          infilePvals = optarg;
          break;
	case 'o':
          outfilename = optarg;
          break;
        case 'h':
          print_help = 1;
          break;
        case 'v':
        case 'V':
          print_version = 1;
          break;
        // no short option needed for --write_pvals
        case 0:
          // long option received, do nothing
          break;
        default:
          print_help = 1;
        }
    }

  if (!print_help && !print_version && infileChromNames.empty())
    {
      std::cerr << "Error:  Required file containing mapping from integers to chromosome names was not supplied."
	   << std::endl << std::endl;
      print_help = 1;
    }
  
  if (!print_help && !print_version && infilePvals.empty())
    {
      std::cerr << "Error:  Required file containing scaled -log10(P) values and their occurrence counts was not supplied."
	   << std::endl << std::endl;
      print_help = 1;
    }
 
  // Print usage and exit if necessary
  if (print_help)
    {
      std::cerr << "Usage:  " << argv[0] << " [options] < in.PvalueData.txt > out.FDR.bed\n"
           << "\n"
           << "Options: \n"
           << "  -n, --num_entries=INT          The number of lines of input (required for optimization)\n"
           << "  -p, --inputPvalData=FILE       A file of scaled -log10(P) values and their occurrence counts\n"
	   << "  -c, --inputChromNames=FILE     A file containing the mapping from integers to chromosome names\n"
           << "  --write_pvals                  Output P-values in column 6 (P-values are not output by default)\n"
           << "  -f, --fdr_threshold=THRESHOLD  Do not output sites with FDR > THRESHOLD (1.00)\n"
           << "  -i, --input=FILE               A file to read input from (STDIN)\n"
           << "  -o, --output=FILE              A file to write output to (STDOUT)\n"
           << "  -v, --version                  Print the version information and exit\n"
           << "  -h, --help                     Display this helpful help\n"
           << "\n"
           << " output (sent to stdout) will be a .bed5 file with FDR in field 5\n"
           << "\tor, if --write_pvals is specified, a .bed6 file with P-values appended in field 6\n"
           << " input (via \"-i FILE\" or piped in from stdin) consists of many rows of values in 4 columns.\n"
           << std::endl
           << std::endl;
      return -1;
    }

  if (print_version)
    {
      std::cout << argv[0] << " version " << hotspot2_VERSION_MAJOR
           << '.' << hotspot2_VERSION_MINOR << std::endl;
      return 0;
    }

  std::ios_base::sync_with_stdio(false); // calling this static method in this way turns off checks, speeds up I/O

  if (!infilename.empty() && infilename != "-")
    {
      if (freopen(infilename.c_str(), "r", stdin) == NULL)
        {
          std::cerr << "Error: Couldn't open input file " << infilename << std::endl;
          return 1;
        }
    }
  if (!outfilename.empty() && outfilename != "-")
    {
      if (freopen(outfilename.c_str(), "w", stdout) == NULL)
        {
          std::cerr << "Error: Couldn't open output file " << outfilename << " for writing" << std::endl;
          return 1;
        }
    }

  std::ifstream ifsChromNames(infileChromNames.c_str()), ifsPvals(infilePvals.c_str());
  if (!ifsChromNames)
    {
      std::cerr << "Error:  Unable to open file \"" << infileChromNames << "\" for read." << std::endl;
      return 1;
    }
  if (!ifsPvals)
    {
      std::cerr << "Error:  Unable to open file \"" << infilePvals << "\" for read." << std::endl;
      return 1;
    }

  std::map<int, std::string*> IntToChromNameMap;
  if (!buildIntToChromNameMap(ifsChromNames, IntToChromNameMap))
    return -1;

  std::vector<std::pair<int, long double> > PvalToFDRmapping;
  if (!buildFDRmapping(ifsPvals, PvalToFDRmapping))
    return -1;

  if (!parseAndProcessInput(IntToChromNameMap, PvalToFDRmapping, numEntries, fdr_threshold, write_pvals ? true : false))
    return -1;

  return 0;
}
