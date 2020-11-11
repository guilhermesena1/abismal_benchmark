#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <fstream>
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::runtime_error;
using std::ostringstream;
using std::istream;
using std::ostream;
using std::ifstream;

struct bed_entry {
  string chrom;
  size_t start, end;
  string name;
  size_t score;
  char strand;
  size_t a,b,c;
  double cov;
};

template <class T>
T& operator >>(T &in, bed_entry &rhs) {
  in >> rhs.chrom >> rhs.start >> rhs.end >> rhs.name >> rhs.score
     >> rhs.strand >> rhs.a >> rhs.b >> rhs.c >> rhs.cov;
  return in;
}

template <class T>
T& operator <<(T &out, const bed_entry &rhs) {
  out << rhs.chrom << '\t' << rhs.start << '\t'
      << rhs.end << '\t' << rhs.name << '\t' << rhs.score
      << '\t' << rhs.strand << '\t' << rhs.a << '\t' << rhs.b
      << '\t' << rhs.c << '\t' << rhs.cov;
  return out;
}

static bool
equivalent(const bed_entry &a, const bed_entry &b) {
  return a.chrom == b.chrom && a.start == b.start && a.end == b.end &&
         a.strand == b.strand;
}

int
main(int argc, const char **argv) {

  bool verbose = false;
  OptionParser opt_parse(strip_path(argv[0]),
                        "find entries with more coverage in <1> than in <2>",
                        "<good-bed> <bad-bed>");
  opt_parse.add_opt("verbose", 'v', "print more run info",false, verbose);
  vector<string> leftover_args;
  opt_parse.parse(argc, argv, leftover_args);
  bed_entry bed_good, bed_bad;

  if (leftover_args.size() != 2) {
    cerr << "Need two arguments: <good-bed> and <bad-bed>";
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }

  ifstream in_good(leftover_args[0]);
  ifstream in_bad(leftover_args[1]);

  if (!in_good.good())
    throw runtime_error("problem opening file: " + leftover_args[0]);

  if (!in_bad.good())
    throw runtime_error("problem opening file: " + leftover_args[1]);

  while((in_good >> bed_good) && (in_bad >> bed_bad)) {
    if (!equivalent(bed_good, bed_bad)) {
      ostringstream fail;
      fail << "bed entries in same line not equivalent\n";
      fail << bed_good << '\n';
      fail << bed_bad << '\n';
      throw runtime_error(fail.str());
    }
    if (bed_good.cov > bed_bad.cov) {
      cout << bed_good << '\n';
    }
  }

  return EXIT_SUCCESS;
}

