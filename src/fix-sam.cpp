#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <sstream>
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "htslib_wrapper.hpp"
#include "dna_four_bit.hpp"
#include "cigar_utils.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::ostringstream;
using std::stoi;
using std::runtime_error;

size_t
get_tag_value(const sam_rec &aln, const string value) {
  const string to_find = value + ":";
  auto the_nm_tag = find_if(begin(aln.tags), end(aln.tags),
                            [&](const string &t) {
                              return t.compare (0, 3, to_find) == 0;
                            });

  if (the_nm_tag == end(aln.tags))
    return std::numeric_limits<size_t>::max();

  // GS: NM:i:xxxx, get xxxx
  assert(the_nm_tag->size() > 5);
  return stoi(the_nm_tag->substr(5));
}

int
main(int argc, const char **argv) {
  bool verbose = false;
  string outfile = "";
  OptionParser opt_parse(strip_path(argv[0]),
                         "removes malformatted SAM entries",
                        "<sam-input>");
  vector<string> leftover_args;

  opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
  opt_parse.add_opt("output", 'o', "output file", false, outfile);
  opt_parse.parse(argc, argv, leftover_args);
  if (leftover_args.size() != 1) {
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }
  SAMReader in(leftover_args[0]);
  sam_rec aln, mate;
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out << in.get_header();
  while (in >> aln) {
    if (!aln.empty())
      out << aln << '\n'; 
  }
  return EXIT_SUCCESS;
}
