#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <sstream>
#include <unordered_map>
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "htslib_wrapper.hpp"
#include "dna_four_bit.hpp"

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
using std::unordered_map;

struct aln_pos {
  bool strand;
  uint32_t pos;
  string chrom;
  aln_pos() {
    strand = false;
    chrom = "";
    pos = 0;
  }
  aln_pos(const sam_rec &aln) :
    strand(check_flag(aln, samflags::read_rc)), pos(aln.pos), chrom(aln.rname) {};
};

bool
is_valid_aln(const aln_pos &truth, const aln_pos &input) {
  static const long long int max_diffs_for_correct = 50;
  return ((truth.strand == input.strand) && 
           (truth.chrom == input.chrom) &&
           (abs(static_cast<long long int>(truth.pos) - static_cast<long long int>(input.pos)) <=
            max_diffs_for_correct));
}

int main(int argc, const char **argv) {
  bool verbose = false;
  string incorrect_output = "";
  OptionParser opt_parse(strip_path(argv[0]),
                         "finds better matches between sam files sorted by name",
                        "<sam-truth> <sam-input>");
  vector<string> leftover_args;

  opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
  opt_parse.add_opt("output", 'o', "print incorrect reads to file", false,
                    incorrect_output);
  opt_parse.parse(argc, argv, leftover_args);
  if (leftover_args.size() != 2) {
    cerr << "Need two arguments: <sam-truth> <sam-input>\n";
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }

  if (verbose)
    cerr << "[reading truth SAM file and storing results]\n";
  SAMReader in_truth(leftover_args[0]), in_input(leftover_args[1]);
  sam_rec aln;
  unordered_map<string, aln_pos> the_truth;

  size_t tot_reads = 0;
  size_t correct_reads = 0;
  size_t incorrect_reads = 0;
  while (in_truth >> aln) {
    ++tot_reads;
    the_truth[aln.qname] = aln_pos(aln);
  }

  const bool print_incorrect = !incorrect_output.empty();
  std::ofstream out_incorrect;
  if (print_incorrect)
    out_incorrect.open(incorrect_output.c_str(), std::ios::binary);

  while (in_input >> aln) {
    const auto the_truth_ptr = the_truth.find(aln.qname);
    if (the_truth_ptr == end(the_truth))
      throw runtime_error("mapped read does not exist in truth: " + aln.qname);

    const aln_pos the_input = aln_pos(aln);
    const bool valid = is_valid_aln(the_truth_ptr->second, the_input);
    correct_reads += valid;
    incorrect_reads += !valid;
    if (print_incorrect & !valid)
      out_incorrect << aln << "\n";
  }

  const size_t reported_reads = correct_reads + incorrect_reads;
  cout << "tot_reads: " << tot_reads << "\n";
  cout << "correct_reads: " << correct_reads << "\n";
  cout << "incorrect_reads: " << incorrect_reads << "\n";
  cout << "sensitivity: " << correct_reads /static_cast<double>(tot_reads) << "\n";
  cout << "specificity: " << correct_reads /static_cast<double>(reported_reads) << "\n";
}
