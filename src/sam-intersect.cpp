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
  string cigar;
  aln_pos() {
    cigar = "";
    strand = false;
    chrom = "";
    pos = 0;
  }
  aln_pos(const sam_rec &aln) :
    strand(check_flag(aln, samflags::read_rc)), pos(aln.pos), chrom(aln.rname), cigar(aln.cigar) {};
  bool operator==(const aln_pos &rhs) const {
    return (strand == rhs.strand && pos == rhs.pos && chrom == rhs.chrom);
  }
};

int main(int argc, const char **argv) {
  bool verbose = false;
  string same_output_a = "";
  string same_output_b = "";
  string different_output_a = "";
  string different_output_b = "";
  OptionParser opt_parse(strip_path(argv[0]),
                         "finds better matches between sam files sorted by name",
                        "<sam-truth> <sam-input>");
  vector<string> leftover_args;

  opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
  opt_parse.add_opt("same-a", '\0', "print same reads to file", true,
                    same_output_a);
  opt_parse.add_opt("same-b", '\0', "print same reads to file", true,
                    same_output_b);
  opt_parse.add_opt("diff-a", '\0', "print different reads in a to file", true,
                    different_output_a);
  opt_parse.add_opt("diff-b", '\0', "print different reads in b to file", true,
                    different_output_b);

  opt_parse.parse(argc, argv, leftover_args);
  if (leftover_args.size() != 2) {
    cerr << "Need two arguments: <sam-a> <sam-b>\n";
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }

  if (verbose)
    cerr << "[reading SAM file a and storing results]\n";

  SAMReader in_a(leftover_args[0]), in_b(leftover_args[1]);
  sam_rec aln;
  unordered_map<string, sam_rec> a_lookup;

  size_t reads_a = 0;
  size_t reads_b = 0;
  size_t same_reads = 0;
  size_t diff_reads = 0;
  while (in_a >> aln) {
    ++reads_a;
    a_lookup[aln.qname] = aln;
  }

  std::ofstream out_same_a, out_same_b, out_diff_a, out_diff_b;
  out_same_a.open(same_output_a.c_str(), std::ios::binary);
  out_same_b.open(same_output_b.c_str(), std::ios::binary);
  out_diff_a.open(different_output_a.c_str(), std::ios::binary);
  out_diff_b.open(different_output_b.c_str(), std::ios::binary);

  if (verbose)
    cerr << "[reading SAM file b and comparing to a]\n";

  while (in_b >> aln) {
    ++reads_b;
    const auto a_ptr = a_lookup.find(aln.qname);
    if (a_ptr != end(a_lookup)) {
      const aln_pos a = aln_pos(a_ptr->second);
      const aln_pos b = aln_pos(aln);
      if (a == b) {
        ++same_reads;
        out_same_a << a_ptr->second << "\n";
        out_same_b << aln << "\n";
      }
      else {
        ++diff_reads;
        out_diff_a << a_ptr->second << "\n";
        out_diff_b << aln << "\n";
      }
    }
  }

  cout << "reads_a: " << reads_a << "\n";
  cout << "reads_b: " << reads_b << "\n";
  cout << "same_reads: " << same_reads << "\n";
  cout << "diff_reads: " << diff_reads << "\n";
}
