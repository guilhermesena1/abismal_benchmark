#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <sstream>
#include <unordered_map>
#include <algorithm>
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
using std::unordered_map;
using std::min;
using std::max;

static size_t
erase_frag_from_name(sam_rec &aln) {
  if (aln.qname.find("FRAG:") == 0) {
    aln.qname.erase(0, 5);
    return 2;
  }
  return 1;
}


void
get_chrom_and_pos(const sam_rec &aln, 
                  string &chrom, uint32_t &the_start, uint32_t &the_end) {
  const string the_qname = aln.qname;

  size_t current, previous = 0;

  // get the read number
  current = the_qname.find_first_of("_");
  previous = current + 1;

  // get the chrom
  current = the_qname.find_first_of(":", previous);
  chrom = the_qname.substr(previous, current - previous);
  previous = current + 1;

  // get the start
  current = the_qname.find_first_of("-", previous);
  const string start_str = the_qname.substr(previous, current - previous);

  const bool is_start_numeric = 
    (std::count_if(begin(start_str), end(start_str),
        [](const char c) { return !isdigit(c); }) == 0);

  if (is_start_numeric)
    the_start = stoi(start_str);
  previous = current + 1;

  // get the end
  current = the_qname.find_first_of("_", previous);
  const string end_str = the_qname.substr(previous, current - previous);

  const bool is_end_numeric = 
    (std::count_if(begin(end_str), end(end_str),
        [](const char c) { return !isdigit(c); }) == 0);

  if (is_end_numeric)
    the_end = stoi(end_str) + 1;
}

size_t
get_edit_distance(const sam_rec &aln) {
  auto the_nm_tag = find_if(begin(aln.tags), end(aln.tags),
                            [](const string &t) {
                              return t.compare (0, 3, "NM:") == 0;
                            });

  // GS: NM:i:xxxx, get xxxx
  assert(the_nm_tag->size() > 5);
  return stoi(the_nm_tag->substr(5));
}

bool
check_valid(const sam_rec &aln) {
  /*
  if (get_edit_distance(aln) == 0)
    return true;*/
  try {
    string chrom;
    uint32_t the_start, the_end;
    const uint32_t readlen = aln.seq.size();
    get_chrom_and_pos(aln, chrom, the_start, the_end);

    the_start = min(the_start, the_end - readlen);
    the_end = max(the_end, the_start + readlen);

    const uint32_t read_st = aln.pos - 1;
    const uint32_t read_end = read_st + cigar_rseq_ops(aln.cigar);

    const bool same_chrom = (aln.rname == chrom);
    const bool same_pos =
      ((read_st >= the_start && read_st <= the_end) ||
      (read_end >= the_start && read_end <= the_end));
    return same_chrom && same_pos;
  } catch (const std::invalid_argument &e) {
    cerr << "problem with line : " << aln << "\n";
    return true;
  }
  return false;
}

int main(int argc, const char **argv) {
  bool verbose = false;
  size_t total_reads = 0;
  string incorrect_output = "";
  OptionParser opt_parse(strip_path(argv[0]),
                         "finds better matches between sam files sorted by name",
                        "<sam-truth> <sam-input>");
  vector<string> leftover_args;

  opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
  opt_parse.add_opt("num-reads", 'n', "number of reads simulated", true, total_reads);
  opt_parse.add_opt("output", 'o', "print incorrect reads to file", false,
                    incorrect_output);
  opt_parse.parse(argc, argv, leftover_args);
  if (leftover_args.size() != 1) {
    cerr << "Need two arguments: <sam-truth> <sam-input>\n";
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }
  if (verbose)
    cerr << "[reading truth SAM file and storing results]\n";
  SAMReader in(leftover_args[0]);
  sam_rec aln;

  size_t correct_reads = 0;
  size_t incorrect_reads = 0;

  const bool print_incorrect = !incorrect_output.empty();
  std::ofstream out_incorrect;
  if (print_incorrect)
    out_incorrect.open(incorrect_output.c_str(), std::ios::binary);

  while (in >> aln) {
    if (aln.rname != "*") {
      const bool valid = check_valid(aln);
      const size_t mult = erase_frag_from_name(aln);
      correct_reads += mult*valid;
      incorrect_reads += mult*(!valid);
      if (print_incorrect && !valid)
        out_incorrect << aln << "\n";
    }
  }

  const size_t missed_reads = total_reads - correct_reads;
  cout << "correct_reads: " << correct_reads << "\n";
  cout << "incorrect_reads: " << incorrect_reads << "\n";
  cout << "sensitivity: " << correct_reads/static_cast<double>(correct_reads + missed_reads) << endl;
  cout << "specificity: " << correct_reads/static_cast<double>(correct_reads + incorrect_reads) << endl;
  cout << "F1: " << correct_reads/static_cast<double>(correct_reads + (incorrect_reads + missed_reads)/2.0) << endl;
}
