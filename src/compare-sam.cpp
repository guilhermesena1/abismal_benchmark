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

void
print_debug(const sam_rec &input, const sam_rec &truth) {
  cerr << "==========================================\n";
  cerr << "input: " << input << '\n';
  cerr << "truth: " << truth << '\n';
  cerr << "==========================================\n";
}

template <typename T> T&
write_read_for_debug(const sam_rec &aln, T &out) {
  out << "@" << aln.qname << "_input_" << aln.rname << '_' << aln.pos << '\n'
  // bsmap puts read in the direction of the reference
      << (check_flag(aln, samflags::read_rc) ? revcomp(aln.seq) : aln.seq)
      << "\n+\n" << std::string(aln.seq.size(), 'B') << '\n';
  return out;
}

// goal is to quantify:
// (2) superoptimal_uniqs: truth worse than input
// (1) optimal_uniqs: equal to the truth and map to same positions
// (2) unknown uniqs: equal to truth and map to different position
// (2) suboptimal_uniqs: truth better, uniq in truth
// (3) false_uniqs: truth better, ambig in truth
//
struct mapstats {
  size_t truth_entries;
  size_t input_entries;
  size_t both;
  // better in input than in truth
  size_t superoptimal_uniq;

  // uniq in both, equal score, same pos
  size_t optimal_uniq;

  // uniq in both, equal score, differnet pos
  size_t unknown_uniq;

  // uniq in both, better in truth than input
  size_t suboptimal_uniq;

  // ambig in truth, better in truth than input
  size_t false_uniq;

  // uniq in truth, not reported by mapper
  size_t missed_uniq;

  // uniq in mapper, not reported by truth
  size_t dangling_uniq;

  // paired in truth, unpaired in input
  size_t unpaired;

  // paired in input, unpaired in truth
  size_t paired_truth_ambig;
  size_t paired_unknown;

  mapstats();
  string tostring();
  void update(const sam_rec &input, const sam_rec &truth);
};

mapstats::mapstats() {
  superoptimal_uniq = optimal_uniq = unknown_uniq = suboptimal_uniq
                    = false_uniq = missed_uniq = dangling_uniq
                    = paired_truth_ambig = paired_unknown = unpaired
                    = truth_entries = input_entries = both = 0;
}

string
mapstats::tostring() {
  static const string tab = "    ";
  const size_t pos =
    superoptimal_uniq + optimal_uniq + unknown_uniq + dangling_uniq;
  const size_t neg = 
    false_uniq + missed_uniq + suboptimal_uniq;

  const double accuracy = pos / static_cast<double>(pos + neg);
  ostringstream out;
  out << "total:\n"
      << tab << "both: " << both << '\n'
      << tab << "input_only: " << input_entries - both << '\n'
      << tab << "truth_only: " << truth_entries - both << '\n'
      << tab << "accuracy: " << accuracy << '\n'
      << "single:\n"
      << tab << "superoptimal_uniq: " << superoptimal_uniq << '\n'
      << tab << "optimal_uniq: " << optimal_uniq << '\n'
      << tab << "unknown_uniq: " << unknown_uniq << '\n'
      << tab << "suboptimal_uniq: " << suboptimal_uniq << '\n'
      << tab << "false_uniq: " << false_uniq << '\n'
      << tab << "missed_uniq: " << missed_uniq << '\n'
      << tab << "dangling_uniq: " << dangling_uniq << '\n'
      << "pairs:\n"
      << tab << "paired_truth_ambig: " << paired_truth_ambig << '\n'
      << tab << "paired_unknown: " << paired_unknown << '\n'
      << tab << "unpaired: " << unpaired << '\n';
  return out.str();
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
is_same_strand(const sam_rec &a, const sam_rec &b) {
  return check_flag(a, samflags::read_rc) ==
         check_flag(b, samflags::read_rc);
}

bool
is_same_pos(const sam_rec &a, const sam_rec &b) {
  return (a.rname == b.rname) &&
         // GS: account for possible soft clipping differences
         (abs(static_cast<long long int>(a.pos) - static_cast<long long int>(b.pos) < 30)) &&
         is_same_strand(a, b);
}

bool
is_paired(const sam_rec &a) {
  return a.rnext == "=";
}

void
mapstats::update(const sam_rec &input, const sam_rec &truth) {
  // get score from NM:i: tag
  const size_t diffs_input = get_edit_distance(input);
  const size_t diffs_truth = get_edit_distance(truth);

  // check for same chrom, pos and strand
  const bool same_pos = is_same_pos(input, truth);

  // abismal prints secondary aln on non-uniqs
  const bool is_truth_uniq = !check_flag(truth, samflags::secondary_aln);

  const bool is_truth_paired = is_paired(truth);
  const bool is_input_paired = is_paired(input);

  // both paired or both unpaired
  if (is_truth_paired == is_input_paired) {
    // if truth is better than input, either suboptimal or false uniq
    if (diffs_truth < diffs_input) {

      // these are disjoint
      optimal_uniq += same_pos;
      suboptimal_uniq += (!same_pos && is_truth_uniq);
      false_uniq += (!same_pos && !is_truth_uniq);
    }

    // if truth is equal to input, either optimal or unknown uniq
    else if (diffs_truth == diffs_input) {
      false_uniq += !is_truth_uniq;
      optimal_uniq += same_pos && is_truth_uniq;
      unknown_uniq += !same_pos && is_truth_uniq;
    } else {
      // GS enforcing not same pos just in case the edit distance is
      // somehow calculated differently
      superoptimal_uniq += !same_pos;
      optimal_uniq += same_pos;
      /*
      if (!same_pos) {
        print_debug(input, truth);
        write_read_for_debug(input, cout);
      }*/
    }
  }

  // update stats on reads with different pairing than truth
  else {
    // separate cases where mate is unpaired because ambiguous or
    // because the other mate is unmapped for an unknown reason
    if (is_input_paired && !is_truth_paired) {
      paired_truth_ambig += !is_truth_uniq;
      paired_unknown += is_truth_uniq;
    }
    unpaired += (!is_input_paired && is_truth_paired);
  }
}

int
main(int argc, const char **argv) {
  bool verbose = false;
  OptionParser opt_parse(strip_path(argv[0]),
                         "finds better matches between sam files sorted by name",
                        "<sam-truth> <sam-input>");
  vector<string> leftover_args;

  opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
  opt_parse.parse(argc, argv, leftover_args);
  if (leftover_args.size() != 2) {
    cerr << "Need two arguments: <sam-truth> <sam-input>\n";
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }
  SAMReader in_truth(leftover_args[0]), in_input(leftover_args[1]);

  sam_rec aln_truth, aln_input;
  mapstats stats;

  string last_input_qname = "";
  string last_truth_qname = "";
  while (in_input >> aln_input) {
    if (aln_input.qname < last_input_qname)
      throw runtime_error("input not sorted by qname\n");

    ++stats.input_entries;
    bool stop_searching = false;
    bool found_matching_read = false;

    if (aln_truth.qname == aln_input.qname) {
      ++stats.both;
      stats.update(aln_input, aln_truth);
      found_matching_read = true;
      stop_searching = true;
    }

    // truth ahead, we need to keep reading input
    else if (aln_truth.qname > aln_input.qname)
      stop_searching = true;

    // truth behind, need to keep searching
    while (!stop_searching) {
      if (in_truth >> aln_truth) {
        if (aln_truth.qname < last_truth_qname)
          throw runtime_error("truth not sorted by qname\n");

        ++stats.truth_entries;
        // truth hasn't reached out to the current read. Reads don't
        // exist in input, missed if unique
        if (aln_truth.qname < aln_input.qname)
          stats.missed_uniq += !check_flag(aln_truth, samflags::secondary_aln);

        // potentially the same read, process their differences
        else if (aln_truth.qname == aln_input.qname) {
          ++stats.both;
          found_matching_read = true;
          stop_searching = true;
          stats.update(aln_input, aln_truth);
        }

        // truth ahead, stop searching for input and read a new line
        else
          stop_searching = true;

        last_input_qname = aln_input.qname;
      }
      else stop_searching = true; //EOF
    }
    stats.dangling_uniq += (!found_matching_read);
    last_truth_qname = aln_truth.qname;
  }

  cout << stats.tostring() << endl;
  return EXIT_SUCCESS;
}
