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
get_edit_distance(const sam_rec &aln) {
  auto the_nm_tag = find_if(begin(aln.tags), end(aln.tags),
                            [](const string &t) {
                              return t.compare (0, 3, "NM:") == 0;
                            });

  // GS: NM:i:xxxx, get xxxx
  assert(the_nm_tag->size() > 5);
  return stoi(the_nm_tag->substr(5));
}

static int
fix_aln(sam_rec &aln) {
  auto the_nm_tag = find_if(begin(aln.tags), end(aln.tags),
                            [](const string &t) {
                              return t.compare (0, 3, "NM:") == 0;
                            });

  auto the_xm_tag = find_if(begin(aln.tags), end(aln.tags),
                            [](const string &t) {
                              return t.compare (0, 3, "XM:") == 0;
                            });
  auto the_xr_tag = find_if(begin(aln.tags), end(aln.tags),
                            [](const string &t) {
                              return t.compare (0, 3, "XR:") == 0;
                            });

  int edit_distance = stoi(the_nm_tag->substr(5));
  const string meth_call = the_xm_tag->substr(5);
  const string xr = *the_xr_tag;
  const int unmeth_c = std::count_if(begin(meth_call), end(meth_call),
       [](const char c) { return (c == 'x' || c == 'h' || c == 'z' || c == 'u'); });

  if (unmeth_c > edit_distance) {
    cerr << "this read has edit_distance " << edit_distance
         << " and unmeth_c " << unmeth_c << "\n";

    cerr << aln << endl;
    throw runtime_error("");
  }
  edit_distance -= unmeth_c;
  aln.tags.clear();
  aln.tags.push_back("NM:i:" + std::to_string(edit_distance));

  // need this for format_reads
  aln.tags.push_back(xr);

  return edit_distance;
}

static bool
has_mate(const sam_rec &aln) {
  return aln.rnext == "=" ||
         aln.rnext != "*";
}

static bool
are_mates(const sam_rec &one, const sam_rec &two) {
  return ((one.rnext == "=" && two.rnext == "=") ||
          (one.rnext == two.qname)) &&
           one.pnext == two.pos &&
           two.pnext == one.pos;
}

static bool
valid_edit_distance(const int edit_distance, const size_t read_size) {
  static double MAX_EDIT_DISTANCE_FRAC = 1.0;
  return (static_cast<double>(edit_distance) <= 
          MAX_EDIT_DISTANCE_FRAC * static_cast<double>(read_size));
}

static size_t readlen(const sam_rec &aln) {
  const int ans = cigar_qseq_ops(aln.cigar) -
                  get_soft_clip_size(aln.cigar);

  assert(ans >= 0);
  return static_cast<size_t>(ans);
}

static bool
report_aln(const sam_rec &aln, const int edit_distance) {
  return valid_edit_distance(edit_distance, readlen(aln));
}

static bool
report_mates(const sam_rec &a, const sam_rec &b,
             const int ea, const int eb) {
  return valid_edit_distance(ea + eb, readlen(a) + readlen(b));
}

static bool
report_se_read(const sam_rec &a, const int e) {
  return valid_edit_distance(e, readlen(a));
}


static void
fix_mate_name(sam_rec &mate) {
  assert(mate.qname.back() == '1');
  mate.qname.back() = '2';
}

static void
convert_to_se_aln(sam_rec &aln) {
  aln.rnext = "*";
  aln.pnext = 0;
  aln.tlen = 0;
}

static void
flip_xr_tag(sam_rec &aln) {
  auto xr_tag_itr = find_if(begin(aln.tags), end(aln.tags),
                               [](const string &t) {
                                 return t.compare (0, 3, "XR:") == 0;
                               });
  const bool a_rich = (xr_tag_itr->back() == 'A');
  xr_tag_itr->pop_back();
  xr_tag_itr->pop_back();

  *xr_tag_itr += (a_rich) ? "CT" : "GA";

}


int
main(int argc, const char **argv) {
  bool verbose = false;
  bool pbat = false;
  bool single_end = false;
  string outfile = "";
  OptionParser opt_parse(strip_path(argv[0]),
                         "finds better matches between sam files sorted by name",
                        "<sam-input>");
  vector<string> leftover_args;

  opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
  opt_parse.add_opt("pbat", 'p', "data is PBAT", false, pbat);
  opt_parse.add_opt("single-end", 's', "filter reads in SE mode", false, single_end);
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
  if (!single_end) {
    while (in >> aln) {
      const int d_aln = fix_aln(aln);
      if (has_mate(aln)) {
        if (in >> mate) {
          fix_mate_name(mate);
          if (!are_mates(aln, mate)) {
            cerr << "reads below are not mates\n";
            cerr << aln << '\n';
            cerr << mate << '\n';
            throw runtime_error("");
          }
        }  
        const int d_mate = fix_aln(mate);
        if (report_mates(aln, mate, d_aln, d_mate)) {
          out << aln << '\n';
          out << mate << '\n';
        } else {
          if (report_aln(aln, d_aln)) {
            convert_to_se_aln(aln);
            out << aln << '\n';
          }
          if (report_aln(mate, d_mate)) {
            convert_to_se_aln(mate);
            out << mate << '\n';
          }
        }
      } else {
        cerr << "found unpaired mate, which bismark does not do by default.\n";
        cerr << "Did you mean to run it on SE mode? (-s)\n";
        cerr << aln << endl;
        throw runtime_error("");
      }
    }
  } else {
    while (in >> aln) {
      const int d_aln = fix_aln(aln);
      if (report_se_read(aln, d_aln)) {
        if (pbat)
          flip_xr_tag(aln);
        out << aln << '\n';
      }
    }
  }

  return EXIT_SUCCESS;
}
