#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <random>
#include <algorithm>

#include <omp.h>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::ostream;
using std::runtime_error;
using std::transform;

typedef vector<uint8_t> Genome;
typedef vector<uint8_t> Read;
typedef int16_t score_t;

static const char encode_8bit = {};
static const char encode_4bit = {};

inline char random_base() {return "ACGT"[rand() % 4];}
inline double rand_double() { return ((double) rand() / (RAND_MAX)); }

static void
make_random_genome(Genome &g, const size_t genome_size) {
  g.resize(genome_size);
  for (auto it(begin(g)); it != end(g); ++it)
    *it = random_base();
}

static void
make_random_read(Read &r, const size_t read_size) {
  r.resize(read_size);
  for(auto it(begin(r)); it != end(r); ++it)
    *it = random_base();
}

static void
sample_positions(vector<uint32_t> &positions, const size_t num_pos,
                 const size_t lim) {
  positions.resize(num_pos);
  for (auto it(begin(positions)); it != end(positions); ++it)
    *it = static_cast<uint32_t>(rand()) % lim;
}

template <typename T>
T& operator << (T &out, const vector<uint8_t> &r) {
  for (auto it(begin(r)); it != end(r);  ++it)
    out << *it;
  return out;
}

template <typename T>
T& operator << (T &out, const vector<uint32_t> &pos) {
  for (auto it(begin(pos)); it != end(pos);  ++it)
    out << *it << ' ';
  return out;
}

inline score_t
full_compare_unencoded(Read::const_iterator read_st,
             const Read::const_iterator read_end,
             const score_t cutoff,
             Genome::const_iterator g) {
  score_t d = 0;
  for (; (read_st != read_end) && (d <= cutoff);) {
    d += (*read_st != *g);
    ++g, ++read_st;
  }
  return d;
}


static void
time_mismatch_unencoded(const Read &r, const Genome &g,
                        const vector<uint32_t> &pos) {
  score_t lim = r.size();
  uint32_t best_pos = 0;
  const auto genome_st(begin(g));
  const auto read_st(begin(r));
  const auto read_end(end(r));
  for (auto it(begin(pos)); it != end(pos); ++it) {
    const score_t diffs =
      full_compare_unencoded(read_st, read_end, lim, genome_st + *it);
    if (diffs < lim) {
      lim = diffs;
      best_pos = *it;
    }
  }
  cerr << "[best mismatch: " << lim << " best pos: " << best_pos << "]\n";
}

inline score_t
full_compare_eight_bit(Read::const_iterator read_st,
             const Read::const_iterator read_end,
             const score_t cutoff,
             Genome::const_iterator g) {
  score_t d = 0;
  for (; (read_st != read_end) && (d <= cutoff);) {
    d += ((*read_st & *g) == 0);
    ++g, ++read_st;
  }
  return d;
}


static void
time_mismatch_eight_bit(const Read &r, const Genome &g,
                        const vector<uint32_t> &pos) {
  score_t lim = r.size();
  uint32_t best_pos = 0;
  const auto genome_st(begin(g));
  const auto read_st(begin(r));
  const auto read_end(end(r));
  for (auto it(begin(pos)); it != end(pos); ++it) {
    const score_t diffs =
      full_compare_eight_bit(read_st, read_end, lim, genome_st + *it);
    if (diffs < lim) {
      lim = diffs;
      best_pos = *it;
    }
  }
  cerr << "[best mismatch: " << lim << " best pos: " << best_pos << "]\n";
}




static void
encode_eight_bit(vector<uint8_t> &v) {
  transform(begin(v), end(v), begin(v), [](const uint8_t c) {
      if (c == 'A') return 1;
      if (c == 'C') return 2;
      if (c == 'G') return 4;
      if (c == 'T') return 8;
      return 0;
  });
}

static void
pack_genome(const Genome &g, vector<size_t> &gp) {
  const size_t packed_size = g.size() / 16;
  gp.resize(packed_size);
  auto ig(begin(g));
  for (auto it(begin(gp)); it != end(gp); ++it) {
    *it = ((*g) | ((*++g) << 4) | ((*++g) << 8) | ((*++g) << 12) |
          ((*++g) << 16) | ((*++g) << 20) | ((*++g) << 24) | ((*++g) << 28) |
          ((*++g) << 32) | ((*++g) << 36) | ((*++g) << 40) | ((*++g) << 44) |
          ((*++g) << 48) | ((*++g) << 52) | ((*++g) << 56) | ((*++g) << 60));
    ++it;
  }
}

int
main(int argc, const char **argv) {
  bool verbose = false;
  size_t genome_size = 1e7;
  size_t read_size = 100;
  size_t num_hits = 10000;
  OptionParser opt_parse(strip_path(argv[0]),
                         "finds better matches between sam files sorted by name",
                        "<sam-input>");
  vector<string> leftover_args;
  opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
  opt_parse.add_opt("genome-size", 'N', "genome size", false, genome_size);
  opt_parse.add_opt("read-size", 'n', "read size", false, read_size);
  opt_parse.add_opt("num-hits", 'h', "number of hits", false, num_hits);
  opt_parse.parse(argc, argv, leftover_args);
  
  if (leftover_args.size() != 0) {
    cerr << "Need two arguments: <sam-truth> <sam-input>\n";
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }
  if (opt_parse.help_requested()) {
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }

  Read r;
  Genome genome;
  vector<uint32_t> positions;
  if (verbose)
    cerr << "[sampling genome, read and positions]\n";
  make_random_genome(genome, genome_size);
  make_random_read(r, read_size);
  sample_positions(positions, num_hits, genome.size() - r.size());

  /***************** UNCONVERTED *******************/
  if (verbose)
    cerr << "[simulating mismatch count with characters]\n";
  double time = omp_get_wtime();
  time_mismatch_unencoded(r, genome, positions); 
  cerr << "[time elapsed: " << (omp_get_wtime() - time) << "]\n";

  /***************** ENCODED WITH 8 BITS ************/
  encode_eight_bit(genome);
  encode_eight_bit(r);
  time = omp_get_wtime();
  time_mismatch_eight_bit(r, genome, positions);
  cerr << "[time elapsed: " << (omp_get_wtime() - time) << "]\n";

  /***************** ENCODED WITH 8 BITS ************/
  vector<size_t> genome_packed;
  pack_genome(genome, genome_packed);
}
