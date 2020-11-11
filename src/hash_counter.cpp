#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <cmath>
#include <fstream>
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

using std::cerr;
using std::cin;
using std::cout;
using std::string;
using std::vector;
using std::runtime_error;
using std::endl;
using std::min;
using std::ifstream;
using std::ofstream;

enum base {a,c,g,t,n};
inline double
rand_double() {
  return ((double)rand())/RAND_MAX;
}
inline uint8_t randbase(const vector<double> &probs) {
  double runif = rand_double();
  if (runif <= probs[a]) return 'A';
  if (runif <= probs[c]) return 'C';
  if (runif <= probs[g]) return 'G';
  return 'T'; //unknown base
}

inline uint8_t
to_twobit(const uint8_t base) {
  if(base == 'A') return 0;
  if(base == 'C') return 1;
  if(base == 'G') return 2;
  return 3;
}

template <class G>
static void
load_genome(const std::string &genome_file, G &genome) {
  std::ifstream in(genome_file);
  if (!in)
    throw std::runtime_error("bad genome file: " + genome_file);

  const auto begin_pos = in.tellg();
  in.seekg(0, std::ios_base::end);
  const size_t file_size = in.tellg() - begin_pos;
  in.seekg(0, std::ios_base::beg);

  genome.clear();
  genome.reserve(file_size); // pad on at start; the

  std::string line;
  while (getline(in, line))
    if (line[0] != '>') {
      line.erase(std::remove_if(std::begin(line), std::end(line),
          [](unsigned char c) {
            return c != 'A' && c != 'C' && c != 'T' && c != 'G';
          }), std::end(line));
      copy(std::begin(line), std::end(line), std::back_inserter(genome));
    }
}

inline void
revcomp_genome(vector<uint8_t>& s) {
  std::transform(s.begin(), s.end(), s.begin(), complement);
  std::reverse(s.begin(), s.end());
}


namespace hash {
  size_t available_bits = 8;
  size_t k_two,
  k_three_poly,
  k_three_bit;

  size_t hash_two_mod,
         hash_three_mod_bit,
         hash_three_mod_poly;

  inline uint8_t
  hash_two(const uint8_t base) {
    return (to_twobit(base) & 1);
  } // A/G=0, C/T=1

  inline uint8_t
  hash_three(const uint8_t base) {
    if (base == 'C' || base == 'T') return 2;
    if (base == 'A') return 0;
    return 1;
  }
  double get_theoretical_two(const double pw, const double ps) {
    return pow(pw + ps, k_two);
  }
  
  double get_theoretical_three_poly(const double pw, const double ps) {
    return pow(0.25 + pw*pw + ps*ps, k_three_poly);
  }

  static size_t
  shift_hash_three_poly(const size_t hash, const uint8_t base) {
    return (((3*hash)%hash_three_mod_poly) +
             hash_three(base)) % hash_three_mod_poly;
  }

  static size_t
  shift_hash_two(const size_t hash, const uint8_t base) {
    return ((hash << 1) | hash_two(base)) & hash_two_mod;
  }
}

template <size_t (*shift_hash)(const size_t, const uint8_t)>
static void
count_buckets(const vector<uint8_t> &genome,
              const size_t key_weight,
              vector <size_t> &counter) {
  auto it = begin(genome);
  const auto end_first = it + min(genome.size(), key_weight);

  size_t hash = 0;
  size_t last_n = 0;
  size_t cnt = 0;
  for (; it != end_first; ++it) {
    ++cnt;
    //if (*it == 'N') last_n = cnt;
    hash = shift_hash(hash, *it);
  }
  if (cnt - last_n >= key_weight)
    ++counter[hash];

  for (; it != end(genome); ++it) {
    ++cnt;
   // if (*it == 'N') last_n = cnt;
    hash = shift_hash(hash, *it);
    if (cnt - last_n >= key_weight)
      ++counter[hash];
  }
}

static void
count_all(const vector<uint8_t> &genome,
          vector<size_t> &counter_two_bit,
          vector<size_t> &counter_three_poly,
          vector<size_t> &counter_three_bit) {
  cerr << "[counting two bit]\n";
  count_buckets<hash::shift_hash_two>(genome, hash::k_two, counter_two_bit);
  // hash three bit

  //count_buckets<hash::shift_hash_three_bit>(genome, hash::k_three_bit,
  //  counter_three_bit);

  cerr << "[counting three bits]\n";

  // hash three poly
  count_buckets<hash::shift_hash_three_poly>(genome, hash::k_three_poly,
    counter_three_poly);
  cerr << "[counting three poly]\n";
}

static double
get_hit_prob(const vector<size_t> &vals) {
  size_t den = 0;
  double ans = 0.0;
  for (auto it(begin(vals)); it != end(vals); ++it)
    den += *it;

  double den_d = double(den);
  for (auto it(begin(vals)); it != end(vals); ++it) {
    const double prob = *it / den_d;
    ans += prob*prob;
  }
  return ans;
}

int main(const int argc, const char **argv) {
  // genome size
  bool verbose = false;
  bool revcomp = false;
  size_t N = 100;
  string genome_file, outfile = "";
  /**************COMMAND ARGS**********************/
  OptionParser opt_parse(strip_path(argv[0]),
                        "analysis of collision probabilities", "<genome-file>");
  vector<string> leftover_args;
  opt_parse.add_opt("bits", 'b', "number of bits", true, hash::available_bits);
  opt_parse.add_opt("size", 'n', "genome size",false, N);
  opt_parse.add_opt("revcomp", 'r', "count buckets in forward+reverse strand",
                    false, revcomp);
  opt_parse.add_opt("verbose", 'v', "print more run info",false, verbose);
  opt_parse.add_opt("output", 'o', "output file for statistics", false,
                    outfile);

  opt_parse.parse(argc, argv, leftover_args);
  if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
  }
  if (opt_parse.option_missing()) {
    cerr << opt_parse.option_missing_message() << endl;
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }

  if (leftover_args.size() != 1) {
    cerr << "Please provide a genome file as leftover arg\n";
    cerr << opt_parse.option_missing_message() << endl;
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }
  genome_file = leftover_args[0];

  /**************END COMMAND ARGS**********************/

  // number of bits in each case
  hash::k_two = hash::available_bits;
  hash::k_three_poly = floor(hash::available_bits*log(2)/log(3));

  // maximum hash value in each case
  hash::hash_two_mod = (1 << hash::k_two) - 1ll;
  hash::hash_three_mod_poly = (size_t) (pow(3, hash::k_three_poly) + 0.5);

  // random genome
  vector<uint8_t> genome(N);
  if (verbose)
    cerr << "[loading custom genome]\n";
  load_genome(genome_file, genome);

  if (verbose)
    cerr << "[calculating pw and ps]\n";

  double pw =
    count_if(begin(genome), end(genome), 
      [](const uint8_t c) { return c == 'A' || c == 'T'; } ) /
    static_cast<double>(2 * genome.size());

  double ps = 0.5 - pw;

  if (verbose)
    cerr << "[pw: " << pw << " ps: " << ps << "]\n";

  // base frequency probabilities
  vector<double> probs(4, 0.0);
  probs[a] = probs[t] = pw;
  probs[c] = probs[g] = ps;

  cerr << "[available bits: " << hash::available_bits << "]\n";
  cerr << "[k2: " << hash::k_two << "]\n";
  //cerr << "[k three bit: " << hash::k_three_bit << "]\n";
  cerr << "[k3: " << hash::k_three_poly << "]\n";
  
  /************ COUNTING ****************************/
  vector <size_t> counter_two_bit,
                  counter_three_poly,
                  counter_three_bit;

  // allocate counter space
  counter_two_bit.resize(hash::hash_two_mod + 1);
  counter_three_poly.resize(hash::hash_three_mod_poly);

  // fill with zeros
  std::fill(begin(counter_two_bit), end(counter_two_bit), 0);
  std::fill(begin(counter_three_poly), end(counter_three_poly), 0);

  cerr << "[counting forward strand]\n";
  count_all(genome, counter_two_bit, counter_three_poly, counter_three_bit);
  if (revcomp) {
    cerr << "[reverse complimenting genome]\n";
    revcomp_genome(genome);
    cerr << "[counting reverse strand]\n";
    count_all(genome, counter_two_bit, counter_three_poly, counter_three_bit);
  }
  string k_str = std::to_string(hash::available_bits);

  if (!outfile.empty()) {
    ofstream os_two(outfile + ".twobit");
    for (size_t i = 0; i < counter_two_bit.size(); ++i)
      os_two << i << " " << counter_two_bit[i] << "\n";

    ofstream os_three_poly(outfile + ".threebit");
    for (size_t i = 0; i < counter_three_poly.size(); ++i)
      os_three_poly << i << " " << counter_three_poly[i] << "\n";

    ofstream stats(outfile);
    stats << "k\tk3\tp2\tp2_theo\tp3\tp3_theo\tp2/p3\tp2_theo/p3_theo"
          << "\tp2/p2_theo\tp3/p3_theo\t(p2/p2_theo)/(p3/p3_theo)\t0.93^b\n";
    double p2b = get_hit_prob(counter_two_bit),
           p3p = get_hit_prob(counter_three_poly);

    double p2b_theo = hash::get_theoretical_two(pw, ps);
    double p3p_theo = hash::get_theoretical_three_poly(pw, ps);

    // double p3b = get_hit_prob(counter_three_bit);
    stats << hash::available_bits << "\t"
          << hash::k_three_poly << "\t"
          << p2b << "\t"
          << p2b_theo << "\t"
          << p3p << "\t"
          << p3p_theo << "\t"
          << p2b/p3p << "\t"
          << p2b_theo / p3p_theo << "\t"
          << p2b / p2b_theo << "\t"
          << p3p / p3p_theo << "\t"
          << (p2b / p2b_theo) / (p3p / p3p_theo) << "\t"
          << pow(0.93, hash::available_bits);

  }
  return 0;
}
