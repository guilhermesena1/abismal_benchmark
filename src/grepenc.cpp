#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "dna_four_bit.hpp"

using std::string;
using std::endl;
using std::cerr;
using std::cout;
using std::vector;
using std::transform;
using std::copy;
using std::begin;
using std::end;
using std::distance;

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
    if (line[0] != '>')
      copy(std::begin(line), std::end(line), std::back_inserter(genome));
}

char
encode_twobit(const char c) {
  if (c == 'A' || c == 'G') return '0';
  if (c == 'C' || c == 'T') return  '1';
  return c;
}

char
encode_threebit(const char c) {
  if (c == 'A') return '0';
  if (c == 'G') return '1';
  if (c == 'C' || c == 'T') return  '2';
  return c;
}

void
encode_genome(const string &genome, string &encoded_genome,
              const size_t alphabet_size) {
  encoded_genome.clear();
  encoded_genome.reserve(genome.size());
  if (alphabet_size == 2)
    transform(begin(genome), end(genome), std::back_inserter(encoded_genome),
              encode_twobit);
  else if (alphabet_size == 3)
    transform(begin(genome), end(genome), std::back_inserter(encoded_genome),
              encode_threebit);
  else copy(begin(genome), end(genome), begin(encoded_genome));

}

void
find_all_occurrences(
                     const string &pattern,
                     const string &encoded_genome,
                     const string &genome,
                     const char strand) {
  cerr << "pattern: " << pattern << "\n";
  size_t found = encoded_genome.find(pattern, 0);
  while (found != string::npos) {
    cout << std::string(begin(genome) + found,
                        begin(genome) + found + pattern.size())
         << '\t' << found << '\t' << strand << '\n';
    found = encoded_genome.find(pattern, found + 1);
  }
}
int
main(int argc, const char **argv) {
  size_t alphabet_size = 4;
  bool verbose = false;
  bool skip_revcomp = false;
  OptionParser opt_parse(strip_path(argv[0]),
                        "find sequences with given encoding",
                        "<pattern> <genome-file>");
  vector<string> leftover_args;

  opt_parse.add_opt("verbose", 'v', "print more run info",false, verbose);
  opt_parse.add_opt("skip-revcomp", 's', "do not count in reverse complement",
                    false, skip_revcomp);
  opt_parse.add_opt("alphabet", 'a', "alphabet encoding (2 3 or 4)",
                    false, alphabet_size);

  opt_parse.parse(argc, argv, leftover_args);
  if (alphabet_size < 2 || alphabet_size > 4) {
    cerr << "invalid alphabet size: " << alphabet_size
         << ". Must be 2, 3 or 4\n";
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }
  if (leftover_args.size() != 2) {
    cerr << "Need two arguments: <pattern> <genome file>\n";
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }
  const string genome_file = leftover_args[1],
               pattern = leftover_args[0];
  string genome, encoded_genome;
  if (verbose) cerr << "[loading genome]\n";
  load_genome(genome_file, genome);
  transform(begin(genome), end(genome), begin(genome), toupper);

  if (verbose) cerr << "[genome size: " << genome.size() << "]\n";

  if (verbose) cerr << "[encoding genome in " << alphabet_size
                    << " letters]\n";
  encode_genome(genome, encoded_genome, alphabet_size);

  if (verbose) cerr << "[finding pattern in + strand]\n";
  find_all_occurrences(pattern, encoded_genome, genome, '+');

  if (!skip_revcomp) {
    if (verbose) cerr << "[revcomping genome]\n";
    revcomp_inplace(genome);

    if (verbose) cerr << "[encoding revcomped genome]\n";
    encode_genome(genome, encoded_genome, alphabet_size);

    if (verbose) cerr << "[finding pattern in - strand]\n";
    find_all_occurrences(pattern, encoded_genome, genome, '-');
  }
  return EXIT_SUCCESS;
}
