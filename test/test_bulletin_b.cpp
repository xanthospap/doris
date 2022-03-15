#include "iers_bulletin.hpp"
#include <cassert>
#include <cstdio>
#include <random>

std::random_device rd;  // obtain a random number from hardware
std::mt19937 gen(rd()); // seed the generator
constexpr const int max_block_size = 100;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Error. Usage %s <BULLETIN B>\n", argv[0]);
    return 1;
  }

  // make an array of sectin1 blocks, to hold file records
  dso::IersBulletinB_Section1Block blocks[max_block_size];

  // initialize a Bulletin B instance
  dso::IersBulletinB bb(argv[1]);

  // collect all values ...
  int num_records = bb.parse_section1(blocks);
  assert(num_records > 0);
  printf("Number of Section 1 records collected: %d\n", num_records);

  // random mjd within the range of collected records
  std::uniform_int_distribution<int> distr(0, num_records - 1);
  int index = distr(gen);
  int error = bb.get_section1_at(blocks[index].mjd, blocks[max_block_size - 1]);
  assert(!error);
  assert(blocks[index].mjd == blocks[max_block_size - 1].mjd);
  assert(blocks[index].x == blocks[max_block_size - 1].x);
  assert(blocks[index].y == blocks[max_block_size - 1].y);
  assert(blocks[index].dut1 == blocks[max_block_size - 1].dut1);
  assert(blocks[index].dX == blocks[max_block_size - 1].dX);
  assert(blocks[index].dY == blocks[max_block_size - 1].dY);
  assert(blocks[index].xerr == blocks[max_block_size - 1].xerr);
  assert(blocks[index].yerr == blocks[max_block_size - 1].yerr);
  assert(blocks[index].dut1err == blocks[max_block_size - 1].dut1err);
  assert(blocks[index].dXerr == blocks[max_block_size - 1].dXerr);
  assert(blocks[index].dYerr == blocks[max_block_size - 1].dYerr);

  return 0;
}