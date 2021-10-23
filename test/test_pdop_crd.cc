#include "doris_rinex.hpp"
#include "doris_utils.hpp"

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage %s [DORS RINEX] [SINEX]\n", argv[0]);
    return 1;
  }

  // aim ...
  printf("Extrapolating coordinates for all beacons found in RINEX %s using "
         "the SINEX file %s\n",
         argv[1], argv[2]);

  // initialize the rinex file (get beacon list)
  ids::DorisObsRinex rnx(argv[1]);

  // get a list of all beacons by their site id
  auto beacons = rnx.stations();
  printf(
      "\nDORIS RINEX read through; it contains observations from %lu beacons\n",
      beacons.size());

  // allocate and fill the list of site id's
  int num_sites = beacons.size();
  char **sites = new char *[num_sites];
  for (int i = 0; i < num_sites; i++) {
    sites[i] = new char[5];
    std::memset(sites[i], 0, 5);
    std::strncpy(sites[i], beacons[i].m_station_id, 4);
  }

  // etrapolate coordinates for the reference time of RINEX
  dso::datetime<dso::nanoseconds> t = rnx.ref_datetime();
  int error = ids::extrapolate_sinex_coordinates(argv[2], sites, num_sites, t);
  if (error)
    fprintf(stderr, "[ERROR] Failed to extrapolate coordinates! (error=%d)\n",
            error);

  // deallocate memmory
  for (int i = 0; i < num_sites; i++)
    delete[] sites[i];
  delete[] sites;

  return error;
}