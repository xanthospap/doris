#include "doris_rinex.hpp"

int dso::RinexDataBlockIterator::next() noexcept {
  char line[DorisObsRinex::MAX_RECORD_CHARS];

  // try getting the next line from the RINEX stream ...
  if (!rnx->stream().getline(line, DorisObsRinex::MAX_RECORD_CHARS)) {
    if (rnx->stream().eof())
      return -1;
    fprintf(stderr,
            "[ERROR] Failed to read lines from RINEX file (traceback: %s)\n",
            __func__);
    return 1;
  }

  // ... this should be a data record line; resolve it
  if (int error = rnx->resolve_data_epoch(line, cheader); error) {
    fprintf(
        stderr,
        "[ERROR] Failed parsing data header line, error=%d (traceback: %s)\n",
        error, __func__);
    fprintf(stderr, "Line is: \"%s\"\n", line);
    return 2;
  }

  // get block of observations that follow ...
  if (int error = rnx->read_data_block(cheader, cblock); error) {
    fprintf(stderr, "[ERROR] Failed parsing data block line! (traceback: %s)\n",
            __func__);
    return 3;
  }

  return 0;
}

dso::datetime<dso::nanoseconds> dso::RinexDataBlockIterator::corrected_l1_epoch(
    dso::datetime<dso::nanoseconds> &tl2) const noexcept {
  dso::datetime<dso::nanoseconds> t = cheader.m_epoch;
  if (!rnx->receiver_clock_offsets_applied()) [[likely]] {
    // clock offset in nanoseconds (from seconds)
    dso::nanoseconds noff(static_cast<dso::nanoseconds::underlying_type>(
        cheader.m_clock_offset * 1e9));
    t.add_seconds(noff);
  }

  // get the L1/L2 offset in nanoseconds
  const auto l12_offset = rnx->l21_date_offset();

  // correct tl2
  tl2 = t;
  tl2.add_seconds(l12_offset);

  // reference epoch for L1
  return t;
}
