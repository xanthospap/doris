#include "tides.hpp"

std::vector<dso::DoodsonOceanTideConstituent> dso::filter_waves(
    const std::vector<dso::DoodsonNumber> &waves,
    const std::vector<dso::DoodsonOceanTideConstituent> &freqs) noexcept {
  /* to be retuned ... */
  std::vector<dso::DoodsonOceanTideConstituent> vecnew;
  /* for every constintuent given */
  for (const auto &wave : waves) {
    /* we have coefficients for the constintuent */
    if (auto it = std::find_if(freqs.begin(), freqs.end(),
                               [&](const dso::DoodsonOceanTideConstituent &f) {
                                 return wave == f.doodson_number();
                               });
                               it != freqs.end()) {
      vecnew.emplace_back(*it);
    }
  }

  /* all done */
  return vecnew;
}
