#include "tides.hpp"
#include <stdexcept>

dso::DoodsonOceanTideConstituent::DoodsonOceanTideConstituent(
    const dso::DoodsonNumber d, int max_degree, int max_order)
    : doodson(d), maxl(max_degree), maxm(max_order),
      DelCpl(new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
          max_degree+1, max_degree+1)),
      DelSpl(new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
          max_degree+1, max_degree+1)),
      DelCmi(new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
          max_degree+1, max_degree+1)),
      DelSmi(new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
          max_degree+1, max_degree+1)) {
  if (max_order > max_degree)
    throw std::runtime_error("[ERROR] Invalid degree/order sizes for "
                             "dso::DoodsonOceanTideConstituent\n");
}

void dso::DoodsonOceanTideConstituent::set_null() noexcept {
  if (DelCpl) {
    delete DelCpl;
  }
  DelCpl = nullptr;
  if (DelSpl) {
    delete DelSpl;
  }
  DelSpl = nullptr;
  if (DelCmi) {
    delete DelCmi;
  }
  DelCmi = nullptr;
  if (DelSmi) {
    delete DelSmi;
  }
  DelSmi = nullptr;
  maxl = maxm = 0;
}

dso::DoodsonOceanTideConstituent::~DoodsonOceanTideConstituent() noexcept {
  set_null();
}

dso::DoodsonOceanTideConstituent::DoodsonOceanTideConstituent(
    const dso::DoodsonOceanTideConstituent &other) noexcept {
  if (other.maxl != 0) {
    // use one of the matrices to check sizes. Asserts that all matrices are
    // of equal size!
    if (this->DelCpl->num_elements() != other.DelCpl->num_elements()) {
      this->set_null();
      assert(this->resize(other.maxl) == dso::iStatus::ok());
    }
    *DelCpl = *(other.DelCpl);
    *DelSpl = *(other.DelSpl);
    *DelCmi = *(other.DelCmi);
    *DelSmi = *(other.DelSmi);
  } else {
    assert(other.maxm == 0);
    this->set_null();
  }
  doodson = other.doodson;
}

dso::DoodsonOceanTideConstituent &dso::DoodsonOceanTideConstituent::operator=(
    const dso::DoodsonOceanTideConstituent &other) noexcept {
  if (this != &other) {
    if (other.maxl != 0) {
      // use one of the matrices to check sizes. Asserts that all matrices are
      // of equal size!
      if (this->DelCpl->num_elements() != other.DelCpl->num_elements()) {
        this->set_null();
        assert(this->resize(other.maxl) == dso::iStatus::ok());
      }
      *DelCpl = *(other.DelCpl);
      *DelSpl = *(other.DelSpl);
      *DelCmi = *(other.DelCmi);
      *DelSmi = *(other.DelSmi);
    } else {
      assert(other.maxm == 0);
      this->set_null();
    }
    doodson = other.doodson;
  }
  return *this;
}

dso::DoodsonOceanTideConstituent::DoodsonOceanTideConstituent(
    dso::DoodsonOceanTideConstituent &&other) noexcept
    : doodson(std::move(other.doodson)), maxl(other.maxl), maxm(other.maxm) {
  this->set_null();
  // just copy the pointers!
  DelCpl = other.DelCpl;
  DelSpl = other.DelSpl;
  DelCmi = other.DelCmi;
  DelSmi = other.DelSmi;
  other.set_null();
}

dso::DoodsonOceanTideConstituent &dso::DoodsonOceanTideConstituent::operator=(
    dso::DoodsonOceanTideConstituent &&other) noexcept {
  doodson = std::move(other.doodson);
  maxl = other.maxl;
  maxm = other.maxm;
  this->set_null();
  DelCpl = other.DelCpl;
  DelSpl = other.DelSpl;
  DelCmi = other.DelCmi;
  DelSmi = other.DelSmi;
  other.set_null();
  return *this;
}

dso::iStatus dso::DoodsonOceanTideConstituent::resize(int maxDegree) noexcept {
  assert(maxDegree > 0);
  ++maxDegree;
  try {
    if (DelCpl) {
      DelCpl->resize(maxDegree, maxDegree);
    } else {
      DelCpl = new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
          maxDegree, maxDegree);
    }
    if (DelSpl) {
      DelSpl->resize(maxDegree, maxDegree);
    } else {
      DelSpl = new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
          maxDegree, maxDegree);
    }
    if (DelCmi) {
      DelCmi->resize(maxDegree, maxDegree);
    } else {
      DelCmi = new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
          maxDegree, maxDegree);
    }
    if (DelSmi) {
      DelSmi->resize(maxDegree, maxDegree);
    } else {
      DelSmi = new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
          maxDegree, maxDegree);
    }
  } catch (std::exception &) {
    return dso::iStatus(1);
  }
  return dso::iStatus::ok();
}
