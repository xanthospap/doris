#include "tides.hpp"
#include <stdexcept>

dso::DoodsonOceanTideConstituent::DoodsonOceanTideConstituent(
    const dso::DoodsonNumber d, int max_degree, int max_order)
    : doodson(d), maxl(max_degree), maxm(max_order),
      DelCpl(new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
          max_degree + 1, max_degree + 1)),
      DelSpl(new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
          max_degree + 1, max_degree + 1)),
      DelCmi(new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
          max_degree + 1, max_degree + 1)),
      DelSmi(new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
          max_degree + 1, max_degree + 1)) {
  if (max_order > max_degree || max_degree <= 1)
    throw std::runtime_error("[ERROR] Invalid degree/order sizes for "
                             "dso::DoodsonOceanTideConstituent\n");
}

#ifdef DEBUG
void dso::DoodsonOceanTideConstituent::print_matrix_sizes() const noexcept {
  if (DelCpl) {
    printf("\tDeltaC+ size: %dx%d, elements: %lu\n", DelCpl->rows(),
           DelCpl->cols(), DelCpl->num_elements());
  } else {
    printf("\tDeltaC+ is NULL!\n");
  }
  if (DelSpl) {
    printf("\tDeltaS+ size: %dx%d, elements: %lu\n", DelSpl->rows(),
           DelSpl->cols(), DelSpl->num_elements());
  } else {
    printf("\tDeltaS+ is NULL!\n");
  }
  if (DelCmi) {
    printf("\tDeltaC- size: %dx%d, elements: %lu\n", DelCmi->rows(),
           DelCmi->cols(), DelCmi->num_elements());
  } else {
    printf("\tDeltaC- is NULL!\n");
  }
  if (DelSmi) {
    printf("\tDeltaS- size: %dx%d, elements: %lu\n", DelSmi->rows(),
           DelSmi->cols(), DelSmi->num_elements());
  } else {
    printf("\tDeltaS- is NULL!\n");
  }
}
#endif

void dso::DoodsonOceanTideConstituent::deallocate() noexcept {
  if (DelCpl)
    delete DelCpl;
  if (DelSpl)
    delete DelSpl;
  if (DelCmi)
    delete DelCmi;
  if (DelSmi)
    delete DelSmi;
}

void dso::DoodsonOceanTideConstituent::set_null() noexcept {
  DelCpl = nullptr;
  DelSpl = nullptr;
  DelCmi = nullptr;
  DelSmi = nullptr;
  maxl = maxm = 0;
}

dso::DoodsonOceanTideConstituent::~DoodsonOceanTideConstituent() noexcept {
  deallocate();
  set_null();
}

dso::DoodsonOceanTideConstituent::DoodsonOceanTideConstituent(
    const dso::DoodsonOceanTideConstituent &other) noexcept {
  if (other.maxl != 0) {
    // use one of the matrices to check sizes. Asserts that all matrices are
    // of equal size!
    if (this->DelCpl->num_elements() != other.DelCpl->num_elements()) {
      this->deallocate();
      this->set_null();
      assert(this->resize(other.maxl) == dso::iStatus::ok());
    }
    *DelCpl = *(other.DelCpl);
    *DelSpl = *(other.DelSpl);
    *DelCmi = *(other.DelCmi);
    *DelSmi = *(other.DelSmi);
  } else {
    assert(other.maxm == 0);
    this->deallocate();
    this->set_null();
  }
  maxl = other.maxl;
  maxm = other.maxm;
  doodson = other.doodson;
}

dso::DoodsonOceanTideConstituent &dso::DoodsonOceanTideConstituent::operator=(
    const dso::DoodsonOceanTideConstituent &other) noexcept {
  if (this != &other) {
    if (other.maxl != 0) {
      // use one of the matrices to check sizes. Asserts that all matrices are
      // of equal size!
      if (this->DelCpl->num_elements() != other.DelCpl->num_elements()) {
        this->deallocate();
        this->set_null();
        assert(this->resize(other.maxl) == dso::iStatus::ok());
      }
      *DelCpl = *(other.DelCpl);
      *DelSpl = *(other.DelSpl);
      *DelCmi = *(other.DelCmi);
      *DelSmi = *(other.DelSmi);
    } else {
      assert(other.maxm == 0);
      this->deallocate();
      this->set_null();
    }
    maxl = other.maxl;
    maxm = other.maxm;
    doodson = other.doodson;
  }
  return *this;
}

dso::DoodsonOceanTideConstituent::DoodsonOceanTideConstituent(
    dso::DoodsonOceanTideConstituent &&other) noexcept
    : doodson(std::move(other.doodson)), maxl(other.maxl), maxm(other.maxm) {
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
  this->deallocate();
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
