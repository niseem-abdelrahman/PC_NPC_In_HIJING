#pragma once

#include "GeneralQCorrelator.h"

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <utility>
#include <vector>

namespace wrapper_qcorr {

class WrapperQCorrEvent {
public:
  using Base = qcorr::GeneralQCorrelator;
  using Result = Base::Result;
  using Subevent = Base::Subevent;
  using Mask = qcorr::Mask;

  static constexpr int kMaxPtBins = 27;
  static constexpr Mask kRefBit = 1ULL << 0;

  explicit WrapperQCorrEvent(double eta_gap = 0.7, int n_pt_bins = kMaxPtBins)
    : calc_(eta_gap),
      n_pt_bins_(n_pt_bins)
  {
    if (n_pt_bins_ < 1 || n_pt_bins_ > kMaxPtBins) {
      throw std::invalid_argument("WrapperQCorrEvent: n_pt_bins must satisfy 1 <= n_pt_bins <= 27.");
    }
  }

  void clear_event() { calc_.clear_event(); }
  void reserve(std::size_t n) { calc_.reserve(n); }

  // Reserve room for event-level block sums. This avoids repeated rehashing in large analyses.
  void reserve_block_cache(std::size_t n) const { calc_.reserve_block_cache(n); }

  void set_eta_gap(double eta_gap) { calc_.set_eta_gap(eta_gap); }
  [[nodiscard]] double eta_gap() const { return calc_.eta_gap(); }

  [[nodiscard]] int n_pt_bins() const { return n_pt_bins_; }
  [[nodiscard]] std::size_t size() const { return calc_.size(); }

  static constexpr Mask ref_bit() { return kRefBit; }

  static constexpr Mask poiA_bit(int pt_bin) {
    return 1ULL << (1 + pt_bin);
  }

  static constexpr Mask poiB_bit(int pt_bin) {
    return 1ULL << (1 + kMaxPtBins + pt_bin);
  }

  static std::vector<bool> gmc_none(int k) {
    return std::vector<bool>(static_cast<std::size_t>(k), false);
  }

  static std::vector<Subevent> all_any(int k) {
    return std::vector<Subevent>(static_cast<std::size_t>(k), Subevent::Any);
  }

  void add_track(double phi,
                 double eta,
                 double pt,
                 int pt_bin,
                 double acc_weight,
                 double gmc_weight,
                 bool in_reference,
                 bool in_poiA = true,
                 bool in_poiB = true) {
    if (pt_bin < 0 || pt_bin >= n_pt_bins_) {
      throw std::out_of_range("WrapperQCorrEvent::add_track: pt_bin out of range.");
    }

    Base::Track t;
    t.phi = phi;
    t.eta = eta;
    t.pt = pt;
    t.acc_weight = acc_weight;
    t.gmc_weight = gmc_weight;
    t.categories = 0;

    if (in_reference) {
      t.categories |= kRefBit;
    }
    if (in_poiA) {
      t.categories |= poiA_bit(pt_bin);
    }
    if (in_poiB) {
      t.categories |= poiB_bit(pt_bin);
    }

    calc_.add_track(t);
  }

  Result integrated_same(const std::vector<int>& harmonics,
                         const std::vector<bool>& use_gmc_weight = {}) const {
    return calc_.same_event(build_integrated_legs(harmonics, use_gmc_weight,
                                                  all_any(static_cast<int>(harmonics.size()))));
  }

  Result integrated_two_subevent(const std::vector<int>& harmonics,
                                 const std::vector<Subevent>& subevents,
                                 const std::vector<bool>& use_gmc_weight = {}) const {
    return calc_.two_subevent(build_integrated_legs(harmonics, use_gmc_weight, subevents));
  }

  Result differential_one_same(const std::vector<int>& harmonics,
                               const std::vector<int>& poi_legs,
                               int pt_bin,
                               const std::vector<bool>& use_gmc_weight = {}) const {
    return calc_.same_event(build_differential_one_legs(harmonics, poi_legs, pt_bin,
                                                        use_gmc_weight,
                                                        all_any(static_cast<int>(harmonics.size()))));
  }

  Result differential_one_two_subevent(const std::vector<int>& harmonics,
                                       const std::vector<int>& poi_legs,
                                       int pt_bin,
                                       const std::vector<Subevent>& subevents,
                                       const std::vector<bool>& use_gmc_weight = {}) const {
    return calc_.two_subevent(build_differential_one_legs(harmonics, poi_legs, pt_bin,
                                                          use_gmc_weight, subevents));
  }

  Result differential_two_same(const std::vector<int>& harmonics,
                               const std::vector<int>& poiA_legs,
                               int pt_bin_a,
                               const std::vector<int>& poiB_legs,
                               int pt_bin_b,
                               const std::vector<bool>& use_gmc_weight = {}) const {
    return calc_.same_event(build_differential_two_legs(harmonics, poiA_legs, pt_bin_a,
                                                        poiB_legs, pt_bin_b, use_gmc_weight,
                                                        all_any(static_cast<int>(harmonics.size()))));
  }

  Result differential_two_two_subevent(const std::vector<int>& harmonics,
                                       const std::vector<int>& poiA_legs,
                                       int pt_bin_a,
                                       const std::vector<int>& poiB_legs,
                                       int pt_bin_b,
                                       const std::vector<Subevent>& subevents,
                                       const std::vector<bool>& use_gmc_weight = {}) const {
    return calc_.two_subevent(build_differential_two_legs(harmonics, poiA_legs, pt_bin_a,
                                                          poiB_legs, pt_bin_b,
                                                          use_gmc_weight, subevents));
  }

  // ------------------------------------------------------
  // Common paper-style wrappers
  // ------------------------------------------------------
  Result d00_same(int n1, int n2,
                  bool gmc_leg1 = false, bool gmc_leg2 = false) const {
    return integrated_same({n1, n2}, {gmc_leg1, gmc_leg2});
  }

  Result d00_gap(int n1, int n2,
                 bool gmc_leg1 = false, bool gmc_leg2 = false) const {
    return integrated_two_subevent({n1, n2}, {Subevent::A, Subevent::B}, {gmc_leg1, gmc_leg2});
  }

  Result d100(int n1, int n2, int n3, int pt_bin,
              bool g1 = false, bool g2 = false, bool g3 = false) const {
    return differential_one_same({n1, n2, n3}, {0}, pt_bin, {g1, g2, g3});
  }

  Result d100_gap(int n1, int n2, int n3, int pt_bin,
                  bool g1 = false, bool g2 = false, bool g3 = false) const {
    return differential_one_two_subevent({n1, n2, n3}, {0}, pt_bin,
                                         {Subevent::A, Subevent::A, Subevent::B},
                                         {g1, g2, g3});
  }

  Result d1200_gap(int n1, int n2, int n3, int n4,
                   int pt_bin_a, int pt_bin_b,
                   bool g1 = false, bool g2 = false,
                   bool g3 = false, bool g4 = false) const {
    return differential_two_two_subevent({n1, n2, n3, n4},
                                         {0}, pt_bin_a,
                                         {1}, pt_bin_b,
                                         {Subevent::A, Subevent::A, Subevent::B, Subevent::B},
                                         {g1, g2, g3, g4});
  }

private:
  Base calc_;
  int n_pt_bins_;

  void check_shapes(const std::vector<int>& harmonics,
                    const std::vector<bool>& use_gmc_weight,
                    const std::vector<Subevent>& subevents) const {
    const std::size_t k = harmonics.size();
    if (k < 2 || k > 6) {
      throw std::invalid_argument("WrapperQCorrEvent: only k=2..6 is supported.");
    }
    if (!use_gmc_weight.empty() && use_gmc_weight.size() != k) {
      throw std::invalid_argument("WrapperQCorrEvent: use_gmc_weight size mismatch.");
    }
    if (!subevents.empty() && subevents.size() != k) {
      throw std::invalid_argument("WrapperQCorrEvent: subevents size mismatch.");
    }
  }

  static bool contains_index(const std::vector<int>& where, int idx) {
    return std::find(where.begin(), where.end(), idx) != where.end();
  }

  std::vector<Base::Leg> build_integrated_legs(const std::vector<int>& harmonics,
                                               const std::vector<bool>& use_gmc_weight,
                                               const std::vector<Subevent>& subevents) const {
    check_shapes(harmonics, use_gmc_weight, subevents);
    std::vector<Base::Leg> legs;
    legs.reserve(harmonics.size());

    for (std::size_t i = 0; i < harmonics.size(); ++i) {
      Base::Leg leg;
      leg.harmonic = harmonics[i];
      leg.required_categories = kRefBit;
      leg.subevent = subevents.empty() ? Subevent::Any : subevents[i];
      leg.use_gmc_weight = use_gmc_weight.empty() ? false : use_gmc_weight[i];
      legs.push_back(leg);
    }
    return legs;
  }

  std::vector<Base::Leg> build_differential_one_legs(const std::vector<int>& harmonics,
                                                     const std::vector<int>& poi_legs,
                                                     int pt_bin,
                                                     const std::vector<bool>& use_gmc_weight,
                                                     const std::vector<Subevent>& subevents) const {
    check_shapes(harmonics, use_gmc_weight, subevents);
    if (pt_bin < 0 || pt_bin >= n_pt_bins_) {
      throw std::out_of_range("WrapperQCorrEvent: pt_bin out of range in differential_one.");
    }

    const Mask poi_bit = poiA_bit(pt_bin);
    std::vector<Base::Leg> legs;
    legs.reserve(harmonics.size());

    for (std::size_t i = 0; i < harmonics.size(); ++i) {
      Base::Leg leg;
      leg.harmonic = harmonics[i];
      leg.required_categories = contains_index(poi_legs, static_cast<int>(i)) ? poi_bit : kRefBit;
      leg.subevent = subevents.empty() ? Subevent::Any : subevents[i];
      leg.use_gmc_weight = use_gmc_weight.empty() ? false : use_gmc_weight[i];
      legs.push_back(leg);
    }
    return legs;
  }

  std::vector<Base::Leg> build_differential_two_legs(const std::vector<int>& harmonics,
                                                     const std::vector<int>& poiA_legs,
                                                     int pt_bin_a,
                                                     const std::vector<int>& poiB_legs,
                                                     int pt_bin_b,
                                                     const std::vector<bool>& use_gmc_weight,
                                                     const std::vector<Subevent>& subevents) const {
    check_shapes(harmonics, use_gmc_weight, subevents);
    if (pt_bin_a < 0 || pt_bin_a >= n_pt_bins_ || pt_bin_b < 0 || pt_bin_b >= n_pt_bins_) {
      throw std::out_of_range("WrapperQCorrEvent: pt_bin out of range in differential_two.");
    }

    const Mask bit_a = poiA_bit(pt_bin_a);
    const Mask bit_b = poiB_bit(pt_bin_b);

    std::vector<Base::Leg> legs;
    legs.reserve(harmonics.size());

    for (std::size_t i = 0; i < harmonics.size(); ++i) {
      const bool is_a = contains_index(poiA_legs, static_cast<int>(i));
      const bool is_b = contains_index(poiB_legs, static_cast<int>(i));
      if (is_a && is_b) {
        throw std::invalid_argument("WrapperQCorrEvent: a leg cannot belong to both poiA_legs and poiB_legs.");
      }

      Base::Leg leg;
      leg.harmonic = harmonics[i];
      if (is_a) {
        leg.required_categories = bit_a;
      } else if (is_b) {
        leg.required_categories = bit_b;
      } else {
        leg.required_categories = kRefBit;
      }
      leg.subevent = subevents.empty() ? Subevent::Any : subevents[i];
      leg.use_gmc_weight = use_gmc_weight.empty() ? false : use_gmc_weight[i];
      legs.push_back(leg);
    }
    return legs;
  }
};

}  // namespace wrapper_qcorr
