#pragma once
#include <array>
#include <climits>

template<int N>
class Bitmask {
private:
    using SubBitmask = uint64_t;
    static constexpr SubBitmask kOne = 1;
    static constexpr uint8_t kSubSize = sizeof(SubBitmask) * CHAR_BIT;
    static constexpr uint8_t kLogSubSize = 6;
    static_assert(kSubSize == (1 << kLogSubSize), "Inconsistent Bitmask");
    
    static inline SubBitmask bitFor(int i) {
        while (i > kSubSize) {
            i -= kSubSize;
        }
        return kOne << i;
    }
    
    static inline size_t nonzeroBits(SubBitmask subMask) {
        return __builtin_popcountll(subMask);
    }

public:
    Bitmask() {
        mask.fill(0);
    }
    
    bool get(int i) const {
        debug_assert(i < kSubSize * N && i >= 0);
        return mask[i >> kLogSubSize] & bitFor(i);
    }

    void set(int i) {
        mask[i >> kLogSubSize] |= bitFor(i);
    }

    void clear(int i) {
        mask[i >> kLogSubSize] &= ~bitFor(i);
    }

    int countBefore(int i) const {
        int count = 0;
        for (int j = 0; j < i >> kLogSubSize; ++j) {
            count += nonzeroBits(mask[j]);
        }
        auto firstDigitsMask = bitFor(i) - 1;
        count += nonzeroBits(mask[i >> kLogSubSize] & firstDigitsMask);
        return count;
    }
    
    void print() const {
        for (int i = 0; i < N; ++i) {
            std::cout << mask[i] << " ";
        }
        std::cout << std::endl;
    }

private:
    std::array<SubBitmask, N> mask;
};
