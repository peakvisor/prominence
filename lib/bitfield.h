#pragma once

#include <deque>
#include <vector>
#include <cstdint>
#include <fstream>

template <uint32_t Side>
class BoolArrayField {
public:
    BoolArrayField() {
        v.fill(false);
    }
    
    bool get(uint32_t i, uint32_t j) const {
        return v[i + j * Side];
    }
    
    uint32_t getCount() const {
        return count;
    }
    
    void set(uint32_t i, uint32_t j, bool value) {
        bool previousValue = v[i + j * Side];
        if (previousValue != value) {
            count += value ? 1 : -1;
            v[i + j * Side] = value;
        }
    }
    
    std::deque<uint32_t> encodeRunLength() const {
        bool seeking = true;
        std::deque<uint32_t> result;
        uint32_t runLength = 0;
        for (const auto x : v) {
            if (seeking == x) {
                result.emplace_back(runLength);
                runLength = 0;
                seeking = !seeking;
            } else {
                ++runLength;
            }
        }
        return result;
    }
    
    void clear() {
        v.fill(false);
    }
    
private:
    std::array<bool, Side * Side> v;
    uint32_t count = 0;
};

class BoolDequeField {
public:
    BoolDequeField(PlanarIndex from, PlanarIndex to) : origin(from), size(to - from) {
        std::cout << static_cast<size_t>(size.x) * static_cast<size_t>(size.y) << endl;
        v.resize(static_cast<size_t>(size.x) * static_cast<size_t>(size.y));
    }
    
    inline bool get(const PlanarIndex &planar) const {
        auto offset = planar - origin;
        return getFromOffset(offset);
    }
    
    inline bool getFromOffset(const PlanarIndex &offset) const {
        return v[offset.linearizeBig(size.x)];
    }
    
    inline void set(const PlanarIndex &planar, bool value) {
        auto offset = planar - origin;
        v[offset.linearizeBig(size.x)] = value;
    }
    
    void clear() {
        v.clear();
    }
    
private:
    const PlanarIndex origin; // TODO: optimize for origin = 0, 0 if nec
    const PlanarIndex size;
    std::deque<bool> v;
};

template <uint8_t radius>
class VicinityBitfield {
public:
    void set(const PlanarIndex &offset) {
        assert(isInside(offset));
        _field.set(offset.x + radius, offset.y + radius, true);
    }
    
    bool get(const PlanarIndex &offset) {
        assert(isInside(offset));
        return _field.get(offset.x + radius, offset.y + radius);
    }
    
    bool isInside(const PlanarIndex &offset) const {
        return std::abs(offset.x) <= radius && std::abs(offset.y) <= radius ;
    }
    
    void clear() {
        _field.clear();
    }
    
private:
    BoolArrayField<2 * radius + 1> _field;
};
